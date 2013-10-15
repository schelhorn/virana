#!/usr/bin/env python
#from __future__ import print_function

#import cProfile
#import line_profiler
import numpy

import sys
import re

import tempfile
import subprocess
import shutil
import os
import os.path
import logging
import bz2
import zlib

import math
import string

from collections import defaultdict

from subprocess import PIPE
import time



NON_ID = ''.join(c for c in map(chr, range(256)) if not c.isalnum())
NON_ID = NON_ID.replace('_', '').replace('-', '')

try:
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    message = 'This script requires the BioPython python package\n'
    sys.stderr.write(message)
    sys.exit(1)

try:
    from plumbum import cli
except ImportError:
    message = 'This script requires the plumbum python package\n'
    sys.stderr.write(message)
    sys.exit(1)


try:
    import HTSeq
except ImportError:
    message = 'This script requires the HTSeq python package\n'
    sys.stderr.write(message)
    sys.exit(1)

KHMER_AVAILABLE = True
try:
    import khmer
except ImportError:
    KHMER_AVAILABLE = False

#from io import BufferedRandom

logging.basicConfig(level=logging.INFO, format='%(message)s')


# def profile_this(fn):
#     def profiled_fn(*args, **kwargs):
#         fpath = fn.__name__ + ".profile"
#         prof = cProfile.Profile()
#         ret = prof.runcall(fn, *args, **kwargs)
#         prof.dump_stats(fpath)
#         return ret
#     return profiled_fn

from operator import itemgetter
from heapq import nlargest
from itertools import repeat, ifilter

class Counter(dict):
    '''Dict subclass for counting hashable objects.  Sometimes called a bag
    or multiset.  Elements are stored as dictionary keys and their counts
    are stored as dictionary values.

    >>> Counter('zyzygy')
    Counter({'y': 3, 'z': 2, 'g': 1})

    '''

    def __init__(self, iterable=None, **kwds):
        '''Create a new, empty Counter object.  And if given, count elements
        from an input iterable.  Or, initialize the count from another mapping
        of elements to their counts.

        >>> c = Counter()                           # a new, empty counter
        >>> c = Counter('gallahad')                 # a new counter from an iterable
        >>> c = Counter({'a': 4, 'b': 2})           # a new counter from a mapping
        >>> c = Counter(a=4, b=2)                   # a new counter from keyword args

        '''
        self.update(iterable, **kwds)

    def __missing__(self, key):
        return 0

    def most_common(self, n=None):
        '''List the n most common elements and their counts from the most
        common to the least.  If n is None, then list all element counts.

        >>> Counter('abracadabra').most_common(3)
        [('a', 5), ('r', 2), ('b', 2)]

        '''
        if n is None:
            return sorted(self.iteritems(), key=itemgetter(1), reverse=True)
        return nlargest(n, self.iteritems(), key=itemgetter(1))

    def elements(self):
        '''Iterator over elements repeating each as many times as its count.

        >>> c = Counter('ABCABC')
        >>> sorted(c.elements())
        ['A', 'A', 'B', 'B', 'C', 'C']

        If an element's count has been set to zero or is a negative number,
        elements() will ignore it.

        '''
        for elem, count in self.iteritems():
            for _ in repeat(None, count):
                yield elem

    # Override dict methods where the meaning changes for Counter objects.

    @classmethod
    def fromkeys(cls, iterable, v=None):
        raise NotImplementedError(
            'Counter.fromkeys() is undefined.  Use Counter(iterable) instead.')

    def update(self, iterable=None, **kwds):
        '''Like dict.update() but add counts instead of replacing them.

        Source can be an iterable, a dictionary, or another Counter instance.

        >>> c = Counter('which')
        >>> c.update('witch')           # add elements from another iterable
        >>> d = Counter('watch')
        >>> c.update(d)                 # add elements from another counter
        >>> c['h']                      # four 'h' in which, witch, and watch
        4

        '''
        if iterable is not None:
            if hasattr(iterable, 'iteritems'):
                if self:
                    self_get = self.get
                    for elem, count in iterable.iteritems():
                        self[elem] = self_get(elem, 0) + count
                else:
                    dict.update(self, iterable) # fast path when counter is empty
            else:
                self_get = self.get
                for elem in iterable:
                    self[elem] = self_get(elem, 0) + 1
        if kwds:
            self.update(kwds)

    def copy(self):
        'Like dict.copy() but returns a Counter instance instead of a dict.'
        return Counter(self)

    def __delitem__(self, elem):
        'Like dict.__delitem__() but does not raise KeyError for missing values.'
        if elem in self:
            dict.__delitem__(self, elem)

    def __repr__(self):
        if not self:
            return '%s()' % self.__class__.__name__
        items = ', '.join(map('%r: %r'.__mod__, self.most_common()))
        return '%s({%s})' % (self.__class__.__name__, items)

    # Multiset-style mathematical operations discussed in:
    #       Knuth TAOCP Volume II section 4.6.3 exercise 19
    #       and at http://en.wikipedia.org/wiki/Multiset
    #
    # Outputs guaranteed to only include positive counts.
    #
    # To strip negative and zero counts, add-in an empty counter:
    #       c += Counter()

    def __add__(self, other):
        '''Add counts from two counters.

        >>> Counter('abbb') + Counter('bcc')
        Counter({'b': 4, 'c': 2, 'a': 1})


        '''
        if not isinstance(other, Counter):
            return NotImplemented
        result = Counter()
        for elem in set(self) | set(other):
            newcount = self[elem] + other[elem]
            if newcount > 0:
                result[elem] = newcount
        return result

    def __sub__(self, other):
        ''' Subtract count, but keep only results with positive counts.

        >>> Counter('abbbc') - Counter('bccd')
        Counter({'b': 2, 'a': 1})

        '''
        if not isinstance(other, Counter):
            return NotImplemented
        result = Counter()
        for elem in set(self) | set(other):
            newcount = self[elem] - other[elem]
            if newcount > 0:
                result[elem] = newcount
        return result

    def __or__(self, other):
        '''Union is the maximum of value in either of the input counters.

        >>> Counter('abbb') | Counter('bcc')
        Counter({'b': 3, 'c': 2, 'a': 1})

        '''
        if not isinstance(other, Counter):
            return NotImplemented
        _max = max
        result = Counter()
        for elem in set(self) | set(other):
            newcount = _max(self[elem], other[elem])
            if newcount > 0:
                result[elem] = newcount
        return result

    def __and__(self, other):
        ''' Intersection is the minimum of corresponding counts.

        >>> Counter('abbb') & Counter('bcc')
        Counter({'b': 1})

        '''
        if not isinstance(other, Counter):
            return NotImplemented
        _min = min
        result = Counter()
        if len(self) < len(other):
            self, other = other, self
        for elem in ifilter(self.__contains__, other):
            newcount = _min(self[elem], other[elem])
            if newcount > 0:
                result[elem] = newcount
        return result

class CLI(cli.Application):
    """RNA-Seq and DNA-Seq short read analysis by mapping to known reference sequences."""
    PROGNAME = "vmap"
    VERSION = "1.0.0"
    DESCRIPTION = \
"""DESCRIPTION: virana vmap - short read mapping for clinical metagenomics.

The virana mapping utility ('vmap') is a wrapper around indexing and mapping
fascilities of three short read alignments tools, STAR, BWA-MEM, and SMALT.
vmap is able to generate reference indexes for each of these mappers from a
single FASTA file and then map short reads in FASTQ format against these indexes.
During mappings, alignments are output as SAM-files, unsorted BAM-files,
taxonomic bins (by taxonomic families), and a special virana format that
summarizes reads that align to specific taxonomic reference databases such as
viruses and additionaly capture multimapping information. In addition, some mappers
also supoort output of unmapped reads in FASTQ format for later assembly or
special treatment of chimeric alignments and their output in SAM files.

https://github.com/schelhorn/virana

Schelhorn S-E, Fischer M, Tolosi L, Altmueller J, Nuernberg P, et al. (2013)
Sensitive Detection of Viral Transcripts in Human Tumor Transcriptomes.
PLoS Comput Biol 9(10): e1003228. doi:10.1371/journal.pcbi.1003228"""

    USAGE = """USAGE: The program has four modes that can be accessed by
       [vmap | python vmap.py] [rnaindex | dnaindex | rnamap | dnamap]"""

    def main(self, *args):

        print 'CLI main'

        if args:
            print self.USAGE
            print("ERROR: Unknown command %r" % (args[0]))
            return 1

        if not self.nested_command:
            print self.USAGE
            print("ERROR : No command given")
            return 1


class SAMHits:
    """ Converts SAM output of mappers into bzipped HIT files. """

    def __init__(self, output_file, sample_id, refseq_filter=None, min_mapping_score=None,\
                     min_alignment_score=None, max_mismatches=None,\
                     max_relative_mismatches=None, min_continiously_matching=None,\
                     filter_complexity=False):

        self.output_file                = bz2.BZ2File(output_file, 'wb', buffering=100 * 1024 * 1024)
        self.sample_id                  = sample_id.translate(None, NON_ID)
        self.refseq_filter              = refseq_filter
        self.max_mismatches             = max_mismatches
        self.max_relative_mismatches    = max_relative_mismatches
        self.current_group              = []
        self.min_mapping_score          = min_mapping_score
        self.min_alignment_score        = min_alignment_score
        self.min_continiously_matching  = min_continiously_matching
        self.filter_complexity          = filter_complexity

        self.re_matches = re.compile(r'(\d+)M')
        self.re_dels    = re.compile(r'(\d+)D')

    def count(self, parsed_line):

        if parsed_line is None:
            return

        read_key, read_name, flag, ref_name, ref_position, mapping_score,\
                    cigar, mate_ref_name, mate_ref_position, insert_size, seq, qual,\
                    is_end1, is_end2, number_mismatches, alignment_score,\
                    number_hits, is_reverse, is_primary, is_mapped, is_mate_mapped,\
                    is_paired, number_matches, read_end_pos, max_match = parsed_line

        if not is_mapped:
            return

        if self.filter_complexity:

            avg_compression = float(len(zlib.compress(seq)))/len(seq)

            if avg_compression < 0.5:
                return

            # length = len(seq)
            # counts = [seq.count(nuc) for nuc in 'ACGT']
            # min_count = length * 0.10
            # max_count = length * 0.50
            # for count in counts:
            #     if count < min_count or count > max_count:
            #         return None

            # counter = Counter()
            # for i in range(length - 2):
            #     counter[seq[i: i + 3]] += 1
            # maximal = length - 4

            # highest = sum([v for k, v in counter.most_common(2)])
            # if highest > (maximal / 3.0):
            #     return None

            # self.passed.append(avg_compression)

        pair_id = ''
        if is_end1:
            pair_id = '/1'

        elif is_end2:
            pair_id = '/2'

        read_name = self.sample_id + ';' + read_name + pair_id

        # Initialize new current group
        if len(self.current_group) == 0:
            self.current_group = [read_name, seq, []]

        # Write old current group to file
        if read_name != self.current_group[0]:
            self._write_group()
            self.current_group = [read_name, seq, []]

        try:
            refseq_group, family, organism, identifier = ref_name.split(';')[:4]
        except ValueError:
            sys.stderr.write('Read mapped to malformed reference sequence %s, skipping\n' % ref_name)
            return

        if self.min_continiously_matching:

            if self.min_continiously_matching > max_match:
                return

        if self.max_mismatches\
            and int(number_mismatches) > self.max_mismatches:
            return

        if self.max_relative_mismatches\
            and int(number_mismatches) / float(len(seq))\
            > self.max_relative_mismatches:
            return

        if self.min_mapping_score\
            and self.min_mapping_score > mapping_score:
            return

        if self.min_alignment_score\
            and self.min_alignment_score > alignment_score:
            return

        start = int(ref_position) + 1

        self.current_group[2].append([refseq_group, family, organism, identifier, str(start), str(read_end_pos)])

    def _write_group(self):
        passed = True

        if self.refseq_filter:
            passed = False
            for refseq_group, family, organism, identifier, start, end in self.current_group[2]:
                if passed:
                    break
                for f in self.refseq_filter:
                    if refseq_group == f:
                        passed = True
                        break
        if passed:
            description = []
            for identifier in self.current_group[2]:
                description.append(';'.join(identifier))
            description = '|'.join(description)

            record = SeqRecord(Seq(self.current_group[1]))
            record.id = 'Read;' + self.current_group[0]
            record.description = description

            SeqIO.write([record], self.output_file, "fasta")

    def write(self):

        self._write_group()
        self.output_file.close()

class SAMParser:

    def parse(self, line):

        if line[0] == '@':
            return None

        alignment   = HTSeq._HTSeq.SAM_Alignment.from_SAM_line(line)
        read_name   = alignment.read.name
        seq         = alignment.read.seq
        qual        = numpy.array(alignment.read.qual)
        flag        = alignment.flag
        cigar       = None

        is_paired       = (flag & 1)
        is_mapped       = not (flag & 4)
        is_mate_mapped  = alignment.mate_aligned is not None #not (flag & 8)
        is_reverse      = (flag & 16)
        is_end1         = (flag & 64)
        is_end2         = (flag & 128)
        is_primary      = not (flag & 256)

        read_key = (read_name, is_end1)

        ref_name            = None
        ref_position        = None
        mapping_score       = 0
        mate_ref_name       = None
        mate_ref_position   = None
        insert_size         = None
        alignment_score     = 0
        read_end_pos        = None

        if is_mate_mapped and alignment.mate_start:
            mate_ref_name       = alignment.mate_start.chrom
            mate_ref_position   = alignment.mate_start.start

        number_hits         = 0
        alignment_score     = 0
        number_mismatches   = 0

        number_matches      = 0
        max_match           = 0

        if is_mapped:

            ref_name            = alignment.iv.chrom
            ref_position        = alignment.iv.start
            read_end_pos        = alignment.iv.end
            alignment_score     = alignment.aQual
            cigar               = alignment.cigar

            if is_mate_mapped:
                insert_size         = alignment.inferred_insert_size

            for c in cigar:
                if c.type == 'M':
                    number_matches += c.size
                    max_match = max(max_match, c.size)

            for tag, value in alignment.optional_fields:
                if tag == 'NM':
                    number_hits = value
                elif tag == 'AS':
                    alignment_score = value
                elif tag == 'NH':
                    number_mismatches = value

        return read_key, read_name, flag, ref_name, ref_position, mapping_score,\
            cigar, mate_ref_name, mate_ref_position, insert_size, seq, qual,\
            is_end1, is_end2, number_mismatches, alignment_score,\
            number_hits, is_reverse, is_primary, is_mapped, is_mate_mapped,\
            is_paired, number_matches, read_end_pos, max_match


class SAMQuality:

    def __init__(self, file_path):

        self.file_path          = file_path

        self.stored             = defaultdict(Counter)
        self.all_references     = defaultdict(int)
        self.primary_references = defaultdict(int)
        self.complement         = string.maketrans('ATCGN', 'TAGCN')

        if KHMER_AVAILABLE:
            self.ktable = khmer.new_ktable(10)

    def _get_complement(self, sequence):

        return sequence.translate(self.complement)[::-1]

    def _get_summary(self, counter):
        """"Returns five numbers (sum, extrema, mean, and std)
            for a max_frequency counter """

        maximum  = 0
        minimum  = sys.maxint
        thesum   = 0
        allcount = 0
        mode     = [0, None]

        items       = 0.0
        mean        = 0.0
        m2          = 0.0
        variance    = 0.0

        for item in counter:

            count   = counter[item]
            if count > mode[0]:
                mode = [count, item]

            allcount += count
            maximum = max(maximum, item)
            minimum = min(minimum, item)
            thesum  += (count * item)

            x = 1
            while x <= count:
                items       += 1
                delta       = item - mean
                mean        = mean + delta / items
                m2          = m2 + delta * (item - mean)
                variance    = m2 / items
                x += 1

        std = math.sqrt(variance)

        return allcount, thesum, minimum, maximum, mode[1], mean, std


    def _to_unit(self, item, is_percentage=False):
        """ Convert a numeric to a string with metric units """

        if is_percentage:
            return ('%-.3f' % (item * 100)) + '%'
        converted = None
        try:
            item = float(item)
            if item > 10**12:
                converted = str(round(item / 10**9,3))+'P'
            elif item > 10**9:
                converted = str(round(item / 10**9,3))+'G'
            elif item > 10**6:
                converted = str(round(item / 10**6,3))+'M'
            elif item > 10**3:
                converted = str(round(item / 10**3,3))+'K'
            else:
                converted = str(round(item,3))
        except:
            converted = str(item)

        return converted

    def _str_metrics(self, data):

        str_metrics = []

        for (item, metric) in sorted(data.keys()):
            counter = data[(item, metric)]
            if not hasattr(counter.iterkeys().next(), 'real'):
                for element, count in counter.most_common():
                    str_metrics.append(self._str_metric(item, metric, element, count, no_units=True))
            else:
                summary = self._get_summary(counter)
                str_metrics.append(self._str_metric(item, metric, *summary))

        return str_metrics

    def _str_metric(self, item, metric, count, thesum='', minimum='',\
                        maximum='', mode='', mean='', std='', no_units=False):

        counters = [count, thesum, minimum, maximum, mode, mean, std]
        counters = map(str, counters)

        if no_units:
            items = [item, metric] +  counters
        else:
            units = map(self._to_unit, counters)
            items = [item, metric] + units

        return '%-15s\t%-60s\t%12s\t%12s\t%12s\t%12s\t%12s\t%12s\t%12s\n' \
                    % tuple(items)


    def _count_read(self, metric, data, sample):

        item = 'read'

        (insert_size, alignment_score, mapping_score, length,\
            q20_length, avg_phred_quality, number_hits, is_reverse) = data

        self.stored[(item, metric + ' mappings')][number_hits] += sample
        self.stored[(item, metric + ' insert')][insert_size] += sample


    def _count_segment(self, metric, data, sample):

        item = 'segment'

        (insert_size, alignment_score, mapping_score, length,\
            q20_length, avg_phred_quality, number_hits, is_reverse) = data

        self.stored[(item, metric + ' algq')][alignment_score] += sample
        self.stored[(item, metric + ' mapq')][mapping_score] += sample
        self.stored[(item, metric + ' length')][length] += sample
        self.stored[(item, metric + ' q20length')][q20_length] += sample
        self.stored[(item, metric + ' meanbasequal')][avg_phred_quality] += sample
        self.stored[(item, metric + ' reverse')][is_reverse] += sample

    def count(self, parsed_line):

        if parsed_line is None:
            return

        #print_metric('Item' , 'Metric', 'Count', 'Sum', 'Min', 'Max', 'Mode', 'Mean', 'STD')

        read_key, read_name, flag, ref_name, ref_position, mapping_score,\
                    cigar, mate_ref_name, mate_ref_position, insert_size, seq, qual,\
                    is_end1, is_end2, number_mismatches, alignment_score,\
                    number_hits, is_reverse, is_primary, is_mapped, is_mate_mapped,\
                    is_paired, number_matches, read_end_pos, max_match = parsed_line

        phred_quality       = qual - 33
        avg_phred_quality   = numpy.mean(phred_quality)
        length              = len(seq)
        mate_reference_id   = mate_ref_name
        reference_id        = ref_name
        reference           = reference_id is not None and reference_id != '*'
        insert_size         = insert_size and abs(insert_size) or insert_size
        is_segment1         = not is_paired or (is_paired and is_end1)
        is_reverse          = is_reverse
        is_unique           = is_primary and number_hits == 1
        is_translocation    = is_paired and is_mapped and is_mate_mapped\
                                    and (mate_reference_id != '=' and reference_id != mate_reference_id)
        is_part_of_doublemap    = is_paired and is_mapped and is_mate_mapped
        is_part_of_halfmap      = is_paired and (is_mapped != is_mate_mapped)
        is_part_of_nomap        = is_paired and not is_mapped and not is_mate_mapped


        # Count length until first low quality base call
        q20_length = 0
        for q in phred_quality:
            if q < 20:
                break
            q20_length += 1

        # Count kmers
        if KHMER_AVAILABLE:
            if not is_reverse:
                self.ktable.consume(seq)
            else:
                self.ktable.consume(self._get_complement(seq))

        if reference:

            self.all_references[reference_id] += 1

            if is_primary:
                self.primary_references[reference_id] += 1


        data   = (insert_size, alignment_score, mapping_score, length,\
                        q20_length, avg_phred_quality, number_hits, is_reverse)

        sample = 1

        self._count_segment('sequenced', data, sample)
        if is_mapped:
            self._count_segment('sequenced mapped multi', data, sample)
            if is_primary:
                self._count_segment('sequenced mapped primary', data, sample)
                if number_hits and is_unique:
                    self._count_segment('sequenced mapped primary unique', data, sample)

            if is_segment1:
                self._count_read('sequenced mapped multi', data, sample)
                if is_primary:
                    self._count_read('sequenced mapped primary', data, sample)

        if is_paired:
            self._count_segment('sequenced paired', data, sample)

            if is_part_of_doublemap:
                self._count_segment('sequenced paired doublemap', data, sample)

                if is_primary:
                    self._count_segment('sequenced paired doublemap primary', data, sample)

                if is_segment1:
                    self._count_read('sequenced paired doublemap multi', data, sample)

                    if is_primary:
                        self._count_read('sequenced paired doublemap primary', data, sample)

                        if number_hits and is_unique:
                            self._count_read('sequenced paired doublemap primary unique', data, sample)

                            if is_translocation:
                                self._count_read('sequenced paired doublemap primary unique translocation', data, sample)

            elif is_part_of_halfmap:

                self._count_segment('sequenced paired halfmap', data, sample)

                # The mapped segment
                if is_mapped:

                    self._count_segment('sequenced paired halfmap mapped', data, sample)

                    if is_primary:

                        self._count_read('sequenced paired halfmap mapped primary', data, sample)

                        if number_hits and is_unique:
                            self._count_read('sequenced paired halfmap mapped primary unique', data, sample)

                    elif not is_primary:

                        self._count_read('sequenced unpaired mapped multi', data, sample)

                # The unmapped segment
                elif not is_mapped:
                    self._count_segment('sequenced paired halfmap unmapped', data, sample)

            elif is_part_of_nomap:

                self._count_segment('sequenced paired nomap', data, sample)

                if is_segment1:

                    self._count_read('sequenced paired nomap', data, sample)

        elif not is_paired:

            self._count_segment('sequenced unpaired', data, sample)

            if is_mapped:

                self._count_segment('sequenced unpaired mapped', data, sample)

                if is_primary:

                        self._count_read('sequenced unpaired mapped primary', data, sample)

                        if number_hits and is_unique:

                            self._count_read('sequenced paired unpaired mapped primary unique', data, sample)

                elif not is_primary:

                    self._count_read('sequenced unpaired mapped multi', data, sample)


            elif not is_mapped:

                self._count_segment('sequenced unpaired unmapped', data, sample)

                if is_segment1:
                    self._count_read('sequenced unpaired unmapped', data, sample)

    def write(self):

        with open(self.file_path, 'w') as output_file:

            all_references = sorted([(count, reference) for reference, count\
                                    in self.all_references.iteritems()], reverse=True)

            for j, (count, reference) in enumerate(all_references[:30]):
                self.stored[('segment', 'multireference_' + str(j+1))][reference] = count

            primary_references = sorted([(count, reference) for reference, count\
                                            in self.primary_references.iteritems()], reverse=True)

            for j, (count, reference) in enumerate(primary_references[:30]):
                self.stored[('segment', 'primaryreference_' + str(j+1))][reference] = count

            # Extract top-ranking kmers
            if KHMER_AVAILABLE:
                kmer_frequencies = []
                for i in range(0, self.ktable.n_entries()):
                    n = self.ktable.get(i)
                    if n > 0:
                        kmer_frequencies.append((n, self.ktable.reverse_hash(i)))
                kmer_frequencies = sorted(kmer_frequencies, reverse=True)
                for j, (frequency, kmer) in enumerate(kmer_frequencies[:10]):
                    self.stored[('segment', 'kmer_' + str(j+1))][kmer] = frequency

            output_file.writelines(self._str_metrics(self.stored))


class SAMTaxonomy:
    """ Provides taxonomic summary information from a SAM file stream. """

    def __init__(self, file_path):

        self.file_path = file_path

        self.count_primaries = Counter()
        self.detailed_information = {}

        self._last_read                 = (None, None)
        self._last_read_human_prim      = 0
        self._last_read_human_sec       = 0
        self._last_organisms            = set()

    def count(self, parsed_line):

        if parsed_line is None:
            return

        read_key, read_name, flag, ref_name, ref_position, mapping_score,\
                    cigar, mate_ref_name, mate_ref_position, insert_size, seq, qual,\
                    is_end1, is_end2, number_mismatches, alignment_score,\
                    number_hits, is_reverse, is_primary, is_mapped, is_mate_mapped,\
                    is_paired, number_matches, read_end_pos, max_match = parsed_line

        if is_mapped:

            refseq_group, family, organism, gi = ref_name.split(';')[:4]

            if is_primary:
                self.count_primaries[organism] += 1

            if organism not in self.detailed_information:
                # refseq_group. family, gis, avg_mapping_score,
                # avg_seq_length, avg_number_hits, avg_alignment_score, avg_nr_mismatches
                initial = [refseq_group,
                           family,
                           set([gi]),
                           [int(mapping_score), 1],
                           [len(seq), 1],
                           [0, 0],
                           [alignment_score, 1],
                           [number_mismatches, 1],
                           0,
                           0]

                self.detailed_information[organism] = initial

            else:
                entry = self.detailed_information[organism]
                entry[2].add(gi)
                entry[3][0] += int(mapping_score)
                entry[3][1] += 1
                entry[4][0] += len(seq)
                entry[4][1] += 1
                entry[6][0] += alignment_score
                entry[6][1] += 1
                entry[7][0] += number_mismatches
                entry[7][1] += 1

            if is_primary:
                entry = self.detailed_information[organism]
                entry[5][0] += number_hits
                entry[5][1] += 1

            if self._last_read == (None, None):
                self._last_read = read_key

            if self._last_read != read_key:

                for last_organism in self._last_organisms:

                    self.detailed_information[last_organism][8]\
                        += self._last_read_human_prim

                    self.detailed_information[last_organism][9]\
                        += self._last_read_human_sec

                self._last_read = read_key
                self._last_organisms = set()
                self._last_read_human_prim  = 0
                self._last_read_human_sec   = 0

            self._last_organisms.add(organism)

            if organism == 'Homo_sapiens':
                if is_primary:
                    self._last_read_human_prim += 1
                else:
                    self._last_read_human_sec += 1

    def get_summary(self, top=100):

        lines = []

        lines.append('%10s\t%20s\t%20s\t%-20s\t%10s\t%10s\t%10s\t%5s\t%5s\t%5s\t%10s\t%10s\n'\
            % ('Count', 'Group', 'Family', 'Organism', 'Targets', 'ReadLen', 'Hits', 'Map', 'Algn', 'Mism', 'HuP', 'HuS'))

        top_organisms = self.count_primaries.most_common(top)

        for organism, count in top_organisms:

            refseq_group, family, identifiers,\
                avg_mapping_score, avg_seq_length, avg_number_hits,\
                avg_alignment_score, avg_nr_mismatches, human_prim, human_sec\
                = self.detailed_information[organism]

            avg_len = int(avg_seq_length[0] / float(avg_seq_length[1]))
            if avg_number_hits[1] == 0:
                avg_hits = 0
            else:
                avg_hits = int(avg_number_hits[
                               0] / float(avg_number_hits[1]))

            avg_mapping_score = int(avg_mapping_score[
                                    0] / float(avg_mapping_score[1]))

            avg_alignment_score = int(avg_alignment_score[
                                      0] / float(avg_alignment_score[1]))

            avg_nr_mismatches = int(avg_nr_mismatches[
                                    0] / float(avg_nr_mismatches[1]))

            nr_ids = len(identifiers)

            if count > 10**6:
                count = str(round(count / float(10**6), 3)) + 'M'
            if human_prim > 10**6:
                human_prim = str(round(human_prim / float(10**6), 3)) + 'M'
            if human_sec > 10**6:
                human_sec = str(round(human_sec / float(10**6), 3)) + 'M'
            if nr_ids > 10**6:
                nr_ids = str(round(nr_ids / float(10**6), 3)) + 'M'

            lines.append('%10s\t%20s\t%20s\t%-20s\t%10s\t%10i\t%10i\t%5i\t%5i\t%5i\t%10s\t%10s\n'\
                % (str(count), refseq_group[:20], family[:20], organism[:20],\
                    str(nr_ids), avg_len, avg_hits, avg_mapping_score,\
                    avg_alignment_score, avg_nr_mismatches, str(human_prim),\
                    str(human_sec)))

        return lines

    def write(self):

        with open(self.file_path, 'w') as output_file:
            output_file.writelines(self.get_summary())


class Index(cli.Application):

    def setup_logging(self):
        if self.debug:
            logging.getLogger().setLevel(logging.DEBUG)

        return self.debug

    def validate_paths(self):

        # Check if genome directory is existing
        if not os.path.exists(self.reference_file):
            sys.stderr.write(
                'Reference file %s nor existing, exiting' % self.reference_file)
            sys.exit(1)

        # Check if output directory is existing
        if not os.path.exists(self.index_dir):
            logging.debug(
                'Making output directory for index at %s' % self.index_dir)
            os.makedirs(self.index_dir)

    def get_command_line(self):

        # # Make named pipe to extract genomes
        # pipe_path = os.path.abspath(os.path.join(self.genome_dir, 'pipe.fa'))
        # if os.path.exists(pipe_path):
        #     os.unlink(pipe_path)
        # os.mkfifo(pipe_path)

        pass

    def run_index_process(self):

        # Run index generation process
        command_line = self.get_command_line()
        process = subprocess.Popen(' '.join(command_line), shell=True, stdout=PIPE, stderr=PIPE)

        # Block until streams are closed by the process
        stdout, stderr = process.communicate()

        if stderr:
            sys.stderr.write(stderr)

        if self.debug and stdout:
            print stdout

    def main(self, *args):

        self.setup_logging()
        self.validate_paths()
        self.run_index_process()

@CLI.subcommand("rnaindex")
class RNAIndex(Index):
    """ Creates a STAR index from a FASTA genome reference """

    reference_file = cli.SwitchAttr(
        ['-r', '--reference_file'], str, mandatory=True,
        help="Sets the reference genome FASTA file.")

    index_dir = cli.SwitchAttr(['-i', '--index_dir'], str, mandatory=True,
                               help="Sets the index output directory." +
                               " Directory will be generated if not existing." +
                               " Directory will be filled with several index files.")
    threads = cli.SwitchAttr(
        ['-t', '--threads'], cli.Range(1, 512), mandatory=False,
        help="Sets the number of threads to use",
        default=1)

    max_ram = cli.SwitchAttr(
        ['-m'], cli.Range(1, 400000000000), mandatory=False,
        help="Sets the maximum amount of memory (RAM) to use (in bytes)",
        default=400000000000)

    path = cli.SwitchAttr(['-p', '--path'], str, mandatory=False,
                          help="Path to STAR executable",
                          default='STAR')
    sparse = cli.Flag(
        ["-s", "--sparse"], help="If given, a sparse index that requires less " +
        " RAM in the mapping phase will be constructed")

    debug = cli.Flag(["-d", "--debug"], help="Enable debug output")

    def get_command_line(self):

        # Make star command line
        cline = [self.path] + ['--runMode', 'genomeGenerate',
                        '--genomeDir', self.index_dir,
                        '--limitGenomeGenerateRAM', str(self.max_ram),
                        '--runThreadN', str(self.threads),
                        '--genomeFastaFiles'] + self.reference_files

        # Add parameters for sparse (memory-saving) index generation
        if self.sparse:
            cline += ['--genomeSAsparseD', '2',
                      '--genomeChrBinNbits', '12',
                      '--genomeSAindexNbases', '13']

        else:
            cline += ['--genomeSAsparseD', '1',
                      '--genomeChrBinNbits', '18',
                      '--genomeSAindexNbases', '15']

        if self.debug:
            print ' '.join(cline)

        return cline

@CLI.subcommand("dnaindex")
class DNAIndex(Index):
    """ Creates a BWA index from a FASTA reference file """

    reference_file  = cli.SwitchAttr(['-r', '--reference_file'], str, mandatory=True,
                                 help="Sets the input reference FASTA file.")
    index_dir   = cli.SwitchAttr(['-i', '--index_dir'], str, mandatory=True,
                                 help="Sets the index output directory." +
                                 " Directory will be generated if not existing." +
                                 " Directory will be filled with several index files.")
    path        = cli.SwitchAttr(['-p', '--path'], str, mandatory=False,
                                 help="Path to BWA executable",
                                 default='bwa')
    debug       = cli.Flag(["-d", "--debug"], help="Enable debug output")

    def get_command_line(self):

        # Make star command line
        cline = [self.path] + ['index', '-a', 'bwtsw', '-p', os.path.join(self.index_dir, 'index'), self.reference_file]

        return cline

@CLI.subcommand("varindex")
class VARIndex(Index):
    """ Creates a SMALT index from a FASTA reference file """

    reference_file  = cli.SwitchAttr(['-r', '--reference_file'], str, mandatory=True,
                                 help="Sets the input reference FASTA file.")
    index_dir   = cli.SwitchAttr(['-i', '--index_dir'], str, mandatory=True,
                                 help="Sets the index output directory." +
                                 " Directory will be generated if not existing." +
                                 " Directory will be filled with several index files.")
    path        = cli.SwitchAttr(['-p', '--path'], str, mandatory=False,
                                 help="Path to SMALT executable",
                                 default='smalt_x86_64')
    debug       = cli.Flag(["-d", "--debug"], help="Enable debug output")

    def get_command_line(self):

        # Make star command line
        cline = [self.path] + ['index', '-k', '20', '-s', '2', os.path.join(self.index_dir, 'index'), self.reference_file]

        return cline

class Mapper(cli.Application):

    def setup_logging(self):
        if self.debug:
            logging.getLogger().setLevel(logging.DEBUG)

        return self.debug

    def validate_paths(self, mapper_executable):

        mapper_path      = self.mapper_path and self.mapper_path or mapper_executable
        samtools_path    = self.samtools_path and self.samtools_path or 'samtools'

        # Check if genome directory is existing
        if not os.path.exists(self.index_dir):
            message = 'Index directory %s not existing, exiting' % self.index_dir
            sys.stderr.write(message)
            sys.exit(1)

        # Check for number of read input files
        if len(self.reads) not in (1, 2):
            message = 'Invalid number of FASTQ files; supply either one (single end) or two (paired end)\n'
            sys.stderr.write(message)
            sys.exit(1)

        # Make temporary directories
        if self.temp_path:
            temp_path = tempfile.mkdtemp(dir=self.temp_path)
        else:
            temp_path = tempfile.mkdtemp()

        # Try if we can make the relevant output files
        outputs = ['unmapped1', 'unmapped2', 'taxonomy', 'qual', 'hits', 'sam', 'bam', 'chimeric_mappings']
        for output in outputs:
            try:
                attribute = getattr(self, output)
            except AttributeError:
                continue

            if attribute is None or attribute == '':
                continue
            try:
                with file(attribute, 'a'):
                    os.utime(attribute, None)
            except IOError:
                sys.stderr.write('Could not write output file %s\n' % attribute)
                sys.exit(1)

        return mapper_path, samtools_path, temp_path

    def get_command_line(self, mapper_path, temp_path):
        pass

    def run_mapper_process(self, mapper_path, temp_path):

        command_line = ' '.join(self.get_command_line(mapper_path, temp_path))
        process = subprocess.Popen(command_line, shell=True, stdout=PIPE)

        if self.debug:
            print 'Running: %s' % command_line

        return process

    def setup_outputs(self, samtools_path):

        taxonomy = None
        if self.taxonomy:
            taxonomy = SAMTaxonomy(self.taxonomy)

        quality = None
        if self.qual:
            quality = SAMQuality(self.qual)

        hits = None
        if self.hits:
            hits = SAMHits(self.hits, self.sample_id, self.hit_filter,
                        self.min_mapping_score,
                        self.min_alignment_score,
                        self.max_mismatches,
                        self.max_relative_mismatches,
                        self.min_continiously_matching,
                        self.filter_complexity)

        sam_file = None
        if self.sam:
            sam_file = open(self.sam, 'w', buffering=100 * 1024 * 1024)

        bam_process = None
        if self.bam:
            with open(self.bam, 'wb', buffering=100 * 1024 * 1024) as bam_file:

                command_line = samtools_path + ['view', '-b', '-1', '-S', '-@', '4', '/dev/stdin']
                command_line = ' '.join(command_line)

                if self.debug:
                    print 'Running: %s' % command_line

                bam_process = subprocess.Popen(command_line, shell=True, stdout=bam_file, stdin=PIPE)

        parser = None
        if taxonomy or quality or hits:
            parser = SAMParser()

        return parser, taxonomy, quality, hits, sam_file, bam_process

    def post_process(self, bam_process, sam_file, hits, taxonomy, quality, temp_path):

        if bam_process:
            bam_process.stdin.close()

        if sam_file:
            sam_file.close()

        if hits:
            hits.write()

        if taxonomy:
            print ''.join(taxonomy.get_summary(10))
            taxonomy.write()

        if quality:
            quality.write()

        shutil.rmtree(temp_path)

    def main(self, *args):

        debug = self.setup_logging()
        mapper_path, samtools_path, temp_path = self.validate_paths(self.mapper_path)
        mapper_process = self.run_mapper_process(mapper_path, temp_path)

        parser, taxonomy, quality, hits, sam_file, bam_process = self.setup_outputs(samtools_path)

        last_time = time.time()
        i = 0
        for line in iter(mapper_process.stdout.readline, ''):

            i += 1

            if debug and i > 0 and (i % 100000) == 0:
                this_time = time.time()
                print 'Mapping at %i reads per hour' % ((100000 * 3600) / (this_time - last_time))
                last_time = this_time

            if sam_file:
                sam_file.write(line)

            if bam_process:
                bam_process.stdin.write(line)

            if line[0] == '@':
                continue

            if parser:
                parsed_line = parser.parse(line)

            if taxonomy:
                taxonomy.count(parsed_line)

                # if i > 0 and (i % 100000) == 0:
                #     print ''.join(taxonomy.get_summary(10))

            if quality:
                quality.count(parsed_line)

            if hits:
                hits.count(parsed_line)

        self.post_process(bam_process, sam_file, hits, taxonomy, quality, temp_path)


@CLI.subcommand("rnamap")
class RNAmap(Mapper):
    """ Map input reads against a STAR index """

    index_dir = cli.SwitchAttr(['-i', '--index_dir'], str, mandatory=True,
                               help="Sets the index output directory")

    threads = cli.SwitchAttr(
        ['-t', '--threads'], cli.Range(1, 512), mandatory=False,
        help="Sets the number of threads to use",
        default=1)

    taxonomy = cli.SwitchAttr(
        ['-x', '--taxonomy'], str, mandatory=False,
        help="Output path for the taxonomy file; setting this option will also enable regular taxonomy output to stdout during mapping",
        default='')

    mapper_path = cli.SwitchAttr(['--star_path'], str, mandatory=False,
                          help="Path to STAR executable",
                          default='')

    samtools_path = cli.SwitchAttr(['--samtools_path'], str, mandatory=False,
                          help="Path to samtools executable",
                          default='')

    temp_path = cli.SwitchAttr(['--temporary_path'], str, mandatory=False,
                          help="Path to temporary directory in which to generate temp files. All temp files with be automatically deleted after execution is complete.",
                          default='')

    min_mapping_score = cli.SwitchAttr(['--min_mapping_score'], cli.Range(1, 255), mandatory=False,
                          help="Mimimum mapping score for saved hits (only applied to -v/--virana_hits)",
                          default=None)

    min_alignment_score = cli.SwitchAttr(['--min_alignment_score'], cli.Range(1, 255), mandatory=False,
                          help="Mimimum alignment score for saved hits (only applied to -v/--virana_hits)",
                          default=None)

    max_mismatches = cli.SwitchAttr(['--max_mismatches'],  cli.Range(0, 10000000), mandatory=False,
                          help="Maximum number of mismatches for saved hits (only applied to -v/--virana_hits)",
                          default=None)

    max_relative_mismatches = cli.SwitchAttr(['--max_relative_mismatches'], float, mandatory=False,
                          help="Maximum number of mismatches relative to read length for saved hits (only applied to -v/--virana_hits)",
                          default=None)

    min_continiously_matching = cli.SwitchAttr(['--min_continiously_matching'], cli.Range(0, 10000000), mandatory=False,
                          help="Minimum number of continious matches for saved hits (only applied to -v/--virana_hits)",
                          default=None)

    filter_complexity = cli.Flag(['--filter_complexity'],
                          help="Discard low-complexity reads (only applied to -v/--virana_hits). Adds some extra processing load to the mapping and may discard important information. Applies to all output files, including quality files (!)",
                          default=False)

    bam = cli.SwitchAttr(['-b', '--bam'], str, mandatory=False,
                         help="Path to unsorted, unindexed output BAM file",
                         default='')

    sam = cli.SwitchAttr(['-s', '--sam'], str, mandatory=False,
                         help="Path to output SAM file",
                         default='')

    qual = cli.SwitchAttr(['-q', '--qual'], str, mandatory=False,
                         help="Path to output quality file",
                         default='')

    hits = cli.SwitchAttr(['-v', '--virana_hits'], str, mandatory=False,
                          help="Path to bzip2-compressed tab-delimited output hit file",
                          default='')

    sample_id = cli.SwitchAttr(['--sample_id'], str, mandatory=False,
                          help="Alphanumeric string ([0-9a-zA-Z_-]*) used to designate sample information within the hit file",
                          default='no_sample_id')

    unmapped1 = cli.SwitchAttr(['--unmapped_end_1'], str, mandatory=False,
                         help="Output path to uncompressed fastq file containing unmapped reads, first ends only for paired ends.",
                         default='')

    unmapped2 = cli.SwitchAttr(['--unmapped_end_2'], str, mandatory=False,
                         help="Output path to uncompressed fastq file containing unmapped reads, second ends only for paired ends.",
                         default='')

    splice_junctions = cli.SwitchAttr(['--splice_junctions'], str, mandatory=False,
                         help="Input path to splice junction file (currently not implemented)",
                         default='')

    chimeric_mappings = cli.SwitchAttr(['--chimeric_mappings'], str, mandatory=False,
                         help="Ouput path to SAM file containing chimeric mappings",
                         default='')

    hit_filter = cli.SwitchAttr(
        ['-f', '--virana_hit_filter'], str, list=True, mandatory=False,
        help="Only generate hit groups that include at last one read mapping to a reference of this reference group.",
        default=[])

    debug = cli.Flag(["-d", "--debug"], help="Enable debug information")

    reads = cli.SwitchAttr(
        ['-r', '--reads'], str, list=True, mandatory=True,
        help="Sets the input reads. Add this parameter twice for paired end reads.")

    zipped = cli.Flag(["-z", "--zipped"], help="Input reads are zipped")

    sensitive = cli.Flag(
        ["--sensitive"], help="If given, mapping will process slower and more sensitive")


    def get_command_line(self, mapper_path, temp_path):

        first_ends      = []
        second_ends     = []
        single_ends     = []

        if len(self.reads) == 2:
            first, second = self.reads
            first_ends.append(first)
            second_ends.append(second)

        elif len(self.reads) == 1:
            single_ends.append(self.reads[0])

        if single_ends and not first_ends and not second_ends:
            reads = [','.join(single_ends)]

        elif first_ends and second_ends:
            reads = [','.join(first_ends), ','.join(second_ends)]

        command_line = [mapper_path] + ['--runMode', 'alignReads',
                        '--genomeDir', self.index_dir,
                        '--runThreadN', str(self.threads),
                        '--readMatesLengthsIn', 'NotEqual',
                        '--outFileNamePrefix', os.path.join(
                            temp_path, 'out'),
                        '--outSAMmode', 'Full',
                        '--outSAMstrandField', 'None',
                        '--outSAMattributes', 'Standard',
                        '--outSAMunmapped', 'Within',
                        '--outStd', 'SAM',
                        '--outFilterMultimapNmax', '1000',
                        '--outSAMprimaryFlag', 'AllBestScore',
                        '--outSAMorder', 'PairedKeepInputOrder']

        if self.unmapped1 or self.unmapped2:
            command_line += ['--outReadsUnmapped', 'Fastx']
        else:
            command_line += ['--outReadsUnmapped', 'None']


        if self.zipped:
            command_line += ['--readFilesCommand', 'zcat']

        if self.sensitive:
            command_line += ['--outFilterMultimapScoreRange', '10',
                          '--outFilterMismatchNmax', '60',
                          '--outFilterMismatchNoverLmax', '0.3',
                          '--outFilterScoreMin', '0',
                          '--outFilterScoreMinOverLread', '0.3',
                          '--outFilterMatchNmin', '0',
                          '--outFilterMatchNminOverLread', '0.66',
                          '--seedSearchStartLmax', '12',
                          '--winAnchorMultimapNmax', '50']

        command_line += ['--readFilesIn'] + reads

        return command_line

    def post_process(self, bam_process, sam_file, hits, taxonomy, quality, temp_path):

        try:
            if self.unmapped1:
                shutil.move(os.path.join(temp_path, 'out' + 'Unmapped.out.mate1'),\
                self.unmapped1)
        except IOError:
            pass

        try:
            if self.unmapped2:
                shutil.move(os.path.join(temp_path, 'out' + 'Unmapped.out.mate2'),\
                 self.unmapped2)
        except IOError:
            pass

        try:
            if self.chimeric_mappings:
                shutil.move(os.path.join(temp_path, 'out' + 'Chimeric.out.sam'),\
                 self.chimeric_mappings)
        except IOError:
            pass

        super(RNAmap, self).post_process(bam_process, sam_file, hits, taxonomy, quality, temp_path)



@CLI.subcommand("dnamap")
class DNAmap(Mapper):
    """ Map input reads with BWA-MEM against a BWA index """

    index_dir = cli.SwitchAttr(['-i', '--index_dir'], str, mandatory=True,
                               help="Sets the index output directory")

    threads = cli.SwitchAttr(
        ['-t', '--threads'], cli.Range(1, 512), mandatory=False,
        help="Sets the number of threads to use",
        default=1)

    taxonomy = cli.SwitchAttr(
        ['-x', '--taxonomy'], str, mandatory=False,
        help="Output path for the taxonomy file; setting this option will also enable regular taxonomy output to stdout during mapping",
        default='')

    samtools_path = cli.SwitchAttr(['--samtools_path'], str, mandatory=False,
                          help="Path to samtools executable",
                          default='')

    temp_path = cli.SwitchAttr(['--temporary_path'], str, mandatory=False,
                          help="Path to temporary directory in which to generate temp files. All temp files with be automatically deleted after execution is complete.",
                          default='')

    min_mapping_score = cli.SwitchAttr(['--min_mapping_score'], cli.Range(1, 255), mandatory=False,
                          help="Mimimum mapping score for saved hits (only applied to -v/--virana_hits)",
                          default=None)

    min_alignment_score = cli.SwitchAttr(['--min_alignment_score'], cli.Range(1, 255), mandatory=False,
                          help="Mimimum alignment score for saved hits (only applied to -v/--virana_hits)",
                          default=None)

    max_mismatches = cli.SwitchAttr(['--max_mismatches'],  cli.Range(0, 10000000), mandatory=False,
                          help="Maximum number of mismatches for saved hits (only applied to -v/--virana_hits)",
                          default=None)

    max_relative_mismatches = cli.SwitchAttr(['--max_relative_mismatches'], float, mandatory=False,
                          help="Maximum number of mismatches relative to read length for saved hits (only applied to -v/--virana_hits)",
                          default=None)

    min_continiously_matching = cli.SwitchAttr(['--min_continiously_matching'], cli.Range(0, 10000000), mandatory=False,
                          help="Minimum number of continious matches for saved hits (only applied to -v/--virana_hits)",
                          default=None)

    filter_complexity = cli.Flag(['--filter_complexity'],
                          help="Discard low-complexity reads (only applied to -v/--virana_hits). Adds some extra processing load to the mapping and may discard important information. Applies to all output files, including quality files (!)",
                          default=False)

    sample_id = cli.SwitchAttr(['--sample_id'], str, mandatory=False,
                          help="Alphanumeric string ([0-9a-zA-Z_-]*) used to designate sample information within the hit file",
                          default='no_sample_id')

    bam = cli.SwitchAttr(['-b', '--bam'], str, mandatory=False,
                         help="Path to unsorted, unindexed output BAM file",
                         default='')

    sam = cli.SwitchAttr(['-s', '--sam'], str, mandatory=False,
                         help="Path to output SAM file",
                         default='')

    qual = cli.SwitchAttr(['-q', '--qual'], str, mandatory=False,
                         help="Path to output quality file",
                         default='')

    hits = cli.SwitchAttr(['-v', '--virana_hits'], str, mandatory=False,
                          help="Path to bzip2-compressed tab-delimited output virana hit file",
                          default='')

    interleaved = cli.Flag(['--interleaved'],
                          help="Inputs FASTQ is an interleaved paired end file. ",
                          default=False)

    hit_filter = cli.SwitchAttr(
        ['-f', '--virana_hit_filter'], str, list=True, mandatory=False,
        help="Only generate hit groups that include at last one read mapping to a reference of this reference group.",
        default=[])

    debug = cli.Flag(["-d", "--debug"], help="Enable debug information")

    zipped = cli.Flag(["-z", "--zipped"], help="Input reads are zipped")

    sensitive = cli.Flag(
        ["--sensitive"], help="If given, mapping will process slower and more sensitive")

    mapper_path = cli.SwitchAttr(['--bwa_path'], str, mandatory=False,
                          help="Path to BWA executable",
                          default='')

    reads = cli.SwitchAttr(
        ['-r', '--reads'], str, list=True, mandatory=True,
        help="Sets the input reads. Add this parameter twice for paired end reads.")

    def get_command_line(self, mapper_path, temp_path):

        command_line = [mapper_path] + ['mem', '-t', str(self.threads), '-M', os.path.join(self.index_dir, 'index')]
        if self.interleaved:
            command_line += ['-p']
        command_line += self.reads

        return command_line


@CLI.subcommand("varmap")
class VARmap(Mapper):
    """ Map input reads with SMALT against a SMALT index """

    index_dir = cli.SwitchAttr(['-i', '--index_dir'], str, mandatory=True,
                               help="Sets the index output directory")

    threads = cli.SwitchAttr(
        ['-t', '--threads'], cli.Range(1, 512), mandatory=False,
        help="Sets the number of threads to use",
        default=1)

    taxonomy = cli.SwitchAttr(
        ['-x', '--taxonomy'], str, mandatory=False,
        help="Output path for the taxonomy file; setting this option will also enable regular taxonomy output to stdout during mapping",
        default='')

    samtools_path = cli.SwitchAttr(['--samtools_path'], str, mandatory=False,
                          help="Path to samtools executable",
                          default='')

    temp_path = cli.SwitchAttr(['--temporary_path'], str, mandatory=False,
                          help="Path to temporary directory in which to generate temp files. All temp files with be automatically deleted after execution is complete.",
                          default='')

    min_mapping_score = cli.SwitchAttr(['--min_mapping_score'], cli.Range(1, 255), mandatory=False,
                          help="Mimimum mapping score for saved hits (only applied to -v/--virana_hits)",
                          default=None)

    min_alignment_score = cli.SwitchAttr(['--min_alignment_score'], cli.Range(1, 255), mandatory=False,
                          help="Mimimum alignment score for saved hits (only applied to -v/--virana_hits)",
                          default=None)

    max_mismatches = cli.SwitchAttr(['--max_mismatches'],  cli.Range(0, 10000000), mandatory=False,
                          help="Maximum number of mismatches for saved hits (only applied to -v/--virana_hits)",
                          default=None)

    max_relative_mismatches = cli.SwitchAttr(['--max_relative_mismatches'], float, mandatory=False,
                          help="Maximum number of mismatches relative to read length for saved hits (only applied to -v/--virana_hits)",
                          default=None)

    min_continiously_matching = cli.SwitchAttr(['--min_continiously_matching'], cli.Range(0, 10000000), mandatory=False,
                          help="Minimum number of continious matches for saved hits (only applied to -v/--virana_hits)",
                          default=None)

    filter_complexity = cli.Flag(['--filter_complexity'],
                          help="Discard low-complexity reads (only applied to -v/--virana_hits). Adds some extra processing load to the mapping and may discard important information. Applies to all output files, including quality files (!)",
                          default=False)

    sample_id = cli.SwitchAttr(['--sample_id'], str, mandatory=False,
                          help="Alphanumeric string ([0-9a-zA-Z_-]*) used to designate sample information within the hit file",
                          default='no_sample_id')

    bam_input = cli.Flag(['--bam_input'],
                          help="Input is a bam file",
                          default=False)


    bam = cli.SwitchAttr(['-b', '--bam'], str, mandatory=False,
                         help="Path to unsorted, unindexed output BAM file",
                         default='')

    sam = cli.SwitchAttr(['-s', '--sam'], str, mandatory=False,
                         help="Path to output SAM file",
                         default='')

    qual = cli.SwitchAttr(['-q', '--qual'], str, mandatory=False,
                         help="Path to output quality file",
                         default='')

    hits = cli.SwitchAttr(['-v', '--virana_hits'], str, mandatory=False,
                          help="Path to bzip2-compressed tab-delimited output virana hit file",
                          default='')

    hit_filter = cli.SwitchAttr(
        ['-f', '--virana_hit_filter'], str, list=True, mandatory=False,
        help="Only generate hit groups that include at last one read mapping to a reference of this reference group.",
        default=[])

    debug = cli.Flag(["-d", "--debug"], help="Enable debug information")

    zipped = cli.Flag(["-z", "--zipped"], help="Input reads are zipped")

    sensitive = cli.Flag(
        ["--sensitive"], help="If given, mapping will process slower and more sensitive")

    mapper_path = cli.SwitchAttr(['--smalt_path'], str, mandatory=False,
                          help="Path to SMALT executable",
                          default='')

    debug = cli.Flag(["-d", "--debug"], help="Enable debug information")

    reads = cli.SwitchAttr(
        ['-r', '--reads'], str, list=True, mandatory=True,
        help="Sets the input reads. Add this parameter twice for paired end reads.")

    def get_command_line(self, mapper_path, temp_path):

        command_line = [mapper_path] + ['map', '-f', 'sam', '-l', '-n', self.threads, '-O', '-p']
        if self.sensitive:
            command_line += ['-x']
        command_line += self.reads

        return command_line


if __name__ == "__main__":
    CLI.run()
