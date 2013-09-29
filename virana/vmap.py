#!/usr/bin/env python
#from __future__ import print_function

import cProfile

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

from collections import defaultdict, Counter

from subprocess import PIPE

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


def profile_this(fn):
    def profiled_fn(*args, **kwargs):
        fpath = fn.__name__ + ".profile"
        prof = cProfile.Profile()
        ret = prof.runcall(fn, *args, **kwargs)
        prof.dump_stats(fpath)
        return ret
    return profiled_fn


class CLI(cli.Application):
    """RNA-Seq and DNA-Seq short read analysis by mapping to known reference sequences"""
    PROGNAME = "vmap"
    VERSION = "1.0.0"
    DESCRIPTION = """Virana vmap is an interface to the NCBI and ensembl reference databases that can
                     generate reference indexes for the short read mappers STAR (RNA-Seq) and
                     BWA-MEM (DNA-Seq). Short reads can be mapped to arbitrary combinations of
                     reference databases and the results can be summarized by taxonomic family
                     as well as stored as SAM file, unsorted BAM file, or as a HIT file that
                     models multimapping reads between specific reference databases."""
    USAGE = """The program has four modes that can be accessed by `vmap rnaindex`, `vmap dnaindex`, `vmap rnamap`, and `vmap dnamap.`"""

    def main(self, *args):

        if args:
            print("Unknown command %r" % (args[0]))
            print self.USAGE
            return 1

        if not self.nested_command:
            print("No command given")
            print self.USAGE
            return 1


@CLI.subcommand("rnaindex")
class RNAIndex(cli.Application):
    """ Creates a STAR index from a FASTA genome reference """

    reference_files = cli.SwitchAttr(
        ['-r', '--reference_file'], str, list=True, mandatory=True,
        help="Sets the reference genome(s) FASTA file." +
        " Multiple occurrences of this parameter are allowed.")
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
                          default='')
    sparse = cli.Flag(
        ["-s", "--sparse"], help="If given, a sparse index that requires less " +
        " RAM in the mapping phase will be constructed")

    debug = cli.Flag(["-d", "--debug"], help="Enable debug output")

    def main(self):

        if self.debug:
            logging.getLogger().setLevel(logging.DEBUG)

        # Obtain star executable
        star = [self.path and self.path or 'STAR']

        # Check if genome directory is existing
        for reference_file in self.reference_files:
            if not os.path.exists(reference_file):
                sys.stdout.write(
                    'Reference file %s nor existing, exiting' % reference_file)
                sys.exit(1)

        # Check if output directory is existing
        if not os.path.exists(self.index_dir):
            logging.debug(
                'Making output directory for index at %s' % self.index_dir)
            os.makedirs(self.index_dir)

        # # Make named pipe to extract genomes
        # pipe_path = os.path.abspath(os.path.join(self.genome_dir, 'pipe.fa'))
        # if os.path.exists(pipe_path):
        #     os.unlink(pipe_path)
        # os.mkfifo(pipe_path)

        # Make star command line
        cline = star + ['--runMode', 'genomeGenerate',
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

        # Run STAR reference generation process
        star_process = subprocess.Popen(' '.join(cline), shell=True, stdout=PIPE, stderr=PIPE)

        # Block until streams are closed by the process
        stdout, stderr = star_process.communicate()

        if stderr:
            sys.stderr.write(stderr)

        if self.debug and stdout:
            print stdout


@CLI.subcommand("dnaindex")
class DNAIndex(cli.Application):
    """ Creates a BWA index from a FASTA reference file """

    reference_file  = cli.SwitchAttr(['-r', '--reference_file'], str, mandatory=True,
                                 help="Sets the input reference FASTA file.")
    index_dir   = cli.SwitchAttr(['-i', '--index_dir'], str, mandatory=True,
                                 help="Sets the index output directory." +
                                 " Directory will be generated if not existing." +
                                 " Directory will be filled with several index files.")
    path        = cli.SwitchAttr(['-p', '--path'], str, mandatory=False,
                                 help="Path to BWA executable",
                                 default='')
    debug       = cli.Flag(["-d", "--debug"], help="Enable debug output")

    if debug:
        logging.getLogger().setLevel(logging.DEBUG)

    def main(self):

        # Obtain star executable
        bwa = [self.path and self.path or 'bwa']

        # Check if genome directory is existing
        if not os.path.exists(self.reference_file):
            sys.stdout.write('Genome file %s nor existing, exiting' % self.reference_file)
            sys.exit(1)

        # Check if output directory is existing
        if not os.path.exists(self.index_dir):
            logging.debug('Making output directory %s' % self.index_dir)
            os.makedirs(self.index_dir)

        # Make star command line
        cline = bwa + ['index', '-a', 'bwtsw', '-p', os.path.join(self.index_dir, 'index'), self.reference_file]

        if self.debug:
            print ' '.join(cline)

        # Run BWA index generation process
        bwa_process = subprocess.Popen(' '.join(cline), shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = bwa_process.communicate()

        if stderr:
            sys.stderr.write(stderr)

        if self.debug and stdout:
            print stdout





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
        qual        = alignment.read.qual
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

        phred_quality       = [q - 33 for q in qual]
        avg_phred_quality   = sum(phred_quality) / float(len(phred_quality))
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

@CLI.subcommand("rnamap")
class RNAmap(cli.Application):
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

    star_path = cli.SwitchAttr(['--star_path'], str, mandatory=False,
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

    def main(self):

        if self.debug:
            logging.getLogger().setLevel(logging.DEBUG)

        # Obtain star executable
        star        = [self.star_path and self.star_path or 'STAR']
        samtools    = [self.samtools_path and self.samtools_path or 'samtools']

        # Check if genome directory is existing
        if not os.path.exists(self.index_dir):
            sys.stdout.write('Index directory %s not existing, exiting' % self.genome_dir)
            sys.exit(1)

        if self.temp_path:
            temp_path = tempfile.mkdtemp(dir=self.temp_path)
        else:
            temp_path = tempfile.mkdtemp()

        first_ends      = []
        second_ends     = []
        single_ends     = []

        if len(self.reads) == 2:
            first, second = self.reads
            first_ends.append(first)
            second_ends.append(second)

        elif len(self.reads) == 1:
            single_ends.append(self.reads[0])

        else:
            sys.stdout.write('Invalid number of fastq files; provide either one (single end) or two (paired end)')
            sys.exit(1)

        if single_ends and not first_ends and not second_ends:
            reads = [','.join(single_ends)]

        elif first_ends and second_ends:
            reads = [','.join(first_ends), ','.join(second_ends)]


        star_cline = star + ['--runMode', 'alignReads',
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
            star_cline += ['--outReadsUnmapped', 'Fastx']
        else:
            star_cline += ['--outReadsUnmapped', 'None']


        if self.zipped:
            star_cline += ['--readFilesCommand', 'zcat']

        if self.sensitive:
            star_cline += ['--outFilterMultimapScoreRange', '10',
                          '--outFilterMismatchNmax', '60',
                          '--outFilterMismatchNoverLmax', '0.3',
                          '--outFilterScoreMin', '0',
                          '--outFilterScoreMinOverLread', '0.3',
                          '--outFilterMatchNmin', '0',
                          '--outFilterMatchNminOverLread', '0.66',
                          '--seedSearchStartLmax', '12',
                          '--winAnchorMultimapNmax', '50']

        star_cline += ['--readFilesIn'] + reads

        if self.debug:
            print ' '.join(star_cline)

        # Try if we can make the relevant files
        touch_files = [self.unmapped1, self.unmapped2, self.taxonomy, self.qual, self.hits, self.sam, self.bam]
        for file_path in touch_files:
            if file_path is None or file_path == '':
                continue
            try:
                with file(file_path, 'a'):
                    os.utime(file_path, None)
            except IOError:
                sys.stderr.write('Could not write output file %s\n' % file_path)
                sys.exit(1)

        star_process = subprocess.Popen(' '.join(
                                        star_cline), shell=True, stdout=PIPE)

        parser      = SAMParser()

        if self.taxonomy:
            taxonomy    = SAMTaxonomy(self.taxonomy)

        if self.qual:
            quality     = SAMQuality(self.qual)

        if self.hits:
            hits        = SAMHits(self.hits, self.sample_id, self.hit_filter,
                        self.min_mapping_score,
                        self.min_alignment_score,
                        self.max_mismatches,
                        self.max_relative_mismatches,
                        self.min_continiously_matching,
                        self.filter_complexity)

        if self.sam:
            sam_file = open(self.sam, 'w')

        if self.bam:
            with open(self.bam, 'wb', buffering=100 * 1024 * 1024) as bam_file:
                samtools_cline = samtools + [
                    'view', '-b', '-1', '-S', '-@', '4', '/dev/stdin']
                if self.debug:
                    print ' '.join(samtools_cline)
                samtools_process = subprocess.Popen(' '.join(samtools_cline), shell=True, stdout=bam_file, stdin=PIPE)


        do_sam = self.sam
        do_bam = self.bam
        do_taxonomy = self.taxonomy
        do_qual = self.qual
        do_hits = self.hits
        do_parse = do_taxonomy or do_qual or do_hits

        for i, line in enumerate(iter(star_process.stdout.readline, '')):

            if do_sam:
                sam_file.write(line)

            if do_bam:
                samtools_process.stdin.write(line)

            if line[0] == '@':
                continue

            if do_parse:
                parsed_line = parser.parse(line)

            if  do_taxonomy:
                taxonomy.count(parsed_line)
                if i > 0 and (i % 50000) == 0:
                    print ''.join(taxonomy.get_summary(10))

            if do_qual:
                quality.count(parsed_line)

            if do_hits:
                hits.count(parsed_line)

        if do_bam:
            samtools_process.stdin.close()

        if do_sam:
            sam_file.close()

        if do_hits:
            hits.write()

        if do_taxonomy:
            print ''.join(taxonomy.get_summary(10))
            taxonomy.write()

        if do_qual:
            quality.write()

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

        shutil.rmtree(temp_path)

@CLI.subcommand("dnamap")
class DNAmap(cli.Application):
    """ Map input reads against a BWA index """

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

    hit_filter = cli.SwitchAttr(
        ['-f', '--virana_hit_filter'], str, list=True, mandatory=False,
        help="Only generate hit groups that include at last one read mapping to a reference of this reference group.",
        default=[])

    debug = cli.Flag(["-d", "--debug"], help="Enable debug information")

    zipped = cli.Flag(["-z", "--zipped"], help="Input reads are zipped")

    sensitive = cli.Flag(
        ["--sensitive"], help="If given, mapping will process slower and more sensitive")


    bwa_path = cli.SwitchAttr(['--bwa_path'], str, mandatory=False,
                          help="Path to BWA executable",
                          default='')

    debug = cli.Flag(["-d", "--debug"], help="Enable debug information")

    if debug:
        logging.getLogger().setLevel(logging.DEBUG)

    reads = cli.SwitchAttr(
        ['-r', '--reads'], str, list=True, mandatory=True,
        help="Sets the input reads. Add this parameter twice for paired end reads.")

    def main(self):

        if self.debug:
            logging.getLogger().setLevel(logging.DEBUG)

        # Obtain star executable
        bwa         = [self.bwa_path and self.bwa_path or 'bwa']
        samtools    = [self.samtools_path and self.samtools_path or 'samtools']

        # Check if genome directory is existing
        if not os.path.exists(self.index_dir):
            sys.stdout.write('Index directory %s not existing, exiting'\
                % self.genome_dir)
            sys.exit(1)

        if len(self.reads) not in (1, 2):
            message = 'Invalid number of FASTQ files; supply either one (single end) or two (paired end)\n'
            sys.stderr.write(message)
            sys.exit(1)

        bwa_cline = bwa + ['mem', '-t', str(self.threads), '-M', os.path.join(self.index_dir, 'index')]

        bwa_cline += self.reads

        if self.debug:
            print ' '.join(bwa_cline)

        bwa_process = subprocess.Popen(' '.join(bwa_cline), shell=True, stdout=PIPE)

        parser      = SAMParser()

        if self.taxonomy:
            taxonomy    = SAMTaxonomy(self.taxonomy)

        if self.qual:
            quality     = SAMQuality(self.qual)

        if self.hits:
            hits        = SAMHits(self.hits, self.sample_id, self.hit_filter,
                        self.min_mapping_score,
                        self.min_alignment_score,
                        self.max_mismatches,
                        self.max_relative_mismatches,
                        self.min_continiously_matching,
                        self.filter_complexity)

        if self.sam:
            sam_file = open(self.sam, 'w', buffering=100 * 1024 * 1024)

        if self.bam:
            with open(self.bam, 'wb', buffering=100 * 1024 * 1024) as bam_file:
                samtools_cline = samtools + [
                    'view', '-b', '-1', '-S', '-@', '4', '/dev/stdin']
                if self.debug:
                    print ' '.join(samtools_cline)
                samtools_process = subprocess.Popen(' '.join(samtools_cline), shell=True, stdout=bam_file, stdin=PIPE)

        do_sam = self.sam
        do_bam = self.bam
        do_taxonomy = self.taxonomy
        do_qual = self.qual
        do_hits = self.hits
        do_parse = do_taxonomy or do_qual or do_hits

        for i, line in enumerate(iter(bwa_process.stdout.readline, '')):

            if do_sam:
                sam_file.write(line)

            if do_bam:
                samtools_process.stdin.write(line)

            if do_parse:
                parsed_line = parser.parse(line)

            if do_taxonomy:
                taxonomy.count(parsed_line)

                if i > 0 and (i % 10000) == 0:
                    print ''.join(taxonomy.get_summary(10))

            if do_qual:
                quality.count(parsed_line)

            if do_hits:
                hits.count(parsed_line)

        if do_bam:
            samtools_process.stdin.close()

        if do_sam:
            sam_file.close()

        if do_hits:
            hits.write()

        if do_taxonomy:
            print ''.join(taxonomy.get_summary(10))
            taxonomy.write()

        if do_qual:
            quality.write()


if __name__ == "__main__":
    CLI.run()
