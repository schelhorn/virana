#!/usr/bin/env python

import sys
import os
import logging
import tempfile
import shutil
import glob

import subprocess
import bz2
import re
import glob

from shutil import rmtree
from subprocess import PIPE
from collections import defaultdict, Counter

try:
    import pysam
except ImportError:
    message = 'This script requires the pysam python package\n'
    sys.stderr.write(message)
    sys.exit(1)

try:
    from plumbum import cli
except ImportError:
    message = 'This script requires the plumbum python package\n'
    sys.stderr.write(message)
    sys.exit(1)

try:
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC, Gapped
except ImportError:
    message = 'This script requires the BioPython python package\n'
    sys.stderr.write(message)
    sys.exit(1)

try:
    import HTSeq
except ImportError:
    message = 'This script requires the HTSeq python package\n'
    sys.stderr.write(message)
    sys.exit(1)

logging.basicConfig(level=logging.INFO, format='%(message)s')


def which(program):

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

class CLI(cli.Application):
    """Identification and visualization of homologous regions based on multi-mapping reads."""
    PROGNAME = "vhom"
    VERSION = "1.0.0"
    DESCRIPTION = \
    """DESCRIPTION: virana vhom - identification of homologous regions.

    The virana homology utility ('vhom') analyzes the homology relationships
    within mapped reads in order to extract homologous regions, i.e.
    nucleotide stretches that display high sequence similarity to a pathogen
    and, optionally, also to human factors. These regions can be analyzed with
    regard to the homologous (i.e., transcriptomic and genomic) contexts
    of the read data. This greatly facilitates delineation of microbial from
    human sequence regions (i.e., syntenic sequence stretches consisting of
    assembled reads that align well to several references).

    https://github.com/schelhorn/virana

    Schelhorn S-E, Fischer M, Tolosi L, Altmueller J, Nuernberg P, et al. (2013)
    Sensitive Detection of Viral Transcripts in Human Tumor Transcriptomes.
    PLoS Comput Biol 9(10): e1003228. doi:10.1371/journal.pcbi.1003228"""

    USAGE = """USAGE: The program has one mode that can be accessed by
       [vhom | python vhom.py] regions
       """

    def main(self, *args):

        if args:
            print self.DESCRIPTION
            print
            print self.USAGE
            print("ERROR: Unknown command %r" % (args[0]))
            return 1

        if not self.nested_command:
            print self.DESCRIPTION
            print
            print self.USAGE
            print("ERROR : No command given")
            return 1

class SequenceProxy:

    def __init__(self, references_path, human_transcripts_path=None):

        self.references_path = references_path
        self.human_transcripts_path = human_transcripts_path

        self.human_transcript_records = None
        self.reference_records = None
        self.hit_records = {}

        self.transcript_ranges = None

    def add_hit_file(self, hit_file_path):

        logging.debug('SequenceProxy: adding new hit file %s' % hit_file_path)

        hit_handle = bz2.BZ2File(hit_file_path, 'rb')
        try:
            hit_handle.read(10)
        except IOError:
            hit_handle = open(hit_file_path, 'rb')
        except:
            sys.stderr.write('SequenceProxy: could not process hit file %s due to errors, skipping.\n' % hit_file_path)
            return

        added = 0
        for hit_record in SeqIO.parse(hit_handle, "fasta"):
            identifier  = hit_record.id
            description = hit_record.description.strip().split(' ')[1]
            if identifier in self.hit_records:
                self.hit_records[identifier].description += ('|' + description)
            else:
                self.hit_records[identifier] = hit_record
            added += 1

        logging.debug('SequenceProxy: hit file with %i entries added' % added)

    def get_read_record(self, read_id):

        try:
            record = self.hit_records[read_id]
            return SeqRecord(record.seq, record.id, '', '')

        except KeyError:
            logging.debug('SequenceProxy: unknown hit record %s, skipping' % read_id)
            return None

    def get_reference_record(self, identifier):

        if not self.reference_records:
            logging.debug('SequenceProxy: indexing reference records in order to provide viral and human genomic context')
            self.reference_records = SeqIO.index(self.references_path, 'fasta')

        try:
            record = self.reference_records[identifier]
            return SeqRecord(record.seq, record.id, '', '')
        except KeyError:
            return None

    def get_transcript_record(self, identifier):

        if not self.human_transcripts_path:
            return None

        if not self.human_transcript_records:
            logging.debug('SequenceProxy: indexing human transcripts')
            self.human_transcript_records = SeqIO.index(
                self.human_transcripts_path, 'fasta')

        try:
            record = self.human_transcript_records[identifier]
            return SeqRecord(record.seq, record.id, '', '')

        except KeyError:
            return None

    def get_overlapping_human_transcript_identifiers(self, chromosome, start, end):

        if not self.human_transcripts_path:
            return set()

        if not self.transcript_ranges:
            self._set_transcript_ranges()

        interval = HTSeq.GenomicInterval(chromosome, start, end, '.')
        all_transcript_ids = set()
        for interval, transcript_ids in self.transcript_ranges[interval].steps():
            all_transcript_ids = all_transcript_ids.union(transcript_ids)

        return all_transcript_ids

    def _set_transcript_ranges(self):

        logging.debug(
            'SequenceProxy: Preparing transcript ranges, please wait...')

        self.transcript_ranges = HTSeq.GenomicArrayOfSets(
            'auto', stranded=False)

        for record in SeqIO.parse(self.human_transcripts_path, 'fasta'):
            transcript_id = record.id
            for exon in record.description.strip().split(' ')[1].split('|'):
                chromosome, start, end = exon.split(';')
                start, end = sorted((int(start), int(end)))
                if start == end:
                    continue
                interval = HTSeq.GenomicInterval(chromosome, start, end, ".")
                self.transcript_ranges[interval] += transcript_id


class JalviewRunner:

    def __init__(self, jalview_jar_dir, tmp_dir=None):

        self.jalview_jar_dir = jalview_jar_dir

        if tmp_dir:
            self.tmp_dir = tempfile.mkdtemp(
                dir=os.path.abspath(os.path.expandvars(tmp_dir)))
        else:
            self.tmp_dir = tempfile.mkdtemp()

        self.config = ['#---JalviewX Properties File---',
                       '#Fri Nov 04 14:23:53 CET 2011',
                       'FIGURE_AUTOIDWIDTH=true',
                       'ID_ITALICS=false',
                       'SORT_ALIGNMENT=No sort',
                       'SHOW_IDENTITY=true',
                       'FONT_NAME=SansSerif',
                       'GAP_SYMBOL=.',
                       'SHOW_QUALITY=false',
                       'SHOW_GROUP_CONSERVATION=false',
                       'FONT_STYLE=plain',
                       'ANTI_ALIAS=false',
                       'SORT_BY_TREE=false',
                       'SHOW_CONSENSUS_HISTOGRAM=false',
                       'SHOW_OVERVIEW=false',
                       'DEFAULT_COLOUR=Nucleotide',
                       'SHOW_CONSENSUS_LOGO=false',
                       'SHOW_ANNOTATIONS=true',
                       'SHOW_UNCONSERVED=true',
                       'AUTO_CALC_CONSENSUS=true',
                       'PAD_GAPS=true',
                       'FONT_SIZE=10',
                       'RIGHT_ALIGN_IDS=false',
                       'WRAP_ALIGNMENT=true',
                       'FASTA_JVSUFFIX=false',
                       'PILEUP_JVSUFFIX=false',
                       'SHOW_JVSUFFIX=false']

    def _make_configuration(self, sub_tmp_dir):

        config_path = os.path.join(sub_tmp_dir, 'config.txt')
        logging.debug(
            'JalviewRunner: writing configuration file %s' % config_path)
        with open(config_path, 'w') as config_file:
            for c in self.config:
                config_file.write(c + '\n')

        return config_path

    def _delete_tmp_dir(self):

        try:
            rmtree(self.tmp_dir)
        except OSError:
            pass

    def _delete_sub_tmp_dir(self, sub_tmp_dir):

        try:
            rmtree(sub_tmp_dir)
        except OSError:
            pass

    def run(self, fasta_path, png_output_path):

        fasta_path = os.path.abspath(os.path.expandvars(fasta_path))

        png_output_path = os.path.abspath(os.path.expandvars(png_output_path))

        sub_tmp_dir = tempfile.mkdtemp(dir=self.tmp_dir)
        config_path = self._make_configuration(sub_tmp_dir)

        logging.debug(
            'JalviewRunner: running Jalview on fasta %s with temp dir %s' %
            (fasta_path, sub_tmp_dir))

        cline = ['java', '-Djava.ext.dirs=%s/lib' % self.jalview_jar_dir,
                 '-cp %s/jalview.jar' % self.jalview_jar_dir, 'jalview.bin.Jalview',
                 '-open', fasta_path, '-nodisplay', '-png', png_output_path,
                 '-colour', 'nucleotide', '-props', config_path]

        # Run sibelia process
        java_process = subprocess.Popen(
            ' '.join(cline), shell=True, stdout=PIPE, stderr=PIPE)

        # Block until streams are closed by the process
        stdout, stderr = java_process.communicate()

        if stderr:
            sys.stderr.write('JalviewRunner: ' + stderr.replace('\n', ' '))

        self._delete_sub_tmp_dir(sub_tmp_dir)

        return png_output_path

    def __del__(self):

        self._delete_tmp_dir()


class LastzRunner:

    def __init__(self, lastz_path=None, word_length=7):

        self.lastz_path = lastz_path
        self.word_length = word_length

    def align_to_bam_file(self, reference_fasta_path, query_fasta_path, output_bam_path, multiple=False, assert_record=None):

        logging.debug('LastzRunner: running on reference %s and query %s' %
                     (reference_fasta_path, query_fasta_path))
        output_sam_path = os.path.abspath(
            os.path.expandvars(output_bam_path.replace('.bam', '.sam')))
        output_bam_unsorted_path = os.path.abspath(
            os.path.expandvars(output_bam_path + '.unsorted'))

        logging.debug(
            'LastzRunner: aligning with output in temporary sam file %s' %
            output_sam_path)
        with open(output_sam_path, 'w') as output_sam_handler:
            for line in self._align(reference_fasta_path, query_fasta_path, multiple):
                output_sam_handler.write(line)

        logging.debug(
            'LastzRunner: transforming sam into unsorted bam file %s' %
            output_bam_unsorted_path)
        input_sam_handler = pysam.Samfile(output_sam_path, "r")
        output_bam_file = pysam.Samfile(
            output_bam_unsorted_path, "wb", template=input_sam_handler)

        logging.debug(
            'LastzRunner: copying from sam file to bam file')
        for s in input_sam_handler:
            output_bam_file.write(s)
        output_bam_file.close()

        logging.debug('LastzRunner: sorting and indexing bam file %s' %
                      output_bam_path)
        pysam.sort(output_bam_unsorted_path,
                   output_bam_path.replace('.bam', ''))

        pysam.index(output_bam_path)

    def align_to_samlines(self, reference_fasta_path, query_fasta_path, multiple=False):

        logging.debug('LastzRunner: running on reference %s and query %s' %
                      (reference_fasta_path, query_fasta_path))

        for line in self._align(reference_fasta_path, query_fasta_path, multiple):
            yield line

    def _align(self, reference_fasta_path, query_fasta_path, multiple):

        reference_fasta_path = os.path.abspath(
            os.path.expandvars(reference_fasta_path))
        query_fasta_path = os.path.abspath(
            os.path.expandvars(query_fasta_path))

        if multiple:
            cline = [self.lastz_path] + [
                reference_fasta_path +
                '[multiple]', query_fasta_path + '[nameparse=darkspace]',
                '--format=sam', '--strand=both', '--gapped', '--chain',
                '--nogfextend', '--seed=match%i' % self.word_length]

        else:
            cline = [self.lastz_path] + [
                reference_fasta_path, query_fasta_path +
                '[nameparse=darkspace]',
                '--format=sam', '--strand=both', '--gapped',
                '--chain', '--nogfextend', '--seed=match%i' % self.word_length]

        # Run lastz process
        logging.debug('Running lastz with command line ' + ' '.join(cline))
        lastz_process = subprocess.Popen(
            ' '.join(cline), shell=True, stdout=PIPE, stderr=PIPE)

        # Filter multiple read mappings (one to each strand) by only keeping
        # the one with the most matches
        alignments = []
        for i, line in enumerate(iter(lastz_process.stdout.readline, '')):
            if line.startswith('@'):
                yield line
            elif line.lower().startswith('read'):

                fields = line.split('\t')
                read_name = fields[0]
                read_cigar = fields[5]
                if not alignments:
                    alignments = [read_name, line, read_cigar]
                    continue
                else:
                    if alignments[0] == read_name:
                        this_mapping = sum(
                            [int(c[:-1]) for c in re.findall('[0-9]*M', read_cigar)])
                        other_mapping = sum(
                            [int(c[:-1]) for c in re.findall('[0-9]*M', alignments[2])])

                        if other_mapping > this_mapping:
                            yield alignments[1]
                        else:
                            yield line
                    else:
                        yield line
                    alignments = []
            else:
                yield line


class ConsensusRunner:

    def __init__(self, ambiguity_cutoff=0.3):

        self.ambiguity_cutoff = ambiguity_cutoff
        self.ambiguity = {
            'GT': 'K', 'AC': 'M', 'AG': 'R', 'CT': 'Y', 'CG': 'S',
            'AT': 'W', 'CGT': 'B', 'ACG': 'V', 'ACT': 'H', 'AGT': 'D',
            'ACGT': 'N'}

        self.alphabet = Gapped(IUPAC.unambiguous_dna, "-")

    def run(self, input_sam_path, output_fasta_path):

        samfile = pysam.Samfile(input_sam_path, "rb")
        references = samfile.references

        logging.debug(
            'ConsensusRunner: writing consensus of bam file %s' % input_sam_path)

        assert len(references) == 1

        consensus = defaultdict(list)
        conditions = {}

        base_mapper = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4, '-': 5}
        ambiguity = self.ambiguity

        alphabet = Gapped(IUPAC.ambiguous_dna)
        frequency_cutoff = self.ambiguity_cutoff

        # Get all conditions of aligned reads
        for alignedread in samfile.fetch(references[0]):

            name = alignedread.qname
            length = len(alignedread.seq)
            condition = None

            if name.lower().startswith('read'):
                condition = ';'.join(name.split(';')[:2])
            else:
                condition = ';'.join(name.split(';')[:-1])

            if condition in conditions:
                if name not in conditions[condition][0]:
                    conditions[condition][0].add(name)
                    conditions[condition][1] += 1
                    conditions[condition][2] += length
            else:
                conditions[condition] = [set([name]), 1, length]

        # Prepare count dictionaries
        consensus = defaultdict(list)

        for i, pileupcolumn in enumerate(samfile.pileup(references[0], max_depth=100000)):

            if i % 1000 == 0:
                logging.debug(
                    'ConsensusRunner: making consensus at position %i' % i)

            counters = {}
            for condition in conditions:
                counters[condition] = [0, 0, 0, 0, 0, 0]  # ACGTN-

            for pileupread in pileupcolumn.pileups:

                alignment = pileupread.alignment
                name = alignment.qname

                if pileupread.is_del:
                    base = '-'
                else:
                    base = alignment.seq[pileupread.qpos].upper()

                if name.lower().startswith('read'):
                    condition = ';'.join(name.split(';')[:2])
                    counters[condition][base_mapper[base]] += 1
                else:
                    condition = ';'.join(name.split(';')[:-1])
                    counters[condition][base_mapper[base]] += 1

            for condition, counts in counters.iteritems():
                depth = float(sum(counts))
                bases = []
                a, c, g, t, n, gap = counts

                if depth > 0:
                    if (n / depth) >= frequency_cutoff:
                        consensus[condition].append('N')
                    else:
                        if (a / depth) >= frequency_cutoff:
                            bases.append('A')
                        if (c / depth) >= frequency_cutoff:
                            bases.append('C')
                        if (g / depth) >= frequency_cutoff:
                            bases.append('G')
                        if (t / depth) >= frequency_cutoff:
                            bases.append('T')

                if len(bases) > 1:
                    consensus[condition].append(
                        ambiguity[''.join(sorted(bases))])
                elif len(bases) == 1:
                    consensus[condition].append(bases[0])
                else:
                    consensus[condition].append('-')

        # Split consensuses by type
        reference_consensuses = []
        read_consensuses = []
        for condition, sequence_list in consensus.iteritems():
            if condition.lower().startswith('read'):
                read_consensuses.append((''.join(sequence_list), condition))
            else:
                reference_consensuses.append(
                    (''.join(sequence_list), condition))

        reference_consensuses.sort(reverse=True)
        read_consensuses.sort(reverse=True)

        # Get start and end of reads (remove fully gapped positions)
        start = None
        end = None
        for read_consensus, condition in read_consensuses:

            # Forward direction
            for i, char in enumerate(read_consensus):
                if char != '-':
                    if start is None:
                        start = i
                    else:
                        start = min(start, i)
                    break

            # Reverse direction
            for i, char in enumerate(read_consensus[::-1]):
                if char != '-':
                    if end is None:
                        end = i
                    else:
                        end = min(end, i)
                    break

        if end == 0 or end is None:
            end = None
        else:
            end = -end

        # Filter all records by start and end of reads
        reference_consensuses =\
            [(r[start:end], c) for r, c in reference_consensuses]

        read_consensuses =\
            [(r[start:end], c) for r, c in read_consensuses]

        # Write consensus records to file
        with open(output_fasta_path, 'w',  buffering=100 * 1024 * 1024) as output_handler:
            for consensus, condition in read_consensuses + reference_consensuses:
                stats = ';'.join(map(str, conditions[condition][1:]))
                record = SeqRecord(
                    Seq(consensus, alphabet=alphabet), id=condition + ';' + stats, name='', description='')
                SeqIO.write([record], output_handler, "fasta")


class Group:

    def __init__(self, family_name, qualified_family_name, sequence_proxy, tmp_dir=None):

        self.family_name = family_name
        self.qualified_family_name = qualified_family_name
        self.sequence_proxy = sequence_proxy

        self.regions = []
        self.fasta_path = None

        if tmp_dir:
            self.tmp_dir = tempfile.mkdtemp(
                dir=os.path.abspath(os.path.expandvars(tmp_dir)))
        else:
            self.tmp_dir = tempfile.mkdtemp()

    def add_region(self, read_id, pathogen_mapping_locations, human_mapping_locations):

        region = Region(self.sequence_proxy, self.tmp_dir)

        region.add_read(read_id)

        for mapping_location in list(pathogen_mapping_locations) + list(human_mapping_locations):
            reference = ';'.join(mapping_location[:-2])
            region.add_reference(
                reference, int(mapping_location[-2]), int(mapping_location[-1]))

        self.regions.append(region)

    def write_outputs(self, output_dir, lastz_runner, consensus_builder, jalview_runner):

        if not self.regions:
            logging.debug('Group %s: no regions for family' % self.family_name)
            return

        made_output_dir = False
        for i, region in enumerate(self.regions):

            if not made_output_dir:
                output_dir = os.path.abspath(
                    os.path.expandvars(os.path.join(output_dir, self.qualified_family_name)))
                logging.debug('Group %s: writing output files to %s' %
                              (self.family_name, output_dir))
                try:
                    os.makedirs(output_dir)
                except OSError:
                    pass
                made_output_dir = True

            logging.debug('Group %s: writing region number %i files to %s' %
                          (self.family_name, i + 1, output_dir))

            unaligned_path = os.path.join(
                output_dir, 'region_%i_unaligned.fa.bzip2' % (i + 1))
            compressed_unaligned = bz2.BZ2File(unaligned_path, 'wb')
            with open(region.get_unaligned_fasta_path()) as unaligned_fasta:
                for line in unaligned_fasta:
                    compressed_unaligned.write(line)

            bam_path = os.path.join(
                output_dir, 'region_%i_alignment.bam' % (i + 1))
            temporary_alignment_bam_path = region.get_alignment_bam_path(
                lastz_runner, assert_record='Read')
            shutil.copy(temporary_alignment_bam_path, bam_path)

            shutil.copy(
                temporary_alignment_bam_path + '.bai', bam_path + '.bai')

            consensus_path = os.path.join(
                output_dir, 'region_%i_consensus.fa' % (i + 1))
            temporary_consensus_path = region.get_consensus_fasta_path(
                consensus_builder, temporary_alignment_bam_path)
            shutil.copy(temporary_consensus_path, consensus_path)

            if jalview_runner:

                jalview_path = os.path.join(
                    output_dir, 'region_%i_consensus.png' % (i + 1))
                temporary_jalview_path = region.get_consensus_figure_path(
                    jalview_runner, temporary_consensus_path)
                shutil.copy(temporary_jalview_path, jalview_path)

    def _delete_temporary_dir(self):

        self.fasta_path = None
        try:
            rmtree(self.tmp_dir)
        except OSError:
            pass

    def filter_regions(self, min_region_length=50, min_read_number=5):

        logging.debug('Group %s: filtering %i candidate regions' %
                      (self.family_name, len(self.regions)))

        filtered = [r for r in self.regions if len(
            r) >= min_read_number and r.get_longest_reference_length() >= min_region_length]

        ordered = sorted(filtered, reverse=True)

        if ordered:
            length = ordered[0].length
        else:
            length = 0

        logging.debug(
            'Group %s: %i candidate regions remained after filtering, the longest is %i bp long' %
            (self.family_name, len(ordered), length))

        self.regions = ordered

        return ordered

    def merge_regions(self, max_gap_length):
        """ Merges candidate regions into homologous regions. """

        logging.debug('Group %s: merging %i candidate regions' %
                      (self.family_name, len(self.regions)))

        if len(self.regions) > 1:

            potentially_mergable = self.regions
            not_mergable = []

            while len(potentially_mergable) > 1:

                merged = False
                current = potentially_mergable[0]
                compared_to = potentially_mergable[1:]

                for region in compared_to:
                    if region.overlaps(current, max_gap_length):
                        region.merge(current)
                        region.clean_references(max_gap_length)
                        #logging.debug('Group %s: merged a region. %i potentially mergable candidate regions remaining' % (self.family_name, len(potentially_mergable)))
                        potentially_mergable = compared_to
                        merged = True
                        break

                if not merged:
                    not_mergable.append(current)
                    potentially_mergable = compared_to
                    #logging.debug('Group %s: not merged a region. %i potentially mergable candidate regions remaining' % (self.family_name, len(potentially_mergable)))

            results = not_mergable + potentially_mergable

            logging.debug('Group %s: merged into %i regions' %
                          (self.family_name, len(results)))

            self.regions = results

        else:
            logging.debug(
                'Group %s: found only 1 region, no mergin necessary' % self.family_name)

    def add_transcripts_to_regions(self):

        logging.debug('Group %s: enriching %i regions with human transcript information'
                      % (self.family_name, len(self.regions)))

        for region in self.regions:

            added_transcripts = 0

            for identifier, locations in region.references.iteritems():
                chromosome = identifier.split(';')[-1]

                for start, end in locations:
                    transcript_identifiers = self.sequence_proxy.get_overlapping_human_transcript_identifiers(
                        chromosome, start, end)

                    for transcript_identifier in transcript_identifiers:
                        region.add_transcript(transcript_identifier)
                        added_transcripts += 1

            logging.debug('Group %s: added %i transcripts to a region' %
                          (self.family_name, added_transcripts))

    def __str__(self):

        s = "Group %s with %i regions" % (self.family_name, len(self.regions))

        return s


class GroupGenerator:

    """ Produces Groups from Hitfiles. Reads are assigned to groups based on
        taxonomic families. Transcripts overlapping human read mappings as well
        as genome sequences of pathogens the reads map to are added to the
        Groups.
    """

    def __init__(self, sequence_proxy, reference_database_filter=None, pathogen_family_filter=None, tmp_dir=None):

        self.sequence_proxy = sequence_proxy

        self.tmp_path = tmp_dir
        self.reference_database_filter = reference_database_filter
        self.pathogen_family_filter = pathogen_family_filter

        self.groups = {}

        if tmp_dir:
            self.tmp_dir = tempfile.mkdtemp(
                dir=os.path.abspath(os.path.expandvars(tmp_dir)))
        else:
            self.tmp_dir = tempfile.mkdtemp()

    def get_groups(self, min_read_number=5):

        logging.debug(
            'GroupGenerator: generating groups from %i hit entries, please wait...' %
            len(self.sequence_proxy.hit_records))


        for read_id, record in self.sequence_proxy.hit_records.iteritems():

            # Extract mapping locations to human DNA, human transcript, and
            # pathogen genomes
            mapping_locations = record.description.strip().split(
                ' ')[1].split('|')

            pathogen_families = defaultdict(set)
            human_mapping_locations = set()

            for mapping_location in mapping_locations:
                fields = mapping_location.split(';')

                # Convert start and end to integers
                fields[-2] = int(fields[-2])  # start
                fields[-1] = int(fields[-1])  # end

                # Add mapping locations to human cDNA or human DNA
                if fields[2] == 'Homo_sapiens':
                    human_mapping_locations.add(tuple(fields))

                # Add mapping locations to non-human references
                else:
                    reference_database, family = fields[:2]
                    if self.reference_database_filter and reference_database\
                            not in self.reference_database_filter:
                        continue
                    if self.pathogen_family_filter and family\
                            not in self.pathogen_family_filter:
                        continue
                    pathogen_families[
                        (reference_database, family)].add(tuple(fields))

            if not pathogen_families:
                continue

            # Process the hits to different pathogen families
            for (reference_database, family), pathogen_locations in pathogen_families.iteritems():

                qualified_family_name = reference_database + '_' + family

                # Obtain existing group or make new group
                if qualified_family_name in self.groups:
                    group = self.groups[qualified_family_name]
                else:
                    group = Group(family, qualified_family_name,
                                  self.sequence_proxy, self.tmp_dir)
                    self.groups[qualified_family_name] = group

                group.add_region(
                    read_id, pathogen_locations, human_mapping_locations)

        # Exclude groups that have collected to few reads
        for qualified_family_name, group in self.groups.items():
            if len(group.regions) < min_read_number:
                del self.groups[qualified_family_name]
            else:
                logging.debug(
                    'GroupGenerator: made new group for family %s' % qualified_family_name)

        return self.groups


class Region:

    def __init__(self, sequence_proxy, tmp_dir=None):

        self.sequence_proxy = sequence_proxy

        self.references = defaultdict(set)
        self.reads = set()
        self.transcripts = set()
        self.length = None

        self.sorted_reference_positions = None

        self.master_tmp = tmp_dir
        self.tmp_dir = None

        self.unaligned_fasta_path = None
        self.alignment_sam_path = None
        self.consensus_fasta_path = None

    def add_reference(self, name, start, end):

        assert start < end
        self.length = None
        self.longest_reference_id = None
        self.references[name].add((start, end))

    def add_read(self, name):

        self.reads.add(name)

    def add_transcript(self, name):

        self.transcripts.add(name)

    def get_tmp_path(self):

        if self.tmp_dir:
            return self.tmp_dir

        if self.master_tmp:
            self.tmp_dir = tempfile.mkdtemp(
                dir=os.path.abspath(os.path.expandvars(self.master_tmp)))
        else:
            self.tmp_dir = tempfile.mkdtemp()

        return self.tmp_dir

    def _distance(self, range1, range2):

        first, second = sorted((range1, range2))

        if first[0] <= first[1] < second[0]:
            return (second[0] - first[1])
        else:
            return 0

    def _delete_temporary_dir(self):

        self.unaligned_fasta_path = None
        if self.tmp_dir:
            try:
                rmtree(self.tmp_dir)
            except OSError:
                pass

    def _get_unaligned_fasta_path(self):

        output_path = os.path.join(self.get_tmp_path(), 'unaligned_region.fa')
        logging.debug('Region: writing unaligned fasta to %s' % output_path)

        assert self.reads

        with open(output_path, 'w',  buffering=100 * 1024 * 1024) as output_handler:

            # Cut and write references
            human_positions = []
            pathogen_positions = []
            for length, identifier, start, end in self.get_sorted_reference_positions():
                if ';Homo_sapiens;' in identifier:
                    human_positions.append((identifier, start, end))
                else:
                    pathogen_positions.append((identifier, start, end))

            for identifier, start, end in pathogen_positions + human_positions:
                record = self.sequence_proxy.get_reference_record(identifier)
                if record:
                    record = SeqRecord(record.seq, record.id, '', '')
                    #record.id += ';%i-%i' % (start, end)
                    SeqIO.write(
                        [record[start - 1:end]], output_handler, "fasta")
                else:
                    pass
                    #logging.debug('Region: could not retrieve reference %s' % identifier)

            # Write full-length transcripts
            for identifier in self.transcripts:
                record = self.sequence_proxy.get_transcript_record(identifier)
                if record:
                    record = SeqRecord(record.seq, record.id, '', '')
                    SeqIO.write([record], output_handler, "fasta")
                else:
                    logging.debug(
                        'Region: could not retrieve reference %s' % identifier)

            # Write reads
            for read_id in sorted(self.reads):
                record = self.sequence_proxy.get_read_record(read_id)
                if record:
                    # Transform read ID so that LASTZ can work with it (it
                    # makes nicknames based on some characters )
                    identifier = record.id.replace(
                        ':', '_').replace('|', '_').replace('/', '_')
                    record = SeqRecord(record.seq, identifier, '', '')
                    SeqIO.write([record], output_handler, "fasta")
                else:
                    logging.debug(
                        'Region: could not retrieve read %s' % read_id)

            self.unaligned_fasta_path = output_path

        return self.unaligned_fasta_path

    def get_unaligned_fasta_path(self):

        if not self.unaligned_fasta_path:
            self._get_unaligned_fasta_path()

        return self.unaligned_fasta_path

    def _align_fastas(self, aligner, assert_record=None):

        # Write reference fasta
        queries_path = self.get_unaligned_fasta_path()
        first_record = SeqIO.parse(open(queries_path, "rU"), "fasta").next()
        reference_path = os.path.join(self.get_tmp_path(), 'reference.fa')

        logging.debug(
            'Region: writing longest reference to %s' % reference_path)

        SeqIO.write([first_record], reference_path, "fasta")

        # Do alignment
        bam_path = os.path.join(self.get_tmp_path(), 'aligned.bam')
        logging.debug(
            'Region: starting alignment using queries %s to sam path %s' %
            (queries_path, bam_path))

        aligner.align_to_bam_file(
            reference_path, queries_path, bam_path, assert_record=assert_record)

        return bam_path

    def _build_alignment_consensus(self, consensus_builder, alignment_bam_path):

        consensus_path = os.path.join(self.get_tmp_path(), 'consensus.fa')
        logging.debug('Region: writing consensus to %s' % consensus_path)

        consensus_builder.run(alignment_bam_path, consensus_path)

        return consensus_path

    def get_alignment_bam_path(self, aligner, assert_record=None):

        return self._align_fastas(aligner, assert_record)

    def get_consensus_fasta_path(self, consensus_builder, alignment_bam_path):

        return self._build_alignment_consensus(consensus_builder, alignment_bam_path)

    def get_consensus_figure_path(self, jalview_runner, consensus_fasta_path):

        png_path = os.path.join(self.get_tmp_path(), 'consensus_figure.png')

        return jalview_runner.run(consensus_fasta_path, png_path)

    def overlaps(self, other, max_gap_length):

        overlap = (set(self.references) & set(other.references))

        # Overlap is solely based on non-human references
        overlap = [ref for ref in overlap if not ';Homo_sapiens;' in ref]

        if not overlap:
            return False

        for reference in overlap:
            reflist = self.references[reference]
            for start, end in reflist:
                for other_start, other_end in other.references[reference]:
                    distance = self._distance((start, end),
                                             (other_start, other_end))
                    if distance <= max_gap_length:
                        return True
        return False

    def merge(self, other):

        for reference, reflist in other.references.iteritems():
            for (start, end) in reflist:
                self.add_reference(reference, start, end)

        for read in other.reads:
            self.add_read(read)

    def clean_references(self, max_gap_length):

        for reference, reflist in self.references.iteritems():

            start_list = sorted(reflist)
            saved = list(start_list[0])
            result = set()

            for item in start_list:
                if self._distance(saved, item) <= max_gap_length:
                    saved[1] = max(saved[1], item[1])
                else:
                    result.add(tuple(saved))
                    saved[0] = item[0]
                    saved[1] = item[1]
            result.add(tuple(saved))
            self.references[reference] = result

    def to_sorted_record_ids(self):

        references = []
        for name, locations in self.references.iteritems():
            for start, end in locations:
                references.append(end - start, name)

        reads = defaultdict(list)
        for name, locations in self.references.iteritems():
            condition = name.split(';')[1]
            for start, end in locations:
                reads[condition].append(end - start, name)

        all_record_ids = sorted(references, reverse=True)
        for condition, the_reads in reads.iteritems():
            all_record_ids += sorted(the_reads, reverse=True)

        for length, record_id in all_record_ids:
            yield record_id

    def __len__(self):

        return len(self.reads)

    def __cmp__(self, other):

        if self.get_longest_reference_length() <\
                other.get_longest_reference_length():
            return -1
        elif self.get_longest_reference_length() ==\
                other.get_longest_reference_length():
            return 0
        return 1

    def get_longest_reference_length(self):

        if self.length is not None:
            return self.length

        positions = self.get_sorted_reference_positions()
        if positions:
            self.length = self.get_sorted_reference_positions()[0][0]
        else:
            self.length = 0

        return self.length

    def get_sorted_reference_positions(self):

        if self.sorted_reference_positions is not None:
            return self.sorted_reference_positions

        lengths = []
        for reference, the_set in self.references.iteritems():
            for start, end in the_set:
                lengths.append([end - start, reference, start, end])

        self.sorted_reference_positions = sorted(lengths, reverse=True)

        return self.sorted_reference_positions

    def __str__(self):
        return '<Region> of length %10i with %10i reads and %5i references' %\
            (self.get_longest_reference_length(),
             len(self.reads), len(self.references))


class RegionStatistics:

    def __init__(self, input_dir):

        self.input_dir = input_dir
        self.stats = []

    def _add_sequence(self, qualified_name, region_id,
                      longest_human_reference, longest_pathogen_reference, record):

        state = 'pathogen'
        if longest_human_reference:
            state = 'ambiguous'

        reference_type = qualified_name.split('_')[0]
        family_name = '_'.join(qualified_name.split('_')[1:])

        read_type, sample_id, number_reads, basepairs = record.id.split(';')

        coverage = int(basepairs) / float(len(longest_pathogen_reference))

        self.stats.append(
            [reference_type, family_name, region_id, sample_id, state,
             number_reads, str(len(longest_pathogen_reference)), basepairs, '%.5f' % coverage])

    def run(self, output_file):

        output_file = os.path.abspath(os.path.expandvars(output_file))

        if os.path.dirname(output_file) and not os.path.exists(os.path.dirname(output_file)):
            logging.debug('Making directories for output file %s' % output_file)
            os.makedirs(os.path.dirname(output_file))

        for qualified_family_name in os.listdir(self.input_dir):

            family_dir = os.path.join(self.input_dir, qualified_family_name)

            if not os.path.isdir(family_dir):
                continue

            consensus_regions = glob.glob(os.path.join(
                family_dir, 'region_*_consensus.fa'))

            for consensus_region in consensus_regions:

                base_bame = os.path.basename(consensus_region)

                region_id = base_bame.split('_')[1]

                reads = []
                longest_human_reference = None
                longest_pathogen_reference = None

                for consensus_record in SeqIO.parse(consensus_region, 'fasta'):
                    if consensus_record.id.lower().startswith('read'):
                        # ends with read number; cumulative bases
                        reads.append(consensus_record)
                    elif ';Homo_sapiens;' in consensus_record.id:
                        if longest_human_reference is None or\
                                len(consensus_record) > longest_human_reference:
                            longest_human_reference = consensus_record
                    else:
                        if longest_pathogen_reference is None or\
                                len(consensus_record) > longest_pathogen_reference:
                            longest_pathogen_reference = consensus_record

                for read in reads:
                    self._add_sequence(qualified_family_name, region_id,
                                       longest_human_reference, longest_pathogen_reference, read)

        with open(os.path.abspath(os.path.expandvars(output_file)), 'w') as output_handler:
            output_handler.write(
                '\t'.join(['reference_type', 'family_name', 'region_id',
                        'sample_id', 'taxonomic_origin', 'number_reads',
                        'region_length', 'basepairs', 'coverage']) + '\n')
            for entries in sorted(self.stats):
                line = '\t'.join(entries) + '\n'
                output_handler.write(line)


@CLI.subcommand("regions")
class RegionRunner(cli.Application):

    """ Derive homologous groups and regions from hit files """

    hit_files = cli.SwitchAttr(
        ['-v', '--virana_hits'], str, list=True, mandatory=True,
        help="Add hit file for analysis. May be supplied multiple times. May contain globbing characters that are expanded to matchin file names.")

    lastz_path = cli.SwitchAttr(['-z', '--lastz_path'], str, mandatory=False,
                                help="Path to lastz executable",
                                default='lastz')

    jalview_jar_dir = cli.SwitchAttr(
        ['-j', '--jalview_jar_dir'], str, mandatory=False,
        help="Directory containing the jalview.jar file",
        default=None)

    references_path = cli.SwitchAttr(
        ['-r', '--references'], str, mandatory=True,
        help="Input fasta file containing genomic references of pathogens and possibly of human chromosomes as well (the latter is experimental)",
        default='')

    cdna_path = cli.SwitchAttr(['-c', '--cdna'], str, mandatory=False,
                               help="Input fasta file containing human cDNA records.",
                               default=None)

    output_dir = cli.SwitchAttr(['-o', '--output_dir'], str, mandatory=True,
                                help="Output directory that will be filled with subdirectories corresponding to homologous regions. The directory will be generated if it is not existing.",
                                default='')

    tmp_dir = cli.SwitchAttr(['-t', '--tmp_dir'], str, mandatory=False,
                             help="Directory in which temorary files are stored. Temporary files will be stored in subdirectories of the provided directory and will be deleted after completion.",
                             default=None)

    stat_path = cli.SwitchAttr(['-s', '--region_stats'], str, mandatory=False,
                               help="Path to output statistics tab-delimited text file.",
                               default='')

    reference_database_filter = cli.SwitchAttr(
        ['-b', '--reference_database_filter'], str, list=True, mandatory=False,
        help="Specifies which kind of reference databases are considered when extracting hits from hit files. May be specified multiple times. If any are specified, all reference databases not specified are filtered out. By default, this parameter is empty.",
        default=[])

    pathogen_family_filter = cli.SwitchAttr(
        ['-f', '--pathogen_family_filter'], str, list=True, mandatory=False,
        help="Specifies which kind of pathogen families are considered when extracting hits from hit files. May be specified multiple times. If any are specified, all families not specified are filtered out. By default, this parameter is empty.",
        default=[])

    min_read_number = cli.SwitchAttr(
        ['-m', '--min_read_number'], cli.Range(1, 1000), mandatory=False,
        help="Minimum number of reads that are required to be present in homologous region. Regions with fewer reads will be omitted from the results.",
        default=5)

    max_gap_length = cli.SwitchAttr(
        ['-l', '--max_gap_length'], cli.Range(1, 1000), mandatory=False,
        help="Maximum number bases that two candidate homologous regions are distant from each other with regard to their positions on a common reference sequence in order for being eligable for merging.",
        default=50)

    min_region_length = cli.SwitchAttr(
        ['-x', '--min_region_length'], cli.Range(1, 1000), mandatory=False,
        help="Minimum number bases of the longest reference sequence of each homologous region that is generated. Shoer regions will be omitted from the results.",
        default=100)

    word_length = cli.SwitchAttr(
        ['-w', '--word_length'], cli.Range(1, 21), mandatory=False,
        help="Word length of the lastz alignment process. Shorter word lengths allow more sensitive but less specific alignments.",
        default=7)

    ambiguity_cutoff = cli.SwitchAttr(
        ['-a', '--ambiguity_cutoff'], float, mandatory=False,
        help="Ratio of variant positions within a column of the homologous region sequence alignment so that the corresponding consensus sequence at this positions show a ambiguous base instead of the majority base.",
        default=0.3)

    debug = cli.Flag(
        ["-d", "--debug"], help="Enable debug information")

    def main(self):

        if self.debug:
            logging.getLogger().setLevel(logging.DEBUG)

        # Make sequence proxy for managing hit files, references, and cdna
        proxy = SequenceProxy(self.references_path, self.cdna_path)

        for hit_file in self.hit_files:
            hit_file = os.path.expandvars(hit_file)
            for expanded in glob.glob(hit_file):
                proxy.add_hit_file(expanded)

        # Generate homologous groups
        generator = GroupGenerator(proxy,
                                   reference_database_filter=self.reference_database_filter,
                                   pathogen_family_filter=self.pathogen_family_filter,
                                   tmp_dir=self.tmp_dir)
        groups = generator.get_groups(min_read_number=self.min_read_number)

        # Prepare analysis modules ('runners') for postprocessing homologous
        # regions
        consensus = ConsensusRunner(ambiguity_cutoff=self.ambiguity_cutoff)

        if not which(self.lastz_path) or not which(self.lastz_path):
            logging.error('Invalid path to lastz: %s' % self.lastz_path)
            sys.exit(1)

        lastz = LastzRunner(self.lastz_path, word_length=self.word_length)

        if self.jalview_jar_dir:

            jalview_path = os.path.join(self.jalview_jar_dir, 'jalview.jar')
            if not os.path.isfile(jalview_path):
                logging.error('Invalid path to jalview: %s' % jalview_path)
                sys.exit(1)

            jalview = JalviewRunner(
                jalview_jar_dir=self.jalview_jar_dir, tmp_dir=self.tmp_dir)

        else:
            jalview = None

        # Make homologous regions within each homologous group
        for name, group in groups.iteritems():
            group.merge_regions(max_gap_length=self.max_gap_length)
            group.filter_regions(
                min_region_length=self.min_region_length, min_read_number=self.min_read_number)

            if self.cdna_path:
                group.add_transcripts_to_regions()

            group.write_outputs(self.output_dir, lastz, consensus, jalview)
            group._delete_temporary_dir()

        # Run statistics on all homologous regions
        if self.stat_path:
            statistics = RegionStatistics(self.output_dir)
            statistics.run(self.stat_path)

if __name__ == "__main__":
    CLI.run()
