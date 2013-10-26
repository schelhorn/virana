#!/usr/bin/env python

import sys
import urllib2
import hashlib
import random
import logging
import string
import gzip
import os

from random import choice, randint, shuffle

try:
    from plumbum import cli
except ImportError:
    message = 'This script requires the plumbum python package\n'
    sys.stderr.write(message)
    sys.exit(1)

logging.basicConfig(level=logging.INFO, format='%(message)s')

class Sequence:

    def __init__(self, bases='', accession=None):

        self.bases = bases

        if accession:

              logging.debug('Obtaining sequence %s from genbank' % accession)
              url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' +\
                      'db=nucleotide&rettype=fasta&id=%s&email=no@reply.com&tool=python' % accession
              try:
                  response = urllib2.urlopen(url)
                  fasta_lines = response.read().split('\n')
                  self.bases = ''.join(fasta_lines[1:])
              except urllib2.HTTPError:
                  logging.warning(
                      'Warning: Unknown identifier %s, bases not stored' % accession)

    def crop_to_random_region(self, region_length):

        start = randint(0, len(self.bases) - region_length - 1)
        end = start + region_length
        logging.debug('Cropping sequence to positions %i-%i' % (start, end))

        self.bases = self.bases[start:end]

    def insert_sequence(self, dna_insert, at_position=None):

        if at_position is None:
            start = randint(len(dna_insert.bases), len(
                self.bases) - len(dna_insert.bases) - 1)
        else:
            start = at_position
        logging.debug('Inserting sequence at position %i' % start)

        self.bases = self.bases[:start] + dna_insert.bases + self.bases[start:]

    def __add__(self, other):

        return Sequence(self.bases + other.bases)

    def __len__(self):

        return len(self.bases)

    def copy(self):

        return Sequence(self.bases[:])


class ReadSimulator:

    def __init__(self, read_length, inner_distance_m, inner_distance_s,\
                    substitution_rate=0.01, insertion_rate=0.001, deletion_rate=0.001):

        self.read_length = read_length
        self.inner_distance_m = inner_distance_m
        self.inner_distance_s = inner_distance_s
        self.substitution_rate = substitution_rate
        self.insertion_rate = insertion_rate
        self.deletion_rate = deletion_rate

        self.trans_complement = string.maketrans('ATCGN', 'TAGCN')

    def get_complement(self, dna_string):

        return dna_string.translate(self.trans_complement)[::-1]

    def get_reads_iterator(self, genome, coverage):

        len_genome = len(genome)
        number_bases = len_genome * coverage
        pairs_needed = number_bases / (2 * self.read_length)

        pairs_generated = 0

        logging.debug(
            'Simulating %i paired end reads of length %i to obtain coverage %i' %
            (pairs_needed, self.read_length, coverage))

        while pairs_generated < pairs_needed:

            end1_start = randint(0, len_genome - 1)
            end1_end = end1_start + self.read_length

            inner_distance = int(random.normalvariate(
                self.inner_distance_m, self.inner_distance_s))

            end2_start = max(end1_start, end1_end + inner_distance)
            end2_end = end2_start + self.read_length

            if end2_end > (len_genome - 1):
                continue

            end_1 = genome.bases[end1_start:end1_end + 1 + self.read_length]
            end_2 = self.get_complement(
                genome.bases[end2_start - self.read_length:end2_end + 1])

            end_1_mutated = []
            for base in end_1:
                base_call = base
                if random.random() <= self.substitution_rate:
                    base_call = random.choice('ACGT'.replace(base, ''))
                if random.random() <= self.insertion_rate:
                    base_call += random.choice('ACGT')
                if random.random() <= self.deletion_rate:
                    base_call = base_call[:-1]
                end_1_mutated.append(base_call)

            end_2_mutated = []
            for base in end_2:
                base_call = base
                if random.random() < self.substitution_rate:
                    base_call = random.choice('ACGT'.replace(base, ''))
                if random.random() < self.insertion_rate:
                    base_call += random.choice('ACGT')
                if random.random() < self.deletion_rate:
                    base_call = base_call[:-1]

                end_2_mutated.append(base_call)

            finished_end_1 = ''.join(end_1_mutated)[:self.read_length]
            finished_end_2 = ''.join(end_2_mutated)[:self.read_length]

            pairs_generated += 1

            yield finished_end_1, finished_end_2


class ReadWriter:

    def __init__(self, pathname_end1, pathname_end2):

        self.pathname_end1 = pathname_end1
        self.pathname_end2 = pathname_end2

        self.sequence_iterators = []

    def add_reads_iterator(self, iterator):

        self.sequence_iterators.append(iterator)

    def write_fastqs(self):

        all_read_pairs = []

        for sequence_iterator in self.sequence_iterators:

            # Make random library identifier
            library_identifier = ''.join(
                choice(string.lowercase) for x in range(20))

            for i, (end_1, end_2) in enumerate(sequence_iterator):

                # Make random read identifier
                read_identifier = hashlib.sha1(
                    library_identifier + str(i)).hexdigest()

                all_read_pairs.append((read_identifier, end_1, end_2))

        logging.debug('Shuffling reads')
        # Suffle reads
        shuffle(all_read_pairs)

        # Write reads to file
        with gzip.GzipFile(self.pathname_end1, 'w') as handler_end1:
            with gzip.GzipFile(self.pathname_end2, 'w') as handler_end2:

                logging.debug('Writing reads')
                for identifier, end_1, end_2 in all_read_pairs:
                    handler_end1.write(
                        '@%s/1\n%s\n+\n%s\n' % (identifier, end_1, 'I' * len(end_1)))
                    handler_end2.write(
                        '@%s/2\n%s\n+\n%s\n' % (identifier, end_2, 'I' * len(end_2)))


class CLI(cli.Application):
    """Simulates metagenomic short sequence reads in FASTQ format."""
    PROGNAME = "vsim"
    VERSION = "1.0.0"
    DESCRIPTION = \
"""DESCRIPTION: Virana vsim - simulates metagenomic short read data.

The Virana simulation utility ('vsim') downloads reference genomes based on
provided gene identifiers (GIs) and generates simulated, paired-end FASTQ reads
from these genomes at user-specified coverages and inner distances. The output
of this module can be employed as test data sets for later Virana analysis stages.

https://github.com/schelhorn/virana

Schelhorn S-E, Fischer M, Tolosi L, Altmueller J, Nuernberg P, et al. (2013)
Sensitive Detection of Viral Transcripts in Human Tumor Transcriptomes.
PLoS Comput Biol 9(10): e1003228. doi:10.1371/journal.pcbi.1003228"""

    USAGE = """USAGE: The program has one mode that can be accessed by
       [vsim | python vsim.py] reads
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


@CLI.subcommand("reads")
class Reads(cli.Application):
    """ Simulates metagenomic short read data."""

    references = cli.SwitchAttr(['-r', '--reference_identifiers'],
                                    str, list=True, mandatory=True,
                                    help="GI:coverage that specifies the reference sequence from which reads should be generated as well as the sequencing coverage. Argument may be supplied multiple times. Example: '9626196:100' will generate 100-fold coverage of the rous sarcoma virus")

    output_path  = cli.SwitchAttr(['-o', '--output_path'], str, mandatory=True,
                                    help="Basename of output fastq file. Note that two zipped fastq files will be generated, one for each end. For example, if the basename (this argument) is /home/user/myfastqs, the the two files /home/user/myfastqs_1.fq.gz and /home/user/myfastqs_2.fq.gz will be generated. If not present, all necessary subdirectories will also be generated")

    read_length     = cli.SwitchAttr(['-l', '--read_length'],
                                    int, default=150,
                                    help="Specifies read length the paired read ends, in basepairs. Read length here denotes the lengths of the single segments, so that a value of 150 would result in 2x150 bp reads")

    inner_distance_m = cli.SwitchAttr(['-m', '--mean_inner_distance'],
                                    int, default=200,
                                    help="Specifies the mean (mu) of inner distance of the paired read ends, in basepairs. Note that this is not identical to the insert size, which commonly includes the lengths of the reads")

    inner_distance_s = cli.SwitchAttr(['-s', '--sdev_inner_distance'],
                                    float, default=20.0,
                                    help="Specifies the standard deviation (sigma) of inner distance of the paired read ends, in basepairs")

    substitution_rate = cli.SwitchAttr(['--substitution_rate'],
                                    float, default=0.01,
                                    help="Specifies the substitution error rate of the read simulation process. A value of 0.01 corresponds to 1% of randomly substituted bases")

    insertion_rate = cli.SwitchAttr(['--insertion_rate'],
                                    float, default=0.001,
                                    help="Specifies the insertion error rate of the read simulation process. A value of 0.001 corresponds to 0.1% chance of a randomly inserted base per base of the read")

    deletion_rate = cli.SwitchAttr(['--deletion_rate'],
                                    float, default=0.001,
                                    help="Specifies the deletion error rate of the read simulation process. A value of 0.001 corresponds to 0.1% chance of a randomly deleted base per base of the read")

    debug       = cli.Flag(["-d", "--debug"], help="Enable debug messages")

    def main(self):
        """ Randomly generated paired end reads """

        random.seed(0)

        if self.debug:
            logging.getLogger().setLevel(logging.DEBUG)

        read_1_path = self.output_path + '_1.fq.gz'
        read_2_path = self.output_path + '_2.fq.gz'
        dirname = os.path.dirname(read_1_path)

        if dirname and not os.path.exists(dirname):
            logging.debug('Making subdirectory %s' % dirname)
            os.makedirs(dirname)

        writer      = ReadWriter(read_1_path, read_2_path)
        simulator   =  ReadSimulator(read_length=self.read_length, inner_distance_m=self.inner_distance_m,\
                                             inner_distance_s=self.inner_distance_s, substitution_rate=self.substitution_rate,\
                                             insertion_rate=self.insertion_rate, deletion_rate=self.deletion_rate)

        for reference in self.references:
            try:
                gi, coverage = map(int, reference.split(':'))
            except ValueError:
                logging.error('Invalid reference definition %s' % reference)
                sys.exit(1)

            logging.debug('Obtaining reference sequence for gi:%i' % gi)

            sequence    = Sequence(accession=str(gi))

            logging.debug('Reference sequence has length %i' % len(sequence))

            reads       = simulator.get_reads_iterator(sequence, coverage=int(coverage))

            writer.add_reads_iterator(reads)

        logging.debug('Generating simulated reads and writing to %s/%s' % (read_1_path, read_2_path))

        writer.write_fastqs()

if __name__ == "__main__":
    CLI.run()


