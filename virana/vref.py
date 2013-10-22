#!/usr/bin/env python

""" Virana reference tool for downlaoding taxonomically annotated
    reference genomes. Part of the Virana package.

    (c) 2013, Sven-Eric Schelhorn, MPI for Informatics.
"""


import sys
import logging
import gzip
import string
import os
import tarfile

from io import BytesIO
from collections import defaultdict

try:
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
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
    import ftputil
except ImportError:
    message = 'This script requires the ftputil python package\n'
    sys.stderr.write(message)
    sys.exit(1)


NON_ID = ''.join(c for c in map(chr, range(256)) if not c.isalnum())
NON_ID = NON_ID.replace('_', '').replace('-', '')
NUC_TRANS = string.maketrans('uryswkmbdhvURYSWKMBDHV', 'tnnnnnnnnnnTNNNNNNNNNN')

logging.basicConfig(level=logging.INFO, format='%(message)s')


class DataDownloader(object):
    """ Generic class for downloading fasta records via FTP. """

    def __init__(self, host_name, base_path, user='anonymous', password='anonymous'):

        self.group     = 'NoGroup'
        self.host_name = host_name
        self.base_path = base_path

        self.host = ftputil.FTPHost(host_name, user, password)
        self.host.chdir(base_path)

    @classmethod
    def _get_fasta_record(self, record, group, family, organism, identifier, sub_identifier=None, description=''):

        record.seq._data = record.seq._data.translate(NUC_TRANS)

        group       = group.translate(None, NON_ID)
        family      = family.translate(None, NON_ID)
        organism    = organism.translate(None, NON_ID)
        identifier  = identifier.translate(None, NON_ID)

        if sub_identifier:
            sub_identifier = sub_identifier.translate(None, NON_ID)
            record.id  = ';'.join([group, family, organism, identifier, sub_identifier])
        else:
            record.id           = ';'.join([group, family, organism, identifier])

        record.description  = description
        record.name         = ''

        return record

    def get_fasta_records(self):

        raise NotImplementedError

class SilvaDownloader(DataDownloader):
    """ Downlaods rRNA transcripts from the Silva database. """

    def __init__(self, silva_release='111'):

        self.silva_release      = silva_release
        self.silva_host_name    = 'ftp.arb-silva.de'
        self.base_path          = '/release_%s/Exports/' % self.silva_release

        super(SilvaDownloader, self).__init__(self.silva_host_name, self.base_path)

        self.group = 'rRNA'

    def get_fasta_records(self):

        remote_path = self.base_path + 'SSURef_%s_NR_tax_silva_trunc.fasta.tgz' % self.silva_release
        if self.host.path.isfile(remote_path):

            with self.host.file(remote_path, 'rb') as remote_handle:
                remote_content = BytesIO(remote_handle.read())

                with tarfile.open(fileobj=remote_content) as tar:

                    for subentry in tar.getnames():
                        if subentry.endswith('.fasta'):
                            logging.debug('Processing rRNA file %s' % subentry)
                            subhandle = tar.extractfile(subentry)

                            for record in SeqIO.parse(subhandle, "fasta"):
                                identifier  = record.id.split('.')[0]
                                description = '_'.join(record.description.split(' ')[1:])
                                taxonomy    = description.split(';')
                                organism    = taxonomy[-1]
                                family      = 'Unassigned'

                                for level in taxonomy:
                                    if level.endswith('eae'):
                                        family = level
                                        break

                                yield self._get_fasta_record(record, self.group, family, organism, identifier)

class HumanTranscriptDownloader(DataDownloader):
    """ Downloads human transcripts fasta records that contain the locations
        of their exons on the human genome within their description line.
    """

    def __init__(self):

        self.ensemble_host_name  = 'ftp.ensembl.org'
        self.gtf_path            = '/pub/current_gtf/homo_sapiens'
        self.cdna_path           = '/pub/current_fasta/homo_sapiens/cdna'
        self.ncrna_path          = '/pub/current_fasta/homo_sapiens/ncrna'

        super(HumanTranscriptDownloader, self).__init__(self.ensemble_host_name, self.gtf_path)

        self.chromosomes    = set()
        self.genes          = defaultdict(set)
        self.regions        = defaultdict(set)
        self.group          = 'Homo_sapiens_cDNA'

    def _cache_annotations(self):
        """ Obtains gtf annotations for human transcripts from ensemble and
            caches them in memory.
        """

        self.host.chdir(self.gtf_path)

        entries = self.host.listdir(self.host.curdir)

        for entry in entries:

            if self.host.path.isfile(entry):
                if entry.startswith('Homo_sapiens.') and entry.endswith('.gtf.gz'):
                    logging.debug('Processing cDNA annotations %s' % entry)

                    with self.host.file(self.gtf_path + '/' + entry, 'rb') as zipped_handle:
                        remote_content = BytesIO(zipped_handle.read())

                        gzipfile = gzip.GzipFile(fileobj=remote_content)
                        for i, gtf_line in enumerate(gzipfile):
                            if i % 100000 == 0:
                                logging.debug('Cached %i human transcript annotations' % i)

                            self._add_annotation(gtf_line)
                        gzipfile.close()

    def _get_raw_transcripts(self):
        """ Obtains coding and noncoding human transcripts from ensemble and
            yields them as fasta records if they have exons that have a known location
            on the human genome. Fasta records contain the exon locations in their header
            lines.
        """

        for subpath in [self.cdna_path, self.ncrna_path]:

            self.host.chdir(subpath)

            entries = self.host.listdir(self.host.curdir)

            for entry in entries:
                if self.host.path.isfile(entry):
                    if entry.endswith('cdna.all.fa.gz') or entry.endswith('ncrna.fa.gz'):

                        logging.debug('Processing human transcript file %s' % entry)

                        with self.host.file(subpath + '/' + entry, 'rb') as zipped_handle:
                            remote_content = BytesIO(zipped_handle.read())

                            gzipfile = gzip.GzipFile(fileobj=remote_content)
                            for i, record in enumerate(SeqIO.parse(gzipfile, "fasta")):

                                if record.id.startswith('ENST'):

                                    if i % 100000 == 0:
                                        logging.debug('Retrieved %i human transcripts' % i)

                                    record.description = ''

                                    yield record
                            gzipfile.close()

    def _add_annotation(self, gtf_line):
        """ Parses gtf annotations and stores them in memory. """

        if gtf_line[0] == '#':
            return

        fields = gtf_line.strip().split('\t')
        chromosome, source, feature, start, end, score, strands, frame, attributes\
                                                                    = fields

        start, end = sorted([int(start), int(end)])

        if feature != 'exon':
            return

        attributes = attributes.replace(' ', '').split(';')
        transcript_ids = [identifier.split('"')[1] for identifier in attributes\
                                if identifier.startswith('transcript_id')]

        gene_names = [identifier.split('"')[1] for identifier in attributes\
                                if identifier.startswith('gene_name')]

        for transcript_id in transcript_ids:
            self.genes[transcript_id] = self.genes[transcript_id].union(gene_names)
            self.regions[transcript_id].add((chromosome, start, end))

    def get_fasta_records(self):
        """ Yields annotated fasta records of human transcripts. """

        logging.debug('Caching annotations')
        self._cache_annotations()

        logging.debug('Downloading and annotating transcripts')
        for record in self._get_raw_transcripts():
            transcript_id   = record.id

            if transcript_id in self.regions:

                mapping_locations = self.regions[transcript_id]

                description = '|'.join([';'.join(map(str, location))\
                        for location in mapping_locations])

                # Obtain gene name for transcript
                gene_names = self.genes[transcript_id]
                if gene_names:
                    gene_name = list(gene_names)[0]
                else:
                    gene_name = transcript_id

                yield self._get_fasta_record(record, self.group , 'Hominidae', 'Homo_sapiens', gene_name, sub_identifier=transcript_id, description=description)

class RefSeqDownloader(DataDownloader):
    """ Generic downloader that known how to interpret and parse genbank records """

    def __init__(self, group):

        self.refseq_host_name  = 'ftp.ncbi.nih.gov'
        self.base_path         = '/genomes/' + group

        super(RefSeqDownloader, self).__init__(self.refseq_host_name, self.base_path)

        self.group              = group


    @classmethod
    def _get_project_from_genbank(self, gb_record):
        """ Parses genbank project identifier from genbank record. """

        project = ''
        for dbxref in gb_record.dbxrefs:
            for entry in dbxref.split(' '):
                if entry.startswith('BioProject'):
                    project = entry.split(':')[1]
                    return project

        if not project:
            for dbxref in gb_record.dbxrefs:
                for entry in dbxref.split(' '):
                    if entry.startswith('Project'):
                        project = entry.split(':')[1]
                        return project

        return project

    @classmethod
    def _get_annotations_from_genbank(self, gb_record):
        """ Parses taxonomic annotation from genbank record. """

        accession   = gb_record.annotations['accessions'][0]
        organism    = gb_record.annotations['organism'].replace(' ', '_')
        organism    = '_'.join(organism.split('_')[:2])
        taxonomy    = ("; ").join(gb_record.annotations['taxonomy'])
        gi          = gb_record.annotations['gi']

        family      = 'Unassigned'
        if len(gb_record.annotations['taxonomy']) > 0:
            for level in gb_record.annotations['taxonomy']:
                if level.endswith('dae') or level.endswith('aceae'):
                    family = level
                    break

        if family == 'Unassigned':
            logging.debug('Unassigned family found based on %s' % taxonomy)

        return organism, accession, gi, taxonomy, family

    @classmethod
    def _get_record_from_genbank(self, gb_record):
        """ Parses sequence from genbank record. """

        return SeqRecord(gb_record.seq)

    @classmethod
    def _parse_genbank_records(self, file_handle):

        for gb_record in SeqIO.parse(file_handle, 'genbank'):

            project = self._get_project_from_genbank(gb_record)
            organism, accession, gi, taxonomy, family =\
                self._get_annotations_from_genbank(gb_record)
            record = self._get_record_from_genbank(gb_record)

            if len(record) > 0:
                yield record, project, taxonomy, family, organism, gi

    def _get_nc(self, entry_path):

        if self.host.path.isfile(entry_path):
            logging.debug('Processing refseq entry %s' % entry_path)
            with self.host.file(entry_path) as remote_handle:

                for record, project, taxonomy, family, organism, gi\
                    in self._parse_genbank_records(remote_handle):

                    cleaned = self._get_fasta_record(record, self.group , family, organism, gi)

                    yield cleaned

    def _get_nz(self, entry_path):

        scaffold_path = entry_path.replace('.gbk', '.scaffold.gbk.tgz')
        if self.host.path.isfile(scaffold_path):
            logging.debug('Processing refseq entry %s' % scaffold_path)
            remote_handle = self.host.file(scaffold_path, 'rb')
            remote_content = BytesIO(remote_handle.read())
            tar = tarfile.open(fileobj=remote_content)
            for subentry in tar.getnames():
                if subentry.endswith('.gbk'):
                    logging.debug('Processing subaccession %s' % subentry)
                    subhandle = tar.extractfile(subentry)
                    for record, project, taxonomy, family, organism, gi in self._parse_genbank_records(subhandle):
                        yield self._get_fasta_record(record, self.group , family, organism, gi)
            tar.close()
            remote_handle.close()

    def get_fasta_records(self):

        species_dirs = self.host.listdir(self.host.curdir)

        for species_dir in species_dirs:

            species_path = self.base_path + '/' + species_dir

            if self.host.path.isdir(species_path):

                self.host.chdir(species_path)

                logging.debug('Processing species directory %s' % species_path)

                if self.host.path.isdir(species_path):

                    self.host.chdir(species_path)

                    entries = self.host.listdir(self.host.curdir)

                    for entry in entries:

                        if self.host.path.isfile(self.host.curdir + '/' + entry):

                            if entry.endswith('.gbk') or entry.endswith('.gbk.gz'):

                                logging.debug('Processing entry %s' % entry)

                                entry_path = species_path + '/' + entry

                                if entry.startswith('NC_'):
                                    for record in self._get_nc(entry_path):
                                        yield record

                                elif entry.startswith('NZ_'):
                                    for record in self._get_nz(entry_path):
                                        yield record

class HumanGenomeDownloader(DataDownloader):

    def __init__(self):

        self.ensemble_host_name  = 'ftp.ensembl.org'
        self.base_path           = '/pub/current_fasta/homo_sapiens/dna'

        super(HumanGenomeDownloader, self).__init__(self.ensemble_host_name, self.base_path)

        self.group          = 'Homo_sapiens'

    def get_fasta_records(self):

        entries = self.host.listdir(self.host.curdir)
        for entry in entries:

            if self.host.path.isfile(entry):

                # You may want to change your filtering criteria here if you would like to include patches and haplotypes
                if entry.endswith('fa.gz') and 'dna_sm.primary_assembly' in entry:

                    logging.debug('Processing human entry %s' % entry)

                    with self.host.file(self.host.getcwd() + '/' + entry, 'rb') as zipped_handle:
                        remote_content = BytesIO(zipped_handle.read())

                        gzipfile = gzip.GzipFile(fileobj=remote_content)

                        for record in SeqIO.parse(gzipfile, "fasta"):

                            identifier      = record.id
                            organism        = 'Homo_sapiens'
                            family          = 'Hominidae'

                            yield self._get_fasta_record(record, self.group, family, organism, identifier)

                        gzipfile.close()

class UniVecDownloader(DataDownloader):

    def __init__(self):

        self.univec_host_name  = 'ftp.ncbi.nih.gov'
        self.base_path         = '/pub/UniVec'

        super(UniVecDownloader, self).__init__(self.univec_host_name, self.base_path)

        self.group          = 'UniVec'

    def get_fasta_records(self):

        logging.debug('Processing Univec entries')

        with self.host.file(self.host.getcwd() + '/' + 'UniVec') as remote_handle:
            for record in SeqIO.parse(remote_handle, "fasta"):
                organism    = '_'.join(record.description.split(' ')[1:])
                family      = 'Unassigned'
                identifier  = '000000000'

                yield self._get_fasta_record(record, self.group, family, organism, identifier)

class CLI(cli.Application):
    """Downloads taxonomically annotated genomic and trasncriptomic reference sequences."""
    PROGNAME = "vref"
    VERSION = "1.0.0"
    DESCRIPTION = \
"""DESCRIPTION: Virana vref - downloads taxonomically annotated reference sequences.

The Virana reference utility ('vref') downloads up-to-date human and microbial
reference sequences for analysis of metagenomic short read data. In addition
to allowing convenient download and pooling of reference genomes and transcriptomes,
vref employs a simple but effective taxonomic annotation scheme that is used
by later stages of the Virana pipeline.

https://github.com/schelhorn/virana

Schelhorn S-E, Fischer M, Tolosi L, Altmueller J, Nuernberg P, et al. (2013)
Sensitive Detection of Viral Transcripts in Human Tumor Transcriptomes.
PLoS Comput Biol 9(10): e1003228. doi:10.1371/journal.pcbi.1003228"""

    USAGE = """USAGE: The program has four modes that can be accessed by
       [vref | python vref.py] [fasta, blast] """

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

@CLI.subcommand("fasta")
class References(cli.Application):
    """ Obtains NCBI RefSeq and Ensembl reference genomes and transcriptomes."""

    valid_references = ['Fungi', 'Fungi_DRAFT', 'Bacteria',
              'Bacteria_DRAFT', 'Homo_sapiens',
              'Viruses', 'UniVec', 'Plasmids',
              'Protozoa', 'rRNA', 'Homo_sapiens_cDNA']

    references = cli.SwitchAttr(['-r', '--reference'],
                                   cli.Set(*valid_references, case_sensitive=True),
                                   list=True, default=valid_references,
                                   mandatory=False,
                                   help="Sets the kind of references to obtain; microbial ('Fungi', 'Fungi_DRAFT', 'Bacteria', 'Bacteria_DRAFT', 'Protozoa', 'Viruses'), 'Plasmids'), human ('Homo_sapiens', 'Homo_sapiens_cDNA'), as well as taxonomically mixed ('rRNA' (from SILVA), 'UniVec') references are available.")

    fasta_path  = cli.SwitchAttr(['-o', '--output_file'], str, mandatory=True,
                                 help="Sets the fasta output file. Note that all downloaded fasta records are stored within a single fasta file.")

    zipped      = cli.Flag(["-z", "--zipped"], help="Write gzipped output")
    debug       = cli.Flag(["-d", "--debug"], help="Enable debug messages")


    def main(self):
        """ Downloads genome fasta files from reference databases"""

        if self.debug:
            logging.getLogger().setLevel(logging.DEBUG)

        if os.path.dirname(self.fasta_path) and not os.path.exists(os.path.dirname(self.fasta_path)):
            logging.debug('Making directories for output file %s' % self.fasta_path)
            os.makedirs(os.path.dirname(self.fasta_path))

        if self.zipped:
            logging.debug('Making zipped output file %s' % self.fasta_path)
            output_file = gzip.open(self.fasta_path, 'wb')

        else:
            logging.debug('Making regular output file %s' % self.fasta_path)
            output_file = open(self.fasta_path, 'w', buffering=100 * 1024 * 1024)


        for reference in self.references:

            if reference == 'UniVec':
                downloader = UniVecDownloader()
            elif reference == 'Homo_sapiens':
                downloader = HumanGenomeDownloader()
            elif reference == 'rRNA':
                downloader = SilvaDownloader()
            elif reference == 'Homo_sapiens_cDNA':
                downloader = HumanTranscriptDownloader()
            else:
                downloader = RefSeqDownloader(group=reference)

            for record in downloader.get_fasta_records():
                SeqIO.write([record], output_file, "fasta")

        output_file.close()

@CLI.subcommand("blast")
class Blast(cli.Application):
    """ Obtains blast databases from NCBI."""

    valid_references = ['nt', 'nr']

    references = cli.SwitchAttr(['-r', '--reference_database'],
                                   cli.Set(*valid_references, case_sensitive=True),
                                   list=True, default=valid_references,
                                   mandatory=False,
                                   help="Sets the kind of database to obtain; argument can be supplied multiple times. Common choices for reference databases are 'nr' and 'nt'.")

    database_path  = cli.SwitchAttr(['-o', '--output_path'], str, mandatory=True,
                                 help="Sets the output directory for the blast database. ")

    debug       = cli.Flag(["-d", "--debug"], help="Enable debug messages")

    def main(self):
        """ Downloads blast database files """

        if self.debug:
            logging.getLogger().setLevel(logging.DEBUG)

        if not os.path.exists(self.database_path):
            logging.debug('Making output directory %s' % self.database_path)
            os.makedirs(self.database_path)

        host_name = 'ftp.ncbi.nih.gov'
        host = ftputil.FTPHost(host_name, 'anonymous', 'anonymous')

        for reference in self.references:

            path = '/blast/db'
            host.chdir(path)

            entries = host.listdir(host.curdir)

            for entry in entries:

                if host.path.isfile(entry):

                    if entry.startswith(reference) and entry.endswith('.tar.gz'):
                        logging.debug('Downloading %s' % entry)
                        local_path = os.path.join(self.database_path, entry)
                        host.download(entry, local_path, 'b')
                        logging.debug('Unpacking %s' % entry)
                        tfile = tarfile.open(local_path, 'r:gz')
                        tfile.extractall(self.database_path)
                        os.remove(local_path)

if __name__ == "__main__":
    CLI.run()

