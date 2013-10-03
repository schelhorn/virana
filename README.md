Virana
======

The Virus Analysis Toolkit

Note: We are currently expanding this Readme in order to include analysis recommendations for various scenarios as well as a detailled installation guide. Please visit this page again in a week for an illustrated guide and additional command line help.  A version of this tool for the workflow engine `Galaxy` that automatically satisfies all dependencies during installation is also available at http://testtoolshed.g2.bx.psu.edu/view/mzeidler/virana_main

### Summary

Virana, _the virus analysis toolkit_, is a Python-based software toolkit for analyzing metatranscriptomic (and, to a degree, metagenomic) sequence data in a context of clinical metagenomics in order to:

- identify pathogen nucleotide sequences in transcriptomic (and, while less tested, genomic) short read data
- identify pathogen transcripts diverged from known references of various microbial species
- identify transcripts expressed at low abundances across multiple samples
- identify transcripts homologous to human factors such as certain viral oncogenes

A recent analysis of human tumor transcriptomes displays the ability of Virana to detect diverged viral nucleotide sequences and human-viral homologs with high sensitivity. The manuscript is available at doi:10.1371/journal.pcbi.1003228 and http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1003228 

### Components

Depending on the branch you choose, virana is available either in a locally executable version or as a tool for the workflow engine  `Galaxy`. In the following, we assume that you want to use the locally executable, command-line version. Virana is implemented in Python and requires the executables `samtools`, `STAR` in the path or provided as additional command line arguments. The jar file of the visualization tool `jalview` may be supplied to automatically generate alignment diagrams.

Virana consists of three components:

1. `vref` ( _virana reference_ ), a tool that obtains viral, human, fungi, and bacterial references sequences as well as their taxonomic annotation from reference archives
2. `vmap` ( _virana mapper_ ), a tool that employs the short read mappers `STAR` and `bwa mem` to map transcriptomic reads and genomic reads against annotated references, respectively
3. `vhom` ( _virana homology_ ), a tool that analyzes the homology relationships within mapped reads in order to extract _homologous regions_, i.e. nucleotide stretches that display high sequence similarity to a pathogen and, optionally, also to human factors

### Code example

A typical analysis may consists of three steps. First, a reference database is obtained in fasta format that contains the human reference sequences, all known viral references, as well as known contaminants from the UniVec database.

```shell
  vref fasta -o my_reference_db.fa -r Homo_sapiens -r UniVec -r Viruses
```
Other options for references are `rRNA`, `Fungi`, `Plasmids`, `Protozoa`, `Homo_sapiens_cDNA`, `Bacteria`

For later analyses, we may also want to obtain human transcripts (coding as well as noncoding) in order to compare them with metatranscriptomic reads

```shell
  vref fasta -o my_transcriptome_db.fa -r Homo_sapiens_cDNA
```

Next, we map short read RNA-Seq data against the reference database by first generating an index for the short read mapper `STAR`, using multiple compute threads, and then map reads in fastq format. As an output, we obtain both a standard `bam` file as well as a special file use in later steps that better captures homology relationships between multi-mapping reads. In this example, we restrict our output to reads that are alignable to viral reference sequences while conserving additional mappings of this reads to other (e.g. human) references. 

```shell
  vmap rnaindex -i my_index_dir -r my_reference_db.fa --threads 16
```

```shell
  vmap rnamap -i my_index_dir -f Viruses -r reads_end_1.fq.gz -r reads_end_2.fq.gz --zipped -b bam_output.bam --threads 16 -v virana_output.bz2
```

Here, we may also mention additional outputs; for example, providing `-x taxonomy.txt` enables output of taxonomic information of read hits during the mapping `--unmapped_end_1` allows for outputting unmapped reads as fastq files, and `--chimeric_mappings` activates output of chimerically mapping reads in sam format.

Last, the mappings are analyzed within their homologous context and regions that correspond to factors that are either of clearly pathogen origin OR that are homologs to both human and pathogen are reconstructed from the mappings in a greedy micro-assembly. These regions allow for either automatic or manual analysis in order to determine the true origin of the reads.

```shell
  vhom regions --references my_reference_db.fa --cdna my_transcriptome_db.fa --output_dir viral_families --region_stats output_statistics.txt --virana_hits virana_output.bz2
```

Within the output directory `viral_families`, subdirectories are generated for each viral taxonomic family that the reads you specified earlier map to. Each family directory, in turn, contains a number of homologous regions that are represented by aligned, multi-fasta consensus sequences and visualizations of these consensus sequences using `jalview`.

### Availability of other data

Validation data that has been used in the PLOS Computational Biology manuscript that introduces Virana are not available on GitHub due to size limitations. Many of these files are available in public archives as indicated in the manuscript. Please contact us at sven@mpi-inf.mpg.de for access to the remaining files not available elsewhere.
