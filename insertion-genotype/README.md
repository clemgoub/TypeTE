# insertion-genotype
Genotyping Mobile Element Insertions based on remapping reads

This is a streamlined implementation of pipeline for genotyping mobile element insertions
described in Wildschutte et al. 2015 PMID:26503250:

Wildschutte JH, Baron A, Diroff NM, Kidd JM. Discovery and characterization of Alu repeat sequences via precise local read
assembly. Nucleic Acids Res. 2015 Dec 2;43(21):10292-307. doi: 10.1093/nar/gkv1089. Epub 2015 Oct 25.

https://www.ncbi.nlm.nih.gov/pubmed/26503250

## Required Components

Requires python2.X, pysam module, samtools and bwa

samtools must be in your path


The path to the version of bwa to be used is given on the command line. The pipeline is designed for bwa version 0.5.9 but should also work with more recent versions


## Running the pipeline

Step 1: Build mapping index

The genotyping works by mapping reads against reference and alternative haplotype
sequences at each locus.  Thus, for each locus a new reference for mapping is constructed
that contains just those alleles.

The sequence information is stored in a simple table format.  An example file for dog
elements is found in input-files/combined-TE-EXAMPLE.txt

This file consists of 8 columns and includes a header row:

**siteID**	unique identifier for site, this is used for naming directories and loci

**insBP**	insertion position.  This information is not used in genotyping

**gTSDl**  size of target site duplication.  This information is also not used

**alignFragBegin, alignFragEnd**  the 1-based coordinates of the listed sequence flanking an insertion. This is typically the element of interest plus 500 bp.

**genomeSeq**  the full DNA sequence of the reference alelle

**altSeq** the full DNA sequence of the alternative allele

The python script create-alternative-alleles.py does the appropriate setup, making
a set of directories for use in mapping and making bwa-indexed files for each locus.
If there are many thousands of loci this is make take many hours per sample.
For troubleshooting, you can first perform a test on smaller set of loci.


Example cmds for running step 1:

First, make a base directory for organizing files

mkdir genotyping-base

Then run:
python insertion-genotype/create-alternative-alleles.py \
--allelefile  input-files/combined-TE.txt \
--allelebase genotyping-base \
--bwa bwa-0.5.9/bwa

This makes a separate directory for each locus under genotyping-base/locusAlleles/.  The option given to --bwa should point to the executable for the version of bwa you would like to use.

Step 2: Genotyping a BAM file

For each locus, the script process-sample.py extracts read pairs that map nearby and remaps them to the
reference and alternative alleles.  Reads with MAPQ0 (map equally well to each allele) are discarded
and the genotype is determined based on the number of read pairs in support of each allele.

The program performs this in order for all loci by first extracting reads for all loci,
then mapping all loci, then filtering and calculating genotype likelihoods, and finally
writing out a VCF for the sample.  For large sets of loci this can take many hours to complete.

example cmd for step 2:

python insertion-genotype/process-sample.py \
--allelefile  input-files/combined-TE.txt \
--allelebase genotyping-base \
--bwa bwa-0.5.9/bwa \
--samplename SAMPLENAME \
--bam PATH/TO/BAM/FILE.bam


The resulting per-sample VCF file will be in genotyping-base/samples/SAMPLENAME/SAMPLENAME.vcf
and other associated information will be in the same directory.


The script process-sample-remap.py can be used to repeat only the mapping portion.

## Update November 2017:
Two important updates have been added:
The default max insert size for bwa sampe has been changed to 1000, this allows proper processing if seuqnicng libraries with expected insert sizes > 500bp

An optional additional filtering step using an exclusion file has been added (--excludefile).
This is a (1-based) set of intervals on the targeted sequences (chr1_xxx_genome and chr1_xxx_alternative, etc) for exclusion.  Read pairs that map nearly entirely within the exclusion region are removed from consideration in genotyping.  This is an important
improvement for genotyping large insertions (full length LINE or ERV proviruses) as some reads may map entirely within the element and not having any anchoring information unique to that region.

## Update May 2019:
Added new module to genotype using only split reads that cross breakpoint/TSD region. Just takes one end when there is sequence deleted also from the empty allele

Step one:
python2 setup-for-splitread.py 

requires that the program age be in the path, use a recent version of bwa, such as 0.7.15 or newer.
This step is only required to be run once.


Step two:
python process-sample-remap-splitread.py

Ran with proper files for each sample.  Uses extracted reads from previous run






 
