# reGenotypeTE
This pipeline is currenly under construction by Jainy Thomas (University of Utah) and Clement Goubert (Cornell University).
Elaborated with the collaboration of Jeffrey M. Kidd (University of Michigan)

## Purpose

reGenotypeTE is a pipeline dedicated to re-genotype Mobile Element Insertion (MEI) previously scored with MELT (Mobile Element Locator Tool, ). reGenotypeTE extracts reads from previously detected polymorphic MEI and performs using populations data a local re-assembly of both presence and absence loci. Eventually, remapping of the reads at the infividual level allow to score the genotype of the MEI using a modified version of Li's et al. genotype likelihood. This method drammatically improves the quality of the genotypes of reported MEI and can be directly used after a MELT run on both de-novo (non-reference) and deletion (reference) insertions.

## Installation

### Dependencies

These programs need to be installed and their path reported in the file "parameterfile.init"
Perl executable must be in the user path

* PERL https://www.perl.org/
* PARALLEL https://www.gnu.org/software/parallel/
* PICARD https://broadinstitute.github.io/picard/
* BEDTOOLS http://bedtools.readthedocs.io/en/latest/
* SEQTK https://github.com/lh3/seqtk
* BAMUTILS https://genome.sph.umich.edu/wiki/BamUtil
* SPADES http://cab.spbu.ru/software/spades/
* MINIA http://minia.genouest.org/
* CAP3 http://seq.cs.iastate.edu/cap3.html
* BLAST ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
* BWA http://bio-bwa.sourceforge.net/bwa.shtml
* BGZIP http://www.htslib.org/doc/bgzip.html
* TABIX http://www.htslib.org/doc/tabix.html

### Download and install

Clone from git repository:

```sh
git clone https://github.com/clemgoub/reGenotypeTE
cd reGenotypeTE
```

Then edit the file "parameterfile.init" following the indications within

## Running reGenotypeTE

### Input

### De-novo insertion (non-reference)

### Deletions (reference-insertions)
