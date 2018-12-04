#! /bin/bash

python create-alternative-alleles.py --allelefile automatic_input_test_1/input_test_2 --allelebase automatic_input_test_1 --bwa ~/bin/bwa-0.7.16a/bwa

list=$(cut -f 1 CEU_to_genotype) # need the list of the individuals to genotype
for ind in $list
do

python process-sample.py --allelefile deletion_genotypes_3/new_input_one_TSD --allelebase deletion_genotypes_3 --samplename $ind --bwa ~/bin/bwa-0.7.16a/bwa --bam /vbod2/cgoub$

done