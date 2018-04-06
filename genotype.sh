#! /bin/bash

### dep: BGZIP, TABIX, vcf-merge
### usage: ./genotype.sh <input_from_de_novo_create_input.sh> projectname names_and_bam

BWA="~/bin/bwa-0.7.16a/bwa"

mkdir -p $2 # adjust with path of working directory

####################################
# STEP 1: Create alternate alleles #
####################################

python insertion-genotype/create-alternative-alleles.py --allelefile $2/$1 --allelebase $2 --bwa $BWA


########################################
# STEP 2: Genotype loci per individual #
########################################

### GO IN LINE MODE

IFS_BAK=$IFS
IFS=$'\n'

list=$(cat $3) ### Need list of individuals
for line in $list
do

ind=$(printf '%s\n' "$line" | awk '{print $1}')
bam=$(printf '%s\n' "$line" | awk '{print $2}')

python insertion-genotype/process-sample.py --allelefile $2/$1 --allelebase $2 --samplename $ind --bwa $BWA --bam $bam

done

### QUIT LINE MODE
IFS=$IFSq_BAK
IFS_BAK=

############################
# STEP 3: CONCATENATE VCFs #
############################

##### 3a: tabiz/gz-compress

### has to be done in */samples/ folder
folders=$(cut -f 1 $3) ## take the ind ID's, should be folders names 
for fol in $folders
do

bgzip -c $2/samples/$fol/$fol.vcf > $2/samples/$fol/$fol.vcf.gz
tabix -p vcf $2/samples/$fol/$fol.vcf.gz

##### 3b: merge vcf

vcf-merge $2/samples/*/*.vcf.gz | bgzip -c > $2.REgenotypeTE.vcf.gz # FINAL VCF (do not include wanrings about mappability and flags if corrected full length)

done