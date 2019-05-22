#! /bin/bash

# TEST-SCRIPT only performing genotyping
# Author Clement Goubert
#
# Usage: ./Genotype.sh <allelefile> <allelebase> <parameterfile>

source $3

paste <(date | awk '{print $4}') <(echo "Genotyping...")

### create alternatives alleles
python2.7 insertion-genotype/create-alternative-alleles.py --allelefile $1 --allelebase $2 --bwa $BWA

### genotype per individual
cat $BAMFILE | $PARALLEL -j $CPU --colsep '\t' --results $2/genotyping_logs "python2.7 insertion-genotype/process-sample.py --allelefile $1 --allelebase $2 --samplename {1} --bwa $BWA --bam $BAMPATH/{2}"

paste <(date | awk '{print $4}') <(echo "Genotyping... Done")
paste <(date | awk '{print $4}') <(echo "Generating outputs...")

### tabix the individuals vcfs
cd $2/samples/
folders=$(ls -d */ | sed 's/\///g')
for fol in $folders
do

$BGZIP -c $fol/$fol.vcf > $fol/$fol.vcf.gz
$TABIX -p vcf $fol/$fol.vcf.gz

done

### merging final VCF
vcf-merge ./*/*.vcf.gz | $BGZIP -c > $2/out.TypeTE.vcf.gz
