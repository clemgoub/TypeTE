#! /bin/bash

###########################################
# TypeTE - run_TypeTE_NoRef-SVA-L1.sh     #
#                                         #
# This is the main script of the pipeline #
# to genotype Non-Reference insertions    #
# SVA and L1 based en MELT2.1.4 outputs   #
#                                         #
# Author: Clement Goubert                 #
# Date: 04/16/2019                        #
# Version: 1.0                            #
###########################################

# Usage: ./run_TypeTE_NoRef-SVA-L1.sh <SVA/L1>

#load the user options, outdir path and dependencies paths
source parameterfile_NoRef.init

#START
echo "#####################"
echo "#      TypeTE       #"
echo "#####################"

################################################
# 0: Setup and create input from MELT files ####
################################################

#locate working directoty
whereamI=$(pwd)

#Creates the $OUTDIR

{
if [ $OUTDIR == "" ]; then
    echo "OUTDIR not set, don't know where to create output folder..."
    exit 0
fi
}

mkdir -p $OUTDIR/$PROJECT

#Creates the fasta file with sequences of each insertion using Julie's script

testfile=$(file $VCF | awk '{print $2}')
if [[ $testfile == "gzip" ]]
	then
	python2.7 ./MEI-VCF-to-FASTA/Convert_MELT_vcf_to_fasta.py <(zcat $VCF) $1 > $OUTDIR/$PROJECT/$1"".consensi.fa
	else
	python2.7 ./MEI-VCF-to-FASTA/Convert_MELT_vcf_to_fasta.py $VCF $1 > $OUTDIR/$PROJECT/$1"".consensi.fa
fi

#Creates the <project>.input in $OUTDIR/$PROJECT

paste <(date | awk '{print $4}') <(echo "preparing input from MELT vcf...")

./input_from_melt.sh $VCF $PROJECT

join -11 -21 <(sort -k1,1 $OUTDIR/$PROJECT/$PROJECT.input) <(grep '>' $OUTDIR/$PROJECT/$1"".consensi.fa | sed 's/_MEINFO/\tMEINFO/g;s/>//g' | awk '{print $1"\t"$1"_"$2}' | sort -k1,1) > $OUTDIR/$PROJECT/$PROJECT.genoinput

paste <(date | awk '{print $4}') <(echo "DONE.")

#######################################
# 6: Generate input for genotyping ####
#######################################

paste <(date | awk '{print $4}') <(echo "Generating input table for genotyping...")

rm $OUTDIR/$PROJECT/$PROJECT.allele

./de_novo_create_input_SVA-L1.sh $OUTDIR/$PROJECT/$PROJECT.genoinput $OUTDIR/$PROJECT/$PROJECT.allele

paste <(date | awk '{print $4}') <(echo "Generating input table for genotyping...Done.")

#####################
# 7: re-Genotype ####
#####################

paste <(date | awk '{print $4}') <(echo "Genotyping...")

### create alternatives alleles
python2.7 insertion-genotype/create-alternative-alleles.py --allelefile $OUTDIR/$PROJECT/$PROJECT.allele --allelebase $OUTDIR/$PROJECT --bwa $BWA

### genotype per individual
cat $BAMFILE | $PARALLEL -j $CPU --colsep '\t' --results $OUTDIR/$PROJECT/genotyping_logs "python2.7 insertion-genotype/process-sample.py --allelefile $OUTDIR/$PROJECT/$PROJECT.allele --allelebase $OUTDIR/$PROJECT --samplename {1} --bwa $BWA --bam $BAMPATH/{2}"


paste <(date | awk '{print $4}') <(echo "Genotyping... Done")
paste <(date | awk '{print $4}') <(echo "Generating outputs...")

### tabix the individuals vcfs
cd $OUTDIR/$PROJECT/samples/
folders=$(ls -d */ | sed 's/\///g')
for fol in $folders
do

$BGZIP -c $fol/$fol.vcf > $fol/$fol.vcf.gz
$TABIX -p vcf $fol/$fol.vcf.gz

done

### merging final VCF
vcf-merge ./*/*.vcf.gz | $BGZIP -c > $OUTDIR/$PROJECT/$PROJECT.TypeTE.vcf.gz

cd $whereamI

paste <(date | awk '{print $4}') <(echo "Generating outputs... Done")

paste <(date | awk '{print $4}') <(echo "TypeTE completed!!!")


