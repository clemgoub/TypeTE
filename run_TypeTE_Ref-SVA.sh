#! /bin/bash

###########################################
# TypeTE - run_TypeTE_Del-SVA             #
#                                         #
# This is the main script of the pipeline #
# This script genotype the reference SVA  #
# insertions                              #
# Author: Clement Goubert                 #
# Date: 04/15/2019                        #
# Version: 1.0                            #
###########################################

#load the user options, outdir path and dependencies paths
source parameterfile_Ref.init # load the parameterfile as argument

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
 rm -r $OUTDIR/$PROJECT/* #cleanup previous failed run

#Creates the <project>.input in $OUTDIR/$PROJECT

paste <(date) <(echo "preparing input from MELT vcf...")

./input_from_melt_Del.sh $VCF $PROJECT

paste <(date) <(echo "DONE.")

 #########################################
 # 1:  #### Match MEI with RM insertions #
 #########################################

paste <(date) <(echo "Finding corresponding Repeat Masker insertions on reference genome...")

# if 'chr' in the vcf keep like that, otherwise, remove the 'chr' from the RM track for the bedtool intersct!

if grep "chr" $OUTDIR/$PROJECT/$PROJECT.input > /dev/null
then
    perl 01_DelP_findcorrespondinginsertion_v3.2-SVA.pl -t $RM_TRACK -f $OUTDIR/$PROJECT/$PROJECT.input -p $OUTDIR/$PROJECT/RM_intervals.out
else
    perl 01_DelP_findcorrespondinginsertion_v3.2-SVA.pl -t <(sed 's/chr//g' $RM_TRACK) -f $OUTDIR/$PROJECT/$PROJECT.input -p $OUTDIR/$PROJECT/RM_intervals.out
fi

#########################################
# 2:  #### Find TSD and TE coordinates  #
#########################################

paste <(date) <(echo "Finding TSD and TE coordinates...")

perl 03_DelP_findTSD_forRMTEcordinates_v3.4.pl -t $OUTDIR/$PROJECT/RM_intervals.out/file.correspondingRepeatMaskerTEs.txt -p $OUTDIR/$PROJECT/output_TSD_Intervals.out -g $GENOME

#########################################
# 3:  #### Create input for genotyping  #
#########################################

paste <(date) <(echo "generating inputs for genotyping...")

paste <(sort -k1,1 $OUTDIR/$PROJECT/output_TSD_Intervals.out/TEcordinates_with_bothtsd_cordinates.v.3.4.txt) <(sort -k1,1 $OUTDIR/$PROJECT/RM_intervals.out/file.correspondingRepeatMaskerTEs.txt) | cut -f 1,2,3,4,11 > $OUTDIR/$PROJECT/RM_insertions_TSD_strands

./deletion_create_input.sh $OUTDIR/$PROJECT/RM_insertions_TSD_strands > $OUTDIR/$PROJECT/$PROJECT.allele

#####################
# 4: re-Genotype ####
#####################

paste <(date | awk '{print $4}') <(echo "Genotyping...")

# remove older files in case of re-run
rm -r $OUTDIR/$PROJECT/locusAlleles
rm -r $OUTDIR/$PROJECT/samples

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


