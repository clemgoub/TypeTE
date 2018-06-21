#! /bin/bash

###########################################
# reGenotypeTE - run_reGonotypeTE_Del     #
#                                         #
# This is the main script of the pipeline #
# This script genotype the non-reference  #
# insertions                              #
# Author: Clement Goubert                 #
# Date: 06/2018                           #
# Version: 1.0                            #
###########################################

#load the user options, outdir path and dependencies paths
source parameterfile_Del.init # load the parameterfile as argument

#START
echo "###########################"
echo "#      reGenotypeTE       #"
echo "###########################"

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

#Creates the <project>.input in $OUTDIR/$PROJECT

paste <(date | awk '{print $4}') <(echo "preparing input from MELT vcf...")

./input_from_melt_Del.sh $VCF $PROJECT

paste <(date | awk '{print $4}') <(echo "DONE.")

#########################################
# 1:  #### Match MEI with RM insertions #
#########################################

paste <(date | awk '{print $4}') <(echo "Finding corresponding Repeat Masker insertions on reference genome...")

perl 01_DelP_findcorrespondinginsertion_v3.0.pl -t <(sed 's/chr//g' $RM_TRACK) -f $OUTDIR/$PROJECT/$PROJECT.input -p $OUTDIR/$PROJECT/RM_intervals.out

#########################################
# 2:  #### Find Mappability Intervals   #
#########################################

paste <(date | awk '{print $4}') <(echo "Computing Mappability score of MEI...")

perl 02_DelP_findmappabilityscores_genomicintervalsv2.0.pl -t hg19wgEncodeCrgMapabilityAlign100mer_index -f $OUTDIR/$PROJECT/RM_intervals.out/file.correspondingRepeatMaskerTEs.txt -db jainys_db -u jainy -pd wysql123 -p $OUTDIR/$PROJECT


#####################
# 7: re-Genotype ####
#####################

# paste <(date | awk '{print $4}') <(echo "Genotyping...")

# ### create alternatives alleles
# python insertion-genotype/create-alternative-alleles.py --allelefile $OUTDIR/$PROJECT/$PROJECT.allele --allelebase $OUTDIR/$PROJECT --bwa $BWA

# ### genotype per individual
# list=$(cut -f 1 $BAMFILE)
# for ind in $list
# do

# bamf=$(grep "$ind" $BAMFILE | cut -f 2)
# python insertion-genotype/process-sample.py --allelefile $OUTDIR/$PROJECT/$PROJECT.allele --allelebase $OUTDIR/$PROJECT --samplename $ind --bwa $BWA --bam $BAMPATH/$bamf

# done

# paste <(date | awk '{print $4}') <(echo "Genotyping... Done")
# paste <(date | awk '{print $4}') <(echo "Generating outputs...")

# ### tabix the individuals vcfs
# cd $OUTDIR/$PROJECT/samples/
# folders=$(ls -d */ | sed 's/\///g')
# for fol in $folders
# do

# $BGZIP -c $fol/$fol.vcf > $fol/$fol.vcf.gz
# $TABIX -p vcf $fol/$fol.vcf.gz

# done

# ### merging final VCF
# vcf-merge ./*/*.vcf.gz | $BGZIP -c > $OUTDIR/$PROJECT/$PROJECT.reGenotypeTE.vcf.gz

# cd $whereamI

# paste <(date | awk '{print $4}') <(echo "Generating outputs... Done")

# paste <(date | awk '{print $4}') <(echo "reGenotypeTE completed!!!")

