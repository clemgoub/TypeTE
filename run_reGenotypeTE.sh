#! /bin/bash

###########################################
# reGenotypeTE - run_reGonotypeTE         #
#                                         #
# This is the main script of the pipeline #
#                                         #
# Author: Clement Goubert                 #
# Date: 04/2018                           #
# Version: 1.0                            #
###########################################

#load the user options, outdir path and dependencies paths
source parameterfile.init

#START
echo "###########################"
echo "#      reGenotypeTE       #"
echo "###########################"

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

./input_from_melt.sh $VCF $PROJECT


#Join ind names to coordinates and generate the list of locus/individuals ("$OUTDIR/$PROJECT/sample_file.txt.list.txt")

paste <(date | awk '{print $4}') <(echo "DONE.")
echo "Joining individuals and MEI coordinates..."

perl makelist_v1.0.pl -t $BAMFILE -f $OUTDIR/$PROJECT/$PROJECT.input -p $OUTDIR/$PROJECT

# split the input in order to paralellize read extraction

paste <(date | awk '{print $4}') <(echo "DONE.")
echo "Splitting individuals for paralellization of read extraction..."

perl 02_splitfile_jt_v3.0_pipeline.pl -f $OUTDIR/$PROJECT/file.list.txt -s yes -n $individual_nb -p $OUTDIR/$PROJECT

# Process bams: extract reads from bam files and extract mappability

#cd in splitfile directory
cd  $OUTDIR/$PROJECT/splitbyindividuals

#perl 03_processbam_extract_GM_scoresv15.0.pl \
#-t $BAMFILE \ #bam id list
#-f MEI_1KGP_data_10av \
#-g /vbod2/jainy/SGDP/Project2/hg19.refFIX.fa \
#-bl /vbod2/cgoubert/Correct_Genotypes/1KGP_bams \
#-p version15.0 -pt /home/jainy/software/picard-2.9.2 \
#-m $MAP -db jainys_db -u jainy -pd wysql123 \
#-mt hg19wgEncodeCrgMapabilityAlign100mer_index


#comes back to working dir
cd $whereamI

# $SCRIPT_PATH/process_bam........pl $BAMFILE $BAMPATH -p $OUTDIR

# orienTE_extracTE.pl -d $OUTDIR/processbamout -t TE_directory_from_indetify_mei_from_RM.sh -l list_outout_from_indetify_mei_from_RM.sh -g ExtractGenomicSequences

# de_novo_create_input.sh (correct it, and modify Jainy's table before)
