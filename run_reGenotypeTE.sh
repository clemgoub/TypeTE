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
#Creates the $OUTDIR

{
if [ $OUTDIR == "" ]; then
    echo "OUTDIR not set, don't know where to create output folder..."
    exit 0
fi
}

mkdir -p $OUTDIR/$PROJECT

#Creates the <project>.input in $OUTDIR/$PROJECT

./input_from_melt.sh $VCF $PROJECT


#Join ind names to coordinates and generate the list of locus/individuals ("$OUTDIR/$PROJECT/sample_file.txt.list.txt")

perl makelist_v1.0.pl -t $BAMFILE -f $OUTDIR/$PROJECT/$PROJECT.input -p $OUTDIR/$PROJECT

# split the input in order to paralellize read extraction

perl 02_splitfile_jt_v3.0_pipeline.pl -f $OUTDIR/$PROJECT/sample_file.txt.list.txt -s yes -n $individual_nb

#splitfile.pl -out $OUTDIR uses file.list.text based on -individuals --> splitbyindividuals folder created with one file per individual_nb + list_of_splitted_files_names in $OUTDIR

#cd  $OUTDIR/$PROJECT/splitbyindividuals

# $SCRIPT_PATH/process_bam........pl $BAMFILE $BAMPATH -p $OUTDIR

# orienTE_extracTE.pl -d $OUTDIR/processbamout -t TE_directory_from_indetify_mei_from_RM.sh -l list_outout_from_indetify_mei_from_RM.sh -g ExtractGenomicSequences

# de_novo_create_input.sh (correct it, and modify Jainy's table before)