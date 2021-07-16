#! /bin/bash

###########################################
# TypeTE - run_TypeTE_NRef                #
#                                         #
# This is the main script of the pipeline #
# to genotype Non-Reference insertions    #
#                                         #
# Author: Clement Goubert                 #
# Date: 07/2021                           #
# Version: 1.1                            #
###########################################

# ***changelog v1.0 --> v1.1***
# - Output vcf:
#    - Cleanup output vcfs from irrelevant info fields in header
#    - Reference alleles are now printed in the traditionnal (REF/ALT) format, with REF = TE present = 0, and ALT = TE absent (deletion) = 1. 
# - Hard code python2.7 in assembly script to match Spades requirements
# - Improve Non-Reference allele reconstruction script at TSD
# - Clean bugs and silence non-threatening error messages
# - Change parameterfile_NoRef.ini to parameterfile_NRef.ini to match regular script naming
# - Create tutorial section (upcoming manuscript)

#load the user options, outdir path and dependencies paths
source ./parameterfile_NRef.init

#START
echo "###########################"
echo "# TypeTE - Non-Reference  #"
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

./input_from_melt.sh $VCF $PROJECT

paste <(date | awk '{print $4}') <(echo "DONE.")

#################################################
# 1: Joining individuals and MEI coordinates ####
#################################################

echo "Joining individuals and MEI coordinates..."
perl makelist_v1.0.pl -t $BAMFILE -f $OUTDIR/$PROJECT/$PROJECT.input -p $OUTDIR/$PROJECT

#####################################
# 2: Split input per individuals ####
#####################################

paste <(date | awk '{print $4}') <(echo "DONE.")
paste <(date | awk '{print $4}') <(echo "Splitting individuals for paralellization of read extraction...")

perl 02_splitfile_jt_v3.0_pipeline.pl -f $OUTDIR/$PROJECT/file.list.txt -s yes -n $individual_nb -p $OUTDIR/$PROJECT

#######################
# 3: Extract reads ####
#######################

paste <(date | awk '{print $4}') <(echo "Extracting reads...")

cd  $OUTDIR/$PROJECT/splitbyindividuals #cd in the splitfile directory

cat ../List_of_split_files.txt | $PARALLEL -j $CPU --results $OUTDIR/$PROJECT/Process_bams "perl $whereamI/03_processbam_forreadextract_v15.0.pl -g $GENOME -t $BAMFILE -f {} -p $OUTDIR/$PROJECT -bl $BAMPATH -pt $PICARD -sq $SEQTK -bu $BAMUTILS -bt $BEDTOOLS" &> $OUTDIR/$PROJECT/extractreads.log

cd $whereamI #comes back to the working dir

# paste <(date | awk '{print $4}') <(echo "Extracting mappability scores...")

# perl denovo_extract_GM_scoresv1.0.pl -t hg19wgEncodeCrgMapabilityAlign100mer_index -f $OUTDIR/$PROJECT/$PROJECT.input -p $OUTDIR/$PROJECT/gmscore_all -db $DATABASE -u $USER -pd $PSWD


####################################################################
# 4: Find TE annotations and consensus using RepeatMasker track ####
####################################################################

paste <(date | awk '{print $4}') <(echo "Finding Repbase consensus for each MEI...")

# counting loci and dividing into subfiles for paralellization

total_locus=$(ls -lh $OUTDIR/$PROJECT/IGV | awk ' NR > 1 {print $NF}')
nb_locus=$(echo "total_locus" | wc -l)
per_file=$( echo "$nb_locus/$CPU" | bc)
if(($per_file < 1))
then 
	splitnb=$((1))
else
	splitnb=$(((nb_locus+1)/$CPU))
fi

mkdir -p $OUTDIR/$PROJECT/splitbylocus
split -l $splitnb <(echo "$total_locus") --additional-suffix .part $OUTDIR/$PROJECT/splitbylocus/locus_
ls $OUTDIR/$PROJECT/splitbylocus/locus_* > $OUTDIR/$PROJECT/splitbylocus/List_of_loci_files.txt

# Run in parallel find TE from Repbase

mkdir -p  $OUTDIR/$PROJECT/Repbase_intersect # creates the output folder if inexistent
rm -r $OUTDIR/$PROJECT/Repbase_intersect/position_and_TE &>/dev/null # remove the output table for safety if one already exist
ls $OUTDIR/$PROJECT/splitbylocus/*.part | $PARALLEL -j $CPU --results $OUTDIR/$PROJECT/Repbase_intersect "./identify_mei_from_RM.sh {} $OUTDIR/$PROJECT/Repbase_intersect"

# Extract the TE sequence
paste <(date | awk '{print $4}') <(echo "Extracting the TE sequences in fasta format...")

mkdir -p $OUTDIR/$PROJECT/Repbase_intersect/TE_sequences
awk '{print $3}' $OUTDIR/$PROJECT/Repbase_intersect/position_and_TE | sort | uniq | awk '{print $1"#SINE/Alu"}' > $OUTDIR/$PROJECT/Repbase_intersect/TE_headers

TEheads=$(cat $OUTDIR/$PROJECT/Repbase_intersect/TE_headers)
for head in $TEheads
do
	name=$(echo "$head" | sed 's/\#SINE\/Alu//g')
	echo "$name"
	perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$head") $RM_FASTA > $OUTDIR/$PROJECT/Repbase_intersect/TE_sequences/$name"".fasta
done

cat $OUTDIR/$PROJECT/Repbase_intersect/TE_sequences/*.fasta > $OUTDIR/$PROJECT/Repbase_intersect/TE_sequences/RM_consensus.fa
rm $OUTDIR/$PROJECT/Repbase_intersect/TE_headers &>/dev/null

paste <(date | awk '{print $4}') <(echo "Done! Results in $2")

############################################################
# 5: de novo Assembly, orientation and find TSDs of MEI ####
############################################################

paste <(date | awk '{print $4}') <(echo "Assembling MEI, retreiving orientation and TSDs...")

rm -r $OUTDIR/$PROJECT/Assembled_TEreads &>/dev/null
perl 04_orientTE_extractTE_v5.0_pipeline.pl -p $OUTDIR/$PROJECT -d $OUTDIR/$PROJECT/orientTE -g $OUTDIR/$PROJECT/ExtractGenomicsequences -t $OUTDIR/$PROJECT/Repbase_intersect/TE_sequences -l $OUTDIR/$PROJECT/Repbase_intersect/position_and_TE -sp $SPADES -mn $MINIA -cp $CAP3 -bp $BLAST -cu $CPU &> $OUTDIR/$PROJECT/assembly.log

#######################################
# 6: Generate input for genotyping ####
#######################################

paste <(date | awk '{print $4}') <(echo "Generating input table for genotyping...")

rm $OUTDIR/$PROJECT/$PROJECT.allele &>/dev/null
FullLength="$OUTDIR/$PROJECT/Assembled_TEsequences.*.txt.full_len.fasta"
./de_novo_create_input_v3.sh $OUTDIR/$PROJECT/$PROJECT.input $OUTDIR/$PROJECT/genomeloc.strand.prediction.5.1.txt $OUTDIR/$PROJECT/Repbase_intersect $FullLength $OUTDIR/$PROJECT/$PROJECT.allele

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


