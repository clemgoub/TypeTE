#! /bin/bash

####################################################################
# 4: Find TE annotations and consensus using RepeatMasker track ####
####################################################################

source parameterfile_tri.init
whereamI=$(pwd)

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
rm -r $OUTDIR/$PROJECT/Repbase_intersect/position_and_TE # remove the output table for safety if one already exist
ls $OUTDIR/$PROJECT/splitbylocus/*.part | $PARALLEL -j $CPU --results $OUTDIR/$PROJECT/Repbase_intersect "./identify_mei_from_RM_tri.sh {} $OUTDIR/$PROJECT/Repbase_intersect"


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
rm $OUTDIR/$PROJECT/Repbase_intersect/TE_headers

paste <(date | awk '{print $4}') <(echo "Done! Results in $2")

