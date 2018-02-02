#! /bin/bash

##########################################
# reGenotypeTE - identify_mei_from_RM.sh #
#                                        #
# This script intersects the             #
# discordant reads of each individals    #
# with the RepeatMasker track in order   #
# to identify the TE sequence            #
#                                        #
# This script is run for the loci for    #
# which the assembly has failed          #
#                                        #
# Author: Clement Goubert                #
# Date: 02/01/2018                       #
##########################################

### USAGE: ./identify_mei_from_RM.sh loci_bam_folder output_folder #folders must be in full path (no relative) without final "/"

### Set paths (CHANGE IN FINAL VERSION TO MATCH ALL SCRIPTS)
BEDTOOLS="/home/cgoubert/bin/bedtools2/bin"
RM_TRACK="/home/cgoubert/CorrectHet/RepeatMasker_Alu_hg38_lift_hg19.bed"
RM_FASTA="/home/cgoubert/bin/RepeatMasker/Libraries/dc20170127-rb20170127/homo_sapiens/refinelib"

rm -r $2 # erase previous output if same name
mkdir $2 # creates the output folder if inexistent

# Pick up the discordant mates from the bam files, generates .bed per individual / locus
paste <(date | awk '{print $4}') <(echo "generates individual .bed file per locus...")

cd $1
listeofposition=$(ls)
for pos in $listeofposition
do
	echo "$pos..."
	mkdir -p $2/$pos # create one folder per position in the outfolder
	cd $1/$pos # switch inside the position folder in the bam folder
	listofbams=*.bam # list the individual bams. One bam / individual / locus without full TE assembly.
	for bam in $listofbams
	do
		samtools view -F 14 $bam | awk '($3!=$7 && $7!="=")' | awk '{print "chr"$7"\t"$8"\t"($8+100)}' > $2/$pos/$bam"".bed # generates individual bed from bam
	done

done
# Concatenate individual .bed per locus
paste <(date | awk '{print $4}') <(echo "concatenate individual discordant reads .bed per locus...")

cd $2 # switch to output folder
mkdir -p concatenated_bed # creates the folder that will hold the concatenated bed per position
listeofposition2=$(ls | grep -v 'bed') # make the list of position folders
for pos2 in $listeofposition2 # loop over every position to concatenate the beds
do
	echo "$pos2..."
	cd $2/$pos2
	#positions=$(ls *.bed | sed $'s/\./\t/g' | cut -f 2 | sort | uniq) # creates the list of unique positions (locus)
	#for posi in $positions
	#do
cat *$pos2*.bed | grep -v "chr\*" > $2/concatenated_bed/$pos2""_concat.bed # concatenates the individual bed per locus

done


# Intersect the locus .bed with the RM track
paste <(date | awk '{print $4}') <(echo "intersect concatenated locus .bed with RepeatMasker track...")

cd $2/concatenated_bed
filelist=*_concat.bed # list the locus bedfile in preparation of the intersection with the RM track

for file in $filelist
do
	echo "$file..."
	# intersect for each locus the RepeatMasker track and select the best overlap
	TE=$($BEDTOOLS/bedtools intersect -wao -a <(awk '{print $0"\t"NR}' $file) -b $RM_TRACK | awk '{print $1"_"$2"_"$3"_"$4"\t"$0}' | sort -k1,1 -k9,9nr | sort -u -k1,1 | sed $'s/repName=//g;s/;/\t/g' | cut -f 9 | sort | uniq -c | sort -k1,1nr | grep -v "\." | head -n 1)
	echo "$TE"
	paste <(echo "$file") <(echo $TE) | sed $'s/_concat\.bed//g' | sed 's/\.\.\///g' >> $2/position_and_TE # print the locus | number of discordant mates that support | TE (one line per locus)
	#paste <(echo "$file") <(echo $TE) | sed $'s/_concat\.bed//g' | sed 's/\.\.\///g;s/\t/xXx/g;s/ /xXx/g' >> $2/position_and_TE # print the locus | number of discordant mates that support | TE (one line per locus)
done

# Extract the TE sequence
paste <(date | awk '{print $4}') <(echo "Extracting the TE sequences in fasta format...")

mkdir $2/TE_sequences
awk '{print $3}' $2/position_and_TE | sort | uniq | awk '{print $1"#SINE/Alu"}' > $2/TE_headers

TEheads=$(cat $2/TE_headers)
for head in $TEheads
do
	name=$(echo "$head" | sed 's/\#SINE\/Alu//g')
	echo "$name"
	perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$head") $RM_FASTA > $2/TE_sequences/$name"".fasta
done

rm $2/TE_headers

paste <(date | awk '{print $4}') <(echo "Done! Results in $2")