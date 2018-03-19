#! /bin/bash

##########################################
# reGenotypeTE - intersect_RM.sh         #
#                                        #
# This sub-script intersects the         #
# discordant reads of each individals    #
# with the RepeatMasker track in order   #
# to identify the TE sequence            #
#                                        #
# Run for every sub list  of its main    #
# scrit "identify_mei_from_RM.sh"        #
#                                        #
# Author: Clement Goubert                #
# Date: 02/07/2018                       #
##########################################

### Set paths (CHANGE IN FINAL VERSION TO MATCH ALL SCRIPTS) 
BEDTOOLS="/home/cgoubert/bin/bedtools2/bin"
RM_TRACK="/home/cgoubert/CorrectHet/RepeatMasker_Alu_hg38_lift_hg19_corrected.bed"
RM_FASTA="/home/cgoubert/bin/RepeatMasker/Libraries/dc20170127-rb20170127/homo_sapiens/refinelib"

listeofposition=$(cat $1)
subfile=$(echo "$1")
mkdir -p $2/fol_$subfile

for pos in $listeofposition
do
	mkdir -p $2/fol_$subfile/$pos # create one folder per position in the outfolder
	cd $3/$pos # switch inside the position folder in the bam folder
	listofbams=*.bam # list the individual bams. One bam / individual / locus without full TE assembly.
	for bam in $listofbams
	do
		samtools view -F 14 $bam | awk '($3!=$7 && $7!="=")' | awk '{print "chr"$7"\t"$8"\t"($8+100)}' > $2/fol_$subfile/$pos/$bam"".bed # generates individual bed from bam
	done

done
# Concatenate individual .bed per locus
#paste <(date | awk '{print $4}') <(echo "concatenate individual discordant reads .bed per locus...")

cd $2/fol_$subfile # switch to output folder
mkdir -p concatenated_bed # creates the folder that will hold the concatenated bed per position
listeofposition2=$(ls | grep -v 'bed') # make the list of position folders
for pos2 in $listeofposition2 # loop over every position to concatenate the beds
do
	echo "$pos2..."
	cd $2/fol_$subfile/$pos2
	#positions=$(ls *.bed | sed $'s/\./\t/g' | cut -f 2 | sort | uniq) # creates the list of unique positions (locus)
	#for posi in $positions
	#do
cat *$pos2*.bed | grep -v "chr\*" > $2/fol_$subfile/concatenated_bed/$pos2""_concat.bed # concatenates the individual bed per locus

done


# Intersect the locus .bed with the RM track
#paste <(date | awk '{print $4}') <(echo "intersect concatenated locus .bed with RepeatMasker track...")

cd $2/fol_$subfile/concatenated_bed
filelist=*_concat.bed # list the locus bedfile in preparation of the intersection with the RM track

for file in $filelist
do
	echo "$file..."
	# intersect for each locus the RepeatMasker track and select the best overlap
	TE=$($BEDTOOLS/bedtools intersect -wao -a <(awk '{print $0"\t"NR}' $file) -b $RM_TRACK | awk '{print $1"_"$2"_"$3"_"$4"\t"$0}' | sort -k1,1 -k9,9nr | sort -u -k1,1 | sed $'s/repName=//g;s/;/\t/g' | cut -f 9 | sort | uniq -c | sort -k1,1nr | grep -v "\." | head -n 1)
	echo "$TE"
	paste <(echo "$file") <(echo $TE) | sed $'s/_concat\.bed//g' | sed 's/\.\.\///g' >> $2/fol_$subfile/position_and_TE # print the locus | number of discordant mates that support | TE (one line per locus)
	#paste <(echo "$file") <(echo $TE) | sed $'s/_concat\.bed//g' | sed 's/\.\.\///g;s/\t/xXx/g;s/ /xXx/g' >> $2/position_and_TE # print the locus | number of discordant mates that support | TE (one line per locus)
done