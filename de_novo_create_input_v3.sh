#! /bin/bash

##########################################
# reGenotypeTE - de_novo_create_input.sh #
#                                        #
# This script creates the input file     #
# of the de-novo genotyping script       #
#                                        #
# Author: Clement Goubert                #
# Date: 03/2018                          #
##########################################

# usage: ./de_novo_create_input.sh <input>.input <orientTE_extractTE.pl>.outputtable <indentify_mei_from_RM.sh>.outfolder <Assembled.TE> <output.file>

# SET PATH FOR FILES AND PROGRAMS
source ./parameterfile_NoRef.init
rm -r $5 &> /dev/null

## Add line to process Jainy's output
join -a1 -11 -21 <(sort -k1,1 <(join -11 -21 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.2,2.3,2.4,2.5 <( paste <(sed 's/ /\t/g' $1 | cut -f 1 | sed 's/_/\t/g') <(sed 's/ /\t/g' $1| cut -f 2-) | awk '{print $1":"($2-250)"-"($2+250),$3,$4,$5,$6,$7}' | sort -k1,1 ) <(sort -k1,1 $2))) <(sort -k1,1 $3/position_and_TE) \
| awk '{if (NF == 10) {print $1,$2,$3,$4,$5,$6,$7,$8,$10"#SINE/Alu",$6} else {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' > $OUTDIR/$PROJECT/pre_input

input=$(cat $OUTDIR/$PROJECT/pre_input) # change it so $1 is the output of previous line

#input=$(join -a1 -11 -21 <(sort -k1,1 <(join -11 -21 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.2,2.3,2.4,2.5 <(sort -k1,1 $1) <(sort -k1,1 $2))) <(sort -k1,1 $3/position_and_TE) \
#| awk '{if (NF == 10) {print $1,$2,$3,$4,$5,$6,$7,$8,$10"#SINE/Alu",$6} else {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}')


# change delimiter (IFS) to new line.
IFS_BAK=$IFS
IFS=$'\n'

for list in $input
do


pos=$(printf '%s\n' "$list" | awk '{print $1}')
assembled=$(printf '%s\n' "$list" | awk '{print $7}')
TSD=$(printf '%s\n' "$list" | awk '{print $NF}' | sed 's/null//g')
direction=$(printf '%s\n' "$list" | awk '{if ($8 == "?") {print $5} else {print $8}}') # if not assembled use MELT orientation and TSD
header=$(printf '%s\n' "$list" | awk '{print $9}')

# create bed files -500bp ---- BP (left)  BP ---- +500bp (right)
sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print $1"\t"($2-250)"\t"($3+250)"\t"$4"\t"$5"-"$6}' > $OUTDIR/$PROJECT/region.bed

echo ""
echo "$pos"
echo "$assembled"
echo "$TSD"
echo "$direction"
echo "$header"

$BEDTOOLS getfasta -fi $GENOME -bed $OUTDIR/$PROJECT/region.bed > $OUTDIR/$PROJECT/region''$pos''.fasta
makeblastdb -in $OUTDIR/$PROJECT/region''$pos''.fasta -out $OUTDIR/$PROJECT/region''$pos''.fasta -dbtype 'nucl'
blastn -query <(echo "$TSD") -db $OUTDIR/$PROJECT/region''$pos''.fasta -word_size 6 -outfmt 6 | sort -k12,12nr -k1,1 | head -n 1 | awk '{if ($9 < $10) {print $2"\t"$9"\t"$10} else {print $2"\t"$10"\t"$9}}' > $OUTDIR/$PROJECT/blast''$pos''.TSD.bed

# Check if we found the expected TSD
if [[ -s "$OUTDIR/$PROJECT/blast''$pos''.TSD.bed" ]]
then # left = 0-TSD_start right = TSD_end-1000
	awk '{print $1"\t1\t500"}' $OUTDIR/$PROJECT/blast''$pos''.TSD.bed > $OUTDIR/$PROJECT/left.bed
	awk '{print $1"\t501\t1000"}' $OUTDIR/$PROJECT/blast''$pos''.TSD.bed > $OUTDIR/$PROJECT/right.bed
else # cut at the middle!
	awk '{print $1"\t1\t"$2}' $OUTDIR/$PROJECT/blast''$pos''.TSD.bed > $OUTDIR/$PROJECT/left.bed
	awk '{print $1"\t"$3"\t1000"}' $OUTDIR/$PROJECT/blast''$pos''.TSD.bed > $OUTDIR/$PROJECT/right.bed
fi

# Check if Alu is fully assembled
if [[ $assembled == "yes" ]] ### For assembled Alu
then #fully assembled TE have TSDs in 5' and 3' within the assembled sequence

#extract the 5' flank up to TSD
$BEDTOOLS getfasta -fi $OUTDIR/$PROJECT/region''$pos''.fasta -bed $OUTDIR/$PROJECT/left.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $OUTDIR/$PROJECT/$pos"".left.seq
#extract the 3' flank after TSD
$BEDTOOLS getfasta -fi $OUTDIR/$PROJECT/region''$pos''.fasta -bed $OUTDIR/$PROJECT/right.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $OUTDIR/$PROJECT/$pos"".right.seq

 	if [[ $direction == "-" ]] # if the TE is reverse...
 	then
 		# ...get assemble TE, Reverse and Complement 
		perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$header") $4 | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2 | rev | tr "ATGCatgc" "TACGtacg" > $OUTDIR/$PROJECT/$pos"".TE.seq
 	else
 		# ...esle if + paste TE as assembled
		perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$header") $4 | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' |  cut -f 2 | sed 's/,//g;s/null//g' > $OUTDIR/$PROJECT/$pos"".TE.seq
	fi
	# write the .allele file line for this locus
	paste <(echo "$pos") <(sed 's/:/\t/g' <(echo "$pos") | cut -f 1) <(sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print ($3-250)}') <(echo ".") <(sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print ($2-250)"\t"($3+250)}') <(awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  $OUTDIR/$PROJECT/region''$pos''.fasta | cut -f 2)  <(paste -d, $OUTDIR/$PROJECT/$pos"".left.seq $OUTDIR/$PROJECT/$pos"".TE.seq $OUTDIR/$PROJECT/$pos"".right.seq | sed 's/,//g') >> $5

else ### For non-assembled
#The TE taken from repbase don't have TSD, so the TSD is to add.

#extract the 5' flank up to TSD
$BEDTOOLS getfasta -fi $OUTDIR/$PROJECT/region''$pos''.fasta -bed $OUTDIR/$PROJECT/left.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $OUTDIR/$PROJECT/$pos"".left.seq
#extract the 3' flank after TSD
$BEDTOOLS getfasta -fi $OUTDIR/$PROJECT/region''$pos''.fasta -bed $OUTDIR/$PROJECT/right.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $OUTDIR/$PROJECT/$pos"".right.seq

 	if [[ $direction == "-" ]]
 	then
 		# ...get reference TE, Reverse and Complement  
		perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$header") $RM_FASTA | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2 | rev | tr "ATGCatgc" "TACGtacg" > $OUTDIR/$PROJECT/$pos"".TE.seq
 	else
 		# ...esle if + get reference TE as in repbase (+)
		perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$header") $RM_FASTA | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' |  cut -f 2 | sed 's/,//g;s/null//g' > $OUTDIR/$PROJECT/$pos"".TE.seq
	fi
	# write the .allele file line for this locus
 	paste <(echo "$pos") <(sed 's/:/\t/g' <(echo "$pos") | cut -f 1) <(sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print ($3-250)}') <(echo ".") <(sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print ($2-250)"\t"($3+250)}') <(awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  $OUTDIR/$PROJECT/region''$pos''.fasta | cut -f 2)  <(paste -d, $OUTDIR/$PROJECT/$pos"".left.seq <(echo "$TSD") $OUTDIR/$PROJECT/$pos"".TE.seq <(echo "$TSD") $OUTDIR/$PROJECT/$pos"".right.seq | sed 's/,//g') >> $5

fi

# ### cleaning
# rm left.bed
# rm right.bed
# rm *.seq

done

rm $OUTDIR/$PROJECT/pre_input
rm $OUTDIR/$PROJECT/right.bed
rm $OUTDIR/$PROJECT/*.seq
rm $OUTDIR/$PROJECT/region*
rm $OUTDIR/$PROJECT/blast*.bed

# return delimiter to previous value
IFS=$IFSq_BAK
IFS_BAK=
