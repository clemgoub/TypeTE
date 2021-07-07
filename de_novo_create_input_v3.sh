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

## Add line to process Jainy's output
join -a1 -11 -21 <(sort -k1,1 <(join -11 -21 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.2,2.3,2.4,2.5 <( paste <(sed 's/ /\t/g' $1 | cut -f 1 | sed 's/_/\t/g') <(sed 's/ /\t/g' $1| cut -f 2-) | awk '{print $1":"($2-250)"-"($2+250),$3,$4,$5,$6,$7}' | sort -k1,1 ) <(sort -k1,1 $2))) <(sort -k1,1 $3/position_and_TE) \
| awk '{if (NF == 10) {print $1,$2,$3,$4,$5,$6,$7,$8,$10"#SINE/Alu",$6} else {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' > pre_input

input=$(cat pre_input) # change it so $1 is the output of previous line

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
sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print $1"\t"($2-250)"\t"($3+250)"\t"$4"\t"$5"-"$6}' > $OUTPUT/$PROJECT/region.bed
# sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print $1"\t"($2-250)"\t"($3-250)"\t"$4"\t"$5"-"$6}'  > left.bed # -250 from left boundary that is already -250 = BP-250bp; -250 from right boundary = BP
# sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print $1"\t"($2+250)"\t"($3+250)"\t"$4"\t"$5"-"$6}'  > right.bed # +250 from right boundary that is already +250 = BP+250bp; +250 from left boundary = BP

echo "$pos"
echo "$assembled"
echo "$TSD"
echo "$direction"
echo "$header"

$BEDTOOLS/bedtools getfasta -fi $GENOME -bed $OUTPUT/$PROJECT/region.bed > $OUTPUT/$PROJECT/region''$pos''.fasta
makeblastdb -in $OUTPUT/$PROJECT/region''$pos''.fasta -out $OUTPUT/$PROJECT/region''$pos''.fasta -dbtype 'nucl'
blastn -query <(echo "$TSD") -db $OUTPUT/$PROJECT/region''$pos''.fasta -word_size 6 -outfmt 6 | sort -k12,12nr -k1,1 | awk '{if ($9 < $10) {print $2"\t"$9"\t"$10} else {print $2"\t"$10"\t"$9}}' > $OUTPUT/$PROJECT/blast''$pos''.TSD.bed
awk '{print $1"\t1\t"$2}' > $OUTPUT/$PROJECT/left.bed
awk '{print $1"\t3\t1000"}' > $OUTPUT/$PROJECT/right.bed

#awk '{print $1"\t1\t1000"}' $OUTPUT/$PROJECT/blast''$pos''.TSD.bed > $OUTPUT/$PROJECT/fullregion''$pos''.bed
#$BEDTOOLS/bedtools substract -a $OUTPUT/$PROJECT/fullregion''$pos''.bed -b $OUTPUT/$PROJECT/blast''$pos''.TSD.bed > $OUTPUT/$PROJECT/leftright_''$pos''.bed

if [[ $assembled == "yes" ]] ### For assembled Alu
#fully assembled TE have TSD in 5' and 3' (or it is assumed)
then

#echo "left"
$BEDTOOLS/bedtools getfasta -fi $GENOME -bed $OUTPUT/$PROJECT/left.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $OUTPUT/$PROJECT/$pos"".left.seq
#echo "right"
$BEDTOOLS/bedtools getfasta -fi $GENOME -bed $OUTPUT/$PROJECT/right.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $OUTPUT/$PROJECT/$pos"".right.seq

 	if [[ $direction == "-" ]]
 	then
 		# if - get assemble TE, Reverse and Complement 
		perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$header") $4 | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2 | rev | tr "ATGCatgc" "TACGtacg") > $OUTPUT/$PROJECT/$pos"".TE.seq
 	else
 		# esle if + paste TE in + 
		perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$header") $4 | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' |  cut -f 2) | sed 's/,//g;s/null//g' > $OUTPUT/$PROJECT/$pos"".TE.seq
	fi

 	paste <(echo "$pos") <(sed 's/:/\t/g' <(echo "$pos") | cut -f 1) <(sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print ($3-250)}') <(echo ".") <(sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print ($2-250)"\t"($3+250)}') <(sed 's/,//g' $OUTPUT/$PROJECT/region''$pos''.fasta)  <(paste -d, $OUTPUT/$PROJECT/$pos"".left.seq $OUTPUT/$PROJECT/$pos"".TE.seq $OUTPUT/$PROJECT/$pos"".right.seq | sed 's/,//g') >> $5



else ### For non-assembled
#The TE taken from repbase don't have TSD, so the TSD is to add.


#echo "left"
$BEDTOOLS/bedtools getfasta -fi $GENOME -bed left.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $pos"".left.seq
#echo "right"
$BEDTOOLS/bedtools getfasta -fi $GENOME -bed right.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $pos"".right.seq

 	if [[ $direction == "-" ]]
 	then
 		# if - get assemble TE, Reverse and Complement 
		perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$header") $RM_FASTA | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2 | rev | tr "ATGCatgc" "TACGtacg") > $OUTPUT/$PROJECT/$pos"".TE.seq
 	else
 		# esle if + paste TE in + 
		perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$header") $RM_FASTA | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' |  cut -f 2) | sed 's/,//g;s/null//g' > $OUTPUT/$PROJECT/$pos"".TE.seq
	fi

 	paste <(echo "$pos") <(sed 's/:/\t/g' <(echo "$pos") | cut -f 1) <(sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print ($3-250)}') <(echo ".") <(sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print ($2-250)"\t"($3+250)}') <(sed 's/,//g' $OUTPUT/$PROJECT/region''$pos''.fasta)  <(paste -d, $OUTPUT/$PROJECT/$pos"".left.seq <(echo "$TSD") $OUTPUT/$PROJECT/$pos"".TE.seq <(echo "$TSD") $OUTPUT/$PROJECT/$pos"".right.seq | sed 's/,//g') >> $5

fi

# ### cleaning
# rm left.bed
# rm right.bed
# rm *.seq

done

#rm pre_input

# return delimiter to previous value
IFS=$IFSq_BAK
IFS_BAK=


#### NEED TO FIX --> turn RC the assembled TSD if on neg strand to print a user-friendly table
#awk '{if ($7 == "yes" && $8 == "-") {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10; printf system("echo " $10 " | rev | tr \"ATGCatgc\" \"TACGtacg\"")} else {if ($7 == "yes" && $8 == "?" && $5 == "-") {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10, printf system("echo " $10 " | rev | tr \"ATGCatgc\" \"TACGtacg\"")} else {print $0}}}' pre_input

