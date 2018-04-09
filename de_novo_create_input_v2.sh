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

# usage: ./de_novo_create_input.sh <input>.input <orientTE_extractTE.pl>.outputtable <indentify_mei_from_RM.sh>.outfolder

# SET PATH FOR FILES AND PROGRAMS
BEDTOOLS="/home/cgoubert/bin/bedtools2/bin"
FullLength="/vbod2/jainy/SGDP/Project2/findgenotype/Assembled_TEsequences.4.5.txt.full_len.fasta"
RefGenome="/vbod2/cgoubert/Correct_Genotypes/Ref_Genome/genome/hg19.fa"
RM_FASTA="/home/cgoubert/bin/RepeatMasker/Libraries/dc20170127-rb20170127/homo_sapiens/refinelib"

### Add line to process Jainy's output
join -a1 -11 -21 <(sort -k1,1 <(join -11 -21 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.2,2.3,2.4,2.5 <(sort -k1,1 $1) <(sort -k1,1 $2))) <(sort -k1,1 $3/position_and_TE) \
| awk '{if (NF == 10) {print $1,$2,$3,$4,$5,$6,$7,$8,$10"#SINE/Alu",$6} else {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' > pre_input

input=$(cat pre_input) # change it so $1 is the output of previous line
# change delimiter (IFS) to new line.
IFS_BAK=$IFS
IFS=$'\n'

for list in $input
do


pos=$(printf '%s\n' "$list" | awk '{print $1}')
#TE=$()
assembled=$(printf '%s\n' "$list" | awk '{print $7}')
TSD=$(printf '%s\n' "$list" | awk '{print $NF}' | sed 's/null//g')
direction=$(printf '%s\n' "$list" | awk '{if ($8 == "?") {print $5} else {print $8}}')
header=$(printf '%s\n' "$list" | awk '{print $9}')

# create bed files -500bp ---- BP (left)  BP ---- +500bp (right)
sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print $1"\t"($2-250)"\t"($3-250)"\t"$4"\t"$5"-"$6}'  > left.bed # -250 from left boundary that is already -250 = BP-250bp; -250 from right boundary = BP
sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print $1"\t"($2+250)"\t"($3+250)"\t"$4"\t"$5"-"$6}'  > right.bed # +250 from right boundary that is already +250 = BP+250bp; +250 from left boundary = BP

echo "$pos"
echo "$assembled"
echo "$TSD"
echo "$direction"
echo "$header"

if [[ $assembled == "yes" ]] ### For assembled
then

	#echo "left"
	$BEDTOOLS/bedtools getfasta -fi $RefGenome -bed left.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $pos"".left.seq
	#echo "right"
	$BEDTOOLS/bedtools getfasta -fi $RefGenome -bed right.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $pos"".right.seq

	if [[ $direction == "-" ]]
	then
		# if reverse paste TE RC | TSD RC
		paste -d, <(perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$header") $FullLength | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | sed "s/$TSD//g" | cut -f 2 | rev | tr "ATGCatgc" "TACGtacg") <(echo "$TSD" | rev | tr "ATGCatgc" "TACGtacg") | sed 's/,//g;s/null//g' > $pos"".TE.seq
	else
		# esle if + paste TE in + | TSD
		paste -d, <(perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$header") $FullLength | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | sed "s/$TSD//g" | cut -f 2) <(echo "$TSD") | sed 's/,//g;s/null//g' > $pos"".TE.seq
	fi

	paste <(echo "$pos") <(sed 's/:/\t/g' <(echo "$pos") | cut -f 1) <(sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print ($3-250)}') <(echo ".") <(sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print ($2-250)"\t"($3+250)}') <(paste -d, $pos"".left.seq $pos"".right.seq | sed 's/,//g')  <(paste -d, $pos"".left.seq $pos"".TE.seq $pos"".right.seq | sed 's/,//g') >> input_test_2



else ### For non-assembled


	#echo "left"
	$BEDTOOLS/bedtools getfasta -fi $RefGenome -bed left.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $pos"".left.seq
	#echo "right"
	$BEDTOOLS/bedtools getfasta -fi $RefGenome -bed right.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $pos"".right.seq

	if [[ $direction == "-" ]]
	then
		# if reverse paste TE RC | TSD from MELT (RC)
		paste -d, <(perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$header") $RM_FASTA | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2 | rev | tr "ATGCatgc" "TACGtacg" ) <(echo "$TSD") | sed 's/,//g;s/null//g' > $pos"".TE.seq
	else
		# esle if + paste TE in + | TSD from MELT
		paste -d, <(perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$header") $RM_FASTA | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2) <(echo "$TSD") | sed 's/,//g;s/null//g' > $pos"".TE.seq
	fi

	paste <(echo "$pos") <(sed 's/:/\t/g' <(echo "$pos") | cut -f 1) <(sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print ($3-250)}') <(echo ".") <(sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print ($2-250)"\t"($3+250)}') <(paste -d, $pos"".left.seq $pos"".right.seq | sed 's/,//g')  <(paste -d, $pos"".left.seq $pos"".TE.seq $pos"".right.seq | sed 's/,//g') >> input_test_2

fi

done
# return delimiter to previous value
IFS=$IFSq_BAK
IFS_BAK=

#### NEED TO FIX --> turn RC the assembled TSD if on neg strand to print a user-friendly table
#awk '{if ($7 == "yes" && $8 == "-") {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10; printf system("echo " $10 " | rev | tr \"ATGCatgc\" \"TACGtacg\"")} else {if ($7 == "yes" && $8 == "?" && $5 == "-") {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10, printf system("echo " $10 " | rev | tr \"ATGCatgc\" \"TACGtacg\"")} else {print $0}}}' pre_input

