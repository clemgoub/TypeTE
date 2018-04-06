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

# SET PATH FOR FILES AND PROGRAMS
BEDTOOLS="/home/cgoubert/bin/bedtools2/bin"
FullLength="/vbod2/jainy/SGDP/Project2/findgenotype/Assembled_TEsequences.4.5.txt.full_len.fasta"
RefGenome="/vbod2/cgoubert/Correct_Genotypes/Ref_Genome/genome/hg19.fa"
RM_FASTA="/home/cgoubert/bin/RepeatMasker/Libraries/dc20170127-rb20170127/homo_sapiens/refinelib"

### Add line to process Jainy's output
join -a1 -11 -21 <(sort -k1,1 <(join -11 -21 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.2,2.3,2.4,2.5 <(sort -k1,1 Alu_denovo.input) <(sort -k1,1 Jainy_expected_input))) <(sort -k1,1 ../300CEU/position_and_TE) \
| awk '{if (NF)print $0"#SINE/Alu"}' > pre_input

input=$(cat pre_input) # change it so $1 is the output of previous line
# change delimiter (IFS) to new line.
IFS_BAK=$IFS
IFS=$'\n'

for list in $input
do


pos=$(printf '%s\n' "$list" | awk '{print $1}')
#TE=$()
assembled=$(printf '%s\n' "$list" | awk '{print $7}')
TSD=$(printf '%s\n' "$list" | awk '{print $10}' | sed 's/null//g')
direction=$(printf '%s\n' "$list" | awk '{if ($8 == "?") {print $5} else {print $8}}')
header=$(printf '%s\n' "$list" | awk '{print $9}')

sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print $1"\t"($2-500)"\t"($3-500)"\t"$4"\t"$5"-"$6}'  > left.bed
sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print $1"\t"($2+500)"\t"($3+500)"\t"$4"\t"$5"-"$6}'  > right.bed

echo "$pos"
echo "$assembled"
echo "$TSD"
echo "$direction"
echo "$header"

if [[ $assembled == "yes" ]] ### For assembled
then

	#echo "left"
	$BEDTOOLS/bedtools getfasta -fi $Ref_Genome -bed left.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $pos"".left.seq
	#echo "right"
	$BEDTOOLS/bedtools getfasta -fi $Ref_Genome -bed right.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $pos"".right.seq

	if [[ $direction == "-" ]]
	then
		# if reverse paste TE RC | TSD RC
		paste -d, <(perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$header") $FullLength | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | sed "s/$TSD//g" | cut -f 2 | rev | tr "ATGCatgc" "TACGtacg") <(echo "$TSD" | rev | tr "ATGCatgc" "TACGtacg") | sed 's/,//g;s/null//g' > $pos"".TE.seq
	else
		# esle if + paste TE in + | TSD
		paste -d, <(perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$header") $FullLength | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | sed "s/$TSD//g" | cut -f 2) <(echo "$TSD") | sed 's/,//g;s/null//g' > $pos"".TE.seq
	fi

	paste <(echo "$pos") <(sed 's/:/\t/g' <(echo "$pos") | cut -f 1) <(sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print ($3-500)}') <(echo ".") <(sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print ($2-500)"\t"($3+500)}') <(paste -d, $pos"".left.seq $pos"".right.seq | sed 's/,//g')  <(paste -d, $pos"".left.seq $pos"".TE.seq $pos"".right.seq | sed 's/,//g') >> input_test_2



else ### For non-assembled


	#echo "left"
	$BEDTOOLS/bedtools getfasta -fi $Ref_Genome -bed left.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $pos"".left.seq
	#echo "right"
	$BEDTOOLS/bedtools getfasta -fi $Ref_Genome -bed right.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $pos"".right.seq

	if [[ $direction == "-" ]]
	then
		# if reverse paste TE RC | TSD from MELT (RC)
		paste -d, <(perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$header") $RM_FASTA | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2 | rev | tr "ATGCatgc" "TACGtacg" ) <(echo "$TSD") | sed 's/,//g;s/null//g' > $pos"".TE.seq
	else
		# esle if + paste TE in + | TSD from MELT
		paste -d, <(perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$header") $RM_FASTA | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2) <(echo "$TSD") | sed 's/,//g;s/null//g' > $pos"".TE.seq
	fi

	paste <(echo "$pos") <(sed 's/:/\t/g' <(echo "$pos") | cut -f 1) <(sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print ($3-500)}') <(echo ".") <(sed 's/:/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print ($2-500)"\t"($3+500)}') <(paste -d, $pos"".left.seq $pos"".right.seq | sed 's/,//g')  <(paste -d, $pos"".left.seq $pos"".TE.seq $pos"".right.seq | sed 's/,//g') >> input_test_2

fi

done
# return delimiter to previous value
IFS=$IFSq_BAK
IFS_BAK=