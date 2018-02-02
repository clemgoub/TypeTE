#! /bin/bash

##########################################
# reGenotypeTE - de_novo_create_input.sh #
#                                        #
# This script creates the input file     #
# of the de-novo genotyping script       #
#                                        #
# Author: Clement Goubert                #
# Date: 12/2017                          #
##########################################


list=$(cat position_and_TE_2)
for position in $list
do

pos=$(sed 's/xXx/\t/g' <(echo "$position") | cut -f 1)
TE=$(sed 's/xXx/\t/g' <(echo "$position") | cut -f 3)
TSD=$(sed 's/xXx/\t/g' <(echo "$position") | cut -f 4)
direction=$(sed 's/xXx/\t/g' <(echo "$position") | cut -f 5)

sed 's/xXx/\t/g;s/:/\t/g;s/-/\t/g' <(echo "$position") | awk '{print $1"\t"($2-250)"\t"($3-250)"\t"$4"\t"$5"-"$6}' | sed 's/1KG_7143-/1KG_7143/g' > left.bed
sed 's/xXx/\t/g;s/:/\t/g;s/-/\t/g' <(echo "$position") | awk '{print $1"\t"($2+250)"\t"($3+250)"\t"$4"\t"$5"-"$6}' | sed 's/1KG_7143-/1KG_7143/g' > right.bed

#echo "left"
~/bin/bedtools2/bin/bedtools getfasta -fi /vbod2/cgoubert/Correct_Genotypes/Ref_Genome/genome/hg19.fa -bed left.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $pos"".left.seq
#echo "right"
~/bin/bedtools2/bin/bedtools getfasta -fi /vbod2/cgoubert/Correct_Genotypes/Ref_Genome/genome/hg19.fa -bed right.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $pos"".right.seq
echo "$TE"

if [[ $direction == "-" ]]
then
	# if reverse paste TE RC | TSD 
	paste -d, <(perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$TE#SINE/Alu") ~/bin/RepeatMasker/Libraries/dc20170127-rb20170127/homo_sapiens/refinelib | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2 | rev | tr "ATGCatgc" "TACGtacg" | cut -f2) <(echo "$TSD") | sed 's/,//g;s/null//g' > $pos"".TE.seq
else
	# esle if + paste TE in + | TSD
	paste -d, <(perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$TE#SINE/Alu") ~/bin/RepeatMasker/Libraries/dc20170127-rb20170127/homo_sapiens/refinelib | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2) <(echo "$TSD") | sed 's/,//g;s/null//g' > $pos"".TE.seq
fi

paste <(echo "$pos") <(sed 's/:/\t/g' <(echo "$pos") | cut -f 1) <(sed 's/xXx/\t/g;s/:/\t/g;s/-/\t/g' <(echo "$position") | awk '{print ($3-250)}') <(echo ".") <(sed 's/xXx/\t/g;s/:/\t/g;s/-/\t/g' <(echo "$position") | awk '{print ($2-250)"\t"($3+250)}') <(paste -d, $pos"".left.seq $pos"".right.seq | sed 's/,//g')  <(paste -d, $pos"".left.seq $pos"".TE.seq $pos"".right.seq | sed 's/,//g') >> input_test_2

done
