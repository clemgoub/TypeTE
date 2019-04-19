#! /bin/bash

##########################################
# reGenotypeTE - de_novo_create_input.sh #
# for SVA and L1, without assembly       #
#                                        #
# This script creates the input file     #
# of the de-novo genotyping script using #
# MELT info only                         #
#                                        #
# Author: Clement Goubert                #
# Date v1: 04/16/2019                    #
##########################################

# usage: ./de_novo_create_input_SVA-L1.sh $OUTDIR/$PROJECT.input <output.file>

# SET PATH FOR FILES AND PROGRAMS 
source parameterfile_NoRef.init
input=$(cat $1)

# change delimiter (IFS) to new line.
IFS_BAK=$IFS
IFS=$'\n'

for list in $input
do


pos=$(printf '%s\n' "$list" | awk '{print $1}')
TSD=$(printf '%s\n' "$list" | awk '{print $6}' | sed 's/null//g;s/d//g;s/0tsd//g')
direction=$(printf '%s\n' "$list" | awk '{print $5}')
header=$(printf '%s\n' "$list" | awk '{print $7}')

##DEBUG
#echo "$pos"
#echo "$TSD"
#echo "$direction"
#echo "$header"


# create bed files -500bp ---- BP (left)  BP ---- +500bp (right)
sed 's/_/\t/g' <(echo "$pos") | awk '{print $1"\t"($2-500)"\t"$2}'  > left.bed # -500 from BP
sed 's/_/\t/g' <(echo "$pos") | awk '{print $1"\t"$2"\t"($2+500)}'  > right.bed # +500 from BP

	$BEDTOOLS/bedtools getfasta -fi $GENOME -bed left.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $pos"".left.seq
	$BEDTOOLS/bedtools getfasta -fi $GENOME -bed right.bed | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2  > $pos"".right.seq

	if [[ $direction == "-" ]]
	then
		# if reverse paste TE RC | TSD from MELT (RC)
		paste -d, <(perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$header") $OUTDIR/$PROJECT""/*.consensi.fa | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2 | rev | tr "ATGCatgc" "TACGtacg" ) <(echo "$TSD") | sed 's/,//g;s/null//g' > $pos"".TE.seq
##DEBUG
#		echo "$tsd"
#		echo "$direction"
	else
		# esle if + paste TE in + | TSD from MELT
		paste -d, <(perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$header") $OUTDIR/$PROJECT""/*.consensi.fa | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | cut -f 2) <(echo "$TSD") | sed 's/,//g;s/null//g' > $pos"".TE.seq
##DEBUG
#		echo "$tsd"
#              	echo "$direction"
	fi

	paste <(echo "$pos") <(sed 's/_/\t/g' <(echo "$pos")) <(echo ".") <(sed 's/_/\t/g;s/-/\t/g' <(echo "$pos") | awk '{print ($2-500)}') <(sed 's/_/\t/g' <(echo "$pos") | awk '{print ($2+500)}') <(paste -d, $pos"".left.seq $pos"".right.seq | sed 's/,//g')  <(paste -d, $pos"".left.seq $pos"".TE.seq $pos"".right.seq | sed 's/,//g') >> $2
	#tail -n 1 $2

### cleaning
rm left.bed
rm right.bed
rm *.seq

done

#rm pre_input

# return delimiter to previous value
IFS=$IFSq_BAK
IFS_BAK=

