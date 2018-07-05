#! /bin/bash

# usage deletion_create_input.sh TEcordinates_with_bothtsd_cordinates.v.3.3.txt

source parameterfile_Del.init

# find the 1000 bp up and downstream of the TE

insertions=$(cat $1)
for locus in $insertions
do

locus=$(awk '{print $1}' $insertions)
left=$(awk '{print $2}' $insertions)
right=$(awk '{print $3}' $insertions)
TSD=$(awk '{print $4}' $insertions)

if [[ TSD == "noTSDs" ]]
	
	then

	left_TE_del=$(sed 's/:/\t/g;s/-/\t/g' <(echo "$left") | awk '{print $1,($2-500),$2}'  > left.bed)
	right_TE_del=$(sed 's/:/\t/g;s/-/\t/g' <(echo "$right") | awk '{print $1,$3,($3+250)}'  > right.bed)
	TE_del=$(sed 's/:/\t/g;s/-/\t/g' <(echo "$left")) # the left and right TSD coordinates in case of "noTSDs" are actually the start and stop of the TE of the input file


	# extract the sequence of them | change the refgenome path by a user input
	left=$($BEDTOOLS/bedtools getfasta -fi $GENOME -bed <(echo "$left_TE_del"))
	right=$($BEDTOOLS/bedtools getfasta -fi $GENOME -bed <(echo "$right_TE_del"))
	noTE=$(paste -d "," <(echo "$left") <(echo "$right") | awk 'getline seq {print seq}' | sed $'s/,//g')
	TE=$($BEDTOOLS/bedtools getfasta -fi $GENOME -bed <(echo "$TE_del") | awk 'getline seq {print seq}')


	else
		
fi
left_TE_del=$(awk '{print $9"\t"($10-500)"\t"$10}' $1)
right_TE_del=$(awk '{print $9"\t"($11-1)"\t"($11+499)}' $1)
TE_del=$(awk '{print $9"\t"($10-500)"\t"($11+500)}' $1)


# extract the sequence of them | change the refgenome path by a user input
left=$($BEDTOOLS/bedtools getfasta -fi $GENOME -bed <(echo "$left_TE_del"))
right=$($BEDTOOLS/bedtools getfasta -fi $GENOME -bed <(echo "$right_TE_del"))
noTE=$(paste -d "," <(echo "$left") <(echo "$right") | awk 'getline seq {print seq}' | sed $'s/,//g')
TE=$($BEDTOOLS/bedtools getfasta -fi $GENOME -bed <(echo "$TE_del") | awk 'getline seq {print seq}')

# generate the table
infos=$(awk '{print $7":"$8"-"($8+1)"\t"$1"\t"$2"\t\.\t"($2-500)"\t"($3+500)}' $1)
paste <(echo "$infos") <(echo "$noTE") <(echo "$TE")

done