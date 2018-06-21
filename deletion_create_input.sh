#! /bin/bash

# usage deletion_create_input.sh TEcordinates_with_bothtsd_cordinates.v.3.3.txt

# find the 1000 bp up and downstream of the TE
left_TE_del=$(awk '{print $9"\t"($10-500)"\t"$10}' $1)
right_TE_del=$(awk '{print $9"\t"($11-1)"\t"($11+499)}' $1)
TE_del=$(awk '{print $9"\t"($10-500)"\t"($11+500)}' $1)


# extract the sequence of them | change the refgenome path by a user input
left=$(bedtools getfasta -fi /Users/clementgoubert/Documents/Postdoc/Lynn/Hetero_correct_project/Ref_Genome/genome/hg19.fa -bed <(echo "$left_TE_del"))
right=$(bedtools getfasta -fi /Users/clementgoubert/Documents/Postdoc/Lynn/Hetero_correct_project/Ref_Genome/genome/hg19.fa -bed <(echo "$right_TE_del"))
noTE=$(paste -d "," <(echo "$left") <(echo "$right") | awk 'getline seq {print seq}' | sed $'s/,//g')
TE=$(bedtools getfasta -fi /Users/clementgoubert/Documents/Postdoc/Lynn/Hetero_correct_project/Ref_Genome/genome/hg19.fa -bed <(echo "$TE_del") | awk 'getline seq {print seq}')

# generate the table
infos=$(awk '{print $7":"$8"-"($8+1)"\t"$1"\t"$2"\t\.\t"($2-500)"\t"($3+500)}' $1)
paste <(echo "$infos") <(echo "$noTE") <(echo "$TE")
