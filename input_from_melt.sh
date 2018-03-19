#! /bin/bash

# make the input table for DE NOVO ALU insertion detected by MELT
#
# to fix: 
# 1) make it to detect the header line (here #71)
# 2) detect if input vcf is in .gz or not
# 3) Check if <INS:ME:ALU> is the call in the new MELT. Otherwise just grab all lines for standard MELT run
#
# USAGE: ./input_from_melt.sh genotypes.vcf(.gz) inputname

zcat $1 | awk -F $'\t' '/<INS:ME:ALU>/ && NR > 71 {print $1,$2,$8}' | sed 's/TSD=/\t/g;s/;AN/\t/g;s/MEINFO=/\t/g;s/,/\t/g;s/;NS/\t/g' | awk '{print $1"_"$2,$1,$2,$5,$8,$10}' > $2"".input