#! /bin/bash

# make the input table for reference ALU insertion detected by MELT
#
# USAGE: ./input_from_melt_Del.sh genotypes.vcf(.gz) inputname

source parameterfile_Ref.init

testfile=$(file $1 | awk '{print $2}')
if [[ $testfile == "gzip" ]]
then

headline=$(zcat $1 | grep -n "#CHROM" | cut -d : -f 1)
MEINFO=$(zcat $1 | awk -v varline="$headline" -F $'\t' 'NR > varline {print $1,$2,$8}' | head -n 1 | sed 's/;/\t/g' |  awk '{for (i=1; i<=NF; ++i) { if ($i ~ "SVTYPE") printf i } } ')
zcat $1 | awk -v varline="$headline" -F $'\t' 'NR > varline {print $1,$2,$8}' | sed 's/;/\t/g' | awk -v tsd="$TSD" -v meifo="$MEINFO" '{print $1"_"$2,$1,$2, $meifo}' | sed 's/,/\t/g;s/TSD=//g;s/SVTYPE=//g' | awk '{print $1,$2,$3,$4}' > $OUTDIR/$PROJECT/$2"".input

else

headline=$(cat $1 | grep -n "#CHROM" | cut -d : -f 1)
MEINFO=$(cat $1 | awk -v varline="$headline" -F $'\t' 'NR > varline {print $1,$2,$8}' | head -n 1 | sed 's/;/\t/g' |  awk '{for (i=1; i<=NF; ++i) { if ($i ~ "SVTYPE") printf i } } ')
cat $1 | awk -v varline="$headline" -F $'\t' 'NR > varline {print $1,$2,$8}' | sed 's/;/\t/g' | awk -v tsd="$TSD" -v meifo="$MEINFO" '{print $1"_"$2,$1,$2, $meifo}' | sed 's/,/\t/g;s/TSD=//g;s/SVTYPE=//g' | awk '{print $1,$2,$3,$4}' > $OUTDIR/$PROJECT/$2"".input


fi
