#! /bin/bash

# make the input table for DE NOVO ALU insertion detected by MELT
#
# to fix: 
# 1) make it to detect the header line (here #71) --> OK
# 2) detect if input vcf is in .gz or not --> OK
# 3) Check if <INS:ME:ALU> is the call in the new MELT. Otherwise just grab all lines for standard MELT run --> YES OK
#
# USAGE: ./input_from_melt.sh genotypes.vcf(.gz) inputname

source parameterfile.init

testfile=$(file $1 | awk '{print $2}')
if [[ $testfile == "gzip" ]]
then

headline=$(zcat $1 | grep -n "#CHROM" | cut -d : -f 1)
TSD=$(zcat $1 | awk -v varline="$headline" -F $'\t' '/<INS:ME:ALU>/ && NR > varline {print $1,$2,$8}' | head -n 1 | sed 's/;/\t/g' |  awk '{for (i=1; i<=NF; ++i) { if ($i ~ "TSD") printf i } } ')
MEINFO=$(zcat $1 | awk -v varline="$headline" -F $'\t' '/<INS:ME:ALU>/ && NR > varline {print $1,$2,$8}' | head -n 1 | sed 's/;/\t/g' |  awk '{for (i=1; i<=NF; ++i) { if ($i ~ "MEINFO") printf i } } ')
zcat $1 | awk -v varline="$headline" -F $'\t' '/<INS:ME:ALU>/ && NR > varline {print $1,$2,$8}' | sed 's/;/\t/g' | awk -v tsd="$TSD" -v meifo="$MEINFO" '{print $1"_"$2,$1,$2, $tsd, $meifo}' | sed 's/,/\t/g;s/TSD=//g;s/MEINFO=//g;s/:/\t/g' | awk '{print $1,$2,$3,$5,$NF,$4}' > $OUTDIR/$PROJECT/$2"".input

else

headline=$(cat $1 | grep -n "#CHROM" | cut -d : -f 1)
TSD=$(cat $1 | awk -v varline="$headline" -F $'\t' '/<INS:ME:ALU>/ && NR > varline {print $1,$2,$8}' | head -n 1 | sed 's/;/\t/g' |  awk '{for (i=1; i<=NF; ++i) { if ($i ~ "TSD") printf i } } ')
MEINFO=$(cat $1 | awk -v varline="$headline" -F $'\t' '/<INS:ME:ALU>/ && NR > varline {print $1,$2,$8}' | head -n 1 | sed 's/;/\t/g' |  awk '{for (i=1; i<=NF; ++i) { if ($i ~ "MEINFO") printf i } } ')
cat $1 | awk -v varline="$headline" -F $'\t' '/<INS:ME:ALU>/ && NR > varline {print $1,$2,$8}' | sed 's/;/\t/g' | awk -v tsd="$TSD" -v meifo="$MEINFO" '{print $1"_"$2,$1,$2, $tsd, $meifo}' | sed 's/,/\t/g;s/TSD=//g;s/MEINFO=//g;s/:/\t/g' | awk '{print $1,$2,$3,$5,$NF,$4}' > $OUTDIR/$PROJECT/$2"".input


fi
