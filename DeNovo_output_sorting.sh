#/bin/bash

# $1 is the output folder (full path)
# $2 is the folder that contains the file produced by the first pipeline (Jainy)
wd=$(pwd)
mkdir -p $1 # creates the output folder
cd $1 # moves to the output folder

echo "#################################################################"
echo "# MODULE 1: Sort empty locus assembly blast against empty locus #"
echo "#################################################################"

# sorting of the assembly including only the non-discordant pairs, non discordant mates and non-split reads
if [ -s $1/names ]; then rm $1/names; fi
if [ -s $1/bp ]; then rm $1/bp; fi
paste <(date | awk '{print $4}') <(echo "sorting blast outputs...")

### count the output with hit on the empty locus (0 < bp < n)

cd $2/genomenosoftblast # moves to the assembly pipeline output folder



filelist=$(ls -lh | grep 'out.tabular.out' | awk '$5 != 0' | awk '{print $9}') # generates the filelist of blast outputs which have an hit on the empty locus

for file in $filelist
do

sort -k12,12nr $file | sort -u -k1,1| cut -f 4 >> $1/bp 
echo "$file" | sed $'s/\.genomenosoftblast\.out\.tabular\.out//g;s/\./\t/g' >> $1/names

done

paste $1/names $1/bp > $1/emptylocus_hit # the loop and this command generates the files of indiv/locus/hit_in_bp


### count the outputs with no hit (0 bp)
if [ -s $1/emptylocusNohit ]; then rm $1/emptylocusNohit; fi

filelist2=$(ls -lh | grep 'out.tabular.out'| awk '$5 == 0' | awk '{print $9}') # generates the filelist of blast outputs which don't have an hit on the empty locus

for file in $filelist2
do

names2=$(echo "$file" | sed $'s/\.genomenosoftblast\.out\.tabular\.out//g;s/\./\t/g')
paste <(echo "$names2") <(echo "0") >> $1/emptylocusNohit 
done # the loop and this command generates the files of indiv/locus/hit_in_bp

cd $1

cat emptylocus_hit emptylocusNohit | sort -k1,1 -k2,2 > $1/no_softclipped_besthit_empty_blast.out # concatenates hits and no hits to get the full list of hit size for empty locus assembly attempt

paste <(date | awk '{print $4}') <(echo "done")

echo "##################################################################################################"
echo "# MODULE 2: Intersect the discordant mates with the RM track to find the polymorphic TE sequence #"
echo "##################################################################################################"

# will pick up the discordant mates from the bam files, generates bed per individual then per locus; interesct locus bed with RepeatMasker track to find TE.

paste <(date | awk '{print $4}') <(echo "exporting discordant mates position in .bed format...")

cd $2/IGV # switch to the IGV folder of the pipeline to find the full bams
mkdir -p $1/disco_mates # creates the subfolder in this output directory

listofbams=*.bam 
for bam in $listofbams
do
samtools view -F 14 $bam | awk '($3!=$7 && $7!="=")' | awk '{print "chr"$7"\t"$8"\t"($8+100)}' > $1/disco_mates/$bam"".bed # generates individual bed from bam
done

paste <(date | awk '{print $4}') <(echo "concatenate individual discordant reads .bed per locus...") 

cd $1/disco_mates 

positions=$(ls | sed $'s/\./\t/g' | cut -f 2 | sort | uniq) # creates the list of unique positions (locus)
for pos in $positions 
do
cat *$pos*.bed | grep -v "chr\*" > $pos""_concat.bed # concatenates the individual bed per locus
done

paste <(date | awk '{print $4}') <(echo "intersect concatenated locus .bed with RepeatMasker track...")

filelist=*_concat.bed # list the locus bedfile in preparation of the intersection with the RM track

for file in $filelist
do
echo "$file..."
# intersect for each locus the RepeatMasker track and select the best overlap
TE=$(~/bin/bedtools2/bin/bedtools intersect -wao -a <(awk '{print $0"\t"NR}' $file) -b $wd/RepeatMasker_Alu_hg38_lift_hg19.bed | awk '{print $1"_"$2"_"$3"_"$4"\t"$0}' | sort -k1,1 -k9,9nr | sort -u -k1,1 | sed $'s/repName=//g;s/;/\t/g' | cut -f 9 | sort | uniq -c | sort -k1,1nr | grep -v "\." | head -n 1)
echo "$TE"
paste <(echo "$file") <(echo $TE) | sed $'s/_concat\.bed//g' | sed 's/\.\.\///g;s/\t/xXx/g;s/ /xXx/g' >> $1/position_and_TE # print the position | number of discordant mates that support | TE (one line per positio)
done


echo "########################################################################################################"
echo "# MODULE 3: Blast of the discordand reads and the contigs (including split reads) again polymorphic TE #"
echo "########################################################################################################"

cd $1

mkdir -p TE_db
mkdir -p contig_blast
mkdir -p disread_blast

filelist=$(cat position_and_TE)
for file in $filelist
do

cd $1

#extract the info from the position_and_TE file
position=$(echo "$file" | sed 's/xXx/\t/g' | awk '{print $1}')
contig=$(echo "$file" | sed 's/xXx/\t/g' | awk '{print $3"#SINE/Alu"}')
TE=$(echo "$file" | sed 's/xXx/\t/g' | awk '{print $3}')
cov=$(echo "$file" | sed 's/xXx/\t/g' | awk '{print $2}')

paste <(date | awk '{print $4}') <(echo "$position $TE being processed...")

#extract and rename with coverage the predicted polymorphic Alu contig corresponding to the position
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$contig") ~/bin/RepeatMasker/Libraries/dc20170127-rb20170127/homo_sapiens/refinelib > TE_db/$position"".$TE"".$cov"".fasta

#make blast db out of the contig
~/bin/ncbi-blast-2.6.0+-src/rmblast/bin/makeblastdb -in TE_db/$position"".$TE"".$cov"".fasta -out TE_db/$position"".$TE"".$cov"".fasta -dbtype 'nucl'


#################################
#disco reads blast per individual
paste <(date | awk '{print $4}') <(echo "Running discordant reads blast... ")

list_of_files=$(ls /vbod2/jainy/SGDP/Project2/findgenotype/findgenotype_withconcat/reconstructTE/$position/)
	for disco in $list_of_files
	do

	# THE PATH TO THE QUERY SHOULD BE CHANGED TO $2/ WHEN THIS OUTPUT WILL BE INTEGRATED TO THE FINAL PIPELINE [$2/reconstructTE ?]
	~/bin/ncbi-blast-2.6.0+-src/rmblast/bin/blastn -query /vbod2/jainy/SGDP/Project2/findgenotype/findgenotype_withconcat/reconstructTE/$position/$disco -db TE_db/$position"".$TE"".$cov"".fasta -outfmt 6 -out disread_blast/$disco"".out
	sort -k1,1 -k12,12nr disread_blast/$disco"".out | sort -u -k1,1 > disread_blast/sorted_$disco"".out

	done

###############################
#contig blast per individual
## so far need to copy the different contigs in one folder

paste <(date | awk '{print $4}') <(echo " Running contig blast...")

mkdir -p $1/assembled_contigs # creates a folder to gather the assembled contigs / may be removed later

cd $2

rsync ExtractedReads/*/Renamedcontigs/*.rename.fasta $1/assembled_contigs

###loop each individual

cd $1/assembled_contigs

list_of_files2=$(ls *$position"".rename.fasta)

	for empty in $list_of_files2

	do

	nam=$(echo "$empty" | sed 's/\//\t/g' | cut -f 3)
	~/bin/ncbi-blast-2.6.0+-src/rmblast/bin/blastn -query $empty -db $1/TE_db/$position"".$TE"".$cov"".fasta -outfmt 6 -out $1/contig_blast/$nam"".out
	sort -k1,1 -k12,12nr $1/contig_blast/$nam"".out | sort -u -k1,1 > $1/contig_blast/sorted_$nam"".out

	done #for individual loop

done #for position loop

cd $wd

echo "#####################################################"
echo "# MODULE 4: Collect the discordant reads statistics #"
echo "#####################################################"

# Goes back to the IGV folder to get the statistics on the total number of pared reads

paste <(date | awk '{print $4}') <(echo "Collecting stats...")

cd $2/IGV # switch to the IGV folder of the pipeline to find the full bams
mkdir -p $1/disco_stats # creates the subfolder in this output directory
if [ -s $1/disco_stats/disco_stats ]; then rm $1/disco_stats/disco_stats; fi # removes the output if already present

filelist=*.bam
for file in $filelist


do

indiv=$(echo "$file" | sed 's/\./\t/g' | cut -f 1)
locus=$(echo "$file" | sed 's/\./\t/g' | cut -f 2)
stats=$(samtools flagstat $file | awk '{print $1}' | perl -anF'\t|\n' -e'$n=@F-1if!$n;for(0..$n){push@{$$m[$_]},$F[$_]}''END{print map{join"\t",@$_,"\n"}@$m}' | cut -f 1,9,12,13)
echo "stats for bam $file..."
paste <(echo "$indiv") <(echo "$locus") <(echo "$stats")  >> $1/disco_stats/disco_stats

done


echo "######################################"
echo "# MODULE 5: Filter the blast outputs #"
echo "######################################"

paste <(date | awk '{print $4}') <(echo "Filtering blast outputs...")

####### Filter the discordant read blast on polymorphic TE according to a threshold set as $3 in the bash script (ex 98 for 98% identity)
cd $1/disread_blast

paste <(date | awk '{print $4}') <(echo "Filtering discordant reads blast...")

if [ -s $1/disco_blast_TE_hit ]; then rm $1/disco_blast_TE_hit; fi # removes the output if already present

hits=$(ls -lh | grep -v ' 0' | awk '{print $NF}' | grep 'sorted') ### filter the hit output

for hit in $hits
do

thresh=$(awk -v th=$3 '$3 >= th' $hit | wc -l)

key=$(echo "$hit" | sed 's/sorted_//g;s/\.discordant/\t/g;s/\./_/g' | cut -f 1)

paste <(echo "$key") <(echo "$thresh")  >> $1/disco_blast_TE_hit

done

no_hits=$(ls -lh | grep ' 0' | awk '{print $NF}' | grep 'sorted') ### and created "0" lines on the no hit, empty outputs

for no in $no_hits
do

key=$(echo "$no" | sed 's/sorted_//g;s/\.discordant/\t/g;s/\./_/g' | cut -f 1)

paste <(echo "$key") <(echo "0") >> $1/disco_blast_TE_hit

done


####### Filter the contigs blast against the polymorrphic TE according to $4 on the bash script

cd $1/contig_blast
paste <(date | awk '{print $4}') <(echo "Filtering contig blast...")
if [ -s $1/contig_blast_TE_hit ]; then rm $1/contig_blast_TE_hit; fi

hits=$(ls -lh | grep -v ' 0' | awk '{print $NF}' | grep 'sorted')

for hit in $hits
do

thresh=$(awk -v th=$3 '$4 >= th' $hit | wc -l)

key=$(echo "$hit" | sed 's/sorted_//g;s/\.rename/\t/g;s/\./_/g' | cut -f 1)

paste <(echo "$key") <(echo "$thresh") >> $1/contig_blast_TE_hit

done

no_hits=$(ls -lh | grep ' 0' | awk '{print $NF}' | grep 'sorted')

for no in $no_hits
do

key=$(echo "$no" | sed 's/sorted_//g;s/\.rename/\t/g;s/\./_/g' | cut -f 1)

paste <(echo "$key") <(echo "0") >> $1/contig_blast_TE_hit

done


###

#sorting the two datasets 
join -11 -21 -a1 <(sort -k1,1 disco_blast_TE_hit) <(sort -k1,1 contig_blast_TE_hit) | sed $'s/ /\t/g' | awk '{if ($3 == "") {print $0"\t0"} else {print $0}}' > $1/blast_vs_mei_disco_and_contig


echo "################################################"
echo "# MODULE 6: Generating summary output for test #"
echo "################################################"

cd $1

#join the empty locus assembly blast against the empty reference with the total number of paires reads per extracted region
bp_and_disco=$(join -11 -21 <(sort -k1,1 <(awk '{print $1"_"$2"\t"$3}' no_softclipped_besthit_empty_blast.out)) <(sort -k1,1 <(awk '{print $1"_"$2"\t"$4}' disco_stats/disco_stats)))
#join this with the number of discordant mates mapping the polymorphic TE as well as the number of uniques hits of the blast of the contigs vs the polymorphic TE (to assess presence of the polymorphic TE)
cat <(echo "ind pos bp paired disco contig") <(join -11 -21 <(sort -k1,1 <(echo "$bp_and_disco")) <(sort -k1,1 blast_vs_mei_disco_and_contig)) | sed 's/ /\t/g;s/_/\t/g' > statistics_table




