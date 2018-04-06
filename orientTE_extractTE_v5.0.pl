#!/usr/bin/perl -w
#######################################################
# Author :  Jainy Thomas
# date   :  November 2017
# email  :  jainythomas1@gmail.com
# Pupose :  to concatenate all reads from a locus from all individuals and to find orientation of the TE (plus or minus)
#			and to extract TE sequences
#           
#####################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Cwd;
use Data::Dumper;
use Bio::SeqIO;
use Bio::DB::Fasta;
use File::Copy;
use List::MoreUtils qw(uniq);
use feature 'fc';


my $version = "5.0";
my $scriptname = "orientTE_extractTE.pl";
my $changelog = "
#   - v1.0 = 3 November 2017 
#	- v2.5 = 8 November 2017 
#			   extract full length TE from the assembled sequences and prints to a file (didnt finish it)
#   - v3.0 = 25 January 2018
#				changed the assembly method from cap3 to other next gene seq assemblers 
#	- v4.0 = 31 January 2018
#				extraction of TEsequences also made possible
#	- v4.1 = 31 January 2018
#				reverse compliment the flank for if the TE is in the minus orientation when blasted
#	- v4.2 = 2 February 2018
#				changed how the TEsequences read by the script, from the table and from folder with fasta sequences	
#	- v4.5 = 5 February 2018
#				generates a output table with prediction and if a full lengthTE was able to extract (yes or no), output the fasta 
#				sequence with both TSDs in one file and the others in another file (header locus,queryTE,orientation, TE)
#	- v5.0 = 4 April 2018
#				modifies the output table with TSD and name of the sequence

\n";

my $usage = "\nUsage [$version]: 
    perl $scriptname  -d <directory> -g <directory> -t <TEdirectory> -l <list of TEinsertion table> [-p <path of the outputdirectory>][-o <output file>] [-v] [-c] [-h] 
	
    MANDATORY ARGUMENTS:	
    -d,--readdir    (STRING) input directory with all reads for all loci <path of the directory needs to be given>
	-g,--genomedir  (STRING) input directory with the corresponding reference genome sequence   
    -t,--TEdir	  	(STRING) input directory with TE sequences (separate fasta files with location name)
    -l,--list		(STRING) input file containing the information on TE insertion
    
    
    OPTIONAL ARGUMENTS:
    -p,--path   (STRING) output directory name (path)
                         Default = <current working directory>
    -o,--output (STRING) output file with strand prediction
    -te,--TEout (STRING) output file with Assembled TE sequences only if nearly full length
    -c,--chlog  (BOOL)   Print updates
    -v,--v      (BOOL)   Print version if only option
    -s,--verbose(BOOL)   The script will talk to you
    -h,--help>  (BOOL)   Print this usage\n\n";

#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($rdir,$gdir,$TEdir,$listte,$teout,$path,$fineout,$verbose,$help,$v,$chlog);
GetOptions ('d=s' => \$rdir,
			'g=s' => \$gdir,
			't=s' => \$TEdir,
			'l=s' => \$listte,
            'p=s' => \$path,
            'o=s' => \$fineout,
           'te=s' => \$teout,
            'c'   => \$chlog, 
            'h'   => \$help,
            's'   => \$verbose, 
            'v'   => \$v);

#check step to see if mandatory argument is provided + if help/changelog
die "\n Script $scriptname version $version\n\n" if ((! $rdir)  && (! $help) && (! $rdir)   && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ((! $rdir) ||(! $gdir) ||(! $TEdir) || ($help));
my $cwd = getcwd();
$path = $cwd if (!$path) ;
$fineout = "$path/genomeloc.strand.prediction.$version.txt" if (! $fineout);
$teout = "$path/Assembled_TEsequences.$version.txt" if (! $teout);
die "\n -d $rdir does not exist?\n\n"  if (! -e $rdir);
die "\n -g $gdir does not exist?\n\n"  if (! -e $gdir);
die "\n -t $TEdir does not exist?\n\n"  if (! -e $TEdir);
$rdir = $1 if ($rdir =~ /^(.*)\/$/);
$gdir = $1 if ($gdir =~ /^(.*)\/$/);
$TEdir = $1 if ($TEdir =~ /^(.*)\/$/);

#my $logfile = "$path/$genomloc.log";
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $CAP3pro = "/home/jainy/software/CAP3";#vader server
my $BLASTpro = "/home/jainy/software/ncbi-blast-2.6.0+/bin";
my $miniapro = "/home/jainy/software/minia-v2.0.7-Source/build/bin";
my $SPAdepro = "/home/jainy/software/SPAdes-3.11.1-Linux/bin";#vaderserver,Yodaserver

#-d "$path/Assembled_TEreads"?die "$path/Assembled_TEreads already exist? please delete/rename the existing folder\n":make_path ("$path/Assembled_TEreads");
my $genomloc;
my $directory;
my @zerofiles;
my $R1out;
my %TEinfo;
my $contigmaxname; 
my $contigreadnum;
my $totalreadnum; 
my $outpath;
my %strandprediction;
my %extractcontiginfo;
my %alldirectory;
#load TE info
%TEinfo = &load_file($listte);
#print Dumper %TEinfo;
#concatenate all the files in the directory 

my @dir =`ls $rdir`;
foreach $directory (@dir) {
	chomp ($directory);
	$alldirectory{$directory} =1;
	print STDERR "the directory now analysing is $directory\n";
	
	make_path ("$path/Assembled_TEreads/$directory") ;
	$outpath = "$path/Assembled_TEreads/$directory";
	#my @files = `ls $rdir/$directory`;
	my @R1files = glob ("'$rdir/$directory/*R1.fasta'");
	my @R2files = glob ("'$rdir/$directory/*R2.fasta'");
	my @UPfiles = glob ("'$rdir/$directory/*Up.fasta'");
	#print STDERR " files loaded are @files\n";
	#%contiginfo = ();
	my $R1out = "$path/Assembled_TEreads/$directory/$directory.concatenated.R1.fasta";
	my $R2out = "$path/Assembled_TEreads/$directory/$directory.concatenated.R2.fasta";
	my $Upout = "$path/Assembled_TEreads/$directory/$directory.concatenated.Up.fasta";
	#my $renamepath = "$path/Assembled_TEreads/$directory";
	&concatenatefiles($R1out,@R1files);
	&concatenatefiles($R2out,@R2files);
	&concatenatefiles($Upout,@UPfiles);
	system ("cat $path/Assembled_TEreads/$directory/$directory.concatenated.R1.fasta $path/Assembled_TEreads/$directory/$directory.concatenated.R2.fasta $path/Assembled_TEreads/$directory/$directory.concatenated.Up.fasta > $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta") == 0 or die ("unable to concatenate discordantsplitreads $directory \n");

	my $R1filesize = (-s "$R1out");
	my $R2filesize = (-s "$R2out");
	
	#assemble the concatenated sequence
	
	# print STDERR "assemble the concatenated sequence\n";
# 	if (($R1filesize > 0) || ($R1filesize > 0)) {
# 		#Assemble only the mapped reads		   
# 		system ("$SPAdepro/spades.py -1 $R1out -2 $R2out -s $Upout --careful --only-assembler -o $path/Assembled_TEreads/$directory/$directory.allreadsSPAdeout") == 0 or 
# 		system ("$SPAdepro/dipspades.py -1 $R1out -2 $R2out -s $Upout --only-assembler --expect-rearrangements -o $path/Assembled_TEreads/$directory/$directory.allreadsdipSPAdeout") == 0 or 
# 		system ("$miniapro/minia -in $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta -kmer-size 45 -abundance-min 3 -out $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta_k45_ma3") == 0 or 
# 		system ("$CAP3pro/cap3 $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta > $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.asmbl.fasta") == 0 or die ("unable to assemble fasta $directory \n");
# 		if (-e "$path/Assembled_TEreads/$directory/$directory.allreadsSPAdeout/scaffolds.fasta") {
# 			#Rename scaffolds.fasta 
# 			copy("$path/Assembled_TEreads/$directory/$directory.allreadsSPAdeout/scaffolds.fasta", "$path/Assembled_TEreads/$directory/$directory.allreads.scaffolds.fasta") or die "Copy failed scaffolds.fasta $directory:$!";
# 			#rename the fasta sequence with its filename
# 			my $assembledfile = "$directory.allreads.scaffolds.fasta";
# 			&renameseq_filename($assembledfile,$outpath,$directory);
# 		} elsif (-e "$path/Assembled_TEreads/$directory/$directory.allreadsdipSPAdeout/dipspades/consensus_contigs.fasta") { 
# 			copy("$path/Assembled_TEreads/$directory/$directory.allreadsdipSPAdeout/dipspades/consensus_contigs.fasta", "$path/Assembled_TEreads/$directory/$directory.allreads.consensus_contigs.fasta") or die "Copy failed consensus_contigs.fasta $directory:$!";
# 			my $assembledfile = "$directory.allreads.consensus_contigs.fasta";
# 			&renameseq_filename($assembledfile,$outpath,$directory);
# 		} elsif (-e "$path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta_k45_ma3.contigs.fa") { 
# 			my $assembledfile = "$directory.concatenated.allreads.fasta_k45_ma3.contigs.fa";
# 			&renameseq_filename($assembledfile,$outpath,$directory);
# 		} elsif (-e "$path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta.cap.contigs") {
# 			my $assembledfile = "$directory.concatenated.allreads.fasta.cap.contigs";
# 			&renameseq_filename($assembledfile,$outpath,$directory);
# 		}
# 	} else {
# 	 push (@zerofiles,$R1out);
# 	}
	#blast the contig with TE sequence to find the strand
	
	my $queryte = &find_TEseq($directory);
	print STDERR "the query file is $queryte.fasta\n";
	&blast_file ($TEdir,"$queryte.fasta",$directory);	#blast with the query TEsequence
	print STDERR "the query is $TEdir/$queryte.fasta\n";
	print STDERR "the directory is $directory\n";
	my $blastfile = "$outpath/$directory.$queryte.fasta.tabular.blast.out";
	my $TEblast = &parse_blast($blastfile);
	#print Dumper %$TEblast;
	#parse the blastoutput for max score
	&blast_file ($gdir,"$directory.extract.seq.fa",$directory);
	my $blastgenomefile = "$outpath/$directory.extract.seq.fa.tabular.blast.out";
	my $genomeblast = &parse_blast($blastgenomefile);
	#print Dumper %$genomeblast;
	#predict the orientation of TE insertion
	&compare_twoHOH($directory,$TEblast,$genomeblast);	
}

print Dumper %extractcontiginfo;
print Dumper %alldirectory;

print STDERR "DONE with prediction of orientation.... \n";

print STDERR "\n\nstarting to extract TEsequences.... \n";

my $left = 5;
my $right = 3;
my @matchingtsdlist;
my @full_lengTEs;
my @partial_teseq;
#my $teseq;
&loadseq_tohash();
print Dumper %extractcontiginfo;
#find the contigs with TEsequence and flanking sequence and extract the flanking sequence
foreach my $locus (sort keys %extractcontiginfo) {
	print "locus is $locus\n";
	my ($genomeloc,$rest) = split (/\./,$locus,2);
	my ($leftkmer);
	my ($rightkmer);
	my $leftflankseq;
	my $rightflankseq;
	my $start;
	my $end;
	my $strand;
	my $seq;
	my $qTE;
	my $newlocus;
	my $teseq;
	my $testart;
	my $teend;
	my $tsdseq;
	my $tsdsucc;
	my $contlen;
	foreach my $info (keys %{$extractcontiginfo{$locus}}) {
		#print "locus is $info\n";
		$start = $extractcontiginfo{$locus}{$info} if ($info eq 'sstart');
		$end = $extractcontiginfo{$locus}{$info} if ($info eq 'send');
		$strand = $extractcontiginfo{$locus}{$info} if ($info eq 'strand');
		#print STDERR " strand is $strand\n"; 
		$seq = $extractcontiginfo{$locus}{$info} if ($info eq 'seq');
		$qTE = $extractcontiginfo{$locus}{$info} if ($info eq 'qid');
		$contlen = $extractcontiginfo{$locus}{$info} if ($info eq 'slen');
		
	}
	$newlocus = $locus.".".$start.".".$end.".".$strand if ($strand eq 'plus');
	$newlocus = $locus.".".$end.".".$start.".".$strand if ($strand eq 'minus');
	my $seq_obj = Bio::Seq->new( -display_id => $newlocus, -seq => $seq); # create object with seq
	if ($strand eq "plus") {# here I am extracting flanking sequence of the TE
		my $leftstart = $start - 40;
		my $leftend = $start + 5;
		my $rightend = $end + 40;
		my $rightstart = $end - 5;
		print "for a plus $leftstart\t$leftend\t$rightstart\t$rightend\n";
		$leftstart = 1 if ($leftstart <= 0);
		$rightend = $contlen if ($rightend > $contlen);
		$leftflankseq = $seq_obj->subseq($leftstart,$leftend) ;
		$rightflankseq = $seq_obj->subseq($rightstart,$rightend);
	} elsif ($strand eq "minus") {
		my $leftstart = $end - 40;
		my $leftend = $end + 5;
		my $rightend = $start + 40;
		my $rightstart = $start - 5;
		
		print "for a minus $leftstart\t$leftend\t$rightstart\t$rightend\n";
		$leftstart = 1 if ($leftstart <= 0);
		$rightend = $contlen if ($rightend > $contlen);
		my $leflankseq = $seq_obj->subseq($leftstart,$leftend);
		my $riflankseq = $seq_obj->subseq($rightstart,$rightend);
		$rightflankseq = &revcom_seq($leflankseq);
		$leftflankseq = &revcom_seq($riflankseq);
	}
	print "$locus left seq is $leftflankseq\n";
	print "$locus right seq is $rightflankseq\n";
	#identifies TSD sequence
	$leftkmer = &generate_kmer_pos($leftflankseq,$left,$newlocus) if defined ($leftflankseq);
	$rightkmer = &generate_kmer_pos($rightflankseq,$right,$newlocus) if defined ($rightflankseq);
	#print Dumper %$leftkmer;
	#print Dumper %$rightkmer;
	&compare_kmers($leftkmer,$rightkmer);
	my $tsdinfo = $matchingtsdlist[0];
	@matchingtsdlist = (); 
	#print "tsd for $newlocus is $tsdinfo\n";
	if (defined ($tsdinfo)) {
		($testart,$teend,$tsdseq,$tsdsucc) = &findcord_bothtsds($tsdinfo);
		print "$locus output from findcordbothtsd is $testart, $teend\n";
		$teseq = $seq_obj->subseq($testart,$teend);
		$teseq = &revcom_seq($teseq) if ($strand eq 'minus');
		if ($tsdsucc eq "yes") {
			
			my $newID = join (".", $genomeloc,$tsdseq,$qTE);
			&checkprediction($genomeloc,$tsdsucc,$newID);
			my $seqout_obj = Bio::Seq->new( -display_id => $newID, -seq => $teseq);
			push (@full_lengTEs,$seqout_obj);
			
		} elsif ($tsdsucc eq "no") {
			
			my $newid =  join (".", $locus,$qTE,"noTSD");
			&checkprediction($genomeloc,$tsdsucc,$newid);
			my $seqout_obj = Bio::Seq->new( -display_id => $newid, -seq => $teseq);
			push (@partial_teseq,$seqout_obj);
		}
	} else {
		$teseq = $seq_obj->subseq($start,$end) if ($strand eq 'plus');
		$teseq = $seq_obj->subseq($end,$start) if ($strand eq 'minus');
		$teseq = &revcom_seq($teseq) if ($strand eq 'minus');
		my $new_id =  join (".", $locus,$qTE,"noTSD");
		my $seqout_obj = Bio::Seq->new( -display_id => $new_id, -seq => $teseq);
		push (@partial_teseq,$seqout_obj);
		&checkprediction($genomeloc,"no",$new_id);
	}
}
#compare two hashes of alldirectory strand prediction
#print the strandprediction hash
&printarray_fasta ("$teout.full_len.fasta",@full_lengTEs);
&printarray_fasta ("$teout.partial.fasta",@partial_teseq);
&find_no_te_cand();
print Dumper %strandprediction;
&print_keyvalue_HOH($fineout,\%strandprediction);

exit;
#-----------------------------------------------------------------------------
#----------------------------------- SUB -------------------------------------
#-----------------------------------------------------------------------------
sub checkprediction {
	my ($gloca,$fTEstatus,$ID) = @_;
	if (exists ($strandprediction{$gloca}))  {
		$strandprediction{$gloca}{'fulllengTE'} = $fTEstatus;
		$strandprediction{$gloca}{'seq_id'} = $ID;
	} else {
		$strandprediction{$gloca}{'fulllengTE'} = $fTEstatus;
		$strandprediction{$gloca}{'strand'} = '?';
		$strandprediction{$gloca}{'seq_id'} = $ID;
	}
}
sub load_file {
	my ($file1) = @_;
	my %stefile;
	open (my $th, "<", $file1) or confess "\n ERROR (main): could not open to read $file1 $!\n";
		while (my $data = <$th>) { #reading the table
			chomp $data; #storing the table values to the $data and removing the extraline
			my @name = split(/\s+/,$data); # splitting the data based on tab and storing into the arrray
			my $dirc = $name[0];#first column is genome location
			$stefile{$dirc} = $name[2]; #name of the TE inserted is loaded 
		}
	return (%stefile);
	close $th;	
}
sub blast_file {
	my ($qpath,$qloc,$directory) = @_;
	unless (-e "$outpath/Renamed_Assembledseq/$directory.rename.fasta.nhr") {
		system ("$BLASTpro/makeblastdb -in $outpath/Renamed_Assembledseq/$directory.rename.fasta -input_type fasta -dbtype nucl") == 0 or die ("unable to makeblastdb on $outpath/Renamed_Assembledseq/$directory.rename.fasta \n");
	}
	if ($qloc eq "$directory.extract.seq.fa") {
		system ("$BLASTpro/blastn -db $outpath/Renamed_Assembledseq/$directory.rename.fasta -query $qpath/$qloc -evalue 0.0001 -outfmt \"6 qseqid sseqid pident qlen length slen sstrand qstart qend sstart send qcovs mismatch gapopen evalue bitscore \" -out $outpath/$qloc.tabular.blast.out");
	} else {
		system ("$BLASTpro/blastn -db $outpath/Renamed_Assembledseq/$directory.rename.fasta -query $qpath/$qloc -evalue 0.0001 -outfmt \"6 qseqid sseqid pident qlen length slen sstrand qstart qend sstart send qcovs mismatch gapopen evalue bitscore \" -out $outpath/$directory.$qloc.tabular.blast.out");
	}
}	
sub concatenatefiles {
	my ($outputfile, @files) = @_;
	foreach $genomloc (@files){
		chomp ($genomloc);
		#print STDERR "the file now analysing is $genomloc\n";
		#system ("cat $rdir/$directory/$genomloc >> $out") == 0 or die ("unable to concatenate $genomloc \n");
		system ("cat $genomloc >> $outputfile") == 0 or die ("unable to concatenate $genomloc \n");#when using glob path is included so no need to give path
	}
}
sub find_TEseq {
	my ($gloc) = @_;
	#print STDERR "the sequence TE is $gloc";
	my $TEseqfile;
	if (exists ($TEinfo{$gloc}))  {
		$TEseqfile = $TEinfo{$gloc};
	}
	else {
		die  "Not able to get the TEseq file \n";
	}
	return ($TEseqfile);
}
sub renameseq_filename {
	my ($contigfile,$fpath,$directory) = @_;
	print STDERR "$fpath and $contigfile\n";
	make_path  ("$fpath/Renamed_Assembledseq");
	open (my $bhout, ">","$fpath/Renamed_Assembledseq/$directory.rename.fasta") or die "\n ERROR (main): could not open to write $fpath/Renamed_Assembledseq/$directory.rename.fasta $!\n";
	open (my $bh, "<", "$fpath/$contigfile") or confess "\n ERROR (main): could not open to read $contigfile $!\n";
		while(my $dataline = <$bh>) {
			chomp($dataline);
			#print STDERR "$line\n";
				if ($dataline =~ m/^\>\w+\d+/) {
					$dataline =~ s/^\>(\w+\d+)/\>$directory\.$1/;
					print $bhout "$dataline\n";
				}
				else {
					print $bhout "$dataline\n";
				}
		}
	close $bh;
	close $bhout;
}
sub parse_blast {
	#parse blast tabular output unique hit with highest blast score
	my $blast = shift;
	my %d = ();
	open (my $bl, "<", $blast) or confess "\n ERROR (sub parse_blast): could not open to read $blast $!\n";
	LINE: while (my $line = <$bl>) {
		chomp $line;
		next LINE if ($line !~ /\w/);#  0      1        2    3     4         5     6      7       8      9   10   11    12       13      14     15
 		my @col = split(/\t/,$line); # qseqid sseqid pident qlen aln_length slen sstrand qstart qend sstart send qcovs mismatch gapopen evalue bitscore
		my $contig = $col[1];
		my $score = $col[15];
		my $strand = $col[6];
		my $pident = $col[2];
		if 	($pident > 95) {	
			$d{$contig}{'score'}  = $score;
			$d{$contig}{'strand'}  = $strand;
			$d{$contig}{'sstart'}  = $col[9];
			$d{$contig}{'send'}  = $col[10];
			$d{$contig}{'qid'} = $col[0];
			$d{$contig}{'slen'} = $col[5];
			if (exists ($d{$contig})) {
				my $currentscore = $d{$contig}{'score'};
				if ($currentscore > $score){
					next;
				} else {
					$d{$contig}{'score'} = $score;
					$d{$contig}{'strand'} = $strand;
					$d{$contig}{'sstart'}  = $col[9];
					$d{$contig}{'send'}  = $col[10];
					$d{$contig}{'qid'} = $col[0];
					$d{$contig}{'slen'} = $col[5];
				}
			}
		}
	}
	#my $noofdisread = keys (%d);
	return (\%d);
}	
sub compare_twoHOH {
	my ($genomloc,$TEhash,$ghash) = @_;
	my %TEhash = %$TEhash;
	
	foreach my $contigkey (sort {$TEhash{$b}->{'score'} <=> $TEhash{$a}->{'score'} } keys (%TEhash)) {#sort on values of HoH
		#print STDERR "key in TEH is $contigkey\n";
		my $TEstrand;
		my $genomestrand;
		if (exists (${$ghash}{$contigkey})) {
			#print STDERR " $contigkey exists in ghash\n";
			#%extractcontiginfo = %{ $TEhash{$contigkey} };# copying only next level of keys and values
			$extractcontiginfo{$contigkey} = $TEhash{$contigkey};
			foreach my $features (keys ${$TEhash}{$contigkey}) {#derefence second level
				$TEstrand = ${$TEhash}{$contigkey} {$features} if ($features eq 'strand');
				$genomestrand = ${$ghash}{$contigkey} {$features} if ($features eq 'strand');
				#print STDERR "TEstrand is $TEstrand and gstrand is $genomestrand \n";
			}
			if ((($TEstrand eq 'plus') && ($genomestrand eq 'minus')) || (($TEstrand eq 'minus') && ($genomestrand eq 'plus'))) {
				$strandprediction{$genomloc}{'strand'} = '-';	
			} elsif ((($TEstrand eq 'plus') && ($genomestrand eq 'plus')) || (($TEstrand eq 'minus') && ($genomestrand eq 'minus'))) {	
				$strandprediction{$genomloc}{'strand'} = '+';	
			}
			last;		
		} else {
			next;
		}
	}
	return (\%strandprediction,\%extractcontiginfo);
}
sub find_no_te_cand {
		foreach my $locus (sort keys %alldirectory) {
		if (exists ($strandprediction{$locus})) {
			next;
		} else {
			$strandprediction{$locus}{'fulllengTE'} = 'no';
			$strandprediction{$locus}{'strand'} = '?';
			$strandprediction{$locus}{'seq_id'} = 'noTE';
		}	
	}
	return (\%strandprediction);
}
sub print_keyvalue_HOH {
	my ($filepath,$hashtoprint) = @_;
	my %hashtoprint =%$hashtoprint;
	open (my $kp,">","$filepath") || die ("failed to open file to write $filepath $!\n");
	foreach my $mate (sort keys %hashtoprint) {
		print $kp "$mate\t$hashtoprint{$mate}{'fulllengTE'}\t$hashtoprint{$mate}{'strand'}\t$hashtoprint{$mate}{'seq_id'}\n";
	}
	close $kp;
}
sub printarray_fasta {
	my ($filetoprint,@seqarray) = @_;		
	my $seqarray_obj = Bio::SeqIO->new(-file   => ">$filetoprint",
										      -format => 'fasta') 
									         or die "\t    ERROR - Failed to create SeqIO FH object from $filetoprint  $!\n";  
	foreach my $seque (@seqarray) {
		$seqarray_obj->write_seq($seque);	
	}
}
sub loadseq_tohash {
	foreach  my $locus (keys %extractcontiginfo) {
		my ($dirinfo,$contigname) = split(/\./,$locus,2); 
		my $readio_obj = Bio::SeqIO->new(-file 	=> "$path/Assembled_TEreads/$dirinfo/Renamed_Assembledseq/$dirinfo.rename.fasta", 
									-format     => "fasta" 
									) or die "\t    ERROR - Failed to create SeqIO FH object from $outpath/Renamed_Assembledseq/$dirinfo.rename.fasta $!\n";  
		while (my $seque = $readio_obj->next_seq() ){
			my $header = $seque->display_id;
			my $sequence = $seque->seq;
			if ($header =~ /$locus/) {
				$extractcontiginfo{$locus}{'seq'} = $sequence;
			} else {
				next;
			}
		}
	}
}
sub generate_kmer_pos {
	my $kmersize = 10;
	my ($seq,$direction,$genomeloc) = @_;
	my $lenseq = length $seq;	
	#print "the length of the sequence $lenseq\n";
	#print "the length of the left right $lenleft\n";    
    my %kmers;    
    for (my $j = $kmersize; $j <= 24; $j++) {#index starts from 0
    	my $klen = $j;
		if ($direction == 3) {
			for (my $i = 0; ($i + $j) <= $lenseq; $i++) {
				$kmers{$genomeloc}{$klen}{$i} = substr($seq, $i, $j);
			}
		}
		if ($direction == 5) {# index starts from the kmerlength number
			for (my $i = $j; $i <= $lenseq; $i++) {#to extract kmers from the end to the beginning#check the lenght 
				#my $s = ($i+$j); 
				#print "i is $i\n";
				#print "j is $j\n";
				#print "s is $s\n";
				$kmers{$genomeloc}{$klen}{$i} = substr($seq, -$i, $j ); #This is the logic used #-4,4 #-5,4 #-6,4
			}
		}
	}
    return (\%kmers);
}
sub compare_kmers {
	my ($l, $r) = @_;
	my %lefthash = %$l;
	my %righthash = %$r;
	my %tsdlist;
	my @potentialtsdlist;
	for my $loc (keys %lefthash) {
		for my $klen ( sort { $a <=> $b } keys $lefthash{$loc} ) {#sort keys numerically
			#print "$klen is the first level key\n"; 
			if (exists ($righthash{$loc}{$klen})) {
				for my $lindex (sort keys $lefthash{$loc}{$klen}) {
					#print "$lindex is the second level key for secondhash\n"; 
					for my $rindex (sort keys $righthash{$loc}{$klen}) {
					     my $lefttsd = $lefthash{$loc}{$klen}{$lindex};
					     my $righttsd = $righthash{$loc}{$klen}{$rindex};
						 #if ( $lefthash{$loc}{$klen}{$lindex} eq $righthash{$loc}{$klen}{$rindex} ) {
						 if ( fc($lefttsd) eq fc($righttsd) ) {#to check equality irrespective of the case
							 #print "for $klen is the matching TSD  $lefthash{$loc}{$klen}{$lindex} is at $lindex, $righthash{$loc}{$klen}{$rindex} is at $rindex\n";
							 my $tsdinfo = "$loc.$klen.$lefthash{$loc}{$klen}{$lindex}.$lindex.$rindex";
							 push (@potentialtsdlist,$tsdinfo);
							 #print "$tsdinfo\n";
						 } #elsif (amatch ($lefthash{$g},["i 9%"],$righthash{$m}) ) {
							#print "the approximate matching TSD $lefthash{$g} is at $g, $righthash{$m} at  $m\n";
						 #}
					}			
				}
				#print "@potentialtsdlist\n";
				my $count = @potentialtsdlist;
				if ($count == 1) {
					push (@matchingtsdlist,$potentialtsdlist[0]);
					@matchingtsdlist = uniq @matchingtsdlist ;
					#print STDERR "@matchingtsdlist\n";
				} else {
					@potentialtsdlist = ();
				}		
			} else {
				print STDERR "not identified in righthash\n"
			}
		}
	}
	return (@matchingtsdlist);
}
sub findcord_bothtsds {
	my ($element) = @_;
	my $genomloc;
	my $tsd1start;
	my $tsd1end;
	my $tsd2start;
	my $tsd2end;
	my $tsdsuccess;
							#		0					1					 2  3   4   5    6		7          8  9	
	print "$element\n";#8:128356417-128356917.NODE_2_length_316_cov_59.472803.96.239.minus.15.AGCCAGGCAATTTTT.39.6
	my @details = split (/\./,$element);
	$genomloc = $details[0];
	my $te_start = $details[3];
	my $te_end = $details[4];
	my $tstrand = $details[5];
	my $tkmerlen = $details[6];
	my $tsd = $details[7];
	my $lefindex = $details[8];
	my $rigindex = $details[9];
	#my $gID = $genomloc.".".$details[1].".".$details[2];
	my $diff = ($lefindex - $tkmerlen);
	if ($diff < 10) {
		if ($tstrand eq 'plus') {
			#my $gchr = $tchr;
			$tsd1start = (($te_start+5) - $lefindex)+1;
			$tsd1end = ($tsd1start+$tkmerlen);
			$tsd2start = ($te_end-5) + $rigindex;
			$tsd2end = ($tsd2start + $tkmerlen);
			$tsdsuccess = "yes";
			#$gloc = $gchr.":".$tsd1start."-".$tsd1end."\t".$gchr.":".$tsd2start."-".$tsd2end;
			#$TEcor_bothtsd{$gloc} = $tsd;
			#$exactmatch{$gID} = 'exactmatch';
		} elsif ($tstrand eq 'minus') {#leftflank and rightflank is reverse complemented to identify the cordinates need to be modified
			#my $gchr = $tchr;
			
			$tsd1end = (($te_start+5) - $rigindex)+1;#  the sequence is reverse complemented so it is polyTend
			$tsd1start = ($tsd1end - $tkmerlen);
			$tsd2end = (($te_end-5) + $lefindex);
			$tsd2start = ($tsd2end - $tkmerlen)-1 ;#+1 and -1 is to get the TE cordinates after and before TSD
			$tsdsuccess = "yes";
			#$gloc = $gchr.":".$tsd1start."-".$tsd1end."\t".$gchr.":".$tsd2start."-".$tsd2end;
			#my $revtsd = &revcom_seq($tsd);      
			#$TEcor_bothtsd{$gloc} = $revtsd;
			#$exactmatch{$gID} = 'exactmatch';
		}
	} else {
		print STDERR "TSD with exact match cannot be identified for $element\n";
		$tsdsuccess = "no";
		$tsd1start = $te_start;
		$tsd2end  = $te_end;
	}
	return ($tsd1start,$tsd2end,$tsd,$tsdsuccess);
}
sub revcom_seq {
	my ($inputseq) = @_; 
	my $seqobj = Bio::Seq->new(-seq => "$inputseq", #to obtain reverse compliment of tsd when the strand is minus as I am extracting the reverse compliment of the sequence 
                         	   -alphabet => 'dna' );
	$seqobj = $seqobj->revcom; #reverse compliment the sequence
	my $revtsd = $seqobj->seq;
	return ($revtsd);
}