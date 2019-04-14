#!/usr/bin/perl -w
#######################################################
# Author :  Jainy Thomas
# date   :  July 2017
# email  :  jainythomas1@gmail.com
# Pupose :  find target site duplication for deletion pipeline if the TEcordinates are provided
#           
#####################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Cwd;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Data::Dumper;
use String::Approx 'amatch';
use List::MoreUtils qw(uniq);
use feature 'fc';


my $version = "3.4";
my $scriptname = "findtsdinfasta.pl";
my $changelog = "
#   - v1.0 = 27 November 2017 
#	 
# 	- v2.0 = 1 December 2017
#				to do fuzzy matching for the ones no exact matches have identified
# 				to failed ones can do one thing, extract 10 bp more polyA and see if that help (RM coordinates are off at poly A region)	
#	- v2.1 = 17 January 2018
#				fixed a bug- TSDs were reported in the antisense orientation where TE insertion is minus
#  				output both TSDs and their start and stop positions
#				output file with approximate match TSDs
#				output file with no TSDs identified
#	- v2.2 = 23 January 2018
#				changed the length of bp that is extracted to 45 bp, 5bp internal of 5'end of TE and 40bp outside and 10bp internal of 3' end of TE and 35 bp outside RM cordinates of TE
#	- v3.0 = 10 April 2018			
#				verified and corrected the cordinates of the TSDs generated 
#	- v3.3 = 12 April 2018
#				changed input format for RM, so changed script accordingly
#				added TE locus to the output, also added approximate TSDs and no TSD data has the same cordinates as the RMout file
#				TSD is in proper orientation as that of the insertion
#				To do : the cordinates of  appmatch TSD will have to be adjusted, if I find more cases it would be easy especially + strand cases
#
#	-v3.4 = 14 November 2018
#				Previus version outputted all possible matches for approximate match. Now changed in such a way the first one with the lowest left index and lowest right index will be outputted
#				also changed the max. kmer to 26
#				compare_kmers_approx criteria changed from 3% to 5%

\n";

my $usage = "\nUsage [$version]: 
    perl $scriptname -t <table RM out in bed format> -f <file containing break points> [-p <path of the outputdirectory>][-o <output file>] [-v] [-c] [-h] 
	
	
	
    MANDATORY ARGUMENT:	
    -t,--table (string) repeatmasker output in bed format (locusname,chr, start,end,name,family,strand,)
    or 
    -f,--file  (string) file containing TE cordinates with flanking sequences annotated as left flank and right flank
    -g,--genome(string) change path in the script if not providing	  
    OPTIONAL ARGUMENTS:
    -p,--path   (STRING) output directory name (path)
                         Default = <current working directory>
    -o,--output (STRING) output file
    -c,--chlog  (BOOL)   Print updates
    -v,--v      (BOOL)   Print version if only option
    -s,--verbose(BOOL)   The script will talk to you
    -h,--help>  (BOOL)   Print this usage\n\n";
   


#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($file,$path,$out,$table,$GENOME,$verbose,$help,$v,$chlog);
GetOptions ('f=s' => \$file,
            'p=s' => \$path,
            'o=s' => \$out,
            't=s' => \$table,
            'g=s' => \$GENOME,
            'c'   => \$chlog, 
            'h'   => \$help,
            's'   => \$verbose, 
            'v'   => \$v);

#check step to see if mandatory argument is provided + if help/changelog
die "\n Script $scriptname version $version\n\n" if ((! $table)  &&  (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if  ($help);
my $cwd = getcwd();
$path = $cwd if (!$path) ;
$file = "$path/Genomicflankseq.extract.$version.fasta" if (! $file);# can be given in teh commandline if the file is already generated and would like to reuse
$GENOME = "/vbod2/jainy/SGDP/Project2/hg19.refFIX.fa" if (! $GENOME);
$out = "$path/tsdinfo.v.$version.txt" if (!$out) ;
my $appout = "$path/appmatch_tsdinfo.v.$version.txt";
my $notsdlist = "$path/no_tsds.v.$version.txt";

#my $logfile = "$file.log";
my $wholeout = "$path/TEcordinates_with_bothtsd_cordinates.v.$version.txt";#should be removed before doing reruns
#my $appwholeout = "$path/appmatch_TEcordinates_with_bothtsd_cordinates.v.$version.txt";
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

=pod
extract sequence 40 bp flanking the insertion using bioperl and rename RMTEcordnates_leftflank and RMTEcordinates_rightflank
always the flank of the plus insertion (include 3 bp inside of the element)
Step 1.load into hash key as header and array as the sequence
       or read the file and process it
Step 2 split the sequence before loading to array so each letter has an index
		Access the left flank and right flank 
step 3 generate kmers of for left flank and right flank
step 4 compare the hashes for kmers, do approximate check as well, regenerated kmers

=cut
#extracting the sequence based on the RM cordinates
  
make_path("$path") if ($path);
#my %gloci;
#extracted flanking genomic sequence
&extract_genomicseqflank() if ($table);
#print Dumper %gloci,"\n"; 
my %flankseq = load_seq ();
#print Dumper  %flankseq;
my $left = 5;
my $right = 3;
my @matchingtsdlist;
my %TEcor_bothtsd;
my %TEcor_notsd;
my %exactmatch;
my @approxtsdlist;
my %TEcor_bothtsd_approx;
my @likelynotsds;
my %apptsdfound;
#my %duplicatesappmatch;
#my @multipelTSDinfo;
my %approxtsdhash;

foreach my $locus (keys %flankseq) {
	my ($leftkmer);
	my ($rightkmer);
	
	foreach my $flanks (keys %{$flankseq{$locus}}) {
		#print "locus is $locus\n";
		my $leftseq = $flankseq{$locus}{$flanks} if ($flanks eq 'leftflank');		
		my $rightseq = $flankseq{$locus}{$flanks} if ($flanks eq 'rightflank');
		$leftkmer = &generate_kmer_pos($leftseq,$left,$locus) if defined ($leftseq);
		$rightkmer = &generate_kmer_pos($rightseq,$right,$locus) if defined ($rightseq);
		#print Dumper %$leftkmer;
		#print Dumper %$rightkmer;
	}
	&compare_kmers($leftkmer,$rightkmer);
}

&print_array($out,@matchingtsdlist);#print the perfect match TSD output#output is tsdinfo
#&findcord_afteronetsd(@matchingtsdlist);
&findcord_bothtsds(@matchingtsdlist);
&print_keyvalue($wholeout,%TEcor_bothtsd);
#print Dumper %exactmatch;
my %appromatch = &find_notsd_cand(%exactmatch);#takes the locus for which no exact matches has been identified or no proper match identified

#print "appromatch\n";

#print Dumper %appromatch;
foreach my $locus (keys %appromatch) {# same process is repeated and checked for approximate matches
	my ($leftkmer);
	my ($rightkmer);
	foreach my $flanks (keys %{$appromatch{$locus}}) {
		#print "locus is $locus\n";
		my $leftseq = $appromatch{$locus}{$flanks} if ($flanks eq 'leftflank');		
		my $rightseq = $appromatch{$locus}{$flanks} if ($flanks eq 'rightflank');
		$leftkmer = &generate_kmer_pos($leftseq,$left,$locus) if defined ($leftseq);
		$rightkmer = &generate_kmer_pos($rightseq,$right,$locus) if defined ($rightseq);
		#print Dumper %$leftkmer;
		#print Dumper %$rightkmer;
	}
	&compare_kmers_approx($leftkmer,$rightkmer);
}
print Dumper %approxtsdhash,"\n";
&print_array($appout,@approxtsdlist);
#&findcord_afteronetsd_approx(@approxtsdlist);
#&findcord_bothTSDs_approx(@approxtsdlist);
&findcord_bothTSDs_approxhASH(%approxtsdhash);
#&print_keyvalue($appwholeout,%TEcor_bothtsd_approx);
&print_keyvalue($wholeout,%TEcor_bothtsd_approx);
#print Dumper %duplicatesappmatch, "\n";
#my $uniqueapprox = &choose_bestapproxi(\%duplicatesappmatch);
#print Dumper %$uniqueapprox, "\n";
#finding the elements for which no approximate TSDs have been identified
my @noTSDswithapp = &findinarray(\@likelynotsds,\%apptsdfound);
print "@noTSDswithapp\n";
if (@noTSDswithapp) {
	print STDERR "failed to identify TSDs for some. Please check file $notsdlist\n";
	&print_array($notsdlist,@noTSDswithapp);
	&print_keyvalue($wholeout,%TEcor_notsd);
	
} else {
	print STDERR "putative TSDs have been identified for all cases\n";
}
#print Dumper %TEcor_notsd_approx;

#-----------------------------------------------------------------------------
#------------------------------ SUBROUTINES ----------------------------------
#-----------------------------------------------------------------------------

#if sequence is provided with leftflank and rightflank notations it is loaded into a hash

sub load_seq {
	my %fseq;
	my $readio_obj = Bio::SeqIO->new(-file 	=> "$file", 
								-format     => "fasta" 
								) or die "\t    ERROR - Failed to create SeqIO FH object from $file $!\n";  
	while (my $seq = $readio_obj->next_seq() ){
		my $header = $seq->display_id;
		my $sequence = $seq->seq;
		#my @seqarray = split (//,$sequence);
		if ($header =~ /(.*)\_leftflank$/) {
			my $head = "$1";
			#$fseq{$head}{'leftflank'} = [@seqarray];
			$fseq{$head}{'leftflank'} = $sequence;
		} elsif ($header =~ /(.*)\_rightflank$/) {
			my $head = "$1";
			#$fseq{$head}{'rightflank'} = [@seqarray];
			$fseq{$head}{'rightflank'} = $sequence;
		}
	}
	return %fseq;
}
#kmers are generated 
sub generate_kmer_pos {
	my $kmersize = 3;
	my ($seq,$direction,$genomeloc) = @_;
	my $lenseq = length $seq;	
	#print "the length of the sequence $lenseq\n";
	#print "the length of the left right $lenleft\n";    
    my %kmers;    
    for (my $j = $kmersize; $j <= 26; $j++) {#index starts from 0
    	my $klen = $j;
		if ($direction == 3) {
			for (my $i = 0; $i + $j <= $lenseq; $i++) {
				$kmers{$genomeloc}{$klen}{$i} = substr($seq, $i, $j);
			}
		}
		if ($direction == 5) {# index starts from the kmerlength number
			for (my $i = $j; $i <= $lenseq; $i++) {#to extract kmers from the end to the beginning
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
		for my $klen ( sort { $a <=> $b } keys %{$lefthash{$loc}} ) {#sort keys numerically
			#print "$klen is the first level key\n"; 
			if (exists ($righthash{$loc}{$klen})) {
				for my $lindex (sort keys %{$lefthash{$loc}{$klen}}) {
					#print "$lindex is the second level key for secondhash\n"; 
					for my $rindex (sort keys %{$righthash{$loc}{$klen}}) {
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
sub extract_genomicseqflank {
	# index the genome and connect to the fasta file
	my $reindex;
	my $indexfile = "$GENOME.index";
		if (-e $indexfile) {
			$reindex = 0;
			print "\t Genome previously indexed - Skipping indexing...\n\t (if you want to reindex the genome, delete $indexfile)\n";
		} else {
			$reindex = 1;
			print "\t  Genome not indexed - Indexing $GENOME...\n";
		}
	my $db = Bio::DB::Fasta->new( $GENOME, -reindex=>$reindex) or die "\t    ERROR - Failed to create Bio::DB::Fasta object from $GENOME $!\n";
	#create list of the ID of the genome file  
	my @dbIDs = $db->get_all_ids();
	#make_path  ("$path/ExtractGenomicFlanks");
	my $extractgenomicseqout = "$path/Genomicflankseq.extract.$version.fasta";
	
	my $gout = Bio::SeqIO->new(-format => 'Fasta',
								-file   => ">>$extractgenomicseqout") 
								or 	die "\t    ERROR - Failed to create SeqIO FH object from $extractgenomicseqout $!\n";
	my $log = "$extractgenomicseqout.log";
	open (LOG, ">>$log") or die "\t    ERROR - can not create log file $log $!\n";
	open (my $th, "<", $table) or confess "\n ERROR (main): could not open to read $table $!\n";
		while (my $data = <$th>) { #reading the table
			chomp $data; #storing the table values to the $data and removing the extraline
			my @col = split(/\t/,$data); # splitting the data based on tab and storing into the arrray
			my $ulocus = $col[0];
			my $chr = $col[1];
			my $start = $col[2];
			my $end = $col[3];
			my $strand = $col[6];
			#$gloci{$ulocus} =1;
			if ($strand eq "+") {
				my $leftstart = $start - 40;
				my $leftend = $start + 5;
				my $rightend = $end + 35;
				my $rightstart = $end - 10;
				my $uniqueid = $ulocus.".".$chr.".".$leftend.".".$rightstart.".".$strand;#chr, same start as the element RM start, right start is 5 bp less than the RM end
				my $newLId = join('_',$uniqueid,"leftflank");
				my 	$newRId = join('_',$uniqueid,"rightflank");
				if ("@dbIDs" =~ m/(\S*)($chr)(\S*)/) {
					my $subSeq = $db->seq($chr,$leftstart,$leftend);	# extract target sequence
					my $seqobj = Bio::Seq->new( -display_id => $newLId,
													   -seq => $subSeq); # create object with target
					$gout->write_seq($seqobj);
					#print $gout $seqobj;	# print it out (in fasta format)
					$subSeq = $db->seq($chr,$rightstart,$rightend);	# extract target sequence
					$seqobj = Bio::Seq->new( -display_id => $newRId,
													   -seq => $subSeq); # create object with target
					$gout->write_seq($seqobj);
					#print $gout $seqobj;	# print it out (in fasta format)
				
				} else {
					print LOG "$chr was not found in $GENOME\n";
				}
			} elsif ($strand eq "C") {
				my $leftend = $end + 40;
				my $leftstart = $end - 5;
				my $rightstart = $start - 35;
				my $rightend = $start + 10;
				my $uniqueid = $ulocus.".".$chr.".".$rightend.".".$leftstart.".".$strand;# chr, RM end +5 on the polyT side, RM start of the element in the minus orienation
				my $newLId = join('_',$uniqueid,"leftflank");
				my 	$newRId = join('_',$uniqueid,"rightflank");
				if ("@dbIDs" =~ m/(\S*)($chr)(\S*)/) {# the extracted sequence is reverse complimented when strand is C but the id is modified slight compared to the RM cordinates
					my $subSeq = $db->seq($chr,$leftstart,$leftend);	# extracted target sequence 5' region of the TE in the antisense orientation
					my $seqobj = Bio::Seq->new( -display_id => $newLId,
													   -seq => $subSeq); # create object with target
					$seqobj = $seqobj->revcom;
					$gout->write_seq($seqobj);
					#print $gout $seqobj;	# print it out (in fasta format)
					$subSeq = $db->seq($chr,$rightstart,$rightend);	# extracted target sequence 3' region of the TE in the antisense orientation
					$seqobj = Bio::Seq->new( -display_id => $newRId,
													-seq => $subSeq); # create object with target
					$seqobj = $seqobj->revcom;
					$gout->write_seq($seqobj);
					#print $gout $seqobj;	# print it out (in fasta format)
				} else {
					print LOG "$chr was not found in $GENOME\n";
				}
			}
		}
	close LOG;
	close $th;
}
#also finds the cases that does not have  exact TSD match to another hash for running approx match
sub findcord_bothtsds {
	my (@arraywithtsdinfo) = @_;
	my $gloc;
	foreach my $element (@arraywithtsdinfo) {
		#print "$element\n";
		my @details = split (/\./,$element);
		my $uniqlocus = $details[0];
		my $tchr = $details[1];
		my $tstart = $details[2];
		my $tend = $details[3];
		my $tstrand = $details[4];
		my $tkmerlen = $details[5];
		my $tsd = $details[6];
		my $lefindex = $details[7];
		my $rigindex = $details[8];
		my $gID = $uniqlocus.".".$tchr.".".$tstart.".".$tend.".".$tstrand;
		my $diff = ($lefindex - $tkmerlen);
		if ($diff < 10) {
			if ($tstrand eq '+') {
				my $gchr = $tchr;
				my $tsd1start = ($tstart - $lefindex)+1;
				my $tsd1end = ($tsd1start + $tkmerlen)-1;
				my $tsd2start = $tend + $rigindex;
				my $tsd2end = ($tsd2start + $tkmerlen) - 1;
				$gloc = $uniqlocus."\t".$gchr.":".$tsd1start."-".$tsd1end."\t".$gchr.":".$tsd2start."-".$tsd2end;
				$TEcor_bothtsd{$gloc} = $tsd;
				$exactmatch{$gID} = 'exactmatch';
			} elsif ($tstrand eq 'C') {
				my $gchr = $tchr;
				my $tsd1end = $tstart - $rigindex;
				my $tsd1start = ($tsd1end - $tkmerlen) + 1;
				#  because the tsd is identifed on reverse complimented sequence
				my $tsd2end = ($tend + $lefindex) - 1;
				my $tsd2start = ($tsd2end - $tkmerlen)+1;
				$gloc = $uniqlocus."\t".$gchr.":".$tsd1start."-".$tsd1end."\t".$gchr.":".$tsd2start."-".$tsd2end;
				my $revtsd = &revcom_seq($tsd);  #reverse compliment of the tsd is printed    
				$TEcor_bothtsd{$gloc} = $revtsd;
				$exactmatch{$gID} = 'exactmatch';
			}
		} else {
			print STDERR "TSD with exact match cannot be identified for $gloc\n";
		}
	}
	return (\%TEcor_bothtsd,\%exactmatch);
}
sub print_array {
	my ($filetoprint,@array)= @_;
	my @sortedarray = sort @array;
	open (my $ah,">","$filetoprint") || die ("cannot open file (sub print_array) $filetoprint to write $!\n");
		foreach my $dataline (@sortedarray) {
			print $ah "$dataline\n";
		}
	close ($ah);
}
sub print_keyvalue {
	my ($filepath,%hashtoprint) = @_;
	open (my $kp,">>","$filepath") || die ("failed to open file to write $filepath $!\n");
	foreach my $mate (sort keys %hashtoprint) {
		print $kp "$mate\t$hashtoprint{$mate}\n";
	}
	close $kp;
}
sub find_notsd_cand {
	my (%exactcandi) = @_;
	my %approcandi;
	foreach my $locus (sort keys %flankseq) {
		#my $nu = keys %flankseq;
		#my $c = keys %exactcandi;
		#print "$nu and $c are number of keys in two hash\n";
		if (exists ($exactcandi{$locus})) {
			next;
		} else {
			$approcandi{$locus} = $flankseq{$locus};
			#print STDERR "locus is $locus and value is $flankseq{$locus}\n";
		}
	}
	return (%approcandi);
}
sub compare_kmers_approx {
	my ($l, $r) = @_;
	my %lefthash = %$l;
	my %righthash = %$r;
	my @approx_poten_tsdlist;
	
	for my $loc (keys %lefthash) {
	my @listapprox=();
		print "the locus is $loc\n";
		$approxtsdhash{$loc}=[@listapprox];
		for my $klen ( sort { $b <=> $a } keys %{$lefthash{$loc}} ) {#sort keys numerically in the descending order
			#print "$klen is the first level key\n"; 
			if (exists ($righthash{$loc}{$klen})) {
				for my $lindex (sort keys %{$lefthash{$loc}{$klen}}) {
					#print "$lindex is the second level key for secondhash\n"; 
					for my $rindex (sort keys %{$righthash{$loc}{$klen}}) {
					     my $lefttsd = $lefthash{$loc}{$klen}{$lindex};
					     my $righttsd = $righthash{$loc}{$klen}{$rindex};
						 #if ( $lefthash{$loc}{$klen}{$lindex} eq $righthash{$loc}{$klen}{$rindex} ) {
						 #if ( fc($lefttsd) eq fc($righttsd) ) {
						 if (amatch ($lefthash{$loc}{$klen}{$lindex},["i 3%"],$righthash{$loc}{$klen}{$rindex}) ) {
							 #print "for $klen is the approximate matching TSD  $lefthash{$loc}{$klen}{$lindex} is at $lindex, $righthash{$loc}{$klen}{$rindex} is at $rindex\n";
							 my $tsdinfo = "$loc.$klen.$lefthash{$loc}{$klen}{$lindex}.$lindex.$rindex.$righthash{$loc}{$klen}{$rindex}";
							 
							 push (@approx_poten_tsdlist,$tsdinfo);
							 
							 #print "$tsdinfo\n";
						 } #elsif (amatch ($lefthash{$g},["i 9%"],$righthash{$m}) ) {
							#print "the approximate matching TSD $lefthash{$g} is at $g, $righthash{$m} at  $m\n";
						 #}
					}			
				}
				#print "@potentialtsdlist\n";
				my $count = @approx_poten_tsdlist;
				if ($count == 0) {
					next;
				} else {
					push (@approxtsdlist,@approx_poten_tsdlist);
					push @{$approxtsdhash{$loc}},@approx_poten_tsdlist;
					@{$approxtsdhash{$loc}} = uniq @{$approxtsdhash{$loc}};
					@approxtsdlist = uniq @approxtsdlist ;
					
					#print STDERR "@approxtsdlist\n";
					@approx_poten_tsdlist = ();
					last;
				}		
			} else {
				print STDERR "not identified in righthash\n"
			}
		}
	}
	return (@approxtsdlist);
}
sub findcord_bothTSDs_approx {
	my (@arraywithtsdinfo) = @_;
	my $gloca;
	
	foreach my $element (@arraywithtsdinfo) {
		#print "$element\n";
		my @details = split (/\./,$element);
		my $uniqlocus = $details[0];
		#$duplicatesappmatch{$uniqlocus} = [@multipelTSDinfo];
		my $tchr = $details[1];
		my $tstart = $details[2];
		my $tend = $details[3];
		my $tstrand = $details[4];
		my $tkmerlen = $details[5];
		my $ltsd = $details[6];
		my $lefindex = $details[7];
		my $rigindex = $details[8];
		my $rtsd = $details[9];
		my $gID = $uniqlocus.".".$tchr.".".$tstart.".".$tend.".".$tstrand;
		my $diff = ($lefindex - $tkmerlen);
		if ($diff < 10) {
			
			if ($tstrand eq '+') {
				my $gchr = $tchr;
				my $tsd1start = ($tstart - $lefindex)+1;
				my $tsd1end = ($tsd1start + $tkmerlen)-1;
				my $tsd2start = $tend + $rigindex;
				my $tsd2end = ($tsd2start + $tkmerlen)-1;
				$gloca = $uniqlocus."\t".$gchr.":".$tsd1start."-".$tsd1end."\t".$gchr.":".$tsd2start."-".$tsd2end;
				#print "$gloca\n";
				$TEcor_bothtsd_approx{$gloca} = $ltsd.".".$rtsd;
				#my $eledetails= $gloca ."=".$ltsd.".".$rtsd;
				$apptsdfound{$gID} = 1;
				#push (@{$duplicatesappmatch{$uniqlocus}},{$element => $eledetails});
				#push (@multipelTSDinfo, {$element => $eledetails}) ;
			} elsif ($tstrand eq 'C') {
				my $gchr = $tchr;
				my $tsd1end = $tstart - $rigindex;
				my $tsd1start = ($tsd1end - $tkmerlen) + 1;#  because the tsd is identifed on reverse complimented sequence
				
				my $tsd2end = ($tend + $lefindex) - 1;
				my $tsd2start = ($tsd2end - $tkmerlen)+1;
				$gloca = $uniqlocus."\t".$gchr.":".$tsd1start."-".$tsd1end."\t".$gchr.":".$tsd2start."-".$tsd2end;
				my $revltsd = &revcom_seq($ltsd); 
				my $revrtsd = &revcom_seq($rtsd);     
				$TEcor_bothtsd_approx{$gloca} = $revrtsd.".".$revltsd;
				#my $eledetails= $gloca ."=".$ltsd.".".$rtsd;
				$apptsdfound{$gID} = 1;
				#push (@{$duplicatesappmatch{$uniqlocus}},{$element => $eledetails});
				#push (@multipelTSDinfo, {$element => $eledetails}) ;
			}
		} else {
			#print STDERR "TSD with approximate match cannot be identified for $element\n";
			push (@likelynotsds,$gID);
			@likelynotsds = uniq @likelynotsds;
		}
	}
	return (\%TEcor_bothtsd_approx,\%apptsdfound,\@likelynotsds);
}
sub findcord_bothTSDs_approxhASH {
	my (%arraywithtsdinfo) = @_;
	my $gloca;
	
	foreach my $element (%arraywithtsdinfo) {
		foreach my $hash_ref (@{$arraywithtsdinfo{$element}}) {
			
			#print "$element\n";
			my @details = split (/\./,$hash_ref);
			my $uniqlocus = $details[0];
			my $tchr = $details[1];
			my $tstart = $details[2];
			my $tend = $details[3];
			my $tstrand = $details[4];
			my $tkmerlen = $details[5];
			my $ltsd = $details[6];
			my $lefindex = $details[7];
			my $rigindex = $details[8];
			my $rtsd = $details[9];
			my $gID = $uniqlocus.".".$tchr.".".$tstart.".".$tend.".".$tstrand;
			my $diff = ($lefindex - $tkmerlen);
			if ($diff < 10) {			
				if ($tstrand eq '+') {
					my $gchr = $tchr;
					my $tsd1start = ($tstart - $lefindex)+1;
					my $tsd1end = ($tsd1start + $tkmerlen)-1;
					my $tsd2start = $tend + $rigindex;
					my $tsd2end = ($tsd2start + $tkmerlen)-1;
					$gloca = $uniqlocus."\t".$gchr.":".$tsd1start."-".$tsd1end."\t".$gchr.":".$tsd2start."-".$tsd2end;
					#print "$gloca\n";
					$TEcor_bothtsd_approx{$gloca} = $ltsd.".".$rtsd;
					#my $eledetails= $gloca ."=".$ltsd.".".$rtsd;
					$apptsdfound{$gID} = 1;
					#push (@{$duplicatesappmatch{$uniqlocus}},{$element => $eledetails});
					#push (@multipelTSDinfo, {$element => $eledetails}) ;
				} elsif ($tstrand eq 'C') {
					my $gchr = $tchr;
					my $tsd1end = $tstart - $rigindex;
					my $tsd1start = ($tsd1end - $tkmerlen) + 1;#  because the tsd is identifed on reverse complimented sequence
				
					my $tsd2end = ($tend + $lefindex) - 1;
					my $tsd2start = ($tsd2end - $tkmerlen)+1;
					$gloca = $uniqlocus."\t".$gchr.":".$tsd1start."-".$tsd1end."\t".$gchr.":".$tsd2start."-".$tsd2end;
					my $revltsd = &revcom_seq($ltsd); 
					my $revrtsd = &revcom_seq($rtsd);     
					$TEcor_bothtsd_approx{$gloca} = $revrtsd.".".$revltsd;
					#my $eledetails= $gloca ."=".$ltsd.".".$rtsd;
					$apptsdfound{$gID} = 1;
					#push (@{$duplicatesappmatch{$uniqlocus}},{$element => $eledetails});
					#push (@multipelTSDinfo, {$element => $eledetails}) ;
				}
				last;
			} else {
				#print STDERR "TSD with approximate match cannot be identified for $element\n";
				push (@likelynotsds,$gID);
				@likelynotsds = uniq @likelynotsds;
			}
		}
	}
	return (\%TEcor_bothtsd_approx,\%apptsdfound,\@likelynotsds);
}
sub revcom_seq {
	my ($inputseq) = @_; 
	my $seqobj = Bio::Seq->new(-seq => "$inputseq", #to obtain reverse compliment of tsd when the strand is minus as I am extracting the reverse compliment of the sequence 
                         	   -alphabet => 'dna' );
	$seqobj = $seqobj->revcom; #reverse compliment the sequence
	my $revtsd = $seqobj->seq;
	return ($revtsd);
}
sub findinarray {
	my ($array,$hash) = @_;
	my %hash2 = %$hash;
	#print Dumper %hash2;
	#print "@likelynotsds\n"; 
	my @uniquevalues;
	foreach my $data (@$array) {
		chomp $data;		
		if (exists ($hash2{$data})) {
			print "$data\n";
			next;	
		} else {
			#print "$data\n";
			push (@uniquevalues,$data);
			my @position = split /\./,$data;#$gID = $uniqlocus.".".$tchr.".".$tstart.".".$tend.".".$tstrand;
			my $start_te = ($position[2]-5) if ($position[4] eq "+");
			my $end_te = ($position[3] + 10) if ($position[4] eq "+");
			$end_te = ($position[3]+5) if ($position[4] eq "C");#5 should be added to the 3'end
			$start_te = ($position[2] - 10) if ($position[4] eq "C");
			my $glocat = $position[0]."\t".$position[1].":".$start_te."-".$end_te."\t".$position[1].":".$start_te."-".$end_te;
			$TEcor_notsd{$glocat} = "noTSDs";
		}
	}
	#print "@uniquevalues";
	return @uniquevalues;
}
