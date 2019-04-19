#!/usr/bin/perl -w
#######################################################
# Author :  Jainy Thomas
# date   :  July 2017
# email  :  jainythomas1@gmail.com
# Pupose :  find corresponding TE insertions from RM output data in bed format if you are given a breakpoint
#           
#####################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Cwd;
use Data::Dumper;


my $version = "3.0";
my $scriptname = "findcorrespondingTEinsertions.pl";
my $changelog = "
#   - v1.0 = 30 July 2017 
#	- v2.0 = 4 August 2017
#				print the corresponding breakpoints as well in another file
#	- v3.0 = 12 April 2018
#				-not done that-> check the breakpoint with the end of the RM cordinate as in addition to the start
#				only output is with unique locus name
# 				the script is modified to find only Alu
#	- v3.1 = 10 January 2019
#				introduced another round of finding RM cordinates for which breakpoints were more 100 bp far, the new cut off is 220bp.
#	-v3.2-LINE1 = 19 April 2019			 
#				modified for LINE1 only
#				
#
\n";

my $usage = "\nUsage [$version]: 
    perl $scriptname -t <table RM out in bed format> -f <file containing break points> [-p <path of the outputdirectory>][-o <output file>] [-v] [-c] [-h] 
	
	
	
    MANDATORY ARGUMENT:	
    -t,--table (string) repeatmasker output in bed format (chr, start,end,name,length,strand,)
    -f,--file  (string) file containing break points (locus,chr,brkpoint)    	  
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
my ($file,$path,$rawout,$table,$verbose,$help,$v,$chlog);
GetOptions ('f=s' => \$file,
            'p=s' => \$path,
            'o=s' => \$rawout,
            't=s' => \$table,
            'c'   => \$chlog, 
            'h'   => \$help,
            's'   => \$verbose, 
            'v'   => \$v);

#check step to see if mandatory argument is provided + if help/changelog
die "\n Script $scriptname version $version\n\n" if ((! $table) && (! $file) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ((! $table) || (! $file) ||  ($help));
my $cwd = getcwd();
$path = $cwd if (!$path) ;
$rawout = "$path/file.correspondingRepeatMaskerTEs.txt" if (!$rawout);
my $logfile = "$path/fileRMfailed.log";
my $loginter = "$path/fileRMfailed.inter.log";
#my $wholeout = "$file.withbrkpointsRMcordinates.txt";
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
#my $bamlocation = "/vbod2/cgoubert/Correct_Genotypes/1KGP_bams";#vader server

my %repeatmaskerfile = ();
#my $chro;
#my $brkpont;
#my $individual;
my @valuestoprint;
my @logvalues;
my @valueswithbrkpnt;
my $found;
my $rip;
my %brkpointcorre;

make_path("$path") if ($path);
print STDERR "loading the Repeatmasker data....\n";
%repeatmaskerfile = load_file ($table);
#print Dumper %repeatmaskerfile;
print STDERR "finding the corresponding elements....\n";
open (my $fh, "<", $file) or confess "\n ERROR (main): could not open to read $file $!\n";
	while(<$fh>) {
		chomp (my $line = $_);
		my @col = split(/\s+/,$line);
		#$individual = $col[0];
		my $chro = $col[1] ;
		my $brkpont = $col[2];
		$rip = $col[0];
		#$found = $rip;
		push (@logvalues,$rip);
		my $num = 50;
		&find_element($num,$chro,$brkpont,$rip);
	}
#&print_array($rawout,@valuestoprint);
#&print_array($loginter,@logvalues)
print (@logvalues);
if (@logvalues) {
	my (@remvalues) = @logvalues ;
	my @final_log;
	my $lastindex = $#remvalues;
	#print STDERR "$lastindex is the last index\n";
	print STDERR "now checking the ones that are failed in the first round\n";
	for(my $i = $lastindex; $i >= 0 ; $i--)  {
		#print "the element is $remvalues[$i]\n";
		my $element =$remvalues[$i];
		my ($chromo,$breakp) = split(/\_/,$element);
		#print "$chromo,$breakp chromo and breakpoint from logvalues\n";
		my $num = 110;
		my $present = &find_element($num,$chromo,$breakp,$element);
		push (@final_log, $element) unless ($present == 2);
	}
	print STDERR "Please check logfile for values for which Reference element not identified...\n";
	&print_array($logfile,@final_log);
}
&print_array($rawout,@valueswithbrkpnt);
exit;
#-----------------------------------------------------------------------------
#----------------------------------- SUB -------------------------------------
#-----------------------------------------------------------------------------
sub load_file {
	my ($file1) = @_;
	my %sbamfile;
	open (my $th, "<", $file1) or confess "\n ERROR (main): could not open to read $file1 $!\n";
		while (my $data = <$th>) { #reading the table
			chomp $data; #storing the table values to the $data and removing the extraline
			my @namebam = split(/\s+/,$data); # splitting the data based on tab and storing into the arrray
			my $chr = $namebam[0];
			$sbamfile{$chr}{$data} = [($namebam[0],$namebam[1],$namebam[2],$namebam[3])]; #loading to the hash
		}
	return (%sbamfile);
	close $th;	
}

sub find_element {
	my ($number,$chrom,$brekpont,$marker) = @_;
	#print "$number,$chrom,$brekpont are number chromo and brekpoint\n";
	my $found = 1;
	foreach my $chromosome (sort keys %repeatmaskerfile) {
		my $chromo;
		my $start;
		my $startminus;
		my $startplus;
		my $te;
		my $len;
		
		
		if ($chrom eq $chromosome) {
			foreach my $whline (sort keys %{ $repeatmaskerfile{$chromosome}}) {#sort the next level of keys
				#print "$whline\n";
				$chromo = $repeatmaskerfile{$chromosome}{$whline}->[0];
				$start = $repeatmaskerfile{$chromosome}{$whline}->[1];
				$end = $repeatmaskerfile{$chromosome}{$whline}->[2];
				$te = $repeatmaskerfile{$chromosome}{$whline}->[3];
				$len = $end-$start;
				$startminus = ($start - $number);
				$startplus = ($start + $number);
				if (($brekpont >= $startminus) && ($brekpont <= $startplus) && ($te =~ /LINE1.*/)) {
					print STDERR "$marker found!\n";
					#print STDERR "$whline\n";
					my $withbrkp = "$marker\t$whline"; 
					#push (@valuestoprint,$whline);
					if (exists $brkpointcorre{$marker}) {
						my $valueline= $brkpointcorre{$marker};
						my $colums = split(/\t/,$valueline);
						my $telen= ($columns[2]-$columns[1])+1;
						if ($len > $telen) {
							$brkpointcorre{$marker}=$whline;
						}
					} else {
						$brkpointcorre{$marker}=$whline;
					}
					
					push (@valueswithbrkpnt,$withbrkp);				
					pop (@logvalues);
					$found = 2;
				} 
			}
		} else {
			next;
		}
	}
	return ($found);
}

sub print_array {
	my ($outfile, @array)= @_;
	my @sortedarray = sort @array;
	open (my $fhout,">",$outfile) || die ("cannot open file $outfile to write $!\n");
		foreach my $dataline (@sortedarray) {
			print $fhout "$dataline\n";
		}
	close ($fhout);
}

