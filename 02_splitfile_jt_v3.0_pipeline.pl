#!/usr/bin/perl -w
#######################################################
# Author :  Jainy Thomas
# date   :  May 2017
# email  :  jainythomas1@gmail.com
# Pupose :  split a file to multiple files and make the list in another file
#####################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use File::Copy;
use File::chdir;
use File::Path qw(make_path remove_tree);
use Data::Dumper;
use Cwd;

my $version = "3.0";
my $scriptname = "splitfile_jt_v3.0.pl";
my $changelog = "
#   - v1.0 = 21 May 2017
#	- v2.0 = 17 June 2017
#  				added the option sort the files based on the individual name
#	- v3.0 = 29 August 2017
#				split based on the number of individuals is the only option
#
#
#
#
#
\n";
my $usage = "\nUsage [$version]: 
    perl $scriptname -f <file that needs to be split>  -n <number of individuals> [-s <print the sorted  file or not, give yes or no>] [-o name of outputfile with the list of files][-v] [-c] [-h] 
	
	MANDATORY ARGUMENT:	
 
    -f,--file  (string) file 
    -n, --number of individuals (STRING) split by the number of individuals	  
    
    OPTIONAL ARGUMENTS:
    
    -p, --path   (STRING)  path of the file                
    -s,--presort (STRING) print the sorted file based on the individual name before splitting files 
    -c,--chlog   (BOOL)   Print updates
    -v,--v       (BOOL)   Print version if only option
    -h,--help>   (BOOL)   Print this usage\n\n";
   
#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($file,$presort,$path,$indinumber,$out,$help,$v,$chlog);
GetOptions ('f=s' => \$file,
            's=s' => \$presort,
            'o=s' => \$out,
            'n=s' => \$indinumber,
            'p=s' => \$path,	
            'c'   => \$chlog, 
            'h'   => \$help,
            'v'   => \$v);

#check step to see if mandatory argument is provided + if help/changelog
die "\n Script $scriptname version $version\n\n" if ((! $file) && (! $indinumber)&&  (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ((! $file) ||  ($help));
my $cwd = getcwd();
$path = $cwd if (!$path) ;
$out = "$path/List_of_split_files.txt" if (! $out);
#$indinumber = 10 if (! $indinumber);
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my %sortindividual;
my $sortedfile = "$file.sorted.txt";;
my @listofsplitfiles;

&load_hash();

#print STDERR "the files is sorting now based on the individual name\n";
&sorthashbyvalue() if ($presort eq "yes") ;
print STDERR "Starting to split based on $indinumber individuals ....\n";
&split_byindivi($indinumber) if ($indinumber);
print STDERR "DONE..\n";
exit;
#-------------------------------------------------------------------------------
#----------------------------------- SUB ---------------------------------------
#-------------------------------------------------------------------------------

sub split_byindivi {
	
	my ($indinu) = shift;
	my $currentindividual;
	my %keyindividual;
	my $numofindi;
	my @elementstoprint;
	my $count =0;
	my $lastindivi;
	
	make_path ("$path/splitbyindividuals");
	foreach my $name (sort { $sortindividual{$a} cmp  $sortindividual{$b}} keys %sortindividual) {
		$currentindividual = $sortindividual{$name};
		
		push (@elementstoprint,$name);
		#print $fp "$name\n";
		$keyindividual{$currentindividual} = 1;
		$numofindi = keys (%keyindividual);
		
		if ($numofindi > "$indinu") {
			$count++;
			$lastindivi = pop (@elementstoprint);
			%keyindividual= ();
			&print_array($count,@elementstoprint);
			@elementstoprint =();
			push (@elementstoprint,$lastindivi);
		} else {
			next;
			print Dumper %keyindividual;
		}	
	}
	$count++;
	&print_array($count,@elementstoprint);
	&print_filename(@listofsplitfiles);
}

sub load_hash {
	open (my $fh,"<",$file) or die ("cannot open file $file to read $!\n");
	while (<$fh>) {
		chomp (my $data = $_);
		my @col = split (/\s+/,$data);
		#my @firstcol = split (/\./,$col[0]);
		my $individual = $col[0];
		$sortindividual{$data} = $individual;
	}
	return (%sortindividual);
}

sub sorthashbyvalue {
	open (my $hp,">","$file.sorted.txt") || die ("failed to open file to write (sub sorthashbyvalue) sorted file $!\n");
	foreach my $name (sort { $sortindividual{$a} cmp  $sortindividual{$b}} keys %sortindividual) {
		print $hp "$name\n";
	}
	close $hp;
}

sub print_array {
	my ($num,@array) = @_;
	
	open (my $ah,">","$path/splitbyindividuals/$file.$num.individuals.sorted.txt") || die ("cannot open file (sub print_array) to write $!\n");
	push (@listofsplitfiles,"$file.$num.individuals.sorted.txt");
		foreach my $dataline (@array) {
			print $ah "$dataline\n";
		}
	close ($ah);
}

sub print_filename {
	my @list = @_;
	my @sortedsplitfiles = sort (@list);
	#print "the array contains @sortedsplitfiles\n";	
	open (my $fhout,">",$out) || die ("cannot open file $out to write (sub print_filename) $!\n");
		foreach my $filename (@sortedsplitfiles) {
			print $fhout "$filename\n";
		}
	close $fhout;
}