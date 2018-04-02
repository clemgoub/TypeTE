#!/usr/bin/perl -w
#######################################################
# Author :  Jainy Thomas
# date   :  May 2017
# email  :  jainythomas1@gmail.com
# Pupose :  split a file to multiple files based on individuals or line number of bytes and make the list in another file 
#####################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use File::Copy;
use File::chdir;
use File::Path qw(make_path remove_tree);
use Data::Dumper;

my $version = "3.0";
my $scriptname = "splitfile_jt_v3.0.pl";
my $changelog = "
#   - v1.0 = 21 May 2017
#	- v2.0 = 17 June 2017
#  				added the option sort the files based on the individual name
#	- v3.0 = 29 August 2017
#				split based on the number of individuals is an option
#
#
#
#
#
\n";
my $usage = "\nUsage [$version]: 
    perl $scriptname -f <file that needs to be split> -s <sort the file or not, give yes or no> -l <line number> or-b <bytes> [-o name of outputfile with the list of files][-p name of the splitfiles] [-v] [-c] [-h] 
	
	MANDATORY ARGUMENT:	
 
    -f,--file  (string) file 
    -l,--line   (STRING) line number that needs to be split
    or
    -b,--bytes  (STRING) bytes 
    or 
    -n, --number of individuals (STRING) split by the number of individuals	  
    OPTIONAL ARGUMENTS:
    
    -pr,--prefix (STRING) name of the split files(aa,ab,ac will be added to the prefix) available only for lines/bytes splitting  
    -p, --path   (STRING)  path of the file                
    -s,--presort (STRING) sort the file based on the individual name before splitting files based on line or bytes also prints the sorted files
    -c,--chlog   (BOOL)   Print updates
    -v,--v       (BOOL)   Print version if only option
    -h,--help>   (BOOL)   Print this usage\n\n";
   
#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($file,$line,$bytes,$prefix,$presort,$path,$indinumber,$out,$help,$v,$chlog);
GetOptions ('f=s' => \$file,
            'l=s' => \$line,
            'b=s' => \$bytes,
           'pr=s' => \$prefix,
            's=s' => \$presort,
            'o=s' => \$out,
            'n=s' => \$indinumber,
            'p=s' => \$path,	
            'c'   => \$chlog, 
            'h'   => \$help,
            'v'   => \$v);

#check step to see if mandatory argument is provided + if help/changelog
die "\n Script $scriptname version $version\n\n" if ((! $file)&&  (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ((! $file) ||  ($help));
$prefix = "$file.sorted.split" if ((! $prefix) && ($presort));
$out = "$prefix.sorted.list.txt" if ((! $out) && ($presort)) ;
$prefix = "$file.split" if (! $prefix);
$out = "$file.split.list.txt" if (! $out);
$path = "./splitbyindividuals" if (! $path);
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my %sortindividual;
my $sortedfile;
my @listofsplitfiles;
if ($presort eq "yes") {
	print STDERR "the files is sorting now based on the file name\n";
	$sortedfile = "$file.sorted.txt";
	&load_hash();
	&sorthashbyvalue();
	&split_file ('l',$line,$sortedfile) if ($line);
	&split_file ('b',$bytes,$sortedfile) if ($bytes);
	print STDERR "Starting to split based on individuals....\n";
	&split_byindivi($indinumber) if ($indinumber);
}
if ($presort eq "no") {
	if ($line) {
		print STDERR "splitting files based on line number\n";
		&split_file('l',$line,$file) ;
	}
	if ($bytes) {
		print STDERR "splitting files based on bytes\n";
		&split_file('b',$bytes,$file) ;
	}
	
	if ($indinumber) {
		print STDERR "splitting based on individuals....\n";
		&load_hash();
		&split_byindivi($indinumber);
	}
}
print STDERR "DONE..\n";
exit;
#-------------------------------------------------------------------------------
#----------------------------------- SUB ---------------------------------------
#-------------------------------------------------------------------------------
sub split_file {
	my ($criteria,$value,$file1) = @_;
	my @splitfiles;
	my $splitfiledir = "splitfiles";
	mkdir $splitfiledir;
	copy ("$file1", "$splitfiledir/$file1") or die "Copy failed $!\n";
	{ #consider the current directory as following for the steps within the double brackets
		(local $CWD) = "$splitfiledir/";
			if ($criteria eq 'l') {
				  `split -l $value $file1 $prefix`;
			} elsif ($criteria eq 'b') {
				 `split -b $value $file1 $prefix`;
			}
	}
	unlink "$splitfiledir/$file1";
	@splitfiles = `ls $splitfiledir`; 
	my @sortedsplitfiles = sort (@splitfiles);
	print "the array contains @splitfiles\n";	
	open (my $fhout,">",$out) || die ("cannot open file $out to write $!\n");
		foreach my $filename (@sortedsplitfiles) {
			print $fhout "$filename";
		}
	close $fhout;
}

sub split_byindivi {
	
	my ($indinu) = shift;
	my $currentindividual;
	my %keyindividual;
	my $numofindi;
	my @elementstoprint;
	my $count =0;
	my $lastindivi;
	
	make_path ("$path") if ($path);
	foreach my $name (sort { $sortindividual{$a} cmp  $sortindividual{$b}} keys %sortindividual) {
		$currentindividual = $sortindividual{$name};
		
		push (@elementstoprint,$name);
		#print $fp "$name\n";
		$keyindividual{$currentindividual} = 1;
		$numofindi = keys (%keyindividual);
		
		if ($numofindi == "$indinu") {
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
	open (my $ah,">","$path/$file.$num.individuals.sorted.txt") || die ("cannot open file (sub print_array) to write $!\n");
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