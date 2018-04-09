#!/usr/bin/perl -w
#######################################################
# Author 	:  Jainy Thomas
# date   	:  Jan 2017
# email  	:  jainythomas1@gmail.com
# Purpose 	:  to attach the name of the individual to the list of coordinates
#####################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use File::Path qw(make_path remove_tree);

my $version = "1.5";
my $scriptname = "makelist.pl";
my $changelog = "
#   - v1.0 = 17 July 2017
#	- v1.5 = 6 April 2018
#				added path of the outputdirectory
#
#
#
\n";
my $usage = "\nUsage [$version]: 
    perl $scriptname -t <BAM ID table> -f <file with TE breakpoints> [-p path of output directory] [-o name of outputfile] [-v] [-c] [-h] [-b]
	
    MANDATORY ARGUMENT:	
    -t,--table 	(STRING) file contain accession information first column needs to be the IDs, second column BAMIDs
    -f,--file  	(STRING) file with list of coordinates
    	  
    OPTIONAL ARGUMENTS:
  	-p --path 	(STRING) path of output directory 
    -o --output (STRING) name of the outputfile    
    -c,--chlog  (BOOL)   Print updates
    -v,--v      (BOOL)   Print version if only option
    -s,--verbose(BOOL)   The script will talk to you
    -h,--help>  (BOOL)   Print this usage\n\n";
   
#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($file,$table,$out,$path,$verbose,$help,$v,$chlog);
GetOptions ('t=s' => \$table,
            'f=s' => \$file,
            'p=s' => \$path,
            'o=s' => \$out,	
            'c'   => \$chlog, 
            'h'   => \$help,
            's'   => \$verbose, 
            'v'   => \$v);

#check step to see if mandatory argument is provided + if help/changelog
die "\n Script $scriptname version $version\n\n" if ((! $table) && (! $file) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ((! $file) || !($table) || ($help));
$out = "$path/$file.list.txt" if (! $out);
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my %individualids;
make_path ("$path");
&load_IDs();
&print_file();
exit;


#-------------------------------------------------------------------------------
#----------------------------------- SUB ---------------------------------------
#-------------------------------------------------------------------------------

sub load_IDs {
	open (my $th,"<",$table) or die ("cannot open file $table to read $!\n");
	while (<$th>) {
		chomp (my $data = $_);
		my @col = split (/\s+/,$data);
		#my @firstcol = split (/\./,$col[0]);
		my $individual = $col[0];
		$individualids{$individual} = 1;
	}
	return (%individualids);
	close ($th);
}

sub print_file {
	open (my $ih, ">",$out ) or die ("cannot write file $out $!\n");
	open (my $fh,"<",$file) or die ("cannot open file $file to read $!\n");
	while (<$fh>) {
		chomp (my $dataline = $_);
		#my @colum = split (/\s+/,$dataline);
		#my $line = join ("\t", @colum);
		foreach my $indi (sort keys %individualids) {
		 	print $ih "$indi\t$dataline\n";
		}
	}
	close ($ih);
	close ($fh);
}