#!/usr/bin/perl -w
#######################################################
# Author :  Jainy Thomas
# date   :  April 2018
# email  :  jainythomas1@gmail.com
# Pupose :  find mappability scores for genomic intervals
#           
#####################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Cwd;
use Data::Dumper;
use DBI;

my $version = "1.0";
my $scriptname = "findmappabilityscores_genomicintervals.pl";
my $changelog = "
#   - v1.0 = April 2018 
#				generates mappability scores specifically to insert in the deletion pipeline
#
#
#
#
#
#		
\n";

my $usage = "\nUsage [$version]: 
    perl $scriptname -t <mysqltable> -f <RMcordinate> -p <path of the outputdirectory> -db<mysql dbinfo> -u<user> -pd <password>  [-v] [-c] [-h] [-s] 
	
	
	
    MANDATORY ARGUMENT:	
    -t,--table 			(STRING) mysql table e.g.hg19wgEncodeCrgMapabilityAlign100mer_index
    -f,--file  			(STRING) file with locus and corresponding Repeatmasker info (10_105097522    10      105097539       105097848       AluYa5  SINE/Alu        +)
    
    -p,--path         	(STRING) output directory name (path)
                            	 Default = <current working directory>
    -db --mysqldbinfo 	(STRING) ex. command: DBI:mysql:database=jainys_db;host=localhost;port=22;     
    -u --user 			(STRING) Username for mysql database
    -pd,--password 		(STRING) password for mysql database
    
    OPTIONAL ARGUMENTS:
    -c,--chlog  		(BOOL)   Print updates
    -v,--v      		(BOOL)   Print version if only option
    -s,--verbose		(BOOL)   The script will talk to you
    -h,--help	  		(BOOL)   Print this usage\n\n";
   


#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($file,$path,$mysqltable,$mysqldb,$user,$password,$verbose,$help,$v,$chlog);
GetOptions ('f=s' => \$file,
            'p=s' => \$path,
            't=s' => \$mysqltable,
            'db=s'=> \$mysqldb,
            'u=s' => \$user,
            'pd=s'=> \$password,
            'c'   => \$chlog, 
            'h'   => \$help,
            's'   => \$verbose, 
            'v'   => \$v);

#check step to see if mandatory argument is provided + if help/changelog
die "\n Script $scriptname version $version\n\n" if ((! $mysqltable) && (! $file) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ((! $mysqltable) || (! $file) ||  ($help));
my $cwd = getcwd();
#die $usage if (($mappability eq "yes") && ((! $mysqldb) || (! $user) ||  ($password))) ;
$path = $cwd if (!$path) ;
my $mapscores = "$path/$file.mappabilityscores.txt";

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

#my $mysqltable = "hg19wgEncodeCrgMapabilityAlign100mer_index";
my %cordinatescore;
my %allindimappingscores;
print STDERR "connecting to mysql database..\n";
my $dbh;
my $dsn = "DBI:mysql:database=$mysqldb;host=localhost;port=22";
$dbh = DBI->connect($dsn, $user, $password,
					{'RaiseError' => 1});

open (my $fh, "<", $file) or confess "\n ERROR (main): could not open to read $file $!\n";
print STDERR "Fetching data from mysql database..\n";
	while(<$fh>) {
		chomp (my $line = $_);
		my @col = split(/\s+/,$line);
		my $uniqueid = $col[0];
		my $chr = $col[1] ;
		my $start = $col[2];
		my $end = $col[3];
		my $genomeloc = $chr.":".$start."-".$end;
		&extractgmscores($uniqueid,$genomeloc);
	}
close $fh;
print STDERR "Printing data from mysql database..\n";
#print Dumper %allindimappingscores,"\n";
make_path ("$path");
&print_AOH (%allindimappingscores);
# Disconnect from the database.
$dbh->disconnect();
#-----------------------------------------------------------------------------
#----------------------------------- SUB -------------------------------------
#-----------------------------------------------------------------------------		
		
sub extractgmscores {
	my ($uid,$location) = @_;
	my @mappingscores = ();
	my @gl = split /[:,-]/, $location;
	#print "@gl\n";
	#my $ind = $gl[0];
	my $chro = $gl[0];
	my $te_start = $gl[1];
	my $te_end = $gl[2];
	my $length = 250;
	my $mstart = $te_start - $length;
	my $mend = $te_end +$length;
	my $givencordinate = $chro.":".$mstart."-".$te_start;#250 flanks upstream of TE
	#print STDERR "$givencordinate\n";
	my $leftscore = &store_mapscores($givencordinate);
	print "left score $length is $leftscore\n" if ($verbose);	#my $left250score = $g
	#push (@mappingscores,$gmscore);
	$givencordinate = $chro.":".$te_end."-".$mend;#250 flanks downstream of TE
	my $rightscore = &store_mapscores($givencordinate);#output only score here
	print "right score $length is $rightscore\n" if ($verbose);
	my $avgscore = ($leftscore+$rightscore)/2;#
	my $gmscore = $length."=>".$avgscore;
	
	push (@mappingscores,$gmscore);
	$length = 500;
	$mstart = $te_start - $length;
	$mend = $te_end +$length;
	$givencordinate = $chro.":".$mstart."-".$te_start;#500 flanks upstream of TE
	$leftscore = &store_mapscores($givencordinate);
	print "left score $length is $leftscore\n" if ($verbose);	
	$givencordinate = $chro.":".$te_end."-".$mend;#50 flanks downstream of TE
	$rightscore = &store_mapscores($givencordinate);#output only score here
	print "right score $length is $rightscore\n" if ($verbose);
	$avgscore = ($leftscore+$rightscore)/2;#
	$gmscore = $length."=>".$avgscore;
	push (@mappingscores,$gmscore);
	$allindimappingscores{$uid} = [@mappingscores];
	#@mappingscores = uniq (@mappingscores);
	return (\%allindimappingscores);
}
sub store_mapscores {
	my ($gcordinate) = @_;
	#print STDERR "$gcordinate\n";
	my $mappingscore;
	if (exists ($cordinatescore{$gcordinate})) {
		my $value = $cordinatescore{$gcordinate};
		$mappingscore = $value;
		#push (@mappingscores,$value);
	} else {
		my ($chro,$start5,$end5) = split /[:,-]/, $gcordinate;
		#print STDERR "$chro,$start5,$end5\n";
		$mappingscore = &extractgmscores_mysql($chro,$start5,$end5);
		#push (@mappingscores,$mappingscore);
	}
	return $mappingscore;
}
sub extractgmscores_mysql {
	my ($rchr,$rstart,$rend) = @_;
	#my $length = $rend - $rstart;
	my $gposition = $rchr.":".$rstart."-".$rend;
	#print "$gposition\n";
	my @startdata =();
	my @enddata=();
	#my $len_avg;
	my $avg;
	#retrieve closest numbers from the table using a command 
	my $sth = $dbh->prepare("(SELECT * FROM $mysqltable WHERE chromo='$rchr' and start <= '$rstart' order by start desc limit 1) 
								union 
							(SELECT * from $mysqltable where chromo='$rchr' and end >= '$rend' order by start asc limit 1)") ;
	$sth->execute();
	my $line = 0;
	while (my $row = $sth->fetchrow_hashref()) {
		my $startd = $row->{'start'};
		push (@startdata,$startd);
		my $endd = $row->{'end'};
		push (@enddata,$endd);
		
		print "Found a row: chr = $row->{'chromo'},$row->{'start'},$row->{'end'},score = $row->{'score'}\n" if ($verbose);
		print "-"x50,"\n"if ($verbose);
		$line++;
	}
	$sth->finish();
	my $tstart = $startdata[0];
	my $tend = $enddata[0] if ($line == 1);
	$tend = $enddata[1] if ($line == 2);
	#retrieve whole raw from mysql table using perlDBI  
	$sth = $dbh->prepare("SELECT * FROM $mysqltable WHERE chromo='$rchr' AND start >= '$tstart' AND end <= '$tend'");
	$sth->execute();
	my $count = 0;
	my $sum = 0;
	
	while (my $row = $sth->fetchrow_hashref()) {
		my $score = $row->{'score'};
		print "Found a row: chr = $row->{'chromo'},$row->{'start'},$row->{'end'},score = $row->{'score'}\n" if ($verbose);
		print "*"x50,"\n"if ($verbose);
		$count++;
		$sum += $score;
	}
	$sth->finish();
	if ($count >= 1) {
		$avg = $sum/$count ;
		$avg = sprintf("%.3f",$avg);
		#$len_avg = $length."=>".$avg;
		#$len_avg = "$tstart-$tend => $avg";
		 
		$cordinatescore{$gposition} = $avg;
	}
	return ($avg);	
}
sub print_AOH {
	#my $fileAOH = @_;
	my %GenomMappAoH = @_;
	my $score;
	open (my $mp,">","$mapscores") || die ("failed to open file to write map scores $!\n");
	print $mp "Locus\tMappability_score_500bp_flank\tMappability_score_250bp_flank\n";
	foreach my $indilocus (sort keys %GenomMappAoH) {
		my $value500;
		my $value250;
		foreach $score (@{$GenomMappAoH{$indilocus}}) {
			my @chrs = split /=>/,$score;
			print "@chrs\n";
			$value500 = $chrs[1] if ($chrs[0] == 500);
			$value250 = $chrs[1] if ($chrs[0] == 250);		
		}
		print $mp "$indilocus\t$value500\t$value250\n";
		#print "$value500\t$value250\n";
	}
	close ($mp);
}