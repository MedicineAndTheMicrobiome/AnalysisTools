#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_n $opt_d $opt_c);

getopts("i:nd:c");

my $usage = "
	usage:
	$0
		[-i <input file name>]
		[-n (number columns)]
		[-d <delimiter, default=tab>];
		[-c (count number of columns only)]

	This script will read in the first line of a delimited column
	text file (i.e. the header) and output the column names on
	multiple lines.

	-i: You can either use STDIN or specify an input file.
	-n: This will number the colums starting from 1.  For example, in 
		case you want you want to pull a column using cut -f.
	-d: You can change this to ,'s for example, for a csv file.
	-c: output the number of columns in the file

	Results goes to STDOUT.
	
";



if(!defined($opt_i) && (-t STDIN) && $ARGV[0] eq ""){
	die $usage;
}

my $File;
if($ARGV[0] eq ""){
	$File=$opt_i;
}else{
	$File=$ARGV[0];
}

my $Number=defined($opt_n);

my $Delim;
if(defined($opt_d)){
	$Delim=$opt_d;
}else{
	$Delim="\t";
}

my $NumCols=defined($opt_c);

###############################################################################

print STDERR "Running $0\n";

my $fh;

if(!(-t STDIN)){
	print STDERR "Reading from STDIN...\n";
	$_=<STDIN>;
}else{
	print STDERR "Opening $File.\n";
	open(IN_FH, "<$File") || die "Could not open $File\n";
	$_=<IN_FH>;
}

chomp $_;
my @header_arr=split /$Delim/, $_, -1;
my $num_cols=$#header_arr+1;

if($NumCols){
	print "$num_cols\n";
	exit;
}

for(my $i=0; $i<$num_cols; $i++){
	if($Number){
		print STDOUT ($i+1) . " ";
	}
	print STDOUT "$header_arr[$i]\n";		

}

print STDERR "Wrote $num_cols columns.\n";
print STDERR "Done.\n";
