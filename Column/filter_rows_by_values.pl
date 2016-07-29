#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_c $opt_v $opt_C $opt_V);

getopts("i:t:c:v:C:V:");

my $usage = "
	usage:
	$0
		-i <input file>
		
		Filter out below:
		[-c <columns, comma separated starting from 0>]
		[-v <cutoffs, comma separated, in same order as columns>]

		Filter out above:
		[-C <columns, comma separated starting from 0>]
		[-V <cutoffs, comma separated, in same order as columns>]


	Filters input files based on the values in the columns by the
	cutoffs.  

	The results are considered ANDs of all the criteria.

";

if(!defined($opt_i)){
	die $usage;
}

my $file=$opt_i;

my @below_col=split /,/, $opt_c;
my @above_col=split /,/, $opt_C;

my @below_val=split /,/, $opt_v;
my @above_val=split /,/, $opt_V;

my $num_below_col=$#below_col+1;
my $num_below_val=$#below_val+1;
my $num_above_col=$#above_col+1;
my $num_above_val=$#above_val+1;

if($num_below_col!=$num_below_val){
	die "Num below columns and cutoff values, don't match.\n";
}

if($num_above_col!=$num_above_val){
	die "Num above columns and cutoff values, don't match.\n";
}

print STDERR "Input File: $file\n";

###############################################################################

open(INPUT_FH, "<$file") || die "Could not open $file\n";
while(<INPUT_FH>){
	chomp;
	my @array=split /\t/, $_;
	

	my $filter=0;
	for(my $i=0; $i<$num_below_col; $i++){
		if($array[$below_col[$i]]<$below_val[$i]){
			$filter=1;	
		}
	}

	for(my $i=0; $i<$num_above_col; $i++){
		if($array[$above_col[$i]]>=$above_val[$i]){
			$filter=1;
		}
	}

	if(!$filter){
		print "$_\n";
	}

}
close(INPUT_FH);

###############################################################################
