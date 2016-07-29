#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_m $opt_k $opt_s $opt_o);

getopts("i:m:k:s:o:");

my $usage = "
	usage:
	$0
		-i <input file>
		-m <map file name>
		-k <key column in input file, starting from 0>
		-s <insert columns in input file, starting from 0>
		-o <output file>

	This script will merge the input file with the map file.
	The map file is read in in its entirety into a hash
	so that it doesn't need to be sorted.

	-k specifies which column in the input file to use as the key.
	-s specifies which column to insert what the key maps to in the output.

";

if(
	!defined($opt_i) || 
	!defined($opt_m) || 
	!defined($opt_k) || 
	!defined($opt_s) || 
	!defined($opt_o) 
){
	die $usage;
}

my $InputFile=$opt_i;
my $MapFile=$opt_m;
my $KeyCol=$opt_k;
my $InsertCol=$opt_s;
my $OutputFile=$opt_o;

###############################################################################

my %map_hash;

open(MAP_FH, "<$MapFile") || die "Could not open map file: $MapFile\n";
while(<MAP_FH>){
	chomp;
	my @array=split /\t/, $_;
	$map_hash{$array[0]}=$array[1];
}
close(KEEP_FH);

###############################################################################

open(INPUT_FH, "<$InputFile") || die "Could not open input file: $InputFile\n";
open(OUTPUT_FH, ">$OutputFile") || die "Could not open output file: $OutputFile\n";

while(<INPUT_FH>){
	chomp;
	my @array=split /\t/, $_;
	my @out_array;

	my $map_val=$map_hash{$array[$KeyCol]};

	if(!defined($map_val)){
		$map_val="UNDEFINED";			
	}

	my $num_col=$#array+1;
	for(my $i=0; $i<($num_col+1); $i++){
		if($i==$InsertCol){
			push @out_array, $map_val;
		}else{
			push @out_array, (shift @array);
		}
	}

	print OUTPUT_FH (join "\t", @out_array) . "\n";

}

close(INPUT_FH);
close(OUTPUT_FH);

###############################################################################

print STDERR "Done.\n";
