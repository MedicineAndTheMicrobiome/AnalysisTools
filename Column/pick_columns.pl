#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_k $opt_d);

getopts("i:k:d:");

my $usage = "
	usage:
	$0
		-i <input file, tab separated columns>

		One of:
		[-k <columns to keep, starting from 1, comma-separated>]
		[-d <columns to delete, starting from 1, comma-separated>]

	Output goes to STDOUT.
		
	Keep example:

		$0 -i inputfile -k 2,1,3

	will change the order of the input file, and eliminate unspecified columns.
	You can also specify a column multiple times if you want it duplicated.

	Delete example:

		$0 -i inputfile -d 2

	will delete column 2 and keep the remaining columns in their original order.


";

if(!defined($opt_i)){
	die $usage;
}

my $file=$opt_i;
my $keep_columns=$opt_k;
my $delete_columns=$opt_d;

my @columns;

my $keep=-1;

if(defined($keep_columns)){
	@columns=split /,/, $keep_columns;
	$keep=1;
}

my %delete_hash;
if(defined($delete_columns)){
	@columns=split /,/, $delete_columns;
	foreach my $ix(@columns){
		$delete_hash{$ix}=1;
	}
	$keep=0;
}

if($keep==-1){
	die "You must specify whether you want columsn deleted (-d) or kept (-k).\n";
}

###############################################################################

open(IN_FH, "<$file") || die "Could not open map file $file\n";

while(<IN_FH>){
	chomp;
	my @in=split /\t/, $_, -1;
	
	my @output;

	if($keep){
		# Keep specified columns
		foreach my $col(@columns){
			push @output, $in[$col-1];
		}
	}else{
		# Exclude specified columns
		for(my $i=0; $i<=$#in; $i++){
			if(!defined($delete_hash{($i+1)})){
				push @output, $in[$i];
			}	
		}
	}

	my $out_str=join "\t", @output;
	print "$out_str\n";
}

close(IN_FH);

