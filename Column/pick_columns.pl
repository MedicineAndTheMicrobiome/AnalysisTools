#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_c);

getopts("i:c:");

my $usage = "
	usage:
	$0
		-i <input file, tab separated columns>
		-c <columns to keep, starting from 1, comma separated>
		
	For example:

		$0 -i inputfile -c 2,1,3

	will change the order of the input file, and eliminate unspecified columns.

	You can also specify a column twice if you want it duplicated.


";

if(!defined($opt_i) || !defined($opt_c)){
	die $usage;
}

my $file=$opt_i;
my $columns=$opt_c;

my @column_order=split /,/, $columns;

###############################################################################

open(IN_FH, "<$file") || die "Could not open map file $file\n";

while(<IN_FH>){
	chomp;
	my @in=split /\t/, $_;
	
	my @output;
	foreach my $col(@column_order){
		push @output, $in[$col-1];
	}

	my $out_str=join "\t", @output;
	print "$out_str\n";
}

close(IN_FH);

