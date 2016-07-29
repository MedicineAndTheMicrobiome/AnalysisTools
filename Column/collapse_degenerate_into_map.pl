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

		$0 -i inputfile -c 1

	This script will read in a file with a degenerate mapping,
	and then produce a distinct mapping file where the keep column
	values are separated by semicolons.  The source column is
	assumed to be column 0.	

	Example input:

		id1	apples
		id1	bananas
		id2	coriander
		id3	penguins
		id3	falcons
		id3	poultry

	Example output:

		id1	apples;bananas
		id2	coriander
		id3	penguins;falcons;poultry
		

	Output goes to STDOUT

";

if(!defined($opt_i) || !defined($opt_c)){
	die $usage;
}

my $file=$opt_i;
my $column=$opt_c;

print STDERR "Input File: $file\n";
print STDERR "Column Index: $column\n";

###############################################################################

open(IN_FH, "<$file") || die "Could not open map file $file\n";

my %double_hash;

while(<IN_FH>){
	chomp;
	my @in=split /\t/, $_;
	
	my $key=$in[0];
	my $val=$in[$column];

	if(!defined($double_hash{$key})){
		my %inner_hash;
		$inner_hash{$val}=1;
		$double_hash{$key}=\%inner_hash;
	}else{
		${$double_hash{$key}}{$val}=1;
	}
}

close(IN_FH);


foreach my $key (sort keys %double_hash){
	
	my @val_arr=sort keys %{$double_hash{$key}};

	my $outstr=join ";", @val_arr;

	print STDOUT "$key\t$outstr\n";


}

