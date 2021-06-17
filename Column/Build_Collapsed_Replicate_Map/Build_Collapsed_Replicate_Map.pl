#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_f $opt_c $opt_o);

my $DEF_FIELDS=5;
my $DEF_COLNAME="collapsed";

getopts("i:f:c:o:");

my $usage = "
	usage:
	$0
		-i <input file name, e.g. summary table, only 1st column is examined.>
		-o <output filename map>

		[-c <collapsed column name, default=", $DEF_COLNAME, ">]
		[-f <num_fields, for replicate_id-free sample ID, def=", $DEF_FIELDS, ">]

	This script will read in the input file and generate a
	mapping for collapsing replicates based on the assumption
	that there are num_fields that identify the underlying
	sample.

	For example (original):

	  0159.AY44.33055.20171011.ST.C02
	  0159.AY44.33066.20171009.ST.2

	Would get collapsed to (collapsed): 

	  0159.AY44.33055.20171011.ST
	  0159.AY44.33066.20171009.ST

	if the -f num_fields option was set to 5

	The output will be two columns:

	<column name>\\t<collapsed column name>\\n
	<original id1>\\t<collapsed id1>\\n
	<original id2>\\t<collapsed id2>\\n
	<original id3>\\t<collapsed id3>\\n
	...
	<original idn>\\t<collapsed idn>\\n


	If the original does not have the specified number of fields,
	then the collapsed id column will match the original.

	
";



if(
	!defined($opt_i) || 
	!defined($opt_o)
)
{
	die $usage;
}

my $input_fn=$opt_i;
my $output_fn=$opt_o;
my $num_fields=$DEF_FIELDS;
my $column_name=$DEF_COLNAME;

if(defined($opt_f)){
	$num_fields=$opt_f;
}

if(defined($opt_c)){
	$column_name=$opt_c;
}

print STDERR "Input File: $input_fn\n";
print STDERR "Output File: $output_fn\n";
print STDERR "Collapsed Column Name: $column_name\n";
print STDERR "Num Fields: $num_fields\n";

###############################################################################

open(IN_FH, "<$input_fn") || die "Could not open $input_fn for reading.\n";
open(OUT_FH, ">$output_fn") || die "Could not open $output_fn for writing.\n";

###############################################################################

my $line_ix=0;
while(<IN_FH>){

	chomp;

	my @fields=split /\t/, $_;

	if($line_ix==0){
		print OUT_FH "$fields[0]\t$column_name\n";
	}else{

		my $original=$fields[0];
		my @arr=split /\./, $original;	

		my $num_found_fields=$#arr+1;
		my $target_fields=($num_found_fields>$num_fields)?$num_fields:$num_found_fields;

		my @keep=();
		for(my $i=0; $i<$target_fields; $i++){
			push @keep, $arr[$i];

		}

		my $trimmed=join ".", @keep;
		print OUT_FH "$original\t$trimmed\n";
	}

	$line_ix++;
}

print STDERR "Done.\n";

