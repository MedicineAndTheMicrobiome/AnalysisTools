#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_l $opt_c $opt_d);

getopts("i:l:c:d");

my $usage = "
	usage:
	$0
		-i <input file>
		-l <target column>
		-c <cutoffs, comma delimited>

		[-d <flag to produce id only list>]

	Reads in input file and based on the column you specify
	partions off the input into multiple outputs, depending
	on the cutoffs you specified.

	For example, if you specify 10,20,30

	You will have 3 files:

	      10 <= x < 20,  in file <input file>.part_10
	      20 <= x < 30,  in file <input file>.part_20
	      30 <= x < inf, in file <input file>.part_30

	If you want to capture all values, you must know the lower bound.
	
";

if(!defined($opt_i) || !defined($opt_l) || !defined($opt_c)){
	die $usage;
}

my $file=$opt_i;
my $column=$opt_l;
my $cutoffs=$opt_c;
my $make_id_list=defined($opt_d);

###############################################################################

my @partition_cutoffs = split /,/, $cutoffs;
@partition_cutoffs=sort {$b<=>$a} @partition_cutoffs;

###############################################################################

open(IN_FH, "<$file") || die "Could not open $file\n";

my @array;
my @acc_arr;
my %partition_info;

while(<IN_FH>){
	chomp;
	my @array=split /\t/, $_;
	foreach my $partition(@partition_cutoffs){
		if($array[$column]>=$partition){
			push @{$partition_info{$partition}}, $_;
			last;
		}
	}
}
close(IN_FH);

###############################################################################

my $outfname_root="$file\.part_";
my $num_length=int(log($partition_cutoffs[0])/log(10))+1;

foreach my $partition(@partition_cutoffs){
	my $part_id=sprintf("%0$num_length"."i", $partition);
	open (OUT_FH, ">$outfname_root$part_id") || die "Could not open $outfname_root$part_id\n";
	foreach my $line(@{$partition_info{$partition}}){
		print OUT_FH "$line\n";
	}
	close(OUT_FH);

	if($make_id_list){
		`cut -f 1 $outfname_root$part_id > $outfname_root$part_id\.id_list`;
	}
}

###############################################################################
