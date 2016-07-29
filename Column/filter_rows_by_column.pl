#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_t);

getopts("i:t:");

my $usage = "
	usage:
	$0
		-i <input file>
		-t <list of targets to keep file>

	Reads in input file and keeps rows where the first column is in the target list.
	It's like grep, but we only look in the first column and the match is exact.

";

if(!defined($opt_i) || !defined($opt_t)){
	die $usage;
}

my $file=$opt_i;
my $keep_list_file=$opt_t;

###############################################################################

my %keep_list_hash;

open(KEEP_FH, "<$keep_list_file") || die "Could not open $keep_list_file\n";
while(<KEEP_FH>){
	chomp;
	my @array=split /\t/, $_;
	$keep_list_hash{$array[0]}=1;
}
close(KEEP_FH);

###############################################################################

open(INPUT_FH, "<$file") || die "Could not open $file\n";
while(<INPUT_FH>){
	chomp;
	my @array=split /\t/, $_;
	if($keep_list_hash{$array[0]}){
		print STDOUT "$_\n";
	}
}
close(INPUT_FH);

###############################################################################
