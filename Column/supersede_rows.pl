#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_t);

getopts("i:t:");

my $usage = "
	usage:
	$0

	Reads in all the files you specify and generates a 
	new mapping of id to data.  The data in the file specified 
	later in the list supersedes the prior one.

	Output is written to STDOUT.
	
	If you place a : after the file name, the number after the colon
	wil specific which column, starting from 0, to use as the data column.

	For example,

	$0 data1.txt:1 data2:2

	will read in the second column from data1 and the 3rd column from data2.
	Anything keys/data pairs in data2 will supersede what is data1.
	

";

###############################################################################

my $file_info;
my %hash;
while($file_info=shift){

	my ($filename, $target_column)=split /:/, $file_info;

	if($target_column eq "" || !defined($target_column)){
		$target_column=1;
	}

	open(FH, "<$filename") || die "Could not open $filename\n";

	while(<FH>){
		chomp;
		my @arr=split /\t/, $_;

		$hash{$arr[0]}=$arr[$target_column];
	}

	close(FH);

}

foreach my $key(sort keys %hash){
	print STDOUT "$key\t$hash{$key}\n";
}

###############################################################################
