#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_n);

getopts("i:n:");

my $usage = "
	usage:
	$0
		-i <input file>
		-n <sub sample size>
	
	This script will randomly subsample (from uniform distribution) from 
	the rows in the input file.  The order of the rows is preserved.
	Sampling is without replacement, so if your sample size is larger
	than the number of rows in the file, you will get back the original
	file.

";

if(!defined($opt_i) || !defined($opt_n)){
	die $usage;
}

my $file=$opt_i;
my $sample_size=$opt_n;

print STDERR "Input File: $file\n";
print STDERR "Sample Size: $sample_size\n";

###############################################################################

print STDERR "Counting number of rows...\n";

my $i=0;
open(IN_FH, "<$file") || die "Could not open $file\n";
while(<IN_FH>){
	$i++;		
}
close(IN_FH);

my $num_rows=$i;
print STDERR "Num rows: $i\n";
print STDERR "\n";

###############################################################################

# Generate random positions
print STDERR "Generating random samples...\n";
my %random_positions;
my $num_positions=0;
do{
	my $pos=int(rand()*$num_rows);
	if(!defined($random_positions{$pos})){
		$random_positions{$pos}=1;
		$num_positions++;
	}

}while($num_positions<$sample_size && $num_positions<$num_rows); 

# Extract positions are we read in the file
print STDERR "Outputting subsamples..\n";
open(IN_FH, "<$file") || die "Could not open $file again.\n";
my $in;
for(my $i=0; $i<$num_rows; $i++){
	$in=<IN_FH>;
	if(defined($random_positions{$i})){
		print STDOUT $in;
	}	
}

close(IN_FH);

###############################################################################

print STDERR "Done.\n";

