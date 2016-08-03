#!/usr/bin/env perl

###############################################################################

use strict;
use FindBin ();
use Getopt::Std;
use FileHandle;
use File::Basename;
use vars qw($opt_i $opt_o $opt_n);

my $PIPELINE_UTIL_PATH="$FindBin::Bin/pipeline_utilities";
print STDERR "Path of Pipeline Utilities: $PIPELINE_UTIL_PATH\n";

my $SUBSAMPLE_PATH="$PIPELINE_UTIL_PATH/Randomly_Sample_from_FASTA.pl";

getopts("i:o:n:");
my $usage = "usage: 

$0 

	-i <path to input original fasta files>
	-o <path to output subsampled fasta files>
	-n <sample size>

	This script will go through all the fasta files
	in the specified input directory and generate subsamples of
	of sequences and save them into the output directory.

	The subsample will be the minimum of the sample size
	and the actual number of sequences that are available.

";

if(!(
	defined($opt_i) && 
	defined($opt_o) && 
	defined($opt_n))){
	die $usage;
}

my $input_path=$opt_i;
my $output_path=$opt_o;
my $sample_size=$opt_n;

print STDERR "\n";
print STDERR "Input Path: $input_path\n";
print STDERR "Output Path: $output_path\n";
print STDERR "Sample Size: $sample_size\n";

print STDERR `mkdir $output_path`;

###############################################################################

my @filelist=split "\n", `find $input_path`;
my @fastalist;

print STDERR "Found FASTA files: \n";
foreach my $fname(@filelist){
	if($fname=~/\.fasta$/ || $fname=~/\.fa$/){
		push @fastalist, $fname;
		print STDERR "$fname\n";
	}
}

foreach my $fpath(@fastalist){
	my ($name, $path)=fileparse($fpath);
	my $cmd="$SUBSAMPLE_PATH " . 
	"-f $fpath " . 
	"-n $sample_size " .
	"-s 1 " .
	"-r $output_path/$name " . 
	"-L ";

	print STDERR "\n\nExecuting: $cmd\n";

	system $cmd;
	
	
}


###############################################################################

print STDERR "done.\n";

