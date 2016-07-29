#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_f $opt_n $opt_o $opt_s);

getopts("f:n:o:s:");
my $usage = "usage: 

$0 
	-f <fastq file>
	-n <sample size>
	-o <output file name root>
	[-s <random number seed>]

	This script will uniformly subsample from the target fastq file.
	The sampling is WITHOUT REPLACEMENT, so the subsample will converge
	on the input file if the input is not large.

	This subsampling should NOT be used for bootstrapping.
	
	Two files will be output:
		1.) The FASTQ subsample
		2.) An index/id file.  
			You can confirm the subsampling was uniform by looking
			at the first column.  You also generate a list of IDs
			that were sampled by cutting the second column.

	If you want to subsample from forward and reverse reads on two
	different runs, then you will want to set the random numbers seed
	to be the same between runs.

";

if(!(
	defined($opt_f) && 
	defined($opt_n) &&
	defined($opt_o))){
	die $usage;
}

my $input_fname=$opt_f;
my $sample_size=$opt_n;
my $output_fname=$opt_o;
my $random_number_seed=$opt_s;

if(!defined($random_number_seed)){
	$random_number_seed=time;
}
srand($random_number_seed);

print STDERR "\n";
print STDERR "Input FASTQ File: $input_fname\n";
print STDERR "Requested Sample Size: $sample_size\n";
print STDERR "Output Filename Root: $output_fname\n";
print STDERR "Random number seed: $random_number_seed\n";
print STDERR "\n";

###############################################################################
# Make sure files open before wasting any time doing anything

if($input_fname=~/\.gz$/){
	open(FASTQ_FH, "zcat $input_fname| ") || die "Could not open $input_fname\n";
}else{
	open(FASTQ_FH, "<$input_fname") || die "Could not open $input_fname\n";
}

###############################################################################

print STDERR "Counting FASTQ file: $input_fname ...\n";
my ($defline, $prev_defline, $sequence);
my $i;

my $num_recs=0;

while(!eof(FASTQ_FH)){

	my $id=<FASTQ_FH>;
	my $seq=<FASTQ_FH>;
	my $plus=<FASTQ_FH>;
	my $qv=<FASTQ_FH>;

	chomp $plus;
	if($plus ne "+"){
		die "Error: + character not found between sequence and quality.\n";
	}

	$num_recs++;
}

print STDERR "Number of Records: $num_recs\n";

close(FASTQ_FH);

###############################################################################

if($num_recs<$sample_size){
	print STDERR "WARNING: Sample Size ($sample_size) is greater than Number of Records ($num_recs)\n";
	print STDERR "Setting sample size to number of records.\n";
	$sample_size=$num_recs;
}

# Generate a list of unique random indices
my %rand_idx;
while((keys %rand_idx)<$sample_size){
	$rand_idx{int(rand($num_recs))}=1;
}

#my $i=0;
#foreach my $key(sort keys %rand_idx){
#	print "$i $key\n";
#	$i++;
#}

###############################################################################
# Rescan over input fastq file and write out records

if($opt_f=~/\.gz$/){
	open(FASTQ_FH, "zcat $input_fname| ") || die "Could not open $input_fname\n";
}else{
	open(FASTQ_FH, "<$input_fname") || die "Could not open $input_fname\n";
}

$output_fname=~s/\.fastq$//;
open(OUTPUT_FH, ">$output_fname\.fastq") || die "Could not open $output_fname\.fastq for writing.\n";

print STDERR "Reading/Sampling from FASTQ file: $input_fname ...\n";

my @ids;
my $rec_idx=0;

while(!eof(FASTQ_FH)){

        my $id=<FASTQ_FH>;
        my $seq=<FASTQ_FH>;
        my $plus=<FASTQ_FH>;
        my $qv=<FASTQ_FH>;

	if(defined($rand_idx{$rec_idx})){
		print OUTPUT_FH "$id";
		print OUTPUT_FH "$seq";
		print OUTPUT_FH "$plus";
		print OUTPUT_FH "$qv";

		chomp $id;
		push @ids, "$rec_idx\t$id";
	}

        $rec_idx++;
}

close(OUTPUT_FH);

###############################################################################
# Output IDs and index of subsampled records

open(OUTPUT_FH, ">$output_fname\.subsamp_id.txt") || die "Could not open $output_fname\.subsamp_id.txt for writing.\n";

foreach my $id(@ids){
	print OUTPUT_FH "$id\n";
}

close(OUTPUT_FH);

###############################################################################

print STDERR "Done.\n\n";

