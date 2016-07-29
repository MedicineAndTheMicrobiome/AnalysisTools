#!/usr/bin/env perl

use strict;
use Getopt::Std;
use vars qw($opt_f $opt_s $opt_r $opt_o);

getopts("f:s:r:o:");
my $usage = "
usage: 
$0 
	-f <fastq file>

	This script will strip off the /1 or /2 from
	the end of the sequence id by converting it
	into ' 1' or ' 2'.

	The output will go into STDOUT, so you can
	run it through a compressor if you want to.

";

if(!(
	defined($opt_f)
)){
	die $usage;
}

my $fname=$opt_f;

###############################################################################

if($opt_f=~/\.gz$/){
	open(FASTQ_FH, "zcat $fname | ") || die "Could not open $opt_f\n";
}else{
	open(FASTQ_FH, "<$fname") || die "Could not open $opt_f\n";
}

###############################################################################

print STDERR "Reading in FASTQ file: $fname ...\n";

my $num_recs=0;

while(!eof(FASTQ_FH)){

	my $id=<FASTQ_FH>;
	my $seq=<FASTQ_FH>;
	my $plus=<FASTQ_FH>;
	my $qv=<FASTQ_FH>;

	my $defline=$id;
	chomp $defline;
	
	if($defline=~/\/\d$/){
		$defline=~s/\/(\d)$//;

		print STDOUT $defline . " $1\n";
		print STDOUT $seq;
		print STDOUT $plus;
		print STDOUT $qv;
		
	}else{
		die "Error, could not find /1, /2, etc.. at end of defline.\n";
	}
}

close(FASTQ_FH);

print STDERR "Completed.\n";

###############################################################################
