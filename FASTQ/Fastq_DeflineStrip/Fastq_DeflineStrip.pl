#!/usr/bin/env perl

use strict;
use Getopt::Std;
use vars qw($opt_f $opt_c);

getopts("f:c:");
my $usage = "
usage: 
$0 
	-f <fastq file>
	[-c <column number to remove, default = 1 is the first column>]

	This script will read a fastq file, split the 
	defline into multiple columns (based on spaces)
	then remove the column number specified with
	the -c.  

	The output will go into STDOUT, so you can
	run it through a compressor if you want to.

";

if(!(
	defined($opt_f)
)){
	die $usage;
}

my $fname=$opt_f;

my $column=1;
if(defined($opt_c)){
	$column=$opt_c;
}

###############################################################################

if($opt_f=~/\.gz$/){
	open(FASTQ_FH, "zcat $fname | ") || die "Could not open $opt_f\n";
}else{
	open(FASTQ_FH, "<$fname") || die "Could not open $opt_f\n";
}

###############################################################################

print STDERR "Reading in FASTQ file: $fname ...\n";

my $num_recs=0;

$column--;

while(!eof(FASTQ_FH)){

	# Parse FASTQ lines
	my $id=<FASTQ_FH>;
	my $seq=<FASTQ_FH>;
	my $plus=<FASTQ_FH>;
	my $qv=<FASTQ_FH>;
	
	# Extract defline
	my $defline=$id;
	chomp $defline;

	# Split
	$defline=~s/^@//;
	my @components=split / /, $defline;
	
	# Remove specified column
	my @new_comp;
	my $num_comp=$#components+1;
	for(my $i=0; $i<$num_comp; $i++){
		if($i != $column){
			push @new_comp, $components[$i];
		}
	}

	# Rejoin components
	my $new_defline=join " ", @new_comp;

	# Output FASTQ Record
	print STDOUT "@" . $new_defline . "\n";
	print STDOUT $seq;
	print STDOUT $plus;
	print STDOUT $qv;
}

close(FASTQ_FH);

print STDERR "Completed.\n";

###############################################################################
