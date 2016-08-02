#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_o);
use FileHandle;

getopts("i:o:");

my $usage = "
	usage:
	$0

	-i <list of input files, comma searated>
	-o <results blast out>
	
	This script will merge two blast outputs so that
	the hits/rows are sorted by subject ID and then by
	bit score.

	For Example:
		$0 -i file1,file2 -o results.blast

";

if(!defined($opt_i) || !defined($opt_o)){
	die $usage;
}

my $input_list=$opt_i;
my $output_file=$opt_o;

my @input_files=split ",", $input_list;

print STDERR "Input List:\n";
foreach my $file(@input_files){
	print "\t$file\n";
}
print STDERR "\n";
print STDERR "Output File: $output_file\n";

###############################################################################

my $inputfiles=join " ", @input_files;

my $sort_cmd=
	"sort -k12,12 -r -n $inputfiles -t '\t' | sort -k1,1 -s -t '\t' > $output_file";

system($sort_cmd);

