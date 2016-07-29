#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;

my $usage = "usage: 

$0 < [input FASTA file]  > [tab delimited file]

Takes as input a fasta file, and generated a tab delimited file.
This makes a sequence record that spans multiple lines fit onto a single line.

";

###############################################################################

print STDERR $usage;

print STDERR "Processing FASTA file...\n";

my $record_count=0;

my ($defline, $prev_defline, $sequence);
while(<STDIN>){
	chomp;
	
	if(/^>/){
		$defline=$_;
		if($sequence ne ""){
			process_record($prev_defline, $sequence);
			$sequence="";
		}
		$prev_defline=$defline;
	}else{
		$sequence.=$_;
	}
}
process_record($prev_defline, $sequence);

print STDERR "Completed.\n";

###############################################################################

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	$defline=~s/^>//;
	print STDOUT "$defline\t$sequence\n";

}

#------------------------------------------------------------------------------
