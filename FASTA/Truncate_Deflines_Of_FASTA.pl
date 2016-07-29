#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;

my $usage = "usage: 

$0 < [input FASTA file]  > [output FASTA file]

This program will read in a FASTA file and output a FASTA file with 
the defline truncated.  Only the identifier following the > will
be saved.  All tags will be removed.

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
		$sequence.="$_\n";
	}
}
process_record($prev_defline, $sequence);

print STDERR "Completed.\n";

###############################################################################

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	my @fields=split / /, $defline;
	$defline=$fields[0];

	print STDOUT "$defline\n";
	print STDOUT "$sequence";
}

#------------------------------------------------------------------------------
