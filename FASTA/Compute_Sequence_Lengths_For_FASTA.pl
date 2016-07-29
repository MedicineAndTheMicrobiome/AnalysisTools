#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_i);
getopts("i");

my $usage="
	
	$0 
		[-i Flag for only outputting ID, not entire defline]

	Pipe FASTA file into this script and it will compute the length
	of each record and report it.  First column will be the defline
	of the record, the second column will be the length of the record's
	sequence.  The columns are separated by a tab.  Output goes to
	STDOUT.

";

print STDERR $usage;

###############################################################################

print STDERR "Waiting for input through STDIN...\n";

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

	my $length=length($sequence);
	if(defined($opt_i)){
		$defline=~/^>(\S+)/;
		print STDOUT "$1\t$length\n";
	}else{
		print STDOUT "$defline\t$length\n";
	}
}

#------------------------------------------------------------------------------
