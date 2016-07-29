#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use Sys::Hostname;
use vars qw($opt_i);

getopts("i:");
my $usage = "usage: 
$0 
	-i <input fasta file>

	Reads in a FASTA file and assumes the first tokenizable string in the
	defline is the identifier.  Everything after this identifier
	is the description.  

	The output is:

	<sequence id>\\t<description>\\n

";

if(!defined($opt_i)){
	die $usage;
}

my $input_fasta=$opt_i;

###############################################################################

open(IN_FASTA, "<$input_fasta") || die "Could not open $input_fasta\n";

print STDERR "Processing FASTA file...\n";

my ($defline, $prev_defline, $sequence);
while(<IN_FASTA>){
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
	my @arr=split / /, $defline;

	my $id=shift @arr;
	my $desc=join " ", @arr;

	print "$id\t$desc\n";

}

###############################################################################
