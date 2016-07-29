#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_l);

getopts("l:");
my $usage = "usage: 
$0 
	-l <list of fasta files>

	This script will read in a list of fasta files and generate
	a list of the number of sequences and residues in each fasta file.

";

if(!(
	defined($opt_l))){
	die $usage;
}


my $fasta_file_list=$opt_l;

###############################################################################

open(LIST_FH, "<$fasta_file_list") || die "Could not open $fasta_file_list\n";

my @targets;
while(<LIST_FH>){
	chomp;
	push @targets, $_;	
}

close(LIST_FH);

###############################################################################

my $total_residues=0;
my $total_sequences=0;

foreach my $target(@targets){

	open(FASTA_FH, "<$target") || die "Could not open $target\n";

	print STDERR "Processing FASTA file: $target\n";

	$total_residues=0;
	$total_sequences=0;

	my ($defline, $prev_defline, $sequence);
	while(<FASTA_FH>){
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

	close(FASTA_FH);

	print STDOUT "$target\t$total_residues\t$total_sequences\n";

	print STDERR "Completed.\n";

}

###############################################################################


sub process_record{
	my $defline = shift;
	my $sequence = shift;

	$total_residues+=length($sequence);
	$total_sequences++;

}

###############################################################################
