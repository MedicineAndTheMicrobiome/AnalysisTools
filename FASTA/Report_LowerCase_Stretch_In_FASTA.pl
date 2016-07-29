#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_f);

getopts("f:");
my $usage = "usage: 
$0 
	This program will extract the positions of lower case
	letters in 0-space based coordinates.

	-f <fasta file name>

	For example, for the FASTA record:

	>FLA
	AaAaaaaAaA

	The output will be:

	>FLA
	1       2
	3       7
	8       9

	Where a tab separates the begin/end of each stretch of lower case
	letters, and the number of lower case stretches will equal the 
	number of rows in the output following an echo of the input record's 
	defline.

";

if(!(defined($opt_f))){
	die $usage;
}

###############################################################################

open(FASTA_FH, "<$opt_f") || die "Could not open $opt_f\n";

print STDERR "Processing FASTA file...\n";

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

print STDERR "Completed.\n";

###############################################################################

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	print STDOUT "$defline\n";

	my $length=length($sequence);
	my $i;

	my $was_lc=0;
	my $is_lc;

	for($i=0; $i<$length; $i++){
		my $nuc=substr($sequence,$i,1);		

		if($nuc eq lc($nuc)){
			$is_lc=1;	
		}else{
			$is_lc=0;
		}

		if(!$was_lc && $is_lc){
			print "$i\t";
			$was_lc=1;
		}elsif($was_lc && !$is_lc){
			print "$i\n";
			$was_lc=0;
		}
	}

	if($was_lc){
		print "$i\n";
	}

}

#------------------------------------------------------------------------------
