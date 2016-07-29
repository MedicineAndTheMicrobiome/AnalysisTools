#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_f);

getopts("f:");
my $usage = "usage: 
$0 
	This program will change all the ambiguity codes in a sequence to N.
	
	-f <Input FASTA file>
";

if(!defined($opt_f)){
	die "$usage";
}

###############################################################################

print STDERR "Processing FASTA file...\n";

open(FASTA_FH, "<$opt_f") || die "Could not open $opt_f\n";

my ($defline, $prev_defline, $sequence);
while(<FASTA_FH>){
	
	if($_=~/^>/){
		print $_;
	}else{
		my $seq=uc($_);
		$seq=~s/[MRWSYKVHDBX\*]/N/g;
		print $seq;
	}

}
print STDERR "Completed.\n";

#------------------------------------------------------------------------------
