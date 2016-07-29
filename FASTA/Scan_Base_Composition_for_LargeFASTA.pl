#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw();

getopts("");
my $usage = "usage: 
$0 
	This program will read in a sequence fasta file though STDIN and then
	generate some statistics on the base composition.	
";

print STDERR "$usage";

###############################################################################

print STDERR "Processing FASTA file...\n";

my %nuc_hash;
my $hash_has_content=0;

my ($defline, $prev_defline, $sequence);
while(<STDIN>){
	chomp;
	
	if(/^>/){
		print "$defline\n";
		if($hash_has_content){
			report_nuc_hash();	
		}
		%nuc_hash=();		
	}else{
		$hash_has_content=1;
		my @bases=split //, $_;
		compute_base_comp(\@bases);
	}
}
report_nuc_hash();

print STDERR "Completed.\n";

###############################################################################

sub compute_base_comp{
	my $base_arr_ref=shift;
	foreach my $base(@{$base_arr_ref}){
		$nuc_hash{$base}++;
	}
}

sub report_nuc_hash{
	foreach my $base(sort keys %nuc_hash){
		print "$base: $nuc_hash{$base}\n";	
	}

}


#------------------------------------------------------------------------------
