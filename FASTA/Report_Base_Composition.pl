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

	$sequence=uc($sequence);
	my %composition_hash;
	my @nucs=split //, $sequence;
	my $tot_length=$#nucs+1;
	my $tot_Ns=0;
	my $tot_Amb=0;
	my $tot_NonAmb=0;

	foreach my $nuc(@nucs){
		$composition_hash{$nuc}++;
		if($nuc=~/[CATG]/){
			$tot_NonAmb++;
		}else{
			$tot_Amb++;
			if($nuc=~/[N]/){
				$tot_Ns++;
			}
		}
	}

	print STDOUT "$defline\n";
	foreach my $nuc(sort keys %composition_hash){
		printf("$nuc\t%4i\t%5.2f%%\n", $composition_hash{$nuc}, 100.0*$composition_hash{$nuc}/$tot_length);
	}
	printf("NonAmbg\t%4i\t%5.2f%%\n", $tot_NonAmb, 100.0*$tot_NonAmb/$tot_length);
	printf("Ambg\t%4i\t%5.2f%%\n" , $tot_Amb, 100.0*$tot_Amb/$tot_length);
	printf("Ns\t%4i\t%5.2f%%\n", $tot_Ns, 100.0*$tot_Ns/$tot_length);
}

#------------------------------------------------------------------------------
