#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;

my $usage = "usage: 
$0 
	
	This program will read in a FASTA file from STDIN and report the
	GC content for each record.

	Output will be:

	>id1\\t<GC>
	>id2\\t<GC>
	...
	>idn\\t<GC>
	WeightedAverage\\t<Average GC>
	
	
";

print STDERR "$usage\n";

###############################################################################

my $overall_gc=0;
my $overall_nuc=0;

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

my $weightedAverage=$overall_gc/$overall_nuc;
print "WeightedAverage\t$weightedAverage\n";


print STDERR "Completed.\n";

###############################################################################


sub process_record{
	my $defline = shift;
	my $sequence = shift;

	my $id;
	if($defline=~/^>(\S+)/){
		$id=$1;
	}

	$sequence=uc($sequence);
	my @nuc=split //, $sequence;
	my $length=$#nuc+1;

	my $gc_count=0;
	foreach my $base(@nuc){
		if($base eq "G" || $base eq "C" || $base eq "S" ){
			$gc_count+=1.000;
		}elsif($base eq "R" || $base eq "M"){
			$gc_count+=0.500;
		}elsif($base eq "V" || $base eq "B"){
			$gc_count+=0.666;
		}elsif($base eq "D" || $base eq "H"){
			$gc_count+=0.333;
		}elsif($base eq "N" || $base eq "X"){
			$gc_count+=0.500;
		}
	}

	my $perc_gc=$gc_count/$length;
	
	$overall_gc+=$gc_count;
	$overall_nuc+=$length;

	print "$id\t$perc_gc\n";

}

#------------------------------------------------------------------------------

