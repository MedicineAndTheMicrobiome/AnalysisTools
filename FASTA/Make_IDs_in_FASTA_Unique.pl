#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use Sys::Hostname;
use vars qw($opt_i $opt_o);

getopts("i:o:");
my $usage = "
usage: 

$0 
	-i <input fasta file>
	-o <output fasta file>
	
	Reads in FASTA, and if it detects any IDs that are not unique
	it will append a numeric counter to the ID, the next time it sees the
	same ID.

	The first time an identifier is found, it will not be modified.
	Only subsequent degenerate IDs will be renamed.

	For Example:

	    >example
	    ACTGACTAGCATCAGCTG
	    >example
	    TCGACTGATCGATCT

	Will change into:

	    >example
	    ACTGACTAGCATCAGCTG
	    >example_1
	    TCGACTGATCGATCT

";

if(!defined($opt_i) || !defined($opt_o)){
	die $usage;
}

my $input_fasta=$opt_i;
my $output_fasta=$opt_o;

###############################################################################

my %id_hash;

open(IN_FASTA, "<$input_fasta") || die "Could not open $input_fasta\n";
open(OUT_FASTA, ">$output_fasta") || die "Could not open $output_fasta\n";

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

	if($defline=~/^>(\S+)/){
		my $id=$1;
		my $new_id;
		if(!defined($id_hash{$id})){
			$id_hash{$id}=0;
			$new_id=$id;
		}else{
			$id_hash{$id}++;
			$new_id="$id\_$id_hash{$id}";
		}

		print STDERR "$id -> $new_id\n";
		$defline=~s/^>$id/>$new_id/;
	}else{
		die "Error parsing defline.\n";
	}

	print OUT_FASTA "$defline\n";

	my $length=length($sequence);
	my $width=80;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print OUT_FASTA substr($sequence, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);
}

###############################################################################
