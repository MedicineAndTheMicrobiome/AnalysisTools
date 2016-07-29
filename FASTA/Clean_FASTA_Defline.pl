#!/usr/local/bin/perl

###############################################################################

use strict;
use Getopt::Std;
use Sys::Hostname;
use vars qw($opt_i $opt_o);

getopts("i:o:");
my $usage = "usage: 
$0 
	-i <input fasta file>
	-o <output fasta file>

	Reads in a FASTA file and assumes the first tokenizable string in the
	defline is the identifier.  Everything after this identifier
	is truncated.

";

if(!defined($opt_i) || !defined($opt_o)){
	die $usage;
}

my $input_fasta=$opt_i;
my $output_fasta=$opt_o;

###############################################################################

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

	my @arr=split /\s+/, $defline;
	print OUT_FASTA "$arr[0]\n";

	$sequence=~s/\s+//g;

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
