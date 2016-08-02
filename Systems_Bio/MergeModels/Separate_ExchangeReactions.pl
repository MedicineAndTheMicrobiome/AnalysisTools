#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_i $opt_o);

getopts("i:o:");

my $usage = "
	usage:
	$0
		-i <Input filename>
		-o <Output filename root>

	This script will read in the input reactions file and 
	split out two files:
		1.) Exchange Reactions
		2.) Cellular and Extracellular Reactions

";

if(!defined($opt_i) || !defined($opt_o)){
	die $usage;
}

my $InputFile=$opt_i;
my $OutputFile=$opt_o;

###############################################################################

print STDERR "Input: $InputFile\n";
print STDERR "Output: $OutputFile\n";

###############################################################################

print STDERR "Processing Input File...\n";
my $input_lines=0;

open(IN_FH, "<$InputFile") || die "Could not open input file $InputFile\n";
open(OUT_CELL_FH, ">$OutputFile.cell.tsv") || die "Could not open input file $OutputFile.cell.tsv\n";
open(OUT_EXCH_FH, ">$OutputFile.exchange.tsv") || die "Could not open input file $OutputFile.exchange.tsv\n";

my $fh;
my $temp_filename;

my $REACTION_COLUMN=2;
my $num_cell_rxn=0;
my $num_exch_rxn=0;
while(<IN_FH>){
	chomp;
	my @array=split /\t/, $_;

	$#array=9;
	my $outstr=join "\t", @array;

	# Grab reaction column
	my $reaction=$array[$REACTION_COLUMN];

	if($_=~/<==>/){
		my ($lhs, $rhs)=split /<==>/, $reaction;
		$rhs=~s/\s+//g;

		if($rhs eq ""){
			print OUT_EXCH_FH "$outstr\n";
			$num_exch_rxn++;
		}else{
			print OUT_CELL_FH "$outstr\n";
			$num_cell_rxn++;
		}
	}else{
		print OUT_CELL_FH "$outstr\n";
		$num_cell_rxn++;
	}

	# Output pulse
	$input_lines++;
}

close(IN_FH);
close(OUT_FH);

###############################################################################

print STDERR "\n";
print STDERR "Num Lines Processed: $input_lines\n";
print STDERR "Num Cellular Reactions: $num_cell_rxn\n";
print STDERR "Num Exchange Reactions: $num_exch_rxn\n";
print STDERR "Done.\n";
