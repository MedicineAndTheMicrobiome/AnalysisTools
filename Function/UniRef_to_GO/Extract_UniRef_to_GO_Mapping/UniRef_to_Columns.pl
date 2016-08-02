#!/usr/local/bin/perl

use strict;
use Getopt::Std;
use vars qw ($opt_m $opt_o);

getopts("m:o:");

my $usage = "
	usage:
	$0
		-m <UniProt ID Mapping Tab file>
		-o <Output file name>

	This script will extract UniRef100 to GO mappings.

	If a GO mapping is not found, GO:-1 is used as a placeholder.
";

if(!defined($opt_m) || !defined($opt_o)){ 
	die $usage;
}

my $dat_file=$opt_m;
my $output_fn=$opt_o;

###############################################################################

print STDERR "\n";
print STDERR "UniProt to IDs Mapping Tab File: $dat_file\n";
print STDERR "Output File: $output_fn\n";
print STDERR "\n";

###############################################################################

open(IN_FH, "<$dat_file") || die "Could not open dat file $dat_file\n";
open(OUT_FH, "| sort -u > $output_fn") || die "Could not open output file $output_fn\n";

my $UNIREF100_POS=8;
my $GO_POS=6;

while(<IN_FH>){
	chomp;
	my @fields=split "\t", $_;

	if($fields[$GO_POS] ne "" && $fields[$UNIREF100_POS] ne ""){
		my @gos=split /; /, $fields[$GO_POS];
		my $out=join ";", @gos;
		print OUT_FH "$fields[$UNIREF100_POS]\t$out\n";
	}else{
		#print OUT_FH "$fields[$UNIREF100_POS]\tGO:-1\n";
	}

}

close(IN_FH);

###############################################################################

print STDERR "Done.\n";
