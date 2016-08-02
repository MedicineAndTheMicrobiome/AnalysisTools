#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_a $opt_b $opt_o);

getopts("a:b:o:");

my $usage = "
	usage:
	$0
		-a <Input filename, Primary>
		-b <Input filename, Secondary>
		-o <Output filename root>

	This script will read in two reaction files, and only based
	on the stoichiometry equation, remove redundancies prefering
	the reaction name of the primary (-a) file first.

";

if(!defined($opt_a) || !defined($opt_b) || !defined($opt_o)){
	die $usage;
}

my $Primary=$opt_a;
my $Secondary=$opt_b;
my $Output=$opt_o;

###############################################################################

print STDERR "Primary: $Primary\n";
print STDERR "Secondary: $Secondary\n";
print STDERR "Output: $Output\n";

###############################################################################

print STDERR "Processing Input File...\n";

###############################################################################

my $RXN_COL=2;

sub load_reaction{
	my $filename=shift;
	my $hash_ref=shift;

	open(IN_FH, "<$filename") || die "Could not open input file $filename\n";
	my $i=0;
	while(<IN_FH>){
		chomp;
		my @fields=split "\t", $_;
		my $rxn=$fields[$RXN_COL];
		$rxn=~s/\s+//g;

		${$hash_ref}{$rxn}=$_;
		$i++;
	}
	print STDERR "Num Reactions Loaded: $i\n";
	close(IN_FH);
}

###############################################################################

sub output_reaction{
	my $filename=shift;
	my $hash_ref=shift;

	open(OUT_FH, ">$filename") || die "Could not open output file $filename\n";

	my @reactions=sort keys %{$hash_ref};

	my $i=0;
	foreach my $reaction(@reactions){
		print OUT_FH "${$hash_ref}{$reaction}\n";
		$i++;
	}
	print STDERR "Num Reactions Output: $i\n";
	
	close(OUT_FH);
}


###############################################################################

my %hash;
load_reaction($Secondary, \%hash);
load_reaction($Primary, \%hash);
output_reaction($Output, \%hash);

###############################################################################

print STDERR "\n";
print STDERR "Done.\n";
