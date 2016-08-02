#!/bin/env perl

use strict;
use Getopt::Std;
use vars qw ($opt_b $opt_u);

getopts("b:u:");

my $usage = "
	usage:
	$0

	-b <btab file>

";

if(!defined($opt_b)){
	die $usage;
}

my $alignments_file=$opt_b;

###############################################################################

my $READ_ID_COL=0;
my $QRY_LEN_COL=2;
my $PERC_SIM_COL=11;
my $QRY_ST_COL=6;
my $QRY_EN_COL=7;
my $UNIREF_ID_COL=5;

###############################################################################

open(BTAB_FH, "<$alignments_file") || die "Could not open $alignments_file\n";

while(<BTAB_FH>){
	chomp;
	my @fields=split "\t", $_;

	my $read_id=$fields[$READ_ID_COL];
	my $qry_len=$fields[$QRY_LEN_COL];
	my $perc_sim=$fields[$PERC_SIM_COL]/100.0;
	my $qry_start=$fields[$QRY_ST_COL];
	my $qry_end=$fields[$QRY_EN_COL];
	my $uniref_id=$fields[$UNIREF_ID_COL];

	$uniref_id=~s/^UniRef100_//;

	if($qry_start>$qry_end){
		($qry_start, $qry_end)=($qry_end, $qry_start);
	}

	my $aln_len=$qry_end-$qry_start+1;
	my $perc_len=$aln_len/$qry_len;
	my $perc_comp=$perc_sim*$perc_len;
	$perc_comp=sprintf("%1.4f", $perc_comp);

	print "$read_id\t$perc_comp\t$uniref_id\n";
}

close(BTAB_FH);

###############################################################################

print STDERR "Done.\n";
