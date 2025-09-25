#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i);

getopts("i:c:t:r");

my $usage = "
	usage:
	$0
		-i <input file>

	This is parse the metacyc reactions.dat file so we can from EC to Pathway.

";

if(!defined($opt_i)){
	die $usage;
}

my $file=$opt_i;

###############################################################################

open(IN_FH, "<$file") || die "Could not open $file\n";

open(OUT_FH, ">$file\.ec_pathway.tsv") || die "Could not open $file\.ec_pathway.tsv";

my %EC_to_Pathways_hash;

my $cur_rxn="";
my $cur_EC="";
my $cur_pathway="";

my $ec_list="";
my $pathway_list="";

while(<IN_FH>){
	chomp;
	
	#print STDOUT "$_\n";
	
	if($_=~/^UNIQUE-ID - (.+)/){
		$cur_rxn=$1;
		print "$cur_rxn\n";
	}elsif($_=~/^EC-NUMBER - (.+)/){
		$cur_EC=$1;
		$cur_EC=~s/^EC-//;
	
		my $num_dots=($cur_EC=~tr/\././);
		$cur_EC="$cur_EC" . (".-" x (3-$num_dots));

		if($ec_list eq ""){
			$ec_list=$cur_EC;
		}else{
			$ec_list="$ec_list;$cur_EC";
		}

	}elsif($_=~/^IN-PATHWAY - (.+)/){
		$cur_pathway=$1;

		if($pathway_list eq ""){
			$pathway_list=$cur_pathway;
		}else{
			$pathway_list="$pathway_list;$cur_pathway";	
		}

	}elsif($_=~/^\/\//){
		print "\t\t$ec_list / $pathway_list\n";

		my @ec=split ";", $ec_list;
		my @pw=split ";", $pathway_list;

		foreach my $ecid(@ec){
			foreach my $pwid(@pw){
				print OUT_FH "$ecid\t$pwid\t$cur_rxn\n";
			}
		}

		$ec_list="";
		$pathway_list="";
	}
	
}

close(IN_FH);
close(OUT_FH);

###############################################################################

