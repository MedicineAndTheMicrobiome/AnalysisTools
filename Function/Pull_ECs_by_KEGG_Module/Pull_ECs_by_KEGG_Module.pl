#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_r $opt_m $opt_p);

getopts("r:m:p:");

my $usage = "
usage:
	$0
		-r <KEGG Reaction File (from kegg-ligand directory)>

	One of:	
		[-m <KEGG Module, e.g. M00173>]
		[-p <KEGG Pathway, e.g. rn00620>]

	This script will extract out all the reaction information for
	the specified KEGG reaction module or pathway.

	Output Format:
		1.) E.C. ID
		2.) R##### (Reaction ID)
		3.) Reaction name

	The output file name will be the:
		<KEGG module/pathway ID>.tsv
	
";

if(!defined($opt_r) || (!defined($opt_m) && !defined($opt_p))){
	die $usage;
}

my $reaction_filename=$opt_r;
my $target_module=$opt_m;
my $target_pathway=$opt_p;

print STDERR "Reaction input file: $reaction_filename\n";

my $target_name;
if(defined($target_module)){
	print STDERR "Target Module: $target_module\n";
	$target_name=$target_module;
}

if(defined($target_pathway)){
	print STDERR "Target Pathway: $target_pathway\n";
	$target_name=$target_pathway;
}

print STDERR "Target ID: $target_name\n";

###############################################################################

sub search{
	my $id=shift;
	my $recs_arr=shift;
	
	my $num_recs=$#{$recs_arr}+1;

	for(my $i=0; $i<$num_recs; $i++){
		my $rec=${$recs_arr}[$i];
		if($rec=~/(\S+)/){
			if($1 eq $id){
				return(1);
			}
		}
	}
	return(0);
}

###############################################################################

open(FH, "<$reaction_filename") || die "Could not open $reaction_filename\n";
open(OUT_FH, ">$target_name\.tsv") || die "Could not open $target_name.tsv\n";

my %data_hash;
my $cur_rec_type;

while(<FH>){

	chomp $_;
	
	if($_=~/^ /){
		 # carry over from previous record 
		my $val=$_;
		$val=~s/^\s+//;
		push @{$data_hash{$cur_rec_type}}, $val;	

	}else{
		if($_=~/^([A-Z]+)/){
			$cur_rec_type=$1;

			my $val="";
			if($_=~/^$cur_rec_type\s+(.+)/){
				$val=$1;
			}

			push @{$data_hash{$cur_rec_type}}, $val;	

		}elsif($_=~/^\/\/\//){

			my @key_arr=keys %data_hash;

			# Check if reaction is part of module/pathway
			my $found=0;	
			if(defined($target_module)){
				$found=search($target_module, $data_hash{"MODULE"});		
			}elsif(defined($target_pathway)){
				$found=search($target_pathway, $data_hash{"PATHWAY"});
			}

			# If so, output the record
			if($found){
				my @EC=split /\s+/, $data_hash{"ENZYME"}[0];
				my $rxn_id=$data_hash{"ENTRY"}[0];
				my $name=$data_hash{"NAME"}[0];

				$rxn_id=(split /\s+/, $rxn_id)[0];

				# Some reactions have more than one EC
				foreach my $ec(@EC){
					print OUT_FH "$ec\t$rxn_id\t$name\n";
				}
			}

			%data_hash=();
		}
	}


}
close(FH);
close(OUT_FH);

###############################################################################

print STDERR "Done.\n";
