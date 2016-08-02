#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_p $opt_c $opt_t $opt_o);

getopts("p:c:t:o:");

my $usage = "
	usage:
	$0
		-p <MetaCyc Reaction Info File (Reaction ID, EC, Enzyme ID, Pathway IDs)>
		-c <MetaCyc Pathway Categories (Pathway ID, Type, Description)>
		-t <Target File (Read ID, EC)>
		-o <Output File Root>
	
	This script will assign a pathway ID to each Read ID via the EC number associated with it.
	
	The EC field in the target file may be a single EC or in the format of 'EC1 || EC2 || .. || ECn'.
	
	The output is denomalized, so a single Read ID may have multiple ECs, with each associated
	with multiple pathways.

";

if(!defined($opt_p) || !defined($opt_c) || !defined($opt_t) || !defined($opt_o)){ 
	die $usage;
}

my $reaction_info=$opt_p;
my $pathway_categories=$opt_c;
my $target_file=$opt_t;
my $output_file_root=$opt_o;

###############################################################################

print STDERR "\n";
print STDERR "Reaction Info File: $reaction_info\n";
print STDERR "Pathway Categories File: $pathway_categories\n";
print STDERR "Target File: $target_file\n";
print STDERR "Output File Root: $output_file_root\n";

###############################################################################
# Load EC to Pathway mapping

my %ec_to_pathway_hash;
my %pathway_to_class_hash;;

open(IN_FH, "<$reaction_info") || die "Could not open reaction file $reaction_info\n";

while(<IN_FH>){
	chomp;
	my ($rxn_id, $ec_str, $enzname, $pathway_str)=split /\t/, $_;

	if($ec_str eq ""){
		next;
	}else{
		my @ecs=split ";", $ec_str;
		my @pws=split ";", $pathway_str;

		foreach my $ec (@ecs){
			$ec=~s/^EC-//;

			if(!defined($ec_to_pathway_hash{$ec})){
				@{$ec_to_pathway_hash{$ec}}=();
			}

			push @{$ec_to_pathway_hash{$ec}}, @pws;
		}
	}
}

close(IN_FH);

###############################################################################
# Load Pathway to pathway Type mapping

my %pathway_type_hash;

open(IN_FH, "<$pathway_categories") || die "Could not open pathway file $pathway_categories\n";

while(<IN_FH>){
	chomp;
	my ($pathway_id, $type, $description)=split /\t/, $_;
	@{$pathway_type_hash{$pathway_id}}=split ";", $type;
}

close(IN_FH);

###############################################################################
# Go through the input file and assign a pathway and type to each read

my $num_unspecific_ecs++;

open(IN_FH, "<$target_file") || die "Could not open target file $target_file\n";

open(OUT_VERBOSE_FH, ">$output_file_root\.tsv") || die "Could not open output $output_file_root\.tsv\n";

while(<IN_FH>){
	chomp;
	my ($read_id, $ec_str);
	my @cols=split /\t/, $_;
	my $numCol=scalar(@cols);
#	print STDERR "Num of Cols in target file:$numCol\n";
	if(scalar(@cols) == 2){
		($read_id, $ec_str)=split /\t/, $_;
	}
	else {
		 $ec_str=split /\t/, $_;
		 $read_id="";
	}
#	print STDERR "$read_id\t$ec_str\n";

#	my @ecs=split / \|\| /, $ec_str;
	my @ecs=split /;/, $ec_str;

	foreach my $ec (@ecs){
	#	print STDERR "$ec\n";

		if(defined($ec_to_pathway_hash{$ec})){
			foreach my $pathway(@{$ec_to_pathway_hash{$ec}}){
	#			print STDERR "\t$pathway\n";
				
				foreach my $type(@{$pathway_type_hash{$pathway}}){

					print OUT_VERBOSE_FH
						"$read_id\t$ec\t$pathway\t$type\n";

				}
			}
		}else{
			print OUT_VERBOSE_FH
                                                "$read_id\t$ec\tunknown\tunknown\n";
			$num_unspecific_ecs++;
		}

	}

}

close(OUT_VERBOSE_FH);

###############################################################################

print STDERR "Num unspecific EC's: $num_unspecific_ecs\n";

print STDERR "Done.\n\n";
