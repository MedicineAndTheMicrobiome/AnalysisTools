#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_p $opt_t $opt_o);

getopts("p:t:o:");

my $usage = "
	usage:
	$0
		-p <MetaCyc Pathway (Key,Parent,Description) File>
		-t <Target classification file, where the column to do the lookup is on column 4.>
		-o <Output File>

	This script will map the target classification file up the MetaCyc pathway.

	For example:

		The pathway classification:

			Sugars-And-Polysaccharides      
			
		Will get traced up the MetaCyc tree

			Sugars-And-Polysaccharides;Carbohydrates-Degradation;Degradation;Pathways

";

if(!defined($opt_p) || !defined($opt_t) || !defined($opt_o)){ 
	die $usage;
}

my $pathway_info=$opt_p;
my $targets=$opt_t;
my $output=$opt_o;

###############################################################################

print STDERR "\n";
print STDERR "MetaCyc Pathway File: $pathway_info\n";
print STDERR "Targets File: $targets\n";
print STDERR "Output File: $output\n";

###############################################################################

my %id_to_description_hash;
my %child_to_parent_hash;

open(IN_FH, "<$pathway_info") || die "Could not open dat file $pathway_info\n";

while(<IN_FH>){
	chomp;
	my ($id, $parent, $description)=split /\t/, $_;

	$id_to_description_hash{$id}=$description;
	$child_to_parent_hash{$id}=$parent;
}

close(IN_FH);

###############################################################################

# Truncate meaningless roots

delete($child_to_parent_hash{"Pathways"});

###############################################################################

my %parent_to_child_hash;
my %parent_to_num_children_hash;
# Invert mappings

my @children_ids=keys(%child_to_parent_hash);

foreach my $child(@children_ids){

	my $parent=$child_to_parent_hash{$child};
	
	if(!defined($parent_to_child_hash{$parent})){
		@{$parent_to_child_hash{$parent}}=();
		$parent_to_num_children_hash{$parent}=0;
	}

	push @{$parent_to_child_hash{$parent}}, $child;
	$parent_to_num_children_hash{$parent}++;
}

###############################################################################

sub lookup{
# This function will recursive travel up the classification tree for each
# pathway type of interest
#
	my $in_child=shift;
	my $path=shift;
	my $paths_ref=shift;

	my @children=split ";", $in_child;

	foreach my $child (@children){

		my $parent=$child_to_parent_hash{$child};
		
		if(defined($parent)){

			my $next_path;
			if($path eq ""){
				$next_path=$child;
			}else{
				$next_path="$path;$child";
			}

			lookup($parent, $next_path, $paths_ref);
		}else{
			push @{$paths_ref}, "$path;$child";
			return;		
		}

	}


}

###############################################################################

open(TARGET_FH, "<$targets") || die "Could not open targets file $targets\n";
open(OUTPUT_FH, ">$output") || die "Could not open output file $output\n";

my $FOI=3; # Counting from 0

while(<TARGET_FH>){
	chomp;
	my @fields=split "\t", $_;
	my $child=$fields[$FOI];
	my $read_id=$fields[0];
	
	my @paths=();
	lookup($child, "", \@paths);

	foreach my $path(@paths){
		print OUTPUT_FH "$read_id\t$child\t$path\n";
	}
}

close(TARGET_FH);
close(OUTPUT_FH);

###############################################################################

print STDERR "Done.\n";
