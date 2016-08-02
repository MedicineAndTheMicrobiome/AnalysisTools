#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_r $opt_u $opt_g $opt_o);

getopts("r:u:g:o:");

my $usage = "
	usage:
	$0
		-r <Read to UniRef100 Mapping>
		-u <UniRef100 to GO Mapping>
		-g <GO ID to GO Parent ID Mapping>
		-o <Output File Root>

	For each read, generates a list of GO that it its or 'is a' of.
	
";

if(!defined($opt_r) || !defined($opt_u) || !defined($opt_g) || !defined($opt_o)){ 
	die $usage;
}

my $reads_to_uniref=$opt_r;
my $uniref_to_go=$opt_u;
my $go_to_parent=$opt_g;
my $output_file_root=$opt_o;

###############################################################################

print STDERR "\n";
print STDERR "Reads to UniRef: $reads_to_uniref\n";
print STDERR "UniRef to GO: $uniref_to_go\n";
print STDERR "GO to GO Parent: $go_to_parent\n";
print STDERR "Output File Root: $output_file_root\n";
print STDERr "\n";

###############################################################################

my %uniref_to_go_hash;

print STDERR "Loading UniRef to GO Mapping: $uniref_to_go\n";
open(IN_FH, "<$uniref_to_go") || die "Could not open UniRef to GO Mapping: $uniref_to_go\n";

while(<IN_FH>){
	chomp;
	my ($uniref_id, $go_ids)=split "\t", $_;
	if($go_ids ne "GO:-1"){
		my @go_arr=split /;/, $go_ids;
		push @{$uniref_to_go_hash{$uniref_id}}, @go_arr;	
	}
}

close(IN_FH);
print STDERR "done.\n";

#
#foreach my $uniref_id(keys %uniref_to_go_hash){
#	print STDERR "$uniref_id\n";
#	foreach my $go_map(@{$uniref_to_go_hash{$uniref_id}}){
#		print "\t$go_map\n";
#	}
#}

###############################################################################

print STDERR "Loading GO to GO Parent Mapping: $go_to_parent\n";
open(IN_FH, "<$go_to_parent") || die "Could not open GO to GO Parent Mapping: $go_to_parent\n";

my %go_to_parent_hash;
my %go_description_hash;
while(<IN_FH>){
	chomp;
	my ($go_id, $parent_ids, $description)=split "\t", $_;
	push @{$go_to_parent_hash{$go_id}}, (split /;/, $parent_ids);
	$go_description_hash{$go_id}=$description;
}
print STDERR "done.\n";


#foreach my $go_id(keys %go_to_parent_hash){
#	print STDERR "$go_id\n";
#	foreach my $go_parents(@{$go_to_parent_hash{$go_id}}){
#		print "\t$go_parents\n";
#	}
#}

###############################################################################

sub lookup{
        my $in_child=shift;
	my $GO_hits_hash_ref=shift;

	if($in_child eq "GO:-1"){
		return;
	}

	${$GO_hits_hash_ref}{$in_child}=1;

	my $parent=$go_to_parent_hash{$in_child};

	if(defined($parent)){
		#if($#{$parent}>0){
		#	print STDERR "Several parents.\n";
		#}
		foreach my $parent_id(@{$parent}){			
			lookup($parent_id, $GO_hits_hash_ref);
		}
	}else{
		print STDERR "Error: Undefined parents for $in_child.\n";
		return;
	}
}

###############################################################################

open(OUT_FH, ">$output_file_root") || die "Could not open output file: $output_file_root\n";
open(IN_FH, "<$reads_to_uniref") || die "Could not open reads file: $reads_to_uniref\n";

# Global isa cache
my %go_isas_hash;

while(<IN_FH>){

	chomp;
	my ($read_id, $uniref_str)=split /\t/, $_;
	
	if($uniref_str=~/^UniRef100/){

		# Retrive GO IDs based on UniRef100 ID
		my $go_arr_ref=$uniref_to_go_hash{$uniref_str};

		# Keep track of all hit GO IDs for this read
		my %read_specific_go_isas;

		# Identify all the isa relationships for each GO ID associated with read
		foreach my $go_id(@{$go_arr_ref}){

			# Look up isa IDs for this GO ID in a cache first
			my $go_isas_ref=$go_isas_hash{$go_id};

			# If GO ID is not in cache, look it up and save it
			if(!defined($go_isas_ref)){
				my %GO_hits;

				lookup($go_id, \%GO_hits);
				my @go_isas_arr=(keys %GO_hits);

				$go_isas_hash{$go_id}=\@go_isas_arr;
				$go_isas_ref=\@go_isas_arr;
				
				print STDERR "+";
			}else{
				print STDERR ".";
			}

			# Store all the isas for this read together
			foreach my $isa(@{$go_isas_ref}){
				$read_specific_go_isas{$isa}=1;
			}
		}

		# Pull out all the unique isa from the read specific hash
		my @isas_arr=sort keys %read_specific_go_isas;

		# Build a string of all the isa GO IDs and output with read ID
		my $isa_str=join ";", @isas_arr;
		print OUT_FH "$read_id\t$isa_str\n";

	}	
}

close(OUT_FH);
close(IN_FH);

###############################################################################

print STDERR "Done.\n\n";



















