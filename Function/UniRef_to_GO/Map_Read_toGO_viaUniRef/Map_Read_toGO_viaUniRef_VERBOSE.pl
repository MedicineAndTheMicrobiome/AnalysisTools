#!/usr/local/bin/perl

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

	This script assign a read through all the GO 'is a' relationships
	through the UniRef100 ID it had been previously assigned.
	
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
	push @{$go_to_parent_hash{$go_id}}, (split ";", $parent_ids);
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
        my $path=shift;
        my $paths_ref=shift;

	my $parent=$go_to_parent_hash{$in_child};

	if(defined($parent)){

		foreach my $parent_id(@{$parent}){			

			my $next_path;
			if($path eq ""){
				$next_path=$in_child;
			}else{
				$next_path="$path;$in_child";
			}

			lookup($parent_id, $next_path, $paths_ref);
		}
	}else{
		push @{$paths_ref}, "$path;$in_child";
		return;
	}
}

###############################################################################

open(OUT_FH, ">$output_file_root") || die "Could not open output file: $output_file_root\n";
open(IN_FH, "<$reads_to_uniref") || die "Could not open reads file: $reads_to_uniref\n";

while(<IN_FH>){
	chomp;
	my ($read_id, $uniref_str)=split /\t/, $_;
	
	if($uniref_str=~/^UniRef100/){

		my $go_arr_ref=$uniref_to_go_hash{$uniref_str};

		foreach my $go_id(@{$go_arr_ref}){
			print STDERR "$read_id\t$go_id\n";
			my @paths;
			lookup($go_id, "", \@paths);
			
			foreach my $path(@paths){
				print OUT_FH "$read_id\t$go_id\t$path\n";
			}
		}
	}	
}

close(OUT_FH);
close(IN_FH);

###############################################################################

print STDERR "Done.\n\n";



















