#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_p $opt_t);

getopts("p:t:");

my $usage = "
	usage:
	$0
		-p <MetaCyc Pathway (Key,Parent,Description) File>
		-t <target classification file>

";

if(!defined($opt_p) || !defined($opt_t)){ 
	die $usage;
}

my $pathway_info=$opt_p;
my $targets=$opt_t;

###############################################################################

print STDERR "\n";
print STDERR "MetaCyc Pathway File: $pathway_info\n";

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

#my @parents_with_gt2_children;

#foreach my $parent(keys %parent_to_num_children_hash){
#	if($parent_to_num_children_hash{$parent}>=1){
#		push @parents_with_gt2_children, $parent;
#	}
#}

#foreach my $parent(@parents_with_gt2_children){
#	my $str= join "|", @{$parent_to_child_hash{$parent}};
#	print "$parent\n";
#	print("\t$str\n");
#}

###############################################################################

sub lookup{
	my $child=shift;

	my $cur_child=$child;
	my @links;

	while(defined($child_to_parent_hash{$cur_child})){
		my $parent=$child_to_parent_hash{$cur_child};
		#push @links, $id_to_description_hash{$parent};
		push @links, $parent;
		$cur_child=$parent;
	}

	return(join ";", @links);

}

###############################################################################

open(TARGET_FH, "<$targets") || die "Could not open targets file $targets\n";

my $FOI=1;
while(<TARGET_FH>){
	chomp;
	my @fields=split "\t", $_;
	my $child=$fields[$FOI];
	my $hierstring=lookup($child);
	print "$child\t$hierstring\n";
}

close(TARGET_FH);

###############################################################################

print STDERR "Done.\n";
