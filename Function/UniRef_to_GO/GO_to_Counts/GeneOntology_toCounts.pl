#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_n $opt_o $opt_c $opt_m);
use FileHandle;

getopts("n:o:c:m:");

my $usage = "
	usage:
	$0
		-n <Ontology Assignments (to individual reads)>
		-o <Output Filename Root>
		-c <Category Depth File>
		-m <GO id to description map>

	This script will build a tree based on the GO ID mapping
	and then populate counts based on the read assignments.

";

if(!defined($opt_n) || !defined($opt_o) || !defined($opt_c) || !defined($opt_m)){ 
	die $usage;
}

my $ontologies=$opt_n;
my $output=$opt_o;
my $category_depth=$opt_c;
my $go_id_map_file=$opt_m;


###############################################################################

print STDERR "\n";
print STDERR "GO Mapping: $go_id_map_file\n";
print STDERR "Ontology Assignments: $ontologies\n";
print STDERR "Category Depth File: $category_depth\n";
print STDERR "Output File: $output\n";
print STDERR "\n";

###############################################################################

my %counts_hash;	# Store counts for each GO ID
my %tree_hash;		# Store relationship between parent and child
my %description_hash;	# Store description for each GO ID

###############################################################################
# Load descriptions

open(GO_MAP, "<$go_id_map_file") || die "Could not open map file: $go_id_map_file\n";

while(<GO_MAP>){
	chomp;
	my ($go_id, $parent_id_str, $desc)=split /\t/, $_;

	$description_hash{$go_id}=$desc;
	$counts_hash{$go_id}=0;

	my @parent_ids=split /;/, $parent_id_str;

	foreach my $parent_id(@parent_ids){
		push @{$tree_hash{$parent_id}}, $go_id;
	}
}

print STDERR "Done loading GO mapping.\n\n";


#my @kids=@{$tree_hash{"GO:0006450"}};
#foreach my $kid(@kids){
#	print "$kid\n";
#}

###############################################################################
# Load category depth file

my @category_arr;
my @depth_arr;
my $num_catdep=0;

open(CAT_DEP, "<$category_depth") || die "Could not open category depth file: $category_depth\n";

while(<CAT_DEP>){
	chomp;
	if(substr($_,0,1) eq "#"){
		next;
	}
	my ($category, $depth)=split /\t/, $_;
	my $desc=$description_hash{$category};
	print STDERR "$category [$desc] Depth: $depth\n";
	push @category_arr, $category;
	push @depth_arr, $depth;
	$num_catdep++;
}

close(CAT_DEP);

print STDERR "Done loading category depth file.\n\n";

###############################################################################

sub add_counts{
	my $go_ids_arr_ref=shift;
	foreach my $go_id(@{$go_ids_arr_ref}){
		$counts_hash{$go_id}++;
	}
}

sub draw_tree_counts{
	my $parent_id=shift;
	my $level=shift;
	my $fh=shift;
	my $depth_limit=shift;

	if($depth_limit==0){
		return;
	}
	
	if($counts_hash{$parent_id}>0){
		print {$fh} "$level$description_hash{$parent_id} [$parent_id] ($counts_hash{$parent_id})\n";
		
		if(defined($tree_hash{$parent_id})){
			my @children=sort @{$tree_hash{$parent_id}};
			foreach my $child(@children){
				draw_tree_counts($child, $level . "\t", $fh, $depth_limit-1);
			}
		}
	}

}

###############################################################################

open(IN_FH, "<$ontologies") || die "Could not open ontology assignment file $ontologies\n";

while(<IN_FH>){
	chomp;
	my @fields=split /\t/, $_;

	my $read_id=$fields[0];
	my $go_ids_str=$fields[1];
	
	my @go_ids=split /;/, $go_ids_str;
	add_counts(\@go_ids);
	
}
close(IN_FH);

print STDERR "\nDone summing up counts.\n";

###############################################################################
# Outputs a complete tree, for reference

for(my $i=0; $i<$num_catdep; $i++){
	my $fh=FileHandle->new();
	my $fname="$output\.$category_arr[$i]\.$depth_arr[$i]\.tree.txt";
	$fh->open(">$fname") || die "Could not open $fname\n";
	draw_tree_counts($category_arr[$i], "" , $fh, $depth_arr[$i]);
	$fh->close;
}

###############################################################################

sub sum_of{
	my $arr_ref=shift;
	my $counts=0;
	foreach my $goid(@{$arr_ref}){
		$counts+=$counts_hash{$goid};
	}
	return($counts);
}

sub produce_table{
	my $parent_id=shift;
	my $ancestry=shift;
	my $depth=shift;
	my $accumulation_arr_ref=shift;

	if($depth==0){
		return;
	}

	if($counts_hash{$parent_id}>0){

		my $output_category=0;

		if($depth==1){
			# If we've hit the depth limit;
			#print "*";
			$output_category=1;
		}elsif(!defined($tree_hash{$parent_id})){
			# If the parent has no more children
			#print "+";
			$output_category=1;
		}elsif(sum_of($tree_hash{$parent_id})==0){
			# If the parent has children, but they are all 0 count ghosts
			#print "#";
			$output_category=1;
		}

		#print "$depth $ancestry:$parent_id ($counts_hash{$parent_id})\n";

		if($output_category==1){
			if($ancestry eq ""){
				push @{$accumulation_arr_ref}, "$parent_id#$counts_hash{$parent_id}";
			}else{
				push @{$accumulation_arr_ref}, "$ancestry;$parent_id#$counts_hash{$parent_id}";
			}
		}
	
		# Go into children
		if(defined($tree_hash{$parent_id})){
			my @children=sort @{$tree_hash{$parent_id}};
			foreach my $child(@children){

				my $next_ancestry;
				if($ancestry eq ""){
					$next_ancestry=$parent_id;
				}else{
					$next_ancestry="$ancestry;$parent_id";
				}

				produce_table($child, $next_ancestry, $depth-1, $accumulation_arr_ref);
			}
		}
	}

}

###############################################################################
# Outputs counts for each specified format in the category depth file

for(my $i=0; $i<$num_catdep; $i++){

	my $category=$category_arr[$i];
	my $depth=$depth_arr[$i];

	print STDERR "Preparing counts for $category at depth $depth...\n";
	
	#----------------------------------------------------------------------
	# Extract table

	my @acc_arr=();
	produce_table($category, "", $depth, \@acc_arr);
	
	#----------------------------------------------------------------------
	# Break table down into columns
	
	my @categories;
	my @counts;
	my $total_counts=0;
	my $num_recs=0;
	foreach my $catcount(@acc_arr){
		my ($category, $count)=split /#/, $catcount;
		push @categories, $category;
		push @counts, $count;
		$total_counts+=$count;
		$num_recs++;
	}	

	#----------------------------------------------------------------------
	# Output summary as a list (by rows)

	my $listfh=FileHandle->new();
	my $output_fn="$output\.$category.$depth.summary_list.txt";
	$listfh->open(">$output_fn") || die "Could not open $output_fn\n";

	for(my $i=0; $i<$num_recs;  $i++){

		my @go_id=split /;/, $categories[$i];
		my $desc=$description_hash{(pop @go_id)};
		print {$listfh} "$categories[$i]\t$counts[$i]\t$desc\n";
	}

	$listfh->close;

	#----------------------------------------------------------------------
	# Output summary table (by column)
	my $tabfh=FileHandle->new();
	my $output_fn="$output\.$category.$depth.summary_table.tsv";
	$tabfh->open(">$output_fn") || die "Could not open $output_fn\n";

	# First row
	print {$tabfh} "SampleID\ttotal\t";
	print {$tabfh} (join "\t", @categories);
	print {$tabfh} "\n";

	# Second row
	print {$tabfh} "$output\t$total_counts\t";
	print {$tabfh} (join "\t", @counts);
	print {$tabfh} "\n";

	$tabfh->close;
	
	#----------------------------------------------------------------------

	print STDERR "   ok.\n";
	
}

###############################################################################

print STDERR "\nDone.\n";
