#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_n $opt_o $opt_c);
use FileHandle;

getopts("n:o:c:");

my $usage = "
	usage:
	$0
		-n <Ontology Assignments (to individual reads)>
		-o <Output Filename Root>
		-c <Category Depth File>


	The Ontology Assiments File Format:

	<read id>	Sugars-And-Polysaccharides;Carbohydrates-Degradation;Degradation;Pathways


	The Category Depth Files Format:

	<CategoryID>\\t<Depth>\\n

	CategoryID could be: Pathways, Degradation, Amino-Acid-Degradation, etc.
	Depth is the number of levels to retrieve counts for.  A depth of 0 will
		only retrieve the count for the CategoryID
	

";

if(!defined($opt_n) || !defined($opt_o) || !defined($opt_c)){ 
	die $usage;
}

my $ontologies=$opt_n;
my $output=$opt_o;
my $category_depth=$opt_c;


###############################################################################

print STDERR "\n";
print STDERR "Ontology Assignments: $ontologies\n";
print STDERR "Output File: $output\n";
print STDERR "Category Depth File: $category_depth\n";
print STDERR "\n";

###############################################################################

my %counts_hash;
my %tree_hash;

###############################################################################
# Load category depth file

my @category_arr;
my @depth_arr;
my $num_catdep=0;

open(CAT_DEP, "<$category_depth") || die "Could not open $category_depth\n";

while(<CAT_DEP>){
	chomp;
	my ($category, $depth)=split /\t/, $_;
	push @category_arr, $category;
	push @depth_arr, $depth;
	$num_catdep++;
}

close(CAT_DEP);


###############################################################################

sub update_tree{
	my $comp_ref=shift;

	#print STDERR "Updating Tree:\n";

	for(my $i=0; $i<$#{$comp_ref}; $i++){
		my $child=${$comp_ref}[$i];
		my $parent=${$comp_ref}[$i+1];

		${$tree_hash{$parent}}{$child}=1;
	}
}

sub update_counts{
	my $unique_keys_ref=shift;
	my @keys=keys %{$unique_keys_ref};

	foreach my $key(@keys){
		if(!defined($counts_hash{$key})){
			$counts_hash{$key}=0;
		}
		$counts_hash{$key}++;

		#print STDERR "\t($key)\n";
	}

	#print STDERR "\n\n";
}

sub draw_tree_counts{
	my $parent_id=shift;
	my $level=shift;
	my $fh=shift;
	
	print {$fh} "$level $parent_id ($counts_hash{$parent_id})\n";
	foreach my $child( keys %{$tree_hash{$parent_id}}){
		draw_tree_counts($child, $level . "\t", $fh);
	}

}

sub produce_table{
	my $parent_id=shift;
	my $parent=shift;
	my $depth=shift;
	my $accumulation_arr_ref=shift;

	#print STDERR "Looking at: $parent_id\n";
	my @keys=(keys %{$tree_hash{$parent_id}});
	my $num_kids=$#keys+1;

	if($depth==0 || $num_kids==0){
		push @{$accumulation_arr_ref}, "$parent#$counts_hash{$parent_id}";
	}else{
		foreach my $child(keys %{$tree_hash{$parent_id}}){
			produce_table($child, "$parent:$child", $depth-1, $accumulation_arr_ref); 
		}

	}
}


###############################################################################

sub process_reads{
	my $group_ref=shift;

	my %unique_keys;

	for(my $i=0; $i<=$#{$group_ref}; $i++){
		#print STDERR "${$group_ref}[$i]\n";	
		my @components=split ";", ${$group_ref}[$i];	

		for(my $j=0; $j<=$#components; $j++){
			$unique_keys{$components[$j]}=1;
		}

		update_tree(\@components);
	}

	update_counts(\%unique_keys);	

}

###############################################################################


open(IN_FH, "<$ontologies") || die "Could not open ontology assignment file $ontologies\n";

my @group=();
my $prior_read_id="";

while(<IN_FH>){
	chomp;
	my @fields=split /\t/, $_;

	my $read_id=$fields[0];
	my $ontology=$fields[2];

	if($read_id ne $prior_read_id){
		process_reads(\@group);
		@group=();
	}

	push @group, $ontology;
	$prior_read_id=$read_id;
	
}
process_reads(\@group);

close(IN_FH);

###############################################################################
# Outputs a complete tree, for reference

my $fh=FileHandle->new();
$fh->open(">$output\.tree") || die "Could not open $output\.tree\n";
draw_tree_counts("Pathways", "", $fh);
$fh->close;

###############################################################################
# Outputs counts for each specified format in the category depth file

for(my $i=0; $i<$num_catdep; $i++){

	my $category=$category_arr[$i];
	my $depth=$depth_arr[$i];

	print STDERR "Preparing counts for $category at depth $depth...\n";
	
	#----------------------------------------------------------------------
	# Extract table

	my @acc_arr=();
	produce_table($category, $category, $depth, \@acc_arr);
	
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
		print {$listfh} "$categories[$i] ($counts[$i])\n";
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
	
}

###############################################################################

print STDERR "\nDone.\n";
