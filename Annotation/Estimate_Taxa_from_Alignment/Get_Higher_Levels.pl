#!/usr/bin/env perl

use strict;
use Getopt::Std;
use vars qw ($opt_t $opt_T $opt_m $opt_d $opt_l $opt_o $opt_h);

getopts("t:T:m:d:l:o:h");

my $usage = "
	usage:
	$0
	-t <taxa file (read_id, taxa_id) >
	-T <taxa file columns, eg. \"1,7\", for read_id and taxa_id, respectively, counting from 1>

	-m <taxa names file (taxa_id, name) >
	-d <taxa nodes file (parent_id, child_id) >
	-l <taxa levels file (taxa_id, taxa_level) >

	[-h (suppress header flag)]
	-o <output file name root>

	This script will take in a file which contains a query id (read id) 
	and the taxa id that was assigned to it, and generate the list of
	taxonomic levels above it starting from the superkingdom.

	The output contains the field names in the first row (header)

	The output contains:
		query_id, taxa_id, 
		superkingdom, phylum, class, order, family, genus, and species.

	If the input taxa_id is at a higher taxonomic level than the taxa
	level of the column, then the lowest taxonomic level that could be assigned
	is used a place holder (padding) within parenthesis.  

	Two outputs are produced:
		*.taxa_ids.tsv:	contains the taxa id of each of the levels
		*.taxa_names.tsv:  contains the name of each of the levels
";

###############################################################################

if(!defined($opt_t) || !defined($opt_T) ||  !defined($opt_m) || !defined($opt_d) || 
	!defined($opt_d) || !defined($opt_o)){
	die $usage;
}

my $taxa_ids_file=$opt_t;
my $taxa_ids_file_columns=$opt_T;
my $taxa_names_file=$opt_m;
my $taxa_nodes_file=$opt_d;
my $taxa_levels_file=$opt_l;
my $output_file=$opt_o;
my $suppress_header=($opt_h eq "1");

###############################################################################

my %levels_to_report;

$levels_to_report{"superkingdom"}=0;
$levels_to_report{"kingdom"}=1;
$levels_to_report{"phylum"}=2;
$levels_to_report{"class"}=3;
$levels_to_report{"order"}=4;
$levels_to_report{"family"}=5;
$levels_to_report{"genus"}=6;
$levels_to_report{"species"}=7;

###############################################################################
###############################################################################

# Loads the taxa nodes into a hash.  Each child points to its parent.
sub load_taxa_nodes{
	my $fn=shift;
	print STDERR "Loading taxa nodes: $fn\n";
	open(FH, "<$fn") || die "Could not open $fn for reading.\n";
	my %nodes;
	my $i=0;
	while(<FH>){
		chomp;
		my ($child, $parent)=split "\t", $_;
		$nodes{$child}=$parent;
		$i++;
	}
	close(FH);
	print STDERR "Nodes loaded: $i\n";
	return(\%nodes);
}

###############################################################################

# Loads the levels into a hash
sub load_taxa_levels{
	my $fn=shift;
	print STDERR "Loading taxa levels: $fn\n";
	open(FH, "<$fn") || die "Could not open $fn for reading.\n";
	my %levels;
	my $i=0;
	while(<FH>){
		chomp;
		my ($taxid, $level)=split "\t", $_;
		$levels{$taxid}=$level;
		$i++;
	}
	close(FH);
	print STDERR "levels loaded: $i\n";
	return(\%levels);
}

###############################################################################

# Loads the taxa names into a map.  Spaces are changed into underscores.
sub load_taxa_names{
	my $fn=shift;
	print STDERR "Loading taxa names: $fn\n";
	open(FH, "<$fn") || die "Could not open $fn for reading.\n";
	my %nodes;
	my $i=0;
	while(<FH>){
		chomp;
		my ($taxa_id, $name)=split "\t", $_;
		$name=~s/ /_/g;	# Changes spaces to underscores
		$nodes{$taxa_id}=$name;
		$i++;
	}
	close(FH);

	$nodes{""}="NULL";
	$nodes{"NULL"}="NULL";

	print STDERR "Names loaded: $i\n";
	return(\%nodes);
}

###############################################################################

# Converts ancestory to array of taxa names
sub anc_arr_to_str{
	my $anc_ref=shift;
	my $taxa_names_ref=shift;
	my @anc_str_arr;
	my $anc_str="";
	foreach my $anc (@{$anc_ref}){
		push @anc_str_arr, ${$taxa_names_ref}{$anc};
	}
	$anc_str=join ";", @anc_str_arr;
	return($anc_str);
}

# Keep cache of ancestries for observed taxa, so we don't have to recompute each time
my %ancestry_cache;

# Builds an ancestry for a given taxa id
sub build_ancestry_arr{
	my $in_taxa_id=shift;
	my $nodes_hash_ref=shift;

	if(!defined($ancestry_cache{$in_taxa_id})){

		my @ancestry;
		my $parent_taxa_id;
		my $taxa_id=$in_taxa_id;

		do{
			# Place taxa id in front, so first element is the root
			unshift @ancestry, $taxa_id; 

			# Look up parent
			$parent_taxa_id=${$nodes_hash_ref}{$taxa_id};

			if(!defined($parent_taxa_id)){
				print STDERR "Error: $taxa_id not found in nodes hash (Parent not found).\n";
				$taxa_id=1;
			}else{
				$taxa_id=$parent_taxa_id;
			}

		}while(($taxa_id ne "1") && defined($taxa_id));

		#print STDERR (join ";", @ancestry) . "\n";
		$ancestry_cache{$in_taxa_id}=\@ancestry;
	}

	return($ancestry_cache{$in_taxa_id});
}

###############################################################################

my $nodes_ref;
my $names_ref;
my $levels_ref;

$nodes_ref=load_taxa_nodes($taxa_nodes_file);
$names_ref=load_taxa_names($taxa_names_file);
$levels_ref=load_taxa_levels($taxa_levels_file);

# Compare an array of ancentry arrays
sub get_higher_levels{

	my $taxa_id=shift;
	my $nodes_ref=shift;
	my $levels_ref=shift;
	my @taxa_at_levels;

	my $anc_arr_ref=build_ancestry_arr($taxa_id, $nodes_ref);

	my $num_levels=$#{$anc_arr_ref}+1;

	@taxa_at_levels=("NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL");

	for(my $i=0; $i<$num_levels; $i++){
		my $cur_taxa_id=${$anc_arr_ref}[$i];
		my $level_str=${$levels_ref}{$cur_taxa_id};
		my $level_idx=$levels_to_report{$level_str};
		if(!defined($level_idx)){ 
			next;
		}
		$taxa_at_levels[$level_idx]=$cur_taxa_id;
	}

	return(\@taxa_at_levels);
}


###############################################################################

if(0){
	my $taxa_levels_arr_ref=get_higher_levels("292415", $nodes_ref, $levels_ref);
	print "------\n";
	for(my $i=0; $i<=$#{$taxa_levels_arr_ref}; $i++){
		my $taxa = ${$taxa_levels_arr_ref}[$i];
		print "[$i] $taxa ${$names_ref}{$taxa}\n";
	}
	exit;
}

###############################################################################

open(ID_RESULTS_FH, ">$output_file\.taxa_ids.tsv") || die "Could not open $output_file\.taxa_ids.tsv\n";
open(NAMES_RESULTS_FH, ">$output_file\.taxa_names.tsv") || die "Could not open $output_file\.taxa_names.tsv\n";

if(!$suppress_header){
	my $header=join "\t", (
		"Query_ID",
		"Taxa_ID",
		"Superkingdom",
		"Kingdom",
		"Phylum",
		"Class",
		"Order",
		"Family",
		"Genus",
		"Species"
	);
	print NAMES_RESULTS_FH "# $header\n";
	print ID_RESULTS_FH "# $header\n";
}

###############################################################################

open(FH, "cut -f $taxa_ids_file_columns $taxa_ids_file |") || die "Could not open $taxa_ids_file\n";

my @records=();
my $last_rec="";

# If first column contains #, it's a comment/header line, so skip it.
$_=<FH>;
if(substr($_,0,1) eq "#"){
	$_=<FH>;
}

# Main processing loop
do{
	chomp;

	my ($rec_id, $taxa_id)=split "\t", $_;

        my $taxa_levels_arr_ref=get_higher_levels($taxa_id, $nodes_ref, $levels_ref);

	# Output Ids
	my $outstr=join "\t", @{$taxa_levels_arr_ref};
	print ID_RESULTS_FH "$rec_id\t$taxa_id\t$outstr\n";	

	# Output Names
	my @names;
	my $prev_name="Unknown";

	foreach my $id(@{$taxa_levels_arr_ref}){
		my $cur_name=${$names_ref}{$id};
		if($cur_name eq "NULL"){
			$cur_name=$prev_name;
		}else{
			$prev_name="($cur_name)";
		}
		push @names, $cur_name;
	}
	my $outstr=join "\t", @names;
	print NAMES_RESULTS_FH "$rec_id\t$taxa_id\t$outstr\n";	

}while(<FH>);

close(FH);

close(ID_RESULTS_FH);
close(NAMES_RESULTS_FH);

###############################################################################

print STDERR "Done.\n";

