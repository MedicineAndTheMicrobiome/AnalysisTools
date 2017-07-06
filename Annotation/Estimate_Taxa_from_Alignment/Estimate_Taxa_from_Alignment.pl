#!/usr/bin/env perl

use strict;
use Getopt::Std;
use vars qw ($opt_a $opt_A $opt_t $opt_n $opt_o $opt_h);

getopts("a:A:t:n:o:h");

my $usage = "
	usage:
	$0
	-a <alignment file (read_id, comp_ident, taxa_id) >
	-A <alignment columns for read_id, comp_ident and taxa_id, respectively, eg. \"1,2,6\">
	-t <taxa nodes file (parent_id, child_id) >
	-n <taxa names file (taxa_id, taxa_name) >
	-o <output file>
	[-h (suppress header flag)]

	The alignment file (-a) needs to contain:
		1.) read id
		2.) composite identity
		3.) taxa id	

	The taxa nodes (-t) :
		1.) child taxa id
		2.) parent taxa id

	The taxa names (-n):
		1.) taxa id
		2.) taxa name

	Output looks like (-o):
		1.) read id
		2.) taxa names: best guess ancestor, 90, 75, 60, 45 
		3.) taxa ids: best guess ancestor, 90, 75, 60, 45

";

###############################################################################

if(!defined($opt_a) || !defined($opt_A) || !defined($opt_t) || !defined($opt_n) || !defined($opt_o)){
	die $usage;
}

my $alignments_file=$opt_a;
my $alignments_file_columns=$opt_A;
my $taxa_nodes_file=$opt_t;
my $taxa_names_file=$opt_n;
my $output_file=$opt_o;
my $suppress_header=($opt_h eq "1");

###############################################################################

# Columns in alignment file
my $READ_ID_COL=0;
my $COMPOSIT_ID_COL=1;
my $TAXA_ID=2;


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
				last;
			}else{
				$taxa_id=$parent_taxa_id;
			}

		}while($taxa_id ne "1" && defined($taxa_id));
		unshift @ancestry, "1";

		#print STDERR (join ";", @ancestry) . "\n";
		$ancestry_cache{$in_taxa_id}=\@ancestry;
	}

	return($ancestry_cache{$in_taxa_id});
}

###############################################################################

# Compare an array of ancentry arrays
sub find_highest_common_ancestor{
	my $anc_str_arr_ref=shift;
	
	my @consensus_id;
	my $num_ids=$#{$anc_str_arr_ref}+1;

	if($num_ids==0){
		return("NULL");
	}

	if($num_ids==1){
		my $only_rec=${$anc_str_arr_ref}[0];
		my $lowest_taxa_idx=$#{$only_rec};
		my $lowest_taxa=${$only_rec}[$lowest_taxa_idx];
		return($lowest_taxa);
	}

	# Print input for debugging
	#if(0){
	#	for(my $i=0; $i<=$#{$anc_str_arr_ref}; $i++){
	#		print STDERR (join ";", @{${$anc_str_arr_ref}[$i]}) . "\n";
	#	}
	#}

	my $level_id=0;
	my $all_agree;
	do{
		$all_agree=1;
		my $first_taxa_id=${$anc_str_arr_ref}[0][$level_id];
		for(my $i=1; $i<$num_ids; $i++){
			if(!defined($first_taxa_id) || $first_taxa_id eq "" || 
				($first_taxa_id ne ${$anc_str_arr_ref}[$i][$level_id])){
				$all_agree=0;
			}
		}

		$level_id++;
	}while($all_agree && $level_id<50);

	my $last_shared_id=$level_id-2;
	if($last_shared_id==-1){
		return(0);
	}

	return(${$anc_str_arr_ref}[0][$last_shared_id]);
}


###############################################################################

my $nodes_ref;
my $names_ref;
$nodes_ref=load_taxa_nodes($taxa_nodes_file);
$names_ref=load_taxa_names($taxa_names_file);

if(0){
	print STDERR "\n";
	my $anc_arr1=build_ancestry_arr("61853", $nodes_ref);
	print STDERR anc_arr_to_str($anc_arr1, $names_ref) . "\n\n";
	my $anc_arr2=build_ancestry_arr("6282", $nodes_ref);
	print STDERR anc_arr_to_str($anc_arr2, $names_ref) . "\n\n";
	my $anc_arr3=build_ancestry_arr("408170", $nodes_ref);
	print STDERR anc_arr_to_str($anc_arr3, $names_ref) . "\n\n";

	my @id_arr;
	$id_arr[0]=$anc_arr3;
	$id_arr[1]=$anc_arr2;
	#$id_arr[2]=$anc_arr1;

	my $hca=find_highest_common_ancestor(\@id_arr);

	print STDERR "Common Ancestor: ${$names_ref}{$hca}\n";

	exit;
}

###############################################################################

open(RESULTS_FH, ">$output_file") || die "Could not open $output_file.\n";
open(EXCLUDED_HITS_FH, ">$output_file\.excluded") || die "Could not open $output_file\.excluded.\n";

if(!$suppress_header){
	print STDERR "Header being written.\n";
	my $fields_str;
	$fields_str=join "\t", (
		"Query_ID",
		"Best_Guess_TaxaName",
		"TaxaName90",	
		"TaxaName75",	
		"TaxaName60",	
		"TaxaName45",	
		"Best_Guess_TaxaID",
		"TaxaID90",
		"TaxaID75",
		"TaxaID60",
		"TaxaID45",
		"NumHits90",
		"NumHits75",
		"NumHits60",
		"NumHits45"
	);
	print RESULTS_FH "# $fields_str\n";
}else{
	print STDERR "Header for output file has been suppressed.\n";
}

sub print_arr{
	my $arr_ref=shift;
	for(my $i=0; $i<=$#{$arr_ref}; $i++){
		#print RESULTS_FH "\t${$arr_ref}[$i]\n";
		print RESULTS_FH anc_arr_to_str(${$arr_ref}[$i], $names_ref) . "\n\n";
	}
}

my %no_rank_hash;
%no_rank_hash=(
	"12908" => 1,
	"32644" => 1,
	"37965" => 1,
	"151659" => 1,
	"408169" => 1,
	"684672" => 1,
	"704107" => 1,
	"1306155" => 1,
	"1427524" => 1,
	"1515699" => 1
);



sub process_records{
	my $rec_arr_ref=shift;
	my $nodes_ref=shift;

	my @perc_id_arr;
	my @taxa_id_arr;

	my @p90_arr;
	my @p75_arr;
	my @p60_arr;
	my @p45_arr;

	my @p90_norank_arr;
	my @p75_norank_arr;
	my @p60_norank_arr;
	my @p45_norank_arr;

	if($#{$rec_arr_ref}==-1){
		# Skip if record is empty.
		return;
	}

	my ($rec_id, $perc_id, $taxa_id);

	my $num_records=0;

	foreach my $rec(@{$rec_arr_ref}){
	
		($rec_id, $perc_id, $taxa_id)=split "\t", $rec;

		my $ancestry_arr_ref=build_ancestry_arr($taxa_id, $nodes_ref);

		if($no_rank_hash{${$ancestry_arr_ref}[1]}){
			print EXCLUDED_HITS_FH "$rec_id / $perc_id / ${$names_ref}{$taxa_id} Removed.\n";
			if($perc_id >= .90){
				push @p90_norank_arr, $ancestry_arr_ref;
			}
			if($perc_id >= .75){
				push @p75_norank_arr, $ancestry_arr_ref;
			}
			if($perc_id >= .60){
				push @p60_norank_arr, $ancestry_arr_ref;
			}
			if($perc_id >= .45){
				push @p45_norank_arr, $ancestry_arr_ref;
			}
		}else{

			if($perc_id >= .90){
				push @p90_arr, $ancestry_arr_ref;
			}
			if($perc_id >= .75){
				push @p75_arr, $ancestry_arr_ref;
			}
			if($perc_id >= .60){
				push @p60_arr, $ancestry_arr_ref;
			}
			if($perc_id >= .45){
				push @p45_arr, $ancestry_arr_ref;
			}

		}

		$num_records++;
	}

	# Use no rank taxa assignments if good assignments don't exist
	if($#p90_arr==-1){
		@p90_arr=@p90_norank_arr;
	}
	if($#p75_arr==-1){
		@p75_arr=@p75_norank_arr;
	}
	if($#p60_arr==-1){
		@p60_arr=@p60_norank_arr;
	}
	if($#p45_arr==-1){
		@p45_arr=@p45_norank_arr;
	}

	# Find HCA for each subset
	my $hca90=find_highest_common_ancestor(\@p90_arr);
	my $hca75=find_highest_common_ancestor(\@p75_arr);
	my $hca60=find_highest_common_ancestor(\@p60_arr);
	my $hca45=find_highest_common_ancestor(\@p45_arr);

	# Translate HCA taxa id to name
	my $name90=${$names_ref}{$hca90};
	my $name75=${$names_ref}{$hca75};
	my $name60=${$names_ref}{$hca60};
	my $name45=${$names_ref}{$hca45};

	# Identify best guess
	my @results=($hca90, $hca75, $hca60, $hca45);
	my @suffix=("_90", "_75", "_60", "_45");

	my $i=0;
	for($i=0; $i<=$#results; $i++){
		if($results[$i] ne "NULL"){
			last;	
		}
	}

	my $best_guess_id=$results[$i];
	my $best_guess_name=${$names_ref}{$best_guess_id} . $suffix[$i];

	# Output HCA
	my $outstr=join "\t", (

		$rec_id,

		$best_guess_name,
		$name90 . "_90", 
		$name75 . "_75", 
		$name60 . "_60", 
		$name45 . "_45",

		$best_guess_id,
		$hca90, $hca75, $hca60, $hca45,

		$#p90_arr+1,
		$#p75_arr+1,
		$#p60_arr+1,
		$#p45_arr+1,
	);

	print RESULTS_FH "$outstr\n";

}

###############################################################################
my $ID_COL=0;
my $PERC_COL=1;
my $TAXA_COL=2;

open(FH, "cut -f $alignments_file_columns $alignments_file | ") || die "Could not open $alignments_file\n";

my @records=();
my $last_rec="";
while(<FH>){
	chomp;

	my @fields=split "\t", $_;
	
	if(($fields[$ID_COL] ne $last_rec)){
		process_records(\@records, $nodes_ref);
		@records=();
	}

	$last_rec=$fields[$ID_COL];
	
	my $taxa_id=$fields[$TAXA_COL];
	my $taxa_is_valid=defined(${$nodes_ref}{$taxa_id});

	if($fields[$TAXA_COL] eq ""){
		print STDERR "Warning: $fields[$ID_COL] at $fields[$PERC_COL] is NULL.\n";
	}elsif(!$taxa_is_valid){
		print STDERR "Warning: $fields[$ID_COL] at $fields[$PERC_COL] has invalid taxa id: $taxa_id\n";
	}else{
		push @records, (join "\t", ($fields[$ID_COL], $fields[$PERC_COL], $taxa_id));
	}

}
process_records(\@records, $nodes_ref);

close(FH);

###############################################################################

print STDERR "Done.\n";

