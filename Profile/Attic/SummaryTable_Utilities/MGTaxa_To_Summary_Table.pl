#!/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_o $opt_t);

getopts("i:o:t:");

my $usage = "
usage:
	$0
		-i <MG-Taxa Input File>
		-o <output summary file name>
		[-t <exclude taxon id list file>]

	Builds a summary table entry based on an MG-Taxa
	output file.  Taxon id list specifies which taxa
	will be excluded. 

";

if(!defined($opt_i) || !defined($opt_o)){
	die $usage;
}

my $InputFile=$opt_i;
my $OutputFile=$opt_o;

my $TaxaExclusion;
if(defined($opt_t)){
	$TaxaExclusion=$opt_t
}

###############################################################################

my %taxa_exclusion;
if($TaxaExclusion ne ""){

	print STDERR "Loading taxa exclusion file.\n";

	open(MAP_FH, "<$TaxaExclusion") || die "Could not open $TaxaExclusion\n";
	while(<MAP_FH>){
		chomp;
		my $key=$_;
		$key=~s/ $//g;
		$key=~s/^ //g;
		$taxa_exclusion{$key}=1;
	}
	close(MAP_FH);
	print STDERR "Done.\n\n";

	#---------------------------------------------------------------------

	print STDERR "Taxa to exclude:\n";
	foreach my $key (keys %taxa_exclusion){
		print STDERR "$key\n";
	}
	print STDERR "\n";
}

###############################################################################

open(IN_FH, "<$InputFile") || die "Could not open map file $InputFile\n";

my $STRAIN_IDX=2;
my $SPECIES_IDX=6;
my $GENUS_IDX=9;
my $FAMILY_IDX=12;
my $ORDER_IDX=15;
my $CLASS_IDX=18;
my $PHYLUM_IDX=21;
my $SUPERKINGDOM_IDX=24;

print STDERR "Checking header line of $InputFile.\n";
my $header_line=<IN_FH>;
my @line=split /\t/, $header_line;
if(
	$line[$STRAIN_IDX] ne "taxid" ||
	$line[$SPECIES_IDX] ne "taxid_species" ||
	$line[$GENUS_IDX] ne "taxid_genus" ||
	$line[$FAMILY_IDX] ne "taxid_family" ||
	$line[$ORDER_IDX] ne "taxid_order" ||
	$line[$CLASS_IDX] ne "taxid_class" ||
	$line[$PHYLUM_IDX] ne "taxid_phylum" ||
	$line[$SUPERKINGDOM_IDX] ne "taxid_superkingdom" 
){
	print STDERR "Error:  Unexpected header info:\n";
	print STDERR "$header_line";
	die;
}

my $NUM_CATEGORIES=8;

###############################################################################

my %taxa_counts_hash;
my %taxa_names_hash;

print STDERR "Loading in counts.\n";
while(<IN_FH>){
	chomp;
	@line=split /\t/, $_;

	# Pick out the columns we want
	my @taxa_ids=(
		$line[$STRAIN_IDX],
		$line[$SPECIES_IDX],
		$line[$GENUS_IDX],
		$line[$FAMILY_IDX],
		$line[$ORDER_IDX],
		$line[$CLASS_IDX],
		$line[$PHYLUM_IDX],
		$line[$SUPERKINGDOM_IDX]);

	my @taxa_names=(
		$line[$STRAIN_IDX+1],
		$line[$SPECIES_IDX+1],
		$line[$GENUS_IDX+1],
		$line[$FAMILY_IDX+1],
		$line[$ORDER_IDX+1],
		$line[$CLASS_IDX+1],
		$line[$PHYLUM_IDX+1],
		$line[$SUPERKINGDOM_IDX+1]);


	# Exclude entire line if any of the taxa_id are in the list
	my $taxa_id;
	my $skip=0;
	for(my $i=0; $i<$NUM_CATEGORIES; $i++){
		$taxa_id=$taxa_ids[$i];
		if(defined($taxa_exclusion{$taxa_id})){
			$skip=1;
			last;	
		}
		$taxa_names_hash{$taxa_id}=$taxa_names[$i];
	}

	# Save up counts 
	if(!$skip){
		for(my $i=0; $i<$NUM_CATEGORIES; $i++){
			my $cat=$taxa_ids[$i];
			#print STDERR "Level: $i / Category: $cat\n";
			if(!defined($taxa_counts_hash{$i}{$cat})){
				$taxa_counts_hash{$i}{$cat}=1;
			}else{
				$taxa_counts_hash{$i}{$cat}++;
			}
		}
	}
}

close(IN_FH);

###############################################################################

print STDERR "Outputing counts.\n";
my @category_levels=("strain", "species", "genus", "family", "order", "class", "phylum", "superkingdom");
for(my $i=0; $i<$NUM_CATEGORIES; $i++){

	# Get Taxa IDs and put them in numerical order
	my @taxa=sort {$a <=> $b} keys %{$taxa_counts_hash{$i}};

	my @taxa_names;
	for(my $j=0; $j<=$#taxa; $j++){
		$taxa_names[$j]=$taxa_names_hash{$taxa[$j]};
	}

	# Open new file
	my $output_fn="$OutputFile\.$category_levels[$i]\.summary_table.xls";
	open(OUT_FH, ">$output_fn") || die "Could not open $output_fn\n";

	# Output header info
	my $taxa_names_str=join "\t", @taxa_names;	
	print OUT_FH "sample_id\ttotal\t$taxa_names_str\n";

	# Output counts
	my $total=0;
	my @counts;
	foreach my $taxon(@taxa){
		#print STDERR "Working on $taxa\n";
		my $cat_count=$taxa_counts_hash{$i}{$taxon};
		$total+=$cat_count;
		push @counts, $cat_count;
	}
	my $outline=join "\t", @counts;
	print OUT_FH "$OutputFile\t$total\t$outline\n";

	# Close the file
	close(OUT_FH);
}

###############################################################################

print STDERR "Done.\n\n";
