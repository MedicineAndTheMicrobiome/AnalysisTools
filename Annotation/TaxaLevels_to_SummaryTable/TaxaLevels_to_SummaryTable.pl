#!/usr/bin/env perl

use strict;
use Getopt::Std;
use FileHandle;
use vars qw ($opt_l $opt_o);

getopts("l:o:");

my $usage = "
usage:
    $0

	-l <List of Taxa annotation output files>
	-o <Output filename root>
	
	This script will read in a list of taxa annotation files.
	The files have been processed by Get_Higher_Levels.pl, so
	that each read contains a taxonomic assignment from Domain,
	all the way down to Species.

	The list (-l) file should look like:

		<sample name>, <taxa annotation path>

	The output file name should be the root name of the summary
	tables to generate.

	For example, 
		<output filename root>.domain.summary_table.tsv
		<output filename root>.kingdom.summary_table.tsv
		<output filename root>.phylum.summary_table.tsv
		<output filename root>.class.summary_table.tsv
		<output filename root>.order.summary_table.tsv
		<output filename root>.family.summary_table.tsv
		<output filename root>.genus.summary_table.tsv
		<output filename root>.species.summary_table.tsv

";

###############################################################################

if(
	!defined($opt_l) ||
	!defined($opt_o)
){
	die $usage;
}

my $TaxaFileList=$opt_l;
my $OutputFnameRoot=$opt_o;

###############################################################################

sub read_taxafile_list{
	my $fname=shift;
	
	print STDERR "Loading Taxa File List: $fname\n";

	open(TFH, "<$fname") || die "Could not open $fname.\n";

	my %fn_hash;

	while(<TFH>){
		chomp;
		my ($sample_name, $path)=split "\t", $_;
		$fn_hash{$sample_name}=$path;	
	}

	close(TFH);
	
	return(\%fn_hash);
}

###############################################################################

my $DOMAIN_COL=2;
my $SPECIES_COL=9;
my $TOT_LEVELS=$SPECIES_COL-$DOMAIN_COL+1;

my @level_names=("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species");

sub get_counts_from_taxa_file{
	my $fname=shift;

	my @level_hash_arr;
	
	print STDERR "$fname\n";
	open(TAX_FH, "<$fname") || die "Could not open $fname\n";


	for(my $i=0; $i<$TOT_LEVELS; $i++){
		my $hash;
		$level_hash_arr[$i]=$hash;
	}

	while(<TAX_FH>){
		chomp;
		if($_=~/^#/){
			# Skip comments
			next;
		}
		my @col=split "\t", $_;
		for(my $level=$DOMAIN_COL; $level<=$SPECIES_COL; $level++){
			if(!defined(${$level_hash_arr[$level-$DOMAIN_COL]}{$col[$level]})){
				${$level_hash_arr[$level-$DOMAIN_COL]}{$col[$level]}=1;
			}else{
				${$level_hash_arr[$level-$DOMAIN_COL]}{$col[$level]}++;
			}
		}
	}

	close(TAX_FH);

	print STDERR "Number of categories detected at each level:\n";
	for(my $i=0; $i<$TOT_LEVELS; $i++){
		my $num_categories=keys %{$level_hash_arr[$i]};
		print STDERR "\t$level_names[$i]:  \t$num_categories\n";
	}
	print STDERR "\n\n";

	return(\@level_hash_arr);

}


###############################################################################

my $target_files_hash_ref=read_taxafile_list($TaxaFileList);

my @sample_ids=sort keys %{$target_files_hash_ref};
my %level_counts_by_sample;
foreach my $sample_id(@sample_ids){
	print STDERR "Working on: $sample_id\n";
	my $path=${$target_files_hash_ref}{$sample_id};

	print STDERR "Path: $path\n";

	$level_counts_by_sample{$sample_id}=get_counts_from_taxa_file($path);
}

for(my $lev_ix=0; $lev_ix<$TOT_LEVELS; $lev_ix++){
	print STDERR "Consolidating $level_names[$lev_ix]...\n";

	my %master_hash;

	foreach my $sample_id(@sample_ids){
		my $samp_hash_ref=${$level_counts_by_sample{$sample_id}}[$lev_ix];

		foreach my $key (keys %{$samp_hash_ref}){
			$master_hash{$key}=1;
		}
	}

	my @masterkeys_arr=sort keys %master_hash;
	foreach my $masterkeys (@masterkeys_arr){
		print STDERR "$masterkeys\n";
	}
	print STDERR "\n\n";

	my $outputfname="$OutputFnameRoot.$level_names[$lev_ix].summary_table.tsv";
	open(OUT_FH, ">$outputfname") || die "Coudl not open $outputfname\n";

	my $category_str=join "\t", @masterkeys_arr;
	print OUT_FH "Sample_ID\tTotal\t$category_str\n";
	foreach my $sample_id(@sample_ids){

		my $samp_hash_ref=${$level_counts_by_sample{$sample_id}}[$lev_ix];

		my $total=0;
		my @counts;
		foreach my $masterkeys (@masterkeys_arr){
			my $val=${$samp_hash_ref}{$masterkeys};
			if(!defined($val)){
				$val=0;
			}

			$total+=$val;
			push @counts, $val;
		}

		my $val_str=join "\t", @counts;
		print OUT_FH "$sample_id\t$total\t$val_str\n";
	}


}

print STDERR "Done.\n";

