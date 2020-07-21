#!/usr/bin/env perl

use strict;
use Getopt::Std;
use vars qw($opt_t $opt_b $opt_o);
getopts("t:b:o:");

my $usage = "usage:
$0
	-t <Mothur taxonomy file (read-to-taxonomy map)>
	-b <Mother count table file (read-to-sample map)>
	-o <output file name>

	This script will read in the read taxonomy assignments
	and automatically look up which sample the group belongs
	to and then increment the taxonomy assignment for that group.
	
	When completed, it will output a summary file table for
	each of the 6 taxonomic levels from Kingdom down to Genus.

";

if(!(
	defined($opt_t) && 
	defined($opt_b) && 
	defined($opt_o)
)){
	die $usage;
}

my $taxonomy_file=$opt_t;
my $counttable_file=$opt_b;
my $output_file=$opt_o;

##################################################################################

print STDERR "Taxonomy File: $taxonomy_file\n";
print STDERR "Count Table File: $counttable_file\n";
print STDERR "Output File Name Root: $output_file\n";

##################################################################################

sub load_read_to_taxa_map{
	my $taxa_file=shift;
	
	my %read_to_taxa_hash;

	open(FH, "<$taxa_file") || die "Could not open read-to-taxa file: $taxa_file\n";	
	
	while(<FH>){
		chomp;
		my ($read_id, $assignment)=split "\t", $_;
		
		my @levels=split ";", $assignment;
		my @clean_assignment_arr=();
		my $previous="unknown";
		my $clean_name;
		foreach my $level(@levels){
			if($level eq "unclassified"){
				$clean_name="($previous)";	
			}elsif($level eq "unknown"){
				$clean_name="unknown";
			}else{
				if($level=~/(.+)(\(\d+\))$/){
					$clean_name=$1;
					my $conf=$2;
					$previous=$clean_name;
				}else{
					die "Could not parse taxa string: $assignment.\n";
				}
			}
			push @clean_assignment_arr, $clean_name;
		}

		$read_to_taxa_hash{$read_id}=\@clean_assignment_arr;

		#print "$assignment\n";
		#print (join ";", @clean_assignment_arr) . "\n";
		#print "\n\n";
	}

	close(FH);

	return(\%read_to_taxa_hash);
}

##################################################################################

my $read_to_taxa_hash_ref=load_read_to_taxa_map($taxonomy_file);

my %sample_to_taxa_arr_hash;
my $line_num=0;
my @id_to_samp_ids_arr;

open(FH, "<$counttable_file") || die "Could not open Count Table file: $counttable_file\n";

while(<FH>){

	chomp;
	my $inline=$_;

	if($inline=~/^#/){
		print STDERR "$inline\n";
	}elsif($line_num==0){

		my @col_headers=split "\t", $inline;
		
		my $rep_seq=shift @col_headers;
		my $total=shift @col_headers;

		if($rep_seq ne "Representative_Sequence"){
			print STDERR "ERROR reading column headers: $rep_seq ne Representative_Sequence.\n";
			exit;
		}
		if($total ne "total"){
			print STDERR "ERROR reading column headers: $total ne total.\n";
			exit;
		}
		@id_to_samp_ids_arr=@col_headers;
	
		$line_num++;

	}else{

		my @columns=split "\t", $inline;

		my $rep_seq=shift @columns;
		my $total=shift @columns;

		my $num_samples_to_map=$#columns+1;

		for(my $i=0; $i<$num_samples_to_map; $i++){
			my $redund_info=$columns[$i];
			my ($samp_id, $counts)=split ",", $redund_info;
			my $group_id=$id_to_samp_ids_arr[$samp_id-1];

			my $taxa=${$read_to_taxa_hash_ref}{$rep_seq};

			for(my $times=0; $times<$counts; $times++){
				push @{$sample_to_taxa_arr_hash{$group_id}}, $taxa;
			}
			
		}
		$line_num++;
	}
}
close(FH);

##################################################################################

# Relase memory
%{$read_to_taxa_hash_ref}=();

##################################################################################

sub output_summary_table{
	my $output_filename_root=shift;
	my $taxonomic_level_tag=shift;
	my $counts_hash_ref=shift;
	
	my @group_keys=keys %{$counts_hash_ref};
	my %all_taxa;

	print STDERR "Working on summary table for $taxonomic_level_tag\n";

	# Get all the unique taxa across all samples and per sample totals
	my %group_totals;

	foreach my $group(@group_keys){
		print STDERR "\nSample ID: $group\n";
		my @taxa_keys=keys ${$counts_hash_ref}{$group};
		foreach my $taxa(@taxa_keys){
			$all_taxa{$taxa}=1;
			print STDERR "${$counts_hash_ref}{$group}{$taxa}\t$taxa\n";
			$group_totals{$group}+=${$counts_hash_ref}{$group}{$taxa};
		}
		print STDERR "Total: $group_totals{$group}\n";
	}

	my @all_taxa_arr=sort keys(%all_taxa);

	# Open output file and start writing
	my $fname="$output_filename_root.$taxonomic_level_tag.summary_table.tsv";
	open(FH, ">$fname") || die "Could not open $fname for writing.\n";

	my @out_header=("Sample_ID", "Total", @all_taxa_arr);
	print FH (join "\t", @out_header) . "\n";

	foreach my $group(sort @group_keys){
		my @out_arr=($group, $group_totals{$group});
		foreach my $taxa(@all_taxa_arr){
			my $val=${$counts_hash_ref}{$group}{$taxa};
			if(!defined($val)){
				$val=0;
			}
			push @out_arr, $val;	
		}
		print FH (join "\t", @out_arr) . "\n";
	}

	return;
}

##################################################################################

my @taxa_levels=("kingdom","phylum","class","order","family","genus");
my @group_id_arr=keys %sample_to_taxa_arr_hash;

for(my $lidx=0; $lidx<6; $lidx++){
	
	print STDERR "\n\nWorking on $taxa_levels[$lidx].\n";
	my %counts;
	foreach my $group_id(@group_id_arr){
		foreach my $taxa_arr_ref(@{$sample_to_taxa_arr_hash{$group_id}}){
	
			# Construct taxa name based on taxa level
			my $taxa_name=join ";", (@{$taxa_arr_ref}[0 .. $lidx]);
			#print STDERR "$group_id [$lidx]: $taxa_name\n";
			$counts{$group_id}{$taxa_name}++;
		}
	}

	output_summary_table($output_file, $taxa_levels[$lidx], \%counts);

}

##################################################################################

print STDERR "Done.\n\n";
