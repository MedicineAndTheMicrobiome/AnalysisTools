#!/usr/bin/env perl

use strict;
use Getopt::Std;
use vars qw($opt_t $opt_n $opt_g $opt_o);
getopts("t:n:g:o:");

my $usage = "usage:
$0
	-t <Mothur taxonomy file (read-to-taxonomy map)>
	-n <Mother names file (representative-to-represented read map)>
	-g <Mother group file (read-to-sample map)>
	-o <output file name>

	This script will read in the read taxonomy assignments
	and automatically look up which sample the group belongs
	to and then increment the taxonomy assignment for that group.
	
	When completed, it will output a summary file table for
	each of the 6 taxonomic levels from Kingdom down to Genus.

";

if(!(
	defined($opt_t) && 
	defined($opt_n) && 
	defined($opt_g) && 
	defined($opt_o)
)){
	die $usage;
}

my $taxonomy_file=$opt_t;
my $names_file=$opt_n;
my $group_file=$opt_g;
my $output_file=$opt_o;

##################################################################################

print STDERR "Taxonomy File: $taxonomy_file\n";
print STDERR "Names File: $names_file\n";
print STDERR "Groups File: $group_file\n";
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

#------------------------------------------------------------------------------

sub load_names{
	my $names_file=shift;
	my %read_to_rep_hash;
	open(FH, "<$names_file") || die "Could not open read-to_represented file: $names_file\n";

	my $num_reps=0;
	my $num_reads=0;

	while(<FH>){
		chomp;
		my ($repres, $reads)=split "\t", $_;
		my @reads_arr=split ",", $reads;
		foreach my $read_id (@reads_arr){
			$read_to_rep_hash{$read_id}=$repres;
			$num_reads++;
		}
		$num_reps++;
	}

	close(FH);
	
	print STDERR "Num Representatives: $num_reps\n";
	print STDERR "Num Reads Represented: $num_reads\n";

	return(\%read_to_rep_hash);
	
}

##################################################################################

my $read_to_taxa_hash_ref=load_read_to_taxa_map($taxonomy_file);
my $read_to_rep_hash_ref=load_names($names_file);

my %sample_to_taxa_arr_hash;

open(FH, "<$group_file") || die "Could not open read-to-group file: $group_file\n";
while(<FH>){
	chomp;
	my ($read_id, $group_id)=split "\t", $_;

	# Create a new hash entry if we haven't seen this group before
	if(!defined($sample_to_taxa_arr_hash{$group_id})){
		my @arr;
		$sample_to_taxa_arr_hash{$group_id}=\@arr;
	}
	
	my $rep=${$read_to_rep_hash_ref}{$read_id};

	# If no rep, assume it's by it's its own representative
	if(!defined($rep)){
		$rep=$read_id;
	}

	my $taxa=${$read_to_taxa_hash_ref}{$rep};

	push @{$sample_to_taxa_arr_hash{$group_id}}, $taxa;

	#print STDERR "$read_id ($group_id) --> ", (join ";", @{$taxa}), "\n";
}
close(FH);

##################################################################################

# Relase memory
%{$read_to_taxa_hash_ref}=();
%{$read_to_rep_hash_ref}=();

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
