#!/usr/bin/env perl

use strict;
use Getopt::Std;
use File::Basename;
use Sys::Hostname;
use vars qw ($opt_l $opt_o $opt_k $opt_p);

getopts("l:o:k:p:");

my $KMER_SIZE=30;
my $NUM_KMERS_KEPT=1000;
my $DATABASE_SIZE=1000000000; # 1G

my $usage = "
	usage:
	$0
		-l <fasta file sample list>
		-o <output filename root>
		[-k <kmer size, default=$KMER_SIZE residues>]
		[-p <minimum number of k-mers to keep, default=$NUM_KMERS_KEPT kmers>]	

	This script will take in a list of fasta files, where each
	file contains the reads for a particular biological sample.
	The distribution of k-mers of the specified length (-k) will be identified
	across the samples, then a summary table will be built based on
	top -p proportions across all samples.  

	The (-l) list should be of the format:
		<sample name>\\t<fasta file path>\\n

	For each fasta file, jellyfish is executed to determine the kmers
	across each of the samples.  Kmer counts are then normalized based
	on the total possible kmers found in each sample (to account for differences
	in the input fasta file sizes).  The most frequent kmers across
	all samples are identified, and then the proportions of each of the
	top kmers are output in a summary table.
	
		
";

###############################################################################

if(!defined($opt_l) || !defined($opt_o)){
	die $usage;
}

# Required parameters
my $FastaFileList=$opt_l;
my $OutputFilenameRoot=$opt_o;

# Optional/Default-ready parameters
my $KmerSize=$KMER_SIZE;
my $NumKmersKept=$NUM_KMERS_KEPT;
my $DatabaseSize=$DATABASE_SIZE;

if(defined($opt_k)){
	$KmerSize=$opt_k;
}

if(defined($opt_p)){
	$NumKmersKept=$opt_p;
}

###############################################################################

sub normalize_hash_counts{
	# This subroutine will read in a counts file (kmer \t count \n), then
	# output a normalized proportion file.
	my $in_counts_fname=shift;
	my $out_norm_fname=shift;

	# Count totals
	print STDERR "Counting: $in_counts_fname\n";
	open(IN_FH, "<$in_counts_fname") || die "Could not open $in_counts_fname\n";
	my $sum=0;
	while(<IN_FH>){
		chomp;
		my ($kmer, $count)=split "\t", $_;
		$sum=$sum+$count;
	}
	close(IN_FH);

	print STDERR "Total used for normalization: $sum.\n";	

	# Normalize and output
	print STDERR "Normalizing and Outputing: $out_norm_fname\n";
	open(OUT_FH, ">$out_norm_fname") || die "Could not open $out_norm_fname for writing.\n";
	open(IN_FH, "<$in_counts_fname") || die "Could not open $in_counts_fname again\n";
	while(<IN_FH>){
		chomp;
		my ($kmer, $count)=split "\t", $_;
		my $norm=$count/$sum;
		print OUT_FH "$kmer\t$norm\n";
	}
	close(IN_FH);
	close(OUT_FH);

	return($sum);
}

sub get_kmer_distribution{
	# This subroutine will run jellyfish and extract the counts for the kmers
	my $fasta_fname=shift;
	my $kmer_length=shift;
	my $tmp_fname=shift;

	my %kmer_hash;

	# Temporary file names
	my $stats_fname="$tmp_fname\.stats.txt";
	my $dump_fname="$tmp_fname\.dump.tsv";
	my $hash_fname="$tmp_fname\.hash.bin";
	my $norm_fname="$tmp_fname\.norm.tsv";

	# Jellyfish count command
	my $count_exec_string=
	"jellyfish count --mer-len $kmer_length --size $DatabaseSize --threads=4 " .
	"--output=$hash_fname --stats=$stats_fname --both-strands $fasta_fname";
	print STDERR "*  $count_exec_string\n";
	print STDERR `$count_exec_string`;

	# Confirm that there is only one jellyfish result hash
	if(-e ("$tmp_fname" . "_1")){
		print STDERR "More than one database found...\n";
		die "Error: Please increase the database size (current: $DatabaseSize).\n";
	}
	
	# Jelly fish dump command
	my $dump_exec_string=
	"jellyfish dump --column --tab --output=$dump_fname $hash_fname\_0";
	print STDERR "*  $dump_exec_string\n";
	print STDERR `$dump_exec_string`;

	# Normalize the counts
	print STDERR "Normalizing...\n";
	normalize_hash_counts($dump_fname, $norm_fname);

	# Delete the counts and dump, now that we have the normalized/proportions
	unlink $dump_fname;
	unlink "$hash_fname\_0";

	print STDERR "ok.\n";
	return;
}

sub get_top_kmers{
	# This subroutine will go through all the intermediate normalized kmers
	# and return a list of the top X proportions across all samples
	my $tmp_fname_root=shift;
	my $sample_ids_ref=shift;
	my $num_to_keep=shift;
	
	my $num_samples=$#{$sample_ids_ref}+1;
	my $combined_fname="$tmp_fname_root\.combined.tsv";

	my %hash;	# Keep track of the sum of the normalized counts across all samples

	# Accumulate the frequencies across all the samples
	for(my $i=0; $i<$num_samples; $i++){
		my $cur_samp_id=${$sample_ids_ref}[$i];
		print STDERR "Loading Normalized Counts: $cur_samp_id\n";
		my $fname="$tmp_fname_root\.$cur_samp_id\.norm.tsv";
		open(FH, "<$fname") || die "Could not open normalized counts: $fname\n";
		while(<FH>){
			chomp;
			my ($kmer, $val)=split "\t", $_;	
			if(!defined($hash{$kmer})){
				$hash{$kmer}=$val;
			}else{
				$hash{$kmer}+=$val
			}
		}
		close(FH);
	}

	# Yank out the kmers and find the ones that are the most frequent
	my @kmers=keys %hash;
	my @sorted=sort {$hash{$b} <=> $hash{$a}} @kmers;
	my @tokeep=splice (@sorted, 0, $num_to_keep);

	return(\@tokeep);
}

sub read_file_list{
	# Read in a (sample ID, fasta path) tuple file
	my $fname=shift;
	
	my @sample_ids;
	my @sample_paths;

	open(FH, "<$fname") || die "Could not open $fname for reading.\n";

	while(<FH>){
		chomp;
		my ($sample_id, $sample_path)=split "\t", $_;
		push @sample_ids, $sample_id;
		push @sample_paths, $sample_path;
	}

	return(\@sample_ids, \@sample_paths, ($#sample_ids+1));
}

sub dump_kmer_list{
	# Output a list of Kmers
	my $fname=shift;
	my $kmerlist_ref=shift;
	open(FH, ">$fname") || die "Could not open $fname.\n";
	for(my $i=0; $i<=$#{$kmerlist_ref}; $i++){
		print FH "${$kmerlist_ref}[$i]\n";
	}	
	close(FH);
}

sub extract_to_summary_table{
	# Across all the samples, build a summary table based on the specified kmer list
	my $summary_table_fname=shift;
	my $keep_kmer_list_ref=shift;
	my $tmp_fname_root=shift;
	my $sample_ids_ref=shift;

	my $num_samples=$#{$sample_ids_ref}+1;

	open(FH, ">$summary_table_fname") || die "Could not open $summary_table_fname for writing.\n";
	
	# Output Header
	print FH "sample_id\ttotal\t";
	for(my $i=0; $i<=$#{$keep_kmer_list_ref}; $i++){
		print FH "${$keep_kmer_list_ref}[$i]\t";
	}
	print FH "Remaining\n";

	my %hash;
	my @kept_kmers=@{$keep_kmer_list_ref};

	# Output per sample counts
	for(my $i=0; $i<$num_samples; $i++){
		my $cur_samp_id=${$sample_ids_ref}[$i];

		# initialize hash
		foreach my $kmer (@kept_kmers){
			$hash{$kmer}=0;
		}
		my $remaining=0;
		my $total=0;

		# Read in counts
		open(IN_NORM, "<$tmp_fname_root\.$cur_samp_id\.norm.tsv") || 
			die "Could not open $tmp_fname_root\.$cur_samp_id\.norm.tsv";
	
		while(<IN_NORM>){
			chomp;
			my ($kmer, $val)=split "\t", $_;
			$total+=$val;
			if(defined($hash{$kmer})){
				$hash{$kmer}=$val;
			}else{
				$remaining+=$val;
			}
		}	
		close(IN_NORM);

		# Output line for sample
		print FH "$cur_samp_id\t$total\t";
		foreach my $kmer(@kept_kmers){
			print FH "$hash{$kmer}\t";
		}
		print FH "$remaining\n";

	}

	close(FH);
}

###############################################################################

# Load the target fasta file list
my ($sample_ids_ref, $sample_paths_ref, $num_samples)=read_file_list($FastaFileList);

# Make a temporary filename
my $host_id=hostname();
my $proc_id=$$;
my ($progname, $path)=fileparse($0);
$progname=~s/\.pl$//;
my $tmp_fname="$OutputFilenameRoot\.$host_id\.$proc_id\.$progname";
print STDERR "Temporary file name root: $tmp_fname\n";

# Get the distribution of kmers across al the samples
for(my $i=0; $i<$num_samples; $i++){
	print STDERR "\n";
	my ($cur_sample, $cur_path)=(${$sample_ids_ref}[$i], ${$sample_paths_ref}[$i]);
	print STDERR "Working on: $cur_sample @ $cur_path\n";
	get_kmer_distribution($cur_path, $KmerSize, "$tmp_fname.$cur_sample");
}
print STDERR "\n";

# Identify the top Kmers
print STDERR "Identifying top K-mers across all samples...\n";
my $kept_kmer_ref=get_top_kmers($tmp_fname, $sample_ids_ref, $NumKmersKept);
dump_kmer_list("$tmp_fname.kmers_kept.txt", $kept_kmer_ref);

# Output the summary table with the top kmers
print STDERR "Outputing Summary Table...\n";
extract_to_summary_table("$OutputFilenameRoot\.summary_table.tsv", $kept_kmer_ref, $tmp_fname, $sample_ids_ref);

###############################################################################

print STDERR "done.\n";


