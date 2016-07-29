#!/usr/bin/env perl

use strict;
use Getopt::Std;
use File::Basename;
use Sys::Hostname;
use vars qw($opt_f $opt_r $opt_o $opt_t);

getopts("f:r:o:t:");
my $usage = "
usage: 
$0 
	-f <forward fastq file>
	-r <reverse fastq file>
	-o <output file name root>
	[-t <temporary directory>]

	This script will read in two fastq files,
	identify which read ids are shared between then
	and merge them.

";

if(!(
	defined($opt_f) &&
	defined($opt_r) &&
	defined($opt_o))
){
	die $usage;
}

my $for_name=$opt_f;
my $rev_name=$opt_r;
my $output_name=$opt_o;
my $temp_directory=$opt_t;

my ($output_fname, $output_path)=fileparse($output_name);

if(!defined($temp_directory)){
	$temp_directory=$output_path;
	if($temp_directory eq ""){
		$temp_directory eq "."
	}
}

print STDERR "Temporary Directory: $temp_directory\n";

###############################################################################

sub get_ids{

	my $fname=shift;
	my $hash_ref=shift;

	print STDERR "Getting IDs for $fname.\n";

	if($fname=~/\.gz$/){
		open(FASTQ_FH, "zcat $fname | ") || die "Could not open $fname.\n";
	}else{
		open(FASTQ_FH, "<$fname") || die "Could not open $fname.\n";
	}

	my $num_ids_found=0;
	while(!eof(FASTQ_FH)){

		my $defline=<FASTQ_FH>;
		my $seq=<FASTQ_FH>;
		my $plus=<FASTQ_FH>;
		my $qv=<FASTQ_FH>;

		if($defline=~/^@(\S+)/){
		
			my $id=$1;

			# for SOLEXA reads
			if($id=~/(.+)\/\d$/){
				$id=$1;
			}

			if(!defined(${$hash_ref}{$id})){
				${$hash_ref}{$id}=1;		
			}else{
				${$hash_ref}{$id}++;
			}
			$num_ids_found++;

		}else{
			die "Error: Could not parse ID out of defline: $defline\n";
		}

	}

	print STDERR "Num IDs found: $num_ids_found\n";

	return;
}

###############################################################################

sub extract_ids{

	my $fname=shift;
	my $hash_ref=shift;
	my $out_fn=shift;
	my $remainder_fn=shift;

	print STDERR "Extracting Sequences from $fname.\n";

	if($fname=~/\.gz$/){
		open(FASTQ_FH, "zcat $fname | ") || die "Could not open $fname.\n";
	}else{
		open(FASTQ_FH, "<$fname") || die "Could not open $fname.\n";
	}

	open(OUT_FH, ">$out_fn") || die "Could not open $out_fn for writing.\n";
	open(OUT_REM_FH, ">$remainder_fn") || die "Could not open $remainder_fn or writing.\n";

	my @ids_arr=keys %{$hash_ref};
	my $num_ids=$#ids_arr+1;

	print STDERR "Number of sequences to extract: $num_ids\n";

	my $found=0;
	while(!eof(FASTQ_FH)){

		my $defline=<FASTQ_FH>;
		my $seq=<FASTQ_FH>;
		my $plus=<FASTQ_FH>;
		my $qv=<FASTQ_FH>;

		if($defline=~/^@(\S+)/){
			my $id=$1;

			# for SOLEXA reads
			if($id=~/(.+)\/\d$/){
				$id=$1;
			}

			if(defined(${$hash_ref}{$id})){
				print OUT_FH "$defline$seq$plus$qv";
				$found++;
			}else{
				print OUT_REM_FH "$defline$seq$plus$qv";
			}
		}else{
			die "Error: Could not parse ID out of defline: $defline\n";
		}

	}

	close(FASTQ_FH);
	close(OUT_FH);

	if($found != $num_ids){
		die "Error: Number of IDs targeted ($num_ids) not matching those found ($found).\n";
	}
	print STDERR "Extraction done.  $found sequences pulled.\n";

	return;
}

###############################################################################

sub merge_fastq{

	my $fname1=shift;
	my $fname2=shift;
	my $out_fn=shift;

	print STDERR "Merging FASTQ files: $fname1 and $fname2.\n";

	if($fname1=~/\.gz$/){
		open(FASTQ1_FH, "zcat $fname1 | ") || die "Could not open $fname1.\n";
	}else{
		open(FASTQ1_FH, "<$fname1") || die "Could not open $fname1.\n";
	}

	if($fname2=~/\.gz$/){
		open(FASTQ2_FH, "zcat $fname2 | ") || die "Could not open $fname2.\n";
	}else{
		open(FASTQ2_FH, "<$fname2") || die "Could not open $fname2.\n";
	}

	my $num_sequences_merged=0;

	open(OUT_FH, ">$out_fn") || die "Could not open $out_fn\n";

	while(!eof(FASTQ1_FH)){

		my $defline=<FASTQ1_FH>;
		my $seq=<FASTQ1_FH>;
		my $plus=<FASTQ1_FH>;
		my $qv=<FASTQ1_FH>;
		my $id1;

		if($defline=~/^@(\S+)/){
			$id1=$1;	
			if($id1=~/(.+)\/\d$/){
				$id1=$1;
			}
		}

		print OUT_FH "$defline$seq$plus$qv";

		my $defline=<FASTQ2_FH>;
		my $seq=<FASTQ2_FH>;
		my $plus=<FASTQ2_FH>;
		my $qv=<FASTQ2_FH>;
		my $id2;

		if($defline=~/^@(\S+)/){
			$id2=$1;	
			if($id2=~/(.+)\/\d$/){
				$id2=$1;
			}
		}

		print OUT_FH "$defline$seq$plus$qv";

		if($id1 ne $id2){
			print STDERR "Error:  The two FASTQ files $fname1 and $fname2 are not lining up\n";
			print STDERR "at $id1 / $id2.\n";
			print STDERR "Incompatible FASTQ files.  Entries do not line up exactly.\n\n";
			return(-1);
		}

		$num_sequences_merged++;
	}

	close(FASTQ1_FH);
	close(FASTQ2_FH);
	close(OUT_FH);

	print STDERR "Number of sequences merged: $num_sequences_merged\n";

	return(0);
}

###############################################################################

# Identify read IDs that are common between two fastq files
my %seen_id_hash;
get_ids($for_name, \%seen_id_hash);
get_ids($rev_name, \%seen_id_hash);

# Only keep IDs that have been seen in both files
my %target_id_hash;
foreach my $key(keys %seen_id_hash){
	if($seen_id_hash{$key}==2){
		$target_id_hash{$key}=1;
	}
}

my $hostname=hostname();
my $pid=$$;
my $time=time();
my $tmp_for_out="$temp_directory/$output_fname.tmp.$hostname.$pid.$time.for.fastq";
my $tmp_rev_out="$temp_directory/$output_fname.tmp.$hostname.$pid.$time.rev.fastq";

# Extract sequences from both files
extract_ids($for_name, \%target_id_hash, $tmp_for_out, "$output_name.for.frag.fastq");
extract_ids($rev_name, \%target_id_hash, $tmp_rev_out, "$output_name.rev.frag.fastq");

# Merge sequences from both files
	# Confirm that the ID's are the same
my $err=merge_fastq($tmp_for_out, $tmp_rev_out, "$output_name.paired.fastq");

`rm $tmp_for_out`;
`rm $tmp_rev_out`;

###############################################################################

if($err==-1){
	`rm $output_name.paired.fastq`;
	`rm $output_name.for.frag.fastq`;
	`rm $output_name.rev.frag.fastq`;
	print STDERR "ERROR Creating paired FASTQ file!!!\n\n";
}else{
	print STDERR "Completed.\n\n";
}

###############################################################################
