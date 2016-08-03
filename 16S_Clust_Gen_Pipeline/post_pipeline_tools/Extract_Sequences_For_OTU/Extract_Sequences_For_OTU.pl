#!/usr/bin/env perl

###############################################################################

use strict;
use FindBin ();
use Getopt::Std;
use FileHandle;
use File::Basename;
use vars qw($opt_o $opt_d $opt_f $opt_r $opt_p $opt_g);

my $PIPELINE_UTIL_PATH="$FindBin::Bin/../../pipeline_utilities";
my $EXTRACT_FASTA_PATH="$PIPELINE_UTIL_PATH/Extract_Record_From_FASTA_By_List.pl";



print STDERR "Path of Pipeline Utilities: $PIPELINE_UTIL_PATH\n";

my $OTU_CUTOFF=0.03;

getopts("o:d:f:r:p:g:");
my $usage = "usage: 

$0 

	-o <OTU ID, e.g., Otu00013>
	-d <Output directory>
	-f <FASTA file contain all sequences>
	-r <OTU to (Representative) Read IDs file, i.e., *.an.list file>
	-p <Representative to Represented Read IDs file, i.e. *.names file>
	-g <Read ID to Group ID Mapping file, i.e. *.groups file>

	This script will extract out all the sequences for the
	specified OTU. If you are looking for a specific OTU
	you can probably find it by looking for its taxonomic
	assignment in the *.cons.taxonomy file.

	
	The -p option should be filled with a file like this:
		*.unique.good.filter.unique.precluster.pick.names
	
	This file represents the most collapsed representative-to-
	represented map file.  

	By default, the cutoff is $OTU_CUTOFF for OTU formation.

";

if(!(
	defined($opt_o) && 
	defined($opt_d) && 
	defined($opt_f) && 
	defined($opt_r) && 
	defined($opt_p) && 
	defined($opt_g))){
	die $usage;
}

my $target_otu=$opt_o;
my $output_dir=$opt_d;
my $fasta_file=$opt_f;
my $otu_to_read_file=$opt_r;
my $read_to_rep_file=$opt_p;
my $read_to_grp_file=$opt_g;

print STDERR "\n";
print STDERR "Target OTU: '$target_otu'\n";
print STDERR "Output Dir: $output_dir\n";
print STDERR "FASTA File: $fasta_file\n";
print STDERR "OTU-Read Map (an.list): $otu_to_read_file\n";
print STDERR "Read-Rep Map (names): $read_to_rep_file\n";
print STDERR "Read-Group Map (groups): $read_to_grp_file\n";
print STDERR "\n\n";

print STDERR `mkdir $output_dir/$target_otu`;

if(!(-e "$output_dir/$target_otu")){
	die "Could not make directory: $output_dir/$target_otu\n";
}
my $outdir="$output_dir/$target_otu";

mkdir "$output_dir/$target_otu/ids";
mkdir "$output_dir/$target_otu/fastas";

###############################################################################

my $read_arr_ref=load_anlist_file($otu_to_read_file, $OTU_CUTOFF, $target_otu);
my $represented_reads_arr_ref=load_names_file($read_to_rep_file, $read_arr_ref);
my $read_grp_hash_ref=load_groups_file($read_to_grp_file, $represented_reads_arr_ref);

foreach my $group(sort keys %{$read_grp_hash_ref}){
	write_list(${$read_grp_hash_ref}{$group}, "$outdir/ids/$group\.ids");
	my $cmd_str=
		"$EXTRACT_FASTA_PATH " .
		"-f $fasta_file " .
		"-l $outdir/ids/$group\.ids " .
		"> $outdir/fastas/$group\.fasta";

	system($cmd_str);
}

###############################################################################

sub write_list{
	my $arr_ref=shift;
	my $fname=shift;
	
	print STDERR "Writing list: $fname\n";
	open(FH, ">$fname") || die "Could not open $fname\n";

	foreach my $ids(@{$arr_ref}){
		print FH "$ids\n";
	}

	close(FH);
	return;
}

sub load_anlist_file{
	my $fname=shift;
	my $cutoff=shift;
	my $tar_otu_id=shift;

	# First line is Otu ID list
	# Example:
	# label   numOtus	<OTU ID1>	<OTU ID2> ... <OTU IDn>
	# unique  80810		<read ids,,,>	<read ids,,,> ... <
	# 0.01    41150		<read ids,,,>	<read ids,,,> ...
	# 0.02    17537		<read ids,,,>	<read ids,,,> ...
	# 0.03    8747		<read ids,,,>	<read ids,,,> ...
	#

	# Read first line
	# Read lines until target cutoff found
	# get column for otu ID

	open(FH, "<$fname") || die "Could not open $fname.\n";

	# Extract the header row
	my $header=<FH>;
	my @otu_ids=split "\t", $header;

	# Find column with the correct OTU ID
	my $target_ix=0;
	my $i=0;
	foreach my $otu_id(@otu_ids){
		if($otu_id eq $tar_otu_id){
			$target_ix=$i;
			last;
		}	
		$i++;
	}
	if($target_ix==0){
		die "Error, could not find requested OTU: $tar_otu_id\n";
	}else{
		print STDERR "Targeted OTU found at index: $target_ix\n";
	}

	# Find column with the correct cutoff, then extract that column.
	my $reads_str="";
	while(<FH>){
		chomp;
		my @fields=split "\t", $_;
		if($fields[0] eq $cutoff){
			$reads_str=$fields[$target_ix];
		}
	}
	if($reads_str eq ""){
		die "Could not find requested cutoff: $cutoff\n";
	}

	close(FH);

	my @reads_arr=split ",", $reads_str;
	print STDERR "Num representative reads found: ", ($#reads_arr+1), "\n";

	return(\@reads_arr);

}

#------------------------------------------------------------------------------

sub load_names_file{
	my $fname=shift;
	my $tar_arr_ref=shift;

	my %rep_hash;
	foreach my $id(@{$tar_arr_ref}){
		$rep_hash{$id}=1;
	}

        open(FH, "<$fname") || die "Could not open $fname.\n";
	
	my @all_reads;
	my $cts=0;
	while(<FH>){
		chomp;
		my ($reptv, $reptd)=split "\t", $_;
		if(defined($rep_hash{$reptv})){
			#print STDERR "Rep found: $reptv\n";
			my @reptd_arr=split ",", $reptd;
			foreach my $rep_id(@reptd_arr){
				push @all_reads, $rep_id;
				$cts++;
			}
		}
	}
	
	close(FH);

	print STDERR "Num reads represented: ", ($#all_reads+1), "\n";
	return(\@all_reads);
}

#------------------------------------------------------------------------------

sub load_groups_file{
	my $fname=shift;
	my $tar_reads_ref=shift;
	
	my %tar_hash;
	foreach my $reads(@{$tar_reads_ref}){
		$tar_hash{$reads}=1;
	}
	
	open(FH, "<$fname") || die "Could not open $fname.\n";
		
	my %groups_hash;
	while(<FH>){
		chomp;
		my ($read_id, $group_id)=split "\t", $_;
		if(defined($tar_hash{$read_id})){
			push @{$groups_hash{$group_id}}, $read_id;
		}
	}

	close(FH);
	
	print STDERR "Num Groups: ", (scalar keys %groups_hash), "\n";

	return(\%groups_hash);

}

#------------------------------------------------------------------------------


###############################################################################

print STDERR "done.\n";

