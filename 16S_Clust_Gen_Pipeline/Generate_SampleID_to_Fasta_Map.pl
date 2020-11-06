#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use File::Basename;
use vars qw($opt_p $opt_o);

getopts("p:o:");
my $usage = "usage: 

$0 

	-p <path to fasta files>
	-o <output sample id to fasta mapping file>

	This script will look through all the fasta (.fasta, .fa) files
	in the specified path, and generate a sample id to fasta file
	path.

	This is an optional script, but if you have a directory of
	fasta files and you think the fasta file names are consistently
	named and reflect the name of the sample the fasta file it
	represents, then this script will automatically generate
	a sample id to fasta file name.

	The output file is:

	<generated sample id> \\t <fasta path> \\n

";

if(!(
	defined($opt_p) && 
	defined($opt_o))){
	die $usage;
}

my $target_path=$opt_p;
my $output_fname=$opt_o;

print STDERR "\n";
print STDERR "Target Path: $target_path\n";
print STDERR "Output Filename: $output_fname\n";

###############################################################################

my @filelist=split "\n", `find $target_path`;

my @fastalist;

foreach my $fname(@filelist){
	if($fname=~/\.fasta$/ || $fname=~/\.fa$/){
		push @fastalist, $fname;
	}
}

print STDERR "Found FASTA files: \n";
my %map;
foreach my $fpath(@fastalist){
	print STDERR "$fpath\n";
	my ($name, $path)=fileparse($fpath);
	@{$map{$fpath}}=split /\./, $name;
}

# Remove extensions that are common across all files
my $done=0;
while(!$done){
	my %uniq_hash;
	foreach my $fpath(@fastalist){
		my $ext=pop @{$map{$fpath}};
		$uniq_hash{$ext}=1;
		push @{$map{$fpath}}, $ext;
	}
	my $num_keys=keys %uniq_hash;
	if($num_keys==1){
		foreach my $fpath(@fastalist){
			pop @{$map{$fpath}};
		}
	}else{
		$done=1;
	}
}

# If any sample id's are redundant, try to make it unique
print STDERR "Checking Sample IDs for uniqueness...\n";
my %uniq_hash;
my %cnts_hash;

# Count duplicates
foreach my $fpath(keys %map){
	my $samp_id = join ".", @{$map{$fpath}};
	if(defined($uniq_hash{$samp_id})){
		$uniq_hash{$samp_id}++;
		$cnts_hash{$samp_id}++;
		print STDERR "Duplicated Sample ID found: $samp_id\n";
	}else{
		$uniq_hash{$samp_id}=1;
		$cnts_hash{$samp_id}=1;
	}
}

my %sampid_to_path_hash;
my %samp_to_uniqsamp_hash;
# Append ID with r#
foreach my $fpath(keys %map){
	my $samp_id = join ".", @{$map{$fpath}};
	my $uniq_samp_id=$samp_id;
	if($cnts_hash{$samp_id}>1){
		$uniq_samp_id="$samp_id.r$uniq_hash{$samp_id}";	
		$uniq_hash{$samp_id}--;
	}
	$sampid_to_path_hash{$uniq_samp_id}=$fpath;
	$samp_to_uniqsamp_hash{$uniq_samp_id}=$samp_id;
}

###############################################################################

open(OUT_FH, ">$output_fname") || die "Could not open $output_fname\n";

foreach my $samp_id(sort keys %sampid_to_path_hash){
	print OUT_FH "$samp_id\t$sampid_to_path_hash{$samp_id}\n";
}

close(OUT_FH);

#-----------------------------------------------------------------------------

my $collapse_rep_tsv="$output_fname.clps.tsv";
open(OUT_FH, ">$collapse_rep_tsv") || die "Could not open $collapse_rep_tsv\n";

print OUT_FH "ReplicateID\tSampleID\n";
foreach my $uniq_samp_id(sort keys %samp_to_uniqsamp_hash){
	print OUT_FH "$uniq_samp_id\t$samp_to_uniqsamp_hash{$uniq_samp_id}\n";
}

close(OUT_FH);

###############################################################################

print STDERR "done.\n";

