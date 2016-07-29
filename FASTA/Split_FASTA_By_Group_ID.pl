#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_i $opt_r $opt_x $opt_d);

getopts("i:r:x:d:");
my $usage = "usage: 
$0 
	-i <Input FASTA Filename>
	[-r <Output Root FASTA Filename, eg. clust_>]
	[-x <Output FASTA File extension, eg. .fasta>]
	-d \"<delimitor>\"

	Read in multi-fasta input file and produces a multi-fasta file for each group,
	where a group's id is defined by the string between the > of the
	defline and the delimitor.

	So for example, for the defline:

		>36|*|GQ148864 A/chicken/Israel/182/2008 2008/06/10 2 (PB1)

	The group id is 36 if the delimitor is \"\\|\".

	If the output root FASTA file name is clust_, and the extention is .fasta,
	then the output file that the above sequence will be placed in will be:

		clust_36.fasta

";

if(!(
	defined($opt_i) && 
	defined($opt_d))){
	die $usage;
}

my $in_fasta=$opt_i;
my $output_root=$opt_r;
my $output_extension=$opt_x;
my $delimitor=$opt_d;

print STDERR "Input FASTA: $in_fasta\n";
print STDERR "Delimitor: $delimitor\n";
print STDERR "Output root: $output_root\n";
print STDERR "Output extension: $output_extension\n";

###############################################################################
# Determine groups

open(FASTA_FH, "<$in_fasta") || die "Could not open $in_fasta\n";

my %group_identifiers;
while(<FASTA_FH>){
	chomp;
	if($_=~/^>/){
		my @parts=split /$delimitor/, $_;
		my $id=$parts[0];
		$id=~s/^>//;
		$group_identifiers{$id}=1;
	}
}

close(FASTA_FH);

print STDERR "Identified Groups:\n";
foreach my $id(sort keys %group_identifiers){
	print STDERR "\t$id\n";
}

my @groups=sort keys %group_identifiers;
my $num_groups=$#groups+1;

print STDERR "Num groups: $num_groups\n";

###############################################################################
# Open output files

my %fhandle_hash;
foreach my $group_id(@groups){
	my $outname=$output_root . $group_id . $output_extension;
	$fhandle_hash{$group_id}=new FileHandle ">$outname";
}

###############################################################################
# Read in fasta file in order to write to individual files

print STDERR "Processing FASTA file...\n";
open(FASTA_FH, "<$in_fasta") || die "Could not open $in_fasta\n";
my ($defline, $prev_defline);
my @sequence;
while(<FASTA_FH>){
	chomp;
	
	if(/^>/){
		$defline=$_;
		if($#sequence!=-1){
			process_record($prev_defline, \@sequence);
			@sequence=();
		}
		$prev_defline=$defline;
	}else{
		push @sequence, (split //, $_);
	}
}
process_record($prev_defline, \@sequence);

close(FASTA_FH);

print STDERR "Completed.\n";

###############################################################################

sub process_record{
	my $defline = shift;
	my $sequence_ref = shift;
	
	# Get group id
	my @parts=split /$delimitor/, $defline;
	my $group_id=$parts[0];
	$group_id=~s/^>//;
	
	# Output defline
	#print STDERR "$group_id : $defline\n";
	print {$fhandle_hash{$group_id}} "$defline\n";

	# Output sequence
	my $length=$#sequence+1;
	my $width=80;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		my @out_arr=splice @{$sequence_ref}, 0, $out_width;
		print {$fhandle_hash{$group_id}} (join "", @out_arr) . "\n";

		$length-=$width;
	}while($length>0);

}

###############################################################################

