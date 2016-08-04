#!/usr/local/bin/perl

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_i $opt_k);

getopts("i:k:");

my $usage = "
	usage:
	$0
		-i <Input Cluster Reps File>
		-k <Keep List>	

	Reads in cluster reps file from running
	Collapse_Categories_By_Correlation and then
	only keeps cluster information based on 
	categories in Keep List.

	Output goes to STDOUT.
";

if(!defined($opt_i) || !defined($opt_k)){
	die $usage;
}

my $InputFile=$opt_i;
my $KeepList=$opt_k;

###############################################################################

print STDERR "\n";
print STDERR "Input:    $InputFile\n";
print STDERR "KeepList: $KeepList\n";
print STDERR "\n";

###############################################################################

print STDERR "Loading map file.\n";
open(LIST_FH, "<$KeepList") || die "Could not open $KeepList\n";

my $list_length=0;
my %keep_list;

while(<LIST_FH>){
	chomp;
	my @fields=split /\t/, $_;
	my $key=shift @fields;
	$keep_list{$key}=0;
	$list_length++;
}

close(LIST_FH);

print STDERR "Num Keep List Records: $list_length\n";

###############################################################################

print STDERR "\nProcessing cluster reps...\n";

open(IN_FH, "<$InputFile") || die "Could not open input file $InputFile\n";

my $hdr=<IN_FH>;
print $hdr;

my $num_recs_kept=0;
while(<IN_FH>){
	chomp;

	if($_=~/^Representative:\t(.+)/){
		my $id=$1;
		my $keep_rec;
		if(defined($keep_list{$id})){
			print "$_\n";
			$keep_rec=1;			
			$keep_list{$id}++;
			$num_recs_kept++;
		}else{
			$keep_rec=0;
		}

		while(<IN_FH>){
			chomp;
			if($_ eq ""){
				if($keep_rec){
					print "\n";
				}
				last;
			}else{
				if($keep_rec){
					print "$_\n";
				}else{
					# NOOP
				}
			}
		}

	}else{
		print STDERR "Error: Could not find Representative tag.\n";
		die;
	}

}

close(IN_FH);

print STDERR "Num Recs found and kept: $num_recs_kept\n";

print STDERR "\n";
print STDERR "Missing records:\n";
my $num_missing=0;
foreach my $id(keys %keep_list){
	if($keep_list{$id}==0){
		print STDERR "$id\n";
		$num_missing++;
	}
}
if($num_missing==0){
	print STDERR "None missing.\n";
}
print STDERR "\n";

print STDERR "Done.\n";
