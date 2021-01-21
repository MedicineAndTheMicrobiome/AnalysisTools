#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_i $opt_o $opt_m $opt_d $opt_t $opt_a);

getopts("i:o:m:dta");

my $usage = "
	usage:
	$0
		-i <Input char-delimited file name>
		-o <Output char-delimited file name>
		-m <Output mapping file>
		
		Input delimitor options
		[-d (auto Detected, default)]
		[-t (Tab delimited)]
		[-a (commA delimited)]

	Output goes to STDOUT.

	This script will clean all cells in the input
	and generate a mapping file.

	To keep the mapping consistent between multiple
	files the same version of this script should
	be applied to all files.

	All characters will be converted into periods, underscores, or alpha-numeric.
	If value starts with a number, then X__ will be prepended to it.

";

if(!defined($opt_i) || !defined($opt_o) || !defined($opt_m)){
	die $usage;
}

# Required
my $InputFile=$opt_i;
my $OutputFile=$opt_o;
my $MappingFile=$opt_m;

# Delimitors
my $Delimitor="";
if(defined($opt_d)){
	$Delimitor="";
}
if(defined($opt_t)){
	$Delimitor="\t";
}
if(defined($opt_a)){
	$Delimitor=",";
}

###############################################################################
# Echo selected options/parameters

print STDERR "Input:   $InputFile\n";
print STDERR "Output:  $OutputFile\n";
print STDERR "Mapping: $MappingFile\n";
print STDERR "Delimitor: '$Delimitor'\n";

###############################################################################
# Autodetect field delimiter

if($Delimitor eq ""){
	my $num_commas=0;
	my $num_tabs=0;

	print STDERR "Attempting to autodetected delimitor.\n";
	open(IN_FH, "head -n 20 $InputFile|") || die "Could not open $InputFile.\n";
	while(<IN_FH>){
		chomp;
		my $test=$_;
		$num_commas+=($test=~tr/,/,/);
		$num_tabs+=($test=~tr/\t/\t/);
	}
	close(IN_FH);

	print STDERR "Num Tabs: $num_tabs\n";
	print STDERR "Num Commas: $num_commas\n";
	if($num_tabs>(1.2*$num_commas)){
		$Delimitor="\t";
	}elsif($num_commas>(1.2*$num_tabs)){
		$Delimitor=",";
	}else{
		print STDERR "Could not detect delimitor.\n";
		print STDERR "Maybe it doesn't matter. Assuming <tab>\n";
		$Delimitor="\t";
		
	}
	print STDERR "I believe the delimitor is '$Delimitor'\n";
	print STDERR "\n";
}

###############################################################################

sub clean_name{		
	my $item=shift;
	my $before=$item;

	# Convert characters to underscores first
	$item=~s/ /_/g;
	$item=~s/\(/_/g;
	$item=~s/\)/_/g;
	$item=~s/\]/_/g;
	$item=~s/\[/_/g;

	# Remove excess underscores
	$item=~s/_+/_/g;
	$item=~s/^_//;
	$item=~s/_$//;

	# Remove characters
	$item=~s/"//g;
	$item=~s/'//g;
	
	# Convert characters to .
	$item=~s/\,/\./g;
	$item=~s/\;/\./g;
	$item=~s/\:/\./g;
	$item=~s/\?/\./g;
	$item=~s/\//\./g;
	$item=~s/\-/\./g;

	# Convert characters to string
	$item=~s/\+/p/g;
	$item=~s/\%/pct/g;
	$item=~s/\=/.eq./g;
	$item=~s/\*/.x/g;

	if($item=~/^\d+/){
		$item="X__$item";
	}

	#print STDERR "$before / $item\n";
	
	return($item);
}

###############################################################################

print STDERR "Processing Input File...\n";
my $input_lines=0;
open(IN_FH, "<$InputFile") || die "Could not open input file $InputFile.\n";
open(OUT_FH, ">$OutputFile") || die "Could not open output file $OutputFile.\n";

my %mapping_hash;
my %inv_hash;

while(<IN_FH>){
	chomp;

	if($_ eq ""){
		next;
	}
	
	my @array=split /$Delimitor/, $_, -1;

	my $num_col=$#array+1;
	for(my $i=0; $i<$num_col; $i++){

		my $original=$array[$i];
		my $cleaned=clean_name($original);

		my $prev_orig=$inv_hash{$cleaned};
		if(!defined($prev_orig)){
			$inv_hash{$cleaned}=$original;
		}else{
			if($original ne $prev_orig){
				print STDERR "\n";
				print STDERR "Error: Cleaned name not unique.\n";
				print STDERR "Cleaned:  $cleaned\n";
				print STDERR "Current:  $original\n";
				print STDERR "Previous: $prev_orig\n";
				print STDERR "\n";
				die "Update cleaning algorithm to fix degenerate cleaning.\n";
			} 
		}
		$mapping_hash{$original}=$cleaned;

		$array[$i]=$cleaned;
	}

	my  $outstr=join "$Delimitor", @array;
	print OUT_FH "$outstr\n";

	$input_lines++;
}
close(IN_FH);


print("Writing mapping file...\n");
open(MAP_FH, ">$MappingFile") || die "Could not open output file $MappingFile.\n";
foreach my $orig (sort keys %mapping_hash){
	print MAP_FH "$orig\t$mapping_hash{$orig}\n";
}
close(MAP_FH);



###############################################################################

print STDERR "Num Lines Processed: $input_lines\n";
print STDERR "Done.\n";
