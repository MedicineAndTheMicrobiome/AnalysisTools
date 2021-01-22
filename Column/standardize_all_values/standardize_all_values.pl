#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_i $opt_o $opt_m $opt_h $opt_c $opt_d $opt_t $opt_a);

getopts("i:o:m:hc:dta");

my $usage = "
	usage:
	$0
		-i <Input char-delimited file name>
		-o <Output char-delimited file name>
		-m <Output mapping file>
		
		[-h (only work on header, i.e. first line)]
		[-c <column numbers, comma-separated, starting from 1>]

		Input delimitor options
		[-d (auto Detected, default)]
		[-t (Tab delimited)]
		[-a (commA delimited)]

	By default, this script will clean all cells in the input
	and generate a mapping file.

	If you only want the first line cleaned, use the -h flag.
	If you only want specific columns cleaned, use the -c option.
		e.g. -c 1,3,16
	If you specify both the -h and -c option, then the only
		the specified columsn will be cleaned in the header.

	To keep the mapping consistent between multiple
	files the same version of this script should
	be applied to all files.

	All characters will be converted into periods, underscores, or alpha-numeric.
	If a value starts with a number, then X__ will be prepended to it, so
	it shouldn't be applied to data that is numeric.

";

if(!defined($opt_i) || !defined($opt_o) || !defined($opt_m)){
	die $usage;
}

# Required
my $InputFile=$opt_i;
my $OutputFile=$opt_o;
my $MappingFile=$opt_m;

# Optional
my $HeaderOnly=defined($opt_h);
my $ColumnOnly=defined($opt_c);

my @target_columns=();
if($ColumnOnly){
	@target_columns=split /,/, $opt_c;
}

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
print STDERR "\n";
print STDERR "Delimitor: '$Delimitor'\n";
print STDERR "\n";

if($HeaderOnly){
	print STDERR "Only work on header line.\n";
}else{
	print STDERR "Working on all rows.\n";
}

if($ColumnOnly){
	print STDERR "Only working on columns: ", (join ", ", @target_columns), "\n";
}else{
	print STDERR "Working on all columns.\n";
}	

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
open(IN_FH, "<$InputFile") || die "Could not open input file $InputFile.\n";
open(OUT_FH, ">$OutputFile") || die "Could not open output file $OutputFile.\n";

my %mapping_hash;
my %inv_hash;

my %column_hash;
my $all_col=0;
if($ColumnOnly){
	foreach my $target_col(@target_columns){
		$column_hash{$target_col}=1;
	}
}else{
	$all_col=1;
}

my $input_lines=0;
my $process_row;
my $process_col;
my $num_cells_processed=0;

while(<IN_FH>){
	chomp;

	$process_row=0;

	# Conditions of when to process row
	if($HeaderOnly && $input_lines==0){
		$process_row=1;
	}
	if(!$HeaderOnly){
		$process_row=1;
	}

	my @array=split /$Delimitor/, $_, -1;
	my $num_col=$#array+1;
	for(my $i=0; $i<$num_col; $i++){

		$process_col=0;

		# Conditions of when to process col	
		if($all_col==1){
			$process_col=1;
		}else{
			if($column_hash{$i+1}){
				$process_col=1;
			}
		}

		if($process_col && $process_row){
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

			$num_cells_processed++;
		}
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
print STDERR "Num Cells Cleaned: $num_cells_processed\n";
print STDERR "Done.\n";
