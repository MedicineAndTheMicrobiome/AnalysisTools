#!/usr/bin/env perl

###############################################################################

use strict;
use File::Basename;
use Getopt::Std;
use vars qw ($opt_i $opt_k $opt_I $opt_K $opt_o $opt_s);

getopts("i:k:I:K:o:s");

my $usage = "
	usage:
	$0
		Input File Parameters:
		-i <input file name, i.e. recipient file>
		-k <key column name>
		
		Mapping File Parameters:
		-I <map file name, i.e. donor file>
		-K <key column in the map file>

		Output File Name:
		-o <output file name>

		Options:
		[-s (exclude lines if no mapping, i.e. NAs)]

	This script will merge the input file with the map file.
	The map file is read in in its entirety into a hash
	so that it doesn't need to be sorted.  The column name
	specified should have unique value though.

	This script assumes both files have headers, including
	the first column.

	The order of the input file will be preserved.
	All values from map file (excluding the column key) will
	be appended to end of output file.

	Use the -s option/flag to skip lines

";

if(
	!defined($opt_i) || 
	!defined($opt_k) || 
	!defined($opt_I) || 
	!defined($opt_K) ||
	!defined($opt_o) 
){
	die $usage;
}

my $InputFile=$opt_i;
my $InputColname=$opt_k;

my $MapFile=$opt_I;
my $MapFileColname=$opt_K;

my $OutputFile=$opt_o;

my $SkipNoMap=defined($opt_s);

#------------------------------------------------------------------------------


print STDERR "-----------------------------------------\n";
print STDERR "\n";
print STDERR "Input File: $InputFile\n";
print STDERR "Input File Colname: $InputColname\n";
print STDERR "\n";
print STDERR "Map File: $MapFile\n";
print STDERR "Map File Colname: $MapFileColname\n";
print STDERR "\n";
print STDERR "Skip UnMergeable/UnMapped: $SkipNoMap\n";
print STDERR "Output File Name: $OutputFile\n";
print STDERR "\n";
print STDERR "-----------------------------------------\n";

###############################################################################
# Load the mapping file

my %map_hash;

print STDERR "Loading Map File...\n";
open(MAP_FH, "<$MapFile") || die "Could not open map file: $MapFile\n";

my $line=0;
my $map_header="";
my $num_fields_to_insert;
my $map_key_colnum=-1;

while(<MAP_FH>){
	chomp;
	my @array=split /\t/, $_, -1;
	my @out_arr=@array;
	my $num_col=$#array+1;

	# Read header and find key
	if($line==0){

		my $already_found=0;
		for(my $i=0; $i<$num_col; $i++){
			if($MapFileColname eq $array[$i]){
				print STDERR "$MapFileColname found in column $i (starting from 0).\n";
				$map_key_colnum=$i;

				if($already_found==1){
					print STDERR "Column name found multiple times.\n";
					exit;
				}
				$already_found=1;
			}
		}

		if($map_key_colnum==-1){
			print STDERR "Error: Column name $MapFileColname not found.\n";
			exit;
		}

	}

	# Remove key column
	splice @out_arr, $map_key_colnum, 1;
	$num_fields_to_insert=$#out_arr+1;
	my $out_str=join "\t", @out_arr;
	my $key=$array[$map_key_colnum];

	# Store header line separately
	if($line==0){
		$map_header=$out_str;
	}else{
		$map_hash{$key}=$out_str;
	}

	$line++;
}
close(MAP_FH);

print STDERR "Num columns to insert from map file: $num_fields_to_insert\n";

###############################################################################

print STDERR "-----------------------------------------\n";

my $undef_str="NA" . ("\tNA" x  ($num_fields_to_insert-1));
print STDERR "String used when mapping not possible: '$undef_str'\n";

print STDERR "-----------------------------------------\n";

###############################################################################
# Process the main input file

open(INPUT_FH, "<$InputFile") || die "Could not open input file: $InputFile\n";
open(OUTPUT_FH, ">$OutputFile") || die "Could not open output file: $OutputFile\n";

my %unmapped_hash;
print STDERR "Iterating through input file...\n";
my $line=0;

my $input_key_colnum;
my $num_unmappable_lines=0;
while(<INPUT_FH>){
	chomp;
	my @array=split /\t/, $_, -1;

	my $num_col=$#array+1;
	my @out_array;
	my $skip_line=0;
	if($line==0){

		# Header/Column names

		my $already_found=0;
		for(my $i=0; $i<$num_col; $i++){
			if($InputColname eq $array[$i]){
				print STDERR "$InputColname found in column $i (starting from 0).\n";
				$input_key_colnum=$i;

				if($already_found==1){
					print STDERR "Column name found multiple times.\n";
					exit;
				}
				$already_found=1;
			}
		}

		if($input_key_colnum==-1){
			print STDERR "Error: Column name $InputColname not found.\n";
			exit;
		}

		@out_array=(@array, $map_header);
	}else{
		# Body/Values 
		
		my $key=$array[$input_key_colnum];
		my $map_val=$map_hash{$key};

		if(!defined($map_val)){
			$map_val=$undef_str;
			if(!defined($unmapped_hash{$key})){
				$num_unmappable_lines++;
				print STDERR "No Mapping for: $key\n";
				$skip_line=1;
			}
		}

		@out_array=(@array, $map_val);

	}

	if(!$skip_line){
		print OUTPUT_FH (join "\t", @out_array) . "\n";
	}

	$line++;
}

close(INPUT_FH);
close(OUTPUT_FH);

###############################################################################
print STDERR "\n\n";
print STDERR "Number of unmapping rows in Input/Destination File: $num_unmappable_lines\n";

print STDERR "Done.\n";
