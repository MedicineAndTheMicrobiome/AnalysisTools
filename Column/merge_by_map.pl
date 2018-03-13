#!/usr/bin/env perl

###############################################################################

use strict;
use File::Basename;
use Getopt::Std;
use vars qw ($opt_i $opt_k $opt_h $opt_m $opt_K $opt_H $opt_s $opt_n $opt_o $opt_v);

getopts("i:k:h:m:K:H:s:n:o:v");

my $usage = "
	usage:
	$0
		Input File Parameters:
		-i <input file name, i.e. recipient file>
		-k <key column in the input file, starting from 0>
		-h <input files contain header information, Y=yes, N=no>
		
		Mapping File Parameters:
		-m <map file name, i.e. donor file>
		[-K <key column in the map file, starting from 0, default=0>]
		-H <map file contains header information, Y=yes, N=no>

		Merged Output File Parameters:
		[-s <where to insert columns in input file, starting from 0, default=end>]
		[-n <insert column name if there is header in input file, 
			default: use header name in mapping file or filename]
		-o <output file name>

		Optional:
		[-v (show keys flag)]

	This script will merge the input file with the map file.
	The map file is read in in its entirety into a hash
	so that it doesn't need to be sorted.

	The map file contents (with the exclusion of the key column), will
	be inserted into the input file.  By default, the first column
	of the map file is used as the key and the contents are appended
	to the end of the input file.

	It's important to know if your input files will contain column headers.
	

";

if(
	!defined($opt_i) || 
	!defined($opt_k) || 
	!defined($opt_h) || 
	!defined($opt_m) ||
	!defined($opt_o) 
){
	die $usage;
}

my $InputFile=$opt_i;
my $InputFileKeyCol=$opt_k;
my $InputFileHeader=uc($opt_h);

my $MapFile=$opt_m;
my $MapFileKeyCol=$opt_K;
my $MapFileHeader=uc($opt_H);

my $InsertCol=$opt_s;
my $OuputHeaderName=$opt_n;
my $OutputFile=$opt_o;

my $Verbose=defined($opt_v);

#------------------------------------------------------------------------------

if(!defined($MapFileKeyCol)){
	$MapFileKeyCol=0;
}

if(!defined($InsertCol)){
	$InsertCol=undef;
}

print STDERR "Input File: $InputFile\n";
print STDERR "Input File Key Column: $InputFileKeyCol\n";
print STDERR "Input File Contains Header: $InputFileHeader\n";
print STDERR "\n";
print STDERR "Map File: $MapFile\n";
print STDERR "Map File Key Column: $MapFileKeyCol\n";
print STDERR "Map File Contains Header: $MapFileHeader\n";
print STDERR "\n";
print STDERR "Output Insertion Column: $InsertCol\n";
print STDERR "Output Column Name: $OuputHeaderName\n";
print STDERR "Output File Name: $OutputFile\n";
print STDERR "\n";
print STDERR "Verbose? $Verbose\n";
print STDERR "\n";

###############################################################################
# Load the mapping file

my %map_hash;

print STDERR "Loading Map File...\n";
open(MAP_FH, "<$MapFile") || die "Could not open map file: $MapFile\n";

my $line=0;
my $map_header="";

while(<MAP_FH>){
	chomp;
	my @array=split /\t/, $_, -1;
	my @out_arr=@array;

	# Remove key column
	splice @out_arr, $MapFileKeyCol, 1;
	my $out_str=join "\t", @out_arr;
	my $key=$array[$MapFileKeyCol];

	# Store header line separately
	if($line==0 && $MapFileHeader eq "Y"){
		$map_header=$key;
	}else{
		$map_hash{$key}=$out_str;
	}

	$line++;
}
close(MAP_FH);

if($map_header eq ""){
	$map_header=fileparse($MapFile);
}
if(defined($OuputHeaderName)){
	$map_header=$OuputHeaderName;
}

print STDERR "Output Column Name: $map_header\n";


###############################################################################

if($Verbose){
	sub print_arr{
		my $arr_ref=shift;
		my $num_out=10;
		for(my $i=0; $i<$num_out; $i++){
			print STDERR "${$arr_ref}[$i]\n";
		}
		if(($#{$arr_ref}+1)>$num_out){
			print STDERR "...\n";
		}
	}

	# Output Map Keys
	my @map_keys=keys %map_hash;
	print STDERR "Map Keys:\n";
	print_arr(\@map_keys);
	print STDERR "\n";

	# Output File Keys
	my @in_keys;
	open(IN_FH, "head -n 11 $InputFile |") || die "Could not open $InputFile\n";
	while(<IN_FH>){
		chomp;
		my @arr=split /\t/, $_;
		push @in_keys, $arr[$InputFileKeyCol];
	}
	close(IN_FH);
	
	print STDERR "Input Keys:\n";
	print_arr(\@in_keys);
	print STDERR "\n";

}

###############################################################################

open(INPUT_FH, "<$InputFile") || die "Could not open input file: $InputFile\n";
open(OUTPUT_FH, ">$OutputFile") || die "Could not open output file: $OutputFile\n";

my %unmapped_hash;
print STDERR "Iterating through input file...\n";
my $line=0;
while(<INPUT_FH>){
	chomp;
	my @array=split /\t/, $_, -1;
	my @out_array=@array;

	# By default, insert column at end
	if(!defined($InsertCol)){
		$InsertCol=$#array+1;
	}

	if($line==0 && $InputFileHeader eq "Y"){
		# Output header with new columns inserted in
		print STDERR "Inserting $map_header into column $InsertCol...\n";
		splice @out_array, $InsertCol, 0, $map_header;
	}else{

		my $key=$array[$InputFileKeyCol];
		my $map_val=$map_hash{$key};

		if(!defined($map_val)){
			$map_val="UNDEFINED";			
			if($Verbose && !defined($unmapped_hash{$key})){
				$unmapped_hash{$key}=1;
				print STDERR "No Mapping for: $key\n";
			}
		}

		splice @out_array, $InsertCol, 0, $map_val;

	}
	print OUTPUT_FH (join "\t", @out_array) . "\n";

	$line++;
}

close(INPUT_FH);
close(OUTPUT_FH);

###############################################################################

print STDERR "Done.\n";
