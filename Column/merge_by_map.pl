#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_m $opt_k $opt_h $opt_K $opt_s $opt_o $opt_v);

getopts("i:m:k:h:K:s:o:v");

my $usage = "
	usage:
	$0
		-i <input file name, i.e. recipient file>
		-m <map file name, i.e. donor file>
		-k <key column in the input file, starting from 0>
		-h <files contain header information, Y=yes, N=no>
		-o <output file name>

		Optional:
		[-K <key column in the map file, starting from 0, default=0>]
		[-s <where to insert columns in input file, starting from 0, default=end>]

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
	!defined($opt_m) || 
	!defined($opt_k) || 
	!defined($opt_o) ||
	!defined($opt_h) 
){
	die $usage;
}

my $InputFile=$opt_i;
my $MapFile=$opt_m;
my $InputKeyCol=$opt_k;
my $Header=uc($opt_h);
my $MapKeyCol=$opt_K;
my $InsertCol=$opt_s;
my $OutputFile=$opt_o;
my $Verbose=defined($opt_v);

if(!defined($MapKeyCol)){
	$MapKeyCol=1;
}

if(!defined($InsertCol)){
	$InsertCol=undef;
}

print STDERR "Input File: $InputFile\n";
print STDERR "Map File: $MapFile\n";
print STDERR "Output File: $OutputFile\n";
print STDERR "Contains Header: $Header\n";
print STDERR "\n";
print STDERR "Input Key Column: $InputKeyCol\n";
print STDERR "Map Key Column: $MapKeyCol\n";
print STDERR "Insertion Column: $InsertCol\n";
print STDERR "\n";

###############################################################################

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
	splice @out_arr, $MapKeyCol, 1;
	my $out_str=join "\t", @out_arr;
	my $key=$array[$MapKeyCol];

	# Store header line separately
	if($line==0 && $Header=="Y"){
		$map_header=$out_str;
	}else{
		$map_hash{$key}=$out_str;
	}

	$line++;
}
close(MAP_FH);

###############################################################################

if($Verbose){
	sub print_arr{
		my $arr_ref=shift;
		my $num_out=10;
		for(my $i=0; $i<$num_out; $i++){
			print STDERR "${$arr_ref}[$i]\n";
		}
	}

	# Output Map Keys
	my @map_keys=keys %map_hash;
	print STDERR "Map Keys:\n";
	print_arr(\@map_keys);

	print STDERR "\n";

	# Output File Keys
	my @in_keys;
	open(IN_FH, "head -n 10 $InputFile |") || die "Could not open $InputFile\n";
	while(<IN_FH>){
		chomp;
		my @arr=split /\t/, $_;
		push @in_keys, $arr[$InputKeyCol];
	}
	close(IN_FH);
	
	print STDERR "Input Keys:\n";
	print_arr(\@in_keys);

}

###############################################################################

open(INPUT_FH, "<$InputFile") || die "Could not open input file: $InputFile\n";
open(OUTPUT_FH, ">$OutputFile") || die "Could not open output file: $OutputFile\n";

print STDERR "Iterating through input file...\n";
my $line=0;
while(<INPUT_FH>){
	chomp;
	my @array=split /\t/, $_, -1;
	my @out_array=@array;

	if(!defined($InsertCol)){
		$InsertCol=$#array+1;
	}

	if($line==0 && $Header=="Y"){
		# Output header with new columns inserted in
		splice @out_array, $InsertCol, 0, $map_header;
	}else{

		my $key=$array[$InputKeyCol];
		my $map_val=$map_hash{$key};

		if(!defined($map_val)){
			$map_val="UNDEFINED";			
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
