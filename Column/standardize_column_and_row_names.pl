#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_i $opt_d $opt_t $opt_a);

getopts("i:dta");

my $usage = "
	usage:
	$0
		-i <Input char-delimited file name>
		
		Input delimitor options
		[-d (auto Detected, default)]
		[-t (Tab delimited)]
		[-a (commA delimited)]

	Output goes to STDOUT.

	This script will clean the row and column names.
	It will not modify the data in the cells.

";

if(!defined($opt_i)){
	die $usage;
}

# Required
my $InputFile=$opt_i;

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

print STDERR "Input: $InputFile\n";
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
	$item=~s/,/_/g;
	$item=~s/ /_/g;
	$item=~s/\{/_/g;
	$item=~s/\}/_/g;
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
	$item=~s/\:/\./g;
	$item=~s/\;/\./g;
	$item=~s/\?/\./g;
	$item=~s/\//\./g;
	$item=~s/\-/\./g;

	# Convert characters to string
	$item=~s/\+/p/g;
	$item=~s/\%/pct/g;
	$item=~s/\=/.eq./g;
	$item=~s/\*/.x/g;

	print STDERR "$before / $item\n";
	
	return($item);
}

###############################################################################

print STDERR "Processing Input File...\n";
my $input_lines=0;
open(IN_FH, "<$InputFile") || die "Could not open input file $InputFile\n";

while(<IN_FH>){
	chomp;

	if($_ eq ""){
		next;
	}
	
	my @array=split /$Delimitor/, $_, -1;
	if($input_lines==0){
		my $num_col=$#array+1;
		my @cleaned;
		for(my $i=0; $i<$num_col; $i++){
			$array[$i]=clean_name($array[$i]);
		}
	}else{
		my $num_col=$#array+1;
		$array[0]=clean_name($array[0]);
	}

	my  $outstr=join "$Delimitor", @array;
	print STDOUT "$outstr\n";

	$input_lines++;
}

close(IN_FH);

###############################################################################

print STDERR "Num Lines Processed: $input_lines\n";
print STDERR "Done.\n";
