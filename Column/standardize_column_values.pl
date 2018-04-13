#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_i $opt_b $opt_n $opt_s $opt_e $opt_N $opt_c $opt_l $opt_d $opt_t $opt_a);

getopts("i:bn:seNc:l:dta");

my $usage = "
	usage:
	$0
		-i <Input char-delimited file name>
		
		Conversion options:
		[-b (blanks (\"\") to NA flag)]
		[-n \"<comma-separated list of values to convert to NA>\"]
		[-s (convert spaces to underscore)]
		[-e (convert Excel errors to NA)]
		[-N (convert N/A and n/a to NA)]

		Target Columns:
		[-c <list of Columns to translate, starting from 0, default ALL columns>]

		Split fields into list of items:
		[-l <within field/cell list separator, e.g. \";\">]

		Input delimitor options
		[-d (auto Detected, default)]
		[-t (Tab delimited)]
		[-a (commA delimited)]

	Output goes to STDOUT.

	This script will read in a table and update the values so
	they will be more software friendly. In particular, you can
	use it to update a metadata file where you need to convert some values
	to NA if they were coded using an arbitrary value.

	-b will convert blanks to NA
	-n will convert everything in the specified list to an NA
		e.g. -n \"-2,null!\" will convert the strings -2 and null to NA.
	-s will convert spaces (' ') to underscore ('_')

	-e will convert Excel spreadsheet errors to NA:
		#DIV/0, #N/A, #NAME?, #NULL!, #NUM!, #REF!, #VALUE!
	
";

if(!defined($opt_i)){
	die $usage;
}

# Required
my $InputFile=$opt_i;

# Conversion Options
my $ConvBlanks=defined($opt_b);
my @ConvList=split ",", $opt_n;
my $ConvSpaces=defined($opt_s);
my $ConvExcelErr=defined($opt_e);
my $ConvNsA=defined($opt_N);

my @Columns;
if(defined($opt_c)){
	@Columns=split ",", $opt_c;
}

# Within field list separator
my $ListSep="";
if(defined($opt_l)){
	$ListSep=$opt_l;
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

print STDERR "Input: $InputFile\n";
print STDERR "Delimitor: '$Delimitor'\n";

my %map;
foreach my $conv_target(@ConvList){
	$map{$conv_target}=1;
}

if($ListSep ne ""){
	print STDERR "List Separator: $ListSep\n";
}

print STDERR "Columns to Map:\n";

my %target_columns;
if($#Columns==-1){
	print STDERR "  ALL\n";
}else{
	foreach my $col(@Columns){
		print STDERR "  $col\n";
		$target_columns{$col}=1;
	}
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

my %excel_error_hash;
$excel_error_hash{"#DIV/0"}=1;
$excel_error_hash{"#N/A"}=1;
$excel_error_hash{"#NAME?"}=1;
$excel_error_hash{"#NULL!"}=1;
$excel_error_hash{"#NUM!"}=1;
$excel_error_hash{"#REF!"}=1;
$excel_error_hash{"#VALUE!"}=1;

my %changes_hash;

sub clean_item{
	my $item=shift;

	if($ConvSpaces){
		$item=~s/ /_/g;
		$changes_hash{"Spaces Converted"}++;
	}

	if($ConvBlanks){
		if($item eq ""){
			return("NA");
		}
		$changes_hash{"Blanks"}++;
	}

	if(defined($map{$item})){
		$changes_hash{$item}++;
		return("NA");
	}

	if($ConvExcelErr){
		if(defined($excel_error_hash{$item})){
			$changes_hash{$item}++;
			return("NA");
		}
	}

	if($ConvNsA && (uc($item) eq "N/A")){
		$changes_hash{$item}++;
		return("NA");	
	}
	
	return($item);
}

###############################################################################

print STDERR "Processing Input File...\n";
my $input_lines=0;
open(IN_FH, "<$InputFile") || die "Could not open input file $InputFile\n";

while(<IN_FH>){
	chomp;
	my @array=split /$Delimitor/, $_, -1;

	my $num_col=$#array+1;

	for(my $i=0; $i<$num_col; $i++){

		if($#Columns!=-1){
			if(!defined($target_columns{$i})){
				next;
			}
		}

		my $field=$array[$i];

		if($ListSep ne ""){

			my @list=split /$ListSep/, $field, -1;
			my $num_items=$#list+1;
			for(my $j=0; $j<$num_items; $j++){
				$list[$j]=clean_item($list[$j]);
			}
			$field=join "$ListSep", @list;
		}else{
			$field=clean_item($field);
		}

		$array[$i]=$field;
	}

	my $outstr=join "$Delimitor", @array;

	print STDOUT "$outstr\n";
	$input_lines++;
}

close(IN_FH);

###############################################################################

print STDERR "Changes made:\n";
foreach my $items (sort keys %changes_hash){
	print STDERR "\t$items: $changes_hash{$items}\n";
}
print STDERR "\n";

print STDERR "Num Lines Processed: $input_lines\n";
print STDERR "Done.\n";
