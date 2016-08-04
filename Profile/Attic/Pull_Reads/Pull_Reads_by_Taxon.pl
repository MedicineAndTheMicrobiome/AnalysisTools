#!/usr/bin/env perl

use strict;
use Getopt::Std;
use vars qw($opt_f $opt_t $opt_o);
use FindBin;
getopts("f:t:o:");

my $usage = "usage:
$0
	-f <assignment_detail.txt_tbl.merged file>
	-t \"<taxon string>\"
	[-o <output name>]

	This script will pull out all reads specified by the <taxon string>.

	The taxon string should be in double quotes, so the command line will not 
	break it apart.  If the taxon of interest has spaces in it, please replace
	it with underscores.  Do not use double quotes in the taxon itself.  Spaces
	should separate each taxon level.  The taxon name will not be case sensitive.
	Always start with the domain.  

	The output name is optional.  By default it will concatenate your taxon
	with periods between levels, and produce a file with the extension \".id\".

	For example:
		\"Bacteria Actinobacteria Actinobacteria\"

	Will produce the output file:
		bacteria.actinobacteria.actinobacteria.ids

";

if(!(
	defined($opt_f) &&
	defined($opt_t)
)){
	die $usage;
}


my $input_file=$opt_f;
my $taxon_string=$opt_t;
my $output_name=$opt_o;

########################################################################################

print STDERR "Going to grab out reads for \"$taxon_string\"\n";
$taxon_string=lc($taxon_string);
my @taxon_arr=split / /, $taxon_string;
foreach my $taxon(@taxon_arr){
	print STDERR "  $taxon\n";	
}
print STDERR "\n";

my $taxon_arr_len=$#taxon_arr+1;

########################################################################################

open(IN_FH, "<$input_file") || die "Could not open $input_file\n for reading.";

my @matching_ids;
while(<IN_FH>){
	chomp;
	$_=~s/"//g;

	my @fields=split /\t/, $_;
	my $id=shift @fields;
	
	my $matches=0;
	for(my $i=0; $i<$taxon_arr_len; $i++){
		if(lc($fields[$i]) eq $taxon_arr[$i]){
			$matches++;
		}	
	}

	if($matches==$taxon_arr_len){
		push @matching_ids, $id;
	}

}

close(IN_FH);

########################################################################################

if(!defined($output_name)){
	$output_name=join ".", @taxon_arr;
	$output_name.=".ids";
}
open(OUT_FH, ">$output_name") || die "Could not open $output_name for writing.\n";
foreach my $id(@matching_ids){
	print OUT_FH "$id\n";
}

close(OUT_FH);
