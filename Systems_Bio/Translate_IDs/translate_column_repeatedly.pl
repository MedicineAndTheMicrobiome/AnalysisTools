#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_i $opt_o $opt_m $opt_c);

getopts("i:o:c:m:");

my $usage = "
	usage:
	$0
		-i <Input filename>
		-o <Translated output filename>
		-m <Map file, tab-separated>
		-c <list of Columns to translate, starting from 0>

	Reads in the input file and for each lines splits out columns based on tabs.
	All the IDs that can be translated in that specified column is translated.
	
	This program does the translation brute force so as to avoid parsing and tokenizing.
	It's slow.

	The map file should have the format:
		<key>\\t<new value>\\n

";

if(!defined($opt_i) || !defined($opt_m) || !defined($opt_c)){
	die $usage;
}

my $InputFile=$opt_i;
my $OutputFile=$opt_o;
my $MapFile=$opt_m;
my $ColumnNum=$opt_c;

###############################################################################

print STDERR "Input: $InputFile\n";
print STDERR "Map: $MapFile\n";
print STDERR "Column Number: $ColumnNum\n";

###############################################################################

my %map;
my $num_map_entries=0;
open(MAP_FH, "<$MapFile") || die "Could not open $MapFile\n";
while(<MAP_FH>){
	chomp;
	my ($key, $val)=split /\t/, $_;

	if($val ne ""){
		$map{$key}=$val;
	}else{
		$map{$key}="UNMAPPED";
	}
	$num_map_entries++;
}
close(MAP_FH);

print STDERR "Num Map Entries Read: $num_map_entries\n";

###############################################################################

delete($map{""});
my @keys=keys %map;

sub translate{
	my $in=shift;

	foreach my $key(@keys){
		my $value=$map{$key};
		$in=~s/$key/$value/g;
	}		

	return("$in");
}

###############################################################################

print STDERR "Processing Input File...\n";
my $input_lines=0;

open(IN_FH, "<$InputFile") || die "Could not open input file $InputFile\n";
open(OUT_FH, ">$OutputFile") || die "Could not open input file $OutputFile\n";


my $fh;
my $temp_filename;

while(<IN_FH>){
	chomp;
	my @array=split /\t/, $_;
	my $str_to_translate=$array[$ColumnNum];
	my $translated_str=translate($str_to_translate);
	$array[$ColumnNum]=$translated_str;
	my $joined=join "\t", @array;
	print OUT_FH "$joined\n";

	$input_lines++;
	if(!($input_lines%10)){
		print STDERR ".";
	}
}

close(IN_FH);
close(OUT_FH);

print STDERR "\n";

print STDERR "Num Lines Processed: $input_lines\n";

print STDERR "Done.\n";
