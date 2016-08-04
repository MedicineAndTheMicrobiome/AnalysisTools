#!/bin/env perl

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_m $opt_b $opt_a);

getopts("i:c:m:ba");

my $usage = "
	usage:
	$0
		-i <input file>
		-m <map file>
		[-b bracket annotated complete replace]
		[-a just annotate]

	Reads map file into hash, then for each line in the input file if a token/word/string matches
	one of the strings in the mapping file to be translated, the string will be translated.

	By default, the replacement will just be with whatever it is mapped to.
	With the -b option, square brackets [<replacement>] will be placed around the replacement text.
	With the -a option, the original will be left there with a :<replacement> after it.

	The map file has the format of:

		<search string>\\t<replacement string>\\n
		<search string>\\t<replacement string>\\n
		<search string>\\t<replacement string>\\n
		<search string>\\t<replacement string>\\n


";

if(!defined($opt_i) || !defined($opt_m)){
	die $usage;
}

my $file=$opt_i;
my $map_file=$opt_m;

my $bracket_replacement=0;
my $annotate=0;

if(defined($opt_b)){
	$bracket_replacement=1;
}

if(defined($opt_a)){
	$annotate=1;
}	

###############################################################################

print STDERR "Loading mapping file.\n";
my %map;
open(MAP_FH, "<$map_file") || die "Could not open $map_file\n";
while(<MAP_FH>){
	chomp;
	my ($key, $val)=split /\t/, $_;
	$map{$key}=$val;
}
close(MAP_FH);
print STDERR "Done.\n";

###############################################################################

open(IN_FH, "<$file") || die "Could not open map file $file\n";

while(<IN_FH>){
	chomp;

	# First split by tabs
	my @tab_split_arr = split /\t/, $_;
	for(my $ts_idx=0; $ts_idx<=$#tab_split_arr; $ts_idx++){

		# Then split by spaces
		my @space_split_arr = split / /, $tab_split_arr[$ts_idx];
		for(my $ss_idx=0; $ss_idx<=$#space_split_arr; $ss_idx++){

			# Translate the token by hash lookup
			my $in=$space_split_arr[$ss_idx];
			my $remap=$map{$in};

			# Format replacement
			if(defined($remap)){
				my $out=$remap;
				if($bracket_replacement){
					$out="[$out]";
				}
				if($annotate){
					$out="$in:$out";
				}
				$space_split_arr[$ss_idx]=$out;
			}
		}

		# Rebuild tab split string by joining by spaces
		$tab_split_arr[$ts_idx]=join " ", @space_split_arr;
	}

	# Rebuild entire line by joining by tabs
	my $outline=join "\t", @tab_split_arr;

	# Reattach end line and output
	print "$outline\n";
}

close(IN_FH);


