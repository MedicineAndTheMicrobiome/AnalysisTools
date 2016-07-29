#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_c $opt_m $opt_e);

getopts("c:m:e:");

my $usage = "
	usage:
	$0
		-c <column list, comma separated, starting from 0, that you want to rename>
		-m <map file>
		-e <extension to add as output file name>
		<file1> <file2> ... <filen>

	Reads map file into hash, then for each line in the input file the column
	specified will be translated to the value in the map file.

	Map file should have the format:
		<new value>\\t<old_value1>,<old_value2>,...,<old_valuen>\\n

	Since the map file maybe very large, you can specify a single input file or a
	list of input files.

";

if(!defined($opt_c) || !defined($opt_m) || !defined($opt_e)){
	die $usage;
}

my @columns=split /,/, $opt_c;
my $map_file=$opt_m;
my $extension=$opt_e;

###############################################################################

print STDERR "Loading mapping file.\n";
my $num_new_ids=0;
my %map;
open(MAP_FH, "<$map_file") || die "Could not open $map_file\n";
while(<MAP_FH>){
	chomp;
	my ($new_val, $old_vals_str)=split /\t/, $_;
	my @old_vals=split /,/, $old_vals_str;
	foreach my $old_val(@old_vals){
		$map{$old_val}=$new_val;
	}
	$num_new_ids++;
}
close(MAP_FH);
print STDERR "Num new/destination IDs loaded: $num_new_ids\n";

###############################################################################

my $file;
while($file=shift){
	print STDERR "Working on $file\n";
	open(IN_FH, "<$file") || die "Could not open map file $file\n";
	open(OUT_FH, ">$file\.$extension") || die "Could not open map file $file\.$extension\n";

	while(<IN_FH>){
		chomp;
		my @array=split /\t/, $_;
		foreach my $col(@columns){
			if(defined($map{$array[$col]})){
				$array[$col]=$map{$array[$col]};
			}
		}
		my $outstr=join "\t", @array;
		print OUT_FH "$outstr\n";
	}

	close(IN_FH);
	close(OUT_FH);
}
