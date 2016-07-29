#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_d $opt_e $opt_p);

getopts("d:ep:");

my $usage = "
	usage:
	$0
		[-d <delimitor, default is tab>]
		[-e (flag to exclude files names as first row)]
		[-p <place holder text>]
		<file1> <file2> ... <filen>

	This script will read in multiple text files, and output them in different columns.
	Useful for merging data to read into excel.

";

my $delim="\t";
if(defined($opt_d)){
	$delim=$opt_d;
}

my $filename_headers=$opt_e;
my $place_holder=$opt_p;
my $print_header=!defined($opt_e);

###############################################################################

my $file;
my @buffer;
my $column=0;
my $row=0;
my $max_row=0;
my $max_col=0;
my @filenames;
while($file=shift){

	open(FH, "<$file") || die "Could not open $file\n";
	push @filenames, $file;

	$row=0;
	while(<FH>){
		chomp;
		$buffer[$row][$column]=$_;
		$row++;
		if($max_row<$row){
			$max_row=$row;
		}
	}
	$column++;
	if($max_col<$column){
		$max_col=$column;
	}

	close(FH);
}

if($column==0){
	print $usage;
}

if($print_header){
	my $header_str=join $delim, @filenames;
	print "$header_str\n";
}

for(my $i=0; $i<$max_row; $i++){
	for(my $j=0; $j<$max_col; $j++){

		if($buffer[$i][$j] eq ""){
			print $place_holder;
		}else{
			print $buffer[$i][$j];
		}
	
		if($j+1 != $max_col){
			print "$delim";
		}
	}
	print "\n";
}
