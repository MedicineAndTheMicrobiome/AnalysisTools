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
		<file1> <file2> ... <filen>

	This script will read in multiple text files, and output them in different columns.
	Useful for merging data to read into excel.

	Even if the number of columns per file are not consistent, the
	script will try to pad out the delimitors so everything still lines up.

";

my $delim="\t";
if(defined($opt_d)){
	$delim=$opt_d;
}

my $filename_headers=$opt_e;
my $place_holder=$opt_p;
my $print_header=!defined($opt_e);

###############################################################################

my %buffer; # Contains text
my %delims; # Contins num delims per line

my @filenames;
my $max_lines_across_files=0;
my %max_delim_per_file;
my %lines_per_file;

# First read in all files into memory
while(my $file=shift){

	push @filenames, $file;

	open(FH, "<$file") || die "Could not open $file\n";

	my $num_lines=0;
	my $max_delim=0;
	while(<FH>){
		$_=~s/\n//;
		my $line=$_;
		push @{$buffer{$file}}, $line;
		
		# Count max delim per file
		my @tmp=split "$delim", $line, -1;
		my $num_delim=$#tmp;
		push @{$delims{$file}}, $num_delim;

		if($num_delim>$max_delim){
			$max_delim=$num_delim
		}	
	
		$num_lines++;
	}
	close(FH);

	print STDERR "Max delim/line: $max_delim  Num lines: $num_lines ($file)\n";

	$max_delim_per_file{$file}=$max_delim;
	$lines_per_file{$file}=$num_lines;

	if($num_lines>$max_lines_across_files){
		$max_lines_across_files=$num_lines;
	}
	
}

print STDERR "Max lines across all files: $max_lines_across_files\n";

my @combined_buffer;
my $num_files=$#filenames + 1;


# Padding Loops
for(my $col=0; $col<$num_files; $col++){
		
	my $file=$filenames[$col];

	# Pad lines with text with extra delim if necessary
	for(my $row=0; $row<$lines_per_file{$file}; $row++){
		$combined_buffer[$row][$col]=${$buffer{$file}}[$row] .
			$delim x ($max_delim_per_file{$file}-${$delims{$file}}[$row]);
	}
	
	# Fill out remaining lines with delim
	for(my $row=($lines_per_file{$file}); $row<$max_lines_across_files; $row++){
		$combined_buffer[$row][$col]=$delim x ($max_delim_per_file{$file})
	}
}

# Output Loops
for(my $row=0; $row<$max_lines_across_files; $row++){
	my @outline;
	for(my $col=0; $col<$num_files; $col++){
		push @outline, $combined_buffer[$row][$col]
	}
	print STDOUT (join $delim, @outline) . "\n";

}

print STDERR "Done.\n";

