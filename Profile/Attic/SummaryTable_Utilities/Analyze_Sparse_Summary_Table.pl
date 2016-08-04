#!/usr/bin/env perl

use strict;
use Getopt::Std;
use File::Basename;
use FileHandle;
use vars qw($opt_i $opt_s);
getopts("i:s:");


my $usage = "usage:
$0
	-i <input summary table>

	This script will read in a summary table that is very sparse (i.e.
	mostly 0's) and then for each sample display the the number of 
	singleton categories, doubleton categories, ...

	For example, below, the sample T63R has 360366 items across
	all it's categories.  There are 22,884 singletons, which contributes
	6.35% of that samples abundance.  i.e. 22884*1/360366 *100% = 6.35%.
	For the doubleton, 3500*2/360366 * 100% = 1.94%.

	T63R  Total=360366
	  [0] = 864686 (0.0000%)
	  [1] = 22884 (6.3502%)
	  [2] = 3500 (1.9425%)
	  [3] = 1133 (0.9432%)
	  [4] = 475 (0.5272%)
	  [5] = 274 (0.3802%)
	  [6] = 185 (0.3080%)
	  [7] = 120 (0.2331%)


";

if(!(
	defined($opt_i)
)){
	die $usage;
}


my $input_st=$opt_i;
my @ext=(".summary_table.tsv", ".summary_table.xls");
my ($name, $path, $extension)=fileparse($input_st, @ext);

print STDERR "Input File: $input_st\n";
print STDERR "Output Root: $path/$name\n";

###############################################################################

open(FH, "<$input_st") || die "Could not open $input_st\n";

# Read in header info
my $hdr_line=<FH>;
chomp $hdr_line;

my @hdr_col=split "\t", $hdr_line;

my @categories=@hdr_col;
shift @categories;	# Shift out sample_id
shift @categories;	# Shift out total

my $num_categories=$#categories+1;
print STDERR "Num categories found: $num_categories\n";

# Read in counts, but process 1 line at a time.
while(<FH>){
	chomp;

	my @data_col=split "\t", $_;

	my $sample_id=shift @data_col;
	my $total=shift @data_col;

	my %counts_hash;
	my $sum=0;
	foreach my $count(@data_col){

		if(!defined($counts_hash{$count})){
			$counts_hash{$count}=1;
		}else{
			$counts_hash{$count}++
		}
		$sum+=$count;
	}
	
	print "$sample_id  Total=$total\n";
	for(my $i=0; $i<=20; $i++){
		if(!defined($counts_hash{$i})){
			$counts_hash{$i}=0;
		}
		print "  [$i] = " . $counts_hash{$i} . " (" . sprintf("%4.4f", 100*$counts_hash{$i}*$i/$sum) . "%)\n";
	}

	print "\n\n";

}
