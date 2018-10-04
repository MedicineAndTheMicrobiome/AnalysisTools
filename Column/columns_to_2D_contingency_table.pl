#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_c $opt_r $opt_o);

getopts("i:c:r:o:");

my $usage = "
	usage:
	$0
		-i <input columns file>
		-c <column number in file, for columns in contigency table, from 0>
		-r <column number in file, for rows in contigency table, from 0>	
		-o <output contingency table file name>

	This script will read in a file with at least 2 columns and
	generate a 2D contingency table based on the counts for each pair
	of column combinations.  In other words, if you specify and X
	and Y column (with the -c and -r options), then the script
	will run through all the possible combinations of X and Y
	and count up their occurrences.

	Note that the input file should not have any headers in them.

	For example:
	Input:
		apples	dogs
		apples	dogs
		bananas	monkeys
		bananas	dogs
		oranges	cats
		oranges	dogs


	Output:
		apples	bananas	oranges
	dogs	2	1	1
	cats	0	0	1
	monkeys	0	1	0

	
";

if(!defined($opt_c) || 
	!defined($opt_r) || 
	!defined($opt_o)){
	die $usage;
}

my $column_id=$opt_c;
my $row_id=$opt_r;
my $infile=$opt_i;
my $outfile=$opt_o;

print STDERR "Column Index: $column_id\n";
print STDERR "Row Index: $row_id\n";

###############################################################################

open(IN_FH, "<$infile") || die "Could not open map file $infile\n";

###############################################################################


my %table;
my %column_keys;
my %row_keys;
my $lines=0;

while(<IN_FH>){
	chomp;
	my @cols=split "\t", $_;
	my $ckey=$cols[$column_id];
	my $rkey=$cols[$row_id];

	$column_keys{$ckey}=1;
	$row_keys{$rkey}=1;
	
	if(!defined($table{$ckey}{$rkey})){
		$table{$ckey}{$rkey}=1;
	}else{
		$table{$ckey}{$rkey}++;
	}
	
	$lines++;

}

print STDERR "Lines read: $lines\n";

close(IN_FH);

###############################################################################

my @column_keys=sort keys %column_keys;
my @row_keys=sort keys %row_keys;

print STDERR "Column Keys:\n";
foreach my $key(@column_keys){
	print STDERR "\t$key\n";
}

print STDERR "\n";
print STDERR "Row Keys:\n";
foreach my $key(@row_keys){
	print STDERR "\t$key\n";
}
print STDERR "\n";


###############################################################################


open(OUT_FH, ">$outfile") || die "Could not open output file $outfile\n";

print OUT_FH "\t" . (join "\t", @column_keys) . "\n";

foreach my $rkey (@row_keys){
	print OUT_FH "$rkey";
	foreach my $ckey (@column_keys){
		my $val=$table{$ckey}{$rkey};
		if(!defined($val)){
			$val=0;
		}
		print OUT_FH "\t$val";

	}
	print OUT_FH "\n";
}

close(OUT_FH);


print STDERR "done.\n";

