#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_c $opt_k);

getopts("i:c:k:");

my $usage = "
	usage:
	$0
		-i <input table file>
		-c <key column, first column is 0>
		[-k <columns to keep, default is all>]

	Reads in a tab-separated table and extracts a subset of rows and columns
	based on the unique row values in the specified column.  The key column
	will then be moved to the first column in the output.  Header names
	will be preserved.  If collapsing is required because the key column
	has redundant IDs, then the keep columns are tested to be sure that
	they are redunant.  If the data associated with each row can not be
	collapsed into a single row, and error will thrown.

	Output goes to STDOUT.
	
";

if(!defined($opt_i) || !defined($opt_c)){
	die $usage;
}

my $file=$opt_i;
my $key_column=$opt_c;
my @keep_col;

my $keep_all=0;
if(!defined($opt_k)){
	$keep_all=1;
}else{
	@keep_col=split ",", $opt_k;
}

###############################################################################

print STDERR "Key column: $key_column\n";

if(!$keep_all){
	print STDERR "Keep columns: ", (join ";", @keep_col), "\n";
}

###############################################################################

sub print_arr{
	my $arr_ref=shift;
	foreach my $val (@{$arr_ref}){
		print STDERR "--> $val\n";
	}
}

open(IN_FH, "<$file") || die "Could not open $file\n";

$_=<IN_FH>;
chomp $_;
my @header_arr=split /\t/, $_;

my @data;

my $num_rows=0;
my %key_hash;

while(<IN_FH>){
	chomp;
	$_=~s/\r//g;
	my @columns=split /\t/, $_;

	if(lc($columns[$key_column]) eq "na" || $columns[$key_column] eq ""){
		next;	
	}

	my $keep_string;	
	my @kept;

	if($keep_all){
		for(my $i=0; $i<=$#columns; $i++){
			if($i!=$key_column){
				push @kept, $columns[$i];
			}
		}
	}else{
		for(my $i=0; $i<=$#keep_col; $i++){
			push @kept, $columns[$keep_col[$i]];
		}
	}

	my $keep_info=join "\t", @kept;

	# Make sure if there are duplicate keys, that the saved data is also identical
	if(defined($key_hash{$columns[$key_column]})){
		my $old_info=$key_hash{$columns[$key_column]};
		if($old_info ne $keep_info){
			print STDERR "\n";
			print STDERR "Error!!!:\n";
			print STDERR "\tKey: \"$columns[$key_column]\"\n";
			print STDERR "\tVersion 1: $old_info\n";
			print STDERR "\tVersion 2: $keep_info\n";
			print STDERR "\n";
			die("Redundant keys could not be collapsed because associated data was not redundant. Aborting.");
		}
	}

	$key_hash{$columns[$key_column]}=$keep_info;
	$num_rows++;
}
close(IN_FH);

###############################################################################

# Extract out header information
my @kept;

if($keep_all){
	for(my $i=0; $i<=$#header_arr; $i++){
		if($i!=$key_column){
			push @kept, $header_arr[$i];
		}
	}
}else{
	for(my $i=0; $i<=$#keep_col; $i++){
		push @kept, $header_arr[$keep_col[$i]];
	}
}
my $header_info=join "\t", @kept;

print STDOUT "$header_arr[$key_column]\t$header_info\n";

# Output kept data
foreach my $key (sort keys %key_hash){
	print STDOUT "$key\t$key_hash{$key}\n";
}

###############################################################################

print STDERR "Done.\n";
