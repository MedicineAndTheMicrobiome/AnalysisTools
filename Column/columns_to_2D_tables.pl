#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_c $opt_r $opt_n $opt_m $opt_i $opt_o);

getopts("c:r:n:m:i:o:");

my $usage = "
	usage:
	$0
		-c <column name for columns in output table, counting from 0>
		-r <column name for rows in output table, counting from 0>	

		-n <column names order list, comma separated in quotes>
		-m <row names order list, comma separated in quotes>

		-i <input file name, tab-separated-value>
		-o <output file name root>

	Will read in a multidimension keyed table and extract a 2D table for every column

	For example the input file:

	Site,PatientID,Shannon,Simpson,Evenness
	Subgingival_plaque,158418336,1.89621504536272,0.823471677623173,0.739279720999776
	Subgingival_plaque,159470302,1.88161816575524,0.826349634045539,0.784695723412212
	Supragingival_plaque,158822939,1.86429915856862,0.822572656637802,0.809654837183226
	Subgingival_plaque,159733294,1.86091777931096,0.813378447118344,0.72551833193021

	
";

if(!defined($opt_c) || 
	!defined($opt_r) || 
	!defined($opt_n) || 
	!defined($opt_m) || 
	!defined($opt_i) || 
	!defined($opt_o)){
	die $usage;
}

my $column_colname=$opt_c;
my $row_colname=$opt_r;
my @column_order=split /,/, $opt_n;
my @row_order=split /,/, $opt_m;
my $infile=$opt_i;
my $outfile=$opt_o;

my $num_rows=$#row_order+1;
my $num_cols=$#column_order+1;

print STDERR "Num Expected Rows: $num_rows\n";
print STDERR "Num Expected Columns: $num_cols\n";

###############################################################################

open(IN_FH, "<$infile") || die "Could not open map file $infile\n";

$_=<IN_FH>;
my @headers=split /,/, $_;

###############################################################################

my $column_idx=-1;
my $row_idx=-1;
for(my $i=0; $i<=$#headers; $i++){
	#print "$headers[$i]\n";
	if($headers[$i] eq $column_colname){
		$column_idx=$i;
	}
	if($headers[$i] eq $row_colname){
		$row_idx=$i;
	}
}

if($column_idx==-1){
	die "Could not identify column that matched $column_colname\n";
}else{
	print STDERR "Found $column_colname in column $column_idx\n";
}

if($row_idx==-1){
	die "Could not identify column that matched $row_colname\n";
}else{
	print STDERR "Found $row_colname in column $row_idx\n";
}

###############################################################################

# Read input 
my @table_hash;

while(<IN_FH>){
	chomp;
	my @in=split /,/, $_;

	my $num_tables=0;
	for(my $i=0; $i<=$#headers; $i++){
		my $col=$in[$column_idx];
		my $row=$in[$row_idx];
		my $data=$in[$i];
		${$table_hash[$i]}{$row}{$col}=$data;
		#print STDERR "Saving [$i]{$row}{$col}=$data\n";
	}
}

close(IN_FH);

###############################################################################

# Write output


my $colkey_str=join "\t", @column_order;

for(my $i=0; $i<=$#headers; $i++){

	if($i==$column_idx || $i==$row_idx){
		next;
	}

	open(OUT_FH, ">$outfile\.$headers[$i]") || die "Could not open $outfile\.$headers[$i]\n";

	print STDERR "Working on column $i, ie. $headers[$i]\n";

	print OUT_FH "$colkey_str\n";

	foreach my $rowkey (@row_order){
		my @out_arr;
		push @out_arr, $rowkey;
		foreach my $colkey (@column_order){
			push @out_arr, "${$table_hash[$i]}{$rowkey}{$colkey}";
			#print STDERR "Writing [$i]{$rowkey}{$colkey}=${$table_hash[$i]}{$rowkey}{$colkey}\n";
		}
		my $outstr=join "\t", @out_arr;
		print OUT_FH "$outstr\n";
	} 


}

