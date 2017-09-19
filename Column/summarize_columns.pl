#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_d $opt_m);

getopts("i:d:m:");

my $MAX_UNIQ=5;

my $usage = "
	usage:
	$0
		[-i <input file name>]
		[-d <delimiter, default=tab>];
		[-m <maximum unique to list, default=$MAX_UNIQ>]


	For each column the following summary will be produced:
		1.) column number
		2.) column header name
		3.) categories (counts)

	If there are more categories than $MAX_UNIQ, then ... will be used to summarize
	the remaining.

	For example:
		
		197 smdsc.excoag:  No(318), N/A(47), Yes(3)
		198 smdsc.clinpres:  non-acute(329), acute(39)
		199 smdsc.scadding:  II(108), IV(83), I(80), III(47), 0(46), N/A(4), ...(0)
		200 smdsc.treatsta:  untreated(205), treated(163)
		201 smdsc.multiorg:  No(315), Yes(52), N/A(1)


	Results goes to STDOUT.
	
";



if(!defined($opt_i) && (-t STDIN) && $ARGV[0] eq ""){
	die $usage;
}

my $File;
if($ARGV[0] eq ""){
	$File=$opt_i;
}else{
	$File=$ARGV[0];
}

my $MaxUnique=$MAX_UNIQ;
if(defined($opt_m)){
	$MaxUnique=$opt_m;
}

my $Delim;
if(defined($opt_d)){
	$Delim=$opt_d;
}else{
	$Delim="\t";
}

###############################################################################

if(!(-t STDIN)){
	$_=<STDIN>;
}else{
	open(IN_FH, "<$File") || die "Could not open $File\n";
	$_=<IN_FH>;
}

chomp $_;
my @header_arr=split /$Delim/, $_, -1;
my $num_hdr_cols=$#header_arr+1;

my %col_hash;

my $line=2;
while(<IN_FH>){
	chomp;

	my @cols=split /$Delim/, $_, -1;
	my $num_cols=$#cols+1;

	if($num_cols != $num_hdr_cols){
		print STDERR "WARNING/ERROR: Line $line (starting from 2) has $num_cols columns, but header has $num_hdr_cols columns (line 1)\n";
	}

	for(my $i=0; $i<$num_cols; $i++){
		push @{$col_hash{$header_arr[$i]}}, $cols[$i];
	}

	$line++;
}

for(my $i=0; $i<$num_hdr_cols; $i++){
	my $colname=$header_arr[$i];
		
	# Count categories
	my %cat_hash;
	my @row_arr=@{$col_hash{$colname}};
	foreach my $val (@row_arr){
		$cat_hash{$val}++;
	}

	# Place counts in parallel array
	my @keys_arr=keys %cat_hash;
	my @counts_arr;
	foreach my $key	(@keys_arr){
		push @counts_arr, $cat_hash{$key};	
	}
	my $num_keys=$#keys_arr+1;
	my @seq=0..($num_keys-1);

	# Compute sort order (decreasing from largest counts)
	my @order=sort {$counts_arr[$b] <=> $counts_arr[$a]} @seq;

	# Limit display
	my $disp;
	if($num_keys>=$MaxUnique){
		$disp=$MaxUnique;
	}else{
		$disp=$num_keys;
	}
		
	# Generate output in array
	my @summary;
	my $remaining=0;
	for(my $x=0; $x<$num_keys; $x++){
		my $ix=$order[$x];
		my $key=$keys_arr[$ix];
		if($key eq ""){
			$key="(null)";
		}
		if($x<$disp){
			push @summary, "$key($counts_arr[$ix])";
		}else{
			$remaining+=$counts_arr[$ix];
		}
	}

	# Concatenate with commas
	if($num_keys>$MaxUnique){
		push @summary, "...($remaining)";
	}
	my $out_string= join ", ", @summary;

	# Output string
	my $colnum=$i+1;
	print "$colnum $colname:  $out_string\n";
}

