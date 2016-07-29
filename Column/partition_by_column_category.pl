#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_i $opt_l $opt_d $opt_c $opt_m $opt_s $opt_o);

getopts("i:l:dcm:s:o:");

my $MAX_CATEGORIES=30;

my $usage = "
	usage:
	$0
		-i <input file>
		-l <column, counting from 0>

		[-d <flag to generate id list, named <input file>.<id_list> >]
		[-c <flag to include category as first item of id list requsted in option -d >]
		[-m <max categories allowed, default $MAX_CATEGORIES>]
		[-s <field separator, default is TAB>]
		[-o <output directory>]

	Reads in input file and based on the column you specify
	splits the file into multiple files using the value in that 
	column as a key.

	Output will be written to <input file>.<category>  
	
	Where category depends on what is in the column.  You can
	override the max columns, but it is there just in case you
	specify the wrong column accidently and the column happens
	to contain a large set of possible values.

";

if(!defined($opt_i) || !defined($opt_l)){
	die $usage;
}

my $file=$opt_i;
my $column=$opt_l;
my $max_categories=$opt_m;
my $make_list=defined($opt_d);
my $include_category_in_list=defined($opt_c);
my $field_sep=$opt_s;
my $output_dir=$opt_o;

if(!defined($opt_m)){
	$max_categories=$MAX_CATEGORIES;
}

if(!defined($opt_s)){
	$field_sep="\t";
}

###############################################################################

my $outfname_root="$file";
my %partitions;

open(IN_FH, "<$file") || die "Could not open $file\n";
while(<IN_FH>){
	chomp;
	my @array=split /$field_sep/, $_;
	push @{$partitions{$array[$column]}}, $_;
}
close(IN_FH);

###############################################################################

my @categories=sort keys %partitions;
if(($#categories+1) > $max_categories){
	print STDERR "Max categories exceeded: " . ($#categories+1) . " > $max_categories\n";
	die "Not letting you shoot yourself in the foot by creating too many files.";
}

###############################################################################

foreach my $partition(@categories){

	if(!defined($output_dir)){
		open (OUT_FH, ">$outfname_root\.$partition") || die "Could not open $outfname_root\.$partition\n";
	}else{
		open (OUT_FH, ">$output_dir/$partition") || die "Could not open $output_dir/$partition\n";
	}

	foreach my $line(@{$partitions{$partition}}){
		print OUT_FH "$line\n";
	}
	close(OUT_FH);

	if($make_list){
		if($include_category_in_list){
			`echo $partition > $outfname_root\.$partition\.id_list`;
		}
		`cut -f 1 $outfname_root\.$partition >> $outfname_root\.$partition\.id_list`;
	}
}

###############################################################################
