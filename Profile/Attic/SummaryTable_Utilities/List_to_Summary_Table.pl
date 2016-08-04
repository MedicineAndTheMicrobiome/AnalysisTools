#!/usr/bin/env perl

use strict;
use Getopt::Std;
use vars qw($opt_i $opt_o $opt_c $opt_l $opt_a $opt_n);
getopts("i:o:c:l:a:n:");

my $usage = "usage:
$0
	-i <Input list of items to count>
	-o <output summary_table.tsv name root>
	[-c <column of item, starting from 0, default=0>]
	[-l <column of length, starting from 0>]
	[-a <read length, default=0 bp>]
	[-n <preferred sample name, default is the output file root>]

	Reads in a list of items and generates a count summary table
	and a length normalized abundance table, if lengths are available.

	Use the -a option if you want to use the normalization
	denominator of (subject_length-read_length).

";

if(!(
	defined($opt_i) && 
	defined($opt_o)
)){
	die $usage;
}

my $input_file_name=$opt_i;
my $output_file_name=$opt_o;
my $item_col=$opt_c;
my $preferred_sample_name=$opt_n;

my $length_col=$opt_l;
my $assay_length=$opt_a;

if(!defined($item_col)){
	$item_col=0;
}

if(!defined($length_col)){
	$length_col=undef;
}

if(!defined($assay_length)){
	$assay_length=0;
}

if(!defined($preferred_sample_name)){
	my @path=split "/", $output_file_name;
	$preferred_sample_name=$path[$#path];
}

print STDERR "Input Filename: $input_file_name\n";
print STDERR "Output Filename Root: $output_file_name\n";
print STDERR "\n";
print STDERR "Item column: $item_col\n";
print STDERR "Length column: $length_col\n";
print STDERR "Preferred Sample Name: $preferred_sample_name\n";
print STDERR "\n";
print STDERR "Assay Length: $assay_length\n";
		
##################################################################################

my %item_cts;
my %item_weights;
my $withlen=0;
my $nolen=0;
my $num_items=0;

open(IN_FH, "<$input_file_name") || die "Could not open $input_file_name.\n";

while(<IN_FH>){
	chomp;
	my @fields=split "\t", $_;
	
	my $item=$fields[$item_col];

	my $length=0;
	if(defined($length_col)){
		$length=$fields[$length_col];
	}
	
	if($item eq ""){
		next;
	}else{

		# Compute weight
		if($length>0){

			my $weight=1/($length-$assay_length);

			# Accumulate weighted hits
			my $val=$item_weights{$item};
			if(defined($val)){
				$val+=$weight;
			}else{
				$val=$weight;
			}
			$item_weights{$item}=$val;
			
			$withlen++;
		}

		# Accumulate counts
		my $val=$item_cts{$item};
		if(defined($val)){
			$val++;
		}else{
			$val=1;
		}
		$item_cts{$item}=$val;
		
		# Keep track of total
		$num_items++;
	}

}

close(IN_FH);

##################################################################################

print STDERR "Number items with lengths: $withlen\n";
print STDERR "Number of non-null items: $num_items\n";

my $generate_length_normalized=0;
if($num_items==$withlen && $num_items!=0){
	print STDERR "Length information for all items.  Normalizing...\n";
	$generate_length_normalized=1;
}else{
	print STDERR "Incomplete length information.  Cannot normalize counts by length.\n";
}

# Get categories
my @keys=sort keys %item_cts;
my $num_categories=$#keys+1;
print STDERR "Number of (unique) categories: $num_categories\n";

# Compute totals for both
my $norm_total=0;
my $cts_total=0;
foreach my $key(@keys){
	$norm_total+=$item_weights{$key};	
	$cts_total+=$item_cts{$key};
}

if($norm_total==0){
	print STDERR "WARNING: Zero counts...\n";
}

# Compute values for normalized
my @norms;
if($generate_length_normalized){
	print STDERR "Length normalized total: $norm_total\n";
	foreach my $key(@keys){
		push @norms, ($item_weights{$key}/$norm_total);
	}

	if($norm_total>0){
		$norm_total=1;
	}
}

# Compute values for counts
my @cts;
print STDERR "Count total: $cts_total\n";
foreach my $key(@keys){
	push @cts, $item_cts{$key};
}

##################################################################################
# Output counts and length normalized abundances

my $categories=join "\t", @keys;

# Output length normalized summary table
if($generate_length_normalized){
	my $normalized_fn="$output_file_name.wt_norm.summary_table.tsv";
	open(OUT_NORM_FH, ">$normalized_fn") || die "Could not open $normalized_fn.\n";
	print OUT_NORM_FH "Sample_ID\tTotal\t$categories\n";
	my $normalized_str=join "\t", @norms;
	print OUT_NORM_FH "$output_file_name\t$norm_total\t$normalized_str\n";
	close(OUT_NORM_FH);
}

# Output counts summary table
my $counts_fn="$output_file_name.cts.summary_table.tsv";
open(OUT_CTS_FH, ">$counts_fn") || die "Could not open $counts_fn.\n";
print OUT_CTS_FH "Sample_ID\tTotal\t$categories\n";
my $counts_str=join "\t", @cts;
print OUT_CTS_FH "$preferred_sample_name\t$cts_total\t$counts_str\n";
close(OUT_CTS_FH);

##################################################################################

print STDERR "Done.\n\n";
