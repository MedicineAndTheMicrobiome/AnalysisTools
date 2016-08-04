#!/usr/bin/env perl

use strict;
use Getopt::Std;
use File::Basename;
use FileHandle;
use vars qw($opt_i $opt_o);
getopts("i:o:");

my $usage = "usage:
$0
	-i <input summary_table file>
	-o <output summary_table root filename>

	Reads in a summary file table where the categories are actually semicolon (;) separated lists.
	The lists are then split, and the counts underneath them are redistributed.

	For example, if we have the following categories and counts:

	    total	a;b	a;d
		9	4	5

	Then the output will have the following categories and counts:
	(*.cumulative.summary_table.tsv)

	    total	a	b	d
		18	9	4	5	

	and (even):
	(*.distributed.summary_table.tsv)

	    total	a	b	d
		9	2+2.5	2	2.5


";

if(!(
	defined($opt_i) || 
	defined($opt_o)
)){
	die $usage;
}

my $input_file_name=$opt_i;
my $output_file_name=$opt_o;

print STDERR "Input Filename: $input_file_name\n";
print STDERR "Output Filename Root: $output_file_name\n";
		
################################################################################

open(OUT_CUMU_FH, ">$output_file_name\.cumulative.summary_table.tsv") || 
	die "Could not open > $output_file_name\.cumulative.summary_table.tsv";

open(OUT_DIST_FH, ">$output_file_name\.distributive.summary_table.tsv") || 
	die "Could not open > $output_file_name\.distributive.summary_table.tsv";

#------------------------------------------------------------------------------

open(ST_FH, "<$input_file_name") || die "Could not open < $input_file_name\n";

# Extract header
my $header=<ST_FH>;
chomp $header;

my @header_fields=split "\t", $header;

# Skip over sample_id and counts fields
my $sample_id_field=shift @header_fields;
my $counts_field=shift @header_fields;

my $num_fields=$#header_fields+1;

# Determine which subcategories there will be when categories are split
my %category_hash;
foreach my $category(@header_fields){
	my @subcats=split ";", $category;
	foreach my $subcat(@subcats){
		$category_hash{$subcat}=1;
	}
}
my @subcats_arr=sort keys %category_hash;

# Output header line
my $subheaders_str=join "\t", @subcats_arr;
print OUT_CUMU_FH "Sample_ID\tTotal\t$subheaders_str\n";
print OUT_DIST_FH "Sample_ID\tTotal\t$subheaders_str\n";

# Go through each sample and split out category counts
while(<ST_FH>){
	
	chomp;

	my @count_fields=split "\t", $_;
	my $sample_id=shift @count_fields;
	my $total=shift @count_fields;

	# initialize hash	
	my %cumulative_hash;
	my %distributive_hash;
	foreach my $subcat(@subcats_arr){
		$cumulative_hash{$subcat}=0;
		$distributive_hash{$subcat}=0;
	}

	# add up counts
	for(my $i=0; $i<$num_fields; $i++){
		my @subfields=split ";", $header_fields[$i];
		my $num_subfields=$#subfields+1;
		
		my $dist_count=$count_fields[$i]/$num_subfields;

		foreach my $subfield(@subfields){
			$cumulative_hash{$subfield}+=$count_fields[$i];
			$distributive_hash{$subfield}+=$dist_count;
		}
	}

	# Addup counts and put them in the right order
	my @cum_outcounts;
	my @dis_outcounts;
	my $cum_sum=0;
	my $dis_sum=0;
	foreach my $subcat(@subcats_arr){
		push @cum_outcounts, $cumulative_hash{$subcat};
		push @dis_outcounts, $distributive_hash{$subcat};
		$cum_sum+=$cumulative_hash{$subcat};
		$dis_sum+=$distributive_hash{$subcat};
	}

	# Output counts
	my $cum_outline=join "\t", @cum_outcounts;
	my $dis_outline=join "\t", @dis_outcounts;

	print OUT_CUMU_FH "$sample_id\t$cum_sum\t$cum_outline\n";
	print OUT_DIST_FH "$sample_id\t$dis_sum\t$dis_outline\n";

}

##################################################################################

close(OUT_CUMU_FH);
close(OUT_DIST_FH);
close(ST_FH);

print STDERR "Done.\n";



