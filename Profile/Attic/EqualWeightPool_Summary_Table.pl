#!/usr/bin/env perl

#######################################################################

use strict;
use Getopt::Std;
use File::Basename;
use FileHandle;
use vars qw($opt_i $opt_o $opt_a);
getopts("i:o:a:");


my $usage = "usage:
$0
		-i <input summary_table.xls name>
		[-o <output summary_table.ewp.xls name>]
		[-a <altnernate pool sample name>]

	Reads in a summary table, then then pools all the data together with an equal weighting.
	When completed there will be one line that sums up to 1.

	If -o is not specified, the output file name will be based on the input file name, 
	with the .ewp.xls extension attached.

	If -a is not specified, the sample name used in the EQW summary table will be base on the
	input file name, without the .summary_table.xls extension and without the path.

";

if(!(
	defined($opt_i)
)){
	die $usage;
}

my $input_file_name=$opt_i;

# Get a new sample name
my ($input_file_name_base)=fileparse($input_file_name);
$input_file_name_base=~s/\.summary_table\.xls$//;
my $pooled_sample_name=$input_file_name_base;
if(defined($opt_a)){
	$pooled_sample_name=$opt_a;	
}

# Get output file name
my $output_file_name="$pooled_sample_name.summary_table.ewp.xls";
if(defined($opt_o)){
	$output_file_name=$opt_o;
}

print STDERR "Input Filename: $input_file_name\n";
print STDERR "Output Filename: $output_file_name\n";
		
open(INFH, "<$input_file_name") || die "Could not open $input_file_name\n";

# Read in categories
my $hdr_line=<INFH>;
my @categories=split /\t/, $hdr_line;
my $num_categories=$#categories-1;
print STDERR "Num categories: $num_categories\n";

# Read in counts
my @samples;
while(<INFH>){
	chomp;
	my @counts=split /\t/, $_;
	shift @counts; # Shift out sample nametitle
	shift @counts; # Shift out totals
	push @samples, \@counts;
}
close(INFH);
my $num_samples=$#samples+1;
print STDERR "Num samples: $num_samples\n";

# Normalize counts within sample
my @sum_of_counts_across_samples;
for(my $i=0; $i<$num_categories; $i++){
	$sum_of_counts_across_samples[$i]=0;
}
my $sum_of_counts_across_samples_ref=\@sum_of_counts_across_samples;

for(my $i=0; $i<$num_samples; $i++){
	my $sum=sum($samples[$i]);
	my $norm_ref=divide($samples[$i], $sum);
	$sum_of_counts_across_samples_ref=sum_arrays($sum_of_counts_across_samples_ref, $norm_ref);
}

# Average sample percentages
my $ewp_ref=divide($sum_of_counts_across_samples_ref, $num_samples);

my $sanity_check_sum=sum($ewp_ref);
if(abs($sanity_check_sum - 1)>0.01){
	print STDERR "Warning:  Sum of normalized counts doesn't exactly equal 1.\n";
}


# Output EWP summary table
open(OUTFH, ">$output_file_name") || die "Could not open $output_file_name\n";

print OUTFH $hdr_line;
my $outline=join "\t", @{$ewp_ref};
print OUTFH "$pooled_sample_name\t$sanity_check_sum\t$outline\n";

close(OUTFH);

print STDERR "Done.\n";

##################################################################################

sub print_array{
	my $arr_ref=shift;
	foreach my $val(@{$arr_ref}){
		print "$val\t";		
	}
	print "\n";
}

sub sum_arrays{
	my $arr1_ref=shift;
	my $arr2_ref=shift;
	my @new_arr;
	
	if($#{$arr1_ref} != $#{$arr2_ref}){
		print STDERR ("Arr1: " . ($#{$arr1_ref}+1) . "\n");
		print STDERR ("Arr2: " . ($#{$arr2_ref}+1) . "\n");
		print STDERR  "Arrays are not the same length.\n";
	}
	
	for(my $i=0; $i<=$#{$arr1_ref}; $i++){
		$new_arr[$i]=${$arr1_ref}[$i]+${$arr2_ref}[$i];
	}
	return(\@new_arr);
}

sub sum{
	my $arr_ref=shift;
	my @arr=@{$arr_ref};
	my $sum=0;
	foreach my $val(@arr){
		$sum+=$val;
	}
	return($sum);
}

sub divide{
	my $arr_ref=shift;
	my $divisor=shift;
	my @arr=@{$arr_ref};
	my @new_arr;
	for(my $i=0; $i<=$#arr; $i++){
		push @new_arr, ($arr[$i]/$divisor);
	}
	return(\@new_arr);
}



