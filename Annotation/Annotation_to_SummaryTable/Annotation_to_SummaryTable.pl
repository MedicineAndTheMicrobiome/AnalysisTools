#!/usr/bin/env perl

use strict;
use Getopt::Std;
use FileHandle;
use vars qw ($opt_l $opt_o);

getopts("l:o:");

my $usage = "
usage:
    $0

	-l <List of Taxa annotation output files>
	-o <Output filename root>
	
	This script will read in a list of annotation files.

	The list (-l) file should look like:

		<sample name> \\t <taxa annotation path>

	The output file name should be the root name of the summary
	tables to generate.

	This script will run through each of the samples in the
	list and generate a summary table for each of the columns
	in the annotation file.

	Example output would be:
		
		<output file name root>.Description.summary_table.tsv
		<output file name root>.PFamIDs.summary_table.tsv
		<output file name root>.TIGRFamIDs.summary_table.tsv
		<output file name root>.GOProcIDs.summary_table.tsv
		<output file name root>.GOFuncIDs.summary_table.tsv
		<output file name root>.ECIDs.summary_table.tsv

	The input annotation file should look like:

		1.) Read ID
		2.) PerComp ID
		3.) Length
		4.) Description
		5.) PFamIDs
		6.) TIGRFamIDs
		7.) GOProcIDs
		8.) GOFuncIDs
		9.) ECIDs


";

###############################################################################

if(
	!defined($opt_l) ||
	!defined($opt_o)
){
	die $usage;
}

my $SampleFileList=$opt_l;
my $OutputFnameRoot=$opt_o;

###############################################################################

sub read_sample_file_list{
	my $fname=shift;
	
	print STDERR "Loading Sample-to-File List: $fname\n";

	open(FH, "<$fname") || die "Could not open $fname.\n";

	my %fn_hash;

	while(<FH>){
		chomp;
		my ($sample_name, $path)=split "\t", $_;
		$fn_hash{$sample_name}=$path;	
	}

	close(FH);
	
	return(\%fn_hash);
}

sub print_arr{
	my $arr_ref=shift;
	foreach my $val(@{$arr_ref}){
		print STDERR "  $val\n";
	}
}

sub load_data{
	my $fname=shift;
	my $counts_hash_ref=shift;
	my $start_col=shift;

	open(FH, "<$fname") || die "Could not open annotation file: $fname.\n";

	print STDERR "Starting from column: " . ($start_col+1) . "\n";

	#print STDERR $counts_hash_ref . "\n";

	my $line=0;
	my @function_types;
	my $num_col;
	while(<FH>){

		chomp;
		my @col_arr=split "\t", $_;

		if($line==0){
			@function_types=@col_arr;
			$num_col=($#function_types+1);
			print STDERR "Num columns detected: $num_col\n";
		}else{
			for(my $i=$start_col; $i<$num_col; $i++){

				my $category=$col_arr[$i];

				if(!defined(${${$counts_hash_ref}{$function_types[$i]}}{$category})){
					${${$counts_hash_ref}{$function_types[$i]}}{$category}=1;
				}else{
					${${$counts_hash_ref}{$function_types[$i]}}{$category}++;
				}
			}
		}	
	
		$line++;
	}	

	close(FH);

}

###############################################################################


my $target_files_hash_ref=read_sample_file_list($SampleFileList);
my @sample_ids=sort keys %{$target_files_hash_ref};

my %counts_hash;

my $START_FUNCT_COL=3;
# $counts{$sample_id}{$function}{$category};

foreach my $sample_id(@sample_ids){
	print STDERR "Working on: $sample_id\n";
	my $path=${$target_files_hash_ref}{$sample_id};

	print STDERR "Path: $path\n";

	my %functions_hash;
	$counts_hash{$sample_id}=\%functions_hash;
	load_data($path, $counts_hash{$sample_id},$START_FUNCT_COL);

}

my %function_types_hash;

foreach my $sample_id(@sample_ids){
	my @function_types=keys %{$counts_hash{$sample_id}};
	foreach my $function_type(@function_types){
		$function_types_hash{$function_type}=1;
	}
}

my @unique_function_types=sort keys %function_types_hash;
print STDERR "\nFunctions Types Discovered: \n";
print_arr(\@unique_function_types);
print STDERR "\n";

#${${$counts_hash_ref}{$function_types[$i]}}{$category}++;

foreach my $function_type(@unique_function_types){
	print STDERR "Collecting data across $function_type...\n";
	
	print STDERR "Determining categories across all samples for $function_type\n";
	my %categories_hash;
	foreach my $sample_id (@sample_ids){
		my $function_hash_ref=${$counts_hash{$sample_id}}{$function_type};
		my @categories=keys %{$function_hash_ref};

		foreach my $category(@categories){
			$categories_hash{$category}=1;
		}
	}

	my @categories=sort keys %categories_hash;
	my $num_cat=$#categories+1;
	print STDERR "Found $num_cat categories.\n";

	my $output_fn="$OutputFnameRoot\.$function_type\.summary_table.tsv";
	print STDERR "Writing $output_fn...\n";
	open(FH, ">$output_fn") || die "Could not open $output_fn\n";

	
	my $header=join "\t", @categories;
	print FH "sample_id\ttotal\t" . $header . "\n";

	foreach my $sample_id (@sample_ids){
		my $tot=0;
		my @values;

		foreach my $category(@categories){

			if($category eq "" || $category eq "NA"){
				next;

			}
			my $counts=${${$counts_hash{$sample_id}}{$function_type}}{$category};
			if($counts eq ""){
				$counts=0;
			}
			$tot+=$counts;
			push @values, $counts;
		}
			
		print STDERR "  $sample_id: $tot\n";
		unshift @values, $tot;
		unshift @values, $sample_id;

		my $outline=join "\t", @values;
		print FH "$outline\n";
	}

	close(FH);
	print STDERR "\n";

}


print STDERR "Done.\n";

