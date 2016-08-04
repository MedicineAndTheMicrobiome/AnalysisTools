#!/usr/bin/perl

#use warnings;
use strict;
use Getopt::Std;
use vars qw($opt_m $opt_o);
#use Data::Dumper;


getopts("m:o:");

my $usage = "usage:
$0
	-m <METAREP_file>
	[-o <output name>]

	This script will pull out abundance information from METAREP file.
\n";

if(!(
	defined($opt_m) 
)){
	die $usage;
}

my $input_file=$opt_m;
my $output_name=$opt_o;

if(!defined($output_name)){
	$output_name= $input_file.".summary.xls";
}

###################################################################################
open (MREP, $input_file) or
	die "\n!ERROR!: Could not open METAREP file for reading\n";

open (OUT, ">$output_name") or
	die "\n!ERROR!: Could not open summary file for writing\n";

my $line; #stores a line read from the METAREP data file
my @temp; #stores the contents of $line split by \t
my $num_samples; #stores number of samples present in the METAREP datafile
my @sample_names; #stores sample names from the METAREP data file
my $counter; #counter variable
my %summary; #contains all the information of a summary table except the total count
my %taxa; #stores the taxa names from the METAREP data file
my %sample_total; #stores the total number of reads in a sample

while($line = <MREP>)
{
	chomp($line);
	@temp = split(/\t/, $line);

	if($temp[0] eq "ID")
	{
		$num_samples = scalar(@temp) - 3; #ID, Category, Total
		for($counter = 3; $counter < ($num_samples + 3); $counter++)
		{
			push(@sample_names, $temp[$counter-1]);
		}
	}
	elsif($temp[0] =~ /^[0-9]/)
	{
		for($counter = 3; $counter < ($num_samples + 3); $counter++)
		{
			$temp[1] =~ s/\s+/_/g; #remove spaces in taxon names and replace with underscore
			#$temp[1] = ucfirst($temp[1]); #to capitalize the first letter of a taxon
			$summary{$sample_names[$counter-3]}{$temp[1]} = $temp[$counter-1]; 
			$taxa{$temp[1]} = 0;
			$sample_total{$sample_names[$counter-3]} += $temp[$counter-1];
		}
	}
}

print OUT "Sample\tTotal";

foreach my $taxon (sort keys %taxa)
{
	print OUT "\t$taxon";
}

print OUT "\n";

foreach my $sample (sort keys %summary)
{
	print OUT "$sample\t$sample_total{$sample}\t";

	my $new_sample_flag = 1;

	foreach my $taxon (sort keys %{$summary{$sample}})
	{
		if($new_sample_flag == 1)
		{
			print OUT "$summary{$sample}{$taxon}";
			$new_sample_flag = 0;
		}
		else
		{
			print OUT "\t$summary{$sample}{$taxon}";
		}
		
	}
	print OUT "\n";
}