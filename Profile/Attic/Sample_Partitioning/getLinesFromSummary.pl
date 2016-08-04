#!/usr/local/bin/perl

use strict;
use Getopt::Std;
use File::Basename;
use FileHandle;
use vars qw($opt_i $opt_p);
getopts("i:p:");

my $usage = "usage:
$0
	Split summary table into two files.
		-i <input list file>
		-p <project name>
        reads in files:

                <project>.domain.summary_table.xls
                <project>.phylum.summary_table.xls
                <project>.class.summary_table.xls
                <project>.order.summary_table.xls
                <project>.family.summary_table.xls
                <project>.genus.summary_table.xls

        Generates a table for input file. Maintains column headers and lists all samples in input file.
       Writes out:
                <inputFile>.<level>.summary

";

if(!(
	defined($opt_i) &&
	defined($opt_p)
)){
	die $usage;
}

my $input_file=$opt_i;
my $project_id=$opt_p;

print STDERR "Input Filename: $input_file\n";
print STDERR "Project Filename:$project_id\n";
 
my @levels=("domain", "phylum", "class", "order", "family", "genus");

my @sample_list;

# Read input file without wasting time
open(INFH, "<$input_file") || die "Could not open $input_file\n";

while(<INFH>){
    chomp;
    push @sample_list, $_;
}
close (INFH);


foreach my $level(@levels){

        my $input_file_name="$project_id\.$level\.summary_table.xls";
	my $output_file_name="$input_file\.$level\.summary";

        open(INFH, "<$input_file_name") || die "Could not open $input_file_name\n";
	open(OUTFH, ">$output_file_name")|| die "Could not open $output_file_name\n";

	my $count=0;
	my $num=0;
	while(<INFH>){
	    if ($num==0){
		print OUTFH $_;
	    }
	    else {
		my @fields=split /\t/, $_;
		for(my $i=0; $i<=$#sample_list; $i++){
		    my $search_string=$sample_list[$i];
		    if ($fields[0]=~/$search_string/){
			print OUTFH $_;
			$count++;
		    }
		}
	    }
	    $num++;
	}

	print STDOUT "$count samples found at $level out of '$#sample_list+1'\n";
	close(INFH);
	close(OUTFH);
}


