#!/usr/bin/env perl

use strict;
use Getopt::Std;
use File::Basename;
use FileHandle;
use vars qw($opt_i $opt_l $opt_s);
getopts("i:l:s:");


my $usage = "usage:
$0
	Split summary table into two files.
		-i <input summary table>
		-l <sample list file that need to be included in new summary file>
		-s <'string' to otuput file name >

	Maintains column headers and confirms that all listed sample names are disjoint.

";

if(!(
	defined($opt_i) &&
	defined($opt_l) &&
	defined($opt_s)
)){
	die $usage;
}

my $input_file_name=$opt_i;
my $list_file=$opt_l;
my $string=$opt_s;

print STDERR "Input Filename: $input_file_name\n";
my $basename=$input_file_name;
$basename=~s/\.summary_table\.xls$//;

my @out_names;
my @filehandles;
my $i=0;
open(FILE,"$list_file") || die "Could not open $list_file\n";

my @strings = <FILE>;
my $out_names=$basename . "." . $string .".summary_table.xls";
print "Output File: '$out_names'\n";
my $filehandle=FileHandle->new;
$filehandle->open(">$out_names");
		
my $linenum=0;
my $split_count=0;
open(INFH, "<$input_file_name") || die "Could not open $input_file_name\n";

while(<INFH>){

	if($linenum==0){
			print {$filehandle} $_;
	} else{
		my @fields=split /\t/, $_;

		for(my $i=0; $i<=$#strings; $i++){
			my $search_string=$strings[$i];
	chomp $search_string;
			if($fields[0]=~/$search_string/){
				print {$filehandle} $_;
				$split_count++;
			}
		}
	}
	$linenum++;
}

my $num_samples=$linenum-1;
print STDERR "Number of samples found $split_count out of $num_samples.\n";

close(FILE);
close(INFH);
