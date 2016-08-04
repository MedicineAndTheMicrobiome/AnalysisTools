#!/usr/bin/env perl

#######################################################################

use strict;
use Getopt::Std;
use File::Basename;
use FileHandle;
use vars qw($opt_i $opt_s);
getopts("i:s:");


my $usage = "usage:
$0
	Split summary table into two files.
		-i <input summary table>
		-s <comma separate split string, eg. donorA,donorB>

	Maintains column headers and confirms that all listed split strings are disjoint.

";

if(!(
	defined($opt_i) &&
	defined($opt_s)
)){
	die $usage;
}

my $input_file_name=$opt_i;
my @strings=split /,/, $opt_s;

print STDERR "Input Filename: $input_file_name\n";
my $basename=$input_file_name;
$basename=~s/\.summary_table\.xls$//;

my @out_names;
my @filehandles;
my $i=0;
foreach my $split_str (@strings){
	my $out_names=$basename . "." . $split_str . ".summary_table.xls";
	print "Output Files: '$out_names'\n";
	$filehandles[$i]=FileHandle->new;
	$filehandles[$i]->open(">$out_names");
	$i++;
}
		
my $linenum=0;
my $split_count=0;
open(INFH, "<$input_file_name") || die "Could not open $input_file_name\n";

while(<INFH>){

	if($linenum==0){
		foreach my $fh(@filehandles){
			print {$fh} $_;
		}
	} else{
		my @fields=split /\t/, $_;

		for(my $i=0; $i<=$#strings; $i++){
			my $search_string=$strings[$i];
			if($fields[0]=~/$search_string/){
				print {$filehandles[$i]} $_;
				$split_count++;
			}
		}
	}
	$linenum++;
}

my $num_samples=$linenum-1;
if(($split_count) == $num_samples){
	print STDERR "Split successful.  \n";
}else{
	if($split_count>$num_samples){
		die "Error!  Your split strings are not exclusively separating your samples.\n";
	}else{
		print STDERR "Warning.  Not all of your samples were assignable to a split destination.\n";
	}
}

close(INFH);
