#!/usr/bin/env perl

###############################################################################

use strict;
use File::Basename;
use Getopt::Std;
use vars qw ($opt_i $opt_k $opt_o $opt_n $opt_r $opt_u);

getopts("i:k:o:nru");

my $usage = "
	usage:
	$0
		Input File Parameters:
		-i <input file name>
		-k <targeted column name>
		
		Output File Name:
		-o <output file name>

		Options to pass onto sort:
		[-n (numeric sort)]
		[-r (reverse sort)]
		[-u (keep unique)]

	This script uses linux sort to sort the input file,
	but makes the call easier to follow by allowing the
	user to specify the column name.  This also prevents
	the linux sort from trying to sort the first row of the
	file.

";

if(
	!defined($opt_i) || 
	!defined($opt_k) || 
	!defined($opt_o) 
){
	die $usage;
}

my $InputFile=$opt_i;
my $InputColname=$opt_k;
my $OutputFile=$opt_o;

my $NumericOpt=($opt_n eq "")?0:1;
my $ReverseOpt=($opt_r eq "")?0:1;
my $UniqueOpt=($opt_u eq "")?0:1;

#------------------------------------------------------------------------------


print STDERR "\n";
print STDERR "Running: $0\n";
print STDERR "   Input File: $InputFile\n";
print STDERR "   Targeted Column Name: $InputColname\n";
print STDERR "   Output File Name: $OutputFile\n";
print STDERR "\n";
print STDERR "     Numeric: $NumericOpt\n";
print STDERR "     Reverse: $ReverseOpt\n";
print STDERR "     Unique:  $UniqueOpt\n";
print STDERR "\n";

###############################################################################
# Load the mapping file

open(INPUT_FH, "<$InputFile") || die "Could not open input file: $InputFile\n";

my $tar_fidx=-1;
my $header;

while(<INPUT_FH>){
	chomp;
	$header=$_;
	#print "HEADER: $header\n";
	my @fields=split "\t", $header;	

	my $fidx=1;
	foreach my $fld (@fields){
		if($InputColname eq $fld){
			print "Target Column Found: $fidx (Starting from 1)\n";
			$tar_fidx=$fidx;
		}
		$fidx++;
	}

	if($tar_fidx==-1){
		die "Error: Could not find targeted Column Name.\n";
	}

	last;	
}

close(INPUT_FH);

open(OUT_FH, ">$OutputFile") || die "Could not open output file: $OutputFile for writing.\n";
print OUT_FH "$header\n";
close(OUT_FH);

my $sort_opt;
if($NumericOpt || $ReverseOpt || $UniqueOpt){
	$sort_opt="-" . 
		($NumericOpt?"n":"") .
		($ReverseOpt?"r":"") .
		($UniqueOpt?"u":"");
}else{
	$sort_opt="";
}

my $sort_command="sort $sort_opt -k $tar_fidx";
print STDERR "Sort Command: '$sort_command'\n";

my $sys_res=system("tail -n +2 $InputFile | $sort_command >> $OutputFile");

if($sys_res){
	print STDERR "An error occurred on sort/append: $sys_res.\n";
}

print STDERR "Done.\n";
