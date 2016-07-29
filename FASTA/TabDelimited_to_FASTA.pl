#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_t);

getopts("t:");
my $usage = "usage: 
$0 
	-t <tab delimited file>

	This script will convert a tab delimited file into a fasta file.
	The format of the tab delimited file should be:

	<id1>\\t<sequence1>\\t...\\t<idn>\\t<sequencen>\\n

	You can can 1 to n, pairs of id/sequence in each row.

	Output will go to stdout.

	This script is useful for converting a file that has been created
	with a spreadsheet program, where the id are in adjacent columns/cells
	to the sequence.  
	
";

if(!defined($opt_t)){
	die $usage;
}

my $input_file=$opt_t;

###############################################################################

print STDERR "Reading TabSeparatedValue file...\n";

open(FH, "<$input_file") || die "Could not open $input_file\n";

while(<FH>){

	chomp;
	my @arr=split /\t/, $_;

	my $num_col=$#arr+1;
	if($num_col%2 > 0){
		die "The number of columns is not a multiple of 2 : $num_col\n";
	}

	for(my $i=0; $i<$num_col; $i+=2){
		$arr[$i]=~s/^\s+//;
		$arr[$i]=~s/\s+$//;
		$arr[$i+1]=~s/^\s+//;
		$arr[$i+1]=~s/\s+$//;
		output_fasta_record($arr[$i], $arr[$i+1]);
	}

}

close(FH);

print STDERR "done.\n";


###############################################################################

sub output_fasta_record{

	my $defline=shift;
	my $sequence=shift;

	print ">$defline\n";
	my $length=length($sequence);
	my $width=60;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print substr($sequence, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);

}

#------------------------------------------------------------------------------
