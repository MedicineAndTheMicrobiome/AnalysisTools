#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_o);

getopts("o:");

my $usage = "
	usage:
	$0
		-o <output file name> <range file1> <range file2> ... <range filen>

	takes in list of clear ranges and produces a clear range
	that is the overlap of them all.

	The format of the range file is:

	<sequence identifier>\\tbegin\\tend\\n

	where begin and end should be in 0-space-based coordinates.

	Anything after the 3rd column will be ignored.


";

if(!defined($opt_o)){
	die $usage;
}

my $out_file=$opt_o;

###############################################################################

open(OUT_FH, ">$out_file") || die "Could not open $out_file";

my $in_file;

my %begins;
my %ends;

while($in_file=shift){
	print STDERR "Working on $in_file\n";
	open(IN_FILE_FH, "<$in_file") || die "Could not open $in_file\n";

	while(<IN_FILE_FH>){

		chomp;
		my ($seq_id, $begin, $end)=split /\t/, $_;

		if($begin>$end){
			die "Unexpected coordinate specification.  Begin is greater than end! $begin>$end\n";
		}

		if(!defined($begins{$seq_id})){
			$begins{$seq_id}=$begin;
		}else{
			if($begin>$begins{$seq_id}){
				$begins{$seq_id}=$begin;
			}
		}

		if(!defined($ends{$seq_id})){
			$ends{$seq_id}=$end;
		}else{
			if($end<$ends{$seq_id}){
				$ends{$seq_id}=$end;
			}
		}
	}

	close(IN_FILE_FH);
}

if(length(%begins)!=length(%ends)){
	die "Error, number of begins != number of ends!\n";
}

foreach my $seq_id(sort keys %begins){
	my $outstr=join "\t", ($seq_id, $begins{$seq_id}, $ends{$seq_id});
	print OUT_FH "$outstr\n";
}

print STDERR "done.\n";



