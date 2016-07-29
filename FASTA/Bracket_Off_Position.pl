#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_b $opt_e);

getopts("b:e:");
my $usage = "usage: 
$0 
	-b <begin coordinate>
	-e <end coordinate>

	This program will read in a FASTA file through STDIN and then place
	square brackets around the sequence specified by the begin and end coordinates.
	The begin/end coordinates should be in 0-space based coordinates.
";

if(!(defined($opt_b)) || !defined($opt_e)){
	die $usage;
}

###############################################################################

print STDERR "Processing FASTA file...\n";

my ($defline, $prev_defline, $sequence);
while(<STDIN>){
	chomp;
	
	if(/^>/){
		$defline=$_;
		if($sequence ne ""){
			process_record($prev_defline, $sequence);
			$sequence="";
		}
		$prev_defline=$defline;
	}else{
		$sequence.=$_;
	}
}
process_record($prev_defline, $sequence);

print STDERR "Completed.\n";

###############################################################################

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	print STDOUT "$defline\n";

	substr($sequence, $opt_b,0)="[";
	substr($sequence, $opt_e+1,0)="]";

	my $length=length($sequence);
	my $width=50;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print STDOUT substr($sequence, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);
}

#------------------------------------------------------------------------------
