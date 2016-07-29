#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_w);

getopts("w:");
my $usage = "usage: 
$0 
	-w Width

	This program will read in a FASTA file through STDIN and then
	output the same FASTA file through STDOUT, except with the 
	newly specified width.
";

if(!(defined($opt_w))){
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
	my $length=length($sequence);
	my $width=$opt_w;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print STDOUT substr($sequence, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);
}

#------------------------------------------------------------------------------
