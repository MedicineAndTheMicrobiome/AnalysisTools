#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_m);

my $usage = "usage: 
$0 

	This program will read in a FASTA file through STDIN and
	output the new FASTA where the homopolymers have been compressed 
	throught STDOUT.
	
";

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

sub compress{
	my $sequence=uc(shift);

	while($sequence=~s/A{2,}?/A/g){};
	while($sequence=~s/T{2,}?/T/g){};
	while($sequence=~s/G{2,}?/G/g){};
	while($sequence=~s/C{2,}?/C/g){};

	return $sequence;
}

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	print STDOUT "$defline\n";

	$sequence=compress($sequence);
	
	my $length=length($sequence);
	my $width=80;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print STDOUT substr($sequence, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);
}

#------------------------------------------------------------------------------
