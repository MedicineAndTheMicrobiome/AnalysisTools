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

	This program will extract the coordinates you specify.
	You need to pipe the source FASTA file into STDIN.
	The output file will be sent to STDOUT.  A begin and
	end tag will be appended to the defline to help you
	differentiate that defline from the input FASTA file's
	defline.
	
	If the input file is a multi-fasta file, it will extract
	the same region from each record, which is probably not
	such a useful thing.

	If extraction end is greater than the sequence, then the 
	end of the sequence just be used.  A message will be 
	printed to STDERR.

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

	print STDOUT "$defline /begin=$opt_b /end=$opt_e\n";


	my $orig_len=length($sequence);

	my $extract_end=$opt_e;
	if($extract_end>$orig_len){
		print STDERR "Extraction range ($opt_e) exceeds sequence length ($orig_len)\n";
		$extract_end=$orig_len;
	}

	$sequence=substr($sequence, $opt_b, ($extract_end-$opt_b));

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
