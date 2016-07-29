#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_w $opt_o);

getopts("w:o:");
my $usage = "usage: 

$0 

    This program reads in from STDIN and outputs to STDOUT.

    -w <Width [default 3] zero padded width of the number you want to append>
    -o <Offset [default 0] offset from where to start counting>

	This program will append a .<offset> to the end of the key of the defline entry.
	This is useful if you have defline keys that are not unique.  

	For example, if you set the options -o 0 and -w 3 (the defaults):

	    >key
	    CATGAGTAG
	    >anotherkey /tag=xxx
	    AATGGACCATG

	will become

	    >key.000
	    CATGAGTAG
	    >anotherkey.001 /tag=xxx
	    AATGGACCATG

";

print STDERR $usage;

###############################################################################

print STDERR "Processing FASTA file...\n";

my $record_count;
my $field_width;

if(defined($opt_o)){
	$record_count=$opt_o;
}else{
	$record_count=0;
}

if(defined($opt_w)){
	$field_width=$opt_w;
}else{
	$field_width=3;
}

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

	my @fields=split / /, $defline;
	$fields[0] .= sprintf ".%0" . $field_width . "i", $record_count;
	$defline=join " ", @fields;

	print STDOUT "$defline\n";
	my $length=length($sequence);
	my $width=50;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print STDOUT substr($sequence, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);

	$record_count++;
}

#------------------------------------------------------------------------------
