#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use Sys::Hostname;
use vars qw($opt_i $opt_o $opt_n);

getopts("i:o:n:");
my $usage = "usage: 
$0 
	-i <input fasta file>
	-o <output fasta file>
	-n \"<new name>\"

	Reads in the input fasta file and renames all the records by the new name
	that was specified.  

	If the input and output fasta file name are the same, then the input
	is replaced with the output file.

	Code isn't very interesting, but it's useful for scripting.

	Note that this script isn't very useful for multi-FASTA files
	because the same <new name> is used over and over again.

";

if(!defined($opt_i) || !defined($opt_o) || !defined($opt_n)){
	die $usage;
}

my $input_fasta=$opt_i;
my $output_fasta=$opt_o;
my $new_name=$opt_n;

###############################################################################

# Get unique name for temporary file.
my $host=hostname();
my $pid=$$;
my $temp_fname="$output_fasta\.$host\.$pid\.tmp";

print STDERR "Using $temp_fname as temporary file.\n";

# Make sure we can use this temp file name
if(-e $temp_fname){
	die "Could not open $temp_fname for writing.  It already exists.\n";
}

open(IN_FASTA, "<$input_fasta") || die "Could not open $input_fasta\n";
open(OUT_FASTA, ">$temp_fname") || die "Could not open $temp_fname\n";

print STDERR "Processing FASTA file...\n";

my ($defline, $prev_defline, $sequence);
while(<IN_FASTA>){
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

# Either replace the input file or move the temp to the new output file name
rename($temp_fname, $output_fasta);

print STDERR "Completed.\n";

###############################################################################

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	print OUT_FASTA ">$new_name\n";

	my $length=length($sequence);
	my $width=80;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print OUT_FASTA substr($sequence, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);
}

###############################################################################
