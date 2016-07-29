#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_f $opt_n $opt_o);

getopts("f:n:o:");
my $usage = "usage: 
$0 
	-f <Input FASTA Filename>
	-o <Output Root FASTA Filename>
	-n <number of pieces to split into>

	This program will split the input FASTA file into
	the number of files specified, using the specified
	output root.  There is no balancing of file sizing.
	The placement of records is round robin so that
	consecutive records in the input FASTA file are 
	unlikely to be placed in the same output file.
";

if(!(
	defined($opt_f) && 
	defined($opt_o) && 
	defined($opt_n))){
	die $usage;
}

my $num_files=$opt_n;
my $output_root=$opt_o;

print STDERR "Num Files to break into: $num_files\n";
print STDERR "Output root: $output_root\n";

###############################################################################
# Make sure files open before wasting any time doing anything

open(FASTA_FH, "<$opt_f") || die "Could not open $opt_f\n";

my $i;
my @fhandle_arr;
for($i=0; $i<$num_files; $i++){	
	my $num_ext=sprintf("%03i", $i);
	$fhandle_arr[$i]=new FileHandle ">$output_root\.$num_ext";
}

###############################################################################
# Read in features

print STDERR "Processing FASTA file...\n";
my ($defline, $prev_defline);
my @sequence;
while(<FASTA_FH>){
	chomp;
	
	if(/^>/){
		$defline=$_;
		if($#sequence!=-1){
			process_record($prev_defline, \@sequence);
			@sequence=();
		}
		$prev_defline=$defline;
	}else{
		push @sequence, (split //, $_);
	}
}
process_record($prev_defline, \@sequence);

close(FASTA_FH);

print STDERR "Completed.\n";

###############################################################################

my $record_count=0;

sub process_record{
	my $defline = shift;
	my $sequence_ref = shift;
	
	my $i=$record_count % $num_files;
	
	print {$fhandle_arr[$i]} "$defline\n";

	my $length=$#sequence+1;
	my $width=80;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;

		my @out_arr=splice @{$sequence_ref}, 0, $out_width;
		print {$fhandle_arr[$i]} (join "", @out_arr) . "\n";

		$length-=$width;
	}while($length>0);

	$record_count++;
}

#------------------------------------------------------------------------------

