#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_w);

getopts("w:");
my $usage = "usage: 
$0 
	-w Width

	This program will read in a quality FASTA file through STDIN and then
	output the same quality FASTA file through STDOUT, except with the 
	newly specified quality values per line.
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
		$sequence.=($_ . " ");
	}
}
process_record($prev_defline, $sequence);

print STDERR "Completed.\n";

###############################################################################

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	print STDOUT "$defline\n";

	$sequence=~s/^\s+//;
	my @qv=split /\s+/, $sequence;
	my $num_qv=$#qv;
	my $width=$opt_w;

	my $abs_pos=0;
	while($abs_pos<$num_qv){
	
		for(my $i=0; $i<$width && $abs_pos<$num_qv; $i++){
			if($i>0){
				print STDOUT " $qv[$abs_pos]";
			}else{
				print STDOUT "$qv[$abs_pos]";
			}
			$abs_pos++;
		}		
		print STDOUT "\n";
	}
}

#------------------------------------------------------------------------------
