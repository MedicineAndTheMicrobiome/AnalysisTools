#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_f $opt_e);

getopts("f:e:");
my $usage = "usage: 
$0 
	-f <Input FASTA Filename>
	-e \"expression\"

	Tries to act like grep.  So when the expression is found in the defline,
	then the entire sequence is reported, not just the defline.

";

if(!(
	defined($opt_f) && 
	defined($opt_e))){
	die $usage;
}


my $expression=$opt_e;

###############################################################################
# Read in features

my $num_found=0;

open(FASTA_FH, "<$opt_f") || die "Could not open $opt_f\n";

print STDERR "Processing FASTA file...\n";
my ($defline, $prev_defline, $sequence);
while(<FASTA_FH>){
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

close(FASTA_FH);

print STDERR "$num_found found\n";
print STDERR "Completed.\n";

###############################################################################


sub process_record{
	my $defline = shift;
	my $sequence = shift;

	if($defline=~/$expression/){
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
		
		$num_found++;
	}
}

#------------------------------------------------------------------------------

