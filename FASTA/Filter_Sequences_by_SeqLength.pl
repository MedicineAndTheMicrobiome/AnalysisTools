#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_f $opt_n $opt_x);

getopts("f:n:x:");
my $usage = "usage: 
$0 
	-f <Input FASTA Filename>
	-n <miN sequence length, optional>
	-x <maX sequence length, optional>

	Computes sequence length and outputs sequences within range to
	stdout.

";

if(!(
	defined($opt_f))){
	die $usage;
}

my $min;
my $max;

if(defined($opt_n)){
	$min=$opt_n;
}

if(defined($opt_x)){
	$max=$opt_x;
}

###############################################################################
# Make sure files open before wasting any time doing anything

my $num_kept=0;
my $total=0;

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

my $perc=sprintf("%3.2f", 100*($num_kept/$total));
print STDERR "$num_kept ($perc%) kept out of $total total \n";
print STDERR "Completed.\n";

###############################################################################


sub process_record{
	my $defline = shift;
	my $sequence = shift;

	$total++;
	my $length=length($sequence);

	if(	(!defined($opt_n) || $length>=$min) &&
		(!defined($opt_x) || $length<=$max)){

		    print STDOUT "$defline\n";

		    my $width=60;
		    my $pos=0;
		    do{
			    my $out_width=($width>$length)?$length:$width;
			    print STDOUT substr($sequence, $pos, $width) . "\n";
			    $pos+=$width;
			    $length-=$width;
		    }while($length>0);
		    
		    $num_kept++;
	}
}

#------------------------------------------------------------------------------

