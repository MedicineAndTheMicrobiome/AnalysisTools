#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_f $opt_n $opt_x $opt_k $opt_r);

getopts("f:n:x:lr");
my $usage = "usage: 
$0 
	-f <Input FASTA Filename>
	-n <miN sequence length, optional>
	-x <maX sequence length, optional>
	[-k (generate an keep list, instead of fasta file)]
	[-r (generate an exclusion list, instead of fasta file)]

	Computes sequence length and outputs sequences within range to
	stdout.  If the -l option is specified, then the output file
	will be a list of sequence IDs to excluded.

";

if(!(
	defined($opt_f))){
	die $usage;
}

my $min;
my $max;
my $generate_list;

if(defined($opt_n)){
	$min=$opt_n+0;
}else{
	$min=undef;
}

if(defined($opt_x)){
	$max=$opt_x+0;
	if(uc($max) eq "INF"){
		$max=undef;
	}
}else{
	$max=undef;
}

if(defined($opt_k)){
	$generate_list=1;
}elsif(defined($opt_r)){
	$generate_list=-1;
}else{
	$generate_list=0;
}

###############################################################################
# Make sure files open before wasting any time doing anything

my $num_kept=0;
my $num_excluded=0;
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

	my $keeper=0;
	if(	(!defined($opt_n) || $length>=$min) &&
		(!defined($opt_x) || $length<=$max)){
		$keeper=1;
		$num_kept++;
	}else{
		$keeper=0;
		$num_excluded++;
	}

	$defline=~s/^>//;
	my @defsplit=split /\s+/, $defline;
	my $seq_id=$defsplit[0];

	if($generate_list==0){

	    print STDOUT "$defline\n";

	    my $width=100;
	    my $pos=0;
	    do{
		    my $out_width=($width>$length)?$length:$width;
		    print STDOUT substr($sequence, $pos, $width) . "\n";
		    $pos+=$width;
		    $length-=$width;
	    }while($length>0);
	    
	}elsif($generate_list==-1 && !$keeper){
		print STDOUT "$seq_id\n";
	}elsif($generate_list==1 && $keeper){
		print STDOUT "$seq_id\n";
	}
}

#------------------------------------------------------------------------------

