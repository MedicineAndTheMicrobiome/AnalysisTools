#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_q $opt_o $opt_r);

getopts("q:o:r:");
my $usage = "usage: 
$0 
	-q <gapped quality value file>
	-o <output average qv values file>
	[-r flag to round QVs]

	This script will read in a gapped quality value file and compute
	the average quality value for each position.  The averaging is done
	in base error space.  

	Error_per_base = exp(-LN10*qv/10);

";

if(
	!(defined($opt_q)) ||
	!(defined($opt_o))
){
	die $usage;
}

###############################################################################

my $LN10=log(10);

my $qv_file=$opt_q;
my $out_file=$opt_o;
my $round_qv=defined($opt_r);

open(QV_IN, "<$qv_file") || die "Could not open $qv_file for reading.\n";
open(OUT, ">$out_file") || die "Could not open $out_file for writing.\n";

###############################################################################

my %qv_hash;

my ($defline, $prev_defline, $sequence);
while(<QV_IN>){
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

#------------------------------------------------------------------------------

my $min=undef;
my $max=undef; 


# Check lengths of alignments
foreach my $id(keys %qv_hash){
	my $len=$#{$qv_hash{$id}}+1;
	if(!defined($min) || $len<$min){
		$min=$len;
	}
	if(!defined($max) || $len>$max){
		$max=$len;
	}
}	

my $len;
if($min!=$max){
	die "The min and max are not the same for this set of QV.\n";
}else{
	$len=$max;
	print "Alignment length is: $len\n";
}

# Initialize variables
my @error_rate;
my @coverage;
my $i;

for($i=0; $i<=$len; $i++){
	$error_rate[$i]=0;
	$coverage[$i]=0;
}

# Sum up error rates and coverage
foreach my $id(keys %qv_hash){
	for($i=0; $i<$len; $i++){
		my $qv=${$qv_hash{$id}}[$i];
		if($qv ne "-"){
			$error_rate[$i]+=qv_to_error_rate($qv);	
			$coverage[$i]++;
		}
	}
}

# Perform averaging in error rate, and then convert to qv.
my @average_qv;
for($i=0; $i<$len; $i++){
	if($coverage[$i]==0){
		die "Error coverage at position $i is zero.\n";
	}
	$error_rate[$i]/=$coverage[$i];
	$average_qv[$i]=error_rate_to_qv($error_rate[$i]);
	#print "$error_rate[$i] => $average_qv[$i]\n";
}

# Output QV, round if requested.
for($i=0; $i<$len; $i++){
	if($round_qv){
		$average_qv[$i]=int($average_qv[$i]+0.5);
	}else{
		$average_qv[$i]=sprintf("%5.4f", $average_qv[$i]);
	}
	print OUT "$average_qv[$i]\n"; 
}

print STDERR "Completed.\n";

###############################################################################

# QV/Error rate formalas
sub qv_to_error_rate{
	my $qv=shift;
	return(exp(-$LN10*$qv/10));
}

sub error_rate_to_qv{
	my $err_rate=shift;
	return(-10*log($err_rate)/$LN10);
}

# Load quality values
sub process_record{
	my $defline = shift;
	my $sequence = shift;

	my $id;
	if($defline=~/^(\S+)/){
		$id=$1;	
	}else{
		die "Could not parse out id from defline: $defline\n";
	}

	$sequence=~s/^\s+//;
	my @qv=split /\s+/, $sequence;

	$qv_hash{$id}=\@qv;
}

#------------------------------------------------------------------------------
