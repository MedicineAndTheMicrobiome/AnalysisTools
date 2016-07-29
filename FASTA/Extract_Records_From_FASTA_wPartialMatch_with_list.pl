#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_f $opt_l $opt_i $opt_p);

getopts("f:l:ip");
my $usage = "usage: 
$0 
	-f <Input FASTA Filename>
	-l <List of tokens to find in defline>
	[-i (ignore case)]
	[-p (periods are periods, else they are wildcards)]

	Goes throught the FASTA file and reports all sequences
	which have defline that has any of the tokens specified
	in the list.

";

if(!(
	defined($opt_f) && 
	defined($opt_l))){
	die $usage;
}

my $ignore_case=0;
my $periods_are_periods=0;

if(defined($opt_i)){
	$ignore_case=1;
}

if(defined($opt_p)){
	$periods_are_periods=1;
}

if($ignore_case){
	print STDERR "Ignoring case.\n";
}else{
	print STDERR "Case sensitive search.\n";
}

if($periods_are_periods){
	print STDERR "Periods are escape as periods.\n";
}else{
	print STDERR "Periods are treated as wild cards.\n";
}

###############################################################################

# Read in list of tokens
my @list_of_tokens;
open(LIST_FH, "<$opt_l") || die "Could not open $opt_l\n";
while(<LIST_FH>){
	chomp;

	if($ignore_case){
		$_=lc($_);
	}

	if($periods_are_periods){
		$_=~s/\./\\./g;
	}

	#print STDERR "$_\n";

	push @list_of_tokens, $_;
}
close(LIST_FH);
my $num_tokens=$#list_of_tokens+1;

###############################################################################

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

	my $report=0;

	my $search_defline=$defline;
	if($ignore_case){
		$search_defline=lc($search_defline);
	}

	for(my $i=0; $i<$num_tokens; $i++){
		my $search_term=$list_of_tokens[$i];
		if($search_defline=~/$search_term/){
			$report=1;
			last;
		}
	}

	if($report){
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

