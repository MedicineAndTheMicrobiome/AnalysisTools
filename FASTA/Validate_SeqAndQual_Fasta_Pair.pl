#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_s $opt_q);

getopts("s:q:");
my $usage = "usage: 
$0 
	-s <sequence fasta file>
	-q <quality fasta file>

	This program will verify that your sequence and quality files are compatible.
";

if(!(defined($opt_s)) || !(defined($opt_q))){
	die $usage;
}

###############################################################################

my $sequence_fn=$opt_s;
my $quality_fn=$opt_q;

my @sequence_id_arr;
my %sequence_length_hash;

my @quality_id_arr;
my %quality_length_hash;

print STDERR "Processing sequence file...\n";

open(SEQ_FH, "<$sequence_fn") || die "Could not open $sequence_fn\n";

my ($defline, $prev_defline, $sequence);
while(<SEQ_FH>){
	chomp;
	
	if(/^>/){
		$defline=$_;
		if($sequence ne ""){
			process_sequence_record($prev_defline, $sequence);
			$sequence="";
		}
		$prev_defline=$defline;
	}else{
		$sequence.=$_;
	}
}
process_sequence_record($prev_defline, $sequence);

print STDERR "Completed.\n";

###############################################################################

print STDERR "Processing quality file...\n";

open(SEQ_FH, "<$quality_fn") || die "Could not open $quality_fn\n";

my ($defline, $prev_defline, $sequence);
while(<SEQ_FH>){
	chomp;
	
	if(/^>/){
		$defline=$_;
		if($sequence ne ""){
			process_quality_record($prev_defline, $sequence);
			$sequence="";
		}
		$prev_defline=$defline;
	}else{
		$sequence.=$_ . " ";
	}
}
process_quality_record($prev_defline, $sequence);

print STDERR "Completed.\n";

###############################################################################

my @sequence_id_arr;
my %sequence_length_hash;

my @quality_id_arr;
my %quality_length_hash;

my $num_seq=$#sequence_id_arr+1;
my $num_qual=$#quality_id_arr+1;

if($num_seq != $num_qual){
	print "Sequence counts don't match: $num_seq != $num_qual\n";
}else{
	print "Sequences counts are good: $num_seq == $num_qual\n";
}

my $num_matching_lengths=0;
my $same_order=1;
for(my $i=0; $i<$num_seq; $i++){
	if($sequence_id_arr[$i] ne $quality_id_arr[$i]){
		print "Sequence and Quality files do not line up.\n";
		print "\tsequence record $i: $sequence_id_arr[$i] not equal $quality_id_arr[$i]\n";
		$same_order=0;
	}
	if($sequence_length_hash{$sequence_id_arr[$i]} != $quality_length_hash{$sequence_id_arr[$i]}){
		print "Sequence and Quality lengths do not match up.\n";
		print "\tsequence record $i: $sequence_id_arr[$i] \n";
		print "\tSeqLen=$sequence_length_hash{$sequence_id_arr[$i]} QualLen=$quality_length_hash{$sequence_id_arr[$i]} \n";
	}else{
		$num_matching_lengths++;
	}
}

print "Matches lengths = $num_matching_lengths of $num_seq\n";
if($same_order==1){
	print "Sequences files are in the same order.\n";
}


###############################################################################

sub process_sequence_record{
	my $defline = shift;
	my $sequence = shift;

	my $id;
	if($defline=~/^>(\S+)/){
		$id=$1;
	}else{
		die "Could not parse defline: $defline\n";
	}
	my $length=length($sequence);

	$sequence_length_hash{$id}=$length;
	push @sequence_id_arr, $id;
}

###############################################################################

sub process_quality_record{
	my $defline = shift;
	my $sequence = shift;

	my $id;
	if($defline=~/^>(\S+)/){
		$id=$1;
	}else{
		die "Could not parse defline: $defline\n";
	}

	$sequence=~s/^(\s+)//;
	my @qv=split /\s+/, $sequence;

	#print (join ",", @qv) . "\n";
	$quality_length_hash{$id}=$#qv + 1;
	push @quality_id_arr, $id;

}

#------------------------------------------------------------------------------
