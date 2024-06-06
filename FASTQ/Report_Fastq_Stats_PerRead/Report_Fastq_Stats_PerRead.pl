#!/usr/bin/env perl

use strict;
use Getopt::Std;
use Statistics::Descriptive;
use vars qw($opt_f $opt_o);

getopts("f:o:");
my $usage = "
usage: 
$0 
	-f <fastq file>
	[-o <qv offset, assuming offset 33>]

	This script will read through a FASTQ file and for each record,
	report the length and mean QV.
";

if(!(
	defined($opt_f))){
	die $usage;
}

my $fname=$opt_f;
my $offset=33;

if(defined($opt_o)){
	$offset=$opt_o;
}

print STDERR "QV conversion offset: $offset\n";

###############################################################################

sub calc_qv_stats{
	my $qv_str=shift;
	my $qv_len=shift;
	my $offset=shift;

	my @qv_char_arr=split "", $qv_str;

	my @qv_val_arr;
	my $qv_sum=0;
	for(my $i; $i<$qv_len; $i++){
		my $qv_char=$qv_char_arr[$i];
		my $qv_val=(ord($qv_char)-$offset);
		if($qv_val < 0 || $qv_val >100){
			die "Error: QV out of expected range ($qv_val)\n";
		}

		push @qv_val_arr, $qv_val;
		if($qv_val==0){
			$qv_val=0.01;
		}
		$qv_sum+=log($qv_val);
	}

	my $stat_obj=Statistics::Descriptive::Full->new();
	$stat_obj->add_data(@qv_val_arr);

	my $geo_mean=sprintf("%2.1f", exp($qv_sum/$qv_len));
	my $min=$stat_obj->min();
	my $max=$stat_obj->max();
	my $med=$stat_obj->median();

	return($geo_mean, $min, $max, $med);
}	

###############################################################################

if($opt_f=~/\.gz$/){
	open(FASTQ_FH, "zcat $fname | ") || die "Could not open $opt_f\n";
}else{
	open(FASTQ_FH, "<$fname") || die "Could not open $opt_f\n";
}

print STDERR "Reading in FASTQ file: $fname ...\n";

my @lengths;
my @med_qvs;
my $num_recs=0;
my $num_bases=0;

my @qv_sample;

print STDOUT "#id\tlen\tgeom_mean\tmedian\tmin\tmax\n";

while(!eof(FASTQ_FH)){

	# Read in record
	my $id=<FASTQ_FH>;
	my $seq=<FASTQ_FH>;
	my $plus=<FASTQ_FH>;
	my $qv=<FASTQ_FH>;

	chomp($id);
	chomp($seq);
	chomp($qv);

	my $seq_len=length($seq);
	my $qv_len=length($qv);

	if($seq_len ne $qv_len){
		die "Error: Sequence ($seq_len) and QVs ($qv_len) not the same length.\n";
	}else{
		my ($geo_mean, $min, $max, $med)=calc_qv_stats($qv, $seq_len, $offset);

		($id)=split "\\s+", $id;
		$id=~s/^@//;
		print STDOUT "$id\t$seq_len\t$geo_mean\t$med\t$min\t$max\n";
	}

}

close(FASTQ_FH);

###############################################################################

print STDERR "Completed.\n";

