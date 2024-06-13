#!/usr/bin/env perl

use strict;
use Getopt::Std;
use Statistics::Descriptive;
use vars qw($opt_f $opt_h $opt_o $opt_a);

my $qv_sample_size=1000000;

getopts("f:ho:a:");
my $usage = "
usage: 
$0 
	-f <fastq file>
	[-o <qv offset, assuming offset 33>]
	[-h (no header flag)]
	[-a <alternative name>]

	This script will read through a FASTQ file and report the number of records
	and basic read statistics. It is assumed that there is only one sequence and
	one QV line per file.  So 4 lines per record.

	By default the offset will be set to 33.  If the quality values exceed 41
	an error will be thrown.  

	The following statistics will be reported:
		Name	
		NumRecords
		MedianLength
		MinLength
		MaxLength
		NumBases
		MedianQV
		LB95QV
		UB95QV	

	The QV stats, since there are so many QV to go through are based on 
	randomly sampling a QV from each read.  Then storing up to $qv_sample_size
	values from which to calculate the percentiles.

	The alternative name option (-a) may be used to substitute the filename/path
	with a user specified name.  For example, after stripping off the path
	or just some other name specified in a pipeline for that fastq file.
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
my $nohdr=defined($opt_h);

my $output_name;
if(defined($opt_a)){
	$output_name=$opt_a;
}else{
	$output_name=$fname;
}

print STDERR "QV conversion offset: $offset\n";

###############################################################################

sub sample_qv{
	my $qv_str=shift;
	my $idx=shift;

	my $qv_char=substr($qv_str, $idx, 1);
	my $qv_val=(ord($qv_char)-$offset);

	if($qv_val < 0){
		die "Error: Negative QVs ($qv_val)\n";
	}elsif($qv_val > 55){
		warn "Warning: High QVs ($qv_val)\n";
	}
	
	return($qv_val);
}	

###############################################################################
# Make sure files open before wasting any time doing anything

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

while(!eof(FASTQ_FH)){

	# Read in record
	my $id=<FASTQ_FH>;
	my $seq=<FASTQ_FH>;
	my $plus=<FASTQ_FH>;
	my $qv=<FASTQ_FH>;

	chomp($seq);
	chomp($qv);

	my $seq_len=length($seq);

	if($seq_len ne length($qv)){
		die "Error: Sequence and QVs not the same length.\n";
	}else{
		# Randomly pick position to sample QV from
		my $idx=int(rand() * $seq_len);
		my $qv=sample_qv($qv, $idx);
		#print STDERR "$idx / $seq_len \n";

		if(($#qv_sample+1)<$qv_sample_size){
			# If buffer isn't full keep load it up
			push @qv_sample, $qv;
		}else{
			# If buffer is full, randomly replace a value with a new one
			my $buf_idx=int(rand() * $qv_sample_size);
			$qv_sample[$buf_idx]=$qv;
			#print STDERR "$buf_idx\n";
		}
	
		# Store all lengths
		push @lengths, $seq_len;

		# Keep track of number of records
		$num_recs++;

		# Keep track of total residues
		$num_bases+=$seq_len;
	}

}

close(FASTQ_FH);

###############################################################################

# Compute stats on length
my $len_stat = Statistics::Descriptive::Full->new();
$len_stat->add_data(@lengths);

my $len_min=$len_stat->min();
my $len_max=$len_stat->max();
my $len_lb=$len_stat->percentile(2.5);
my $len_ub=$len_stat->percentile(97.5);
my $len_median=$len_stat->median();

if($len_lb eq "" || $len_ub eq ""){
	$len_lb="NA";
	$len_ub="NA";
}

# Compute stats on QV
print STDERR "QV Subsample size:" . ($#qv_sample+1) . "\n";
#foreach my $qv(@qv_sample){
#	print STDERR "$qv\n";
#}

my $qv_stat = Statistics::Descriptive::Full->new();
$qv_stat->add_data(@qv_sample);

my $qv_median=$qv_stat->median();
my $qv_min=$qv_stat->min();
my $qv_max=$qv_stat->max();
my $qv_lb=$qv_stat->percentile(2.5);
my $qv_ub=$qv_stat->percentile(97.5);

if($qv_lb eq "" || $qv_ub eq ""){
	$qv_lb="NA";
	$qv_ub="NA";
}

# Output statistics
if(!$nohdr){
	my $num_tabs_in_name=($output_name=~s/\t/\t/g);
	my $padding="\t" x $num_tabs_in_name;
	print STDOUT "# Name$padding\tNumRecords\tNumBases\tMedianLength\tMinLength\tMaxLength\tLB95Len\tUB95Len\tMedianQV\tMinQV\tMaxQV\tLB95QV\tUB95QV\n";
}
print STDOUT "$output_name\t$num_recs\t$num_bases\t$len_median\t$len_min\t$len_max\t$len_lb\t$len_ub\t$qv_median\t$qv_min\t$qv_max\t$qv_lb\t$qv_ub\n";

###############################################################################

print STDERR "Completed.\n";

