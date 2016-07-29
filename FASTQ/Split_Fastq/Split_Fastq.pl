#!/usr/bin/env perl

use strict;
use Getopt::Std;
use vars qw($opt_f $opt_a $opt_q $opt_d $opt_l $opt_o);

my $OFFSET=33;

getopts("f:a:q:d:l:o:");
my $usage = "
usage: 
$0 
	-f <fastq file>
	[-l <min length>]
	[-a <output fasta file>]
	[-q <output quality file>]
	[-d <decoded quality file>]
	[-o <offset for decoding quality file, $OFFSET>]

	This script will read in a fastq file and generate
	fasta and quality file. The quality file will not be
	decoded into values unless you use the -d option. 

	The -l flag is the minimum length you want to extract.

";

if(!(
	defined($opt_f))){
	die $usage;
}

my $fname=$opt_f;
my $qual_name=$opt_q;
my $fasta_name=$opt_a;
my $decoded_qual_name=$opt_d;
my $min_length=$opt_l;

if(!defined($min_length)){
	$min_length=0;
}

my $offset;
if(defined($opt_o)){
	$offset=$opt_o;
}else{
	$offset=$OFFSET;
}

###############################################################################

if($opt_f=~/\.gz$/){
	open(FASTQ_FH, "zcat $fname | ") || die "Could not open $opt_f\n";
}else{
	open(FASTQ_FH, "<$fname") || die "Could not open $opt_f\n";
}

my $prod_qual=defined($qual_name);
my $prod_fasta=defined($fasta_name);
my $prod_decoded_qual=defined($decoded_qual_name);

if($prod_fasta){
	open(FASTA_FH, ">$fasta_name") || die "Could not open $fasta_name\n";
}

if($prod_qual){
	open(QUAL_FH, ">$qual_name") || die "Could not open $qual_name\n";
}

if($prod_decoded_qual){
	open(DECODED_QUAL_FH, ">$decoded_qual_name") || die "Could not open $decoded_qual_name\n";
}

###############################################################################

sub decode_qv_to_str{
        my $coded_qv_str=shift;
	my $offset=shift;

	my @qv_chars=split //, $coded_qv_str;
	my @qv_vals;

	for(my $i=0; $i<=$#qv_chars; $i++){
		my $qv_val=(ord($qv_chars[$i])-$offset);	

		if($qv_val < 0 || $qv_val >41){
			die "Error: QV out of expected range ($qv_val).  You may need to change the offset: $offset.\n";
		}
		push @qv_vals, $qv_val;
	}

	my $qv_val_str=join " ", @qv_vals;

        return($qv_val_str);
}

###############################################################################

print STDERR "Reading in FASTQ file: $fname ...\n";

my $num_recs=0;

while(!eof(FASTQ_FH)){

	my $id=<FASTQ_FH>;
	my $seq=<FASTQ_FH>;
	my $plus=<FASTQ_FH>;
	my $qv=<FASTQ_FH>;

	my $defline=$id;
	$defline=~s/^@/>/;	

	my $seq_len=length($seq)-1;

	if($seq_len>=$min_length){

		if($prod_qual){
			print QUAL_FH $defline;
			print QUAL_FH $qv;	
		}

		if($prod_fasta){
			print FASTA_FH $defline;
			print FASTA_FH $seq;
		}

		if($prod_decoded_qual){
			chomp $qv;
			my $d_qv=decode_qv_to_str($qv, $offset);
			print DECODED_QUAL_FH $defline;
			print DECODED_QUAL_FH "$d_qv\n";
		}

	}
}

close(FASTQ_FH);

if($prod_qual){
	close(QUAL_FH);
}
if($prod_fasta){
	close(FASTA_FH);
}

print STDERR "Completed.\n";

###############################################################################
