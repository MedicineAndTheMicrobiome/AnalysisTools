#!/usr/bin/env perl

use strict;
use Getopt::Std;
use POSIX;
use vars qw($opt_f $opt_o $opt_s);

getopts("f:o:s");
my $usage = "
usage: 
$0 
	-f <fastq file>
	[-o <output fastq file root>]
	[-s (skip validation)]

	This script will take a FASTQ file that has forward and reverse
	reads in it and split it into two FASTQ files.

";

if(!(
	defined($opt_f))){
	die $usage;
}

my $fname=$opt_f;
my $outname=$opt_o;

if(!defined($outname)){
	$outname=$fname;
	$outname=~s/\.gz$//;
	$outname=~s/\.fastq$//;
}

my $do_validation;
if(defined($opt_s)){
	$do_validation=0;
}else{
	$do_validation=1;
}

###############################################################################

if($opt_f=~/\.gz$/){
	open(FASTQ_FH, "zcat $fname | ") || die "Could not open $opt_f\n";
}else{
	open(FASTQ_FH, "<$fname") || die "Could not open $opt_f\n";
}


my $outname_f="$outname.R1.fastq";
my $outname_r="$outname.R2.fastq";

open(R1_FH, ">$outname_f") || die "Could not open $outname_f\n";
open(R2_FH, ">$outname_r") || die "Could not open $outname_r\n";

###############################################################################

print STDERR "Input Filename: $fname\n";
print STDERR "Output Filename Root: $outname\n";
print STDERR "Do validation? $do_validation\n";

###############################################################################

print STDERR "Reading in FASTQ file: $fname ...\n";

my $num_recs=0;

my ($r1_id, $r1_seq, $r1_plus, $r1_qv);
my ($r2_id, $r2_seq, $r2_plus, $r2_qv);

while(!eof(FASTQ_FH)){

	if($num_recs%2){
		$r2_id=<FASTQ_FH>;
		$r2_seq=<FASTQ_FH>;
		$r2_plus=<FASTQ_FH>;
		$r2_qv=<FASTQ_FH>;

		if($do_validation){
			my $r1_clean_id=$r1_id;
			my $r2_clean_id=$r2_id;
			chomp $r1_clean_id;
			chomp $r2_clean_id;

			($r1_clean_id)=split /\s+/, $r1_clean_id;
			($r2_clean_id)=split /\s+/, $r2_clean_id;

			my ($r1_read_name, $r1_dir)=split "/", $r1_clean_id;
			my ($r2_read_name, $r2_dir)=split "/", $r2_clean_id;

			my $error=0;
			if($r1_read_name ne $r2_read_name){
				print STDERR "Error $r1_read_name doesn't equal $r2_read_name\n";
				$error=1;
			}
			#if($r1_dir!=1 || $r2_dir!=2){
			#	print STDERR "Error $r1_dir is not /1 or $r2_dir is not /2\n";
			#	$error=1;
			#}

			if($error){
				`mv $outname_f $outname_f\.errored.fastq`;
				`mv $outname_r $outname_r\.errored.fastq`;
				die "Input file read directions not properly paired.\n";
			}
		}

		print R1_FH "$r1_id";
		print R1_FH "$r1_seq";
		print R1_FH "$r1_plus";
		print R1_FH "$r1_qv";

		print R2_FH "$r2_id";
		print R2_FH "$r2_seq";
		print R2_FH "$r2_plus";
		print R2_FH "$r2_qv";


	}else{
		$r1_id=<FASTQ_FH>;
		$r1_seq=<FASTQ_FH>;
		$r1_plus=<FASTQ_FH>;
		$r1_qv=<FASTQ_FH>;
	}

	$num_recs++;

}
print STDERR "\n";
print STDERR "Num records split: $num_recs\n";


close(FASTQ_FH);

close(R1_FH);
close(R2_FH);

print STDERR "Completed.\n";

###############################################################################
