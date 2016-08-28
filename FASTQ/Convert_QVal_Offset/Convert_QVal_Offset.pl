#!/usr/bin/env perl

use strict;
use Getopt::Std;
use vars qw($opt_f $opt_t $opt_o);

my $OFFSET=33;

getopts("f:t:o:");
my $usage = "
usage: 
$0 
	-f <fastq file>
	[-t <target offset, def = 33>]
	-o <converted fastq file>

	This script will read in a fastq file and convert
	the quality values to the offset-33 encoding.

	The offset-33 encoding is used in (PREFERRED):
		Sanger and Illumina 1.8+

	Characters used: !\"#\$%&'()*+,-./0123456789:;<=>?\@ABCDEFGHI


	The offset-64 encoding is used in:
		Solexa, Illumina 1.3+, and Illumina 1.5+

	Characters used: \@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefg

";

if(!(
	defined($opt_f) &&
	defined($opt_o))){
	die $usage;
}

my $in_fastq=$opt_f;
my $target_offset=33;
my $out_fastq=$opt_o;
my $source_offset;

if(defined($opt_t)){
	$target_offset=$opt_t;
}

if($target_offset == 64){
	$source_offset = 33;
}elsif($target_offset == 33){
	$source_offset = 64;
}else{
	die "Target offset is non-standard.  i.e. not 64 or 33.\n";
}

###############################################################################

print STDERR "Input FASTQ: $in_fastq\n";
print STDERR "Target Offset: $target_offset\n";
print STDERR "Source Offset: $source_offset (Assumed)\n";
print STDERR "Output FASTQ: $out_fastq\n";

###############################################################################

print STDERR "Opening FASTQ file: $in_fastq...\n";
if($in_fastq=~/\.gz$/){
	open(FASTQ_FH, "zcat $in_fastq| ") || die "Could not open $in_fastq\n";
}else{
	open(FASTQ_FH, "<$in_fastq") || die "Could not open $in_fastq\n";
}

open(OUT_FQ_FH, ">$out_fastq") || die "Could not open $out_fastq\n";

###############################################################################

sub decode_qv_to_arr{

	# This function convert a coded QV string into an array of QV values.
	# If the values are out of range, an error will be returned.

        my $coded_qv_str=shift;
	my $offset=shift;

	my @qv_chars=split //, $coded_qv_str;
	my @qv_vals;
	my $error_n=0;
	my $error_p=0;

	for(my $i=0; $i<=$#qv_chars; $i++){
		my $qv_val=(ord($qv_chars[$i])-$offset);	

		if($qv_val < 0){
			$error_n--;
		}elsif($qv_val >41){
			$error_p++;
		}

		push @qv_vals, $qv_val;
	}

        return(\@qv_vals, $error_n, $error_p);
}

sub recode_arr_to_qv{

	# This function will convert an array of QV values into a
	# coded QV string based on the specified offset.

	my $qv_arr_ref=shift;
	my $offset=shift;

	my @char_arr;	

	for(my $i=0; $i<=$#{$qv_arr_ref}; $i++){
		push @char_arr, chr(${$qv_arr_ref}[$i]+$offset);
	}

	return(join "", @char_arr);
}

###############################################################################

my $num_recs=0;
my $tot_errn=0;
my $tot_errp=0;

print STDERR "Reading FASTQ file...\n";

while(!eof(FASTQ_FH)){

	my $id=<FASTQ_FH>;
	my $seq=<FASTQ_FH>;
	my $plus=<FASTQ_FH>;
	my $coded_qv_str=<FASTQ_FH>;

	chomp $coded_qv_str;

	my ($qv_val_arr_ref, $errn, $errp)=decode_qv_to_arr($coded_qv_str, $source_offset);

	$tot_errn+=$errn;
	$tot_errp+=$errp;
	
	my $recoded_qv=recode_arr_to_qv($qv_val_arr_ref, $target_offset);

	print OUT_FQ_FH "$id";
	print OUT_FQ_FH "$seq";
	print OUT_FQ_FH "$plus";
	print OUT_FQ_FH "$recoded_qv\n";

}

close(FASTQ_FH);
close(OUT_FQ_FH);

print STDERR "$tot_errn / $tot_errp\n";

if($tot_errn == 0 && $tot_errp == 0){
	print STDERR "QVs in range.\n";
}elsif($tot_errn < 0 && $tot_errp > 0){
	print STDERR "QVs inconsistently out of range.  Cannot confirm offset.\n";
	unlink $out_fastq;
}elsif($tot_errn<0){
	print STDERR "QVs out of range, i.e. negative.  Assumed starting offset was 64, but was 33?\n";
	print STDERR "No conversion performed.\n";
	system("cp $in_fastq $out_fastq");
}elsif($tot_errp>0){
	print STDERR "QVs out of range, i.e. too large.  Assumed started offset was 33, but was 64?\n";
	print STDERR "No conversion performed.\n";
	system("cp $in_fastq $out_fastq");
}

print STDERR "Completed.\n\n";

###############################################################################
