#!/usr/bin/env perl

use strict;
use Getopt::Std;
use vars qw($opt_f $opt_m);

my $MAX_QV=40;
my $OFFSET=33;

getopts("f:m:");
my $usage = "
usage: 
$0 
	-f <fastq file>
	[-m <maximum quality value, def = $MAX_QV>]

	This script will read in a fastq file and limit
	the maximum quality value to the value set by the -m option.  
	The expected quality values should be coded to characters, assuming an offset of $OFFSET. 

	The output will be to STDOUT in case you want to pipe it through compression.

";

if(!(
	defined($opt_f))){
	die $usage;
}

my $in_fastq=$opt_f;
my $max_qv;

if(defined($opt_m)){
	$max_qv=$opt_m;
}else{
	$max_qv=$MAX_QV;
}
$max_qv=$max_qv+0;

###############################################################################

print STDERR "Input FASTQ: $in_fastq\n";
print STDERR "Maximum QV: $max_qv\n";
print STDERR "Source Offset: $OFFSET (Assumed)\n";

###############################################################################

print STDERR "Opening FASTQ file: $in_fastq...\n";
if($in_fastq=~/\.gz$/){
	open(FASTQ_FH, "zcat $in_fastq| ") || die "Could not open $in_fastq\n";
}else{
	open(FASTQ_FH, "<$in_fastq") || die "Could not open $in_fastq\n";
}

###############################################################################

sub decode_qv_to_arr{

	# This function convert a coded QV string into an array of QV values.
	# If the values are out of range, an error will be returned.

        my $coded_qv_str=shift;
	my $offset=shift;

	my @qv_chars=split //, $coded_qv_str;
	my @qv_vals;

	for(my $i=0; $i<=$#qv_chars; $i++){
		my $qv_val=(ord($qv_chars[$i])-$offset);	
		push @qv_vals, $qv_val;
	}

        return(\@qv_vals);
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
my $tot_abv=0;
my $tot_bel=0;
my $tot_residues=0;

print STDERR "Reading FASTQ file...\n";

while(!eof(FASTQ_FH)){

	my $id=<FASTQ_FH>;
	my $seq=<FASTQ_FH>;
	my $plus=<FASTQ_FH>;
	my $coded_qv_str=<FASTQ_FH>;

	chomp $coded_qv_str;

	my $qv_val_arr_ref=decode_qv_to_arr($coded_qv_str, $OFFSET);

	# Limit QV
	my @values=@{$qv_val_arr_ref};
	for(my $i=0; $i <= $#values; $i++){
		if($values[$i] > $max_qv){
			$values[$i] = $max_qv;
			$tot_abv++;
		}elsif($values[$i] < 0){
			$tot_bel++;
		}
	}

	$num_recs++;
	$tot_residues+=(($#values)+1);
	
	my $recoded_qv=recode_arr_to_qv(\@values, $OFFSET);

	print STDOUT "$id";
	print STDOUT "$seq";
	print STDOUT "$plus";
	print STDOUT "$recoded_qv\n";

}

close(FASTQ_FH);

print STDERR "\n";
print STDERR "Total Records: $num_recs\n";
print STDERR "Total Residues: $tot_residues\n"; 
print STDERR "QV below 0: $tot_bel " . sprintf("%6.2f", 100 * $tot_bel/$tot_residues) . "%\n";
print STDERR "QV above $max_qv: $tot_abv " . sprintf("%6.2f", 100 * $tot_abv/$tot_residues) . "%\n";

print STDERR "Completed.\n\n";

###############################################################################
