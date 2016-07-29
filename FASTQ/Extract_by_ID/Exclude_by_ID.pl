#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_f $opt_l $opt_o);

getopts("f:l:o:");
my $usage = "usage: 
$0 
	-f <fastq file>
	-l <read ID list>
	-o <output fastq file>

	This script will read in your read ID list and then
	look for the IDs in your input fastq file, excluding them 
	them as they are found. 

	If the read ID's have @ in front of them, then they will
	be stripped off.  It is assumed you won't have any read ID
	with legitimate @ in them.

";

if(!(
	defined($opt_f) ||
	defined($opt_l) ||
	defined($opt_o))){
	die $usage;
}

my $input_fname=$opt_f;
my $id_list_fname=$opt_l;
my $output_fname=$opt_o;


print STDERR "\n";
print STDERR "Input FASTAQ File: $input_fname\n";
print STDERR "ID List File: $id_list_fname\n";
print STDERR "Output Filename Root: $output_fname\n";
print STDERR "\n";

###############################################################################

print STDERR "Loading ID file.\n";

open(LIST_FH, "<$id_list_fname") || die "Could not open ID list: $id_list_fname\n";

my %id_hash;

while(<LIST_FH>){
	chomp;
	my $id=$_;
	$id=~s/^@//;
	$id_hash{$id}=1;	
}

my $num_ids=keys %id_hash;
print STDERR "Number of Targeted IDs to Exclude: $num_ids\n";

close(LIST_FH);

if(0){
	print STDERR "Records in list (cleaned):\n";
	foreach my $id(keys %id_hash){
		print STDERR "$id\n";
	}
	print STDERR "\n";
}

###############################################################################

if($input_fname=~/\.gz$/){
	open(FASTQ_FH, "zcat $input_fname| ") || die "Could not open $input_fname\n";
}else{
	open(FASTQ_FH, "<$input_fname") || die "Could not open $input_fname\n";
}

open(OUTPUT_FH, ">$output_fname") || die "Could not open $output_fname for writing.\n";

###############################################################################

print STDERR "Reading/Yanking from FASTQ file: $input_fname ...\n";

my $num_included=0;
my $rec_idx=0;

while(!eof(FASTQ_FH)){

        my $defline=<FASTQ_FH>;
        my $seq=<FASTQ_FH>;
        my $plus=<FASTQ_FH>;
        my $qv=<FASTQ_FH>;

	my $id;
	if($defline=~/^@(\S+)/){
		$id=$1;
	}else{
		die "Could not parse ID out of defline: $defline\n";
	}

	if(!defined($id_hash{$id})){
		print OUTPUT_FH "$defline";
		print OUTPUT_FH "$seq";
		print OUTPUT_FH "$plus";
		print OUTPUT_FH "$qv";
		$num_included++;
	}

        $rec_idx++;
}

close(OUTPUT_FH);

print STDERR "$num_included of $rec_idx records included.\n";
print STDERR "$rec_idx records processed.\n";

if(($rec_idx-$num_included)!= $num_ids){
	print STDERR "WARNING: NOT ALL RECORDS REQUESTED FOR EXCLUSION FOUND!\n";
}


###############################################################################

print STDERR "Done.\n\n";

