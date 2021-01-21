#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_r);

my $SEQ_RUN_DIR="/mnt/cmmnas02/SequencingRuns";


getopts("r:");
my $usage = "usage: 

$0 
	-r <sequencing run name> 


	This script will setup a directory for the specified sequencing run in:

	$SEQ_RUN_DIR

	This script will:

	1.) check to confirm directoy of that names does not already exists
	2.) create directory of that name
	3.) create subdirectories for:
		1.) Illumina Samplesheet tsv file
		2.) Sequencing run zip file
		3.) QC30 Run
		4.) Project Index files

";

if(!(defined($opt_r))){
	die $usage;
}

my $seq_run_name=$opt_r;

###############################################################################

print STDERR "\n";
print STDERR "User specified: $seq_run_name\n";
print STDERR "\n";

my $new_run_dir="$SEQ_RUN_DIR/$seq_run_name";

if(-e $new_run_dir){
	print STDERR "Error: $new_run_dir already exists.\n\n";
	print STDERR `ls -ld $new_run_dir`. "\n";
	die;
}


sub mkdir_and_check{
	my $dname=shift;

	print STDERR "Making new directory: $dname\n";
	mkdir $dname;
	`chmod u+wrx $dname`;
	if(!(-e $dname)){
		print STDERR "Could not create new directory: $dname.\n";
		die;
	}
}

mkdir_and_check("$new_run_dir");
mkdir_and_check("$new_run_dir/Samplesheet");
mkdir_and_check("$new_run_dir/Run");
mkdir_and_check("$new_run_dir/QC");
mkdir_and_check("$new_run_dir/QC/QV_35");
mkdir_and_check("$new_run_dir/QC/QV_35/MergeMates");
mkdir_and_check("$new_run_dir/QC/QV_35/QC");
mkdir_and_check("$new_run_dir/QC/QV_30");
mkdir_and_check("$new_run_dir/QC/QV_30/MergeMates");
mkdir_and_check("$new_run_dir/QC/QV_30/QC");
mkdir_and_check("$new_run_dir/QC/QV_25");
mkdir_and_check("$new_run_dir/QC/QV_25/MergeMates");
mkdir_and_check("$new_run_dir/QC/QV_25/QC");
mkdir_and_check("$new_run_dir/ProjectIndices");
mkdir_and_check("$new_run_dir/ProjectIndices/QV_35");
mkdir_and_check("$new_run_dir/ProjectIndices/QV_30");
mkdir_and_check("$new_run_dir/ProjectIndices/QV_25");

print STDERR "\n";
print STDERR "Setting group ownership to CMM.\n";
`chgrp -R CMM $new_run_dir`;

print STDERR "\n\n";
print STDERR "Ok. Done.\n\n";

print STDERR "Remember:\n\n";
print STDERR "1.) Place the Illumina samplesheet into: $new_run_dir/Samplesheet\n";
print STDERR "\n\tFor Example:\n\n";
print STDERR "\t\tmv $seq_run_name\.tsv $new_run_dir/Samplesheet/\n";
print STDERR "\n";
print STDERR "2.) Place the zipped Illumina sequencing run into: $new_run_dir/Run\n";
print STDERR "\n\tFor Example:\n\n";
print STDERR "\t\tmv $seq_run_name\.zip $new_run_dir/Run/\n";
print STDERR "\n";

#------------------------------------------------------------------------------
