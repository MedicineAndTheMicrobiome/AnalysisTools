#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_d $opt_x $opt_s $opt_r $opt_o);
use File::Basename;
use Cwd;

getopts("d:x:s:ro");

my $usage = "
	usage:
	$0

		-d <root data file directory>
		-x <file eXtension to run program on>
		-s <execution Script>
		[-r (rename original)]
		[-o (do not execute if .stderr and .stdout file already exist)]
		
	This script will go into the specified data directory and identify all the files
	with the specificed extension and run the analysis specified in the execution script.

	The -s option will rename the file that the execution script was called on to: <original>.orig
	
	Before running the script, the script will go into the directory, just in case some output
	files are accidentally being created without a path.  Before executing the script specified,
	the environmental variables will be set to:

		TARGET_FILE:  The file identified with the extension specified
		RESULT_FILE:  The file identified as the target, with the extension stripped off.

	The RESULT_FILE is optional if there is a default output file based on the input TARGET_FILE.

	example:
		$0 
			-d \${PROJECT}/data/summary_tables/AllSamples \\
			-x .summary_table.tsv \\
			-s run_malr.csh \\


	where run_malr.csh looks like:
		
		#!/bin/csh
		
		\${COUNT_ANALYSES_PATH}/Fit_Multivarate_ALR_Regression/Fit_Multivarate_ALR_Regression.r \\
			-s \${TARGET_FILE} \\
			-o \${RESULT_FILE} \\
			<other parameters...>

	A run_malr.csh should be written/maintained for each run in order to keep track of which
	parameters were employed in each analysis.  For example, the cutoffs and model would be
	in the <other paramters...>.

	Standard out and standard error will be written to a file next to each input data file that the program
	was run on.
	
	if the -o option is set, then the execution will be skipped if both the stderr and stdout
	file already exist.

";

if(!defined($opt_d) || !defined($opt_x) || !defined($opt_s)){
	die $usage;
}

my $RenameOriginal=0;
if(defined($opt_r)){
	$RenameOriginal=1;	
}

my $SkipRerun=0;
if(defined($opt_o)){
	$SkipRerun=1;
}

my $DataRootDirectory=$opt_d;
my $TargetExtension=$opt_x;
my $ExecutionScript=$opt_s;

###############################################################################

print STDERR "\n";
print STDERR "Data Root Directory: $DataRootDirectory\n";
print STDERR "Target Extension: $TargetExtension\n";
print STDERR "Analysis Execution Script: $ExecutionScript\n";
print STDERR "\n";

###############################################################################

print STDERR "Confirming existance of directories and scripts.\n";

if(-e $DataRootDirectory){
	print STDERR "\tFound $DataRootDirectory\n";
}else{
	die "\tCould not find $DataRootDirectory\n";
}

if(-e $ExecutionScript){
	print STDERR "\tFound $ExecutionScript\n\n";
	if($ExecutionScript=~/^\//){
		print STDERR "Execution script was absolutely specified.  Good.\n";
	}else{
		$ExecutionScript=getcwd() . "/" . $ExecutionScript;
		print STDERR "Execution script was relatively specified.  Modifying to:\n";
		print STDERR "\t$ExecutionScript\n";
	}
}else{
	die "\tCould not find $ExecutionScript\n";
}

print STDERR "\n";

###############################################################################

sub pop_file{
	my $file_path=shift;
	my @file_path_arr=split "/", $file_path;
	my $file=pop @file_path_arr;
	my $popped_file_path=join "/", @file_path_arr;
	return("$popped_file_path/", $file);
}

###############################################################################

sub get_targets{
	my $data_path=shift;
	my $target_extension=shift;

	print STDERR "Looking into $data_path for $target_extension file...\n";
	my $find_res=`find $data_path | sort`;
	my @files=split "\n", $find_res;
	
	my @targeted_files;
	foreach my $file(@files){
		if($file=~/$target_extension$/){
			push @targeted_files, $file;
		}
	}
	return(\@targeted_files);
}

###############################################################################

sub prep_and_execute{
	my $target=shift;
	my $target_ext=shift;
	my $script=shift;
	my $rename=shift;

	my ($target_path, $target_file)=pop_file($target);

	chdir $target_path || die "Could not change into $target_path\n";

	my $output=$target;
	$output=~s/$target_ext$//;

	if($SkipRerun){
		if(-e "$output.stderr" && -e "$output.stdout"){
			print STDERR "Output aleady exists, skipping...\n";
			print STDERR "Not executing: '$script'\n";
			return("Skipped.\n");
		}
	}


	$ENV{"TARGET_FILE"}=$target;
	$ENV{"RESULT_FILE"}=$output;
	my $res=`$script 2> $output.stderr`;

	open(FH, ">$output.stdout") || die "Could not open $output.stdout\n";
	print FH $res;
	close(FH);

	if($rename){
		print STDERR "Renaming $target\n";
		print STDERR "      To $target.orig\n";
		`mv $target $target.orig`
	}

	return($res);
}

###############################################################################

my $target_files_ref=get_targets($DataRootDirectory, $TargetExtension);

foreach my $target(@{$target_files_ref}){
	print STDERR "Working on: $target\n";	

	my $res=prep_and_execute($target, $TargetExtension, $ExecutionScript, $RenameOriginal);

	print "$res\n";
}

print STDERR "Done.\n";

