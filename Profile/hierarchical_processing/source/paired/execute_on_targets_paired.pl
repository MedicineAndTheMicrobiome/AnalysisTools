#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_d $opt_x $opt_s $opt_p $opt_j);
use File::Basename;
use Cwd;

getopts("d:x:s:o:pj:");

my $JOINER="_vs_";

my $usage = "
	usage:
	$0

		-d <root data file directory>
		-x <file eXtension to run program on>
		-s <execution Script>
		[-p (prepend shared prefix to RESULT_FILE)]
		[-j <string used to join names in RESULT_FILE, default = \"$JOINER\">]
		
	This script will go into the specified data directory and for each directory
	identify all the files that match the specified extension and run the specified
	script for each pairwise combination.

	The following environmental variables will be set before the specified
	execution script is initiated:

		TARGET_FILE_1:  The file identified with the extension specified
		TARGET_FILE_2:  The file identified with the extension specified
		RESULT_FILE:  The file identified as the target, with the extension stripped off.

	The RESULT_FILE is optional if there is a default output file based on the input TARGET_FILE.

	The script will not run if there is only one file identified in the directory.

	example:
		$0 
			-d \${PROJECT}/data/summary_tables/AllSamples \\
			-x .summary_table.tsv \\
			-s run_compare.csh \\

	where run_compare.csh looks like:
		
		#!/bin/csh
		
		\${COUNT_ANALYSES_PATH}/Compare_Microbiomes.r \\
			-a \${TARGET_FILE_1} \\
			-b \${TARGET_FILE_2} \\
			-o \${RESULT_FILE} \\
			<other parameters...>

	A run_compare.csh should be written/maintained for each run in order to keep track of which
	parameters were employed in each analysis.  For example, the cutoffs and model would be
	in the <other paramters...>.

	Standard out and standard error will be written to a file next to each input data file that the program
	was run on.

";

if(!defined($opt_d) || !defined($opt_x) || !defined($opt_s)){
	die $usage;
}

my $DataRootDirectory=$opt_d;
my $TargetExtension=$opt_x;
my $ExecutionScript=$opt_s;
my $PrependShared=defined($opt_p);
my $Joiner=defined($opt_j)?$opt_j:$JOINER;

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

sub get_target_directories{
	my $data_path=shift;
	my $target_extension=shift;

	print STDERR "Looking into $data_path for $target_extension file...\n";
	my $find_res=`find $data_path -type d | sort`;
	my @directory_arr=split "\n", $find_res;

	# For each directory, keep track of the ones with 2 or more targets
	my @targetable_dir;
	foreach my $directory(@directory_arr){
		my @ls_res=split "\n", `ls -1 $directory/*$target_extension`;
		my $num_target_files=($#ls_res+1);
		if($num_target_files >= 2){
			print STDERR "$directory: $num_target_files files\n";
			push @targetable_dir, $directory;	
		}
	}

	return(\@targetable_dir);
}

###############################################################################

sub get_targeted_files_from_directory{
	my $directory=shift;
	my $target_extension=shift;

	my @targets=split "\n", `ls $directory/*$target_extension`;
	return(\@targets);
}

###############################################################################

sub remove_common{
	my $str1=shift;
	my $str2=shift;

	#print STDERR "$str1 / $str2\n";

	my @arr1=split "\\.", $str1;
	my @arr2=split "\\.", $str2;

	# Strip common from end
	my $equal=1;
	do{
		my $c1=pop @arr1;
		my $c2=pop @arr2;

		if($c1 ne $c2){
			push @arr1, $c1;
			push @arr2, $c2;
			$equal=0;
		}
	}while($equal==1);

	# Strip common from beginning
	$equal=1;
	do{
		my $c1=shift @arr1;
		my $c2=shift @arr2;

		if($c1 ne $c2){
			unshift @arr1, $c1;
			unshift @arr2, $c2;
			$equal=0;
		}
	}while($equal==1);

	# Join remaining differences
	$str1=join ".", @arr1;
	$str2=join ".", @arr2;

	#print STDERR "$str1 / $str2\n";
	
	return($str1, $str2);
}

sub find_common{
	my $str1=shift;
	my $str2=shift;

	my @arr1=split "\\.", $str1;
        my @arr2=split "\\.", $str2;

	my @shared;

        my $equal=1;
        do{
                my $c1=shift @arr1;
                my $c2=shift @arr2;

                if($c1 eq $c2){
			push @shared, $c1;
                }else{
			$equal=0;
		}	
        }while($equal);

	my $common_str=join ".", @shared;
	return($common_str);

}

###############################################################################

sub prep_and_execute{
	my $target_1=shift;
	my $target_2=shift;
	my $target_ext=shift;
	my $script=shift;

	my ($target_path, $target_file_1)=pop_file($target_1);
	my ($target_path, $target_file_2)=pop_file($target_2);

	chdir $target_path || die "Could not change into $target_path\n";

	my ($clean_1, $clean_2)=remove_common($target_file_1, $target_file_2);
	
	my $prefix="";
	if($PrependShared){
		$prefix=find_common($target_file_1, $target_file_2);
		# If there isn't a . at the end of the prefix, add one
		if(!($prefix=~/\\.$/)){
			$prefix.=".";	
		}
	}

	my $output= $prefix . $clean_1 . $Joiner . $clean_2;

	$ENV{"TARGET_FILE_1"}=$target_1;
	$ENV{"TARGET_FILE_2"}=$target_2;
	$ENV{"RESULT_FILE"}=$output;

	print STDERR "Executing to generate: $output\n";
	my $res=`$script 2> $output.stderr`;

	open(FH, ">$output.stdout") || die "Could not open $output.stdout\n";
	print FH $res;
	close(FH);

	return($res);
}

###############################################################################

my $target_directories_ref=get_target_directories($DataRootDirectory, $TargetExtension);
print STDERR "\n";

foreach my $target_dir(@{$target_directories_ref}){
	print STDERR "Working on files in:\n\t$target_dir\n";	

	my $file_targets_arr_ref=get_targeted_files_from_directory($target_dir, $TargetExtension);
	my $num_targets=$#{$file_targets_arr_ref}+1;
	for(my $i=0; $i<$num_targets; $i++){
		print STDERR "${$file_targets_arr_ref}[$i]\n";
	}

	for(my $i=0; $i<$num_targets; $i++){
		for(my $j=($i+1); $j<$num_targets; $j++){
			prep_and_execute(
				${$file_targets_arr_ref}[$i],
				${$file_targets_arr_ref}[$j],
				$TargetExtension,
				$ExecutionScript
			);
		}
	}
}

print STDERR "Done.\n";

