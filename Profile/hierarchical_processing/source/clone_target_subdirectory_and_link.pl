#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_d $opt_x $opt_o);
use File::Basename;
use Cwd;

getopts("d:x:o:");

my $usage = "
	usage:
	$0

		-d <root Data file directory>
		-x <file eXtension to run program on>
		-o <analysis Output directory>
		
	This script will go into the data directory and identify all the files
	with the specificed extension and run the analysis with the specified
	parameters.  The output will be placed in the output directory, following
	the directory structure of the data file subdirectories.

	example:
		$0 
			-d \${PROJECT}/data/summary_tables/AllSamples \\
			-x .summary_table.tsv \\
			-o \${PROJECT}/analyses/MALR/MALR_run1

	It's recommended that the -o option be a full path, but if you
	just running it locally, you can also specify a '.' for current
	working directory.

";

if(!defined($opt_d) || !defined($opt_x) || !defined($opt_o)){
	die $usage;
}

my $DataRootDirectory=$opt_d;
my $TargetExtension=$opt_x;
my $AnalysisOutputDirectory=$opt_o;

# Remove / from end of paths
$DataRootDirectory=~s/\/$//;
$AnalysisOutputDirectory=~s/\/$//;

###############################################################################

print STDERR "\n";
print STDERR "Data Root Directory: $DataRootDirectory\n";
print STDERR "Target Extension: $TargetExtension\n";
print STDERR "Analysis Output Directory: $AnalysisOutputDirectory\n";
print STDERR "\n";

###############################################################################

print STDERR "Confirming existance of directories and scripts.\n";

if(-e $DataRootDirectory){
	print STDERR "\tFound $DataRootDirectory\n";
}else{
	die "\tCould not find $DataRootDirectory\n";
}

if($AnalysisOutputDirectory eq "." || $AnalysisOutputDirectory eq "./"){
	$AnalysisOutputDirectory=getcwd;
}

if(-e $AnalysisOutputDirectory){
	print STDERR "\tFound $AnalysisOutputDirectory\n";
}else{
	print STDERR "\tCould not find $AnalysisOutputDirectory\n";
	
	my ($dirname, $targetpath)=fileparse($AnalysisOutputDirectory);
	
	print STDERR "\t\tTrying to make $dirname in $targetpath\n";
	mkdir("$targetpath/$dirname");
	if(-e $AnalysisOutputDirectory){
		print STDERR "\t\tSuccess.  Continuing...\n";
	}else{
		die "\nError: Could not make $AnalysisOutputDirectory\n";
	}
}

###############################################################################

sub copy_subdirectories{
	my $src_dir=shift;
	my $dst_dir=shift;

	my @directory_list=split "\n", `find $src_dir -type d | sort`;

	print STDERR "\nCreating empty directories based on template source directory...\n\n";

	foreach my $dir(@directory_list){
		print STDERR "Source: $dir\n";

		# Remove source parent directory
		my $clean_dir=$dir;
		$clean_dir=~s/^$src_dir//;
		$clean_dir=~s/^\///;

		# Template path
		print STDERR "Target: $clean_dir\n";

		# Prepend destination path to template
		my $new_dir="$dst_dir/$clean_dir";
		print STDERR "Destin: $new_dir\n";

		# Make directory
		if(!(-e $new_dir)){
			print STDERR "  creating...";
			mkdir($new_dir);
			if(!(-e $new_dir)){
				die "Error creating $new_dir when copying template subdirectory.\n";
			}else{
				print STDERR "ok.\n";
			}
		}else{
			print STDERR "  exists.\n";
		}

		print STDERR "\n";
	}

	print STDERR "ok.\n";

}

sub find_common_root{
	my $dirA=shift;
	my $dirB=shift;

	my @compA=split "/", $dirA;
	my @compB=split "/", $dirB;

	my $numA=$#compA+1;
	my $numB=$#compB+1;
	
	my $max=($numA>$numB)?$numA:$numB;
	my @common;

	for(my $i=0; $i<$max; $i++){
		if($compA[$i] eq $compB[$i]){
			push @common, $compA[$i];
		}else{
			last;
		}
	}
	
	my $common_str=join "/", @common;

	print STDERR "Shared: $common_str\n";
	return($common_str);
}

sub pop_file{
	my $file_path=shift;
	my @file_path_arr=split "/", $file_path;
	my $file=pop @file_path_arr;
	my $popped_file_path=join "/", @file_path_arr;
	return("$popped_file_path/", $file);
}

sub climb_up_path{
	my $file_path=shift;
	my @file_path_arr=split "/", $file_path;

	for(my $i=0; $i<=$#file_path_arr; $i++){
		$file_path_arr[$i]="..";
	}	

	$file_path=join "/", @file_path_arr;
	return($file_path);
}

sub sym_link_target_files{
	my $src_dir=shift;
	my $dst_dir=shift;
	my $target_extension=shift;

	my $common_root=find_common_root($src_dir, $dst_dir);
	print STDERR "Source: $src_dir\n";
	print STDERR "Destin: $dst_dir\n";
	print STDERR "Common Root: $common_root\n";

	my @file_list=split "\n", `find $src_dir -type f | sort`;

	if($#file_list == -1){
		print STDERR "\n";
		print STDERR "Error: No target files found.\n";
		print STDERR "Note: For robust references, links to links will not be created.\n";
		die "Without targets found, there is no purpose in continuing...";
	}

	print STDERR "\nCreating symbolic links to target files...\n\n";

	foreach my $orig_file(@file_list){

	    	if($orig_file =~/$target_extension$/){
			print "$orig_file\n";

			# This is the location of where we want to create the link
			# This can be specified with an absolute path.
			my $link_path=$orig_file; 
			$link_path=~s/^$src_dir/$dst_dir/;

			# This is the location of the original file relative to the
			# location of the link.  
			my $src_path=$orig_file;
			# Separate filename from location
			my ($link_dir_path, $file)=pop_file($link_path);
			# Remove the beginning of the path, which is common
			# What is left is how many levels we have to climb up
			$link_dir_path=~s/^$common_root\///;
			# Change directories to climb to ..
			my $dotdots=climb_up_path($link_dir_path);
			my $src_from_common=$orig_file;
			$src_from_common=~s/^$common_root//;
			$src_path="$dotdots$src_from_common";

			# Create symbolic link now
			print STDERR "file location:\n\t$src_path\n";
			print STDERR "link location:\n\t$link_path\n";

			if(-l $link_path){
				print STDERR "Link already exists.  Replacing.\n";
				unlink $link_path;	
			}

			symlink($src_path, $link_path);
			print STDERR "\n";
		}
		print STDERR "\n";
	}
}

###############################################################################

copy_subdirectories($DataRootDirectory, $AnalysisOutputDirectory);

sym_link_target_files($DataRootDirectory, $AnalysisOutputDirectory, $TargetExtension);

###############################################################################

print STDERR "Done.\n";
