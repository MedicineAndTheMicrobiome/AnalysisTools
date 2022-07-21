#!/usr/bin/env perl

###############################################################################

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Std;
use vars qw($opt_d $opt_s $opt_p $opt_o);
use File::Basename;

my $SCREENED_MERGE_MATE_DIR="QC/QV_30/MergeMates/screened_fasta";

getopts("d:s:p:o:");
my $usage = "usage: 

$0 

	-d <SequencingRuns Directory>
	-s <'list of sequencing runs to search' file>
	-p <'list of project IDs to include', 0000 and 0001, negative controls automatically included>
	-o <Output directory>

	This script will look into the SequencingRuns directory, finding the sequencing runs
	specified in the sequencing runs list.  For each sequencing run, it will look into
	
	$SCREENED_MERGE_MATE_DIR

	The paired fasta files will be copied into the output directory.

	<output_directory>/<sequencing run>/<project_id>.*.paired.*.fasta

	Note:  The QC and MergeMates must have already been done for this.
	
	Once the script has completed successfully, the output directory should
	be ready for 16S_Clust_Gen_Pipeline/Generate_SampleID_to_Fasta_Map.pl
	
";

if(
	!defined($opt_d)||
	!defined($opt_s)||
	!defined($opt_p)||
	!defined($opt_o)
){
	die $usage;
}

my $SequencingRunsDirectory=$opt_d;
my $SequencingRunsListFname=$opt_s;
my $ProjectIDsListFname=$opt_p;
my $OutputDirectory=$opt_o;

###############################################################################

print STDERR "\n";
print STDERR "Sequencing Runs Directory: $SequencingRunsDirectory\n";
print STDERR "Sequencing Run Target List File Name: $SequencingRunsListFname\n";
print STDERR "Project IDs List File Name: $ProjectIDsListFname\n";
print STDERR "Output Directory: $OutputDirectory\n";
print STDERR "\n";

###############################################################################

sub load_file_list{
	my $list=shift;

	system("dos2unix $list");

	open(IN_FH, "<$list") || die "Could not open $list\n";
	my @load_file;
	while(<IN_FH>){
		chomp;
		if(substr($_, 0, 1) eq "#"){ next;} # Skip comments
		push @load_file, $_;
	}
	close(IN_FH);	
	return(\@load_file);
}

###############################################################################

sub make_dir{
	my $dir=shift;
	print STDERR "Making $dir.\n";
	if(!(-e $dir)){
		mkdir $dir;
	}else{
		print STDERR "$dir exists.\n";
	}
	if(!(-e $dir)){
		die "Could not make or find $dir\n";
	}
}

###############################################################################
# Load targeted sequence runs

my $target_seq_runs_ref=load_file_list($SequencingRunsListFname);
my $target_projects_ids_ref=load_file_list($ProjectIDsListFname);

unshift @{$target_projects_ids_ref}, ("0000", "0001");


# Clean up project ids a little and make sure they are unique
my %pid_hash;
foreach my $pid(@{$target_projects_ids_ref}){
	if(length($pid)>1){
		$pid_hash{$pid}=1;
	}
}
@{$target_projects_ids_ref}=sort keys %pid_hash;

###############################################################################

print STDOUT "\nIdentify Existence of Targeted Sequencing Runs:\n\n";

my $missing_sr_dirs=0;
foreach my $sr(@{$target_seq_runs_ref}){
	#print STDOUT "\t$sr\n";

	if(-e "$SequencingRunsDirectory/$sr"){
		print "$sr found.\n";
	}else{
		print "$sr NOT FOUND!!!\n";
		$missing_sr_dirs=1;
	}
}

if($missing_sr_dirs){
	print STDERR "\n";
	print STDERR "ERROR: Missing sequencing run directories.\n";
	print STDERR "Please correct error(s).\n";
	print STDERR "(You man need to fix typos in the sequencing\n";
	print STDERR " names, or remove them from the list to proceed.)\n";
	print STDERR "\n";
}

###############################################################################

# Try to make output dir
make_dir($OutputDirectory);

foreach my $sr(@{$target_seq_runs_ref}){
	
	my $seqrun_full_dir="$SequencingRunsDirectory/$sr";

	make_dir("$OutputDirectory/$sr");

	print STDERR "Looking for target files in: $SequencingRunsDirectory/$sr\n";

	foreach my $projid (@{$target_projects_ids_ref}){

		print STDERR "Looking for project ID: $projid\n";

		my @files=`ls $seqrun_full_dir/$SCREENED_MERGE_MATE_DIR/$projid.*.paired.*.fasta`;

		foreach my $file(@files){
			chomp $file;
			#print STDERR "Copying file:\n\t$file\nto\n\t$OutputDirectory/$sr\n\n";
			print STDERR `cp -pv $file $OutputDirectory/$sr`;
		}

	}

}

###############################################################################

print STDERR "Done.\n";

