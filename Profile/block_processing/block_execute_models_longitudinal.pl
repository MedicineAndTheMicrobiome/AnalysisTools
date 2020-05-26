#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_s $opt_f $opt_c $opt_t $opt_S $opt_g $opt_p $opt_a $opt_r $opt_o $opt_E);
use File::Basename;
use Cwd;
use Digest::MD5;
use Sys::Hostname;

getopts("s:f:c:t:S:g:p:a:r:o:E");

my $NUM_ALR_VARIABLES=35;

my $usage = "
	usage:
	$0

	-s <summary table file>
	-f <factor file>
	[-c <target variable list file>]
	
	-t <offset/time column name>
	-S <subject ID column name>
	[-g <group column name>]

	ALR (Abundance-Based) Related Options
	[-p <number of ALR variables (for abundance-based analyses), default=$NUM_ALR_VARIABLES>]
	[-a <filename for list of additional categories to include>]

	Factor Levels
	[-r <reference leveling file name (tsv file)>]  (Two columns, 1. Factor Name, 2. Reference Level)
	
	-o <output directory>

	[-E (Do not abort on error.  Keep going on other analyses)]

	This script will run a suite of analyses that use the following
	inputs:

		1.) Summary Table (.summary_table.tsv)
		2.) Factor Table (metadata.tsv)
		3.) Variables of Interest

	The following directories will be created in the output directory:

		abundance_based
		distribution_based
		distance_based



";

if(
	!defined($opt_s) || 
	!defined($opt_f) || 
	!defined($opt_t) ||
	!defined($opt_S) ||
	!defined($opt_o)
){
	die $usage;
}


my $SummaryTable=$opt_s;
my $FactorFile=$opt_f;
my $OffsetsColumn=$opt_t;
my $SubjectIDColumn=$opt_S;
my $OutputDir=$opt_o;

my $GroupColumn="";
if(defined($opt_g)){
	$GroupColumn=$opt_g;
}

my $NumALRVariables=$NUM_ALR_VARIABLES;
my $DontAbort=0;
my $AdditionalALRFile="";;
my $ReferenceRelevelingFile="";
my $TargetVariablesFile="";

if(defined($opt_c)){
	$TargetVariablesFile=$opt_c;
}

if(defined($opt_p)){
	$NumALRVariables=$opt_p;
}

if(defined($opt_a)){
	$AdditionalALRFile=$opt_a;
}

if(defined($opt_r)){
	$ReferenceRelevelingFile=$opt_r;
}

if(defined($opt_E)){
	$DontAbort=1;
}

my $AnalysisName=$TargetVariablesFile;
if($AnalysisName eq ""){
	$AnalysisName="result";
}else{
	$AnalysisName=~s/\.txt$//;
	$AnalysisName=~s/\.lst$//;
	$AnalysisName=File::Basename::fileparse($AnalysisName);
}

###############################################################################

my $ABDNC_DIR="abundance_based";
my $DSTRB_DIR="distribution_based";
my $DSTNC_DIR="distance_based";

###############################################################################

print STDERR "\n";
print STDERR "Summary Table:           $SummaryTable\n";
print STDERR "Factor File:             $FactorFile\n";
print STDERR "Target Variables:        $TargetVariablesFile\n";
print STDERR "Offsets Column Name:     $OffsetsColumn\n";
print STDERR "Subject ID Column Name:  $SubjectIDColumn\n";
print STDERR "Group ID Column Name:    $GroupColumn\n";
print STDERR "Output Directory:        $OutputDir\n";
print STDERR "\n";
print STDERR "Analysis Name:           $AnalysisName\n";
print STDERR "\n";
print STDERR "Reference Leveling File: $ReferenceRelevelingFile";
print STDERR "\n";
print STDERR "Num ALR Variables:       $NumALRVariables\n";
print STDERR "Additional ALR File:     $AdditionalALRFile\n";
print STDERR "\n";
print STDERR "Don't Abort on Errors:   $DontAbort\n";
print STDERR "\n";

###############################################################################

sub run_command{
	my $cmd_name=shift;
	my $cmd_name_shrt=shift;
	my $cmd_str=shift;
	my $chkpt_dir=shift;

	print STDERR "Command Name: $cmd_name\n";
	print STDERR "Command String: $cmd_str\n";
	print STDERR "Checkpoint File's Directory: $chkpt_dir\n";

	my $cmd_str_hashid=Digest::MD5::md5_hex($cmd_str);
	my $checkpt_file="$chkpt_dir/$cmd_name_shrt.$cmd_str_hashid.CHKPT";

	print STDERR "Checkpoint File Name: $checkpt_file\n";

	# Check if command needs to be executed
	if(-e $checkpt_file){
		open(CHK_FH, "tail -n 1 $checkpt_file |") || die "Could not open $checkpt_file.\n";
		my $last_line=<CHK_FH>;
		close(CHK_FH);
		
		$last_line=~s/\s+//g;
		if($last_line eq "Completed."){
			print STDERR "$cmd_name ($checkpt_file) was successfully completed previous. Skipping.\n";
			return;
		}
	}

	# Prep command
	my $clean_cmd=$cmd_str;
	$clean_cmd=~s/\s+/ /g;

	print STDERR "'$clean_cmd'\n";
	
	# Save command to file
	open(CHK_FH, ">$checkpt_file") || die "Could not open $checkpt_file.\n";
	print CHK_FH "$clean_cmd\n";
	close(CHK_FH);

	# Execute command
	my $result=`$clean_cmd 2>&1`;
	my $exit_code=$?;

	# Save output to file
	my $logname="$chkpt_dir/$cmd_name_shrt\.LOG.TXT";
	open(LOG_FH, ">$logname") || die "Could not open $logname.\n";
	print LOG_FH "$result";
	close(LOG_FH);

	# Write status to end of checkpoint file
	if($exit_code==0){
		open(CHK_FH, ">>$checkpt_file") || die "Could not open $checkpt_file.\n";
		print CHK_FH "\nCompleted.\n";
		close(CHK_FH);

		if(-e "$checkpt_file\.failed"){
			unlink "$checkpt_file\.failed"
		}

		return;
	}else{
		open(CHK_FH, ">>$checkpt_file\.failed") || die "Could not open $checkpt_file.\n";

		my $tail_log=`tail -n 10 $logname`;

		print(CHK_FH "$tail_log\n");
		print CHK_FH "\nFailed.\n";

		close(CHK_FH);

		print STDERR "Error: $cmd_name returned with non-zero error code.\n";
		print STDERR "\n";
		print STDERR "$cmd_str\n";
		
		if(!$DontAbort){
			exit;
		}
	}
	
}


###############################################################################

sub run_abundance_based{

	my $output_dir=shift;
	my $summary_table=shift;
	my $factor_file=shift;
	my $target_var_file=shift;
	my $model_name=shift;
	my $offsets_colname=shift;
	my $subjectid_colname=shift;
	my $groups_colname=shift;
	my $reference_leveling_file=shift;
	my $num_alr=shift;
	my $additional_alr_file=shift;

	my $cmd;

	mkdir "$output_dir/abundance";
	
	my $add_alr="";
	if($AdditionalALRFile ne ""){
		$add_alr="-a $AdditionalALRFile";
	}

	my $add_reflev="";
	if($ReferenceRelevelingFile ne ""){
		$add_reflev="-r $ReferenceRelevelingFile";
	}

	my $add_groups="";
	if($groups_colname ne ""){
		$add_groups="";
	}

	my $add_target_var_file="";
	if($target_var_file ne ""){
		$add_target_var_file="-m $target_var_file";
	}

	$cmd=
	"~/git/AnalysisTools/Profile/abundance_based/Plot_Longitudinal_Samples_ALR/Plot_Longitudinal_Samples_ALR.r \
		-s $summary_table \
		-v $num_alr \
		-x \";\" \
		-f $factor_file \
		-t $offsets_colname \
		-i $subjectid_colname \
		-o $output_dir/abundance/$model_name
		$add_alr
		$add_reflev 
		$add_groups 
		$add_target_var_file
	";
	run_command("Longitudinal ALR", "longit_alr", $cmd, 
		"$output_dir/abundance");

}

###############################################################################

sub run_distribution_based{

	my $output_dir=shift;
	my $summary_table=shift;
	my $factor_file=shift;
	my $target_var_file=shift;
	my $model_name=shift;
	my $offsets_colname=shift;
	my $subjectid_colname=shift;
	my $groups_colname=shift;
	my $reference_leveling_file=shift;

	my $cmd;

	mkdir "$output_dir/distribution";

	my $add_reflev="";
	if($reference_leveling_file ne ""){
		$add_reflev="";
	}

	my $add_groups="";
	if($groups_colname ne ""){
		$add_groups="";
	}

	#######################################################################
	# Diversity
	$cmd=
	"~/git/AnalysisTools/Profile/distribution_based/Plot_Longitudinal_Samples_Dist/Plot_Longitudinal_Samples_Dist.r \
                -s $summary_table \
		-f $factor_file \
		-t $offsets_colname \
		-i $subjectid_colname \
             	-o $output_dir/distribution/$model_name \
		-d tail \
		$add_reflev
		$add_groups
	";
	run_command("Plot Diversity Longitudinal", "longit_diversity", 
		$cmd, "$output_dir/distribution");

}


###############################################################################

sub run_distance_based{

	my $output_dir=shift;
	my $summary_table=shift;
	my $factor_file=shift;
	my $target_var_file=shift;
	my $model_name=shift;
	my $offsets_colname=shift;
	my $subjectid_colname=shift;
	my $groups_colname=shift;
	my $reference_leveling_file=shift;

	my $LONGIT_PLOT="longit_plot";

	my $cmd;
	my $dist_type="man";

	mkdir "$output_dir/distance";

	my $add_reflev="";
	if($reference_leveling_file ne ""){
		$add_reflev="";
	}

	my $add_groups="";
	if($groups_colname ne ""){
		$add_groups="";
	}


	$cmd=
	"~/git/AnalysisTools/Profile/distance_based/Plot_Longitudinal_Samples/Plot_Longitudinal_Samples.r \
		-s $summary_table \
		-f $factor_file \
		-t $offsets_colname \
		-i $subjectid_colname \
		-o $output_dir/distance/$model_name \
		-d man \
		$add_reflev \
		$add_groups
	";
	run_command("Plot Distance Longitudinal", "longit_dist_plot",
		$cmd, "$output_dir/distance");


}


###############################################################################

#print STDERR "Summary Table:   $SummaryTable\n";
#print STDERR "Factor File:     $FactorFile\n";
#print STDERR "Covariates List: $Covariates\n";
#print STDERR "Group Variables List: $GroupVar\n";
#print STDERR "Output Directory: $OutputDir\n";
#


if(!(-e $OutputDir)){
	mkdir $OutputDir;
}


run_abundance_based(
	$OutputDir,
	$SummaryTable,
	$FactorFile,
	$TargetVariablesFile,
	$AnalysisName,
	$OffsetsColumn,
	$SubjectIDColumn,
	$GroupColumn,
	$ReferenceRelevelingFile,
	$NumALRVariables,
	$AdditionalALRFile
);


run_distribution_based(
	$OutputDir,
	$SummaryTable,
	$FactorFile,
	$TargetVariablesFile,
	$AnalysisName,
	$OffsetsColumn,
	$SubjectIDColumn,
	$GroupColumn,
	$ReferenceRelevelingFile
);

run_distance_based(
	$OutputDir,
	$SummaryTable,
	$FactorFile,
	$TargetVariablesFile,
	$AnalysisName,
	$OffsetsColumn,
	$SubjectIDColumn,
	$GroupColumn,
	$ReferenceRelevelingFile
);




















































