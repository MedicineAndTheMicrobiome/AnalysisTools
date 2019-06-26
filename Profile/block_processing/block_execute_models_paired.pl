#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_s $opt_S $opt_f $opt_F $opt_c $opt_g $opt_p $opt_A $opt_B $opt_P $opt_o $opt_E);
use File::Basename;
use Cwd;
use Digest::MD5;
use Sys::Hostname;

getopts("s:S:f:F:c:g:p:A:B:P:o:E");

my $NUM_ALR_VARIABLES=15;

my $usage = "
	usage:
	$0

	Summary Table:
	-s <summary table>
	[-S <second summary table, in case pairings were in different files>]

	Factor File:
	-f <factor file>
	-F <column name of sample ids>

	Variable Lists:
	-c <covariates list>
	-g <\"grouped\" variables list>

	Pairing:
	-p <pairings map>
	-A <Sample Group A, (Column name in pairings map file)>
	-B <Sample Group B, (Column name in pairings map file)>

	ALR (Abundance-Based) Related Options:
	[-P <number of ALR variables (for abundance-based analyses), default=$NUM_ALR_VARIABLES>]
	
	Output:
	-o <output directory>

	[-E (Do not abort on error.  Keep going on other analyses)]

	This script will run a suite of analyses that use the following
	inputs:

		1.) Summary Table (.summary_table.tsv)
		2.) Factor Table (metadata.tsv)
		3.) Covariates, to control for  (covariates.txt)
		4.) Predictors/Responses (variable group)

	The following directories will be created in the output directory:

		abundance_based
			A predicts B
			B predicts A
			B-A

		distance_based
			Distance between A and B

";

if(
	!defined($opt_s) || 
	!defined($opt_f) || 
	!defined($opt_F) || 
	!defined($opt_c) || 
	!defined($opt_g) || 
	!defined($opt_p) || 
	!defined($opt_A) || 
	!defined($opt_B) || 
	!defined($opt_o)
){
	die $usage;
}


my $SummaryTable=$opt_s;
my $SummaryTable2=$opt_S;
my $FactorFile=$opt_f;
my $SampID_Colname=$opt_F;
my $Covariates=$opt_c;
my $GroupVar=$opt_g;
my $PairingMap=$opt_p;
my $Aname=$opt_A;
my $Bname=$opt_B;
my $NumALRVariables=$opt_P;
my $OutputDir=$opt_o;
my $DontAbort=$opt_E;

my $AnalysisName=$GroupVar;

if(!defined($SummaryTable2)){
	$SummaryTable2="";
}

if(!defined($opt_P)){
	$NumALRVariables=$NUM_ALR_VARIABLES;
}

if(!defined($opt_E)){
	$DontAbort=0;
}else{
	$DontAbort=1;
}

$AnalysisName=~s/\.txt$//;
$AnalysisName=File::Basename::fileparse($AnalysisName);

my $ABDNC_DIR="abundance_based";
my $DSTRB_DIR="distribution_based";
my $DSTNC_DIR="distance_based";

###############################################################################

print STDERR "\n";
print STDERR "Summary Table:        $SummaryTable\n";
print STDERR "Summary Table 2:      $SummaryTable2\n";
print STDERR "\n";
print STDERR "Factor File:          $FactorFile\n";
print STDERR "Sample ID Colname:    $SampID_Colname\n";
print STDERR "\n";
print STDERR "Covariates List:      $Covariates\n";
print STDERR "Group Variables List: $GroupVar\n";
print STDERR "Analysis Name:        $AnalysisName\n";
print STDERR "\n";
print STDERR "Pairings Map:         $PairingMap\n";
print STDERR "A Name:               $Aname\n";
print STDERR "B Name:               $Bname\n";
print STDERR "\n";
print STDERR "Output Directory:     $OutputDir\n";
print STDERR "\n";
print STDERR "Num ALR Variables:    $NumALRVariables\n";
print STDERR "\n";
print STDERR "Don't Abort on Errors: $DontAbort\n";
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
	my $summary_table2=shift;
	my $factor_file=shift;
	my $covariates=shift;
	my $variable_list=shift;
	my $model_name=shift;
	my $num_alr=shift;
	my $pair_map=shift;
	my $A_colname=shift;
	my $B_colname=shift;

	print STDERR "\n";
	print STDERR "Running Distance Based Analyses:\n";
	print STDERR "  Output Dir: $output_dir\n";
	print STDERR "  Summary Table 1: $summary_table\n";
	print STDERR "  Summary Table 2: $summary_table2\n";
	print STDERR "  Factor File: $factor_file\n";
	print STDERR "  Covariates File: $covariates\n";
	print STDERR "  Grouped Variable Fle: $variable_list\n";
	print STDERR "  Model Name: $model_name\n";
	print STDERR "  Num ALR Variables: $num_alr\n";
	print STDERR "  Pair Map: $pair_map\n";
	print STDERR "  A ColumnName: $A_colname\n";
	print STDERR "  B ColumnName: $B_colname\n";
	print STDERR "\n";

	my $A_pred_B_OUT_DIR="alr_a_resp_b_pred";
	my $B_pred_A_OUT_DIR="alr_b_resp_a_pred";
	my $DIFF_OUT_DIR="alr_diff";
	my $cmd;

	$cmd="cat $covariates $variable_list > $output_dir/cov_var";
	run_command("Concatenate variables into full model list", "concat", $cmd, $output_dir);

	mkdir "$output_dir/abundance";
	mkdir "$output_dir/abundance/$A_pred_B_OUT_DIR";
	mkdir "$output_dir/abundance/$B_pred_A_OUT_DIR";
	mkdir "$output_dir/abundance/$DIFF_OUT_DIR";

	my $sumtabs;
	if($summary_table2 eq ""){
		$sumtabs="-s $summary_table";
	}else{
		$sumtabs="-s $summary_table -S $summary_table2";
	}

	
	$cmd=
	"~/git/AnalysisTools/Profile/abundance_based/Fit_TwoSample_ALR_Regression/Fit_TwoSample_ALR_Regression.r \
		$sumtabs \
		-p $pair_map \
		-f $factor_file \
		-F 1 \
		-M $output_dir/cov_var \
		-e $A_colname \
		-g $B_colname \
		-q $output_dir/cov_var \
		-x \";\" \
		-o $output_dir/abundance/$A_pred_B_OUT_DIR
	";
	run_command("Fit $A_colname as Response", "A_as_resp", $cmd, "$output_dir/abundance/$B_pred_A_OUT_DIR");


	$cmd=
	"~/git/AnalysisTools/Profile/abundance_based/Fit_TwoSample_ALR_Regression/Fit_TwoSample_ALR_Regression.r \
		$sumtabs \
		-p $pair_map \
		-f $factor_file \
		-F 1 \
		-M $output_dir/cov_var \
		-e $B_colname \
		-g $A_colname \
		-q $output_dir/cov_var \
		-x \";\" \
		-o $output_dir/abundance/$A_pred_B_OUT_DIR
	";
	run_command("Fit $B_colname as Response", "B_as_resp", $cmd, "$output_dir/abundance/$A_pred_B_OUT_DIR");


	$cmd=
	"~/git/AnalysisTools/Profile/abundance_based/Fit_TwoSample_ALR_Regression/Fit_PairedDiff_ALR_Regression.r \
		$sumtabs \
		-p $pair_map \
		-f $factor_file \
		-F 1 \
		-M $output_dir/cov_var \
		-B $B_colname \
		-A $A_colname \
		-q $output_dir/cov_var \
		-x \";\" \
		-o $output_dir/abundance/$DIFF_OUT_DIR
	";
	run_command("Fit B-A Difference", "paired_diff", $cmd, "$output_dir/abundance/$DIFF_OUT_DIR");


}

###############################################################################

sub run_distance_based{

	my $output_dir=shift;
	my $summary_table=shift;
	my $summary_table2=shift;
	my $factor_file=shift;
	my $covariates=shift;
	my $variable_list=shift;
	my $model_name=shift;
	my $pair_map=shift;
	my $A_colname=shift;
	my $B_colname=shift;

	print STDERR "\n";
	print STDERR "Running Distance Based Analyses:\n";
	print STDERR "  Output Dir: $output_dir\n";
	print STDERR "  Summary Table 1: $summary_table\n";
	print STDERR "  Summary Table 2: $summary_table2\n";
	print STDERR "  Factor File: $factor_file\n";
	print STDERR "  Covariates File: $covariates\n";
	print STDERR "  Grouped Variable Fle: $variable_list\n";
	print STDERR "  Model Name: $model_name\n";
	print STDERR "  Pair Map: $pair_map\n";
	print STDERR "  A ColumnName: $A_colname\n";
	print STDERR "  B ColumnName: $B_colname\n";
	print STDERR "\n";

	my $cmd;
	my $dist_type="man";
	my $DIST_DIFF="paired_dist_regr";

	$cmd="cat $covariates $variable_list > $output_dir/cov_var";
	run_command("Concatenate variables into full model list", "concat", $cmd, $output_dir);

	mkdir "$output_dir/distance";
	mkdir "$output_dir/distance/$DIST_DIFF";

	my $sumtabs;
	if($summary_table2 eq ""){
		$sumtabs="-s $summary_table";
	}else{
		$sumtabs="-s $summary_table -S $summary_table2";
	}

	my ($st_root_name)=File::Basename::fileparse($summary_table);
	$st_root_name=~s/\.summary_table\.tsv$//;

	#######################################################################

	$cmd=
	"~/git/AnalysisTools/Profile/distance_based/Fit_Paired_Distance_Regression/Fit_Paired_Distance_Regression.r \
		$sumtabs \
		-p $pair_map \
		-f $factor_file \
		-F 1 \
		-M $covariates \
		-B $B_colname \
		-A $A_colname \
		-q $covariates \
		-x \";\" \
		-o $output_dir/distance/$DIST_DIFF
	";
	run_command("Fit Paired Distance Regression", "paired_dist_regr",
		$cmd, "$output_dir/distance/$DIST_DIFF");


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
	$SummaryTable2,
	$FactorFile,
	$Covariates,
	$GroupVar,
	$AnalysisName,
	$NumALRVariables,
	$PairingMap,
	$Aname,
	$Bname
);

run_distance_based(
	$OutputDir,
	$SummaryTable,
	$SummaryTable2,
	$FactorFile,
	$Covariates,
	$GroupVar,
	$AnalysisName,
	$PairingMap,
	$Aname,
	$Bname
);




















































