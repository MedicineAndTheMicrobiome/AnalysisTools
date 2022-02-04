#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_s $opt_S $opt_f $opt_F $opt_r $opt_c $opt_g $opt_p $opt_A $opt_B $opt_P $opt_o $opt_E $opt_t);
use File::Basename;
use Cwd;
use Digest::MD5;
use Sys::Hostname;

getopts("s:S:f:F:r:c:g:p:A:B:P:o:Ec:t:");

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

	-r <reference level filename>

	Variable Lists:
	[-c <covariates list>]
	[-g <\"grouped\" variables list>]

	Pairing:
	-p <pairings map>
	-A <Sample Group A, (Column name in pairings map file)>
	-B <Sample Group B, (Column name in pairings map file)>

	ALR (Abundance-Based) Related Options:
	[-P <number of ALR variables (for abundance-based analyses), default=$NUM_ALR_VARIABLES>]
	
	Output:
	-o <output directory>

	[-t <tag name>]

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
my $ReferenceLevelsFilename=$opt_r;
my $PairingMap=$opt_p;
my $Aname=$opt_A;
my $Bname=$opt_B;
my $NumALRVariables=$opt_P;
my $OutputDir=$opt_o;
my $DontAbort=$opt_E;
my $TagName=$opt_t;

my $Covariates=$opt_c;
my $GroupVar=$opt_g;

if(!defined($Covariates)){
	$Covariates="";
}

if(!defined($GroupVar)){
	$GroupVar="";
}

if(!defined($ReferenceLevelsFilename)){
	$ReferenceLevelsFilename="";
}

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

my $AnalysisName;
if($AnalysisName eq ""){

	if($GroupVar ne ""){
		($AnalysisName)=fileparse($GroupVar);
		$AnalysisName=~s/\.txt$//;
		$AnalysisName=~s/\.lst$//;
		$AnalysisName=~s/\.tsv$//;
	}else{
		$AnalysisName="result";
	}
}

if(!defined($TagName)){
	$TagName=$AnalysisName;
}

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
print STDERR "Tag Name:             $TagName\n";
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
	my $reference_level_file=shift;
	my $covariates=shift;
	my $variable_list=shift;
	my $model_name=shift;
	my $num_alr=shift;
	my $pair_map=shift;
	my $A_colname=shift;
	my $B_colname=shift;
	my $tag_name=shift;

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
	print STDERR "  Tag Name: $tag_name\n";
	print STDERR "\n";

	my $A_pred_B_OUT_DIR="alr_b_resp_a_pred";
	my $B_pred_A_OUT_DIR="alr_a_resp_b_pred";
	my $PRED_RESP_ANALYSIS="pred_resp_analysis";
	my $DIFF_OUT_DIR="alr_diff";
	my $cmd;

	if($covariates eq "" && $variable_list eq ""){
		$cmd="cat /dev/null > $output_dir/cov_var";
	}else{
		$cmd="cat $covariates $variable_list > $output_dir/cov_var";
	}
	run_command("Concatenate variables into full model list", "concat", $cmd, $output_dir);

	mkdir "$output_dir/abundance";
	mkdir "$output_dir/abundance/$A_pred_B_OUT_DIR";
	mkdir "$output_dir/abundance/$B_pred_A_OUT_DIR";
	mkdir "$output_dir/abundance/$PRED_RESP_ANALYSIS";
	mkdir "$output_dir/abundance/$DIFF_OUT_DIR";

	my $sumtabs;
	if($summary_table2 eq ""){
		$sumtabs="-s $summary_table";
	}else{
		$sumtabs="-s $summary_table -S $summary_table2";
	}

	my $reflvl;
	if($reference_level_file eq ""){
		$reflvl="";
	}else{
		$reflvl="-c $reference_level_file";
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
		-u $num_alr \
		-v $num_alr \
		-q $output_dir/cov_var \
		-x \";\" \
		-o $output_dir/abundance/$B_pred_A_OUT_DIR/$model_name \
		-t $tag_name \
		$reflvl
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
		-u $num_alr \
		-v $num_alr \
		-q $output_dir/cov_var \
		-x \";\" \
		-o $output_dir/abundance/$A_pred_B_OUT_DIR/$model_name \
		-t $tag_name \
		$reflvl
	";
	run_command("Fit $B_colname as Response", "B_as_resp", $cmd, "$output_dir/abundance/$A_pred_B_OUT_DIR");

	$cmd=
	"~/git/AnalysisTools/Profile/abundance_based/Compare_ALR_PredResp/Compare_ALR_PredResp.r \
	  -u $output_dir/abundance/$A_pred_B_OUT_DIR/$model_name.p_$A_colname.r_$B_colname.alr_as_pred.tp.coefs.tsv \
	  -v $output_dir/abundance/$B_pred_A_OUT_DIR/$model_name.p_$B_colname.r_$A_colname.alr_as_pred.coefs.tsv \
	  -x $output_dir/abundance/$A_pred_B_OUT_DIR/$model_name.p_$A_colname.r_$B_colname.alr_as_pred.tp.pvals.tsv \
	  -y $output_dir/abundance/$B_pred_A_OUT_DIR/$model_name.p_$B_colname.r_$A_colname.alr_as_pred.pvals.tsv \
	  -a $A_colname \
	  -b $B_colname \
	  -p 0.01 \
	  -d \
	  -t $tag_name \
	  -o $output_dir/abundance/$PRED_RESP_ANALYSIS/$model_name
	";
	run_command("Predictor/Response Analysis", "pred_resp_analysis", $cmd, "$output_dir/abundance/$PRED_RESP_ANALYSIS");

	$cmd=
	"~/git/AnalysisTools/Profile/abundance_based/Compare_ALR_PredResp/Compare_ALR_PredResp.r \
	  -u $output_dir/abundance/$A_pred_B_OUT_DIR/$model_name.p_$A_colname.r_$B_colname.alr_as_pred.tp.coefs.tsv \
	  -v $output_dir/abundance/$B_pred_A_OUT_DIR/$model_name.p_$B_colname.r_$A_colname.alr_as_pred.coefs.tsv \
	  -x $output_dir/abundance/$A_pred_B_OUT_DIR/$model_name.p_$A_colname.r_$B_colname.alr_as_pred.tp.pvals.tsv \
	  -y $output_dir/abundance/$B_pred_A_OUT_DIR/$model_name.p_$B_colname.r_$A_colname.alr_as_pred.pvals.tsv \
	  -a $A_colname \
	  -b $B_colname \
	  -p 0.05 \
	  -d \
	  -t $tag_name \
	  -o $output_dir/abundance/$PRED_RESP_ANALYSIS/$model_name
	";
	run_command("Predictor/Response Analysis", "pred_resp_analysis", $cmd, "$output_dir/abundance/$PRED_RESP_ANALYSIS");


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
		-u $num_alr \
		-x \";\" \
		-o $output_dir/abundance/$DIFF_OUT_DIR/$model_name \
		-t $tag_name \
		$reflvl
	";
	run_command("Fit B-A Difference", "paired_diff", $cmd, "$output_dir/abundance/$DIFF_OUT_DIR");


}

###############################################################################

sub run_distribution_based{

	my $output_dir=shift;
	my $summary_table=shift;
	my $summary_table2=shift;
	my $factor_file=shift;
	my $reference_level_file=shift;
	my $covariates=shift;
	my $variable_list=shift;
	my $model_name=shift;
	my $pair_map=shift;
	my $A_colname=shift;
	my $B_colname=shift;
	my $tag_name=shift;

	print STDERR "\n";
	print STDERR "Running Distribution Based Analyses:\n";
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
	print STDERR "  Tag Name: $tag_name\n";
	print STDERR "\n";

	my $cmd;
	my $DIV_DIFF="paired_div_diff_regr";
	my $STACKED_BP="stacked_barplots";

	if($covariates=="" && $variable_list==""){
		$cmd="touch $output_dir/cov_var"
	}else{
		$cmd="cat $covariates $variable_list > $output_dir/cov_var";
	}
	run_command("Concatenate variables into full model list", "concat", $cmd, $output_dir);

	mkdir "$output_dir/distribution";
	mkdir "$output_dir/distribution/$DIV_DIFF";
	mkdir "$output_dir/distribution/$STACKED_BP";

	my $sumtabs;
	if($summary_table2 eq ""){
		$sumtabs="-s $summary_table";
	}else{
		$sumtabs="-s $summary_table -S $summary_table2";
	}

	my $reflvl;
	if($reference_level_file eq ""){
		$reflvl="";
	}else{
		$reflvl="-c $reference_level_file";
	}

	my ($st_root_name)=File::Basename::fileparse($summary_table);
	$st_root_name=~s/\.summary_table\.tsv$//;

	#######################################################################

	$cmd=
	"~/git/AnalysisTools/Profile/distribution_based/Fit_TwoSample_Diversity_Regression/Fit_PairedDiff_Diversity_Regression.r \
		$sumtabs \
		-p $pair_map \
		-f $factor_file \
		-F 1 \
		-M $output_dir/cov_var \
		-B $B_colname \
		-A $A_colname \
		-q $output_dir/cov_var \
		-o $output_dir/distribution/$DIV_DIFF/$model_name \
		-t $tag_name \
		$reflvl
	";
	run_command("Fit Paired Diversity Difference Regression", "paired_div_diff_regr",
		$cmd, "$output_dir/distribution/$DIV_DIFF");
	
	#######################################################################

	my $CATEGORY_NAME="Category";
	
	$cmd=
	"~/git/AnalysisTools/Column/Make_Pair_Mapping_from_FactorInfo/Paired_to_Metadata.r \
		-p $pair_map \
		-o $output_dir/distribution/$STACKED_BP/paired_as_metadata.tsv \
		-a $A_colname \
		-b $B_colname \
		-c $CATEGORY_NAME
	";

	$cmd=
	run_command("Create Metadata Out of Paired Map", "paired_to_metadata",
		$cmd, "$output_dir/distribution/$STACKED_BP");

	#----------------------------------------------------------------------
	
	`echo $CATEGORY_NAME > $output_dir/distribution/$STACKED_BP/paired_category_name.txt`;

	#----------------------------------------------------------------------

	if($summary_table2 ne ""){
		$cmd=
		"~/git/AnalysisTools/Profile/SummaryTableUtilities/Join_Summary_Tables.r \
			-i $summary_table,$summary_table2 \
			-o $output_dir/distribution/$STACKED_BP/combined.summary_table.tsv
		";
		$cmd=
		run_command("Combine Summary Tables", "combine_summary_tables",
			$cmd, "$output_dir/distribution/$STACKED_BP");
	}else{
		$cmd=
		"cp $summary_table $output_dir/distribution/$STACKED_BP/combined.summary_table.tsv";
		run_command("Renaming Summary Table", "rename_summary_table",
			$cmd, "$output_dir/distribution/$STACKED_BP");
	}

	#----------------------------------------------------------------------
	
	foreach my $t ((2,3,4)){

		$cmd=
		"~/git/AnalysisTools/Profile/distribution_based/Plot_StackedBar/Plot_StackedBar.r \
			-i $output_dir/distribution/$STACKED_BP/combined.summary_table.tsv \
			-f $output_dir/distribution/$STACKED_BP/paired_as_metadata.tsv \
			-M $output_dir/distribution/$STACKED_BP/paired_category_name.txt \
			-o $output_dir/distribution/$STACKED_BP/paired_stacked_bp \
			-T $t \
			-t $tag_name \
			-s \";\"
		";
		$cmd=
		run_command("Create Metadata Out of Paired Map", "paired_to_metadata",
			$cmd, "$output_dir/distribution/$STACKED_BP");

	}


}

###############################################################################

sub run_distance_based{

	my $output_dir=shift;
	my $summary_table=shift;
	my $summary_table2=shift;
	my $factor_file=shift;
	my $reference_level_file=shift;
	my $covariates=shift;
	my $variable_list=shift;
	my $model_name=shift;
	my $pair_map=shift;
	my $A_colname=shift;
	my $B_colname=shift;
	my $tag_name=shift;

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
	print STDERR "  Tag Name: $tag_name\n";
	print STDERR "\n";

	my $cmd;
	my $dist_type="man";
	my $DIST_DIFF="paired_dist_regr";
	my $CLUST_CMP="clust_cmp";

	if($covariates=="" && $variable_list==""){
		$cmd="touch $output_dir/cov_var";
	}else{
		$cmd="cat $covariates $variable_list > $output_dir/cov_var";
	}
	run_command("Concatenate variables into full model list", "concat", $cmd, $output_dir);

	mkdir "$output_dir/distance";
	mkdir "$output_dir/distance/$DIST_DIFF";
	mkdir "$output_dir/distance/$CLUST_CMP";

	my $sumtabs;
	if($summary_table2 eq ""){
		$sumtabs="-s $summary_table";
	}else{
		$sumtabs="-s $summary_table -S $summary_table2";
	}

	my $reflvl;
	if($reference_level_file eq ""){
		$reflvl="";
	}else{
		$reflvl="-c $reference_level_file";
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
		-M $output_dir/cov_var \
		-B $B_colname \
		-A $A_colname \
		-q $output_dir/cov_var \
		-d man \
		-x \";\" \
		-o $output_dir/distance/$DIST_DIFF/$model_name \
		-t $tag_name \
		$reflvl
	";
	run_command("Fit Paired Distance Regression", $DIST_DIFF,
		$cmd, "$output_dir/distance/$DIST_DIFF");
	
	#######################################################################
	
	if($summary_table2 eq ""){
		$summary_table2=$summary_table;
	}

	$cmd=
	"~/git/AnalysisTools/Profile/distance_based/Cluster_Compare_TwoProfiles/Cluster_Compare_TwoProfiles.r \
		-a $summary_table \
		-b $summary_table2 \
		-A $A_colname \
		-B $B_colname \
		-m $pair_map \
		-o $output_dir/distance/$CLUST_CMP/$model_name \
		-d man \
		-t $tag_name
	";
	run_command("Cluster Compare", $CLUST_CMP, $cmd, "$output_dir/distance/$CLUST_CMP");


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
	$ReferenceLevelsFilename,
	$Covariates,
	$GroupVar,
	$AnalysisName,
	$NumALRVariables,
	$PairingMap,
	$Aname,
	$Bname,
	$TagName
);

run_distribution_based(
	$OutputDir,
	$SummaryTable,
	$SummaryTable2,
	$FactorFile,
	$ReferenceLevelsFilename,
	$Covariates,
	$GroupVar,
	$AnalysisName,
	$PairingMap,
	$Aname,
	$Bname,
	$TagName
);

run_distance_based(
	$OutputDir,
	$SummaryTable,
	$SummaryTable2,
	$FactorFile,
	$ReferenceLevelsFilename,
	$Covariates,
	$GroupVar,
	$AnalysisName,
	$PairingMap,
	$Aname,
	$Bname,
	$TagName
);




















































