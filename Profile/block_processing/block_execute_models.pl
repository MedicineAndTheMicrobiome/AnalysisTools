#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_s $opt_f $opt_c $opt_g $opt_p $opt_a $opt_o $opt_E $opt_r $opt_t);
use File::Basename;
use Cwd;
use Digest::MD5;
use Sys::Hostname;

getopts("s:f:c:g:p:a:o:Er:t:");

my $NUM_ALR_VARIABLES=35;

my $usage = "
	usage:
	$0

	-s <summary table file>
	-f <factor file>
	[-c <covariates list file>]
	[-g <\"grouped\" variables list file>]

	ALR (Abundance-Based) Related Options
	[-p <number of ALR variables (for abundance-based analyses), default=$NUM_ALR_VARIABLES>]
	[-a <filename for list of additional categories to include>]

	Factor Levels
	[-r <reference leveling file name (tsv file)>]  (Two columns, 1. Factor Name, 2. Reference Level)
	
	-o <output directory>

	[-E (Do not abort on error.  Keep going on other analyses)]
	[-t <tag name>]

	This script will run a suite of analyses that use the following
	inputs:

		1.) Summary Table (.summary_table.tsv)
		2.) Factor Table (metadata.tsv)
		3.) Covariates, to control for  (covariates.txt)
		4.) Predictors/Responses (variable group)

	The following directories will be created in the output directory:

		abundance_based
			ALR as predictor
			ALR as response
			Pred/Resp analysis

		distribution_based
			Diversity as predictor
			Diversity as response
			Pred/Resp analysis

		distance_based
			Permanova
			Cluster MLL

	Use the -t <tag name> option to help watermark the individual pages in the output pdf file.


";

if(
	!defined($opt_s) || 
	!defined($opt_f) || 
	!defined($opt_o)
){
	die $usage;
}


my $SummaryTable=$opt_s;
my $FactorFile=$opt_f;
my $Covariates=$opt_c;
my $OutputDir=$opt_o;

my $GroupVar="";
my $AnalysisName=$GroupVar;
my $NumALRVariables=$NUM_ALR_VARIABLES;
my $DontAbort;
my $AdditionalALRFile;
my $ReferenceRelevelingFile;
my $TagName="";

if(defined($opt_g)){
	$GroupVar=$opt_g;
}

if(defined($opt_p)){
	$NumALRVariables=$opt_p;
}

if(defined($opt_a)){
	$AdditionalALRFile=$opt_a;
}else{
	$AdditionalALRFile="";
}

if(defined($opt_E)){
	$DontAbort=1;
}else{
	$DontAbort=0;
}

if(defined($opt_r)){
	$ReferenceRelevelingFile=$opt_r;
}else{
	$ReferenceRelevelingFile="";
}

$AnalysisName=~s/\.txt$//;
$AnalysisName=File::Basename::fileparse($AnalysisName);

if($AnalysisName eq ""){
	$AnalysisName="result";
}

if(defined($opt_t)){
	$TagName=$opt_t;
}else{
	$TagName=$AnalysisName;
}

my $ABDNC_DIR="abundance_based";
my $DSTRB_DIR="distribution_based";
my $DSTNC_DIR="distance_based";

my $AS_PRED="as_pred";
my $AS_RESP="as_resp";
my $PVR="pred_vs_resp";

my $PERM="permanova";
my $CLST_MLL="cluster_mll";

###############################################################################

print STDERR "\n";
print STDERR "Summary Table:           $SummaryTable\n";
print STDERR "Factor File:             $FactorFile\n";
print STDERR "Covariates List:         $Covariates\n";
print STDERR "Group Variables List:    $GroupVar\n";
print STDERR "Output Directory:        $OutputDir\n";
print STDERR "Analysis Name:           $AnalysisName\n";
print STDERR "\n";
print STDERR "Reference Leveling File: $ReferenceRelevelingFile";
print STDERR "\n";
print STDERR "Num ALR Variables:       $NumALRVariables\n";
print STDERR "Additional ALR File:     $AdditionalALRFile\n";
print STDERR "\n";
print STDERR "Don't Abort on Errors:   $DontAbort\n";
print STDERR "Tag Name:                $TagName\n";
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
	my $covariates=shift;
	my $variable_list=shift;
	my $model_name=shift;
	my $num_alr=shift;
	my $tag_name=shift;


	my $PRED_OUT_DIR="alr_as_pred";
	my $RESP_OUT_DIR="alr_as_resp";
	my $COMP_DIR="alr_pred_resp_comp";
	my $cmd;

	$cmd="cat $covariates $variable_list > $output_dir/cov_var";
	run_command("Concatenate variables into full model list", "concat", $cmd, $output_dir);

	mkdir "$output_dir/abundance";
	mkdir "$output_dir/abundance/$RESP_OUT_DIR";
	mkdir "$output_dir/abundance/$PRED_OUT_DIR";
	mkdir "$output_dir/abundance/$COMP_DIR";
	
	my $add_alr="";
	if($AdditionalALRFile ne ""){
		$add_alr="-a $AdditionalALRFile";
	}

	my $add_reflev="";
	if($ReferenceRelevelingFile ne ""){
		$add_reflev="-r $ReferenceRelevelingFile";
	}
	
	$cmd=
	"~/git/AnalysisTools/Profile/abundance_based/Fit_ALR_as_Response/Fit_ALR_as_Response.r \
                -s $summary_table \
                -f $factor_file \
                -M $output_dir/cov_var \
                -q $output_dir/cov_var \
                -p $num_alr \
                -o $output_dir/abundance/$RESP_OUT_DIR/$model_name  \
                -x \";\" \
		-t $tag_name \
		$add_alr $add_reflev
	";
	run_command("Fit ALR as Response", "alr_as_resp", $cmd, "$output_dir/abundance/$RESP_OUT_DIR");

	if($variable_list ne ""){

		$cmd=
		"~/git/AnalysisTools/Profile/abundance_based/Fit_ALR_as_Predictor/Fit_ALR_as_Predictor.r \
			-s $summary_table \
			-f $factor_file \
			-y $variable_list \
			-c $covariates \
			-q $output_dir/cov_var \
			-p $num_alr \
			-o $output_dir/abundance/$PRED_OUT_DIR/$model_name \
			-x \";\" \
			-t $tag_name \
			$add_alr $add_reflev
		";
		run_command("Fit ALR as Predictor", "alr_as_pred", $cmd, "$output_dir/abundance/$PRED_OUT_DIR");


		$cmd=
		"~/git/AnalysisTools/Profile/abundance_based/Compare_ALR_PredResp/Compare_ALR_PredResp.r \
			-x $output_dir/abundance/$PRED_OUT_DIR/$model_name.alr_as_pred.pvals.tsv \
			-y $output_dir/abundance/$RESP_OUT_DIR/$model_name.alr_as_resp.pvals.tsv \
			-u $output_dir/abundance/$PRED_OUT_DIR/$model_name.alr_as_pred.coefs.tsv \
			-v $output_dir/abundance/$RESP_OUT_DIR/$model_name.alr_as_resp.coefs.tsv \
			-o $output_dir/abundance/$COMP_DIR/$model_name.alr \
			-p .025 \
			-t $tag_name
		";
		run_command("Compare ALR Pred/Resp", "alr_pred_resp_comp.025", $cmd, "$output_dir/abundance/$COMP_DIR");	

		$cmd=
		"~/git/AnalysisTools/Profile/abundance_based/Compare_ALR_PredResp/Compare_ALR_PredResp.r \
			-x $output_dir/abundance/$PRED_OUT_DIR/$model_name.alr_as_pred.pvals.tsv \
			-y $output_dir/abundance/$RESP_OUT_DIR/$model_name.alr_as_resp.pvals.tsv \
			-u $output_dir/abundance/$PRED_OUT_DIR/$model_name.alr_as_pred.coefs.tsv \
			-v $output_dir/abundance/$RESP_OUT_DIR/$model_name.alr_as_resp.coefs.tsv \
			-o $output_dir/abundance/$COMP_DIR/$model_name.alr \
			-p .050 \
			-t $tag_name
		";
		run_command("Compare ALR Pred/Resp", "alr_pred_resp_comp.050", $cmd, "$output_dir/abundance/$COMP_DIR");	

		$cmd=
		"~/git/AnalysisTools/Profile/abundance_based/Compare_ALR_PredResp/Compare_ALR_PredResp.r \
			-x $output_dir/abundance/$PRED_OUT_DIR/$model_name.alr_as_pred.pvals.tsv \
			-y $output_dir/abundance/$RESP_OUT_DIR/$model_name.alr_as_resp.pvals.tsv \
			-u $output_dir/abundance/$PRED_OUT_DIR/$model_name.alr_as_pred.coefs.tsv \
			-v $output_dir/abundance/$RESP_OUT_DIR/$model_name.alr_as_resp.coefs.tsv \
			-o $output_dir/abundance/$COMP_DIR/$model_name.alr \
			-p .005 \
			-t $tag_name
		";
		run_command("Compare ALR Pred/Resp", "alr_pred_resp_comp.005", $cmd, "$output_dir/abundance/$COMP_DIR");	

	}

}

###############################################################################

sub run_distribution_based{

	my $output_dir=shift;
	my $summary_table=shift;
	my $factor_file=shift;
	my $covariates=shift;
	my $variable_list=shift;
	my $model_name=shift;
	my $tag_name=shift;

	my $PRED_OUT_DIR="div_as_pred";
	my $RESP_OUT_DIR="div_as_resp";
	my $COMP_DIR="div_pred_resp_comp";
	my $STCK_BAR_DIR="stackedbar_plot";
	my $RANK_ABND_DIR="rank_abund_plot";
	my $cmd;


	$cmd="cat $covariates $variable_list > $output_dir/cov_var";
	run_command("Concatenate variables into full model list", "concat", $cmd, $output_dir);

	mkdir "$output_dir/distribution";
	mkdir "$output_dir/distribution/$RESP_OUT_DIR";
	mkdir "$output_dir/distribution/$PRED_OUT_DIR";
	mkdir "$output_dir/distribution/$COMP_DIR";
	mkdir "$output_dir/distribution/$STCK_BAR_DIR";
	mkdir "$output_dir/distribution/$RANK_ABND_DIR";

	my $add_reflev="";
	if($ReferenceRelevelingFile ne ""){
		$add_reflev="-r $ReferenceRelevelingFile";
	}

	#######################################################################
	# Diversity
	$cmd=
	 "~/git/AnalysisTools/Profile/distribution_based/Fit_Diversity_as_Response/Fit_Diversity_as_Response.r \
                -s $summary_table \
                -f $factor_file \
                -M $output_dir/cov_var \
                -q $output_dir/cov_var \
             	-o $output_dir/distribution/$RESP_OUT_DIR/$model_name \
		-t $tag_name \
		$add_reflev
	";
	run_command("Fit Diversity as Response", "div_as_resp", $cmd, "$output_dir/distribution/$RESP_OUT_DIR");


	if($variable_list ne ""){
		$cmd=
		"~/git/AnalysisTools/Profile/distribution_based/Fit_Diversity_as_Predictor/Fit_Diversity_as_Predictor.r \
			-s $summary_table \
			-f $factor_file \
			-y $variable_list \
			-c $covariates \
			-q $output_dir/cov_var \
			-o $output_dir/distribution/$PRED_OUT_DIR/$model_name \
			-t $tag_name \
			$add_reflev
			
		";
		run_command("Fit Diversity as Predictor", "div_as_pred", $cmd, "$output_dir/distribution/$PRED_OUT_DIR");

		$cmd=
		"~/git/AnalysisTools/Profile/abundance_based/Compare_ALR_PredResp/Compare_ALR_PredResp.r \
			-x $output_dir/distribution/$PRED_OUT_DIR/$model_name.div_as_pred.pvals.tsv \
			-y $output_dir/distribution/$RESP_OUT_DIR/$model_name.div_as_resp.pvals.tsv \
			-u $output_dir/distribution/$PRED_OUT_DIR/$model_name.div_as_pred.coefs.tsv \
			-v $output_dir/distribution/$RESP_OUT_DIR/$model_name.div_as_resp.coefs.tsv \
			-o $output_dir/distribution/$COMP_DIR/$model_name.div \
			-a Diversity \
			-p .025 \
			-t $tag_name
		";
		run_command("Compare Diversity Pred/Resp", "div_pred_resp_comp.025", $cmd, 
			"$output_dir/distribution/$COMP_DIR");

		$cmd=
		"~/git/AnalysisTools/Profile/abundance_based/Compare_ALR_PredResp/Compare_ALR_PredResp.r \
			-x $output_dir/distribution/$PRED_OUT_DIR/$model_name.div_as_pred.pvals.tsv \
			-y $output_dir/distribution/$RESP_OUT_DIR/$model_name.div_as_resp.pvals.tsv \
			-u $output_dir/distribution/$PRED_OUT_DIR/$model_name.div_as_pred.coefs.tsv \
			-v $output_dir/distribution/$RESP_OUT_DIR/$model_name.div_as_resp.coefs.tsv \
			-o $output_dir/distribution/$COMP_DIR/$model_name.div \
			-a Diversity \
			-p .050 \
			-t $tag_name
		";
		run_command("Compare Diversity Pred/Resp", "div_pred_resp_comp.050", $cmd, 
			"$output_dir/distribution/$COMP_DIR");
	}

	#######################################################################
	# Plots
	$cmd=
	"~/git/AnalysisTools/Profile/distribution_based/Plot_StackedBar/Plot_StackedBar.r \
		-i $summary_table \
		-f $factor_file \
		-M $output_dir/cov_var \
		-o $output_dir/distribution/$STCK_BAR_DIR/$model_name \
		-s \";\" \
		-t $tag_name
	";
	run_command("Plot Stacked Bar Plots", "stacked_bp", $cmd, "$output_dir/distribution/$STCK_BAR_DIR");

	$cmd=
	"~/git/AnalysisTools/Profile/distribution_based/Plot_RankAbundance_wFactors/Plot_RankAbundance_wFactors.r \
		-i $summary_table \
		-f $factor_file \
		-M $output_dir/cov_var \
		-o $output_dir/distribution/$RANK_ABND_DIR/$model_name \
		-s \";\" \
		-t $tag_name
	";
	run_command("Plot Rank Abundance Plots", "rank_abnd", $cmd, "$output_dir/distribution/$RANK_ABND_DIR");

}


###############################################################################

sub run_distance_based{

	my $output_dir=shift;
	my $summary_table=shift;
	my $factor_file=shift;
	my $covariates=shift;
	my $variable_list=shift;
	my $model_name=shift;
	my $tag_name=shift;

	my $DISTMAT_DIR="dist_mat";
	my $PERMA_DIR="permanova";
	my $CLUST_MLL_DIR="clust_mll";
	my $INFL_BY_FACTOR="infl_by_factors";
	my $cmd;
	my $dist_type="man";

	my $add_reflev="";
	if($ReferenceRelevelingFile ne ""){
		$add_reflev="-r $ReferenceRelevelingFile";
	}

	$cmd="cat $covariates $variable_list > $output_dir/cov_var";
	run_command("Concatenate variables into full model list", "concat", $cmd, $output_dir);

	mkdir "$output_dir/distance";
	mkdir "$output_dir/distance/$DISTMAT_DIR";
	mkdir "$output_dir/distance/$PERMA_DIR";
	mkdir "$output_dir/distance/$CLUST_MLL_DIR";
	mkdir "$output_dir/distance/$INFL_BY_FACTOR";

	my ($st_root_name)=File::Basename::fileparse($summary_table);
	$st_root_name=~s/\.summary_table\.tsv$//;

	#######################################################################
	# Compute distance matrix
	$cmd=
	"~/git/AnalysisTools/Profile/SummaryTableUtilities/Create_Distance_Matrix.r \
                -i $summary_table \
                -d $dist_type \
             	-o $output_dir/distance/$DISTMAT_DIR/$st_root_name
	";
	run_command("Compute Distance Matrix", "dist_mat", $cmd, "$output_dir/distance/$DISTMAT_DIR");

	# Compute permanova
	$cmd=
	"~/git/AnalysisTools/Profile/distance_based/Permanova/Permanova.r \
                -d $output_dir/distance/$DISTMAT_DIR/$st_root_name.man.distmat \
                -f $factor_file \
		-M $output_dir/cov_var \
		-o $output_dir/distance/$PERMA_DIR/$model_name \
		-q $output_dir/cov_var \
		-t $tag_name
	";
	run_command("Permanova", "perma", $cmd, "$output_dir/distance/$PERMA_DIR");
		

	#######################################################################
	# Compute cluster multinomial log linear model
	$cmd=
	"~/git/AnalysisTools/Profile/distance_based/Cluster_MultinomLogLin/Cluster_MultinomLogLin.r \
		-i $summary_table \
		-f $factor_file \
		-M $output_dir/cov_var \
		-q $output_dir/cov_var \
		-o $output_dir/distance/$CLUST_MLL_DIR/$model_name \
		-d $dist_type \
		-t $tag_name \
		$add_reflev
	";
	run_command("Fit Multinomial Log-Linear Model", "clust_mll", $cmd, "$output_dir/distance/$CLUST_MLL_DIR");

	$cmd=
	"~/git/AnalysisTools/Profile/distance_based/Plot_ClusterHeatmap/Plot_ClusterHeatmap.r \
		-i $summary_table \
		-f $factor_file \
		-d $dist_type \
		-M $output_dir/cov_var \
		-c $output_dir/distance/$CLUST_MLL_DIR/$model_name\.$dist_type\.cl_mll.used_samp.lst \
		-o $output_dir/distance/$CLUST_MLL_DIR/$model_name \
		-t $tag_name
	";
	run_command("Plot Metadata Heatmap", "cl_hmp", $cmd, "$output_dir/distance/$CLUST_MLL_DIR");

	$cmd=
	"~/git/AnalysisTools/Profile/distance_based/CH_Cluster_Stop/CH_Cluster_Stop.r \
		-i $summary_table \
		-d $dist_type \
		-c $output_dir/distance/$CLUST_MLL_DIR/$model_name\.$dist_type\.cl_mll.used_samp.lst \
		-o $output_dir/distance/$CLUST_MLL_DIR/$model_name \
		-t $tag_name 
	";
	run_command("Compute Calinski-Harabasz Stopping", "ch_stop", $cmd, "$output_dir/distance/$CLUST_MLL_DIR");

	$cmd=
	"~/git/AnalysisTools/Profile/distance_based/Cluster_Influencers/Cluster_Influencers.r \
		-i $summary_table \
		-d $dist_type \
		-l $output_dir/distance/$CLUST_MLL_DIR/$model_name\.$dist_type\.cl_mll.used_samp.lst \
		-o $output_dir/distance/$CLUST_MLL_DIR/$model_name \
		-t $tag_name
	";
	run_command("Compute Cluster Influencers", "cl_inf", $cmd, "$output_dir/distance/$CLUST_MLL_DIR");
	
	#######################################################################
	# Cluster influencers by factors

	$cmd=
	"~/git/AnalysisTools/Profile/distance_based/Cluster_Influencers/Cluster_Influencers.r \
		-i $summary_table \
		-d $dist_type \
		-o $output_dir/distance/$INFL_BY_FACTOR/$model_name \
		-n $output_dir/cov_var \
		-f $factor_file \
		-t $tag_name
	";
	run_command("Compute Cluster Influencers by Factors", "cl_inf_by_factors",
		$cmd, "$output_dir/distance/$INFL_BY_FACTOR");


}


###############################################################################

print STDERR "Summary Table:   $SummaryTable\n";
print STDERR "Factor File:     $FactorFile\n";
print STDERR "Covariates List: $Covariates\n";
print STDERR "Group Variables List: $GroupVar\n";
print STDERR "Output Directory: $OutputDir\n";

if(!(-e $OutputDir)){
	mkdir $OutputDir;
}


###############################################################################

my $cmd="cat $Covariates $GroupVar > $OutputDir/all_req_var";
run_command("Combine Cov/Grp variables for NA prescreen", "na_prescreen_concat_var",
	$cmd, "$OutputDir");

# Prescreen factors/summary table for NAs
my $cmd=
"~/git/AnalysisTools/Metadata/RemoveNAs/Prescreen_Sample_woFactorNAs/Prescreen_Sample_woFactorNAs.r \
	-s $SummaryTable \
	-f $FactorFile \
	-d $OutputDir \
	-t $OutputDir/all_req_var \
	-q $OutputDir/all_req_var
";
run_command("Run NA prescreen", "na_prescreen", $cmd, "$OutputDir");

print STDERR "\n**************************************************************\n";

my ($stname, $stpath)=fileparse($SummaryTable);
my $screened_summary_table="$OutputDir/$stname";
$screened_summary_table=~s/\.summary_table\.tsv$//;
$screened_summary_table="$screened_summary_table.prescr.summary_table.tsv";
print STDERR "Prescreened Summary Table: $screened_summary_table\n";

my ($fcname, $fcpath)=fileparse($FactorFile);
my $screened_factor_file="$OutputDir/$fcname";
$screened_factor_file=~s/\.tsv$//;
$screened_factor_file="$screened_factor_file.prescr.tsv";
print "Prescreened Factor File: $screened_factor_file\n";

print STDERR "\n**************************************************************\n";

# Remove factors lost due to NA screening
my $cmd=
"cat $screened_factor_file | \
	~/git/AnalysisTools/Column/list_column_headers.pl | \
	grep -v '^\$' > \
	$OutputDir/all_avail_var";
#my $cmd="lch $screened_factor_file | grep -v '^\$' > $OutputDir/all_avail_var; echo done.";
run_command("Get Post-NA Removal Variables", "get_postna_rem_var", $cmd, "$OutputDir");

if($Covariates ne ""){
	my $cmd="grep -f $OutputDir/all_avail_var $Covariates > $OutputDir/covar.screened; echo done.";
	run_command("Screen Covariates", "scr_covar", $cmd, "$OutputDir");
}else{
	my $cmd="touch $OutputDir/covar.screened; echo done.";
	print STDERR "Covariates not specified.  Skipping NA removal Covariate screen.\n";
	run_command("Screen Covariate Variables", "scr_covar", $cmd, "$OutputDir");
}

if ($GroupVar ne  ""){
	my $cmd="grep -f $OutputDir/all_avail_var $GroupVar > $OutputDir/groupvar.screened; echo done.";
	run_command("Screen Grouped Variables", "scr_group_var", $cmd, "$OutputDir");
}else{
	my $cmd="touch $OutputDir/groupvar.screened; echo done.";
	print STDERR "No Group Variables specified...\n";
	run_command("Screen Grouped Variables", "scr_group_var", $cmd, "$OutputDir");
}	

###############################################################################

$SummaryTable=$screened_summary_table;
$FactorFile=$screened_factor_file;
$Covariates="$OutputDir/covar.screened";
$GroupVar="$OutputDir/groupvar.screened";


run_abundance_based(
	$OutputDir,
	$SummaryTable,
	$FactorFile,
	$Covariates,
	$GroupVar,
	$AnalysisName,
	$NumALRVariables,
	$TagName
);


run_distribution_based(
	$OutputDir,
	$SummaryTable,
	$FactorFile,
	$Covariates,
	$GroupVar,
	$AnalysisName,
	$TagName
);

run_distance_based(
	$OutputDir,
	$SummaryTable,
	$FactorFile,
	$Covariates,
	$GroupVar,
	$AnalysisName,
	$TagName
);

###############################################################################

