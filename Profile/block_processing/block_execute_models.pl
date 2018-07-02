#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Temp;
use vars qw ($opt_s $opt_f $opt_c $opt_g $opt_o);
use File::Basename;
use Cwd;
use Digest::MD5;
use Sys::Hostname;

getopts("s:f:c:g:o:");

my $usage = "
	usage:
	$0

	-s <summary table>
	-f <factor file>
	-c <covariates list>
	-g <\"grouped\" variables list>
	
	-o <output directory>

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



";

if(
	!defined($opt_s) || 
	!defined($opt_f) || 
	!defined($opt_c) || 
	!defined($opt_g) || 
	!defined($opt_o)
){
	die $usage;
}


my $SummaryTable=$opt_s;
my $FactorFile=$opt_f;
my $Covariates=$opt_c;
my $GroupVar=$opt_g;
my $OutputDir=$opt_o;
my $AnalysisName=$GroupVar;

$AnalysisName=~s/\.txt$//;
$AnalysisName=File::Basename::fileparse($AnalysisName);

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
print STDERR "Summary Table:        $SummaryTable\n";
print STDERR "Factor File:          $FactorFile\n";
print STDERR "Covariates List:      $Covariates\n";
print STDERR "Group Variables List: $GroupVar\n";
print STDERR "Output Directory:     $OutputDir\n";
print STDERR "Analysis Name:        $AnalysisName\n";
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
		print CHK_FH "\nFailed.\n";
		close(CHK_FH);

		print STDERR "Error: $cmd_name returned with non-zero error code.\n";
		print STDERR "\n";
		print STDERR "$cmd_str\n";
		exit;
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
	
	$cmd=
	"~/git/AnalysisTools/Profile/abundance_based/Fit_ALR_as_Response/Fit_ALR_as_Response.r \
                -s $summary_table \
                -f $factor_file \
                -M $output_dir/cov_var \
                -q $output_dir/cov_var \
                -p $num_alr \
                -o $output_dir/abundance/$RESP_OUT_DIR/$model_name  \
                -x \";\"
	";
	run_command("Fit ALR as Response", "alr_as_resp", $cmd, "$output_dir/abundance/$RESP_OUT_DIR");

	$cmd=
        "~/git/AnalysisTools/Profile/abundance_based/Fit_ALR_as_Predictor/Fit_ALR_as_Predictor.r \
                -s $summary_table \
                -f $factor_file \
                -y $variable_list \
                -c $covariates \
                -q $output_dir/cov_var \
                -p $num_alr \
                -o $output_dir/abundance/$PRED_OUT_DIR/$model_name \
                -x \";\"
	";
	run_command("Fit ALR as Predictor", "alr_as_pred", $cmd, "$output_dir/abundance/$PRED_OUT_DIR");


	$cmd=
      	"~/git/AnalysisTools/Profile/abundance_based/Compare_ALR_PredResp/Compare_ALR_PredResp.r \
                -x $output_dir/abundance/$PRED_OUT_DIR/$model_name.alr_as_pred.pvals.tsv \
                -y $output_dir/abundance/$RESP_OUT_DIR/$model_name.alr_as_resp.pvals.tsv \
                -u $output_dir/abundance/$PRED_OUT_DIR/$model_name.alr_as_pred.coefs.tsv \
                -v $output_dir/abundance/$RESP_OUT_DIR/$model_name.alr_as_resp.coefs.tsv \
                -o $output_dir/abundance/$COMP_DIR/$model_name.alr \
                -p .025
	";
	run_command("Compare ALR Pred/Resp", "alr_pred_resp_comp", $cmd, "$output_dir/abundance/$COMP_DIR");	

}

###############################################################################

sub run_distribution_based{

	my $output_dir=shift;
	my $summary_table=shift;
	my $factor_file=shift;
	my $covariates=shift;
	my $variable_list=shift;
	my $model_name=shift;

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


	#######################################################################
	# Diversity
	$cmd=
	 "~/git/AnalysisTools/Profile/distribution_based/Fit_Diversity_as_Response/Fit_Diversity_as_Response.r \
                -s $summary_table \
                -f $factor_file \
                -M $output_dir/cov_var \
                -q $output_dir/cov_var \
             	-o $output_dir/distribution/$RESP_OUT_DIR/$model_name 
	";
	run_command("Fit Diversity as Response", "div_as_resp", $cmd, "$output_dir/distribution/$RESP_OUT_DIR");

	$cmd=
        "~/git/AnalysisTools/Profile/distribution_based/Fit_Diversity_as_Predictor/Fit_Diversity_as_Predictor.r \
                -s $summary_table \
                -f $factor_file \
                -y $variable_list \
                -c $covariates \
                -q $output_dir/cov_var \
                -o $output_dir/distribution/$PRED_OUT_DIR/$model_name
	";
	run_command("Fit Diversity as Predictor", "div_as_pred", $cmd, "$output_dir/distribution/$PRED_OUT_DIR");

	$cmd=
        "~/git/AnalysisTools/Profile/abundance_based/Compare_ALR_PredResp/Compare_ALR_PredResp.r \
                -x $output_dir/distribution/$PRED_OUT_DIR/$model_name.div_as_pred.pvals.tsv \
                -y $output_dir/distribution/$RESP_OUT_DIR/$model_name.div_as_resp.pvals.tsv \
                -u $output_dir/distribution/$PRED_OUT_DIR/$model_name.div_as_pred.coefs.tsv \
                -v $output_dir/distribution/$RESP_OUT_DIR/$model_name.div_as_resp.coefs.tsv \
                -o $output_dir/distribution/$COMP_DIR/$model_name.div \
                -p .025
	";
	run_command("Compare Diversity Pred/Resp", "div_pred_resp_comp", $cmd, "$output_dir/distribution/$COMP_DIR");

	#######################################################################
	# Plots
	$cmd=
	"~/git/AnalysisTools/Profile/distribution_based/Plot_StackedBar/Plot_StackedBar.r \
		-i $summary_table \
		-f $factor_file \
		-M $output_dir/cov_var \
		-o $output_dir/distribution/$STCK_BAR_DIR/$model_name \
		-s \";\" \
	";
	run_command("Plot Stacked Bar Plots", "stacked_bp", $cmd, "$output_dir/distribution/$STCK_BAR_DIR");

	$cmd=
	"~/git/AnalysisTools/Profile/distribution_based/Plot_RankAbundance_wFactors/Plot_RankAbundance_wFactors.r \
		-i $summary_table \
		-f $factor_file \
		-M $output_dir/cov_var \
		-o $output_dir/distribution/$RANK_ABND_DIR/$model_name \
		-s \";\" \
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

	my $DISTMAT_DIR="dist_mat";
	my $PERMA_DIR="permanova";
	my $CLUST_MLL_DIR="clust_mll";
	my $cmd;
	my $dist_type="man";


	$cmd="cat $covariates $variable_list > $output_dir/cov_var";
	run_command("Concatenate variables into full model list", "concat", $cmd, $output_dir);

	mkdir "$output_dir/distance";
	mkdir "$output_dir/distance/$DISTMAT_DIR";
	mkdir "$output_dir/distance/$PERMA_DIR";
	mkdir "$output_dir/distance/$CLUST_MLL_DIR";

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
		-q $output_dir/cov_var
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
		-d $dist_type
	";
	run_command("Fit Multinomial Log-Linear Model", "clust_mll", $cmd, "$output_dir/distance/$CLUST_MLL_DIR");

	$cmd=
	"~/git/AnalysisTools/Profile/distance_based/Plot_ClusterHeatmap/Plot_ClusterHeatmap.r \
		-i $summary_table \
		-f $factor_file \
		-d $dist_type \
		-M $output_dir/cov_var \
		-c $output_dir/distance/$CLUST_MLL_DIR/$model_name\.$dist_type\.cl_mll.used_samp.lst \
		-o $output_dir/distance/$CLUST_MLL_DIR/$model_name
	";
	run_command("Plot Metadata Heatmap", "cl_hmp", $cmd, "$output_dir/distance/$CLUST_MLL_DIR");

	$cmd=
	"~/git/AnalysisTools/Profile/distance_based/CH_Cluster_Stop/CH_Cluster_Stop.r \
		-i $summary_table \
		-d $dist_type \
		-c $output_dir/distance/$CLUST_MLL_DIR/$model_name\.$dist_type\.cl_mll.used_samp.lst \
		-o $output_dir/distance/$CLUST_MLL_DIR/$model_name
	";
	run_command("Compute Calinski-Harabasz Stopping", "ch_stop", $cmd, "$output_dir/distance/$CLUST_MLL_DIR");

	$cmd=
	"~/git/AnalysisTools/Profile/distance_based/Cluster_Influencers/Cluster_Influencers.r \
		-i $summary_table \
		-d $dist_type \
		-l $output_dir/distance/$CLUST_MLL_DIR/$model_name\.$dist_type\.cl_mll.used_samp.lst \
		-o $output_dir/distance/$CLUST_MLL_DIR/$model_name
	";
	run_command("Compute Cluster Influencers", "cl_inf", $cmd, "$output_dir/distance/$CLUST_MLL_DIR");


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
	$Covariates,
	$GroupVar,
	$AnalysisName,
	35
);


run_distribution_based(
	$OutputDir,
	$SummaryTable,
	$FactorFile,
	$Covariates,
	$GroupVar,
	$AnalysisName
);

run_distance_based(
	$OutputDir,
	$SummaryTable,
	$FactorFile,
	$Covariates,
	$GroupVar,
	$AnalysisName
);




















































