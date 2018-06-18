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


	The following code will be executed:

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
		open(CHK_FH, "<$checkpt_file") || die "Could not open $checkpt_file.\n";
		my $first_line=<CHK_FH>;
		close(CHK_FH);

		$first_line=chomp($first_line);
		if($first_line eq "Completed."){
			print STDERR "$cmd_name ($checkpt_file) was successfully completed previous. Skipping.\n";
			return;
		}
	}

	# Execute command
	my $clean_cmd=$cmd_str;
	print STDERR "'$cmd_str'\n";
	$clean_cmd=~s/\s+/ /g;
	print STDERR "'$clean_cmd'\n";
	my $result=`$clean_cmd`;
	my $exit_code=$?;

	# Write checkpoint file
	if($exit_code==0){
		open(CHK_FH, ">$checkpt_file") || die "Could not open $checkpt_file.\n";
		print CHK_FH "Completed.\n";
		return;
	}else{
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
	my $cmd;

	$cmd="cat $covariates $variable_list > $output_dir/cov_var";

	run_command("Concatenate variables into full model list", "make_full_model", $cmd, $output_dir);

	
	$cmd=
	"~/git/AnalysisTools/Profile/abundance_based/Fit_ALR_as_Response/Fit_ALR_as_Response.r \
                -s $summary_table \
                -f $factor_file \
                -M $output_dir/cov_var \
                -q $output_dir/cov_var \
                -p $num_alr \
                -o $output_dir/$model_name  \
                -x \";\"
	";
	run_command("Fit ALR as Response", "alr_as_resp", $cmd, $output_dir);

	$cmd=
        "~/git/AnalysisTools/Profile/abundance_based/Fit_ALR_as_Predictor/Fit_ALR_as_Predictor.r \
                -s $summary_table \
                -f $factor_file \
                -y $variable_list \
                -c $covariates \
                -q $output_dir/cov_var \
                -p $num_alr \
                -o $output_dir/$model_name \
                -x \";\"
	";
	run_command("Fit ALR as Predictor", "alr_as_pred", $cmd, $output_dir);


	$cmd=
      	"~/git/AnalysisTools/Profile/abundance_based/Compare_ALR_PredResp/Compare_ALR_PredResp.r \
                -x $output_dir/$model_name.alr_as_pred.pvals.tsv \
                -y $output_dir/$model_name.alr_as_resp.pvals.tsv \
                -u $output_dir/$model_name.alr_as_pred.coefs.tsv \
                -v $output_dir/$model_name.alr_as_resp.coefs.tsv \
                -o $output_dir/$model_name \
                -p .025
	";
	run_command("Compare ALR Pred/Resp", "alr_pred_resp_comp", $cmd, $output_dir);	

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
	my $cmd;


	$cmd="cat $covariates $variable_list > $output_dir/cov_var";
	run_command("Concatenate variables into full model list", $cmd, $output_dir);

	$cmd=
	 "~/git/AnalysisTools/Profile/distribution_based/Fit_Diversity_as_Response/Fit_Diversity_as_Response.r \
                -s $summary_table \
                -f $factor_file \
                -M $output_dir/cov_var \
                -q $output_dir/cov_var \
             	-o $output_dir/$model_name 
	";
	run_command("Fit Diversity as Response", "div_as_resp", $cmd, $output_dir);

	$cmd=
        "~/git/AnalysisTools/Profile/distribution_based/Fit_Diversity_as_Predictor/Fit_Diversity_as_Predictor.r \
                -s $summary_table \
                -f $factor_file \
                -y $variable_list \
                -c $covariates \
                -q $output_dir/cov_var \
                -o $output_dir/$model_name
	";
	run_command("Fit Diversity as Predictor", "div_as_pred", $cmd, $output_dir);

	$cmd=
        "~/git/AnalysisTools/Profile/abundance_based/Compare_ALR_PredResp/Compare_ALR_PredResp.r \
                -x $output_dir/$model_name.div_as_pred.pvals.tsv \
                -y $output_dir/$model_name.div_as_resp.pvals.tsv \
                -u $output_dir/$model_name.div_as_pred.coefs.tsv \
                -v $output_dir/$model_name.div_as_resp.coefs.tsv \
                -o $output_dir/$model_name \
                -p .025
	";
	run_command("Compare Diversity Pred/Resp", "div_pred_resp_comp", $cmd, $output_dir);


}


###############################################################################

sub run_distance_based{


	#compute distance matrix

	#permanova

	#cluster_mll

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
	"temp",
	35
);
	



















































