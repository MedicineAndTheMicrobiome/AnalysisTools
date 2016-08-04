#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"cutoff", "c", 1, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"		-i <input summary_table.xls file>\n",
	"		-c <cutoff, ie. .5, or .95, to determine median or 95% upper bound>\n",
	"\n",
	"This script will read in a summary table which contains a set of replicates.\n",
	"\n",
	"For each group of replicates, a chi-squared test of homogeneity is used to\n",
	"determine the 'effective' number of samples that must have been the basis for\n",
	"amplification.\n",
	"\n",
	"The model assumes that a smaller set of DNA was sampled during an early stage\n",
	"of PCR, and that the amplification process only amplified the original sampled proportions.\n",
	"\n",
	"For example, if the proportion of organisms in a sample was 4:3:2:1, a small sample of N=10 may\n",
	"yield a proportion of 3:4:1:2.  This would not be rejected because it is not significantly\n",
	"different than the source distribution.  However if N=100, and your proportion of organisms\n",
	"was 30:40:10:20, then this would be rejected using a chi-squared test.  If you believe the\n",
	"sample replicates should follow a chi-squared difference and an amplification step was used\n",
	"to magnify the original N, then this script will calculate the largest pre-amplification sample\n",
	"size that still is not rejected by chi-squared.\n",
	"\n",
	"Using a cutoff of .5, suggests that you want your input replicates to be close to the median\n",
	"outcome of your effective sample size.  Using a value of 95%, suggests that your replicate\n",
	"is closer to the upperbound of sample dissimilarity within the calculated sample size.\n",
	"\n");

if(!length(opt$input_file) || !length(opt$cutoff)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileNameRoot=gsub("\\.xls$", "", opt$input_file);
Cutoff=opt$cutoff;

###############################################################################
# Load counts from file

cat("Working on ", InputFileName, "\n", sep="");

# Load summary_table.xls
mat<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE))

# Exclude total counts column and convert from string to integer
count_mat=apply((mat[,3:ncol(mat)]), 2, as.integer);
group_names=mat[,1];

#cat("Counts:\n");
#print(count_mat);

#cat("Group Names:\n");
#print(group_names);

#------------------------------------------------------------------------------

NumSamples=nrow(count_mat);
NumCategories=ncol(count_mat);

cat("\n");
cat("Num Samples: ", NumSamples, "\n");
cat("Num Categories: ", NumCategories, "\n");
cat("\n");

#------------------------------------------------------------------------------
# Normalize
prob_mat=matrix(nrow=NumSamples,ncol=NumCategories);
sample_totals=rep(0, NumSamples);
for(i in 1:NumSamples){
    	sample_totals[i]=sum(count_mat[i,]);
	prob_mat[i,]=count_mat[i,]/sample_totals[i];
}

#cat("Normalized matrix:\n");
#print(prob_mat);

###############################################################################

unique_groups=unique(group_names);
num_unique_groups=length(unique_groups);

cat("Num replicates:", num_unique_groups, "\n");

###############################################################################

compute_chisq=function(param, data){
# This function will compute the p-value for the null hypothesis that a set of 
# samples are are from the same distribution.

	sample_size=param;
	mat_info=data;

	# Compute degrees of freedom
	nrow=nrow(mat_info$group_mat);
	ncol=ncol(mat_info$group_mat);
	df=(nrow-1)*(ncol-1);	

	# Compute expected and observed for the stated sample size
	exp=sample_size*mat_info$centroid;
	obs=sample_size*mat_info$group_mat;

	#cat("N = ", sample_size, " df=", df, "\n", sep="");
	#cat("Obs:\n");
	#print(obs);
	#cat("Exp:\n");
	#print(exp);

	# Sum up the chi-squared components
	chisq_stat=0;
	for(obsix in 1:nrow){
		chisq_stat = chisq_stat + sum(((obs[obsix,]-exp)^2)/exp);
	}

	# Get probability of not rejecting this set assuming N
	p=pchisq(chisq_stat, df);
	#cat("n=", sample_size, " Prob(same) ", p, " chisq:", chisq_stat, "\n");
	return(p);
}

#------------------------------------------------------------------------------

comp_eff_SampSize=function(group_mat, centroid, cutoff, max_reads){
# Cycles through effective reads to determine the size necessary to acheive the
# specified cutoff.

	mat_info=list();
	mat_info$group_mat=group_mat;
	mat_info$centroid=centroid;

	p=0;
	effective_reads=1;
	while(p<cutoff){
		p=compute_chisq(effective_reads,mat_info);
		effective_reads = effective_reads + 1;
	}
	
	return(effective_reads);

}

#------------------------------------------------------------------------------

remove_zero_categories=function(mat){
# Remove columns which are consistently zero.
	sum=apply(mat, 2, sum);
	nz=sum>0;
	return(mat[,nz]);
}

#------------------------------------------------------------------------------

effective_sample_size=rep(0, num_unique_groups);

for(i in 1:num_unique_groups){
	cat("Working on: ", unique_groups[i], "\n");
	group_idx=(which(group_names==unique_groups[i]));
	if(length(group_idx)>1){
		nz_mat=remove_zero_categories(prob_mat[group_idx,]);
    		max_reads=max(sample_totals[group_idx]);
		centroid=apply(nz_mat, 2, mean);
		effective_sample_size[i]=comp_eff_SampSize(nz_mat, centroid, Cutoff, max_reads);
	}else{
		cat("\tNo replicate for sample.\n");
	}
	cat("\n");
}

###############################################################################
# Output Results

results_fname=paste(OutputFileNameRoot, ".eff_samp_size.", Cutoff, ".tsv", sep="");
results_fh=file(results_fname, "wt");
for(i in 1:num_unique_groups){
	cat(unique_groups[i], effective_sample_size[i], sep="\t", file=results_fh);
	cat("\n", file=results_fh);
}
close(results_fh);

###############################################################################

cat("Done.\n");
q(status=0);
