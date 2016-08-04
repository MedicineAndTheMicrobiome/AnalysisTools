#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_filename_root", "o", 2, "character"
);

NUM_SUMMARIZED=5;

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	[-o <output filename root>]\n",
	"\n",	
	"This script will compute:\n",
	"	1.) The mean proportion and the 95% confidence interval around the mean\n",
	"		This uses the t-distribution, with (num_samples - 1) degrees of freedom\n",
	"	2.) The 95% prediction interval around the proportion\n",
	"		This models the porportions with a beta distribution using maximum likelihood.\n",
	"\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

###############################################################################

InputFileName=opt$input_file;

if(length(opt$output_filename_root)>0){
	OutputFileNameRoot=opt$output_filename_root;
}else{
	OutputFileNameRoot=gsub(".summary_table.xls", "", opt$input_file);
}

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name Root: ", OutputFileNameRoot, "\n");
cat("\n");

pdf(paste(OutputFileNameRoot, ".beta_fit.pdf", sep=""), height=11, width=8.5);

###############################################################################
###############################################################################

# Load data
inmat=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="*", row.names=1))
#cat("Original Matrix:\n")
#print(inmat);

# Grab columns we need into a vector, ignore totals, we won't trust it.
counts_mat=inmat[,2:(ncol(inmat))];
#print(counts_mat);

num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);

cat("Num Samples:", num_samples, "\n");
cat("Num Categories:", num_categories, "\n");

category_names=colnames(counts_mat);
sample_names=rownames(counts_mat);

###############################################################################
# Normalize

# Sum sample totals
sample_totals=apply(counts_mat, 1, sum);
#print(sample_totals);

# normalize, to compute probabilities
normalized=matrix(0, nrow=num_samples, ncol=num_categories);
for(i in 1:num_samples){
	normalized[i,]=counts_mat[i,]/sample_totals[i];
}
#cat("\nNormalized:\n");
#print(normalized);

prop_means=apply(normalized, 2, mean);
prop_var=apply(normalized, 2, var);

#cat("\nMeans:\n");
#print(prop_means);

#cat("\nVar:\n");
#print(prop_var);

cat("\n");

###############################################################################
# Compute the standard error and 95% confidence intervals for the mean

lb=numeric(num_categories);
ub=numeric(num_categories);

alpha=.05

T=qt(1-alpha/2, num_samples-1);

cat("T at ", alpha, "/2 = ", T, "\n", sep="");
cat("Num samples: ", num_samples, "\n", sep="");

for(i in 1:num_categories){
	bounds=sqrt(prop_var[i]/num_samples)*T;
	lb[i]=prop_means[i]-bounds;	
	ub[i]=prop_means[i]+bounds;	

	if(lb[i]<0){
		lb[i]=0;
	}
	if(ub[i]>1){
		ub[i]=1;
	}
}

###############################################################################
# Compute the prediction intervals

params=list(num_categories);
par(mfrow=c(3,1));

beta_mean=numeric(num_categories);
beta_lb=numeric(num_categories);
beta_ub=numeric(num_categories);
beta_num_samples=numeric(num_categories);

for(i in 1:num_categories){
	#print(normalized[,i]);
	#cat("\n");

	# Extract that that is non zero.  Here we are assuming 0's are not zero, but just not see because
	#   it was beyond the assay's sensitivity.
	gtzero=normalized[,i]>0;
	data=normalized[gtzero,i];
	beta_num_samples[i]=length(data);

	if(length(data)>3){
		# Fit to the beta distribution when there are more than 3 samples
		param=fitdistr(data,"beta", list(shape1=.5, shape2=.5), lower=c(0.0001, 0.0001));
		params[[i]]=param;
		maxplot=max(data)*2;

		# Compute the mean and prediction intervals based on the fitted beta distribution
		beta_mean[i]=param$estimate[1]/(param$estimate[1]+param$estimate[2]);
		beta_lb[i]=qbeta(alpha/2, param$estimate[1], param$estimate[2]);
		beta_ub[i]=qbeta(1-alpha/2, param$estimate[1], param$estimate[2]);

		# Get the distribution for the beta distribution based on the ML estimated parameters
		overlap_pts=seq(0,maxplot, length.out=80);
		beta_density=dbeta(overlap_pts, param$estimate[1], param$estimate[2]);

		# Generate the histogram and overlay the theoretical beta distribution over it
		h=hist(data, plot=F);
		hist(data, freq=F,
			xlim=c(0,maxplot), ylim=c(0, max(beta_density[beta_density!=Inf], h$density)*1.2), 
			main=category_names[i], xlab="Abundance", ylab="Density",
			);
		
		# Label the mean and bounds
		mtext(paste("mean=", sprintf("%3.4f", beta_mean[i]), sep=""), line=.5, cex=.8);
		mtext(paste("(", sprintf("%3.4f", lb[i]), ", ", sprintf("%3.4f",ub[i]), ")", sep=""), line=-.5, cex=.8);
		
		# Mark the bounds
		abline(v=beta_lb[i], lty=2, col="grey");
		abline(v=beta_ub[i], lty=2, col="grey");
		points(overlap_pts, beta_density, type="l", col="blue");
	}else{
		beta_mean[i]=NA;
		beta_lb[i]=NA;
		beta_ub[i]=NA;
	}
}

###############################################################################


fh=file(paste(OutputFileNameRoot, ".confidence_intervals.csv", sep=""), "w");

cat(file=fh, paste(c(
	"#Taxa", "Est. Mean", "LowerBound", "UpperBound", "(Beta)", "Mean", "LowerBound", "UpperBound", "NumSamplesUsed"
	), collapse=","),
	"\n");

for(i in 1:num_categories){
	outline=paste(category_names[i], prop_means[i], lb[i], ub[i],"", beta_mean[i], beta_lb[i], beta_ub[i], beta_num_samples[i], sep=",");
	cat(file=fh, outline, "\n", sep="");

}

###############################################################################

cat("Done.\n")
print(warnings());

q(status=0)
