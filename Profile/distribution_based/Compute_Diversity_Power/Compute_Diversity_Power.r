#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');

DEF_ALPHA=0.05;
DEF_BETA=0.20;

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 2, "character",
	"alpha", "a", 2, "numeric",
	"beta", "b", 2, "numeric",
	"N", "n", 2, "numeric",
	"M", "m", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	[-o <output file root name>]\n",
	"	[-a <alpha, def=", DEF_ALPHA, "]\n",
	"	[-b <beta, def=", DEF_BETA, "]\n",
	"\n",
	"	[-n <sample size 1>]\n",
	"	[-m <sample size 2>]\n",
	"\n",
	"This script will take a summary table and then\n",
	"generate a power analysis at various effect\n",
	"sizes using the Wilcoxon Rank Sum Test statistic.\n",
	"\n",
	"This method assumes that the WRST will be used\n",
	"when performing the actual collected experimental\n",
	"data.\n",
	"\n",
	"It is assumed both sets of samples being compared\n",
	"have the same sample size.\n",
	"The eta^2 value is estimated as an ANOVA effect size\n",
	"but the actual difference in diversity is also reported\n",
	"\n",
	"Based on the specified alpha, the critical value is\n",
	"estimated based on the diversity of the samples in the\n",
	"summary table.  This is in turn used to estimate the\n",
	"gamma=1-beta value in the alternate distribution, which\n",
	"is computed by just adding the effect size to the null\n",
	"distribution.\n",
	"\n",
	"Only the sample sizes that exceed the gamma value specified\n",
	"are reported.  I.e. only experiments that are sufficiently\n",
	"powered are reported.\n",
	"\n",
	"Remember that for Cohen's Eta^2 the effect sizes are:\n",
	"	Small: < 0.02\n",
	"	Medium: ~ 0.13\n",
	"	Large: > 0.26\n",
	"\n",
	"If n and m are specified, then only those sample sizes will be computed.\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;

if(length(opt$output_file)>0){
	OutputFileRoot=opt$output_file;
}else{
	OutputFileRoot=InputFileName;
	OutputFileRoot=gsub("\\.summary_table\\.tsv$", "", OutputFileRoot);
	OutputFileRoot=gsub("\\.summary_table\\.xls$", "", OutputFileRoot);
	cat("No output file root specified.  Using input file name as root.\n");
}

Alpha=ifelse(length(opt$alpha), opt$alpha, DEF_ALPHA);
Beta=ifelse(length(opt$beta), opt$beta, DEF_BETA);

N=ifelse(length(opt$N), opt$N, 0);
M=ifelse(length(opt$M), opt$M, 0);

cat("Input Summary Table: ", InputFileName, "\n");
cat("Output Filename Root: ", OutputFileRoot, "\n");
cat("Target Alpha: ", Alpha, "\n");
cat("Target Beta: ", Beta, "\n");

if(N!=0){
	cat("N: ", N, "\n");
}
if(M!=0){
	cat("M: ", M, "\n");
}

###############################################################################

#OutputPDF = paste(OutputFileRoot, ".div_power_calc.pdf", sep="");
#cat("Output PDF file name: ", OutputPDF, "\n", sep="");
#pdf(OutputPDF,width=8.5,height=11)

###############################################################################

load_summary_file=function(fname){
        cat("Loading Summary Table: ", fname, "\n");
        inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", row.names=1))
        counts_mat=inmat[,2:(ncol(inmat))];
        return(counts_mat);
}

normalize=function(counts){
        totals=apply(counts, 1, sum);
        num_samples=nrow(counts);
        normalized=matrix(0, nrow=nrow(counts), ncol=ncol(counts));

        for(i in 1:num_samples){
                normalized[i,]=counts[i,]/totals[i];
        }

        colnames(normalized)=colnames(counts);
        rownames(normalized)=rownames(counts);
        return(normalized);
}

compute_beta=function(sample_val, effect_size, N1, N2, alpha){
	
	num_bootstraps=1000;

	# Compute NULL distribution, so we can find critical value at alpha
	null_dist=numeric(num_bootstraps);
	for(bs_ix in 1:num_bootstraps){

		x1=sample(sample_val, N1, replace=T);
		x2=sample(sample_val, N2, replace=T);

		wt_res=wilcox.test(x2, x1);
		null_dist[bs_ix]=wt_res$statistic;
	}

	# Compute Alternate distribution, so we can find beta at critical value
	alternate_dist=numeric(num_bootstraps);
	for(bs_ix in 1:num_bootstraps){

		x1=sample(sample_val, N1, replace=T);
		x2=sample(sample_val+effect_size, N2, replace=T);

		wt_res=wilcox.test(x2, x1);
		alternate_dist[bs_ix]=wt_res$statistic;
	}

	#cat("Null: mean=", mean(null_dist), " sd=",sd(null_dist), "\n");
	#cat("Altr: mean=", mean(alternate_dist), " sd=",sd(alternate_dist), "\n");

	# Find critical value at alpha
	sorted_null=sort(null_dist);
	alpha_val=sorted_null[ceiling((1-alpha)*num_bootstraps)];

	# Find beta at critical value
	sorted_alternate=sort(alternate_dist);
	beta_at_val=sum(sorted_alternate<alpha_val)/num_bootstraps;

	#cat("Sample Size:", N1, "/", N2, ", Eff Size: ", effect_size, 
	#	", W at alpha = ", alpha, " is ", alpha_val, ", at W, beta = ", beta_at_val, "\n");

	#dist_ranges=range(c(null_dist, alternate_dist));
	#hist(null_dist, xlim=dist_ranges);
	#abline(v=alpha_val, col="blue");

	#hist(alternate_dist, xlim=dist_ranges);
	#abline(v=alpha_val, col="blue");

	return(beta_at_val);

}

###############################################################################

counts_mat=load_summary_file(InputFileName);
#print(counts_mat);

###############################################################################

counts_mat_samples=rownames(counts_mat);

###############################################################################

normalized_mat=normalize(counts_mat);
#print(normalized_mat);

diversity_arr=diversity(normalized_mat, "shannon");
#print(diversity_arr);

###############################################################################

print(diversity_arr);

#par(mfrow=c(3,1));
#hist(diversity_arr, xlab="diversity", main="Distribution of Diversity");

# Determine which effect sizes to try
range=range(diversity_arr);
span=diff(range);
effect_sizes=seq(0, span/3, length.out=40);

cat("The range of diversity: ", range[1], "-", range[2], "\n");
cat("The span: ", span, "\n");
cat("\n");
cat("Computing power over: \n");
print(effect_sizes);
cat("\n");

N1_range=2:40;

if(N==0){
	N1_range=2:40;
}else{
	N1_range=N;
}

###############################################################################

out_fh=file(paste(OutputFileRoot, ".div_power.csv", sep=""), "w");
cat(file=out_fh,
	paste("eta_sqrd", "effect", "N1", "N2", "alpha", "1-beta", sep=","), "\n", sep="");

for(effect_size in effect_sizes){
	
	# Compute R^2 so we can get ANOVA effect size
	samp1=diversity_arr;
	samp2=diversity_arr+effect_size;
	comb=c(samp1, samp2);
	
	ss_residual=sum((samp1-mean(samp1))^2) + sum((samp2-mean(samp2))^2);
	ss_total=sum((comb-mean(comb))^2);
	r_sqrd=1-(ss_residual/ss_total);

	cat("Cohen's Eta^2: ", r_sqrd, "\n");
	
	for(N1 in N1_range){
		
		N2=ifelse(M==0,N1,M);

		beta_at_alpha=compute_beta(diversity_arr, effect_size, N1, N2, Alpha);

		#cat("Effect Size: ", effect_size, " N1: ", N1, "  N2: ", N2, " (1-beta): ", 1-beta_at_alpha, "\n", sep="");

		if(beta_at_alpha < Beta){
			cat("Effect Size: ", effect_size, " N1: ", N1, "  N2: ", N2, " (1-beta): ", 1-beta_at_alpha, "\n", sep="");
			cat(file=out_fh,
				paste(r_sqrd, effect_size,  N1, N2, Alpha, 1-beta_at_alpha, sep=","), "\n", sep="");
			next;
		}

	}
}

##############################################################################

close(out_fh);

cat("Done.\n")
#dev.off();
warn=warnings();
if(length(warn)){
	#print(warn);
}
q(status=0)
