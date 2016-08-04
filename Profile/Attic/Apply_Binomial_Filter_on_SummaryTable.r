#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"cutoff", "c", 1, "numeric",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	-c <cutoff, eg. 99.9, means count will be non-zero 99.9% of time upon resampling.>\n",
	"	[-o <output summary table file name>\n",
	"\n",	
	"This script will read in a summary_table and filter out all the counts\n",
	"that are not statistically significantly non-zero.  The code uses\n",
	"the binomial distribution to determine the probability a category\n",
	"will be 0 upon resampling.  If the probability is greater than a cutoff\n",
	"eg. .1%, then the count is set to zero.  In other words, if the\n",
	"count will be nonzero 99.9% of the time, then count will not be eliminated.\n",
	"\n",
	"\n");

if(!length(opt$input_file) || !length(opt$cutoff)){
	cat(usage);
	q(status=-1);
}

if(opt$cutoff<0 || opt$cutoff>=100){
	cat("Nonsensical cutoff: ", opt$cutoff, "\n");
	q(status=-1);
}

if(!length(opt$output_file)){
	OutputFileName = paste(opt$input_file, ".bin_filt.", opt$cutoff, sep="");
}else{
	OutputFileName=opt$output_file;
}

###############################################################################

InputFileName=opt$input_file;
Cutoff=opt$cutoff/100.0;

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("Cutoff: ", Cutoff, "\n");

###############################################################################
###############################################################################

# Load data
inmat<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="*", row.names=1))
#cat("Original Matrix:\n")
#print(inmat);

# Grab columns we need into a vector, ignore totals, we won't trust it.
counts_mat=inmat[,2:(ncol(inmat))];
#print(counts_mat);

num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");

###############################################################################
# Normalize

# Sum sample totals
sample_totals=numeric();
for(i in 1:num_samples){
	sample_totals[i]=sum(counts_mat[i,]);
}
#print(sample_totals);

# normalize, to compute probabilities
normalized=matrix(0, nrow=num_samples, ncol=num_categories);
for(i in 1:num_samples){
	normalized[i,]=counts_mat[i,]/sample_totals[i];
}
#print(normalized);

###############################################################################
# Use binomial distribution to determine probability count will not be zero when resampled
prob_cutoff_mat=matrix(logical(), nrow=num_samples, ncol=num_categories);
for(samp_idx in 1:num_samples){
	for(cat_idx in 1:num_categories){
		prob_cutoff=pbinom(0, sample_totals[samp_idx], normalized[samp_idx, cat_idx]);
		prob_cutoff_mat[samp_idx, cat_idx]=prob_cutoff<(1-Cutoff);
	}
}
# Probability of being 0 more than CUTOFF percent of the time
#print(prob_cutoff_mat);

# Filtered results
counts_mat[!prob_cutoff_mat]=0;
#print(counts_mat);

###############################################################################
# Output
fc=file(OutputFileName, "w");

write(paste("sample_id", paste(colnames(inmat), collapse="\t"), sep="\t"), file=fc);
sample_names=rownames(inmat);
for(samp_idx in 1:num_samples){
	total=sum(counts_mat[samp_idx,]);
	outline=paste(sample_names[samp_idx],total,paste(counts_mat[samp_idx,], collapse="\t"), sep="\t");
	write(outline, file=fc);
}
close(fc);	

###############################################################################

writeLines("Done.\n")
print(warnings());

q(status=0)
