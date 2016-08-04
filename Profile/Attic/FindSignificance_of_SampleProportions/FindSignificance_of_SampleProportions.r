#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 2, "character",
	"sample_id", "s", 1, "character",
	"normalize", "n", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	-s <sample id of interest>\n",
	"	[-o <output values>]\n",
	"	[-n (flag whether to normalize values first)]\n",
	"\n",	
	"This script will load the summary table file and determine the rank\n",
	"and p-value of the sample of interest, if you assume the other samples\n",
	"are the null distribution.  This is a one tailed hypothesis test.\n",
	"\n",
	"Note that the p-values may be a little off if there are a lot of ties.\n",
	"\n",
	"\n");

if(!length(opt$input_file) || !length(opt$sample_id)){
	cat(usage);
	q(status=-1);
}

Normalize=FALSE;
if(length(opt$normalize)){
	Normalize=TRUE;
}

# Figure out output filename
if(!length(opt$output_file)){
	OutputFileName = paste(gsub("\\.summary_table\\.xls", "", opt$input_file), ".rank_and_pvalues.txt", sep="");
}else{
	OutputFileName=opt$output_file;
}

###############################################################################

InputFileName=opt$input_file;
SampleOfInterest=opt$sample_id;

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("Sample of Interest: ", SampleOfInterest, "\n");

if(Normalize){
	cat("Normalization is ON\n");
}else{
	cat("Normalization is OFF\n");
}

###############################################################################
###############################################################################

# Load data
cat("Loading Matrix...\n");
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
# Normalize?

if(Normalize){
	# Sum sample totals
	sample_totals=numeric();
	for(i in 1:num_samples){
		sample_totals[i]=sum(counts_mat[i,]);
	}

	# normalize, to compute probabilities
	normalized=matrix(0, nrow=num_samples, ncol=num_categories);
	for(i in 1:num_samples){
		normalized[i,]=counts_mat[i,]/sample_totals[i];
	}

	rownames(normalized)=rownames(counts_mat);
	colnames(normalized)=colnames(counts_mat);
	data_mat=normalized;
}else{
	data_mat=counts_mat;
}

#print(data_mat);

###############################################################################
# Find sample of interest
sample_names=rownames(data_mat);
sample_of_interest_idx=which(sample_names==SampleOfInterest);
if(length(sample_of_interest_idx)>2){
	cat("Error: multiple samples match your sample of interest: ", SampleOfInterest, "\n");
	q(status=-1);
}
cat("Index of Sample of Interest: ", sample_of_interest_idx, "\n");

###############################################################################
# Compute ranks and p-values
soi_ranks=numeric();
pvalues=numeric();
for(i in 1:num_categories){
	ranks=rank(data_mat[,i]);
	soi_rank=ranks[sample_of_interest_idx];

	max_ranks=max(ranks);
	if(soi_rank==max_ranks){
		pvalue=paste("< ", 1/soi_rank);
	}else{
		pvalue=(max_ranks-soi_rank)/max_ranks;
	}

	soi_ranks[i]=soi_rank;
	pvalues[i]=pvalue;

}

###############################################################################

fc=file(OutputFileName, "w");
cat_names=colnames(data_mat);

for(i in 1:num_categories){
	cat(file=fc, cat_names[i], "\t", soi_ranks[i], "\t", pvalues[i], "\n");
}

close(fc);	

###############################################################################

print(warnings());

cat("Done.\n");
q(status=0)
