#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	[-o <output text file root>\n",
	"\n",	
	"This script will read in the summary file and generate summary statistics\n",
	"for the category abundances across the samples.\n",
	"\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$output_file)){
	OutputFileName = paste(gsub("\\.summary_table\\.xls", "", opt$input_file), ".categorical_summary", sep="");
}else{
	OutputFileName = opt$output_file;
}

###############################################################################

InputFileName=opt$input_file;

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");       

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

sample_names=rownames(counts_mat);
category_names=colnames(counts_mat);

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
# Computations

skew=function(x){
	mu=mean(x);
	sigma=sd(x);
	sk=sum(((x-mu)/sigma)^3);
	return(sk);
}

coefvar=function(x){
	mu=mean(x);
	sigma=sd(x);
	cv=sigma/mu;
	return(cv);
}

cat_mean    =apply(normalized, 2, mean);
cat_min     =apply(normalized, 2, min);
cat_max     =apply(normalized, 2, max);
cat_median  =apply(normalized, 2, median);
cat_stdev   =apply(normalized, 2, sd);
cat_skewness=apply(normalized, 2, skew);
cat_coefvar =apply(normalized, 2, coefvar)

###############################################################################
# Output calculations

stat_names=c("Mean", "StDev", "Median", "Min", "Max", "Skew", "CoefVar");
all=cbind(cat_mean, cat_stdev, cat_median, cat_min, cat_max, cat_skewness, cat_coefvar);

sort_info=sort(cat_mean, decreasing=TRUE, index.return=TRUE);
all_sorted=all[sort_info$ix,];

cat("Writing characteristics...\n");

fc=file(OutputFileName, "w");
cat(file=fc, "CatName\t", paste(stat_names, col="\t"), "\n");
for(i in 1:num_categories){
    cat(category_names[i], "\t", paste(all_sorted[i,], col="\t"),  "\n", file=fc, sep="");
}
close(fc);

###############################################################################

writeLines("Done.\n")
if(!is.null(warnings())){
	print(warnings());
}

q(status=0)
