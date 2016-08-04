#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"sample_size", "n", 1, "numeric",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	-n <sample size>\n",
	"	[-o <output summary table file name>\n",
	"\n",	
	"This script will read in a summary table and then randomly sample n (sample size) entries.",
	"\n",
	"\n");

if(!length(opt$input_file) || !length(opt$sample_size)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$output_file)){
	outroot=gsub("\\.xls$", "", opt$input_file);
	OutputFileName = paste(outroot, ".subsample_", opt$sample_size, ".xls", sep="");
}else{
	OutputFileName=opt$output_file;
}

###############################################################################

InputFileName=opt$input_file;
SubsampleSize=opt$sample_size;

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("Subsample Size: ", SubsampleSize, "\n");

###############################################################################
###############################################################################

# Load data
inmat<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="*", row.names=1))
#cat("Original Matrix:\n")
#print(inmat);

# Grab columns we need into a vector, ignore totals, we won't trust it.
counts_mat=inmat[,2:(ncol(inmat))];
#print(counts_mat);

# Summary what we've loaded
num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");

###############################################################################

# Confirm requested sample size is smaller than population
if(SubsampleSize >= num_samples){
	cat("Error:  Subsample Size (", SubsampleSize, ") is >= Number of Samples (", num_samples, ")\n");
	q(status=0);
}

subsample_indices=sort(sample(1:num_samples, SubsampleSize));
#print(subsample_indices);
#print(counts_mat);
subsample_matrix=counts_mat[subsample_indices,];

###############################################################################
# Output
fc=file(OutputFileName, "w");

write(paste("sample_id", paste(colnames(inmat), collapse="\t"), sep="\t"), file=fc);
sample_names=rownames(subsample_matrix);
num_samples=SubsampleSize;
for(samp_idx in 1:num_samples){
	total=sum(subsample_matrix[samp_idx,]);
	outline=paste(sample_names[samp_idx],total,paste(subsample_matrix[samp_idx,], collapse="\t"), sep="\t");
	write(outline, file=fc);
}
close(fc);	

###############################################################################

writeLines("Done.\n")
print(warnings());

q(status=0)
