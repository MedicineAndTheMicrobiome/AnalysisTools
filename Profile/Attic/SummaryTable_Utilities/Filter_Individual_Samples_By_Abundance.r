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
	"	-c <cutoff, eg. 1, if proportion of categories for a sample is less than 1, zero it out.\n",
	"	[-o <output summary table file name>\n",
	"\n",	
	"This script will read in a summary file, and zero out categories where the in-sample proportion\n",
	"is less than the specified cutoff.\n",
	"\n",
	"If a category becomes all zeros across all samples, that category/column is removed.\n",
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
	OutputFileName = paste(gsub("\\.summary_table\\.xls", "", opt$input_file), ".filtered.", sprintf("%02f",opt$cutoff), ".summary_table.xls", sep="");
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

for(sample_idx in 1:num_samples){
	lt_cutoff_idx=normalized[sample_idx,]<Cutoff;
	counts_mat[sample_idx,lt_cutoff_idx]=0;
	normalized[sample_idx,lt_cutoff_idx]=0;
}

#print(normalized);
#print(counts_mat);

category_sums=apply(counts_mat,2,sum);
rem_col=category_sums==0;

outmat=counts_mat[,!rem_col];
#print(outmat);

for(sample_idx in 1:num_samples){
	if(sum(counts_mat[sample_idx,])==0){
		cat("Warning!  After filtering, ", sample_names[sample_idx], " has zero sum.\n", sep="");
	}
}

# Normalize samples so they add back up to one.
total=apply(outmat, 1, sum);
for(sample_idx in 1:num_samples){
	outmat[sample_idx,]=outmat[sample_idx,]/total[sample_idx];
}
outmat=outmat[total>0,];

###############################################################################
# Output
cat("Writing New Matrix...\n");
fc=file(OutputFileName, "w");

write(paste("sample_id", "total", paste(colnames(outmat), collapse="\t"), sep="\t"), file=fc);
sample_names=rownames(inmat);
for(samp_idx in 1:num_samples){
	total=sum(outmat[samp_idx,]);
	outline=paste(sample_names[samp_idx],total,paste(outmat[samp_idx,], collapse="\t"), sep="\t");
	write(outline, file=fc);
}
close(fc);	

###############################################################################

writeLines("Done.\n")
if(!is.null(warnings())){
	print(warnings());
}

q(status=0)
