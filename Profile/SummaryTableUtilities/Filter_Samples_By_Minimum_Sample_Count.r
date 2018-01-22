#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"minimum_total", "c", 1, "numeric",
	"output_file", "o", 2, "character",
	"generate_plot", "p", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	-c <cutoff>\n",
	"	[-o <output summary table file name>\n",
	"	[-p (generate plot)]\n",
	"\n",	
	"This script will read in the summary table, and recompute the total for each sample,\n",
	"then only output the samples with total reads greater than the specified cutoff.\n",
	"\n",
	"\n");

if(!length(opt$input_file) || !length(opt$minimum_total)){
	cat(usage);
	q(status=-1);
}

if(opt$minimum_total<0){
	cat("Nonsensical negative cutoff\n");
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileName=opt$output_file;
MinimumTotalCutoff=opt$minimum_total;
GeneratePlot=opt$generate_plot;

if(length(OutputFileName)==0){
	OutputNameRoot=paste(gsub("\\summary_table\\.tsv$", "", InputFileName), ".min", MinimumTotalCutoff, sep="");
	OutputFileName=paste(OutputNameRoot, ".summary_table.tsv", sep="");
}else{
	OutputNameRoot=OutputFileName;
}

OutputPDFFileName=paste(OutputNameRoot, ".hist.pdf", sep="");

if(length(GeneratePlot)==0){
	GeneratePlot=FALSE;
}else{
	GeneratePlot=TRUE;
}


###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("Minimum Total Cutoff: ", MinimumTotalCutoff, "\n");
cat("Generate Plot: ", GeneratePlot, "\n");
cat("\n");

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

# Compute Totals
total=apply(counts_mat, 1, sum);
#cat("Total:\n");
#print(total);

# Plot Histogram
if(GeneratePlot){
	pdf(OutputPDFFileName, height=8.5, width=11);
	bins=seq(0,max(total)*1.1, MinimumTotalCutoff/4);
	hist(total, xlab="Sample Totals", main="Sample Total Distribution", breaks=bins);
	abline(v=MinimumTotalCutoff, col="blue");
}

keep_idx=total>=MinimumTotalCutoff;
num_samples_to_keep=sum(keep_idx);
cat("Number of Samples to Remove: ", num_samples-num_samples_to_keep, "\n", sep="");
cat("Number of Samples to Keep  : ", num_samples_to_keep, "\n", sep="");

# Subset out rows to keep
outmat=counts_mat[keep_idx,];


###############################################################################
# Output
cat("Writing New Matrix...\n");
fc=file(OutputFileName, "w");

write(paste("sample_id", "total", paste(colnames(outmat), collapse="\t"), sep="\t"), file=fc);
out_num_samples=nrow(outmat);
sample_names=rownames(outmat);
for(samp_idx in 1:out_num_samples){
	total=sum(outmat[samp_idx,]);
	outline=paste(sample_names[samp_idx],total,paste(outmat[samp_idx,], collapse="\t"), sep="\t");
	write(outline, file=fc);
}
close(fc);	

###############################################################################

writeLines("Done.\n")
warns=warnings();
if(length(warns)>0){
	print(warnings());
}

q(status=0)
