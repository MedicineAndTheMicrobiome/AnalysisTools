#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"sample_size", "n", 2, "numeric",
	"num_bootstraps", "b", 2, "numeric",
	"output_file_root", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	[-n <sample size, default=number of input samples >]\n",
	"	[-b <num bootstraps/resamples, default=40>]\n",
	"	[-o <output summary table file name root, default, based on input file name>]\n",
	"\n",	
	"This script will read in a summary table and then randomly sample n (sample size) entries b times, with replacement.\n",
	"The sample ids will be appended with an index to keep them unique.\n",
	"\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

###############################################################################
# Get paramaters

InputFileName=opt$input_file;

if(!length(opt$output_file)){
	OutputFileName=gsub("\\.summary_table\\.xls$", "", opt$input_file);
}else{
	OutputFileName=opt$output_file_root;
}

if(!length(opt$num_bootstraps)){
	NumBootstraps=40;
}else{
	NumBootstraps=opt$num_bootstraps;
}

if(!length(opt$sample_size)){
	SampleSize=-1;	
}else{
	SampleSize=opt$sample_size;
}

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("Num Resamples to Produce:", NumBootstraps, "\n");

###############################################################################
###############################################################################

# Load data
inmat<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="*", row.names=1))
#cat("Original Matrix:\n")
#print(inmat);

# Grab columns we need into a vector, ignore totals, we won't trust it.
counts_mat=inmat[,2:(ncol(inmat))];
#print(counts_mat);

# Summarze what we've loaded
num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");

# Determine default resample size, if one is not specified
if(SampleSize==-1){
	SampleSize=num_samples;
}
cat("Sample Size: ", SampleSize, "\n");

###############################################################################

# Compute padding of numbers with left 0's so names line up
num_width=ceiling(log(NumBootstraps, base=10)) + 1;

# Produce NumBootstrap files each with a different resampling
for(resamples in 1:NumBootstraps){

	# Generate random sample
	subsample_indices=sample(1:num_samples, SampleSize, replace=TRUE); 
	print(subsample_indices);
	#print(counts_mat);

	# Subset out samples from original
	subsample_matrix=counts_mat[subsample_indices,];
	
	# Append sample names with index to ensure unique sample ids
	rownames(subsample_matrix)=paste(rownames(subsample_matrix),".", 1:SampleSize, sep="");
	#print(subsample_matrix);

	###############################################################################
	# Output
	fmt_str=paste("%0", num_width, "i", sep="");
	outname=paste(OutputFileName, ".", sprintf(fmt_str, resamples), ".summary_table.xls", sep="");
	cat("Writing: ", outname, "\n");

	fc=file(outname, "w");

	write(paste("sample_id", paste(colnames(inmat), collapse="\t"), sep="\t"), file=fc);
	sample_names=rownames(subsample_matrix);
	for(samp_idx in 1:SampleSize){
	    total=sum(subsample_matrix[samp_idx,]);
	    outline=paste(sample_names[samp_idx],total,paste(subsample_matrix[samp_idx,], collapse="\t"), sep="\t");
	    write(outline, file=fc);
	}
	close(fc);	
}

###############################################################################

writeLines("Done.\n")
print(warnings());

q(status=0)
