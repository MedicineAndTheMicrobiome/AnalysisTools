#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"num_permutations", "n", 1, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	-n <number of permutations to generate, eg. 40 for 95% CI, 200 for 99% CI.>\n",
	"\n",	
	"This script will read in a summary file table, then for each sample, permute the categories values across the samples.\n",
	"\n",
	"\n");

if(!length(opt$input_file) || !length(opt$num_permutations)){
	cat(usage);
	q(status=-1);
}

###############################################################################

InputFileName=opt$input_file;
NumPermutations=opt$num_permutations;
OutputFilenameRoot=gsub("\\.summary_table\\.xls", "", InputFileName);
OutputDir=paste(OutputFilenameRoot, ".permutations", sep="");

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output Dir: ", OutputDir, "\n");
cat("Output Filename Root: ", OutputFilenameRoot, "\n");       
cat("Num Permutations: ", NumPermutations, "\n");

cat("\n");
cat("Going to make directory: ", OutputDir, "\n", sep="");

create_success=dir.create(OutputDir);
#if(!create_success){
#	cat("Error, could not create directory.\n");
#	q(status=-1);
#}

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

# Get input dimensions
num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");
sample_names=rownames(counts_mat);

# Normalize
sums=apply(counts_mat, 1, sum);
norm_mat=matrix(0, ncol=num_categories, nrow=num_samples);
for(i in 1:num_samples){
	norm_mat[i,]=counts_mat[i,]/sums[i];
}
#print(norm_mat);

###############################################################################

idlen=log(NumPermutations,base=10)+1;
myformat=paste("%0", idlen, "i", sep="");

for(inst in 1:NumPermutations){

	# Create an output file name
	permutation_out_fname=paste(OutputDir, "/permute.", sprintf(myformat, inst), ".summary_table.xls", sep="");
	cat("Generating: ", permutation_out_fname, "\n");

	# Permute each category across all samples
	outmat=matrix(0, nrow=num_samples, ncol=num_categories);
	for(cat_id in 1:num_categories){
		outmat[,cat_id]=norm_mat[sample(1:num_samples, num_samples, replace=FALSE), cat_id];
	}

	# Write out matrix
	fc=file(permutation_out_fname, "w");
	write(paste("sample_id", "total", paste(colnames(counts_mat), collapse="\t"), sep="\t"), file=fc);
	sample_names=rownames(inmat);
	for(samp_idx in 1:num_samples){
		total=sum(outmat[samp_idx,]);
		outline=paste(sample_names[samp_idx],total,paste(outmat[samp_idx,], collapse="\t"), sep="\t");
		write(outline, file=fc);
	}
	close(fc);	

}


q(status=1);

###############################################################################

###############################################################################

writeLines("Done.\n")
if(!is.null(warnings())){
	print(warnings());
}

q(status=0)
