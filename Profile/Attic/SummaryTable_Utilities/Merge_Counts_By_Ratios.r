#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"counts_file", "c", 1, "character",
	"output_file", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary file table>\n",
	"	-c <list of sample names and (original) counts>\n",
	"	-o <output file name root>\n",
	"\n",	
	"This script will read in a summary table and a list of samples with counts.\n",
	"Based on the ratio of the counts, the summary table will be collapsed into\n",
	"a single combined sample using the ratios given.\n",
	"\n",
	"The script will choose the sample with the largest (original) count, and use\n",
	"the remaining samples for down sampling relative to it, in order to preserve\n",
	"as large of a count as possible.\n",
	"\n",
	"Example counts file:\n",
	"	sample1\\t3000\n",
	"	sample2\\t4000\n",
	"	sample3\\t2000\n",
	"	...\n",
	"\n",
	"Because sample2 had the largest count, the multiplier applied to the counts\n",
	"in the summary table will be:  .75 : 1 : .5\n",
	"\n",
	"Two files will be generated:\n",
	"	*.adj_indiv.summary_table.tsv:	will contain the adjusted counts for each sample.\n",
	"	*.adj_comb.summary_table.tsv:   will contain the combined counts in a single new sample.\n",
	"\n",
	"\n");

if(	!length(opt$input_file) || 
	!length(opt$output_file) ||
	!length(opt$output_file)
){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileNameRoot=opt$output_file;
CountsFileName=opt$counts_file;

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("(Original) Counts File Name: ", CountsFileName, "\n");
cat("Output File Name Root: ", OutputFileNameRoot, "\n");       
cat("\n");

###############################################################################

###############################################################################
# Load summary file
cat("Loading: ", InputFileName, "\n", sep="");
inmat=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, 
	check.names=FALSE, comment.char="", row.names=1, quote=""));

#print(inmat);

# Grab columns we need into a vector, ignore totals, we won't trust it.
counts_mat=inmat[,2:(ncol(inmat))];

num_categories=ncol(counts_mat);
num_samples=nrow(counts_mat);
category_names=colnames(counts_mat);
sample_names=rownames(counts_mat);

cat("\n");
cat("Num Categories: ", num_categories, "\n");
cat("Num Samples: ", num_samples, "\n");
cat("\n");

#print(counts_mat);

###############################################################################
# Load counts file

cat("\n");
cat("Loading: ", CountsFileName, "\n", sep="");

ratio_counts=as.vector(read.table(CountsFileName, sep="\t",
        check.names=FALSE, comment.char="", row.names=1, quote=""));
ratio_counts_vect=ratio_counts[,1];
ratio_counts_names=rownames(ratio_counts);
names(ratio_counts_vect)=ratio_counts_names;

print(ratio_counts_vect);

largest_sample_ix=min(which(max(ratio_counts_vect)==ratio_counts_vect));
largest_sample_val=ratio_counts_vect[largest_sample_ix]
cat("\n");
cat("Largest sample is ", ratio_counts_names[largest_sample_ix], " with ", 
	largest_sample_val, " units.\n", sep="");

multiplier=ratio_counts_vect/largest_sample_val;

cat("\n");
cat("Multiplier: \n");
print(multiplier);
cat("\n");

###############################################################################

# Confirm we have counts for all the samples.

if(!setequal(sample_names, ratio_counts_names)){
	x=(setdiff(sample_names, ratio_counts_names));
	y=(setdiff(ratio_counts_names, sample_names));
	cat("Error: Names between summary table and counts file do not match up.\n");
	if(length(x)>0){
		cat("Not in ", CountsFileName, "\n");
		print(x);
	}
	if(length(y)>0){
		cat("Not in ", InputFileName, "\n");
		print(y);
	}
	quit(status=-1);	
}

###############################################################################
# Adjust counts
#print(counts_mat);

adjusted_matrix=matrix(NA, ncol=num_categories, nrow=num_samples);
rownames(adjusted_matrix)=sample_names;
colnames(adjusted_matrix)=category_names;
for(i in 1:num_samples){
	samp_name=sample_names[i];
	adjusted_matrix[samp_name,]=counts_mat[samp_name,]*multiplier[samp_name];
}

#print(adjusted_matrix);

###############################################################################

write_summary_file=function(out_mat, fname){
	fc=file(fname, "w");
	cat(file=fc, paste("sample_id\ttotal", paste(colnames(out_mat), collapse="\t"), sep="\t"));
	cat(file=fc, "\n");
	sample_names=rownames(out_mat);
	num_samples=nrow(out_mat);
	for(samp_idx in 1:num_samples){
		total=sum(out_mat[samp_idx,]);
		outline=paste(sample_names[samp_idx], total, 
			paste(out_mat[samp_idx,], collapse="\t"), sep="\t");
		cat(file=fc, outline);
		cat(file=fc, "\n");
	}
	close(fc);
}

###############################################################################

# Output Individual
write_summary_file(adjusted_matrix, paste(OutputFileNameRoot, ".adj_indiv.summary_table.tsv", sep=""));

# Output Summed
combined=matrix(apply(adjusted_matrix, 2, sum), nrow=1);
rownames(combined)=OutputFileNameRoot;
colnames(combined)=category_names;
write_summary_file(combined, paste(OutputFileNameRoot, ".adj_comb.summary_table.tsv", sep=""));

###############################################################################

cat("Done.\n")
warns=warnings();
if(length(warns)>0){
	print(warnings());
}

q(status=0)
