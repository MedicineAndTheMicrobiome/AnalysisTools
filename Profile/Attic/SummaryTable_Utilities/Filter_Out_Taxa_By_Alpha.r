#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"cutoff", "c", 2, "numeric",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

AbundanceCutoff=1;

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	[-c <cutoff, default = ", AbundanceCutoff, "% >]\n",
	"	[-o <output filtered summary table file name>]\n",
	"\n",	
	"This script will read in a summary file, then compute the average abundance\n",
	"across all the samples (after first normalizing).\n",
	"Taxa are then sorted by abundance and then kept piecewise until the cumulative\n",
	"abundance is greater than 1-(cutoff)\n",
	"A summary table is then generated with those taxa removed.\n",
	"\n",
	"\n", sep="");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

if(length(opt$cutoff)){
	AbundanceCutoff=opt$cutoff;
}

if(AbundanceCutoff<0 || AbundanceCutoff>=100){
	cat("Nonsensical cutoff: ", AbundanceCutoff, "\n");
	q(status=-1);
}else{
	AbundProbCutoff=AbundanceCutoff/100;
}

if(!length(opt$output_file)){
	OutputFileNameRoot = paste(gsub("\\.summary_table\\.xls", "", opt$input_file), sep="");
}else{
	OutputFileNameRoot = opt$output_file;
}

###############################################################################

InputFileName=opt$input_file;

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name Root: ", OutputFileNameRoot, "\n");       
cat("Abundance Cutoff: ", AbundProbCutoff, "\n");
cat("\n");

###############################################################################
###############################################################################

# Load data
cat("Loading Matrix...\n");
inmat=as.matrix(read.delim(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="", quote="", row.names=1))
#cat("Original Matrix:\n")
#print(inmat);

# Grab columns we need into a vector, ignore totals, we won't trust it.
counts_mat=inmat[,2:(ncol(inmat))];
#print(counts_mat);

num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);
sample_names=rownames(counts_mat);

cat("\n");
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");
cat("\n");

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

# Compute abundances across all samples
abundances_across_samples=apply(normalized, 2, mean);
names(abundances_across_samples)=colnames(counts_mat);
#print(abundances_across_samples);

sorted_abundances=sort(abundances_across_samples, decreasing=TRUE);
cumul_abund=cumsum(sorted_abundances);
#print(cumul_abund);
last_taxa_ix_to_keep=min(which(cumul_abund>(1-AbundProbCutoff)));
taxa_to_keep=sorted_abundances[1:last_taxa_ix_to_keep];
#print(taxa_to_keep);

# Subset out couts to keep
if(length(taxa_to_keep)==1){
	outmat=t(t(counts_mat[,names(taxa_to_keep)]));
	colnames(outmat)=c(taxa_to_keep);
}else{
	outmat=counts_mat[,names(taxa_to_keep)];
}

###############################################################################

# Output Filtered Summary Table
cat("Writing New Matrix...\n");
fc=file(paste(OutputFileNameRoot, ".filtered.summary_table.xls", sep=""), "w");

# Write header
cat(file=fc, paste("sample_id", "total", paste(colnames(outmat), collapse="\t"), sep="\t"));
cat(file=fc, "\n");

# Write counts
sample_names=rownames(inmat);
for(samp_idx in 1:num_samples){
	total=sum(outmat[samp_idx,]);
	cat(file=fc, paste(sample_names[samp_idx],total,paste(outmat[samp_idx,], collapse="\t"), sep="\t"));
	cat(file=fc, "\n");
}
close(fc);	

###############################################################################
cat("Writing out taxa that were eliminated...\n");

fc=file(paste(OutputFileNameRoot, ".filtered.removed_taxa.txt", sep=""), "w");
filtered_abundances=sorted_abundances[(last_taxa_ix_to_keep+1):length(sorted_abundances)];
filtered_names=names(filtered_abundances);
for(i in 1:length(filtered_abundances)){
	cat(file=fc, filtered_names[i], "\t", filtered_abundances[i], "\n");
}
close(fc);

###############################################################################

cat("Done.\n")
if(!is.null(warnings())){
	print(warnings());
}

q(status=0)
