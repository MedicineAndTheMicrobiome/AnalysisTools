#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"reads_per_sample", "r", 2, "numeric",
	"num_bootstraps", "b", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"		-i <input summary_table.xls file>\n",
	"		[-r <Reads per sample, default=average reads/sample across all samples>]\n",
	"		[-b <Num bootstrap, default=40>]\n",
	"\n",
	"This script will read in a summary_table.xls file and then generate a series of\n",
	"rarefaction curves through bootstrapping.\n",
	"\n",
	"The output is a line for each bootstrap.\n",
	"\n",
	"\n", sep="");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileNameRoot=gsub(".summary_table.xls$", "", opt$input_file);

###############################################################################
# Load counts from file

cat("Working on ", InputFileName, "\n", sep="");

# Load summary_table.xls
mat<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, row.names=1))

# Exclude total counts column
count_mat=mat[,2:ncol(mat)]
num_cat=ncol(count_mat);

#cat("Input matrix:\n");
#print(count_mat);

#------------------------------------------------------------------------------
# Get sample names
sample_names=rownames(mat);

# Get column/category names
categories=as.vector(colnames(count_mat));
num_categories=length(categories);

#cat("Original category names:\n");
#print(categories);

#------------------------------------------------------------------------------
# Show num samples/categories to be used

NumSamples=nrow(count_mat);
NumCategories=ncol(count_mat);

cat("\n");
cat("Num Samples: ", NumSamples, "\n");
cat("Num Categories: ", NumCategories, "\n");
cat("\n");

#------------------------------------------------------------------------------
# Normalize
prob_mat=matrix(nrow=NumSamples,ncol=NumCategories);
colnames(prob_mat)=colnames(count_mat);
rownames(prob_mat)=rownames(count_mat);
sample_totals=rep(0, NumSamples);
for(i in 1:NumSamples){
    	sample_totals[i]=sum(count_mat[i,]);
	prob_mat[i,]=count_mat[i,]/sample_totals[i];
}
#cat("Normalized matrix:\n");
#print(prob_mat);

mean_sample_counts=as.integer(mean(sample_totals));
cat("\nMean counts per sample: ", mean_sample_counts, "\n", sep="");

#------------------------------------------------------------------------------
# Report rarefaction parameters

# Reads per sample
if(length(opt$reads_per_sample)){
	cat("Reads per sample specified.\n");
	ReadsPerSample=opt$reads_per_sample;
}else{
	cat("Using mean counts (default)\n");
	ReadsPerSample=mean_sample_counts;
}
cat("Reads per sample: ", ReadsPerSample, "\n");

# Num bootstraps
if(length(opt$num_bootstraps)){
	cat("Num bootstraps specified.\n");
	NumBootstraps=opt$num_bootstraps;
}else{
	cat("Using Default Num Bootstraps.\n");
	NumBootstraps=40;
}
cat("Num bootstraps: ", NumBootstraps, "\n");

###############################################################################

values_out_fname=paste(OutputFileNameRoot, ".bs_rarefaction_curves", sep="");
values_out_fh=file(values_out_fname, "wt");

for(bs_id in 1:NumBootstraps){
	cat("Generating Rarefaction Curve Instance: ", bs_id, "\n", sep="");

	# Randomly select donors, with replacement
	random_sequence=sample(1:NumSamples, NumSamples, replace=TRUE);

	rarefaction_curve=c();
	cumulative_unique=c();
	for(sample_idx in 1:NumSamples){
		# Randomly select reads from each donor
		sample_reads=sample(1:NumCategories, ReadsPerSample, prob=prob_mat[random_sequence[sample_idx],], replace=TRUE);
		cumulative_unique=unique(c(cumulative_unique, sample_reads));
		rarefaction_curve=c(rarefaction_curve, length(cumulative_unique));
	}

	cat(paste(OutputFileNameRoot, ".", bs_id, sep=""), rarefaction_curve, file=values_out_fh, sep="\t");
	cat("\n", file=values_out_fh);

}

close(values_out_fh);

###############################################################################

cat("Done.\n");
q(status=0);
