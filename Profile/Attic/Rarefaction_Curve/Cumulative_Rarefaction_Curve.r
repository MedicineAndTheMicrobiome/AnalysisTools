#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"reads_per_sample", "r", 2, "numeric",
	"num_samples", "d", 2, "numeric",
	"num_bootstraps", "b", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n\t\t<Input summary_table.xls FileName>\n\n",
	"	-i <input summary_table.xls file>\n",
	"	[-r <Reads per sample>]\n",
	"	[-d <Num samples to graph>]\n",
	"	[-b <Num bootstrap, default=40>]\n",
	"\n",
	"This script will generate a cumulative rarefaction curve based on the input summary table.\n",
	"\n",
	"The 'Reads per sample' is an option to specify how many reads there are per sample.\n",
	"  If this is not specified, the 'Reads per sample' is estimated by taking the average\n",
	"   number of reads across all samples in the input summary table.\n",
	"\n",
	"The 'Num samples to graph' is an option to specify how many samples to graph along the\n",
	"  the x axis.  If this is not specified, the number of samples to graph will be the number\n",
	"  of samples in the input summary table.\n",
	"\n",
	"By default, the confidence intervals (CI) will be 95% (40 bootstraps per number of samples)\n",
	"If you ask for more or less bootstraps, your CI will be (num bootstraps-2)/(num boostraps).\n",
	"\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileNameRoot=opt$input_file;

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
# Compute max possible taxa

count_mat=apply(count_mat, 2, sum);
num_nonzero_taxa=sum(count_mat>0);
cat("Num Nonzero Taxa:", num_nonzero_taxa, "\n");

#------------------------------------------------------------------------------
# Report rarefaction parameters

# Reads per sample
if(length(opt$reads_per_sample)){
	cat("Reads per sample specified.\n");
	ReadsPerSample=opt$reads_per_sample;
}else{
	ReadsPerSample=mean_sample_counts;
}
cat("Reads per sample: ", ReadsPerSample, "\n");

# Num samples
if(length(opt$num_samples)){
	cat("Num samples specified.\n");
	BS_NumSamples=opt$num_samples;
}else{
	BS_NumSamples=NumSamples;
}
cat("Num samples: ", NumSamples, "\n");

# Num bootstraps
if(length(opt$num_bootstraps)){
	cat("Num bootstraps specified.\n");
	NumBootstraps=opt$num_bootstraps;
}else{
	NumBootstraps=40;
}
cat("Num bootstraps: ", NumBootstraps, "\n");
cat("\n");

###############################################################################

cumulative_rarefaction=function(prob_mat, ReadsPerSample, NumSamples){
# since we are trying to perform cumulative rarefaction, we need to accumulate the donors as well
	
	# Keep track of cumulative discovered taxa
	cumulative_counts=rep(0, NumSamples);

	# Determine available donors and taxa/categories
	num_available_donors=nrow(prob_mat);
	num_categories=ncol(prob_mat);

	# Generate a random ordering of donors
	donor_ordering=sample(1:num_available_donors, NumSamples, replace=TRUE);
	
	#cat("Donor Ordering:\n");
	#print(donor_ordering);

	# Sample reads
	unique_categories=integer(0);
	for(donor_idx in 1:NumSamples){
		cur_donor=donor_ordering[donor_idx];
		sample_reads=sample(1:num_categories, ReadsPerSample, prob=prob_mat[cur_donor,], replace=TRUE);
		unique_categories=unique(c(sample_reads,unique_categories));
		cumulative_counts[donor_idx]=length(unique_categories);
	}

	#cat("Cumulative Counts:\n");
	#print(cumulative_counts);
	return(cumulative_counts);	
}

###############################################################################

cumu_taxa_mat=matrix(0, nrow=NumBootstraps, ncol=BS_NumSamples);
for(bs in 1:NumBootstraps){
	cumu_taxa_mat[bs,]=cumulative_rarefaction(prob_mat, ReadsPerSample, BS_NumSamples);
}

# Compute CI
raref_lb=rep(0,BS_NumSamples);
raref_ub=rep(0,BS_NumSamples);
raref_med=rep(0,BS_NumSamples);

for(bs_num_samples in 1:BS_NumSamples){
	sorted=sort(cumu_taxa_mat[,bs_num_samples]);
	raref_lb[bs_num_samples]=sorted[2];
	raref_ub[bs_num_samples]=sorted[NumBootstraps-1];
	raref_med[bs_num_samples]=median(sorted);
}

cat("LB:\n");
print(raref_lb);
cat("UP:\n");
print(raref_ub);
cat("Median:\n");
print(raref_med);

ci=100*(NumBootstraps-2)/NumBootstraps;

###############################################################################

pdf(paste(OutputFileNameRoot, ".cumulative_rarefaction.pdf", sep=""), height=8.5, width=11);

library(gplots);

plotCI(x=1:BS_NumSamples, y=raref_med, ui=raref_ub, li=raref_lb, ylim=c(0,max(raref_ub)),
	ylab="Num Unique Taxa Discovered",
	xlab="Num Samples Taken", xaxt="n", yaxt="n",
	main=OutputFileNameRoot,
	col=c(rep("black",NumSamples), rep("blue",BS_NumSamples-NumSamples))
);
abline(h=NumCategories);
text(1:BS_NumSamples, raref_med, labels=as.integer(raref_med), pos=4, cex=.7, col="grey25");

print(paste( ReadsPerSample, NumCategories, NumSamples, ci, sep="/"));
margin_text_string=
	sprintf("Cumulative Rarefaction: Reads/Sample = %i, Max Taxa = %i, Num Samples Provided = %i, %2.1f%% C.I.", 
		ReadsPerSample, num_nonzero_taxa, NumSamples, ci);
mtext(margin_text_string, side=3, line=0);
axis(side=1, at=1:BS_NumSamples);
axis(side=2, at=seq(0,max(raref_ub), 25));

###############################################################################

medians_out_fname=paste(OutputFileNameRoot, ".cumulative_rarefaction.medians", sep="");
medians_out_fh=file(medians_out_fname, "wt");
cat(OutputFileNameRoot, raref_med, file=medians_out_fh, sep="\t");
cat("\n", file=medians_out_fh, sep="");
close(medians_out_fh);

###############################################################################

cat("Done.\n");
q(status=0);
