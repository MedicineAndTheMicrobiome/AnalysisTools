#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"reads_per_sample", "r", 2, "numeric",
	"num_samples", "d", 2, "numeric",
	"num_bootstraps", "b", 2, "numeric",
	"low_abund_threshold", "t", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n\t\t<Input summary_table.xls FileName>\n\n",
	"	-i <input summary_table.xls file>\n",
	"	[-r <Reads per sample>]\n",
	"	[-d <Num samples to graph>]\n",
	"	[-b <Num bootstrap, default=40>]\n",
	"	[-t <Low abundancy threshold, default=0, but .01 would be 1%>]\n",
	"\n",
	"This script will generate a rarefaction curve based on the input summary table.\n",
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

# Low abundance threshold
if(length(opt$low_abund_threshold)){
	cat("Low abundance threshold specified.\n");
	LowAbundThresh=opt$low_abund_threshold;
}else{
	LowAbundThresh=0;
}
cat("Low Abundance Threshold: ", LowAbundThresh, "\n");
cat("\n");

###############################################################################

compute_num_organisms=function(num_individuals, ReadsPerSample, prob_mat, cutoff){
	
	num_donors=nrow(prob_mat);
	num_categories=ncol(prob_mat);

	# Select donors
	donor_samples=sample(1:num_donors, num_individuals, replace=TRUE);
	#cat("Selected donors:\n");
	#print(donor_samples);

	# Select reads from each donor
	unique_categories=numeric(0);
	for(donor in donor_samples){
		samples=sample(1:num_categories, ReadsPerSample, replace=TRUE, prob=prob_mat[donor,]);
		t=table(samples);

		categories=as.integer(attr(t,"dimnames")[[1]]);
		probs=as.vector(t)/ReadsPerSample;
		filt_categories=categories[probs>cutoff];

		#print(filt_categories);
		#cat("---\n");
		#print(samples);

		#unique_categories=unique(c(unique_categories, samples));
		unique_categories=unique(c(unique_categories, filt_categories));
	}

	return(length(unique_categories));
} 

###############################################################################

ci=100*(NumBootstraps-2)/NumBootstraps;

raref_min=rep(0, NumSamples);
raref_max=rep(0, NumSamples);
raref_med=rep(0, NumSamples);
raref_mean=rep(0, NumSamples);

for(num_individuals in 1:BS_NumSamples){

	num_org = rep(0, NumBootstraps);
	for(bs in 1:NumBootstraps){
		#cat("Num Individuals: ", num_individuals, " Reads/Sample:", ReadsPerSample, "\n");
		num_org[bs]=compute_num_organisms(num_individuals, ReadsPerSample, prob_mat, LowAbundThresh);
	}

	sorted_num_org=sort(num_org);
	raref_min[num_individuals]=sorted_num_org[2];
	raref_med[num_individuals]=median(sorted_num_org);
	#raref_med[num_individuals]=mean(sorted_num_org);
	raref_max[num_individuals]=sorted_num_org[NumBootstraps-1];
	raref_mean[num_individuals]=mean(sorted_num_org);

	cat(
		"Num Ind = ", num_individuals, " mean = ", raref_mean[num_individuals],
		" med = ", raref_med[num_individuals],
		" (", raref_min[num_individuals], ", " , raref_max[num_individuals],")",
		"\n", sep=""
	);
	
}


###############################################################################

count_mat=apply(count_mat, 2, sum);
num_nonzero_taxa=sum(count_mat>0);
cat("Num Nonzero Taxa:", num_nonzero_taxa, "\n");

###############################################################################

pdf(paste(OutputFileNameRoot, ".rarefaction.pdf", sep=""), height=8.5, width=11);

library(gplots);

print(paste( ReadsPerSample, NumCategories, NumSamples, ci, sep="/"));

###############################################################################
# Plot Medians
plotCI(x=1:BS_NumSamples, y=raref_med, ui=raref_max, li=raref_min, ylim=c(0,max(raref_max)),
	ylab="Num Unique Taxa Discovered",
	xlab="Num Samples Taken", xaxt="n", yaxt="n",
	main=OutputFileNameRoot,
	col=c(rep("black",NumSamples), rep("blue",BS_NumSamples-NumSamples))
);
text(1:BS_NumSamples, raref_med, labels=raref_med, pos=4, cex=.7, col="grey25");
abline(h=NumCategories);

margin_text_string=sprintf("Reads/Sample = %i, Max Taxa = %i, Num Samples Provided = %i, %2.1f%% C.I., Low Abundance Thresh = %4.4f",
			 ReadsPerSample, num_nonzero_taxa, NumSamples, ci, LowAbundThresh);
mtext(margin_text_string, side=3, line=0);
mtext("MEDIAN", side=3, line=-1);
axis(side=1, at=1:BS_NumSamples);
axis(side=2, at=floor(seq(0,max(raref_max), length.out=10)));

###############################################################################
#Plot Means;

plotCI(x=1:BS_NumSamples, y=raref_mean, ui=raref_max, li=raref_min, ylim=c(0,max(raref_max)),
	ylab="Num Unique Taxa Discovered",
	xlab="Num Samples Taken", xaxt="n", yaxt="n",
	main=OutputFileNameRoot,
	col=c(rep("black",NumSamples), rep("blue",BS_NumSamples-NumSamples))
);
text(1:BS_NumSamples, raref_mean, labels=raref_mean, pos=4, cex=.7, col="grey25");
abline(h=NumCategories);

margin_text_string=sprintf("Reads/Sample = %i, Max Taxa = %i, Num Samples Provided = %i, %2.1f%% C.I., Low Abundance Thresh = %4.4f",
			 ReadsPerSample, num_nonzero_taxa, NumSamples, ci, LowAbundThresh);
mtext(margin_text_string, side=3, line=0);
mtext("MEAN", side=3, line=-1);
axis(side=1, at=1:BS_NumSamples);
axis(side=2, at=floor(seq(0,max(raref_max),length.out=10)));

###############################################################################
# Plot both superimposed

plot(x=1:BS_NumSamples, type="n", ylim=c(0,max(raref_max)),
	ylab="Num Unique Taxa Discovered",
	xlab="Num Samples Taken", xaxt="n", yaxt="n",
	main=OutputFileNameRoot
);

lines(1:BS_NumSamples, raref_med, pch=1, cex=1.2, col="blue", type="b");
text(1:BS_NumSamples, raref_med, labels=raref_med, adj=c(2,-2), pos=4, cex=.6, col="blue");

lines(1:BS_NumSamples, raref_mean, pch=1, cex=1.2, col="green", type="b");
text(1:BS_NumSamples, raref_mean, labels=sprintf("%1.1f", raref_mean), adj=c(-2,2), pos=2, cex=.6, col="green");

margin_text_string=sprintf("Reads/Sample = %i, Max Taxa = %i, Num Samples Provided = %i, %2.1f%% C.I., Low Abundance Thresh = %4.4f",
			 ReadsPerSample, num_nonzero_taxa, NumSamples, ci, LowAbundThresh);
mtext(margin_text_string, side=3, line=0);

mtext("Median", side=3, line=-1, col="blue");
mtext("Mean", side=3, line=-2, col="green");

axis(side=1, at=1:BS_NumSamples);
axis(side=2, at=floor(seq(0,max(raref_max),length.out=10)));

###############################################################################

# Compute rate of change of discovering new taxa

rate_of_change=function(rarefaction){
	num_samples_taken=length(rarefaction);
	roc=rep(0, num_samples_taken);
	zero_starting=c(0, rarefaction);
	for(i in 1:num_samples_taken){
		roc[i]=zero_starting[i+1]-zero_starting[i];
	}
	return(roc);
}

roc_med=rate_of_change(raref_med);
roc_mean=rate_of_change(raref_mean);

par(mfrow=c(2,1));

# Plot Median ROC
plot(1:BS_NumSamples, roc_med, type="b", xlab="Num Samples Collected", ylab="New Taxa Discovered", main=OutputFileNameRoot, 
	xaxt="n", yaxt="n");
mtext("Median", side=3, line=0);
axis(side=1, at=1:BS_NumSamples);
axis(side=2, at=floor(seq(0,max(roc_med),length.out=10)));

# Plot Mean ROC
plot(1:BS_NumSamples, roc_mean, type="b", xlab="Num Samples Collected", ylab="New Taxa Discovered", main=OutputFileNameRoot,
	xaxt="n", yaxt="n");
mtext("Mean", side=3, line=0);
axis(side=1, at=1:BS_NumSamples);
axis(side=2, at=floor(seq(0,max(roc_mean),length.out=10)));

###############################################################################
###############################################################################

values_out_fname=paste(OutputFileNameRoot, ".rarefaction.values", sep="");
medians_roc_out_fname=paste(OutputFileNameRoot, ".rarefaction_roc.medians", sep="");
means_roc_out_fname=paste(OutputFileNameRoot, ".rarefaction_roc.means", sep="");

values_out_fh=file(values_out_fname, "wt");
medians_roc_out_fh=file(medians_roc_out_fname, "wt");
means_roc_out_fh=file(means_roc_out_fname, "wt");

cat(paste(OutputFileNameRoot, ".median", sep=""), raref_med, file=values_out_fh, sep="\t");
cat("\n", file=values_out_fh, sep="");
cat(paste(OutputFileNameRoot, ".mean", sep=""), raref_mean, file=values_out_fh, sep="\t");
cat("\n", file=values_out_fh, sep="");
cat(paste(OutputFileNameRoot, ".lo_bound", sep=""), raref_min, file=values_out_fh, sep="\t");
cat("\n", file=values_out_fh, sep="");
cat(paste(OutputFileNameRoot, ".up_bound", sep=""), raref_max, file=values_out_fh, sep="\t");
cat("\n", file=values_out_fh, sep="");

cat(OutputFileNameRoot, roc_med, file=medians_roc_out_fh, sep="\t");
cat("\n", file=medians_roc_out_fh, sep="");

cat(OutputFileNameRoot, roc_mean, file=means_roc_out_fh, sep="\t");
cat("\n", file=means_roc_out_fh, sep="");

close(values_out_fh);
close(medians_roc_out_fh);
close(means_roc_out_fh);

###############################################################################

cat("Done.\n");
q(status=0);
