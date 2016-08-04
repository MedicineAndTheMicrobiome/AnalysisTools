#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"ubiquity", "u", 1, "numeric",
	"noise", "n", 1, "numeric",
	"output_file", "o", 2, "character",
	"generate_plot", "p", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	-u <ubiquity cutoff, if the mean is -u standard deviations greater than 0, remove it, 2 or 3 is good>\n",
	"	-n <noise cutoff, if the mean is is less than -n standard deviations from 0, remove it, 1 is good>\n",
	"	[-o <output summary table file name>\n",
	"	[-p <generate plot>]\n",
	"\n",	
	"This script will read in a summary table and then remove categories that are ubiquitous, ie 'definitely' in \n",
	"every sample and remove categories that are 'close' empty, ie. noise.\n",
	"\n",
	"To turn off the ubiquity filter, use a very large value, ie 100.\n",
	"To turn off the noise filter, set it to 0\n",
	"\n",
	"\n");

if(!length(opt$input_file) || !length(opt$ubiquity) || !length(opt$noise)){
	cat(usage);
	q(status=-1);
}

if(opt$ubiquity<0 || opt$noise<0){
	cat("Nonsensical negative cutoff.\n");
	q(status=-1);
}

InputFileName=opt$input_file;
ubiquity_cutoff=opt$ubiquity;
noise_cutoff=opt$noise;

if(!length(opt$output_file)){
	extension=sprintf("filtered.u%g_n%g", ubiquity_cutoff, noise_cutoff);
	output_root=gsub("\\.summary_table\\.xls", "", opt$input_file);
	output_root=paste(output_root, ".", extension, sep="");

	OutputFileName = paste(output_root, ".summary_table.xls", sep="");
	NoiseFileName = paste(output_root, ".noise.list", sep="");
	UbiquityFileName = paste(output_root, ".ubiquity.list", sep="");
	OutputPDFFileName = paste(output_root, ".pdf", sep="");
}else{
	OutputFileName=opt$output_file;
}

generate_plot=FALSE;
if(length(opt$generate_plot)){
	generate_plot=TRUE;	
}

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("Ubiquity Cutoff: ", ubiquity_cutoff, "\n");
cat("Noise Cutoff: ", noise_cutoff, "\n");
cat("\n");
cat("Generating Plot: ", generate_plot, "\n");

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

column_names=colnames(counts_mat);
noise_category_remove=rep(TRUE,num_categories);
ubiquity_category_remove=rep(TRUE,num_categories);
pdf(OutputPDFFileName, width=11, height=8.5);

max=max(normalized);

par(mfrow=c(2,1));

for(i in 1:num_categories){

	med=median(normalized[,i]);
	mean=mean(normalized[,i]);
	std=sd(normalized[,i]);
	
	lb=mean-std;
	ub=mean+std;

	noise_lb   =mean-(noise_cutoff*std);
	ubiquity_lb=mean-(ubiquity_cutoff*std);

	noise_category_remove[i]=(noise_lb<0);
	ubiquity_category_remove[i]=(ubiquity_lb>0);

	if(generate_plot==TRUE){

		# Plot with respect to max proportion
		hist(normalized[,i], breaks=50, freq=FALSE, xlim=c(0,max*1.1), main=column_names[i], col="grey75", border="grey50", xlab="Propotion");
		lines(density(normalized[,i]), col="blue", lwd=2);
		abline(v=mean, col="black", lwd=1.5);
		abline(v=noise_lb, col="green");
		abline(v=ubiquity_lb, col="green");
		abline(v=ub, col="blue");
		mtext("Global", line=1, cex=.8);
		mtext(sprintf("   Noise Remove: %i", noise_category_remove[i]), line=0, cex=.75);
		mtext(sprintf("Ubiquity Remove: %i", ubiquity_category_remove[i]), line=-.5, cex=.75);

		# Plot to maximize resolution of specific category
		localmax=max(normalized[,i]);
		hist(normalized[,i], breaks=50, freq=FALSE, xlim=c(0,localmax*1.1), main=column_names[i], col="grey75", border="grey50", xlab="Proportion");
		lines(density(normalized[,i]), col="blue", lwd=2);
		abline(v=mean, col="black", lwd=1.5);
		abline(v=noise_lb, col="green");
		abline(v=ubiquity_lb, col="green");
		abline(v=ub, col="blue");
		mtext("Local", line=1, cex=.8);
	}
}

###############################################################################

noise_removed=sum(noise_category_remove);
ubiquity_removed=sum(ubiquity_category_remove);
category_keep=!(noise_category_remove | ubiquity_category_remove);
num_keep=sum(category_keep);

cat("Categories removed:\n");
cat("      Noise: ", noise_removed, "\n", sep="");
cat("   Ubiquity: ", ubiquity_removed, "\n", sep="");
cat("\n");
cat("Num Categories kept: ", num_keep, "\n", sep="");

outmat=counts_mat[,category_keep];

###############################################################################

category_names=colnames(counts_mat);
#print(category_names[noise_category_remove]);
#cat("\n");
#print(category_names[ubiquity_category_remove]);

cat("Writing Removed Noise...\n");
fc=file(NoiseFileName, "w");
write(file=fc, category_names[noise_category_remove]);
close(fc)

cat("Writing Removed Ubiquity...\n");
fc=file(UbiquityFileName, "w");
write(file=fc, category_names[ubiquity_category_remove]);
close(fc);

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
print(warnings());

q(status=0)
