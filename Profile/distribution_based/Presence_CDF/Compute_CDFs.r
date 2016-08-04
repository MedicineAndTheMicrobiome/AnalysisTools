#!/usr/bin/env Rscript

###############################################################################
#                                                                             #
#       Copyright (c) 2013 J. Craig Venter Institute.                         #
#       All rights reserved.                                                  #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                             #
###############################################################################
###############################################################################

library('getopt');

MinExpectedAbundance=1e-4;
LogSteps=1000;
LinearSteps=1000;

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 2, "character",
	"min_exp_abund", "a", 2, "numeric",
	"log_steps", "g", 2, "numeric",
	"linear_steps", "n", 2, "numeric",
	"bs_smooth", "s", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"	-i <input summary_table.xls file>\n",
	"	[-o <output cdf name root>]\n",
	"	[-a <minimum expected abundance, default=", MinExpectedAbundance, "]\n",
	"	[-g <log steps to take, default=", LogSteps, "]\n",
	"	[-n <linear steps to take, default=", LinearSteps, "]\n",
	"	[-s (Flag to enable smoothing with bootstrapping, default=no smoothing)]\n",
	"\n",
	"Computes the CDF for each organism in the specified summary table.\n",
	"\n",
	"This script transforms your summary table into a set of points that underly the CDF\n",
	"curves describing each taxon across the cohort in your summary table.\n",
	"\n",
	"Note: If you are planning to graph the CDF, then leave -g and -n to be at least 1000.\n",
	"If you are going to compare the CDF lines, you should make sure that you take the same\n",
	"log or linear steps, or else you will have mismatching abundances.  The steps\n",
	"refer to the the increments between abundances that will be taken to compute the ubiquities.\n",
	"Because the taxa abundances are very low for some taxa, i.e. the rank abundance curve\n",
	"has a long tail, the log steps will reduce the steps taken (exponentially) between abundances\n",
	"to increase the precision for both graphing and comparing taxa at low abundances.\n",
	"\n",
	"The minimum expected abundances (-a) seems to about right for 454 sequencing depths.\n",
	"This may need to be reduced for higher depth sequencing technologies.\n",
	"\n", sep="");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputRoot=gsub("\\.summary_table\\.xls$", "", InputFileName);
OutputRoot=gsub("\\.summary_table\\.tsv$", "", InputFileName);

if(length(opt$output_file)){
	OutputRoot=opt$output_file;
}

if(length(opt$min_exp_abund)){
	MinExpectedAbundance=opt$min_exp_abund;
}

if(length(opt$log_steps)){
	LogSteps=opt$log_steps;
}

if(length(opt$linear_steps)){
	LinearSteps=opt$linear_steps
}

if(length(opt$bs_smooth)){
	Smooth=T;
}else{
	Smooth=F;
}


cat("Input File: ", InputFileName, "\n");
cat("Output File Root: ", OutputRoot, "\n");
cat("Smoothing?: ", Smooth, "\n");

################################################################################

load_summary_table=function(filename){
	st=as.matrix(read.table(filename, header=TRUE, sep="\t", row.names=1, check.names=F));
	return(st[,2:ncol(st)]);
}

#-------------------------------------------------------------------------------

normalize=function(st){
	sums=apply(st, 1, sum);
	n=matrix(0,nrow=nrow(st), ncol=ncol(st));
	rownames(n)=rownames(st);
	colnames(n)=colnames(st);
	for(i in 1:nrow(st)){
		n[i,]=st[i,]/sums[i];
	}
	return(n);
}

#-------------------------------------------------------------------------------

compute_cdf_on_vector=function(vector, abundances){
	num_samples=length(vector);
	num_abundances=length(abundances);
	ubiquities=rep(0,num_abundances);
	
	for(i in 1:num_abundances){
		ubiquities[i]=sum(abundances[i]<=vector);	
	}

	ubiquities=ubiquities/num_samples;
	return(ubiquities);
}

#-------------------------------------------------------------------------------

write_cdf_to_file=function(matrix, fh){
	ncol=ncol(matrix);
	nrow=nrow(matrix);
	increments=colnames(matrix);
	taxa=rownames(matrix);
	cat(file=fh, ",");
	cat(file=fh, paste(increments, collapse=","), "\n");
	for(i in 1:nrow){
		cat(file=fh, taxa[i], paste(matrix[i,], collapse=","), sep=",");
		cat(file=fh, "\n");
	}	
}

#-------------------------------------------------------------------------------

generate_bs_samples=function(st){

	# Compute probabilities for each sample
	nst=normalize(st);
	
	# Get summary table stats
	num_samples=nrow(nst);	
	num_reads=apply(st, 1, sum);
	num_categories=ncol(nst);
	cat_names=colnames(st);
	samp_names=rownames(st);
	
	cat("Num Samples: ", num_samples, "\n");
	cat("Num Reads: ", num_reads, "\n");

	# Determine how many bootstrap samples to generate
	num_samples_to_generate=max(0, 1000-num_samples);
	cat("Num Samples to BS Generate: ", num_samples_to_generate, "\n");

	# Allocate matrix to store bs samples
	bs_sample=matrix(NA, nrow=num_samples_to_generate, ncol=num_categories);

	# Generate sequence of original samples to bs from
	sample_ix=sample(num_samples, num_samples_to_generate, replace=T);

	# Generate bs samples
	for(bix in 1:num_samples_to_generate){

		# Smoothing across ubiquity
		bs_sample[bix,]=t(rmultinom(1, num_reads[sample_ix[bix]], nst[sample_ix[bix],]));

		# Smoothing across abundance
		gtzero=nst[sample_ix[bix],]>0;
		bs_sample[bix,]=bs_sample[bix,]+runif(num_categories, -.5, .5)*gtzero;

		# Don't let zero counts go to zero
		bs_sample[bix,]=sapply(bs_sample[bix,], function(x){max(0,x);});
	}

	# Generate names for debugging
	bs_names=paste("Bootstrap_", (1:num_samples_to_generate), "_", samp_names[sample_ix], sep="");
	rownames(bs_sample)=bs_names;
	colnames(bs_sample)=cat_names;
	
	# Attach original samples to list
	bs_sample=rbind(bs_sample, st);

	return(bs_sample);
	

}

################################################################################

cat("Loading summary table...\n");
st=load_summary_table(InputFileName);
cat("ok.\n");

if(Smooth){
	st=generate_bs_samples(st);
}

num_taxa=ncol(st);
num_samples=nrow(st);
#print(st);

# Normalize counts
nst=normalize(st);

# Compute min non zero abundance
nonzero_abundance=nst[nst>0];
min_nonzero_abundance=min(nonzero_abundance);
cat("Min nonzero abundance detected: ", min_nonzero_abundance, "\n");

# Compute the abundances for which we will compute the ubiquity for.
min_log_exp_abundance=log(MinExpectedAbundance,10);
log_abundance_steps=10^(seq(min_log_exp_abundance, 0, length.out=LogSteps));
linear_abundance_steps=seq(0,1, length.out=LinearSteps);
abundance_steps=sort(unique(c(log_abundance_steps, linear_abundance_steps)));
num_abundances=length(abundance_steps);
cat("Num Abundances to Measure: ", num_abundances, "\n", sep="");

# Initialize output cdf matrix
output_matrix=matrix(0,nrow=num_taxa, ncol=num_abundances);
colnames(output_matrix)=sprintf("%1.16f",abundance_steps);
rownames(output_matrix)=colnames(nst);

# For each taxa, compute the CDF
cat("Num of Taxa: ", num_taxa, "\n");

progress=as.integer(seq(1,num_taxa, length.out=20));
for(i in 1:num_taxa){
	output_matrix[i,]=compute_cdf_on_vector(nst[,i], abundance_steps)
	if(any(i==progress)){
		cat("\n", sprintf("%3.0f%%", 100*i/num_taxa), " completed ");
	}else{
		if(!i%%10){
			cat(".");
		}
	}
}
cat("\n");
#print(output_matrix);

################################################################################

ext="";
if(LinearSteps>0){
	ext=paste(ext, ".lin", sep="");	
}
if(LogSteps>0){
	ext=paste(ext, ".log", sep="");
}

if(Smooth){
	ext=paste(ext, ".smoothed", sep="");
}

fh=file(paste(OutputRoot, ext, ".cdf", sep=""), "w");
write_cdf_to_file(output_matrix, fh);
close(fh);

################################################################################

writeLines("Done.\n")

q(status=0)
