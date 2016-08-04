#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"num_samples", "n", 1, "numeric",
	"depth", "d", 1, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"	-i <input frequency file>\n",
	"	-n <number of samples>\n",
	"	-d <depth per sample>\n",
	"\n",
	"This script will read in a frequencies file (or a count file)\n",
	"and then generate a rarefaction curve based  on that distribution\n",
	"assuming the number of samples specified and depth.\n",
	"\n", sep="");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
NumSamples=opt$num_samples;
SampleDepth=opt$depth;

cat("Input file: ", InputFileName, "\n");
cat("Num Samples: ", NumSamples, "\n");
cat("Sample Depth: ", SampleDepth, "\n");
cat("\n");

###############################################################################

ci = function(x, alpha){
	n=length(x);
	med=median(x);
	ba=1-(n-2)/n;
	#cat("Best alpha = ", ba, "\n");
	if(ba <= (alpha+.0000000000000001)){
		sorted=sort(x);
		lb=sorted[floor(n*(alpha/2))+1];
		ub=sorted[ceiling(n*(1-(alpha/2)))];
		return(c(med,lb,ub));
	}else{
		return(c(med,NA,NA))
	}
}

###############################################################################
# Load counts from file

cat("Loading counts or frequencies", InputFileName, "\n", sep="");
frequencies=scan(InputFileName);
num_items=length(frequencies);

cat("Num items to sample from: ", num_items, "\n");

total=sum(frequencies);
cat("Sum of frequencies: ", total, "\n");

norm_freq=frequencies/total;

NUM_BOOTSTRAPS=40; # For 95% confidence intervals

# Compute counts
item_ids=1:num_items;
counts=matrix(0, nrow=NumSamples, ncol=NUM_BOOTSTRAPS);
for(bs in 1:NUM_BOOTSTRAPS){
	cat("Working on bootstrap iteration: ", bs, "\n");
	cumulative_unique=c();
	for(n in 1:NumSamples){
		cat("\tn=", n, "\n", sep="");
		samples=sample(item_ids, size=SampleDepth, replace=TRUE, prob=norm_freq);
		cumulative_unique=unique(c(cumulative_unique, samples));
		counts[n, bs]=length(cumulative_unique);
	}
}

# Compute confidence intervals
conf_info=matrix(0, nrow=NumSamples, ncol=3);
for(n in 1:NumSamples){
	conf_info[n,]=ci(counts[n,], .05);
}
	
# Output confidence information
fh=file(paste(InputFileName, ".rarefaction", sep=""), "w");
for(n in 1:NumSamples){
	cat(file=fh, n, ",", paste(conf_info[n,], collapse=","), "\n", sep="");
}
close(fh);

###############################################################################

cat("Done.\n");
q(status=0);
