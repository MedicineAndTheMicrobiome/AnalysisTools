#!/usr/bin/env Rscript

###############################################################################

progname <- commandArgs(FALSE)[4]
args <- commandArgs(TRUE)

#CUTOFF=.05;
CUTOFF=.001;
arg_count=1

if(is.na(args[arg_count])){

	script_name <- unlist(strsplit(progname,"="))[2]
	usage <- paste (
		"\nUsage:\n\t", script_name, "\n\t\t<summary_table.xls>\n\n",
		"This script will read in a summary_table and filter out all the counts\n",
		"that are not statistically significantly non-zero.  The code uses\n",
		"the binomial distribution to determine the probability a category\n",
		"will be 0 upon resampling.  If the probability is greater than a cutoff\n",
		"eg. .1%, then the count is set to zero.  In other words, if the\n",
		"count will be nonzero 99.9% of the time, then count will not be eliminated.\n",
		"\n")

	writeLines(usage)
	writeLines("Input FileName not defined.\n")
	quit(status=0)
}

###############################################################################
# Main program loop

while(!(is.na(args[arg_count]))){
	InputFileName=args[arg_count]
	OutputFileName=sub(".summary_table.xls", ".filt.summary_table.xls", InputFileName);

	if(InputFileName == OutputFileName){
		OutputFileName=paste(OutputFileName, ".filt.summary_table.xls");
	}

	cat("\n")
	cat("Input File Name: ", InputFileName, "\n");
	cat("Output File Name: ", OutputFileName, "\n");       

	###############################################################################
	###############################################################################

	# Load data
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

	###############################################################################
	# Use binomial distribution to determine probability count will not be zero when resampled
	prob_cutoff_mat=matrix(logical(), nrow=num_samples, ncol=num_categories);
	for(samp_idx in 1:num_samples){
		for(cat_idx in 1:num_categories){
			prob_cutoff=pbinom(0, sample_totals[samp_idx], normalized[samp_idx, cat_idx]);
			prob_cutoff_mat[samp_idx, cat_idx]=prob_cutoff<CUTOFF;
		}
	}
	# Probability of being 0 more than CUTOFF percent of the time
	#print(prob_cutoff_mat);
	
	# Filtered results
	counts_mat[!prob_cutoff_mat]=0;
	#print(counts_mat);

	###############################################################################
	# Output
	fc=file(OutputFileName, "w");
	write(paste("sample_id", paste(colnames(inmat), collapse="\t")), file=fc);
	sample_names=rownames(inmat);
	for(samp_idx in 1:num_samples){
		total=sum(counts_mat[samp_idx,]);
		outline=paste(sample_names[samp_idx],total,paste(counts_mat[samp_idx,], collapse="\t"), sep="\t");
		write(outline, file=fc);
	}
	close(fc);	

	###############################################################################

	arg_count=arg_count+1;
}

writeLines("Done.\n")
print(warnings());

q(status=0)
