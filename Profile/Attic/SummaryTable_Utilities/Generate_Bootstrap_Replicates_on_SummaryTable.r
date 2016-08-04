#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"bootstraps_per_sample", "b", 1, "numeric",
	"reads_per_sample", "r", 1, "numeric",
	"unique", "u", 2, "logical",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	-b <num bootstraps per sample>\n", 
	"	[-r <reads per sample]\n",
	"	[-u (flag if you want an index assigned to replicate to make it unique)]\n",
	"	[-o <output summary table file name>\n",
	"\n",	
	"This script will read in a summary table, and generate <num bootstraps per sample> bootstrap replicates for each sample.\n",
	"If -u is not specified, multiple bootstraps will have the same name.\n",
	"If -r <reads per sample> is not specified, then the number of reads for each bootstrap will be the same as the sample count.\n",
	"\n",
	"\n");

if(!length(opt$input_file) || !length(opt$bootstraps_per_sample)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$output_file)){
	OutputFileName = paste(gsub("\\.summary_table\\.xls$", "", opt$input_file), "", sep="");
}else{
	OutputFileName=opt$output_file;
}

###############################################################################

InputFileName=opt$input_file;
NumBootstraps=opt$bootstraps_per_sample;

cat("\n");
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("\n");

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

sample_names=rownames(counts_mat);

###############################################################################
# Normalize

# Sum sample totals
sample_totals=apply(counts_mat, 1, sum);
#print(sample_totals);

# normalize, to compute probabilities
normalized=matrix(0, nrow=num_samples, ncol=num_categories);
for(i in 1:num_samples){
	normalized[i,]=counts_mat[i,]/sample_totals[i];
}
#print(normalized);

###############################################################################
# Compute padding of numbers with left 0's so names line up
num_width=ceiling(log(NumBootstraps, base=10)) + 1;
unique=1:num_samples;

for(bs in 1:NumBootstraps){
	
	fmt_str=paste("%0", num_width, "i", sep="");
	outname=paste(OutputFileName, ".", sprintf(fmt_str, bs), ".summary_table.xls", sep="");
	fc=file(outname, "w");
	write(paste("sample_id", paste(colnames(inmat), collapse="\t"), sep="\t"), file=fc);
	unique_counter=1;

	# Generate replicates
	samples=sample(1:num_samples, num_samples, replace=TRUE);
	for(i in samples){
	
		cat("Bootstrapping: ", sample_names[i], "\n");

		# Determining how many reads to use
		if(is.null(opt$reads_per_sample)){
			reads_per_sample=sample_totals[i];	
		}else{
			reads_per_sample=opt$reads_per_sample;
		}

		obs=sample(1:num_categories, replace=TRUE, reads_per_sample, prob=normalized[i,]);
		counts=as.vector(table(c(obs,1:num_categories))-1);

		if(!is.null(opt$unique)){
			name=paste(sample_names[i], ".", unique[unique_counter], sep="");
			unique_counter= unique_counter + 1;
		}else{
			name=sample_names[i];
		}

		outline=paste(name, reads_per_sample, paste(counts, collapse="\t"), sep="\t");
		write(outline, file=fc);
	}
	close(fc);	
}



###############################################################################

writeLines("Done.\n")
w=warnings();
if(!is.null(w)){
	print(w);
}

q(status=0)
