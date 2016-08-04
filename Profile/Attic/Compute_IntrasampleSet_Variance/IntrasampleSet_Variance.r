#!/usr/local/bin/Rscript

###############################################################################

progname <- commandArgs(FALSE)[4]
args <- commandArgs(TRUE)

arg_count=1

if(is.na(args[arg_count])){

	script_name <- unlist(strsplit(progname,"="))[2]
	usage <- paste (
		"\nUsage:\n\t", script_name, "\n\t\t<Input summary_table.xls FileName>\n\n",
		"Computes variance of samples with in set by averaging the euclidean distance between samples.\n",
		"\n")

	writeLines(usage)
	writeLines("Input FileName not defined.\n")
	quit(status=0)
}

###############################################################################
# Main program loop

while(!(is.na(args[arg_count]))){
	InputFileName=args[arg_count]
	OutputFileName = paste(InputFileName, ".variance", sep="")

	cat("\n")
	cat("Input File Name: ", InputFileName, "\n")

	###############################################################################
	###############################################################################

	# Load data
	mat<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, row.names=1))

	count_mat=mat[,2:ncol(mat)]

	#print(count_mat);

	num_samples=(nrow(count_mat));
	num_cat=(ncol(count_mat));

	cat("   Num samples:", num_samples, "\n", sep="");
	cat("Num categories:", num_cat, "\n", sep="");

	normalized=matrix(rep(0,num_samples*num_cat), ncol=num_cat, nrow=num_samples);
	keep=logical(0);
	for(i in 1:num_samples){
		sample_totals=sum(count_mat[i,]);
		if(sample_totals>0){
			keep[i]=TRUE;
		}else{
			keep[i]=FALSE;
		}
		normalized[i,]=count_mat[i,]/sample_totals;
	}

	normalized=normalized[keep,];
	cat("Samples kept: ", sum(keep), "\n");

	#print(normalized);

	dist=dist(normalized);

	#print(dist);

	sum_of_dist=mean(dist);
	
	fc=file(OutputFileName, "w");
	outline=paste(InputFileName, sum_of_dist, sep="\t");
	write(outline, fc);

	close(fc);
	
	###############################################################################

	arg_count=arg_count+1;
}

writeLines("Done.\n")

q(status=0)
