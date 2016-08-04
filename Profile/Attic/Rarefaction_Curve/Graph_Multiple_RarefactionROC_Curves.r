#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n   ", script_name, "\n",
	"	-i <input rarefaction median file>\n",
	"\n",
	"This script will read in a file of rarefaction rat of change (roc) values and plot them on\n",
	"a single plot.  Colors will be assigned automatically.\n",
	"\n",
	"To generate the input for this script, you can use the UNIX cat command\n",
	"to concatenate all the .medians files together.  You can use the\n",
	"Rarefaction_Curve.r script to generate the roc. If you don't like\n",
	"the sample names, make sure you modify them before running this script.\n",
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
mat<-as.matrix(read.table(InputFileName, sep="\t", check.names=FALSE, row.names=1, fill=TRUE))

numSamples=nrow(mat);
numDonors=ncol(mat);

cat("Num Samples: ", numSamples, "\n", sep="");
cat("Num Donors: ", numDonors, "\n", sep="");

taxa_counts=as.vector(mat);
max=max(taxa_counts[is.finite(taxa_counts)]);
ylimit=c(0, max);
xlimit=c(0, numDonors);

colors=sample(rainbow(numSamples),numSamples);

pdf(paste(OutputFileNameRoot, ".combined_roc.pdf", sep=""), height=8.5, width=11);

reorder=order(mat[,numDonors],decreasing=TRUE);
mat=mat[reorder,];

plot(0,0,type="n", ylim=ylimit, xlim=xlimit,
	xlab="Num Donors",
	ylab="Num New Taxa Discovered",
	main="Rarefaction Rate of Change Curves"
	);
for(sample_idx in 1:numSamples){
	#print(mat[sample_idx,]);
	lines(1:numDonors, mat[sample_idx,],col=colors[sample_idx]);
}

legend(0, ylimit[2], rownames(mat), fill=colors, col=colors, bty="n"); 

cat("Done.\n");
q(status=0);
