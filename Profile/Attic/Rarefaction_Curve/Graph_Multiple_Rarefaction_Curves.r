#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"truncate_paths", "t", 0, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n   ", script_name, "\n",
	"	-i <input rarefaction median file>\n",
	"	[-t (truncate path for sample names)]",
	"\n",
	"This script will read in a file of rarefaction values and plot them on\n",
	"a single plot.  Colors will be assigned automatically.\n",
	"\n",
	"To generate the input for this script, you can use the UNIX cat command\n",
	"to concatenate all the .medians files together.  You can use the\n",
	"Rarefaction_Curve.r script to generate the medians.  If you don't like\n",
	"the sample names, make sure you modify them before running this script.\n",
	"\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileNameRoot=opt$input_file;
TruncatePaths=length(opt$truncate_paths)>0;

###############################################################################
# Load counts from file

cat("Working on ", InputFileName, "\n", sep="");

# Load summary_table.xls
mat<-as.matrix(read.table(InputFileName, sep="\t", check.names=FALSE, row.names=1, fill=TRUE))
oriMat=read.table(InputFileName, sep="\t",fill=TRUE) 
numSamples=nrow(mat);
numDonors=ncol(mat);
mat=cbind(0, mat); #start with zero

cat("Num Samples: ", numSamples, "\n", sep="");
cat("Num Donors: ", numDonors, "\n", sep="");

if(TruncatePaths){
	cat("Truncating paths from sample names.\n");
	sample_names=rownames(mat);
	for(i in 1:numSamples){
		split=strsplit(sample_names[i], "/")[[1]];
		last_token=length(split);
		cat("  ", sample_names[i], " -> ", split[last_token], "\n", sep="");
		sample_names[i]=split[last_token];
		
	}
	rownames(mat)=sample_names;
}


taxa_counts=as.vector(mat);
max=max(taxa_counts[is.finite(taxa_counts)]);
ylimit=c(0, max);
xlimit=c(0, numDonors);
names=c(rownames(mat));
colors=sample(rainbow(numSamples,gamma=.8),numSamples);

pdf(paste(OutputFileNameRoot, ".combined.pdf", sep=""), height=8.5, width=11);

arr=NULL;
for(i in 1:numSamples){
arr[i]=max(mat[i,],na.rm=TRUE)
}
#row.max=as.matrix(cbind(mat[,1],arr[1:length(arr)]));
#row.max<-row.max[rev(order(as.numeric(row.max[,2]))),];
#print(row.max);

mat=cbind(mat,arr[1:length(arr)]);
#print(mat);
reorder=rev(order(as.numeric(mat[,ncol(mat)],decreasing=TRUE,na.last=FALSE)));
mat=mat[reorder,];
mat<-mat[,1:numDonors,drop=FALSE]
#print (rownames(mat));
plot(0,0,type="o", ylim=ylimit, xlim=xlimit,
	xlab="Num Donors",
	ylab="Num Taxa",
	main="Rarefaction Curves"
	);
num=2;
#print (row.max[,1]);
for(sample_idx in 1:numSamples){
	lines(0:(length(mat[sample_idx,])-1), mat[sample_idx,],type="b",col=colors[sample_idx],pch=sample_idx,lty=sample_idx);
}
legend(0, ylimit[2], rownames(mat) , col=colors[1:numSamples],bty="n",pch=c(1:numSamples),lty=c(1:numSamples)); 

cat("Done.\n");
q(status=0);
