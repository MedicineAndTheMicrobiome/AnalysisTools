#!/usr/bin/env Rscript

###############################################################################

library(MASS)
library('getopt');
library('vegan');

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 2, "character",
	"dist_mat_type", "d", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\n\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.xls file>\n",
	"	[-o <output file root name, default is input file base name>]\n",
	"	[-d <euc/wrd/man/bray/horn/bin/gow/tyc, default is euc>]\n",
	"\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}
InputFileName=opt$input_file;

if(length(opt$output_file)>0){
	OutputFileRoot=opt$output_file;
}else{
	OutputFileRoot=InputFileName;
	OutputFileRoot=gsub(".summary_table.xls$", "", OutputFileRoot);
	OutputFileRoot=gsub(".summary_table.tsv$", "", OutputFileRoot);
	cat("No output file root specified.  Using input file name as root.\n");
}

NotDefaultDistanceMatrix=!is.null(opt$dist_mat_type);
MatrixType=opt$dist_mat_type;

if(NotDefaultDistanceMatrix){
	if(MatrixType == "wrd"){
		cat("Using Weighted Rank Difference as distance.\n");
		type="wrd";
	}else if(MatrixType == "man"){
		cat("Using Manhattan as distance.\n");
		type="man";
	}else if(MatrixType == "bray"){
		cat("Using Bray-Curtis as distance\n");
		type="bray";
	}else if(MatrixType == "horn"){
		cat("Using Horn-Morisita as distance\n");
		type="horn";
	}else if(MatrixType == "bin"){
		cat("Using Binomial index as distance\n");
		type="bin";
	}else if(MatrixType == "gow"){
		cat("Using Gower index as distance\n");
		type="gow";
	}else if(MatrixType == "euc"){
		cat("Using Euclidean as distance.\n");
		type="euc";
	}else if(MatrixType == "tyc"){
		cat("Using Yue & Clayton's measure of dissimilarity as a distance.\n");
		type="tyc";
	}

}else{
	cat("Using Euclidean as distance.\n");
	type="euc";
}

###############################################################################
# Load the WRD code

if(type == "wrd"){
	path_comp=strsplit(script_name, "/")[[1]];
	bin_comp_idx=length(path_comp);
	bin_path=paste(path_comp[-bin_comp_idx], collapse="/", sep="");
	if(nchar(bin_path)==0){
	        bin_path=".";
	}
	cat("Binary path: '", bin_path, "'\n", sep="");

	source(paste(bin_path, "WeightedRankDifference.r", sep="/")); #changed path. orig Cluster/WeightedRankDifference.r
}

# Function for ThetaYC
# See http://www.mothur.org/wiki/Thetayc for formula

tyc_fun=function(v1, v2){
	sum_intersect=sum(v1*v2);
	sum_sqrd_diff=sum((v1-v2)^2);
	denominator=sum_sqrd_diff + sum_intersect;
	tyc=1-(sum_intersect/denominator);
	return(tyc);
}

thetaYC=function(matrix){
	
	nsamples=nrow(matrix);
	ycdist=matrix(0, ncol=nsamples, nrow=nsamples);
	for(i in 1:nsamples){
		for(j in 1:nsamples){
			if(i<j){
				ycdist[i,j]=tyc_fun(matrix[i,], matrix[j,]);
			}else{
				ycdist[i,j]=ycdist[j,i];			
			}
		}
	}
	
	as.dist(return(ycdist));
}

###############################################################################

DistanceMatrixTXT = paste(OutputFileRoot, ".", type, ".distmat", sep="");

cat("Output Distance Matrix: ", DistanceMatrixTXT, "\n", sep="");

###############################################################################
# Load data
InMat=as.matrix(read.delim(InputFileName, sep="\t", header=TRUE, row.names=1, check.names=FALSE, comment.char="", quote=""))

#print(dim(InMat));

num_categories=ncol(InMat)-1;
num_samples=nrow(InMat);

cat("Num Categories: ", num_categories, "\n", sep="");
cat("Num Samples: ", num_samples, "\n", sep="");

countsMat=InMat[,2:(num_categories+1)];
#print(countsMat)

Categories=colnames(countsMat);
SampleNames=rownames(countsMat);

#------------------------------------------------------------------------------

# Sum up the number of members in each row
sample_counts=apply(countsMat, 1, sum);
cat("Sample Total:\n");
print(sample_counts);

# Removing zero count samples
nonzero_count_samples=(sample_counts!=0);
if(!all(nonzero_count_samples)){

	cat("Warning: Zero count samples found in summary table.\n");
	cat("Removing zero count samples:\n");
	print(SampleNames[!nonzero_count_samples]);
	cat("\n");

	countsMat=countsMat[nonzero_count_samples,, drop=F];
	sample_counts=apply(countsMat, 1, sum);
	SampleNames=rownames(countsMat);
	num_samples=nrow(countsMat);
	cat("New number of sampels: ", num_samples, "\n");
}

# Normalize counts
cat("Normalizing counts...\n");
normalized=matrix(0, nrow=num_samples, ncol=num_categories);
rownames(normalized)=rownames(countsMat);
colnames(normalized)=colnames(countsMat);
for(i in 1:num_samples){
	normalized[i,]=countsMat[i,]/sample_counts[i];
}	

# Compute distances
cat("Computing distances...\n");
if(type=="euc"){
	dist_mat=dist(normalized);
}else if (type=="wrd"){
	dist_mat=weight_rank_dist_opt(normalized, deg=4);
}else if (type=="man"){
	dist_mat=vegdist(normalized, method="manhattan");
}else if (type=="bray"){
	dist_mat=vegdist(normalized, method="bray");
}else if (type=="horn"){
	dist_mat=vegdist(normalized, method="horn");
}else if (type=="bin"){
	dist_mat=vegdist(normalized, method="bin");
}else if (type=="gow"){
	dist_mat=vegdist(normalized, method="gower");
}else if (type=="tyc"){
	dist_mat=thetaYC(normalized);
}
# Output distance matrix
asFull=as.matrix(dist_mat);
fh=file(DistanceMatrixTXT, "w");

for(i in 1:num_samples){
#	cat(file=fh," ", SampleNames[i]);
	cat(file=fh, " ", SampleNames[i], sep="");
}
cat(file=fh,"\n");
for(i in 1:num_samples){
	cat(file=fh, SampleNames[i], asFull[i,], sep=" ");
#	cat(file=fh, SampleNames[i], asFull[i,]);
	cat(file=fh, "\n");
}
close(fh);
