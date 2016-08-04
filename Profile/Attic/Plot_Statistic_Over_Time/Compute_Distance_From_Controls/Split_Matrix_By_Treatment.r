#!/usr/bin/env Rscript

###############################################################################

library(MASS)
library('getopt');

params=c(
	"treatment_info_file", "t", 1, "character",
	"treatment_column", "c", 1, "numeric",
	"distance_matrix", "d", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-t <treatment information file, tab-seperated>\n",
	"	-c <column to split distance matrix by>\n",
	"\n",
	"	-d <distance matrix>\n",
	"\n",
	"This script will generate a distances from the baseline based on\n",
	"the specified distance matrix.  The samples are grouped.\n",
	"\n");

if(
	!length(opt$treatment_info_file)||
	!length(opt$treatment_column)||
	!length(opt$distance_matrix)
){
	cat(usage);
	q(status=-1);
}

TreatmentFileName=opt$treatment_info_file;
TreatmentColumn=opt$treatment_column;
DistanceMatrixFileName=opt$distance_matrix;

##############################################################################

cat("\n");
cat("Treatment Filename: ", TreatmentFileName, "\n");
cat("Distance Matrix Filename: ", DistanceMatrixFileName, "\n");
cat("\n");

##############################################################################

# Load Treatment Info
treatment_matrix=as.matrix(read.table(TreatmentFileName, sep="\t", header=T, row.names=1));
treatment_col_names=colnames(treatment_matrix);

treatment_colname=treatment_col_names[TreatmentColumn];
treatment_matrix[,TreatmentColumn]=gsub(" ", "_", treatment_matrix[,TreatmentColumn]);
treatment_values=sort(unique(treatment_matrix[,TreatmentColumn]));

cat("Treatments:\n");
print(treatment_values);
cat("\n");

# Load Distance Matrix
distance_matrix=as.matrix(read.table(DistanceMatrixFileName, sep=" ", header=T, row.names=1, check.names=F));
#print(distance_matrix);

##############################################################################

for(trt in treatment_values){
	cat("Working on: ", trt, "\n");
	
	sample_idx=which(treatment_matrix[,TreatmentColumn]==trt);
	sample_names=names(sample_idx);
	print(sample_names);

	trt_distmat=distance_matrix[sample_names, sample_names];
	#print(trt_distmat);

	# Write out distance matrix
	filename=paste(gsub("r_distmat","", DistanceMatrixFileName), trt, ".distmat", sep="");
	cat("Writing: ", filename, "\n");
	fc=file(filename, "w");
	cat(file=fc, " ");
	cat(file=fc, sample_names);
	cat(file=fc, "\n");
	for(samp_name in sample_names){
		cat(file=fc, samp_name, paste(trt_distmat[samp_name,], collapse=" "));
		cat(file=fc, "\n"); 
	}
	close(fc);
	
	cat("\n");
}

##############################################################################

cat("Done.\n")

q(status=0)
