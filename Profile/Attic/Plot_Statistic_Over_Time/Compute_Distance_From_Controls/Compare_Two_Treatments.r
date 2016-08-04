#!/usr/bin/env Rscript

###############################################################################

library(MASS)
library('getopt');

params=c(
	"treatment_info_file", "t", 1, "character",
	"treatment_column", "c", 1, "numeric",
	"distance_matrix", "d", 1, "character",
	"output_filename_root", "o", 1, "character",
	"mds_range", "r", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-t <treatment information file, tab-seperated>\n",
	"	-c <column to split samples by>\n",
	"\n",
	"	-d <distance matrix>\n",
	"	-o <output filename root>\n",
	"\n",
	"	[-r <range for plotting mds, eg. 2 will be +/- 2 in both dimension, default dynamic scaled>]\n",
	"This script will compare two groups of samples, based on the treatment column\n",
	"\n");

if(
	!length(opt$treatment_info_file)||
	!length(opt$treatment_column)||
	!length(opt$distance_matrix) ||
	!length(opt$output_filename_root)
){
	cat(usage);
	q(status=-1);
}

TreatmentFileName=opt$treatment_info_file;
TreatmentColumn=opt$treatment_column;
DistanceMatrixFileName=opt$distance_matrix;
OutputFilenameRoot=opt$output_filename_root;
if(!length(opt$mds_range)){
	MDSRange=0;
}else{
	MDSRange=opt$mds_range;
}

##############################################################################

path_comp=strsplit(script_name, "/")[[1]];
bin_comp_idx=length(path_comp);
bin_path=paste(path_comp[-bin_comp_idx], collapse="/", sep="");
if(nchar(bin_path)==0){
        bin_path=".";
}
cat("Binary path: '", bin_path, "'\n", sep="");

source(paste(bin_path, "wilcoxanova.r", sep="/"));

##############################################################################

cat("\n");
cat("Treatment Filename: ", TreatmentFileName, "\n");
cat("Distance Matrix Filename: ", DistanceMatrixFileName, "\n");
cat("Output Filename Root: ", OutputFilenameRoot, "\n");
cat("Range? ", MDSRange, "\n");
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
samples_in_matrix=colnames(distance_matrix);
#print(distance_matrix);

##############################################################################

pdf(paste(OutputFilenameRoot, ".pdf", sep=""), height=8.5, width=17);

fc=file(paste(OutputFilenameRoot, ".csv", sep=""), "w");
cat(file=fc, "# distmat,comparison,num_trt1,num_trt2,mdn_intra,mdn_inter,pval,med_trt1,med_trt2\n");

num_trtments=length(treatment_values);
for(trt_ix1 in 1:num_trtments){
	trt1=treatment_values[trt_ix1];
	samples1=which(treatment_matrix[,TreatmentColumn]==trt1);
	available_samples1=intersect(names(samples1), samples_in_matrix);

	for(trt_ix2 in 1:num_trtments){

		if(trt_ix1>=trt_ix2){
			next;
		}

		trt2=treatment_values[trt_ix2];
		cat("Trt: ", trt1, " vs ", trt2, "\n", sep="");

		samples2=which(treatment_matrix[,TreatmentColumn]==trt2);
		available_samples2=intersect(names(samples2), samples_in_matrix);

		cat("Samples: ", trt1, "\n", sep="");
		print(available_samples1);
		cat("Samples: ", trt2, "\n", sep="");
		print(available_samples2);	

		#print(distance_matrix);
		all_relevant_samples=c(available_samples1, available_samples2);
		
		sub_distmat=distance_matrix[all_relevant_samples, all_relevant_samples];
		results=wilcoxanova(sub_distmat, available_samples1, available_samples2);
		print(results);

		# Output statistics into text file
		cat(file=fc, c(
			DistanceMatrixFileName,
			paste(trt1, " vs ", trt2, sep=""),
			length(available_samples1),
			length(available_samples2),
			results$intra_med,
			results$inter_med,
			results$pval,
			results$grpA_med,
			results$grpB_med
			), sep=","
		);
		cat(file=fc, "\n");

		# Compute MDS and plot it
		distmat_wcentroid=add_centroids_to_distmat(sub_distmat, available_samples1, available_samples2);
		#print(distmat_wcentroid);

		t1_col="blue";
		t2_col="green";

		mds=isoMDS(as.dist(distmat_wcentroid));
		colors=rep(t1_col, length(all_relevant_samples));
		names(colors)=all_relevant_samples;
		colors[available_samples2]=t2_col;
		colors["cA"]=t1_col;
		colors["cB"]=t2_col;

		par(mfrow=c(1,2))

		num_samples=length(all_relevant_samples);

		rotated=level_rotate(
			c(mds$points["cA",1], mds$points["cA",2]),
			c(mds$points["cB",1], mds$points["cB",2]),
			cbind(mds$points[,1],mds$points[,2])		
		);

		# Plot points with labels scaled
		plot(rotated[1:num_samples,1], rotated[1:num_samples,2], 
			col=colors, cex=.5, main=DistanceMatrixFileName);
		# Plot labels
		text(rotated[1:num_samples,1], rotated[1:num_samples,2], 
			labels=rownames(mds$points)[1:num_samples], pos=3, col=colors, cex=.5);


		# Plot only points and centroid
		if(MDSRange>0){
			plot(rotated[,1], rotated[,2], 
				col=colors, cex=c(rep(.5, num_samples),3,3), main=DistanceMatrixFileName, 
				ylab="Dimension 1", xlab="Dimension 2",
				ylim=c(-MDSRange, MDSRange), xlim=c(-MDSRange, MDSRange));
		}else{
			plot(rotated[,1], rotated[,2], 
				col=colors, cex=c(rep(.5, num_samples),3,3), main=DistanceMatrixFileName, 
				ylab="Dimension 1", xlab="Dimension 2");
		}
		text(rotated["cA",1], rotated["cA",2], labels=trt1, col=t1_col);
		text(rotated["cB",1], rotated["cB",2], labels=trt2, col=t2_col);


		
	}	
}

##############################################################################

cat("Done.\n")

q(status=0)
