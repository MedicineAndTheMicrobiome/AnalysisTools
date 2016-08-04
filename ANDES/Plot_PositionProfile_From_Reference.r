#!/usr/bin/env Rscript

###############################################################################
#                                                                             # 
#       Copyright (c) 2009 J. Craig Venter Institute.                         #     
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

params=c(
	"rmsd_values", "i", 1, "character",
	"annotation", "a", 2, "character",
	"output_root", "o", 2, "character",
	"topn", "n", 2, "numeric",
	"above", "t", 2, "numeric",
	"noise_cutoff", "s", 2, "numeric"
)

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name, 
	"\n",
	"	-i <RMSD values, see below for format.>\n",
	"	[-a <annotation regions file name, 0-space-based coordinates>]\n",
	"	[-o <output filename root]\n",
	"	[-n <Top N RMSD values to label>]\n",
	"	[-t <RMSD threshold to label>]\n",
	"	[-s <Noise cutoff, below which to not draw>]\n",
	"\n",
	"\n",
	"Plots a RMSD from reference graph across the length of the sequence.\n",
	"\n",
	"The format of the input RMSD values file should be:\n\n",
	"1st line: \n",
	"	<reference_name>\\t<dataset1_name>\\t<dataset2_name>\\t...\n",
	"\n2nd to last line: \n",
	"	0.000\\t<dataset1_RMSD_1>\\t<dataset2_RMSD_1>\\t...\n",
	"	0.000\\t<dataset1_RMSD_2>\\t<dataset2_RMSD_2>\\t...\n",
	"	0.000\\t<dataset1_RMSD_3>\\t<dataset2_RMSD_3>\\t...\n",
	"	...\n",
	"	0.000\\t<dataset1_RMSD_L>\\t<dataset2_RMSD_L>\\t...\n",
	"\n",
	"\n",
	"The format of the annotation file (-a option), should be:\n",
	"	<begin>\\t<end>\\n\n",
	"	<begin>\\t<end>\\n\n",
	"	...\n",
	"	<begin>\\t<end>\\n\n",
	"\n\n",
	sep="");

###############################################################################

if(!length(opt$rmsd_values)){
	cat(usage);
	q(status=-1);
}

if(length(opt$rmsd_values)){
	RmsdFile=opt$rmsd_values;	
}

if(length(opt$output_root)){
	OutputRoot=opt$output_root;
}else{
	OutputRoot=opt$rmsd_values;
}

if(length(opt$annotation)){
	AnnotationFile=opt$annotation;
}else{
	AnnotationFile="";
}

if(length(opt$topn)){
	TopNLabel=opt$topn;
}else{
	TopNLabel=0;
}

if(length(opt$above)){
	AboveLabel=opt$above;
}else{
	AboveLabel="";
}

if(length(opt$noise_cutoff)){
	NoiseCutoff=opt$noise_cutoff;
}else{
	NoiseCutoff=0.0;
}

###############################################################################

cat("Input RMSD File: ", RmsdFile, "\n");
cat("Annotation File: ", AnnotationFile, "\n");
cat("Top N Points to Label: ", TopNLabel, "\n");
cat("Points above RMSD Threshold to label: ", AboveLabel, "\n");
cat("Noise Cutoff to not plot: ", NoiseCutoff, "\n");

###############################################################################

# Load the annotation file, if it's specified
if(AnnotationFile!=""){
	annotations=as.matrix(read.table(AnnotationFile, sep="\t"));
}

# Load specified dataset
A=read.table(RmsdFile, header=TRUE, check.names=FALSE)
num_sets=ncol(A[1,])

# Change this if you want a different output format
vsRefProfPDF=paste(OutputRoot, ".pdf", sep="")
pdf(vsRefProfPDF,width=11,height=5.25)

# Computes x axis offsets
len=length(A[,2])
x_axis=1:len;

# Initiate empty plot
y_max=.8;
plot(0,0, xlim=c(0,len), ylim=c(0,y_max), type="n",
		main="RMSD of Samples vs. Reference",
		xlab="Genomic Axis", ylab="RMSD", xaxt="n");
axis(1,at=seq(0,len, 100));

# Plot annotations in grey
if(AnnotationFile!=""){
	num_annotations=nrow(annotations);
	cat("Num annotations to draw: ", num_annotations, "\n");
	for(i in 1:num_annotations){
		rect(
			annotations[i, 1],
			0,
			annotations[i, 2],
			y_max,
			border=NA,
			col="grey85"
		);
	}
}

# Draw reference lines for max and biallelic differences
MAX_RMSD=.632;
BIALLELIC_RMSD=.316;
abline(h=MAX_RMSD,col="gray50",lty="dashed")
abline(h=BIALLELIC_RMSD,col="gray50",lty="dotted")

# Draw line for noise cutoff
abline(h=NoiseCutoff,col="gray50")

# Select a color for each dataset
if(num_sets <= 8){
	# crayola 8 pack
	palette=c("red", "blue", "green", "purple", "orange", "brown", "black", "yellow3");
}else{
	palette=rainbow(num_sets,s=1, v=.80);
}

# Draw the legend based on the headers in the RMSD file.  The first column is the reference.
for(x in 2:num_sets){
	# Plot the legend, so we know what colors go with which dataset
	legend(0, .91-(x/35), names(A)[x], bty="n", text.col=palette[x-1], cex=.8);
}

# Draw the points 
nsamples=ncol(A);
nposition=nrow(A);
for(pos in 1:nposition){
	# Randomly draw among the various samples, so the last sample drawn doesn't always cover the other samples' points
	for(s in sample(2:nsamples)){
		p=A[pos,s];
		if(p>NoiseCutoff){
			points(pos,p, col=palette[s-1], pch=1, cex=.5);
		}
	}
}

# Keep track of the labels
labelled_indices=numeric();
max_vec=apply(A, 1, max);

# Label the top points 
if(TopNLabel>0){
	top_order=order(max_vec, decreasing=T);
	top_idx=(head(top_order, TopNLabel));
	text(top_idx, max_vec[top_idx], labels=top_idx, cex=.4, pos=3, offset=.15 );
	labelled_indices=c(labelled_indices, top_idx);
}

# Label points above the specified threshold
if(AboveLabel!=""){
	top_idx=which(max_vec>AboveLabel);
	text(top_idx, max_vec[top_idx], labels=top_idx, cex=.4, pos=3, offset=.15 );
	labelled_indices=c(labelled_indices, top_idx);
}

###############################################################################

# Output the top points in a file
if(TopNLabel>0 || AboveLabel!=""){

	# Just keep the unique indices
	labelled_indices=sort(unique(labelled_indices));

	# Construct filename
	lab_pts_fname_orb=paste(OutputRoot, ".1-res.txt", sep="");
	lab_pts_fname_zsb=paste(OutputRoot, ".0-spc.txt", sep="");

	# Open files
	fh_orb=file(lab_pts_fname_orb, "w");
	fh_zsb=file(lab_pts_fname_zsb, "w");

	# Output each labelled point 
	num_labelled=length(labelled_indices);
	for(i in 1:num_labelled){
		cat(file=fh_orb, labelled_indices[i], max_vec[labelled_indices[i]], sep="\t");
		cat(file=fh_orb, "\n");
		cat(file=fh_zsb, labelled_indices[i]-1, labelled_indices[i], max_vec[labelled_indices[i]], sep="\t");
		cat(file=fh_zsb, "\n");
	}

	# Close files
	close(fh_orb);
	close(fh_zsb);
}

###############################################################################

# Draw a histogram
for(i in 2:nsamples){
	#hist(A[,i], xlim=c(0, MAX_RMSD), breaks=100, xlab="RMSD", main="Distribution of RMSDs");
	hist(log(A[,i], 10), breaks=100, xlab=expression(Log[10](RMSD)),  main="Distribution of RMSDs");
	if(AboveLabel!=""){
		abline(v=log(AboveLabel, 10), col="grey", lty="dotted");
	}
	mtext(names(A)[i], line=1, cex=.6, col=palette[i-1]);
}

###############################################################################

dev.off();
cat("Done.\n");

###############################################################################
