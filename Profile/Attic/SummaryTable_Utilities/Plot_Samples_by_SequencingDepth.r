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
	"	-i <input summary_table.tsv file>\n",
	"\n",
	"This script will read in a summary table and\n",
	"then generate some plots based on sequencing\n",
	"depth per sample.  Plots include MDS, distance\n",
	"to centroid, and tail statistic.\n",
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

###############################################################################
# Load data
in_mat=as.matrix(read.delim(InputFileName, sep="\t", header=TRUE, row.names=1, check.names=FALSE, comment.char="", quote=""))

#print(in_mat);
counts_mat=in_mat[, 2:ncol(in_mat)];

num_categories=ncol(counts_mat);
num_samples=nrow(counts_mat);

cat("Num Categories: ", num_categories, "\n", sep="");
cat("Num Samples: ", num_samples, "\n", sep="");

#print(countsMat)

category_names=colnames(counts_mat);
sample_names=rownames(counts_mat);

###############################################################################

cat("Normalizing counts...\n");

norm_mat=matrix(0, nrow=num_samples, ncol=num_categories);
rownames(norm_mat)=rownames(counts_mat);
colnames(norm_mat)=colnames(counts_mat);

sample_counts=apply(counts_mat, 1, sum);
for(i in 1:num_samples){
	norm_mat[i,]=counts_mat[i,]/sample_counts[i];
}	

#print(norm_mat);

###############################################################################

# Compute distances
cat("Computing distances...\n");
dist_mat=dist(norm_mat);
#dist_mat=vegdist(norm_mat, "bray");
#print(dist_mat);

pdf(paste(OutputFileRoot, ".count_analysis.pdf", sep=""), width=8.5, height=11);

mds_pts=cmdscale(dist_mat, k=2);
#isoRes=isoMDS(dist_mat, k=2);
#mds_pts=isoRes$points;

overall_centroid=c(mean(mds_pts[,1]), mean(mds_pts[,2]));

cat("\n");
cat("Overall Centroid: \n");
cat("x = ", overall_centroid[1], "\ty = ", overall_centroid[2], "\n");
cat("\n");

#-------------------------------------------------------------------------
# Coloring, don't forget to update the legend if you change this
coloring=rep("black", num_samples);

cutoffs=c(Inf,5000,750,500,250,100);
colors=c("black", "purple", "blue", "green", "orange", "red");
num_cutoffs=length(cutoffs);

coloring=character(num_samples);
for(i in 1:num_cutoffs){
	coloring[sample_counts<cutoffs[i]]=colors[i];
}

plot_legend=function(xrange, yrange){
	xincr=diff(xrange)/10;
	yincr=diff(yrange)/10;
	info=sprintf("< %g", cutoffs);
	fill=colors;
	legend(xrange[2]-xincr, yrange[2]-5*yincr, legend=info, fill=fill, cex=.5);
}
#-------------------------------------------------------------------------

# Points Plot
plot(mds_pts[,1], mds_pts[,2], xlab="Dim 1", ylab="Dim 2", main="All Data Points MDS", col=coloring);
plot_legend(range(mds_pts[,1]), range(mds_pts[,2]));

# Label Plot
plot(mds_pts[,1], mds_pts[,2], xlab="Dim 1", ylab="Dim 2", main="All Data Points MDS", type="n");
plot_legend(range(mds_pts[,1]), range(mds_pts[,2]));
text(mds_pts[,1], mds_pts[,2], labels=sample_names, col=coloring, cex=.4);
plot_legend(range(mds_pts[,1]), range(mds_pts[,2]));

sqr_dist_mat=as.matrix(dist_mat);

#print(sqr_dist_mat);
mean_dist_to_other_points=apply(sqr_dist_mat, 1, function(x){return(sum(x)/(length(x)-1))});

# Plot counts vs distance to centroid
plot(sample_counts, mean_dist_to_other_points, xlab="Reads/Sample", ylab="Mean Distance to Other Samples", 
	main="Difference vs. Sequencing Depth", col=coloring);
plot_legend(range(sample_counts), range(mean_dist_to_other_points));

plot(sample_counts, mean_dist_to_other_points, xlab="Reads/Sample", ylab="Mean Distance to Other Samples",
	main="Difference vs. Sequencing Depth",  type="n");
text(sample_counts, mean_dist_to_other_points, labels=sample_names, col=coloring, cex=.4);
plot_legend(range(sample_counts), range(mean_dist_to_other_points));


###############################################################################


#mem=c("p32d1-H3N2_136", "p2d1-H3N2_58");
#mem=c("p32d7-H3N2_140", "p20d1-H1N1_45");
#mem=c("hc9d28-CTRL_54", "p9d3-H3N2_84"); # artifact of MDS
#mini_mat=counts_mat[mem,];
#tot=apply(mini_mat, 2, sum);
#print(mini_mat[,tot>0]);

tail_statistic=function(x){
        sorted=sort(x, decreasing=TRUE);
        norm=sorted/sum(x);
        n=length(norm);
        tail=0;
        for(i in 1:n){
                tail=tail + norm[i]*((i-1)^2);
        }
        return(sqrt(tail));
}

tail_stat=numeric();
for(i in 1:num_samples){
	tail_stat[i]=tail_statistic(norm_mat[i,]);
}	

plot(sample_counts, tail_stat, xlab="Reads/Sample", ylab="Tail Statistic", main="Diversity vs. Sequencing Depth", col=coloring);
plot_legend(range(sample_counts), range(tail_stat));

plot(sample_counts, tail_stat, xlab="Reads/Sample", ylab="Tail Statistic", main="Diversity vs. Sequencing Depth", type="n");
text(sample_counts, tail_stat, labels=sample_names, col=coloring, cex=.4);
plot_legend(range(sample_counts), range(tail_stat));

par(mfrow=c(3,1));
hist(sample_counts, xlab="Reads/Sample", main="Distribution of Read Depth", breaks=seq(0, max(sample_counts)+50, 50));
abline(v=cutoffs, col=colors);


hist(tail_stat, xlab="Tail Statistic", main="Distribution of Diversity", breaks=20);
hist(mean_dist_to_other_points, xlab="Mean Dist to Other Samples", main="Distribution of Differences", breaks=20);

###############################################################################

dev.off();


