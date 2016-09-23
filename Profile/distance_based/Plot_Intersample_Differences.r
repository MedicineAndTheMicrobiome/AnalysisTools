#!/usr/bin/env Rscript

###############################################################################

library(MASS)
library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 2, "character",
	"use_wrd", "w", 2, "logical",
	"split_char", "s", 2, "character",
	"colors", "c", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	[-o <output file root name>]\n",
	"	[-w (flag: use Weighted Rank Distance as distance, default is Euclidean)]\n",
	"	[-s <short name split character, default is \" \">]\n",
	"	[-c <num colors/groups for preliminary groups, default=4>]\n",
	"\n",
	"	Reads in a summary_table.xls file and generates the follow output files:\n",
	"		1.) MDS Plot\n",
	"		2.) Dendrogram\n",
	"		3.) Heat Map\n",
	"		4.) Distance Matrix\n",
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
	OutputFileRoot=gsub("\\.summary_table\\.tsv$", "", OutputFileRoot);
	OutputFileRoot=gsub("\\.summary_table\\.xls$", "", OutputFileRoot);
	cat("No output file root specified.  Using input file name as root.\n");
}

UseWeightedRankDifference=!is.null(opt$use_wrd);
if(UseWeightedRankDifference){
	cat("Using Weighted Rank Difference as distance.\n");
	type="wrd";
}else{
	cat("Using Euclidean as distance.\n");
	type="euc";
}

SNSplitChar=" ";
if(length(opt$split_char)){
	SNSplitChar=opt$split_char;
}
cat("Spliting long names by: \"", SNSplitChar, "\"\n", sep="");


PrelimColors=4;
if(length(opt$colors)){
	PrelimColors=opt$colors;
}
#cat("Num colors/groups to detect: ", PrelimColors, "\n");

###############################################################################
# Load the WRD code

#path_comp=strsplit(script_name, "/")[[1]];
#bin_comp_idx=length(path_comp);
#bin_path=paste(path_comp[-bin_comp_idx], collapse="/", sep="");
#if(nchar(bin_path)==0){
#        bin_path=".";
#}
#cat("Binary path: '", bin_path, "'\n", sep="");
#
#source(paste(bin_path, "Cluster/WeightedRankDifference.r", sep="/"));

order_dist=function(a, b, deg){
        #asort=sort(a, decreasing=F, index.return=T, method="shell");
        #bsort=sort(b, decreasing=F, index.return=T, method="shell");

        arank=rank(a, ties.method="average");
        brank=rank(b, ties.method="average");

        sort_sqrdiff=sqrt(sum(((arank-brank)^2)*((a+b)/2)^deg));
        #sort_sqrdiff=sqrt(sum(((arank-brank)^2)*(((a-b)/2)^2)));
        return(sort_sqrdiff);

}

###############################################################################

weight_rank_dist=function(M, deg){
        NumSamples=nrow(M);
        order_dist_mat=matrix(0, nrow=NumSamples, ncol=NumSamples);
        for(i in 1:NumSamples){
                for(j in 1:i){
                        order_dist_mat[i,j]=order_dist(M[i,], M[j,], deg);
                }
        }
        rownames(order_dist_mat)=rownames(M);
        return(as.dist(order_dist_mat));
}

###############################################################################

OutputPDF = paste(OutputFileRoot, ".", type, ".pdf", sep="");
DistanceMatrixTXT = paste(InputFileName, ".", type, ".r_distmat", sep="");

cat("Output PDF file name: ", OutputPDF, "\n", sep="");
cat("Output Distance Matrix: ", DistanceMatrixTXT, "\n", sep="");

###############################################################################
# Load data
InMat=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, row.names=1, check.names=FALSE, comment.char=""))

num_categories=ncol(InMat)-1;
num_samples=nrow(InMat);
if(num_samples<PrelimColors){
	PrelimColors=num_samples;
}

cat("Num colors/groups to detect: ", PrelimColors, "\n");

cat("Num Categories: ", num_categories, "\n", sep="");
cat("Num Samples: ", num_samples, "\n", sep="");

countsMat=InMat[,2:(num_categories+1)];
#print(countsMat)

Categories=colnames(countsMat);
SampleNames=rownames(countsMat);

#------------------------------------------------------------------------------

# Generate short names
ShortCat=character(num_categories);
for(i in 1:num_categories){
	taxonomy=unlist(strsplit(Categories[i], SNSplitChar));
	ShortCat[i]=tail(taxonomy,1);
	#print(ShortCat[i]);
}

# Sum up the number of members in each column by doing an outer product
sample_counts=apply(countsMat, 1, sum);
print(sample_counts);

# Normalize counts
cat("Normalizing counts...\n");
normalized=matrix(0, nrow=num_samples, ncol=num_categories);
rownames(normalized)=rownames(countsMat);
colnames(normalized)=colnames(countsMat);
for(i in 1:num_samples){
	normalized[i,]=countsMat[i,]/sample_counts[i];
}	

# Compute distances
if(type=="euc"){
	dist_mat=dist(normalized);
}else if (type=="wrd"){
	dist_mat=weight_rank_dist(normalized, deg=4);
}
print(dist_mat);

# Scale the labels sizes so they don't overlap each other.
# But limit the scaling so it's not ridiculously large.
rlabel_scale=42/num_samples;
if(rlabel_scale > 1){
	rlabel_scale = 1;
}
clabel_scale=42/num_categories;
if(clabel_scale > 1){
	clabel_scale = 1;
}

###############################################################################
###############################################################################

pdf(OutputPDF,width=8.5,height=11)

###############################################################################
# Plot dendrograms
ward_clust=hclust(dist_mat, "ward");
orig_cex=par()$cex.lab;
par(cex=rlabel_scale);
par(cex.main=1.2*(1/rlabel_scale))
par(cex.sub=1/rlabel_scale)
plot(ward_clust, main=InputFileName, xlab="Ward's Minimum Variance");

par(mfrow=c(2,1));

par(cex=rlabel_scale);
par(cex.main=1.2*(1/rlabel_scale))
par(cex.sub=1/rlabel_scale)
plot(hclust(dist_mat, "complete"), xlab="Complete Linkage")

par(cex=rlabel_scale);
par(cex.main=1.2*(1/rlabel_scale))
par(cex.sub=1/rlabel_scale)
plot(hclust(dist_mat, "average"), xlab="Average Linkage")

par(cex=1);
par(cex.main=1.2);
par(cex.sub=1);

###############################################################################
# Draw MDS Plot

# Remove 0 distances with very small number
for(i in 1:length(dist_mat)){
	if(dist_mat[i]==0){
		dist_mat[i]=1e-323;
	}
}

# Assign some colors so we can track the points between dimensions
rainbow_colors=rainbow(PrelimColors, start=0, end=0.65);
clusters=cutree(ward_clust, k=PrelimColors);

#------------------------------------------------------------------------------
# Compute the padding inside of the graph so the text labels don't fall off the screen

# Compute 2D MDS
mds=isoMDS(dist_mat);

xmin=min(mds$points[,1])
xmax=max(mds$points[,1])
xwidth=(xmax-xmin)
xmargin=xwidth*(0.1)

ymin=min(mds$points[,2])
ymax=max(mds$points[,2])
ywidth=(ymax-ymin)
ymargin=ywidth*(0.1)

# Plot with names
par(mfrow=c(1,1));
plot(0, 0, type="n", main=InputFileName, 
	xlim=c(xmin-xmargin, xmax+xmargin),
	ylim=c(ymin-ymargin, ymax+ymargin),
	 xlab="", ylab="")
text(mds$points[,1], mds$points[,2], labels=SampleNames, cex=rlabel_scale, col=rainbow_colors[clusters])

# Plot points
#plot(mds$points[,1], mds$points[,2], main=InputFileName, 
#	xlim=c(xmin-xmargin, xmax+xmargin),
#	ylim=c(ymin-ymargin, ymax+ymargin),
#	xlab="", ylab="",  col=rainbow_colors[clusters])

#------------------------------------------------------------------------------
# Plot 3D MDS

transform_range=function(x, new_min, new_max){

        x_min=min(x);
        x_max=max(x);
        x_toZero=x-x_min;

        x_range=x_max-x_min;
        x_norm=x_toZero/x_range;
        new_range=new_max-new_min;

        new_x=x_norm*new_range;
        new_x=new_x+new_min;

        if(all(is.nan(new_x))){
                new_x=rep(1,length(new_x));
        }

        return(new_x);
}

if(num_samples > 3){

mds3D=isoMDS(dist_mat, k=3);

par(mfrow=c(2,2));
plot(mds$points[,1], mds$points[,2],main="2D X/Y", xlab="x", ylab="y", col=rainbow_colors[clusters])

depth=transform_range( mds3D$points[,3], .5, 2.5);
plot(mds3D$points[,1], mds3D$points[,2], main="3D X/Y", xlab="x", ylab="y", cex=depth, col=rainbow_colors[clusters])
depth=transform_range( mds3D$points[,1], .5, 2.5);
plot(mds3D$points[,2], mds3D$points[,3], main="3D Y/Z", xlab="y", ylab="z", cex=depth, col=rainbow_colors[clusters])
depth=transform_range( mds3D$points[,2], .5, 2.5);
plot(mds3D$points[,1], mds3D$points[,3], main="3D X/Z", xlab="x", ylab="z", cex=depth,  col=rainbow_colors[clusters])

}
###############################################################################
###############################################################################
# Draw Dendrogram Heatmap Plot
 
par(mfrow=c(1,1));

max_sample_name_length=max(length(SampleNames));
max_category_name_length=max(length(ShortCat));

cat("Max Sample Name Length: ", max_sample_name_length, "\n");
cat("Max Category Name Length: ", max_category_name_length, "\n");

# Changing samples margin
#sample_margin=max_sample_name_length*rlabel_scale/25;
sample_margin=7;
category_margin=max_category_name_length*clabel_scale/30;

if(sample_margin<1){sample_margin=1};
if(category_margin<0){category_margin=1};

cat("Sample (Right) Margin: ", sample_margin, "\n");
cat("Category (Bottom) Margin: ", category_margin, "\n");


# Take the top 100 only before generating heat map
TOP=100;
TOP=min(num_categories, TOP);
avg_norm=apply(normalized, 2, mean);
sort_ix=order(avg_norm, decreasing=T);
sort_normalized=normalized[,sort_ix];
top_normalized=sort_normalized[,1:TOP];

hm_clabel_scale=42/TOP;
if(hm_clabel_scale > 1){
	clabel_scale = 1;
}

heatmap(top_normalized, cexRow=rlabel_scale,  cexCol=hm_clabel_scale,
	labCol=colnames(top_normalized),
	#xlab="Taxonomies", ylab="Samples",
	col=rev(rainbow(2^16, start=0, end=0.65)),
	margins=c(category_margin, sample_margin),
	main=sprintf("Top %i Taxa", TOP)
	)

##############################################################################
##############################################################################
# Output distance matrix

asFull=as.matrix(dist_mat);
fh=file(DistanceMatrixTXT, "w");
cat(file=fh, " ", paste(SampleNames, collapse=" "), "\n", sep="");
for(i in 1:num_samples){
	cat(file=fh, SampleNames[i], asFull[i,], sep=" ");
	cat(file=fh, "\n");
}
close(fh);

##############################################################################

cat("Done.\n")
dev.off();

q(status=0)
