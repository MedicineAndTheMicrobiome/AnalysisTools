#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"export_clusters", "c", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"\t\t-i <input summary_table.xls file>\n",
	"\t\t[-c (export clusters)]\n",
	"\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileNameRoot=opt$input_file;
ExportClusters=!is.null(opt$export_clusters);

if(ExportClusters){
        cat("Exporting clusters.\n");
}else{
        cat("Not exporting clusters.\n");
}

###############################################################################
# Load counts from file

cat("Working on ", InputFileName, "\n", sep="");

# Load summary_table.xls
mat<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, row.names=1))

# Exclude total counts column
count_mat=mat[,2:ncol(mat)]
num_cat=ncol(count_mat);

#cat("Input matrix:\n");
#print(count_mat);

#------------------------------------------------------------------------------
# Get sample names
sample_names=rownames(mat);

# Get column/category names
categories=as.vector(colnames(count_mat));
num_categories=length(categories);

#cat("Original category names:\n");
#print(categories);

# Compute shorted names
short_names=character(num_categories);
for(i in 1:num_categories){
	taxonomy=unlist(strsplit(categories[i], " "));
	short_names[i]=taxonomy[length(taxonomy)];
}

#cat("Shortened category names:\n");
#print(short_names);

#------------------------------------------------------------------------------
# Identify 0 count columns
col_sum=apply(count_mat, 2, sum);
col_idx_iszero=(1:num_cat)[col_sum==0];
if(length(col_idx_iszero)>0){
	cat("Zero count categories:\n");
	print(col_idx_iszero);
}else{
	cat("No zero count categories detected.\n");
}

if(length(col_idx_iszero)>0){
	cat("Zero count names.\n");
	print(short_names[col_idx_iszero]);
}

# Remove all zero columns
if(length(col_idx_iszero)>0){
	cat("Removing zero count categories.\n");
	short_names=short_names[-col_idx_iszero];
	count_mat=count_mat[,-col_idx_iszero];
}

#------------------------------------------------------------------------------
# Show num samples/categories to be used

NumSamples=nrow(count_mat);
NumCategories=ncol(count_mat);

cat("\n");
cat("Num Samples: ", NumSamples, "\n");
cat("Num Non-zero Categories: ", NumCategories, "\n");
cat("\n");

#------------------------------------------------------------------------------
# Normalize
sample_totals=numeric(NumSamples);
prob_mat=matrix(nrow=NumSamples,ncol=NumCategories);
colnames(prob_mat)=short_names;
rownames(prob_mat)=sample_names;
for(i in 1:NumSamples){
    	sample_totals[i]=sum(count_mat[i,]);
	prob_mat[i,]=count_mat[i,]/sample_totals[i];
}
#print(prob_mat);

#------------------------------------------------------------------------------
# Sort Matrix by sample names (so they will be in a consistent order later)
sample_names=rownames(count_mat);
sample_name_sort=sort(sample_names, index.return=TRUE);
mat_sort=count_mat[sample_name_sort$ix,];

###############################################################################
###############################################################################

source("WeightedRankDifference.r");

###############################################################################

library(MASS);
library(stats);
library(vegan);

###############################################################################

pdf(paste(OutputFileNameRoot, ".cstop_clustered.pdf", sep=""), height=8.5, width=11);

dist_list=list();

countmat=prob_mat*10000;
intmat=apply(countmat,2,as.integer);
rownames(intmat)=rownames(prob_mat);
#print(intmat);

dist_list[[1]]=weight_rank_dist(prob_mat, 4);
dist_list[[2]]=dist(prob_mat);
dist_list[[3]]=vegdist(intmat, method="morisita");
dist_list[[4]]=vegdist(intmat, method="bray");

dist_name=vector();
dist_name[1]="Weighted Rank Difference, power=4";
dist_name[2]="Euclidean";
dist_name[3]="Morisita-Horn";
dist_name[4]="Bray-Curtis";

num_distances=length(dist_list);
#num_distances=1;


###############################################################################
###############################################################################

source("cstop");

cat("Number of clustering techniques (distances) to try: ", num_distances, "\n");

laymat=matrix(c(1,1,2,3), nrow=2, byrow=T);
layout(laymat);

for(dist_idx in 1:num_distances){

	cat("Working on: ", dist_name[dist_idx], "\n", sep="");
	samp_dist=dist_list[[dist_idx]];
	max_dist=max(samp_dist);

	# Remove 0 distances so isoMDS doesn't freakout
	for(i in 1:length(samp_dist)){
		if(samp_dist[i]==0){
			samp_dist[i]=1e-323;
		}
	}
	#cat("Removed 0 distances.\n");

	# Compute sample clusters
	hcluster=hclust(samp_dist, method="ward");
	clusters=cstop(T=hcluster, D=samp_dist, pvalue=.05);

	if(ExportClusters){
                cluster_out_fn=paste(OutputFileNameRoot, ".recommended_cstop_clusters.", gsub(" ", "_", dist_name[dist_idx]), ".tsv", sep="");
                cluster_out_fh=file(cluster_out_fn, "wt");
                matout=cbind(clusters, names(clusters));
                matout=matout[order(matout[,1]),];
                for(sample in 1:NumSamples){
                        cat(file=cluster_out_fh, matout[sample,], "\n");
                }
        }

	# Sort the results from the clustering by sample name (original order was by cluster id)
	sample_names=names(clusters);
	cluster_assignments=as.vector(clusters);
	sorted_names=sort(sample_names, index.return=TRUE);
	cluster_assignments=cluster_assignments[sorted_names$ix];
	clusters=clusters[sorted_names$ix];

	# Create mapping from name to cluster assignment
	color_list=list();
	for(cl in 1:NumSamples){
		color_list[[sorted_names$x[cl]]]=cluster_assignments[cl];
	}

	#------------------------------------------------------------------------------
	
	# Compute label scales
	label_scale=32/NumSamples;
	label_scale=ifelse(label_scale>1.2, 1.2, label_scale);

	# Plot dendrogram
	denfun=function(n){
		if(is.leaf(n)){
			leaf_attr=attributes(n);
			leaf_name=leaf_attr$label;
			attr(n, "nodePar") = c(leaf_attr$nodePar, list(lab.col=color_list[[leaf_name]]));
		}
		return(n);
	}
	hclust_res=as.dendrogram(hclust(samp_dist, method="ward"));
	DL=dendrapply(hclust_res, denfun);
	plot(DL, cex=label_scale*50, main=dist_name[dist_idx]);

	# Compute MDS 
	#print(samp_dist);
	#mds=isoMDS(samp_dist);

	# Compute display ranges/sizes
	#min=min(mds$points[,1])
	#max=max(mds$points[,1])
	#width=(max-min)
	#margin=width*(0.1)

	# Plot sample name/text
	#plot(mds$points,type="n", xlim=c(min-margin, max+margin), xlab="", ylab="")
	#text(mds$points,labels=sorted_names$x, cex=label_scale, col=clusters)

	# Plot points
	#cat("Plotting points...\n");
	#plot(mds$points,type="p", xlim=c(min-margin, max+margin), xlab="", ylab="", col=clusters, pch=19)

}

#-------------------------------------------------------------------------------

cat("Done.\n");
q(status=0);
