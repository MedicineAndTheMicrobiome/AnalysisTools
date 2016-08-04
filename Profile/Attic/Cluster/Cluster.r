#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"no_bootstrap", "b", 2, "logical",
	"export_clusters", "c", 2, "logical",
	"export_fstats", "f", 2, "logical",
	"export_delta_fstats", "d", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n\t\t<Input summary_table.xls FileName>\n\n",
	"	-i <input summary_table.xls file>\n",
	"	[-b (turn off bootstrapping)]\n",
	"	[-c (export clusters)]\n",
	"	[-f (export fstats)]\n",
	"	[-d (export delta f-stat]\n",
	"\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileNameRoot=opt$input_file;
BootstrapOn=is.null(opt$no_bootstrap);
ExportClusters=!is.null(opt$export_clusters);
ExportFstats=!is.null(opt$export_fstats);
ExportDeltaFstats=!is.null(opt$export_delta_fstats);

if(BootstrapOn){
	cat("Bootstrapping is ON\n");
}else{
	cat("Bootstrapping is OFF\n");
}

if(ExportClusters){
	cat("Exporting clusters.\n");
}else{
	cat("Not exporting clusters.\n");
}

if(ExportFstats){
	cat("Exporting F-statistics.\n");
}else{
	cat("Not exporting F-statistics.\n");
}

if(ExportDeltaFstats){
	cat("Exporting Delta F-statistics.\n");
}else{
	cat("Not exporting Delta F-statistics.\n");
}


# Figure out where we are running, this script so we can find the "WeightedRankDifference.r" code.
path_comp=strsplit(script_name, "/")[[1]];
bin_comp_idx=length(path_comp);
bin_path=paste(path_comp[-bin_comp_idx], collapse="/", sep="");
if(nchar(bin_path)==0){
	bin_path=".";
}
cat("Binary path: '", bin_path, "'\n", sep="");

source(paste(bin_path, "WeightedRankDifference.r", sep="/"));

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
cat("Num Categories: ", NumCategories, "\n");
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

###############################################################################
###############################################################################

library(MASS);
library(stats);
library(vegan);

###############################################################################

pdf(paste(OutputFileNameRoot, ".clustered.pdf", sep=""), height=8.5, width=11);

dist_list=list();

countmat=prob_mat*10000;
intmat=apply(countmat,2,as.integer);
rownames(intmat)=rownames(prob_mat);
#print(intmat);

#dist_list[[1]]=weight_rank_dist(prob_mat, 2);
dist_list[[1]]=weight_rank_dist(prob_mat, 4);
#dist_list[[2]]=dist(prob_mat);
#dist_list[[3]]=vegdist(intmat, method="morisita");

dist_name=vector();
#dist_name[1]="Weighted Rank Difference, power=2";
dist_name[1]="Weighted Rank Difference, power=4";
#dist_name[3]="Weighted Rank Difference, power=6";
#dist_name[2]="Euclidean";
#dist_name[3]="Morisita-Horn";

#num_distances=1;
num_distances=length(dist_list);


###############################################################################
###############################################################################

compute_sum_sqr=function(a_idc, b_idc, distmat){
	#print(a_idc);
	#print(b_idc);
	
	sum_sqr=0;
	num_sums=0;
	for(a_idx in a_idc){
		for(b_idx in b_idc){
			if(a_idx!=b_idx){
				sum_sqr=sum_sqr+distmat[a_idx,b_idx]^2;
				num_sums=num_sums+1;
			}
		}
	}
	
	result=list();
	result$ssqr=sum_sqr;
	result$nsum=num_sums;
	return(result);
}

#-----------------------------------------------------------------------------

computePseudoFstat=function(clusters, samp_dist){

	#print(clusters);
	#print(samp_dist);

	cluster_ids=sort(unique(clusters));
	num_clusters=length(cluster_ids);

	#cat("Cluster Ids: ", paste(cluster_ids, collapse=", "), "\n", sep="");
	#cat("Num Clusters: ", num_clusters, "\n", sep="");

	distmat=as.matrix(samp_dist);
	num_samples=nrow(distmat);

	# Compute between/inter cluster variance 
	ss_sum=0;
	num_sum=0;
	for(i in 1:num_clusters){
		#cat("Working on ", i, " vs all.\n");
		i_idx=which(clusters==i);
		other_cl_idx=(1:num_samples)[-i_idx];
		ss_result=compute_sum_sqr(i_idx,other_cl_idx, distmat);
		ss_sum = ss_sum + ss_result$ssqr;
		num_sum = num_sum + ss_result$nsum;
	}
	intercluster_sqr_dist=ss_sum;
	#cat("Sum Intercluster distance: ", intercluster_sqr_dist, "\n");
	#cat("Num Intercluster distances:", num_sum, "\n");

	# Compute within/intra cluster variance
	ss_sum=0;
	num_sum=0;
	for(i in 1:num_clusters){
		#cat("Working on ", i, " vs self.\n");
		i_idx=which(clusters==i);
		ss_result=compute_sum_sqr(i_idx, i_idx, distmat);
		ss_sum = ss_sum + ss_result$ssqr;
		num_sum = num_sum + ss_result$nsum;
	}
	intracluster_sqr_dist=ss_sum;
	#cat("Sum Intracluster distance: ", intracluster_sqr_dist, "\n");
	#cat("Num Intracluster distances:", num_sum, "\n");

	k=num_clusters;
	n=num_samples;
	BGSS=intercluster_sqr_dist;
	WGSS=intracluster_sqr_dist;
	fstat=(BGSS/(k-1))/(WGSS/(n-k));

	#cat("k = num clusters = ", k, "\n");
	#cat("n = num samples  = ", n, "\n");

	#cat("inter: BGSS/(k-1):", BGSS/(k-1), "\n");
	#cat("intra: WGSS/(n-k):", WGSS/(n-k), "\n");
	#cat("F-stat: ", fstat, "\n");

	#cat("-----------------------------------------------------------------\n");

	return(fstat);

}

###############################################################################

random_clusters=function(cluster_sizes){

	# Get number clusters
	num_samples=sum(cluster_sizes);
	num_cluster_sizes=length(cluster_sizes);
	#cat("Num cluster sizes: ", num_cluster_sizes, "\n");	

	# Randomly generate clusters with specified sizes
	s=sample(num_samples, replace=F);
	cluster_assignments=rep(0,num_samples);
	offset=1;
	for(i in 1:num_cluster_sizes){
		for(fill in 1:cluster_sizes[i]){
			cluster_assignments[s[offset]]=i;
			offset=offset+1;
		}
	}
	#print(cluster_assignments);

	return(cluster_assignments);
}

###############################################################################

bootstrap_cluster=function(cluster_sizes, dist, nbootstrap, test_stat){

	# Make sure sample count matches distance matrix
	num_samples=sum(cluster_sizes);
	num_samples_in_dist=(1+sqrt(1+length(dist)*8))/2;
	if(num_samples != num_samples_in_dist){
		cat("Error: Num samples based on cluster sizes not the same as those in dist matrix.\n");
	}

	# Run bootstrap
	pfstats=rep(0,nbootstrap);
	for(bs_n in 1:nbootstrap){

		# Generate random clusters with same size, different members
		rnd_clst=random_clusters(cluster_sizes);
		
		# Compute pseudo F statistic for each random cluster
		pfstats[bs_n]=computePseudoFstat(rnd_clst, dist);
	}

	# Compute single tailed p-value
	sorted_pfstats=sort(pfstats);
	num_gt_teststat=sum(test_stat<=pfstats);	
	pvalue=num_gt_teststat/nbootstrap;

	#print(pfstats);

	# package results for returning in list
	result=list();
	result$null_dist=pfstats;
	result$pvalue=pvalue;
	result$mean=mean(pfstats);
	result$lb=min(pfstats);
	result$ub=max(pfstats);

	return(result);
}

#######################################################################################

#d=dist(matrix(runif(39),ncol=3));
#bootstrap_cluster(c(5,3,4,1), d, 40,12);
#q(status=0);

#######################################################################################

MAX_NUM_CLUSTERS=20;
#MAX_NUM_CLUSTERS=NumSamples;

MAX_NUM_CLUSTERS=ifelse(MAX_NUM_CLUSTERS>=NumSamples, NumSamples-1, MAX_NUM_CLUSTERS); 
#MAX_NUM_CLUSTERS=NumSamples-1;


cat("Number of clustering techniques (distances) to try: ", num_distances, "\n");

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
	cluster_list=list();
	null_dist_list=list();
	pseudoFstat_vector=rep(0, MAX_NUM_CLUSTERS-1);
	fstat_by_numclusters=rep(0,MAX_NUM_CLUSTERS);
	pvalue_vector=rep(0, MAX_NUM_CLUSTERS-1);
	cat("Running hclust/cutree from k = 2 to ", MAX_NUM_CLUSTERS, "\n", sep="");

	hcluster=hclust(samp_dist, method="ward");
	for(clstr_idx in 1:(MAX_NUM_CLUSTERS-1)){
		
		num_clusters=clstr_idx+1;
		cat("  Clustering with k = ", num_clusters, "\n");

		# Perform clustering
		cluster_asgnts=cutree(hcluster, k=num_clusters);
		pfstat=computePseudoFstat(cluster_asgnts, samp_dist);	

		# Keep track of F-stat and cluster assignments for each cutoff
		pseudoFstat_vector[clstr_idx]=pfstat;
		fstat_by_numclusters[num_clusters]=pfstat;
		cluster_list[[num_clusters]]=cluster_asgnts;

		if(BootstrapOn){
			# cluster_sizes, dist, nbootstrap, test_stat
			bs_result=bootstrap_cluster(as.vector(table(cluster_asgnts)), samp_dist, 200, pseudoFstat_vector[clstr_idx]);
			null_dist_list[[num_clusters]]=bs_result$null_dist;
			pvalue_vector[num_clusters]=bs_result$pvalue;
		}
	}

	# Compute clusters at different heights
	height_cutoffs=seq(max_dist, 0, -(max_dist/1000));
	num_heights_cutoffs=length(height_cutoffs);
	num_clusters_by_height=rep(0, num_heights_cutoffs);
	fstat_by_height=rep(0, num_heights_cutoffs);
	i=1;
	for(heights in height_cutoffs){
		cluster_asgnts=cutree(hcluster, h=heights);	
		#print(cluster_asgnts);
		num_clusters_by_height[i]=length(unique(cluster_asgnts));
		fstat_by_height[i]=fstat_by_numclusters[num_clusters_by_height[i]];
		i=i+1;
	}
	#print(height_cutoffs);
	#print(num_clusters_by_height);

	# Compute change in pfstat over num_clusters
	delta_by_numclusters=rep(0,MAX_NUM_CLUSTERS-2);
	for(i in 2:(MAX_NUM_CLUSTERS-2)){
		delta_by_numclusters[i]=abs(pseudoFstat_vector[i-1]-pseudoFstat_vector[i])/pseudoFstat_vector[i-1];
	}

	delta_by_height=rep(0, num_heights_cutoffs);
	for(i in 2:num_heights_cutoffs){
		delta_by_height[i]=abs(fstat_by_height[i-1]-fstat_by_height[i])/fstat_by_height[i-1];
	}


	#print(delta);
	#print(pseudoFstat_vector);

	par(mfrow=c(2,1));
	par(oma=c(0,0,0,0));
	par(mar=c(5,4,4,2));

	#cat("Pseudo F-stat vector:\n");
	#print(pseudoFstat_vector);

	pseudoFstat_vector[pseudoFstat_vector==Inf]=0;
	log10pseudoFstat_vector=log10(pseudoFstat_vector);
	log10pseudoFstat_vector[!is.finite(log10pseudoFstat_vector)]=0;

	delta_by_numclusters[!is.finite(delta_by_numclusters)]=0;
	ymax=max(log10pseudoFstat_vector);
	
	#------------------------------------------------------------------------------

	cat("Plotting PseudoFstat...\n");
	plot(2:MAX_NUM_CLUSTERS, log10pseudoFstat_vector, type="b", ylab="log(pseudo F Stat)", xlab="Num clusters", main=dist_name[dist_idx],
		ylim=c(0,ymax));
	axis(side=1, at=2:(MAX_NUM_CLUSTERS-1));

	cat("Plotting Delta...\n");
	plot(.5+(2:(MAX_NUM_CLUSTERS-1)), delta_by_numclusters[2:(MAX_NUM_CLUSTERS-1)],
		 col="blue", type="b", ylab="Rate of C&H F-stat Change", xlab="Num clusters");
	axis(side=1, at=2:(MAX_NUM_CLUSTERS-1));

	#------------------------------------------------------------------------------

	cat("Plotting Height/NumClusters...\n");
	plot(height_cutoffs, num_clusters_by_height, ylab="Num Clusters", xlab="Height Cutoff", type="b", col="green", cex=.5, main=dist_name[dist_idx]);

	cat("Plotting Delta Fstat for Heights...\n");
	y_maxfinite=max(delta_by_height[is.finite(delta_by_height)]);
	plot(height_cutoffs, delta_by_height,  col="blue", type="b", ylab="Rate of C&H F-stat Change", xlab="Height Cutoff", cex=.5,
		ylim=c(0, y_maxfinite*1.1));
	
	nonzero_deltas_idx=(delta_by_height>0 & is.finite(delta_by_height));
	print(num_clusters_by_height[nonzero_deltas_idx]);
	labels=sprintf("%i", num_clusters_by_height[nonzero_deltas_idx]);
	text(height_cutoffs[nonzero_deltas_idx], delta_by_height[nonzero_deltas_idx], labels, pos=3, cex=.8);

	#------------------------------------------------------------------------------

	mds=isoMDS(samp_dist);

	# Compute display ranges/sizes
	min=min(mds$points[,1])
	max=max(mds$points[,1])
	width=(max-min)
	margin=width*(0.1)
	label_scale=48/NumSamples;
	label_scale=ifelse(label_scale>1.2, 1.2, label_scale);

	
	par(mfrow=c(1,1));
	par(oma=c(1,1,2,1));
	par(mar=c(1,1,3.5,1));

	# Plot sample name/text
	recommended_cluster_id=which.max(delta_by_numclusters[1:(NumSamples/2)])+1;
	clusters=cluster_list[[recommended_cluster_id]];

	# output cluster members
	if(ExportClusters){
		cluster_out_fn=paste(OutputFileNameRoot, ".recommended_clusters.", gsub(" ", "_", dist_name[dist_idx]), ".tsv", sep="");
		cluster_out_fh=file(cluster_out_fn, "wt");
		matout=cbind(clusters, names(clusters));
		matout=matout[order(matout[,1]),];
		for(sample in 1:NumSamples){
			cat(file=cluster_out_fh, matout[sample,], "\n");
		}
	}

	# output F statistic across k
	if(ExportFstats){
		fstat_out_fn=paste(OutputFileNameRoot, ".fstats.", gsub(" ", "_", dist_name[dist_idx]), ".tsv", sep="");
		fstat_out_fn=file(fstat_out_fn, "wt");
		cat(file=fstat_out_fn, OutputFileNameRoot, pseudoFstat_vector, sep=",");
		cat(file=fstat_out_fn, "\n");
	}

	# output delta F-statistics across k
	if(ExportDeltaFstats){
		fstat_out_fn=paste(OutputFileNameRoot, ".delta_fstats.", gsub(" ", "_", dist_name[dist_idx]), ".tsv", sep="");
		fstat_out_fn=file(fstat_out_fn, "wt");
		cat(file=fstat_out_fn, OutputFileNameRoot, delta_by_numclusters[1:(NumSamples/2)], sep=",");
		cat(file=fstat_out_fn, "\n");
	}

	# Plot recommended clusters
	plot(mds$points,type="n", xlim=c(min-margin, max+margin), xlab="", ylab="", main=dist_name[dist_idx])
	text(mds$points,labels=sample_names, cex=label_scale, col=clusters)
	mtext(sprintf("k = %i", recommended_cluster_id));

	if(BootstrapOn){
		par(mfrow=c(3,3));
	}else{
		par(mfrow=c(3,2));
	}
	par(oma=c(1,1,2,1));
	par(mar=c(1,1,3.5,1));

	# Plot isoMDS for each cluster assignment
	for(num_clusters in 2:MAX_NUM_CLUSTERS){
		clusters=cluster_list[[num_clusters]];

		# Plot sample name/text
		#cat("Plotting names...\n");
		plot(mds$points,type="n", xlim=c(min-margin, max+margin), xlab="", ylab="", main=dist_name[dist_idx])
		text(mds$points,labels=sample_names, cex=label_scale, col=clusters)
		mtext(sprintf("k = %i, pseudo F-stat: %3.4f, p-value: %3.4f", 
			num_clusters, pseudoFstat_vector[num_clusters-1], pvalue_vector[num_clusters]),
			cex=.8);

		# Plot points
		#cat("Plotting points...\n");
		plot(mds$points,type="p", xlim=c(min-margin, max+margin), xlab="", ylab="", col=clusters, pch=19)

		if(BootstrapOn){
			# Plot null distribution	
			#cat("Plotting null distribution...\n");
			hist(null_dist_list[[num_clusters]], xlim=c(0, pseudoFstat_vector[num_clusters-1]*1.25), main="Null Distribution and Pseudo F-Stat");
			abline(v=pseudoFstat_vector[num_clusters-1], col="blue", lty="dashed");
		}
	}

} # End of loop for each distance measure to try clustering with

par(mfrow=c(num_distances,1));

# Plot dendrograms
for(dist_idx in 1:num_distances){
	samp_dist=dist_list[[dist_idx]];
	#print(samp_dist);
	plot(hclust(samp_dist, method="ward"), cex=label_scale, main=dist_name[dist_idx]);
}


wards=function(d){
	return(hclust(d,method="ward"));
}

#-------------------------------------------------------------------------------

# Reorder categories/taxonomies by mean abundance
mean_prob=apply(prob_mat, 2, mean);
ordering=order(mean_prob, decreasing=TRUE);
sorted_prob_mat=prob_mat[,ordering];
num_cat=ncol(sorted_prob_mat);

# Prep top ten matrix
TOP_MAX=20;
top_prob_mat=sorted_prob_mat[,1:TOP_MAX];

color_scheme=rev(grey(seq(0,1, length.out=10)));

#-------------------------------------------------------------------------------
# Draw heatmaps for weighted rank difference

wrd_p4=function(x){
	weight_rank_dist(x,4);
}

heatmap(sorted_prob_mat, cexRow=label_scale,  cexCol=label_scale * 0.70,
	Colv=NA,
	distfun=wrd_p4,
	hclustfun=wards,
	col=color_scheme, 
	scale="none",
	margins=c(7,3),
	main="Weighted Rank Difference: All"
)

heatmap(top_prob_mat, cexRow=label_scale,  cexCol=label_scale * 0.90,
	Colv=NA,
	distfun=wrd_p4,
	hclustfun=wards,
	col=color_scheme, 
	scale="none",
	margins=c(7,3),
	main="Weighted Rank Difference: Top 10"
)

#-------------------------------------------------------------------------------
# Draw heatmaps for euclidean

heatmap(sorted_prob_mat, cexRow=label_scale,  cexCol=label_scale * 0.70,
	Colv=NA,
	distfun=dist,
	hclustfun=wards,
	col=color_scheme, 
	scale="none",
	margins=c(7,3),
	main="Euclidean: All"
)

heatmap(top_prob_mat, cexRow=label_scale,  cexCol=label_scale * 0.70,
	Colv=NA,
	distfun=dist,
	hclustfun=wards,
	col=color_scheme, 
	scale="none",
	margins=c(7,3),
	main="Euclidean: Top 10"
)

#-------------------------------------------------------------------------------



cat("Done.\n");
q(status=0);
