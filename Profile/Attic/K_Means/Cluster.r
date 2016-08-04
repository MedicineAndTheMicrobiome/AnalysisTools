#!/usr/local/packages/R-2.11.1/bin/Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n\t\t<Input summary_table.xls FileName>\n\n",
	"	-i <input summary_table.xls file>\n",
	"\n",
	"Run K-means on summary table.\n",
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

order_dist=function(a, b){
	asort=sort(a, decreasing=T, index.return=T, method="shell");
	bsort=sort(b, decreasing=T, index.return=T, method="shell");

	sort_sqrdiff=sqrt(sum(((asort$ix-bsort$ix)^2)*(asort$x^2+bsort$x^2)));
	#sort_sqrdiff=sqrt(sum(((asort$ix-bsort$ix)^2)*(asort$x*bsort$x)));
	return(sort_sqrdiff);
	
}

###############################################################################

library(MASS);
library(stats);
library(vegan);

###############################################################################

order_dist_mat=matrix(0, nrow=NumSamples, ncol=NumSamples);
for(i in 1:NumSamples){
	for(j in 1:i){
		order_dist_mat[i,j]=order_dist(prob_mat[i,], prob_mat[j,]);
	}
}
rownames(order_dist_mat)=rownames(prob_mat);

pdf(paste(OutputFileNameRoot, ".clustered.pdf", sep=""), height=8.5, width=11);

dist_list=list();

countmat=prob_mat*10000;
intmat=apply(countmat,2,as.integer);
rownames(intmat)=rownames(prob_mat);
#print(intmat);

dist_list[[1]]=as.dist(order_dist_mat);
dist_list[[2]]=dist(prob_mat);
dist_list[[3]]=vegdist(intmat, method="morisita");

dist_name=vector();
dist_name[1]="Weighted Rank Difference";
dist_name[2]="Euclidean";
dist_name[3]="Morisita-Horn";

num_distances=2;
#num_distances=length(dist_list);


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

MAX_NUM_CLUSTERS=NumSamples;
#MAX_NUM_CLUSTERS=NumSamples;

MAX_NUM_CLUSTERS=ifelse(MAX_NUM_CLUSTERS>=NumSamples, NumSamples-1, MAX_NUM_CLUSTERS); 
#MAX_NUM_CLUSTERS=NumSamples-1;


cat("Number of clustering techniques (distances) to try: ", num_distances, "\n");

for(dist_idx in 1:num_distances){

	cat("Working on: ", dist_name[dist_idx], "\n", sep="");
	samp_dist=dist_list[[dist_idx]];


	# Remove 0 distances so isoMDS doesn't freakout
	for(i in 1:length(samp_dist)){
		if(samp_dist[i]==0){
			samp_dist[i]=1e-323;
		}
	}
	#cat("Removed 0 distances.\n");

	# Compute cluster clusters
	cluster_list=list();
	pseudoFstat_vector=rep(0, MAX_NUM_CLUSTERS-1);
	pvalue_vector=rep(0, MAX_NUM_CLUSTERS-1);
	cat("Running hclust/cutree from k = 2 to ", MAX_NUM_CLUSTERS, "\n", sep="");
	for(clstr_idx in 1:(MAX_NUM_CLUSTERS-1)){
		
		num_clusters=clstr_idx+1;
		cat("  Clustering with k = ", num_clusters, "\n");

		# Run k-means multiple times and take the best result
		numTrials=1; # only need to run this once if we're using ward's method
		trial_clusters=matrix(0,nrow=numTrials, ncol=NumSamples);
		trial_pfstat=rep(0,numTrials);
		for(trials in 1:numTrials){
			#km_res=kmeans(as.matrix(samp_dist), num_clusters, algorithm="MacQueen");
			#cluster_asgnts=as.vector(km_res$cluster);

			hcluster=hclust(samp_dist, method="ward");
			cluster_asgnts=cutree(hcluster, k=num_clusters);

			trial_clusters[trials,]=cluster_asgnts;
			trial_pfstat[trials]=computePseudoFstat(cluster_asgnts, samp_dist);	
		}

		#print(trial_pfstat);
		cat("    mean pseudo-fstat: ", mean(trial_pfstat), ", Range=(", min(trial_pfstat), " - ",  max(trial_pfstat), ") \n");

		# Pick the best trials based on highest pseudo-F-stat
		best_cluster=which.max(trial_pfstat);
		pseudoFstat_vector[clstr_idx]=trial_pfstat[best_cluster];
		cluster_list[[num_clusters]]=trial_clusters[best_cluster,];

		# cluster_sizes, dist, nbootstrap, test_stat
		bs_result=bootstrap_cluster(as.vector(table(trial_clusters[best_cluster,])), samp_dist, 100, pseudoFstat_vector[clstr_idx]);
		pvalue_vector[num_clusters]=bs_result$pvalue;
	}

	# Compute change in pfstat
	delta=rep(0,MAX_NUM_CLUSTERS-1);
	for(i in 2:(MAX_NUM_CLUSTERS-2)){
		delta[i]=abs(pseudoFstat_vector[i-1]-pseudoFstat_vector[i])/pseudoFstat_vector[i-1];
	}

par(mfrow=c(1,1));
par(oma=c(0,0,0,0));
par(mar=c(5,4,4,2));
	#cat("Pseudo F-stat vector:\n");
	#print(pseudoFstat_vector);
	plot(2:MAX_NUM_CLUSTERS, log10(pseudoFstat_vector), type="b", ylab="log(pseudo F Stat)", xlab="Num clusters", main=dist_name[dist_idx],
		ylim=c(0,max(log10(pseudoFstat_vector))));
	lines(.5+(2:(MAX_NUM_CLUSTERS-1)), delta[2:(MAX_NUM_CLUSTERS-1)], col="blue", type="b");
	axis(side=1, at=2:(MAX_NUM_CLUSTERS-1));
	mtext("Rate", col="blue", side=2, at=1, line=3);

	cat("Delta:\n");
	print(delta);

par(mfrow=c(3,2));
par(oma=c(1,1,2,1));
par(mar=c(1,1,3.5,1));

	mds=isoMDS(samp_dist);

	# Compute display ranges/sizes
	min=min(mds$points[,1])
	max=max(mds$points[,1])
	width=(max-min)
	margin=width*(0.1)
	label_scale=42/NumSamples;
	label_scale=ifelse(label_scale>1.2, 1.2, label_scale);

	for(num_clusters in 2:MAX_NUM_CLUSTERS){
		clusters=cluster_list[[num_clusters]];

		# Plot texts
		plot(mds$points,type="n", xlim=c(min-margin, max+margin), xlab="", ylab="", main=dist_name[dist_idx])
		text(mds$points,labels=sample_names, cex=label_scale, col=clusters)
		mtext(sprintf("k = %i, pseudo F-stat: %3.4f, p-value: %3.4f", num_clusters, pseudoFstat_vector[num_clusters-1], pvalue_vector[num_clusters]));

		# Plot points
		plot(mds$points,type="p", xlim=c(min-margin, max+margin), xlab="", ylab="", col=clusters, pch=19)
	}
}

par(mfrow=c(num_distances,1));

# Plot dendrograms
for(dist_idx in 1:num_distances){
	samp_dist=dist_list[[dist_idx]];
	#print(samp_dist);
	plot(hclust(samp_dist, method="ward"), cex=label_scale, main=dist_name[dist_idx]);
}

cat("Done.\n");
q(status=0);
