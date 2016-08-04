
sample_data=cbind(
	c( 0, 5,11,11,14,14),
	c( 5, 0,10, 6,13,15),
	c(11,10, 0, 6,17,21),
	c(11, 6, 6, 0,13,15),
	c(14,13,17,13, 0, 6),
	c(14,15,21,15, 6, 0)
	);
sample_data=sqrt(sample_data);
rownames(sample_data)=c("A","B","C","D","E","F");
colnames(sample_data)=c("A","B","C","D","E","F");
print(sample_data);
		

###############################################################################

# Computes centroid based on entire distance matrix and a group.
sum_sqr=function(dist_mat){
	if(nrow(dist_mat)==1){
		return(0);
	}else{
		n=ncol(dist_mat);
		ss=sum(as.dist(dist_mat)^2);
		mean=ss/n;
		return(mean);
	}
}

###############################################################################

# Computes the fstat for a set of clusters and a distance matrix
compute_fstat=function(dist_mat, clusters){
	num_samples=nrow(dist_mat);
	num_clusters=max(clusters);

	cat("Num Samples: ", num_samples, "\n");
	cat("Num Clusters: ", num_clusters, "\n");
	print(clusters);

	#-----------------------------------------------------------------------------
	# Compute SSW
	gss=numeric(num_clusters);
	for(cl in 1:num_clusters){
		#cat("cluster: ", cl, "\n");
		grp_idx=which(clusters==cl);
		group_dist_mat=as.matrix(dist_mat[grp_idx, grp_idx], nrow=length(grp_idx));
		gss[cl]=sum_sqr(group_dist_mat);
	}
	wgss=sum(gss);
	cat("WGSS: ", wgss, "\n");

	#-----------------------------------------------------------------------------
	# Compute SST 
	tss=sum_sqr(dist_mat);
	cat("TSS: ", tss, "\n");

	#-----------------------------------------------------------------------------
	# Compute F-stat
	bgss=tss-wgss;
	k=num_clusters;
	n=num_samples;
	fstat=(bgss/(k-1))/(wgss/(n-k));
	print(fstat);

	#-----------------------------------------------------------------------------
	# Save/Output F-stat and components
	return(fstat);
}

###############################################################################

#sample_data=as.matrix(read.table("/usr/local/projects/NTD/SVN/DAS/ANDES/peptides_wReference.r_distmat", header=TRUE, check.names=FALSE));
#print(sample_data);

num_samples=nrow(sample_data);
example_tree=hclust(as.dist(sample_data),"ward");

num_fstats=num_samples;
fstats=rep(0, num_fstats);

for(k in 2:num_fstats){
	cat("**********************************************************************\n");
	cat("k = ", k, "\n");
	clusters=cutree(example_tree, k);
	#print(clusters);
	result=compute_fstat(sample_data, clusters);
}









