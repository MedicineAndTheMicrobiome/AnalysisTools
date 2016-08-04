###############################################################################

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

cat("\nWeighted Rank Distance Functions Loaded.\n\n");
