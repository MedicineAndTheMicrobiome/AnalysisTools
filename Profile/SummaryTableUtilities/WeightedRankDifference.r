#!/usr/bin/env Rscript

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

weight_rank_dist_opt=function(M, deg){
	NumSamples=nrow(M);
	order_matrix=matrix(0, nrow=nrow(M), ncol=ncol(M));
	for(i in 1:NumSamples){
		order_matrix[i,]=rank(M[i,], ties.method="average");
	}
	
	dist_mat=matrix(0, nrow=NumSamples, ncol=NumSamples);
	colnames(dist_mat)=rownames(M);
	rownames(dist_mat)=rownames(M);
	for(i in 1:NumSamples){
		for(j in 1:i){
			dist_mat[i,j]=
				sqrt(sum((
					(order_matrix[i,]-order_matrix[j,])^2)*
					(((M[i,]+M[j,])/2)^deg)
					)
				);
		}
	}	
	return(as.dist(dist_mat));

}

###############################################################################

test=T;
if(test){
	test_mat=matrix(c(
		1,2,3,4,5,
		100,200,300,400,550,
		80,90,75,100,225,
		5,3,2,1,4,
		500,300,200,100,400), 
		byrow=T,
		ncol=5);
	
	rownames(test_mat)=c("A","B","C","D","E");

	# Normalize counts
	tot=apply(test_mat, 1, sum);
	norm=test_mat;
	for(i in 1:nrow(test_mat)){
		norm[i,]=test_mat[i,]/tot[i];
	}

	cat("\nTest Count Matrix:\n");
	print(test_mat);
	cat("\nTest Normalized Matrix:\n");
	print(norm);

	m=weight_rank_dist(norm, 4);	
	cat("\nWRD Basic Compute:\n");
	print(m);
	cat("\n");

	m=weight_rank_dist_opt(norm, 4);
	cat("WRD Optimized Compute:\n");
	print(m);
	
	# Remove test variables
	rm(test_mat, tot, norm, i, m, test);
	#ls();
}

cat("\nWeighted Rank Distance Functions Loaded.\n\n");
