
#------------------------------------------------------------------------------

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

#------------------------------------------------------------------------------

sample_counts_to_diversity_matrix=function(counts){
	
	cat("Calculating Diversity Indices...\n");
	normalized=normalize(counts);
	num_samples=nrow(counts);
	num_categories=ncol(counts);
	
	cat("Number of Samples: ", num_samples, "\n");
	cat("Number of Categories: ", num_categories, "\n");

	div_names=c("Tail", "Shannon", "Simpson", "Evenness", "SimpsonsRecip");
	num_div_idx=length(div_names);

	div_mat=matrix(0, nrow=num_samples, ncol=num_div_idx);
	colnames(div_mat)=div_names;
	rownames(div_mat)=rownames(normalized);

	for(i in 1:num_samples){
		curNorm=normalized[i,];
		zeroFreeNorm=curNorm[curNorm>0]

		div_mat[i,"Tail"]=tail_statistic(zeroFreeNorm);
		div_mat[i,"Shannon"]=-sum(zeroFreeNorm*log(zeroFreeNorm));
		div_mat[i,"Simpson"]=1-sum(curNorm^2);
		div_mat[i,"Evenness"]=div_mat[i,"Shannon"]/log(length(zeroFreeNorm));
		div_mat[i,"SimpsonsRecip"]=1/sum(curNorm^2);
	}

	# Set evenness to 0 if there is only 1 category.
	evenness=div_mat[,"Evenness"];
	div_mat[is.na(evenness),"Evenness"]=0;

	return(div_mat);
}

#------------------------------------------------------------------------------

