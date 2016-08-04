#!/usr/local/bin/Rscript

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
	"Run principal component analysis on summary table.",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileNameRoot=opt$input_file;

###############################################################################
# Main program loop

cat("Working on ", InputFileName, "\n", sep="");

# Load summary_table.xls
mat<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, row.names=1))

# Exclude total counts column
count_mat=mat[,2:ncol(mat)]
num_orig_cat=ncol(count_mat);
num_orig_samp=nrow(count_mat);
cat("Orig num categories: ", num_orig_cat, "\n");
cat("Orig num samples: ", num_orig_samp, "\n");

# Identify 0 count columns so we can filter them
col_sum=apply(count_mat, 2, sum);
col_idx_iszero=which(col_sum==0);
cat("Zero count categories index:", col_idx_iszero, "\n");

# Remove zero count columns
if(length(col_idx_iszero)>0){
	cat("Removing zero count categories...\n");
	count_mat=count_mat[,-col_idx_iszero];
}

# Get column/category names
cat_names=as.vector(colnames(count_mat));
num_used_categories=ncol(count_mat);

# Compute shortened names
short_names=character(num_used_categories);
for(i in 1:num_used_categories){
	taxonomy=unlist(strsplit(cat_names[i], " "));
	short_names[i]=taxonomy[length(taxonomy)];
}

# Get sample names
sample_names=rownames(count_mat);

##################################################################################

num_used_samples=num_orig_samp;

cat("\n");
cat("Num Used Categories: ", num_used_categories, "\n");
cat("Num Used Samples: ", num_used_samples, "\n");
cat("\n");

##################################################################################

# Normalize
prob_mat=matrix(nrow=num_used_samples, ncol=num_used_categories);
colnames(prob_mat)=short_names;
rownames(prob_mat)=sample_names;
for(i in 1:num_used_samples){
    	sample_totals=sum(count_mat[i,]);
	prob_mat[i,]=count_mat[i,]/sample_totals;
}
print(prob_mat);

mean_prob_vect=apply(prob_mat, 2, mean);
sorted_mean_prob_vect=sort(mean_prob_vect, decreasing=T, index.return=T);
print(sorted_mean_prob_vect);

# reorder matrix
prob_mat=prob_mat[,sorted_mean_prob_vect$ix];
short_names=short_names[sorted_mean_prob_vect$ix];

print(colnames(prob_mat));
PERCENT_CUTOFF=1;
cutoff=max(which(cumsum(sorted_mean_prob_vect$x)<=PERCENT_CUTOFF));
if(cutoff==-Inf){
	cutoff==1;
}
print(prob_mat);


# Compute log ratios

num_ratios=choose(num_used_categories,2);
logratio_matrix=matrix(1, nrow=num_used_samples, ncol=num_ratios);
ratio_names=rep("",num_ratios);
k=1;

num_used_categories=cutoff;
cat("Using ", cutoff, " categories, because they represent <", PERCENT_CUTOFF*100, "% of the data.\n", sep="");

for(i in 1:(num_used_categories-1)){
	for(j in (i+1):num_used_categories){

		logratio=log(prob_mat[,i]/prob_mat[,j]);
		finite_logratios=logratio[is.finite(logratio)];
		
		if(length(finite_logratios)){
			min=-1;
			max=1;
		}

		min=min(finite_logratios);
		max=max(finite_logratios);

		logratio=ifelse(is.nan(logratio), 0, logratio);
		logratio=ifelse(logratio==-Inf, min*2, logratio);
		logratio=ifelse(logratio==Inf, max*2, logratio);

		# If there are still -Inf it's because the max was 0
		logratio=ifelse(logratio==-Inf, -1, logratio);

		if(sum(logratio)!=0){
			ratio_name=paste(short_names[i], "/", short_names[j], sep="");
			ratio_matrix=as.matrix(logratio, ncol=1);
			ratio_names[k]=ratio_name;
			logratio_matrix[,k]=ratio_matrix;
			k=k+1;		
		}
	}
}

colnames(logratio_matrix)=ratio_names;
rownames(logratio_matrix)=sample_names;
print(ratio_names);

# Get rid of all the extra columns
logratio_matrix=logratio_matrix[,1:(k-1)];
#print(logratio_matrix);


###############################################################################

pdf(paste(OutputFileNameRoot, ".pdf", sep=""), height=12, width=12);


pca_result=prcomp(prob_mat, scale=TRUE);
bi=biplot(pca_result, cex=1);
plot(pca_result$x[,1], pca_result$x[,2], "p");
#text(pca_result$x[,1], pca_result$x[,2], sample_names, cex=1);
plot(pca_result$x[,1], pca_result$x[,2], "n");
text(pca_result$x[,1], pca_result$x[,2], sample_names, cex=1);

pca_result=prcomp(logratio_matrix, scale=TRUE);
bi=biplot(pca_result, cex=1);
plot(pca_result$x[,1], pca_result$x[,2], "p");
#text(pca_result$x[,1], pca_result$x[,2], sample_names, cex=1);
plot(pca_result$x[,1], pca_result$x[,2], "n");
text(pca_result$x[,1], pca_result$x[,2], sample_names, cex=1);





###############################################################################

cat("Done.\n");
warnings();
q(status=0);
