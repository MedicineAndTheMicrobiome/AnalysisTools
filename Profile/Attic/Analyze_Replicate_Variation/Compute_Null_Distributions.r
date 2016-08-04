#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
        "input_file", "i", 1, "character",
        "sample_size", "n", 1, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage=paste (
		"\nUsage:\n\t", script_name, "\n",
		"\t\t-i <Input summary_table.xls FileName>\n",
		"\t\t-n <Sample Size>\n",
		"\n",
		"This code will read in a summary table xls file and for each sample compute\n",
		"through bootstrapping the null distribution for an indices.\n",
		"\n");

if(!length(opt$input_file) || !length(opt$sample_size)){
	cat(usage);
	quit(status=-1);
}

InputFileName=opt$input_file;
SampleSize=opt$sample_size;

###############################################################################
# Main program loop

NDPlot= paste(gsub("\\.xls$","",InputFileName), ".null_distributions.pdf", sep="")

cat("\n")
cat("                 Input File Name: ", InputFileName, "\n");
cat("                Output File Name: ", NDPlot, "\n");

###############################################################################
###############################################################################

library(vegan);

# Load data
A=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, row.names=1, check.names=FALSE));
#cat("Original Matrix:\n")
#print(A);

# Extract out only the counts (ignore totals)
count_mat=A[,2:ncol(A)];
num_categories=ncol(count_mat);
num_samples=nrow(count_mat);

cat("Num Samples: ", num_samples, "\n");
cat("Num Taxa: ", num_categories, "\n");

sample_names=rownames(count_mat);
category_names=colnames(count_mat);

###############################################################################

# Compute normalized
total=apply(count_mat, 1, sum);
#print(total);
norm_mat=count_mat/total;
#print(norm_mat);

###############################################################################

order_dist=function(a, b){
        arank=rank(a, ties.method="average");
        brank=rank(b, ties.method="average");

        sort_sqrdiff=sqrt(sum(((arank-brank)^2)*(a^2+b^2)));
        #sort_sqrdiff=sqrt(sum(((arank-brank)^2)*(((a-b)/2)^2)));
        return(sort_sqrdiff);
}

weight_rank_dist=function(M){
        NumSamples=nrow(M);
        order_dist_mat=matrix(0, nrow=NumSamples, ncol=NumSamples);
        for(i in 1:NumSamples){
                for(j in 1:i){
                        order_dist_mat[i,j]=order_dist(M[i,], M[j,]);
                }
        }
        rownames(order_dist_mat)=rownames(M);
        return(as.dist(order_dist_mat));
}

###############################################################################

NUMBOOTSTRAPS=400;

dist_names=character();
dist_names[1]="Weighted Rank Difference";
dist_names[2]="Euclidean";
num_distances=length(dist_names);

pdf(NDPlot, width=11,height=8.5)
laymat=matrix(c(1,2,3,3), nrow=2, ncol=2, byrow=TRUE);
layout(laymat);

for(i in 1:num_samples){

	cat("Working on: ", sample_names[i], "\n");

	sorted_prob=sort(norm_mat[i,], decreasing=TRUE);
	
	# Generate bootstraps
	bootstraps=matrix(0, nrow=NUMBOOTSTRAPS, ncol=num_categories);
	for(b in 1:NUMBOOTSTRAPS){
		samples=sample(1:num_categories, size=SampleSize, replace=TRUE, prob=norm_mat[i,]);
		bootstraps[b,]=as.vector(table(c(c(1:num_categories),samples)))-1;
	}

	#print(bootstraps);
	compare_input=matrix(0, nrow=2, ncol=num_categories);
	compare_input[1,]=norm_mat[i,];
	distances=matrix(0, nrow=NUMBOOTSTRAPS, ncol=num_distances);

	for(b in 1:NUMBOOTSTRAPS){
		compare_input[2,]=bootstraps[b,]/sum(bootstraps[b,]);

		distances[b,1]=weight_rank_dist(compare_input)[1];
		distances[b,2]=dist(compare_input)[1];
	}

	#print(distances);
	for(didx in 1:num_distances){
		hist(distances[,didx], main=dist_names[didx]);
		coef_of_var=sd(distances[,didx])/mean(distances[,didx]);
		mtext(sprintf("Coef of Var: %3.2f", coef_of_var), side=3, line=-1);
	}

	non_zero_idx=sorted_prob>0;
	barplot(sorted_prob[non_zero_idx], las=2);
	
}





dev.off();

###############################################################################

cat("Done.\n")

warnings();
q(status=0)
