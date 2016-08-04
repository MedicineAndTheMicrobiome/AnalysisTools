#!/usr/local/bin/Rscript

###############################################################################

library('getopt');

NumBootstraps=1000;

params=c(
	"input_file", "i", 1, "character",
	"list_of_pairs", "l", 1, "character",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"	-i <input summary_table.xls file>\n",
	"	-l <list of taxa pairs to correlate and identify p-values for.>\n",
	"	[-o <output filename root>]\n",
	"\n",
	"Reads in the summary table and then permutes abundances, and then computes\n",
	"correlation coefficients, in order to generate a null distribution.\n",
	"Then based on the pairs of taxas compute the pearson correlation and p-values.\n",
	"\n",
	"\n",
	"\n", sep="");

if(!length(opt$input_file) || !length(opt$list_of_pairs)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputRoot=gsub(".summary_table.xls", "", InputFileName);
ListOfPairs=opt$list_of_pairs;

if(length(opt$output_file)){
	OutputRoot=opt$output_file;
}

cat("Input File: ", InputFileName, "\n");
cat("Output File Root: ", OutputRoot, "\n");
cat("List of Pairs: ", ListOfPairs, "\n");

################################################################################

load_summary_table=function(filename){
	st=as.matrix(read.table(filename, header=TRUE, sep="\t", row.names=1, check.names=F));
	return(st[,2:ncol(st)]);
}

load_list_pairs=function(filename){
	pairs=as.matrix(read.table(filename, header=F, sep="\t", check.names=F));
	ncol=ncol(pairs);
	cat("Num columns in list: ", ncol, "\n");
	return(pairs);
}

#-------------------------------------------------------------------------------

normalize=function(st){
	sums=apply(st, 1, sum);
	n=matrix(0,nrow=nrow(st), ncol=ncol(st));
	rownames(n)=rownames(st);
	colnames(n)=colnames(st);
	for(i in 1:nrow(st)){
		n[i,]=st[i,]/sums[i];
	}
	return(n);
}

#-------------------------------------------------------------------------------

permute=function(nst){
	nrow=nrow(nst);
	ncol=ncol(nst);
	nabund=nrow*ncol;
	rand_ix=order(runif(nabund));	
	abund_vect=as.vector(nst);
	permuted_vect=abund_vect[rand_ix];
	permuted_mat=matrix(permuted_vect, nrow=nrow, ncol=ncol);
	norm_perm_mat=normalize(permuted_mat);
	return(norm_perm_mat);
}

find_two_tail_pvalue=function(test_stat, sorted_null){
	num_null=length(sorted_null);
	num_greater=sum(sorted_null<test_stat);
	perc_greater=num_greater/num_null;
	if(perc_greater > .5){
		return(2*(1-perc_greater));
	}else{
		return(2*perc_greater);
	}
}


################################################################################

# Load list of pairs
pairs=load_list_pairs(ListOfPairs);

if(ncol(pairs)==1){
	n=nrow(pairs);
	syn_pairs=matrix("", nrow=(n*(n-1)/2),  ncol=2);
	idx=1;
	for(i in 1:n){
		for(j in 1:n){
			if(i<j){
				syn_pairs[idx,]=c(pairs[i,1], pairs[j,1])
				idx=idx+1;
			}
		}
		
	}
	pairs=syn_pairs;
}
#print(pairs);

# Load summary table
st=load_summary_table(InputFileName);
num_taxa=ncol(st);
num_samples=nrow(st);
#print(st);

# Normalize counts
nst=normalize(st);
#print(nst)

################################################################################

pdf(paste(OutputRoot, ".cor.pdf", sep=""), height=11, width=8.5);

################################################################################
# Compute Null distribution
cat("\n");
cat("Computing Null Distribution for this datset.\n");
null_dist=numeric();
null_abund=numeric();
for(bs in 1:NumBootstraps){
	# Permute
	permuted_mat=permute(nst);
	# Compute correlation
	null_dist[bs]=cor(permuted_mat[,1], permuted_mat[,2]);
	null_abund[bs]=mean(c(permuted_mat[,1], permuted_mat[,2]));
		
	cat(".");
}
sorted_null_dist=sort(null_dist);

cat("\n");
null_median=median(null_dist);
cat("Null Distribution Summary Statistics:\n");
summary(null_median);
hist_res=hist(null_dist, breaks=20, plot=FALSE);

################################################################################
# Compute correlation of randomly selected pairs from this summary table
cat("\nComputing correlation distribution for within samples.\n");
in_table_corr=numeric(NumBootstraps);
in_table_abund=numeric(NumBootstraps);
stop=FALSE;
i=1;
while(i<=NumBootstraps){
	pair=sample(1:num_taxa, 2, replace=FALSE);
	taxa_abund_A=(nst[,pair[1]]);
	taxa_abund_B=(nst[,pair[2]]);
	if(!(all(taxa_abund_A==0) || all(taxa_abund_B==0))){
		in_table_corr[i]=cor(nst[,pair[1]], nst[,pair[2]]);
		in_table_abund[i]=mean(c(nst[,pair[1]], nst[,pair[2]]));
		i=i+1;
		cat(".");
	}
}
cat("\n");
summary(in_table_corr);
intable_hist=hist(in_table_corr, breaks=20, plot=FALSE);
cat("\n\n");

################################################################################
# Plot distribution for within sample and null distribution

par(mfrow=c(2,1));

plot(intable_hist, xlab="Correlation Coefficient", ylab="Frequency", main="Input Correlation Distribution", xlim=c(-1,1));
intable_median_cor=median(in_table_corr);
abline(v=intable_median_cor);
text(intable_median_cor, max(intable_hist$counts), sprintf("Median = %3.4f", intable_median_cor), cex=.8, col="blue", pos=4);
plot(in_table_abund, in_table_corr);

plot(hist_res, main="Null Distribution", xlab="Correlation", xlim=c(-1,1));
null_median=median(null_dist);
abline(v=null_median, col="blue");
text(null_median, max(hist_res$counts), sprintf("Median = %3.4f", null_median), cex=.8, col="blue", pos=4);
plot(null_abund, null_dist);

################################################################################

# Compute correlation coefficient between specified pairs
num_pairs=nrow(pairs);
cat("Num Pairs: ", num_pairs, "\n", sep="");
correl=numeric(num_pairs);
pval=numeric(num_pairs);
for(i in 1:num_pairs){
	a=pairs[i,1];
	b=pairs[i,2];
	cat("Comparing: ", a, " vs. ", b, "\n");
	a_val=nst[,a];
	b_val=nst[,b];
	correl[i]=cor(a_val, b_val);
	cat("Correlation: ", correl[i], "\n");

	pval[i]=find_two_tail_pvalue(correl[i], sorted_null_dist);
	cat("p-value: ", pval[i], "\n");
}

################################################################################

# Output correlation and p-values
fh=file(paste(OutputRoot, ".cor.txt", sep=""), "w");
cat(file=fh, paste("InputFileRoot", "TaxaA", "TaxaB", "CorrelationCoefficient", "p-value", sep=","), "\n");
for(i in 1:num_pairs){
	cat(file=fh, paste(OutputRoot, pairs[i,1], pairs[i,2], correl[i], pval[i], sep=","), "\n");
}
close(fh);

# If not too many pairs to compare, output plots
#if(num_pairs<18){
	par(mfrow=c(3,2));
	for(i in 1:num_pairs){
		title=paste(pairs[i,1], "\n vs. \n", pairs[i,2], sep="");
		plot(hist_res, main=title, xlab="Correlation", xlim=c(-1,1));
		abline(v=correl[i], col="blue");
		text(correl[i], max(hist_res$counts), sprintf(" corr = %3.4f", correl[i]), cex=.8, col="blue", pos=4);
		text(correl[i], max(hist_res$counts), sprintf("\n\np-val = %3.4f", pval[i]), cex=.8, col="blue", pos=4);
	}
	dev.off();
#}

################################################################################

cat("Done.\n")

q(status=0)
