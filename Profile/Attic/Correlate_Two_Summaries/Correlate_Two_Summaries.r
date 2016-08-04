#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file_r", "r", 1, "character",
	"input_file_c", "c", 1, "character",
	"output_file", "o", 2, "character",
	"pvalue_cutoff", "p", 2, "numeric",
	"bonferroni_correct", "b", 2, "logical",
	"generate_plot", "g", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-r <input summary table.xls Rows>\n",
	"	-c <input summary table.xls Columns>\n",
	"	-o <output filename root>\n",
	"	[-p <p-value cutoff to report, eg. .05, uncorrected, single test.>]\n",
	"	[-b (flag to bonferroni correct for multiple testing)]\n",
	"	[-g <generate plot>]\n",
	"\n",	
	"This script will read in both summary files, normalize them and then compute a correlation\n",
	"between all the categories in A versus B.\n",
	"\n",
	"The output will be organized such that there will be a row for each category in A and\n",
	"there will be a column for the top N correlations with the categories in B.\n",
	"Both files must share sample names.\n",
	"\n",
	"\n");

if(!length(opt$input_file_r) || !length(opt$input_file_c) || !length(opt$output_file)){
	cat(usage);
	q(status=-1);
}

RowsInputFile=opt$input_file_r;
ColumnsInputFile=opt$input_file_c;
OutputFileRoot=opt$output_file;
PvalueCutoff=opt$pvalue_cutoff;
if(length(PvalueCutoff)==0){
	PvalueCutoff=1;
}
BonferroniCorrect=length(opt$bonferroni_correct)>0;
GeneratePlot=length(opt$generate_plot);

###############################################################################

cat("\n")
cat("Rows Input File Name: ", RowsInputFile, "\n");
cat("Columns Output File Name: ", ColumnsInputFile, "\n");       
cat("Output File Name Root: ", OutputFileRoot, "\n");
cat("Generate Plot: ", GeneratePlot, "\n");
cat("\n");
cat("P-value Cutoff: ", PvalueCutoff, ", correct for multiple testing: ", BonferroniCorrect, "\n");
cat("\n");

###############################################################################
###############################################################################

loadSummaryFileTable=function(filename){
	cat("Loading: '", filename, "'\n", sep="");
	inmat=as.matrix(read.table(filename, sep="\t", header=TRUE, check.names=FALSE, row.names=1))
	counts_mat=inmat[,2:(ncol(inmat))];
	#cat("Counts:\n");
	#print(counts_mat);
	cat("ok.\n");
	return(counts_mat);
}

normalizeCounts=function(counts){
	# Sum sample totals
	sample_totals=numeric();
	sample_totals=apply(counts, 1, sum);
	#print(sample_totals);

	# normalize, to compute probabilities
	normalized=matrix(0, nrow=nrow(counts), ncol=ncol(counts));
	for(i in 1:nrow(counts)){
		normalized[i,]=counts[i,]/sample_totals[i];
	}
	#print(normalized);
	return(normalized);
}
###############################################################################

rows_matrix=loadSummaryFileTable(RowsInputFile);
cols_matrix=loadSummaryFileTable(ColumnsInputFile);

normalized_rows_matrix=normalizeCounts(rows_matrix);
normalized_cols_matrix=normalizeCounts(cols_matrix);
#cat("Normalized rows matrix:\n");
#print(normalized_rows_matrix);
#cat("Normalized cols matrix:\n");
#print(normalized_cols_matrix);

###############################################################################

rows_matrix_samples=rownames(rows_matrix);
cols_matrix_samples=rownames(cols_matrix);

rows_matrix_categories=colnames(rows_matrix);
cols_matrix_categories=colnames(cols_matrix);

rows_matrix_samples_numSamples=length(rows_matrix_samples);
cols_matrix_samples_numSamples=length(cols_matrix_samples);
shared_samples=intersect(rows_matrix_samples, cols_matrix_samples);
num_shared_samples=length(shared_samples);

cat("Num samples in RowsInput: ", rows_matrix_samples_numSamples, "\n", sep="");
cat("Num samples in ColsInput: ", cols_matrix_samples_numSamples, "\n", sep="");
cat("Num shared samples: ", num_shared_samples, "\n", sep="");

# Compute shared samples
sorted_shared_samples=sort(shared_samples);
row_matrix_shrd_idx=numeric(num_shared_samples);
col_matrix_shrd_idx=numeric(num_shared_samples);
for(i in 1:length(sorted_shared_samples)){
	cur_sample=sorted_shared_samples[i];
	row_matrix_shrd_idx[i]=which(cur_sample==rows_matrix_samples);
	col_matrix_shrd_idx[i]=which(cur_sample==cols_matrix_samples);	
}

#print(row_matrix_shrd_idx);
#print(col_matrix_shrd_idx);


# Get number of categories
row_matrix_num_cat=ncol(rows_matrix);
col_matrix_num_cat=ncol(cols_matrix);
correlation_matrix=matrix(0, nrow=row_matrix_num_cat, ncol=col_matrix_num_cat);
pvalue_matrix     =matrix(0, nrow=row_matrix_num_cat, ncol=col_matrix_num_cat);
for(i_cat in 1:row_matrix_num_cat){
	cat("Computing row: ", i_cat, "\n", sep="");
	for(j_cat in 1:col_matrix_num_cat){
		#cat("Working on ", i_cat, " vs ", j_cat, "\n");
		rows_vector=normalized_rows_matrix[row_matrix_shrd_idx,i_cat];
		cols_vector=normalized_cols_matrix[col_matrix_shrd_idx,j_cat];

		#print(rows_vector);
		#print(cols_vector);
		ct=cor.test(rows_vector, cols_vector, method="spearman");
		
		correlation_matrix[i_cat, j_cat]=ct$estimate;
		pvalue_matrix[i_cat, j_cat]=ct$p.value
		#cat("\t", cor, "\n");
	}
}

#print(correlation_matrix);


skewness=function(x){
	mu=mean(x);
	sigma=sd(x);
	gamma=mean(((x-mu)/sigma)^3);
	return(gamma);
}


# Sort correlation matrix
correlation_sort_idx=matrix(0, nrow=row_matrix_num_cat, ncol=col_matrix_num_cat);
correlation_sort_val=matrix(0, nrow=row_matrix_num_cat, ncol=col_matrix_num_cat);
for(i_cat in 1:row_matrix_num_cat){
	sort_rec=sort(abs(correlation_matrix[i_cat,]), decreasing=TRUE, index.return=TRUE);
	correlation_sort_val[i_cat,]=correlation_matrix[i_cat,sort_rec$ix];
	correlation_sort_idx[i_cat,]=sort_rec$ix;
}


if(BonferroniCorrect){
	PvalueCutoff=PvalueCutoff/(col_matrix_num_cat*row_matrix_num_cat);
	cat("Bonferroni Corrected p-value: ", PvalueCutoff, "\n", sep="");
}

#pdf(paste(OutputFileRoot, ".heatmap.pdf", sep=""), height=8.5, width=11);
#image(t(correlation_sort_val),);
fh=file(paste(OutputFileRoot, ".top_correlated.txt", sep=""), "w");
for(i_cat in 1:row_matrix_num_cat){
	cat(rows_matrix_categories[i_cat], "\n", sep="", file=fh);

	for(j_cat in 1:col_matrix_num_cat){

		top_idx=correlation_sort_idx[i_cat,j_cat];
		pval=pvalue_matrix[i_cat,top_idx];
		
		if(pval<PvalueCutoff){
			cat("\t", cols_matrix_categories[top_idx], "\t", 
				correlation_matrix[i_cat,top_idx], "\t", 
				pvalue_matrix[i_cat,top_idx], "\t",
				"\n", file=fh, sep="");

		}else{
			break;
		}
	}
}
close(fh);


# Histogram each column
if(GeneratePlot){
	hist_breaks=seq(-1,1,length.out=12);
	pdf(paste(OutputFileRoot, ".correl_dist_row.pdf", sep=""), height=8.5, width=11);
	for(i_cat in 1:row_matrix_num_cat){
		vals=correlation_matrix[i_cat,];
		skew=skewness(vals);
		med=median(vals);
		hist(vals, xlim=c(-1,1), main=rows_matrix_categories[i_cat], xlab="Correlation Coefficient (rho)", hist_breaks);	
		mtext(sprintf("median: %g", med), line=0, cex=.75);
		mtext(sprintf("skew  : %g", skew), line=-.5, cex=.75);
	}

	pdf(paste(OutputFileRoot, ".correl_dist_col.pdf", sep=""), height=8.5, width=11);
	for(j_cat in 1:col_matrix_num_cat){
		vals=correlation_matrix[,j_cat];
		skew=skewness(vals);
		med=median(vals);
		hist(correlation_matrix[,j_cat], xlim=c(-1,1), main=cols_matrix_categories[j_cat], xlab="Correlation Coefficient (rho)", hist_breaks);	
		mtext(sprintf("median: %g", med), line=0, cex=.75);
		mtext(sprintf("skew  : %g", skew), line=-.5, cex=.75);
	}
}

###############################################################################

writeLines("Done.\n")
print(warnings());

q(status=0)
