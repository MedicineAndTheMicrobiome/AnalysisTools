#!/usr/bin/env Rscript

library('getopt');

REMAINING_TAXA_NAME="Remaining";

params=c(
	"input_file", "i", 1, "character",
	"taxa_list", "t", 1, "character",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

AbundanceCutoff=1;

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	-t <file name for list of taxa to keep>\n",
	"	[-o <output reduced summary table file name>]\n",
	"\n",	
	"This script will read in a summary table, and a list of chosen\n",
	"taxa, and then reduce the number of categories to the number of\n",
	"chosen taxa, plus one for the remaining.\n",
	"\n",
	"\n", sep="");

if(!length(opt$input_file) || !length(opt$taxa_list)){
	cat(usage);
	q(status=-1);
}else{
	InputFileName=opt$input_file;
}

TaxaList=opt$taxa_list;

if(!length(opt$output_file)){
	taxa_list_name = tail(strsplit(TaxaList, "/")[[1]],1);
	
	OutputFileName = gsub("\\.summary_table\\.xls", "", opt$input_file);
	OutputFileName = gsub("\\.summary_table\\.tsv", "", OutputFileName);
	OutputFileName = paste(OutputFileName, ".", taxa_list_name, ".reduced.summary_table.tsv", sep="");
}else{
	OutputFileName = opt$output_file;
	OutputFileName = gsub("\\.summary_table\\.xls", "", OutputFileName);
	OutputFileName = gsub("\\.summary_table\\.tsv", "", OutputFileName);
	OutputFileName = paste(OutputFileName, ".summary_table.tsv", sep="");
}

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("\n");

## Report what we plan on reducing to
#cat("Taxa List: \n");
#for(i in 1:TaxaListLength){
#	cat("  ", TaxaList[i], "\n", sep="");
#}

#cat("\n");

###############################################################################
###############################################################################

load_ids=function(list_fn){
        cat("Loading List: ", list_fn, "\n", sep="");
        list=as.vector(read.table(list_fn, sep="\t")[,1]);
        return(list);
}

# Load the categories to keep
taxa_ids=load_ids(TaxaList);
num_taxa_ids=length(taxa_ids);

# Load data
cat("Loading Matrix...\n");
inmat=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, 
	check.names=FALSE, comment.char="", quote="", row.names=1))
#cat("Original Matrix:\n")
#print(inmat);

# Grab columns we need into a vector, ignore totals, we won't trust it.
counts_mat=inmat[,2:(ncol(inmat))];
#print(counts_mat);

num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);
sample_names=rownames(counts_mat);

cat("\n");
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");
cat("\n");

###############################################################################

categories=colnames(counts_mat);
#print(categories);

shared=sort(intersect(categories, taxa_ids));
# cat("Identified categories: \n");
# print(shared);

if(length(shared) != length(taxa_ids)){
	cat("*************************************************************************************\n");
	cat("* WARNING: Categories specified in extraction list were not found in Summary Table. *\n");
	cat("*************************************************************************************\n");
	cat("These are the missing categories that were requested:\n");
	print(setdiff(taxa_ids, categories));
	cat("\n");
	warning("Categories specified in extraction list were not found in Summary Table.");
}else{
	cat("Great, all targeted categories identified in Summary Table!\n");
}

reduced_mat=t(t(counts_mat[,shared]));
colnames(reduced_mat)=shared;

totals=apply(counts_mat, 1, sum);
kept=apply(reduced_mat, 1, sum);

# Compute what's in the cumulative remainder column
remaining=totals-kept;

# Attached remained to end of matrix
reduced_matrix=cbind(reduced_mat, remaining);

# Rename the remainder matrix column
colnames(reduced_matrix)=c(shared, REMAINING_TAXA_NAME);

#print(totals);
#print(reduced_matrix);
outmat=reduced_matrix;

###############################################################################

# Output Filtered Summary Table
cat("Writing New Matrix...\n");
fc=file(OutputFileName, "w");

# Write header
cat(file=fc, paste("sample_id", "total", paste(colnames(outmat), collapse="\t"), sep="\t"));
cat(file=fc, "\n");

# Write counts
for(samp_idx in 1:num_samples){
	total=sum(outmat[samp_idx,]);
	cat(file=fc, paste(sample_names[samp_idx],total,paste(outmat[samp_idx,], collapse="\t"), sep="\t"));
	cat(file=fc, "\n");
}
close(fc);	

###############################################################################

cat("Done.\n")
if(!is.null(warnings())){
	print(warnings());
}

q(status=0)
