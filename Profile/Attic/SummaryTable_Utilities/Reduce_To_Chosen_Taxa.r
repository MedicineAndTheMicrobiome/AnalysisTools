#!/usr/bin/env Rscript

###############################################################################

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
	"	-t <list of taxa to keep>\n",
	"	[-o <output reduced summary table file name root>]\n",
	"\n",	
	"This script will read in a summary table, and a list of chosen\n",
	"taxa, and then reduce the number of categories to the number of\n",
	"chosen taxa, plus one for the remaining.\n",
	"\n",
	"If the remaining category already exists, then the category\n",
	"counts for those not selected will be added to the remaining column\n",
	"\n", sep="");

if(!length(opt$input_file) || !length(opt$taxa_list)){
	cat(usage);
	q(status=-1);
}else{
	InputFileName=opt$input_file;
	TaxaListName=opt$taxa_list;
}

if(!length(opt$output_file)){
	OutputFileName = gsub("\\.summary_table\\.xls$", "", opt$input_file);
	OutputFileName = gsub("\\.summary_table\\.tsv$", "", OutputFileName);
	OutputFileName = paste(OutputFileName, ".reduced", sep="");
}else{
	OutputFileName = opt$output_file;
}

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Taxa List: ", TaxaListName, "\n");
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

# Load data
cat("Loading Matrix...\n");
inmat=as.matrix(read.delim(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="", quote="", row.names=1))
#cat("Original Matrix:\n")
#print(inmat);

# Grab columns we need into a vectord , ignore totals, we won't trust it.
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

# Load list
TaxaList=scan(TaxaListName, what=character(), sep="\n", quote="", comment.char="");
num_targeted=length(TaxaList);
#print(TaxaList);

###############################################################################

# Confirming remaining column is not in keep list
uc_taxalist=toupper(TaxaList);
keep_list_remainder_ix=which(uc_taxalist=="REMAINDER" | uc_taxalist=="REMAINING");
if(length(keep_list_remainder_ix)>0){
	cat("*********************************************************************************\n");
	cat("*                                                                               *\n");
	cat("*   WARNING: Remainder/Remaining column specified in Category in keep list!!!   *\n");
	cat("*                                                                               *\n");
	cat("*********************************************************************************\n");
	cat("\n");
	cat("Removing:\n");
	print(TaxaList[keep_list_remainder_ix]);
	cat("\n");
	TaxaList=TaxaList[-keep_list_remainder_ix];
}

# Identifying remaining column in counts table
categories=colnames(counts_mat);
uc_categories=toupper(categories);
count_mat_remainder_ix=which(uc_categories=="REMAINDER" | uc_categories=="REMAINING");
if(length(count_mat_remainder_ix)>1){
	cat("***************************************************************************\n");
	cat("*                                                                         *\n");
	cat("*  ERROR:  More than one remainder/remaining category.  Too confusing!!!  *\n");
	cat("*                                                                         *\n");
	cat("***************************************************************************\n");
	quit(status=-1);
}

# Extracting values for remaining column
remaining_counts=counts_mat[,count_mat_remainder_ix, drop=F];
if(!length(count_mat_remainder_ix)){
	cat("Remaining/Remainder counts not found...  Assuming 0's.\n\n");
	remaining_counts=matrix(0, nrow=num_samples, ncol=1);
	colnames(remaining_counts)=REMAINING_TAXA_NAME;
	rownames(remaining_counts)=sample_names;
}
cat("Current remaining counts:\n");
print(remaining_counts);
cat("\n");

# Extract columns that have been requested
shared=sort(intersect(categories, TaxaList));
missing_targets=F;
#print(shared);
num_shared=length(shared);
cat("Num Shared: ", num_shared, "\n");
cat("Num Targeted: ", num_targeted, "\n");
if(num_shared<num_targeted){
	missing_targets=T;
}
#cat("Identified categories: \n");
#print(shared);
#cat("\n");
reduced_mat=counts_mat[,shared, drop=F];
colnames(reduced_mat)=shared;

# Compute total counts (includes remaining) and counts for kept categories (excludes remaining)
totals=apply(counts_mat, 1, sum);
kept=apply(reduced_mat, 1, sum);

# Compute what's in the cumulative remainder column
remaining=totals-kept;

# Attached remained to end of matrix
reduced_matrix=cbind(reduced_mat, remaining);

# Rename the remainder matrix column
colnames(reduced_matrix)=c(shared, REMAINING_TAXA_NAME);
outmat=reduced_matrix;

###############################################################################

# Output Filtered Summary Table
cat("Writing New Matrix...\n");
fc=file(paste(OutputFileName, ".summary_table.tsv", sep=""), "w");

# Write header
cat(file=fc, paste("Sample_ID", "Total", paste(colnames(outmat), collapse="\t"), sep="\t"));
cat(file=fc, "\n");

# Write counts
for(samp_idx in 1:num_samples){
	total=sum(outmat[samp_idx,]);
	cat(file=fc, paste(sample_names[samp_idx],total,paste(outmat[samp_idx,], collapse="\t"), sep="\t"));
	cat(file=fc, "\n");
}
close(fc);	

###############################################################################

if(missing_targets){
	cat("WARNING: Targets specified, but not identified for extraction.\n");
}

if(!is.null(warnings())){
	print(warnings());
}

cat("Done.\n")

q(status=0)
