#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"group_map", "m", 1, "character",
	"output_file", "o", 2, "character", 
	"acc_in_remaining", "r", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

AbundanceCutoff=1;

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.tsv>\n",
	"	-m <grouping file>\n",
	"	[-o <output regrouped summary table file name root>]\n",
	"	[-r <move ungrouped categories into remaining>]\n",
	"\n",	
	"This script will merge the categories in the summary table\n",
	"based on the groups in the grouping file.\n",
	"\n",
	"The contents of the grouping file should be:\n",
	"	<category name>\\t<new group name>\\t[optional comments]\\n\n",
	"\n",
	"There shouldn't be a header in this file.\n",
	"\n",
	"If the -r option is not specified:\n",
	"  If the category is not found in the group file, then it will remain\n",
	"  the same in the output file.\n",
	"\n",
	"Else if the -r is specified, then\n",
	"  Anything remaining (including the counts in Remaining) will be summed up and placed into Remaining\n",
	"\n", sep="");

if(!length(opt$input_file) || !length(opt$group_map)){
	cat(usage);
	q(status=-1);
}else{
	InputFileName=opt$input_file;
	GroupFileName=opt$group_map;
}

if(!length(opt$output_file)){
	OutputFileName = gsub("\\.summary_table\\.xls$", "", opt$input_file);
	OutputFileName = gsub("\\.summary_table\\.tsv$", "", OutputFileName);
	OutputFileName = paste(OutputFileName, ".regrpd", sep="");
}else{
	OutputFileName = opt$output_file;
}

AccumulateInRemaining=logical();
if(!length(opt$acc_in_remaining)){
	AccumulateInRemaining=TRUE;
}else{
	AccumulateInRemaining=FALSE;
}

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Group Map Filename: ", GroupFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("Accumulate In Remaining: ", AccumulateInRemaining, "\n");
cat("\n");

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
category_names=colnames(counts_mat);

cat("\n");
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");
cat("\n");

###############################################################################

# Load map
cat("Loading Map...\n");
map_matrix=as.matrix(read.delim(GroupFileName, sep="\t", header=F, check.names=FALSE, comment.char="", quote=""))[,c(1,2)];
colnames(map_matrix)=c("Original Name", "Group");
print(map_matrix);

unique_groups=sort(unique(map_matrix[,2]));
num_groups=length(unique_groups);
target_categories=sort(map_matrix[,1]);

unmapped_categories=sort(setdiff(category_names, target_categories));
num_unmapped_cat=length(unmapped_categories);

cat("\n");
cat("Number of Target Groups: ", num_groups, "\n");
cat("Number of Unmapped Categories: ", num_unmapped_cat, "\n");


if(AccumulateInRemaining){
	cat("Accumulating unmapped into Remaining...\n");
	new_matrix=matrix(apply(counts_mat[, unmapped_categories, drop=F], 1, sum), ncol=1, nrow=num_samples);	
	colnames(new_matrix)=c("Remaining");
}else{
	cat("Copying unmapped categories into new matrix...\n");
	# Copy categories that won't be touched to new matrix
	new_matrix=counts_mat[, unmapped_categories, drop=F];
}

map_groups=map_matrix[,2];
map_categories=map_matrix[,1];
for(group in unique_groups){
	cat("\nCollapsing ", group, "...\n", sep="");
	categories_in=map_categories[map_groups==group];
	#print(categories_in);
	overlapping=intersect(categories_in, category_names);
	extracted=counts_mat[, overlapping, drop=F];
	cat("Found:\n");
	print(overlapping);
	sums=matrix(apply(extracted, 1, sum), dimnames=list(sample_names, group), ncol=1);
	#print(sums)	
	new_matrix=cbind(sums, new_matrix);
}

outmat=new_matrix;

###############################################################################

# Output Filtered Summary Table
cat("\n\nWriting New Matrix...\n");
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

if(!is.null(warnings())){
	print(warnings());
}

cat("Done.\n")

q(status=0)
