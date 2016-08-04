#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"remove_list", "l", 1, "character",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	-l <list of categories to remove>\n",
	"	[-o <output summary table file name>]\n",
	"\n",	
	"This script will read in a summary table and remove the categories\n",
	"that are specified in the list file.\n",
	"\n",
	"\n");

if(!length(opt$input_file) || !length(opt$remove_list)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$output_file)){
	outputroot=gsub("\\.summary_table\\.xls", "", opt$input_file);
	outputroot=gsub("\\.summary_table\\.tsv", "", opt$input_file);
	OutputFileName = paste(outputroot, ".list_filtered.summary_table.tsv", sep="");
}else{
	OutputFileName=opt$output_file;
}

###############################################################################

InputFileName=opt$input_file;
RemoveList=opt$remove_list;

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Remove List Name: ", RemoveList, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("\n");

###############################################################################
###############################################################################

load_summary_table=function(summary_table_fn){
        # Load data
        cat("Loading Matrix (", summary_table_fn, ") ...\n", sep="");
        inmat=as.matrix(read.table(summary_table_fn, sep="\t", header=TRUE, check.names=FALSE, row.names=1, quote=""))

        #cat("\nOriginal Matrix:\n")
        #print(inmat);

        # Grab columns we need into a vector, ignore totals, we won't trust it.
        counts_mat=inmat[,2:(ncol(inmat))];
        #cat("\nCounts Matrix:\n");
        #print(counts_mat);

        num_samples=nrow(counts_mat);
        num_categories=ncol(counts_mat);
        sample_names=rownames(counts_mat);

        cat("\n");
        cat("Num Samples: ", num_samples, "\n");
        cat("Num Categories: ", num_categories, "\n");
        cat("\n");
        return(counts_mat);
}

###############################################################################

load_ids=function(list_fn){
        cat("Loading List (", list_fn, ") ...\n", sep="");
        list=scan(file=list_fn, what="complex", sep="\t");
        return(list);
}

###############################################################################
# Load counts matrix
counts_mat=load_summary_table(InputFileName);
num_categories=ncol(counts_mat);
num_samples=nrow(counts_mat);

# Load remove list IDs
remove_id_list=load_ids(RemoveList);
num_removal=length(remove_id_list);
cat("Num categories to remove: ", num_removal, "\n");

# Get the category names
category_names=colnames(counts_mat);

# Get indices for the columns we want to remove
rem_idx=num_removal;
cat("\nRemoval Categories:\n");
for(i in 1:num_removal){
	cat("\t", remove_id_list[i]);
	rem_idx[i]=which(remove_id_list[i]==category_names);
	cat(" / ", rem_idx[i], "\n", sep="");
}

# Remove the columns
outmat=counts_mat[,-rem_idx];

# Remove num categories left
num_remaining_categories=ncol(outmat);
cat("Num remaining categories: ", num_remaining_categories, "\n");

###############################################################################
# Output
cat("\nWriting New Matrix...\n");
fc=file(OutputFileName, "w");

write(paste("sample_id", "total", paste(colnames(outmat), collapse="\t"), sep="\t"), file=fc);
sample_names=rownames(counts_mat);
for(samp_idx in 1:num_samples){
	total=sum(outmat[samp_idx,]);
	outline=paste(sample_names[samp_idx],total,paste(outmat[samp_idx,], collapse="\t"), sep="\t");
	write(outline, file=fc);
}
close(fc);	

###############################################################################

cat("Done.\n")
warns=warnings();
if(!is.null(warns)){
	print(warnings());
}
q(status=0)
