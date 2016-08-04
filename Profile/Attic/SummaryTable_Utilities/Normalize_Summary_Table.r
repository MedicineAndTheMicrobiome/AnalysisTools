#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_fname_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table>\n",
	"	[-o <output filename root>]\n",
	"\n",	
	"This script will normalize a summary table's counts into\n",
	"proportions.\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileNameRoot=opt$output_fname_root;

if(length(opt$output_fname_root)){
	OutputFileNameRoot=opt$output_fname_root;
}else{
	OutputFileNameRoot=gsub("\\.summary_table.tsv$", "", InputFileName);
	OutputFileNameRoot=gsub("\\.summary_table.xls$", "", OutputFileNameRoot);
}

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name Root: ", OutputFileNameRoot, "\n");       

###############################################################################
###############################################################################

# Load data
inmat=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, row.names=1, quote=NULL))
#cat("Original Matrix:\n")
#print(inmat);

# Grab columns we need into a vector, ignore totals, we won't trust it.
counts_mat=inmat[,2:(ncol(inmat)), drop=F];

num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);
cat("\n");
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");
cat("\n");


sums=apply(counts_mat, 1, sum);

new_summary_table=matrix(0, nrow=num_samples, ncol=num_categories);
rownames(new_summary_table)=rownames(counts_mat);
colnames(new_summary_table)=colnames(counts_mat);
for(i in 1:num_samples){
	new_summary_table[i,]=counts_mat[i,]/sums[i];
}

###############################################################################

write_summary_file=function(out_mat, fname){
        fc=file(fname, "w");
        cat(file=fc, paste("sample_id\ttotal", paste(colnames(out_mat), collapse="\t"), sep="\t"));
        cat(file=fc, "\n");
        sample_names=rownames(out_mat);
        num_samples=nrow(out_mat);
        for(samp_idx in 1:num_samples){
                total=sum(out_mat[samp_idx,]);
                outline=paste(sample_names[samp_idx], total,
                        paste(out_mat[samp_idx,], collapse="\t"), sep="\t");
                cat(file=fc, outline);
                cat(file=fc, "\n");
        }
        close(fc);
}

# Output new summary table
output_fname=paste(OutputFileNameRoot, ".norm.summary_table.tsv", sep="");

write_summary_file(new_summary_table, output_fname);

###############################################################################

cat("\nDone.\n");
if(!is.null(warnings())){
	print(warnings());
}

q(status=0);
