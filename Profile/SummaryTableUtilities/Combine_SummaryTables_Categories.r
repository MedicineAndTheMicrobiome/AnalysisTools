#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"summary_table_a", "a", 1, "character",
	"summary_table_b", "b", 1, "character",
	"output_file", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-a <input summary_table.tsv file name A>\n",
	"	-b <input summary_table.tsv file name B>\n",
	"	-o <output file name root>\n",
	"\n",	
	"This script will combine the categories of two summary \n",
	"tables.\n",
	"\n",
	"The script will confirm that the sample ID's match.\n", 
	"\n",
	"\n");

if(
	!length(opt$summary_table_a) ||
	!length(opt$summary_table_b) ||
	!length(opt$output_file)
){
	cat(usage);
	q(status=-1);
}

outputroot=gsub("\\.summary_table\\.xls$", "", opt$output_file);
outputroot=gsub("\\.summary_table\\.tsv$", "", outputroot);
OutputFileName=paste(outputroot, ".summary_table.tsv", sep="");

SummaryTableA=opt$summary_table_a;
SummaryTableB=opt$summary_table_b;

###############################################################################

cat("\n")
cat("Summary Table A: ", SummaryTableA, "\n");
cat("Summary Table B: ", SummaryTableB, "\n");
cat("Output File Name Root: ", OutputFileName, "\n");
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
        counts_mat=inmat[,2:(ncol(inmat)), drop=F];
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

write_summary_table=function(out_mat, fname){
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

###############################################################################
# Load counts matrix

counts_mat_a=load_summary_table(SummaryTableA);
num_categories_a=ncol(counts_mat_a);
num_samples_a=nrow(counts_mat_a);
categories_a=colnames(counts_mat_a);

counts_mat_b=load_summary_table(SummaryTableB);
num_categories_b=ncol(counts_mat_b);
num_samples_b=nrow(counts_mat_b);
categories_b=colnames(counts_mat_b);

###############################################################################

cat("Checking that Sample IDs match...\n");

sample_ids_a=rownames(counts_mat_a);
sample_ids_b=rownames(counts_mat_b);

shared_sample_ids=intersect(sample_ids_a, sample_ids_b);
num_shared_ids=length(shared_sample_ids);

cat("Num Shared IDs: ", num_shared_ids, "\n");

if(num_shared_ids<num_samples_a || num_shared_ids<num_samples_b){
	cat("Warning!  Number of shared samples not complete.\n");
}

###############################################################################

cat("Checking that categories don't overlap...\n");
shared_categories=intersect(categories_a, categories_b);
num_shared=length(shared_categories);
if(num_shared){
	cat("Error:  Categories overlapping:\n");
	print(shared_categories);
	quit(-1);
}else{
	cat("Categories are unique...good.\n");
}

###############################################################################

combined_counts_mat=cbind(counts_mat_a[shared_sample_ids,], counts_mat_b[shared_sample_ids,]);
num_new_categories=ncol(combined_counts_mat);
cat("Num Combined Categories: ", num_new_categories, "\n");

# Output
cat("\nWriting New Matrix...\n");
write_summary_table(combined_counts_mat, OutputFileName);

###############################################################################

cat("\nDone.\n")
warns=warnings();
if(!is.null(warns)){
	print(warnings());
}
q(status=0)
