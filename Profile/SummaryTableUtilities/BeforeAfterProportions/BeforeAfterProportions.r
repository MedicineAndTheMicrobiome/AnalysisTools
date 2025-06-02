#!/usr/bin/env Rscript

###############################################################################

library('getopt');
source("~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r");

params=c(
	"input_file_B", "B", 1, "character",
	"input_file_A", "A", 1, "character",
	"output_file", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-B <Before input summary table.xls>\n",
	"	-A <Output input summary table.xls>\n",
	"	-o <output file name\n",
	"\n",	
	"This script will read in both files, match the counts by sample ID,\n",	
	"then for each of the sample IDs that paired properly, calculate the\n",
	"the proportion of reads that remain.\n",
	"For example:\n",	
	"	if Before you had 100 reads, then After you have 30 reads\n",
	"	then the proportion reported will be 0.30.\n",
	"\n",
	"\n");

if(
	!length(opt$input_file_B) || 
	!length(opt$input_file_A) || 
	!length(opt$output_file)){
	cat(usage);
	q(status=-1);
}


InputFileNameBefore=opt$input_file_B;
InputFileNameAfter=opt$input_file_A;
OutputFileName=opt$output_file;

###############################################################################

cat("\n")
cat("Input File Name Before: ", InputFileNameBefore, "\n");
cat("Input File Name After: ", InputFileNameAfter, "\n");
cat("Output File Name Name: ", OutputFileName, "\n");       
cat("\n");

###############################################################################
###############################################################################

before_st=load_summary_file(InputFileNameBefore);
after_st=load_summary_file(InputFileNameAfter);

before_sample_ids=rownames(before_st);
after_sample_ids=rownames(after_st);

num_before_sample_ids=length(before_sample_ids);
num_after_sample_ids=length(after_sample_ids);

cat("Num Before Sample IDs: ", num_before_sample_ids, "\n");
cat("Num After Sample IDs: ", num_after_sample_ids, "\n");

shared_sample_ids=intersect(before_sample_ids, after_sample_ids);
num_shared_sample_ids=length(shared_sample_ids);
cat("Num Shared Sample IDs: ", num_shared_sample_ids, "\n");

shared_before=before_st[shared_sample_ids,,drop=F];
shared_after=after_st[shared_sample_ids,,drop=F];

shared_before_counts=apply(shared_before, 1, sum);
shared_after_counts=apply(shared_after, 1, sum);

counts_mat=cbind(shared_before_counts, shared_after_counts);
print(counts_mat);

proportions=apply(counts_mat, 1, function(x){ x[2]/x[1];});

print(proportions);

out_tab=cbind(names(proportions), proportions);
colnames(out_tab)=c("SampleID", "ProportionKept");
write.table(out_tab, file=OutputFileName, sep="\t", col.names=T, row.names=F, quote=F);


###############################################################################

writeLines("Done.\n")
warns=warnings();
if(length(warns)>0){
	print(warnings());
}

q(status=0)
