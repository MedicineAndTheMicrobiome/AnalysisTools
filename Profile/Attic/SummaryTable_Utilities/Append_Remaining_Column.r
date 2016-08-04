#!/usr/bin/env Rscript

###############################################################################

library('getopt');

REMAINING_CATEGORY_NAME="Remaining";

params=c(
	"input_file", "i", 1, "character",
	"append_file", "t", 1, "character",
	"value", "v", 1, "character",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

AbundanceCutoff=1;

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary_table.tsv>\n",
	"	-t <input append values file (tsv)>\n",
	"	-v <value type: (R)emainder, (T)otal>\n",
	"	[-o <output appended summary table file name root>]\n",
	"\n",	
	"This script will read in a summary table and then a values\n",
	"file.  The values file will be appended to the columns of\n",
	"the summary table.\n",
	"\n",
	"If the value type is R, for Remainder, then the column\n",
	"will just be appended, and the Total column will be updated\n",
	"to reflect the extra counts.\n",
	"\n",
	"If the value type is T, for Total, then this implies that\n",
	"we want the Total to be this value, so a Remainder column\n",
	"will be created to fill in the difference.\n",
	"\n",
	"\n", sep="");

if(!length(opt$input_file) || !length(opt$value) || !length(opt$append_file)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$output_file)){
	OutputFileName = gsub("\\.summary_table\\.xls$", "", opt$input_file);
	OutputFileName = gsub("\\.summary_table\\.tsv$", "", OutputFileName);
}else{
	OutputFileName = opt$output_file;
}

InputFileName=opt$input_file;
AppendFile=opt$append_file
Value=toupper(opt$value);

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Append Values File Name: ", AppendFile, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("Append Value: ", Value, "\n");
cat("\n");

###############################################################################
###############################################################################

# Load summary file data
cat("Loading Matrix...\n");
inmat=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="*", quote="", row.names=1))
#cat("Original Matrix:\n")
#print(inmat);

# Grab columns we need into a vector, ignore totals, we won't trust it.
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

# Load append values
cat("Loading Values...\n");
append_values=as.matrix(read.table(AppendFile, sep="\t", check.names=FALSE, comment.char="*", quote="", row.names=1))
#print(append_values);
colnames(append_values)=c("Append Values");
num_app_samples=nrow(append_values);
app_samp_names=rownames(append_values);
cat("Number of Samples in Append Values File: ", num_app_samples, "\n");

###############################################################################

shared_samples=sort(intersect(app_samp_names, sample_names));

cat("\nShared Samples:\n");
print(shared_samples);

if(length(shared_samples)<num_samples){
	cat("WARNING: Shared samples is less than number of samples in summary table.\n");
}

shared_app_val=append_values[shared_samples,, drop=F];
out_mat=counts_mat[shared_samples,];

cat("\n");
print(shared_app_val);
count_tot=apply(out_mat, 1, sum);

cat("\n");
if(Value=="T"){
	cat("Treating specified values as the Total...\n");
	remainder=shared_app_val-count_tot;
}else if(Value=="R"){
	cat("Treating specified value as the Remaining...\n");
	remainder=shared_app_val;
}else{
	cat("Unrecognized Value Type: ", Value, "\n");
}
cat("\n");

if(any(category_names==REMAINING_CATEGORY_NAME)){
	cat("*********************************************************************\n");
	cat("*  ERROR: ", REMAINING_CATEGORY_NAME, " category already exists!!!  *\n");
	cat("*********************************************************************\n");
	quit(status=-1);
}


colnames(remainder)=REMAINING_CATEGORY_NAME;
out_mat=cbind(out_mat, remainder);

# Output results for sanity check
outcome=matrix(0, nrow=length(shared_samples), ncol=2);
rownames(outcome)=rownames(out_mat);
colnames(outcome)=c("Total", REMAINING_CATEGORY_NAME);
outcome[,1]=apply(out_mat, 1, sum);
outcome[,2]=out_mat[,REMAINING_CATEGORY_NAME];
print(outcome);

###############################################################################

# Output Filtered Summary Table
cat("\nWriting New Matrix...\n");
fc=file(paste(OutputFileName, ".wRem.summary_table.tsv", sep=""), "w");

# Write header
cat(file=fc, paste("Sample_ID", "Total", paste(colnames(out_mat), collapse="\t"), sep="\t"));
cat(file=fc, "\n");

# Write counts
for(samp_idx in 1:num_samples){
	total=sum(out_mat[samp_idx,]);
	cat(file=fc, paste(sample_names[samp_idx],total,paste(out_mat[samp_idx,], collapse="\t"), sep="\t"));
	cat(file=fc, "\n");
}
close(fc);	

###############################################################################

cat("Done.\n")
if(!is.null(warnings())){
	print(warnings());
}

q(status=0)
