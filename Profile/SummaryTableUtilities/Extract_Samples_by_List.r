#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"sample_list", "s", 1, "character",
	"column_name", "c", 2, "character",
	"remove_samples", "R", 2, "logical",
	"output_file", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary_table.tsv>\n",
	"	-s <sample list>\n",
	"	[-c <column name>]\n",
	"	[-R remove samples]\n",
	"	-o <output summary_table file name>\n",
	"\n",	
	"This script will extract the samples on the list from the summary table,\n",
	"writing them to a new summary table.\n",
	"\n",
	"If the -R flag is set, the samples on the list will instead\n",
	"be deleted and the remaining samples written to a the new summary table\n",
	"\n",
	"If the -c column is set, then the sample list file can be a matrix\n",
	"and the column with the sample IDs of interest is specified by name.\n",
	"\n");

if(!length(opt$input_file) || !length(opt$sample_list) || !length(opt$output_file)){
	cat(usage);
	q(status=-1);
}

RemoveSamples=F;
if(length(opt$remove_samples)){
	RemoveSamples=T;
}

ColumnName="";
if(length(opt$column_name)){
	ColumnName=opt$column_name;
}

###############################################################################

InputFileName=opt$input_file;
SampleList=opt$sample_list;
OutputFileName=opt$output_file;

OutputFileName=gsub("\\.summary_table\\.tsv", "", OutputFileName);
OutputFileName=gsub("\\.summary_table\\.xls", "", OutputFileName);
OutputFileName=paste(OutputFileName, ".summary_table.tsv", sep="");


cat("\n");
cat("Input File Name: ", InputFileName, "\n");
cat("Sample List: ", SampleList, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
if(ColumnName!=""){
	cat("Column Name: ", ColumnName, "\n");
}
cat("\n");

if(InputFileName==OutputFileName){
	cat("Error: Input and output summary table file name are the same.\n");
}

###############################################################################
###############################################################################

# Load data
inmat=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="", row.names=1))
#cat("Original Matrix:\n")
#print(inmat);

counts_mat=inmat[,2:(ncol(inmat))];

# Summary what we've loaded
num_samples=nrow(inmat);
num_categories=ncol(counts_mat);
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");

###############################################################################
# Load sample list

if(ColumnName==""){
	# If no column name specified, assume it's just a simple list without a header
	sample_list=as.vector(read.table(SampleList, header=FALSE, comment.char="#", sep="\t")[,1]);
}else{
	# If the column name is specified, then there must be a header.
	samp_list_mat=read.table(SampleList, header=T, comment.char="#", sep="\t");
	sample_list=samp_list_mat[,ColumnName];
}

num_samples_to_extract=length(sample_list);

st_samples=rownames(inmat);

if (RemoveSamples){
	cat("Samples to Remove:\n");
	print(sample_list);
	cat("\n");
        cat("Num Samples to Remove: ", num_samples_to_extract, "\n");
	to_extract=setdiff(st_samples, sample_list);
	
}else{
	cat("Samples to Extract:\n");
	print(sample_list);
	cat("\n");
	cat("Num Samples to Extract: ", num_samples_to_extract, "\n");
	to_extract=intersect(st_samples, sample_list);
}

num_to_extract=length(to_extract);
num_targeted=length(sample_list);
if(num_to_extract<num_targeted){
	cat("Warning: ", num_targeted, " samples targeted.  ", num_to_extract, " samples found.\n", sep="");
}else{
	cat("All samples found...\n");
}

outmat=counts_mat[to_extract,];

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
write_summary_file(outmat, OutputFileName);

###############################################################################

writeLines("Done.\n")
print(warnings());

q(status=0)
