#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

AbundanceCutoff=1;

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input Mothur shared file>\n",
	"	-o <output filename root>\n",
	"\n",	
	"This script will read in the .shared file from Mothur's output of make.shared.\n",
	"The output will be a summary file table\n",
	"\n",
	"\n", sep="");

if(!length(opt$input_file) || !length(opt$output_file)){
	cat(usage);
	q(status=-1);
}else{
	InputFileName=opt$input_file;
	OutputFileNameRoot=opt$output_file;
}

OutputFileNameRoot=gsub("\\.summary_table\\.tsv$", "", OutputFileNameRoot);
OutputFileNameRoot=gsub("\\.summary_table\\.xls$", "", OutputFileNameRoot);

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name Root: ", OutputFileNameRoot, "\n");       
cat("\n");

###############################################################################
###############################################################################

# Load data
cat("Loading Matrix...\n");
inmat=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, 
	check.names=FALSE, comment.char="", quote="")); 

#cat("Original Matrix:\n")
#print(inmat);

# Check for label/cutoff consistency
label=inmat[,"label"];
label=unique(label);
if(length(label)!=1){
	cat("ERROR: There is more than one label/cutoff in this file!!!\n");
	print(label);
	quit(status=-1);
}
cat("Label/Cutoff: ", label, "\n", sep="");

# Get sample IDs
sample_names=inmat[,"Group"]
cat("Sample Names identified: \n");
print(sample_names);
cat("\n");

# Get number of OTUs
num_otus=inmat[,"numOtus"];
num_otus=unique(num_otus);
if(length(num_otus)>1){
	cat("ERROR: There is more than one OTU grouping in this file!!!\n");
	print(num_otus);
	quite(status=-1);
}
num_otus=as.numeric(num_otus);
cat("Number of OTUs: ", num_otus, "\n", sep="");

# Extract OTU counts into count matrix and convert to numeric
counts_char=(inmat[,4:(num_otus+3)]);
rownames(counts_char)=sample_names;

outmat=counts_char;
class(outmat)="numeric";
#print(counts);

###############################################################################

num_samples=nrow(outmat);
num_categories=ncol(outmat);
sample_names=rownames(outmat);

cat("\n");
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");
cat("\n");

###############################################################################

# Output Filtered Summary Table
cat("Writing New Matrix...\n");
out_label=gsub("^0", "", sprintf("%g", (1-as.numeric(label))));
fc=file(paste(OutputFileNameRoot, out_label, ".summary_table.tsv", sep=""), "w");

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

cat("Done.\n")
if(!is.null(warnings())){
	print(warnings());
}

q(status=0)
