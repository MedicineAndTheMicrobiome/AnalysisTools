#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"sample_id", "s", 1, "character",
	"input_file", "i", 1, "character",
	"output_file", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-s <sample id, to use as row/sample name>\n",
	"	-i <bracken input file name, from kraken/bracken>\n",
	"	[-o <output file name e.g. .summary_table.tsv>]\n",
	"\n",	
	"This script will convert a bracken file with reads/abundances\n",
	"into a single entry (row) summary_file.tsv file.\n",
	"\n",
	"This script expects to find 'name' and 'new_est_reads' column names.\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

SampleID=opt$sample_id;
InputFileName=opt$input_file;
OutputFileName=opt$output_file;

cat("\n");
cat("Sample ID: ", SampleID, "\n", sep="");
cat("Input File: ", InputFileName, "\n", sep="");
cat("Output File: ", OutputFileName, "\n", sep="");

###############################################################################

data=read.table(InputFileName, header=T, sep="\t");

# Extract and sort
count_data=data[,c("name", "new_est_reads")];
dec_ord=order(count_data[,"new_est_reads"], decreasing=T);
count_data=count_data[dec_ord,];

cat("\n");
cat("Top 10 Taxa:\n");
print(count_data[1:10,]);

cat("\n");
cat("Converting to summary table...\n");
num_categories=nrow(count_data);
cat_names=count_data[,"name"];
cat_counts=count_data[,"new_est_reads"];

# Create 1 row summary table
out_mat=matrix(NA, nrow=1, ncol=num_categories);
colnames(out_mat)=cat_names;
rownames(out_mat)=SampleID;
out_mat[1,]=cat_counts;

#print(out_mat);

###############################################################################

cat("\nWriting New Matrix...\n");

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

write_summary_file(out_mat, OutputFileName);

###############################################################################

cat("\nDone.\n")
warns=warnings();
if(!is.null(warns)){
        print(warnings());
}
q(status=0)

