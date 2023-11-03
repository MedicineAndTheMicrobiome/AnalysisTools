#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <filename for list of items, eg. .txt, .lst, .tsv>\n",
	"	[-o <output file name e.g. .summary_table.tsv>]\n",
	"\n",	
	"This script will read in a list of items, count them and output\n",
	"a single entry summary table that can be merged with another script.\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileName=InputFileName;

if(!length(opt$output_file)){
	OutputFileName=gsub(".txt$", "", OutputFileName);
	OutputFileName=gsub(".lst$", "", OutputFileName);
	OutputFileName=gsub(".tsv$", "", OutputFileName);
	OutputFileName=paste(OutputFileName, ".summary_table.tsv", sep="");
}else{
	OutputFileName=opt$output_file;
}

cat("Input File: ", InputFileName, "\n", sep="");
cat("Output File: ", OutputFileName, "\n", sep="");

###############################################################################

items=scan(InputFileName, what=character(), sep="\n");

counts=sort(table(items), decreasing=T);

cat_names=names(counts);
cat_names=gsub("^-$", "Unclassified", cat_names);

cat_names=gsub(" ", "_", cat_names);
cat_names=gsub("-", ".", cat_names);
cat_names=gsub("=", ".", cat_names);
cat_names=gsub("\\(", ".", cat_names);
cat_names=gsub("\\)", ".", cat_names);
cat_names=gsub("/", ".", cat_names);
cat_names=gsub(":", ".", cat_names);

names(counts)=cat_names;

###############################################################################

cat("\nWriting New Matrix...\n");

out_mat=matrix(counts, ncol=length(counts), nrow=1);
colnames(out_mat)=names(counts);
rownames(out_mat)=gsub(".summary_table.tsv", "", OutputFileName);

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

