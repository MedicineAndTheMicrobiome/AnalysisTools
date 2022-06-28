#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 2, "character",
	"shorten_category_names", "x", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.tsv>\n",
	"	-x \"<delimitor/separator character>\"\n",
	"	[-o <output summary table root name>]\n",
	"\n",	
	"This script will read in a summary file and \n",
	"shorten all the categories by selecting the\n",
	"right most token once split by the specified\n",
	"separator character.\n",
	"\n");

if(!length(opt$input_file) || !length(opt$shorten_category_names)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$output_file)){
	outputroot=gsub("\\.summary_table\\.xls", "", opt$input_file);
	outputroot=gsub("\\.summary_table\\.tsv", "", opt$input_file);
	OutputFileName = paste(outputroot, ".shortened.summary_table.tsv", sep="");
}else{
	outputroot=gsub("\\.summary_table\\.xls", "", opt$output_file);
	outputroot=gsub("\\.summary_table\\.tsv", "", opt$output_file);
	OutputFileName=opt$output_file;
}

###############################################################################

InputFileName=opt$input_file;
ShortenCategoryNames=opt$shorten_category_names;

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
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

# Load counts matrix
counts_mat=load_summary_table(InputFileName);
num_categories=ncol(counts_mat);
num_samples=nrow(counts_mat);

# Get the category names
category_names=colnames(counts_mat);

splits=strsplit(category_names, ShortenCategoryNames);
short_names=character();

targets=1:num_categories;
shorten=T;

level=1;
while(shorten){

	cat("\n\nLevel :", level, "\n");

	for(i in targets){
		short_names[i]=paste(tail(splits[[i]], level), collapse=ShortenCategoryNames);
		cat(category_names[i], " ->\n\t",  short_names[i], "\n");
	}

	dup_bool=table(short_names)>1;
	dups_names=names(dup_bool[dup_bool==T]);
	cat("\nDuplicated Names after Shortening (tailed at:", level, "):\n");
	print(dups_names);

	targets=c();
	for(dups_ix in 1:length(dups_names)){
		targets=c(targets, which(dups_names[dups_ix]==short_names));
	}

	level=level+1;

	if(length(targets)){
		shorten=T;
	}else{
		shorten=F;
	}
}


cat("\n\nFinal Shortened Name Assignments:\n");
print(short_names);

colnames(counts_mat)=short_names;

outmat=counts_mat;

###############################################################################
# Output
cat("\nWriting New Matrix...\n");
fc=file(paste(OutputFileName), "w");

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
