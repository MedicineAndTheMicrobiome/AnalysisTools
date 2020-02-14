#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "s", 1, "character",
	"keep_regex", "k", 2, "character",
	"remove_regex", "r", 2, "character",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-s <input summary_table.tsv>\n",
	"	[-k \"<keep regular expression\"]\n",
	"	[-r \"<remove regular expression\"]\n",
	"	[-o <output summary_table file name>]\n",
	"\n",	
	"This script will extract the samples from the summary table\n",
	"that are specified in the sample list.\n",
	"\n",
	"Use the -k or -r, to keep or remove sample with that match\n",
	"\n",
	"\n");

if(!length(opt$input_file) || 
	(!length(opt$keep_regex) && !length(opt$remove_regex))){
	cat(usage);
	q(status=-1);
}

###############################################################################

InputFileName=opt$input_file;
OutputFileName=opt$output_file;


OutputFileName=gsub("\\.summary_table\\.tsv", "", OutputFileName);
OutputFileName=gsub("\\.summary_table\\.xls", "", OutputFileName);

OutputFileName=paste(OutputFileName, ".summary_table.tsv", sep="");

regex=character();
keep=logical();

if(length(opt$keep_regex)){
	keep=T;
	regex=opt$keep_regex;	
}
if(length(opt$remove_regex)){
	keep=F;
	regex=opt$remove_regex;	
}

cat("\n");
cat("Input File Name: ", InputFileName, "\n");
cat("Regular Expression: ", regex, "\n");
cat("Keep Matches? : ", keep, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("\n");

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

samp_ids=rownames(counts_mat);

hits=grep(regex, samp_ids);

cat("\nRegex hits:\n");
print(samp_ids[hits]);
num_hits=length(hits);

if(keep){
	cat("Keeping ", num_hits, " samples.\n");

	if(num_hits==0){
		outmat=counts_mat[0,,drop=F];
	}else{
		outmat=counts_mat[hits,,drop=F];
	}
}else{
	cat("Removing ", num_hits, " samples.\n");

	if(num_hits==0){
		outmat=counts_mat;
	}else{
		outmat=counts_mat[-hits,,drop=F];
	}
}

###############################################################################

write_summary_file=function(out_mat, fname){
        fc=file(fname, "w");
        cat(file=fc, paste("sample_id\ttotal", paste(colnames(out_mat), collapse="\t"), sep="\t"));
        cat(file=fc, "\n");

        sample_names=rownames(out_mat);
        num_samples=nrow(out_mat);

	if(num_samples>0){
		for(samp_idx in 1:num_samples){
			total=sum(out_mat[samp_idx,]);
			outline=paste(sample_names[samp_idx], total,
				paste(out_mat[samp_idx,], collapse="\t"), sep="\t");
			cat(file=fc, outline);
			cat(file=fc, "\n");
		}
	}
        close(fc);
}

# Output new summary table
write_summary_file(outmat, OutputFileName);

###############################################################################

writeLines("Done.\n")
print(warnings());

q(status=0)
