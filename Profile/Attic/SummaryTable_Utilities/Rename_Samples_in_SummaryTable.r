#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"mapping", "m", 1, "character",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary_table.tsv>\n",
	"	-m <rename mapping file, tab-separated>\n",
	"	[-o <output summary table file name, e.g. renamed.summary_table.tsv>\n",
	"\n",	
	"This script will rename your sample names in (-i) based on your mapping file (-m).\n",
	"\n",
	"Mapping File format:\n",
	"	<current sample id>\\t<new sample id>\\n\n",
	"\n");

if(!length(opt$input_file) || !length(opt$mapping)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
MappingFileName=opt$mapping;

if(!length(opt$output_file)){
	output_root=gsub("\\.summary_table\\.tsv", "", opt$input_file);
	OutputFileName = paste(output_root, ".renamed.summary_table.tsv", sep="");
}else{
	OutputFileName=opt$output_file;
}

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Mapping File Name: ", MappingFileName, "\n");       
cat("Output File Name: ", OutputFileName, "\n");       
cat("\n");

###############################################################################
###############################################################################

# Load data
cat("Loading Matrix...\n");
inmat=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="*", row.names=1))
#cat("Original Matrix:\n")
#print(inmat);

# Grab columns we need into a vector, ignore totals, we won't trust it.
counts_mat=inmat[,2:(ncol(inmat))];
#print(counts_mat);

num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");

sample_names=rownames(counts_mat);
category_names=colnames(counts_mat);

###############################################################################
# Load map
cat("Loading Map...\n");
inmap=as.matrix(read.table(MappingFileName, sep="\t", header=FALSE));

# Transfer matrix into map
map=list();
for(i in 1:nrow(inmap)){
	map[[inmap[i,1]]]=inmap[i,2];
}
#print(map);

###############################################################################
# Rename

new_sample_names=numeric();
for(i in 1:num_samples){
	old_samp_name=sample_names[i];
	new_samp_name=map[[old_samp_name]];
	if(length(new_samp_name)==0){
		new_sample_names[i]=old_samp_name;
	}else{
		new_sample_names[i]=new_samp_name;
	}
}

#print(new_sample_names);


###############################################################################
# Output
cat("Writing New Matrix...\n");
fc=file(OutputFileName, "w");

write(paste("sample_id", "total", paste(category_names, collapse="\t"), sep="\t"), file=fc);
for(samp_idx in 1:num_samples){
	total=sum(counts_mat[samp_idx,]);
	outline=paste(new_sample_names[samp_idx],total,paste(counts_mat[samp_idx,], collapse="\t"), sep="\t");
	write(outline, file=fc);
}
close(fc);	

###############################################################################

writeLines("Done.\n")
warns=warnings();
if(!is.null(warns)){
	print(warns);
}

q(status=0)
