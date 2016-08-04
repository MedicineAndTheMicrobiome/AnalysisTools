#!/usr/bin/env Rscript

###############################################################################

FIELD_SEP=";";

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"field_separator", "s", 2, "character",
	"output_file", "o", 2, "character",
	"misc_cleanup", "m", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.tsv>\n",
	"	[-s <field separator, default='", FIELD_SEP, "'>]\n",
	"	[-o <output summary table file name>]\n",
	"	[-m <string to remove for extra cleanup, may be a list with the same -f field separator>]\n",
	"\n",	
	"This script will reduce the taxa names in your summary table.\n",
	"\n",
	"If the categories/taxa names are a list, then this\n",
	"script will split the list according to the field\n",
	"separator and keep the one furthest to the right.\n",
	"\n",
	"For example, if the taxa name is: \n",
	"	k__Fungi;p__Basidiomycota;c__unclassified_Basidiomycota;o__unclassified_Basidiomycota;f__unclassified_Basidiomycota\n",
	"\n",
	"Then the kept taxa name will be: \n",
	"	f__unclassified_Basidiomycota\n",
	"\n",
	"\n",
	"The -m option can specify of list of miscellaneous text to also remove.  For example, \n",
	"	-m \"k__;p__;c__;o__;f__;g__;s__\" will remove the f__ from the above example.\n",
	"\n", sep="");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
MappingFileName=opt$mapping;
MiscCleanup=opt$misc_cleanup;

if(!length(opt$field_separator)){
	FieldSeparator=FIELD_SEP;
}else{
	FieldSeparator=opt$field_separator;
}	

if(!length(opt$misc_cleanup)){
	MiscCleanupArr=c();
}else{
	MiscCleanupArr=strsplit(opt$misc_cleanup, split=FieldSeparator)[[1]];
}

if(!length(opt$output_file)){
	output_root=gsub("\\.summary_table\\.xls", "", opt$input_file);
	output_root=gsub("\\.summary_table\\.tsv", "", output_root);
	OutputFileName = paste(output_root, ".rdc.summary_table.tsv", sep="");
}else{
	OutputFileName=opt$output_file;
}

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Field Separator: ", FieldSeparator, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("Miscellaneous Cleanup String: ", MiscCleanup, "\n");
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
# Rename 

cat("Identifying new category names...\n");

new_category_names=category_names;

# Clean up miscellaneous characters 
if(length(MiscCleanupArr)){
	for(patt in MiscCleanupArr){
		for(i in 1:num_categories){
			category_names[i]=gsub(patt, "", category_names[i]);
		}
	}
}

problematic_ids=rep(TRUE, num_categories);
levels_up=1;

while(any(problematic_ids)){

	# Grab last element in list of names as the new name
	for(i in which(problematic_ids)){
		arr=strsplit(category_names[i], split=FieldSeparator)[[1]];
		new_category_names[i]=paste(tail(arr, levels_up), collapse=FieldSeparator);
		if(levels_up>1){
			cat(category_names[i], "  ->\n\t", new_category_names[i], "\n");
		}
	}

	# Confirm category names are still unique
	tbl=table(new_category_names);
	problem_names=names(tbl[tbl>1]);
	if(length(problem_names)){
		cat("\nProblematic names:\n");
		print(tbl[tbl>2]);
		cat("\n");
	}
	
	# Identify categories to include a higher level
	problematic_ids=rep(FALSE, num_categories);
	for(name in problem_names){
		problematic_ids=problematic_ids | (name==new_category_names);
	}
	levels_up=levels_up+1;
}

category_names=new_category_names;


###############################################################################
# Output
cat("Writing New Matrix...\n");
fc=file(OutputFileName, "w");

write(paste("sample_id", "total", paste(category_names, collapse="\t"), sep="\t"), file=fc);
for(samp_idx in 1:num_samples){
	total=sum(counts_mat[samp_idx,]);
	outline=paste(sample_names[samp_idx],total,paste(counts_mat[samp_idx,], collapse="\t"), sep="\t");
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
