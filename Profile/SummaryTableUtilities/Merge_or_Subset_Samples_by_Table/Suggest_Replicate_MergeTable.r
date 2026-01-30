#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"mapping_table", "m", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table>\n",
	"	-m <sample mapping table>\n",
	"\n",	
	"This script will read in both the summary table\n",
	"and suggest a mapping to merge replicates.\n",
	"\n",
	"The script will read in the sample IDs from the summary table\n",
	"look for a numeric extension, and then generate a map\n",
	"that 'Merge_or_Subset_Samples_by_Table.r' needs to perform\n",
	"it's task.\n",
	"\n",
	"For example:\n",
	"	Original                               Collapsed\n",
	"	0275.LS038.TP2.PRE.20250604.ST         0275.LS038.TP2.PRE.20250604.ST\n",
	"	0275.LS038.V0a.TP1.PRE.20241212.OW     0275.LS038.V0a.TP1.PRE.20241212.OW\n",
	"	0275.LS038.V0a.TP1.PRE.20241212.OW.2   0275.LS038.V0a.TP1.PRE.20241212.OW\n",
	"\n",
	"Note: The first sample ID is left alone, but the 2nd and 3rd will be collapsed together.\n",
	"\n");

if(!length(opt$input_file) || !length(opt$mapping_table)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
MappingTable=opt$mapping_table;

options(width=200);

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Mapping Table Name: ", MappingTable, "\n");

###############################################################################
###############################################################################

# Load data
inmat=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, row.names=1, quote=NULL))
#cat("Original Matrix:\n")
#print(inmat);

# Grab columns we need into a vector, ignore totals, we won't trust it.
counts_mat=inmat[,2:(ncol(inmat)), drop=F];

num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);
cat("\n");
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");
cat("\n");

###############################################################################

orig_sample_ids=rownames(inmat);
#cat("Original Sample IDs:\n");
#print(orig_sample_ids);

merged_sample_ids=gsub("\\.\\d+$", "", orig_sample_ids);

outtab=cbind(orig_sample_ids, merged_sample_ids);
colnames(outtab)=c("Original", "Collapsed");

print(outtab);

write.table(outtab, MappingTable, col.names=T, row.names=F, quote=F, sep="\t");

###############################################################################

cat("\nDone.\n");
if(!is.null(warnings())){
	print(warnings());
}

q(status=0);
