#!/usr/bin/env Rscript

###############################################################################

library('getopt');
source("~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r");
source("~/git/AnalysisTools/Metadata/OutputFileLibrary/OutputFileLibrary.r");

params=c(
	"input_file", "i", 1, "character",
	"minimum", "c", 1, "numeric",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	-c <cutoff>\n",
	"	[-o <output root file name>]\n",
	"\n",
	"This script will normalize each sample, then find abundances that are\n",
	"below the specified cutoff, and set them to 0.\n",
	"\n",
	"The sample is then renormalized to sum to 1.\n",
	"\n",
	"\n");

if(!length(opt$input_file) || !length(opt$minimum)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
MinimumCutoff=opt$minimum;
OutputFileRoot=opt$output_file;

if(length(OutputFileRoot)==0){
	OutputFileRoot=InputFileName;
}

OutputFileRoot=gsub("\\.summary_table\\.tsv$", "", OutputFileRoot);
OutputFileRoot=gsub("\\.summary_table\\.xls$", "", OutputFileRoot);

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name Root: ", OutputFileRoot, "\n");       
cat("Minimum Cutoff: ", MinimumCutoff, "\n");
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

sample_ids=rownames(counts_mat);
category_names=colnames(counts_mat);

cat("\n");
cat("Example Sample IDs:\n");
print(head(sample_ids));
cat("...\n\n");
cat("Example Categories:\n");
print(head(category_names));
cat("...\n\n");

###############################################################################

norm=normalize(counts_mat);

kept_stats=matrix(NA, nrow=num_samples, ncol=2);
rownames(kept_stats)=sample_ids;
colnames(kept_stats)=c("Coverage", "NumCatAbvCutoff");

for(samp_id in sample_ids){

	cat("Working on: ", samp_id, "\n");

	abd=norm[samp_id,];

	keep_ix=(abd>=MinimumCutoff);

	kept_prop=sum(abd[keep_ix]);
	num_kept=sum(keep_ix);

	abd[!keep_ix]=0;

	sum_abd=sum(abd);
	abd=abd/sum_abd;

	norm[samp_id,]=signif(abd,5);
	
	#cat("Cumulative above cutoff: ", kept_prop, "\n");
	
	kept_stats[samp_id,]=c(signif(kept_prop,5), num_kept);

}

###############################################################################

outputfn=paste(OutputFileRoot, ".mincutoff_", MinimumCutoff, ".summary_table.tsv", sep="");
write_summary_file(norm, outputfn);

###############################################################################

outmat=cbind(rownames(kept_stats), kept_stats);
colnames(outmat)=c("SampleID", "PreNormCov", "NumCatAbvCutoff");
outputfn=paste(OutputFileRoot, ".mincutoff_", MinimumCutoff, ".stats.tsv", sep="");
write.table(outmat, outputfn, quote=F, row.names=F, col.names=T);

###############################################################################

writeLines("Done.\n")
warns=warnings();
if(length(warns)>0){
	print(warnings());
}

q(status=0)
