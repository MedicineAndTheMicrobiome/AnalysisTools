#!/usr/bin/env Rscript

###############################################################################

library('getopt');

options(useFancyQuotes=F);

params=c(
	"paired_map_file", "p", 1, "character",
	"metadata_output", "o", 1, "character",
	"category_name", "c", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-p <paired map file>\n",
	"	-o <output metdatafile>\n",
	"	[-c <category name, default=Category>]\n",
	"\n",
	"This script will read in a paired map file and regenerate\n",
	"a metadata/factor file.\n",
	"\n", sep="");

if(
	!length(opt$paired_map_file) || 
	!length(opt$metadata_output)
){
	cat(usage);
	q(status=-1);
}

PairedMapFile=opt$paired_map_file;
OutputMetadataFile=opt$metadata_output;

if(length(opt$category_name)){
	CategoryName=opt$category_name;
}else{
	CategoryName="Category";
}

cat("\n");
cat("Paired Map File: ", PairedMapFile, "\n", sep="");
cat("Output Metadata File: ", OutputMetadataFile, "\n", sep="");
cat("Category Name: ", CategoryName, "\n", sep="");
cat("\n");

##############################################################################

load_paired=function(fname){
	cat("Loading Paired Map...\n");
	table=data.frame(read.table(fname,  sep="\t", header=TRUE, 
		row.names=1, check.names=FALSE, comment.char=""));
	return(table);
}

##############################################################################

# Load factors
paired_map=load_paired(PairedMapFile);
subject_ids=rownames(paired_map);
pair_category=colnames(paired_map);

num_subjects=length(subject_ids);

cat("Num Subjects:", num_subjects, "\n", sep="");
cat("Pairing Category: \n");
print(pair_category);

a_samp_ids=as.character(paired_map[,1]);
b_samp_ids=as.character(paired_map[,2]);

a_subj_ids=subject_ids;
b_subj_ids=subject_ids;

names(a_subj_ids)=a_samp_ids;
names(b_subj_ids)=b_samp_ids;

metadata_out=matrix(NA, nrow=2*num_subjects, ncol=3);
colnames(metadata_out)=c(CategoryName, "SubjectID", "SampleID");
rownames(metadata_out)=c(a_samp_ids, b_samp_ids);
metadata_out[,CategoryName]=c(rep(pair_category[1], num_subjects), rep(pair_category[2], num_subjects));
metadata_out[,"SubjectID"]=c(a_subj_ids, b_subj_ids);
metadata_out[,"SampleID"]=c(a_samp_ids, b_samp_ids);

cat("Removing NAs...\n");
nrows_orig=2*num_subjects;
not_na_ix=!is.na(metadata_out[,"SampleID"]);
metadata_out=metadata_out[not_na_ix, c("SampleID", CategoryName, "SubjectID")];
nrows_after_na_remove=nrow(metadata_out);
num_removed=nrows_orig-nrows_after_na_remove;
cat("Number of Rows Removed: ", num_removed, "\n");

cat("Writing New Factor File into: ", OutputMetadataFile, "\n");
fname=paste(OutputMetadataFile);
write.table(metadata_out, file=fname, row.names=F, append=F, quote=F, sep="\t");

###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);