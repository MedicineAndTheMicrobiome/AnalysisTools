#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);
options(width=120);

params=c(
	"input_list_cli", "i", 2, "character",
	"input_list_file", "l", 2, "character",
	"shared_id", "c", 2, "character",
	"output", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	Use one of the input options:\n",
	"	[-i <input, comma separated list of tsv files>]\n",
	"	[-l <input, text file with list of target files>]\n",
	"\n",
	"	-c <shared ID column name.>\n",
	"	-o <output tab_separated column data file>\n",
	"\n",
	"This script will read in the metadata files specified in the comma-separated\n",
	"list of the -i option or file with the -l option, and merge the files based\n",
	"on the shared column names.\n",
	"\n",
	"Only the shared IDs across all the files will be returned.\n",
	"\n");

if(!length(opt$output) || !length(opt$shared_id)){
	cat(usage);
	q(status=-1);
}

InputFNameList=opt$input_list_cli;
InputFNameFile=opt$input_list_file;

Formulas=opt$formulas;
OutputFName=opt$output;
CollapseBySampleID=opt$collapse_samp_id;

use_cli=NA;
if(length(InputFNameList)){
	cat("List specified on commmand line.\n");
	use_cli=T;
}
if(length(InputFNameFile)){
	cat("List specified in file. (", InputFNameFile, ")\n");
	use_cli=F;
}
if(is.na(use_cli)){
	cat("Please specify one of the input options.\n");
}

if(length(opt$shared_id)){
	SharedID=opt$shared_id;
}

cat("\n");
cat("Shared ID: ", SharedID, "\n");
cat("Output Filename: ", OutputFName, "\n");
cat("\n");


target_list=character();
if(use_cli){
	target_list=strsplit(InputFNameList, ",")[[1]];	
}else{
	target_list=read.table(InputFNameFile, as.is=T, sep="\n", header=F, check.names=F,
			comment.char="", quote="")[,1];
}

cat("Target list of files to merge:\n");
print(target_list);

##############################################################################

load_factors=function(fname, prim_key){
	factors=data.frame(read.table(fname,  header=TRUE, 
		row.names=prim_key,
		check.names=FALSE, as.is=T, 
		comment.char="#", quote="", 
		sep="\t"));

	#print(factors);

	dimen=dim(factors);
	cat("Primary Key: ", prim_key, "\n", sep="");
	cat("Rows Loaded: ", dimen[1], "\n");
	cat("Cols Loaded: ", dimen[2], "\n");

	return(factors);
}

write_factors=function(fname, table, prim_key){

	dimen=dim(table);
	cat("Primary Key: ", prim_key, "\n", sep="");
	cat("Rows Exporting: ", dimen[1], "\n");
	cat("Cols Exporting: ", dimen[2], "\n");
	
	out_mat=cbind(rownames(table), table);
	colnames(out_mat)=c(prim_key, colnames(table));
	write.table(out_mat, fname, quote=F, row.names=F, sep="\t");

}

##############################################################################

num_target_files=length(target_list);

cat("Number of Target Files: ", num_target_files, "\n\n", sep="");

factors_list=list();
shared_ids_list=list();

for(fname in target_list){
	cat("Loading: ", fname, "\n", sep="");
	factors_list[[fname]]=load_factors(fname, SharedID);
	shared_ids_list[[fname]]=rownames(factors_list[[fname]]);
}

cat("\n\n");

shared_ids=shared_ids_list[[1]];
for(fname in target_list){
	cat("Intersecting with IDs from: ", fname, "\n");
	shared_ids=intersect(shared_ids, shared_ids_list[[fname]]);
	cat("Num Shared: ", length(shared_ids), "\n\n");	
}

##############################################################################

shared_ids=sort(shared_ids);
combined_factors=matrix(NA, nrow=length(shared_ids), ncol=0);
rownames(combined_factors)=shared_ids;

for(fname in target_list){
	combined_factors=cbind(combined_factors, factors_list[[fname]][shared_ids,]);
}

write_factors(OutputFName, combined_factors, SharedID);

##############################################################################

print(warnings());
q(status=0);
