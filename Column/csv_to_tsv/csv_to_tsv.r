#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);
options(width=120);

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input comma-separated file (.csv)>\n",
	"	[-o <output tab-separated files (.tsv)>]\n",
	"\n",
	"This script will convert the CSV file into a TSV File.\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFile=opt$input_file;
OutputFile=opt$output_file;

if(length(OutputFile)==0){
	OutputFile=gsub("\\.csv$", "\\.tsv", InputFile);
}

if(length(grep("\\.tsv", OutputFile))==0){
	OutputFile=paste(OutputFile, ".tsv", sep="");
}

cat("Input Filename: ", InputFile, "\n");
cat("Output Filename: ", OutputFile, "\n");
cat("\n");

##############################################################################

load_factors=function(fname){
	factors=read.table(fname, header=TRUE, check.names=FALSE, as.is=T, comment.char="", quote="\"", sep=",");

	dimen=dim(factors);
	cat("Rows Loaded: ", dimen[1], "\n");
	cat("Cols Loaded: ", dimen[2], "\n");

	header=colnames(factors);
	if(header[1]==""){
		header[1]="SampleID";
	}
	colnames(factors)=header;

	return(factors);
}

write_factors=function(fname, table){

	dimen=dim(table);
	cat("Rows Exporting: ", dimen[1], "\n");
	cat("Cols Exporting: ", dimen[2], "\n");
	
	write.table(table, fname, quote=F, row.names=F, sep="\t");

}

##############################################################################

data=load_factors(InputFile);
write_factors(OutputFile, data)

##############################################################################

cat("\nDone.\n");

print(warnings());
q(status=0);
