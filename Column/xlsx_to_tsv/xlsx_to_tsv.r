#!/usr/bin/env Rscript

###############################################################################

library(getopt);
library(xlsx); # Depends on: install.packages('xlsx')

options(useFancyQuotes=F);
options(width=120);

params=c(
	"input", "x", 1, "character",
	"output", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-x <Input .xlsx file>\n",
	"	[-o <Output Filename Root>]\n",
	"\n",
	"This script will read in an Excel spreadsheet and\n",
	"then generate a tab-separated text file for each\n",
	"sheet.\n",
	"\n",
	"Output will not need to be dos2unix'd and any\n",
	"excess rows and columns will be trimmed.\n",
	"\n");

if(!length(opt$input)){
	cat(usage);
	q(status=-1);
}

InputFile=opt$input;
OutputFileRoot=opt$output;

if(length(OutputFileRoot)==0){
	OutputFileRoot=gsub(".xlsx", "", InputFile);
	OutputFileRoot=gsub(".xls", "", OutputFileRoot);
}

cat("Input XLSX Filename: ", InputFile, "\n");
cat("Output Filename Root: ", OutputFileRoot, "\n");
cat("\n");

##############################################################################

output_sheet=function(data, outputfname){
	cat("Found Sheet.\n");
	num_rows=nrow(data);
	num_cols=ncol(data);

	cat("Loaded: Rows: ", num_rows, " Cols: ", num_cols, "\n");
	cat("Column Names:\n");
	print(colnames(data));

	empty_rows=apply(data,1,function(x){all(is.na(x))});
	data=data[!empty_rows,,drop=F];

	empty_cols=apply(data,2,function(x){all(is.na(x))});
	data=data[,!empty_cols,drop=F];

	num_rows=nrow(data);
	num_cols=ncol(data);
	cat("Cleaned: Rows: ", num_rows, " Cols: ", num_cols, "\n");

	cat("Writing: ", outputfname, "\n");
	fh=file(outputfname, "w");
	cat(file=fh, paste(colnames(data), collapse="\t"), "\n", sep="");
	for(i in 1:num_rows){
		cat(file=fh, paste(data[i,], collapse="\t"), "\n", sep="");
	}
	close(fh);

	cat("ok.\n\n");

	return(0);
}

i=1;
keep_reading=T;

while(keep_reading){
	cat("Reading sheet: ", i, "\n");
	keep_reading=tryCatch({
		data=read.xlsx(InputFile, i);
		#gc();
		output_sheet(data, paste(OutputFileRoot, ".", i, ".tsv", sep=""));
		keep_reading=T;
	}, error=function(e){
		cat("\nError Detected...\n");
		cat("Probably no more sheets left to read.\n");
		print(e);
		keep_reading=F;	
	}, warning=function(w){
		cat("Warning Detected...\n");
		print(w);
		keep_reading=F;
	});
	i=i+1;
}

##############################################################################

cat("\nDone.\n");

print(warnings());
q(status=0);
