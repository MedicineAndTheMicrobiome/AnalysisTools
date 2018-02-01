#!/usr/bin/env Rscript

###############################################################################

library(getopt);
library(xlsx); # Depends on: install.packages('xlsx')

options(useFancyQuotes=F);
options(width=120);

params=c(
	"input", "x", 1, "character",
	"sheetnum", "s", 2, "character",
	"output", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-x <Input .xlsx file>\n",
	"	[-s <sheet numbers (counting from 1) to extract>, eg. 1,4>]\n",
	"	[-o <Output Filename Root>]\n",
	"\n",
	"This script will read in an Excel spreadsheet and\n",
	"then generate a tab-separated text file for each\n",
	"sheet.\n",
	"\n",
	"Output will not need to be dos2unix'd and any\n",
	"excess rows and columns will be trimmed.\n",
	"\n",
	"Note that if a value is empty, it will be filled in with NAs.\n",
	"\n");

if(!length(opt$input)){
	cat(usage);
	q(status=-1);
}

InputFile=opt$input;
OutputFileRoot=opt$output;
SheetNumList=opt$sheetnum;

if(length(OutputFileRoot)==0){
	OutputFileRoot=gsub(".xlsx", "", InputFile);
	OutputFileRoot=gsub(".xls", "", OutputFileRoot);
}

num_sheets_to_extract=NA;
if(length(SheetNumList)){
	SheetNumArr=as.numeric(strsplit(SheetNumList, ",")[[1]]);
	num_sheets_to_extract=length(SheetNumArr);
}

cat("Input XLSX Filename: ", InputFile, "\n");
cat("Output Filename Root: ", OutputFileRoot, "\n");

if(!is.na(num_sheets_to_extract)){
	cat("Sheets to extract:\n");
	print(SheetNumArr);
}else{
	cat("Extracting all sheets.\n");
}
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

	out_mat=matrix("", nrow=num_rows, ncol=num_cols);;
	for(i in 1:num_rows){
		for(j in 1:num_cols){
			out_mat[i,j]=as.character(data[i,j]);
		}
	}

	cat("Writing: ", outputfname, "\n");
	fh=file(outputfname, "w");
	cat(file=fh, paste(colnames(data), collapse="\t"), "\n", sep="");
	for(i in 1:num_rows){
		cat(file=fh, paste(out_mat[i,], collapse="\t"), "\n", sep="");
	}
	close(fh);

	cat("ok.\n\n");

	return(0);
}

##############################################################################

if(!is.na(num_sheets_to_extract)){

	cat("Trying to extract specified sheets.\n");

	for(i in SheetNumArr){
		data=read.xlsx(InputFile, i, check.names=F);
		output_sheet(data, paste(OutputFileRoot, ".", i, ".tsv", sep=""));
	}

}else{

	cat("Trying to extract all sheets...\n");

	i=1;
	keep_reading=T;
	while(keep_reading){
		
		keep_reading=tryCatch({
			cat("Reading sheet: ", i, "\n");
			data=read.xlsx(InputFile, i, check.names=F);
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
}

##############################################################################

cat("\nDone.\n");

print(warnings());
q(status=0);
