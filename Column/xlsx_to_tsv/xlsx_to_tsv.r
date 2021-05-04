#!/usr/bin/env Rscript

###############################################################################

library(getopt);
library(xlsx); # Depends on: install.packages('xlsx')

options(useFancyQuotes=F);
options(width=120);

params=c(
	"input", "x", 1, "character",
	"sheetnum", "s", 2, "character",
	"output", "o", 2, "character",
	"transpose", "t", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-x <Input .xlsx file>\n",
	"	[-s <sheet numbers (counting from 1) to extract>, eg. 1,4>]\n",
	"	[-o <Output Filename Root>]\n",
	"	[-t <transpose>]\n",
	"\n",
	"This script will read in an Excel spreadsheet and\n",
	"then generate a tab-separated text file for each\n",
	"sheet.\n",
	"\n",
	"Output will not need to be dos2unix'd and any\n",
	"excess rows and columns will be trimmed.\n",
	"\n",
	"Note that if a value is empty, it will be filled in with NAs.\n",
	"\n",
	"Make sure that all columns have a column name, even if it is a sample ID column.\n",
	"The transpose option will automatically transpose the matrix before export.\n",
	"\n");

if(!length(opt$input)){
	cat(usage);
	q(status=-1);
}

InputFile=opt$input;
OutputFileRoot=opt$output;
SheetNumList=opt$sheetnum;
Transpose=length(opt$transpose)>0;

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

cat("Transpose?: ", Transpose, "\n");

##############################################################################

clean_variables=function(in_char_arr){
	ar_len=length(in_char_arr);
	for(i in 1:ar_len){
		tmp=gsub("^\\s+", "", in_char_arr[i]);
		tmp=gsub("\\s+$", "", tmp);
		in_char_arr[i]=tmp;
	}
	return(in_char_arr);
}

##############################################################################

output_sheet=function(data, outputfname, transpose=F){

	cat("Found Sheet.\n");



	num_rows=nrow(data);
	num_cols=ncol(data);

	cat("Loaded: Rows: ", num_rows, " Cols: ", num_cols, "\n");
	cat("Column Names:\n");
	
	cnames=colnames(data);
	cnames=clean_variables(cnames);
	colnames(data)=cnames;

	tab_name=colnames(data)[1];

	empty_rows=apply(data,1,function(x){all(is.na(x))});
	data=data[!empty_rows,,drop=F];

	empty_cols=apply(data,2,function(x){all(is.na(x))});
	data=data[,!empty_cols,drop=F];

	num_rows=nrow(data);
	num_cols=ncol(data);
	cat("Cleaned: Rows: ", num_rows, " Cols: ", num_cols, "\n");

	if(transpose){

		cat("----------------------------------------------\n");
		cat("Before Transpose:\n");
		print(data);

		sample_ids=colnames(data);
		colnames(data)=c();
		if(nrow(data)>1){
			varnames=as.character(data[,1]);
			data=data[,-1, drop=F];
			sample_ids=sample_ids[-1];
		}else{
			varnames=tab_name;	
		}

		cat("\n");
		cat("Identified Variable Names:\n");
		print(varnames);
		cat("Identified Sample IDs:\n");
		print(sample_ids);
		cat("\n");

		cat("Transposing...\n");
		data=t(data);
		data=cbind(sample_ids, data);
		colnames(data)=c("sample_ids", varnames);

		num_rows=nrow(data);
		num_cols=ncol(data);

		cat("After Transpose:\n");
		print(data);
		cat("----------------------------------------------\n");
	}

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
		output_sheet(data, paste(OutputFileRoot, ".", i, ".tsv", sep=""), transpose=Transpose);
	}

}else{

	cat("Trying to extract all sheets...\n");

	i=1;
	keep_reading=T;
	while(keep_reading){
		
		keep_reading=tryCatch({
			cat("Reading sheet: ", i, "\n");
			data=read.xlsx(InputFile, i, check.names=F);
			output_sheet(data, paste(OutputFileRoot, ".", i, ".tsv", sep=""), transpose=Transpose);
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
