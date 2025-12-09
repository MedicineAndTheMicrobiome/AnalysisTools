#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);
options(width=120);

params=c(
	"input_file", "i", 1, "character",
	"output_file_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input tsv file>\n",
	"	-o <output file root>\n",
	"\n",
	"This script will read in the input file, apply R's make.names\n",
	"and write the cleaned up results to the output file.\n",
	"\n",
	"The following will be generated:\n",
	"	<output file root>.rfriendly.tsv\n",
	"	<output file root>.rfriendly.map\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFile=opt$input_file;
OutputFileRoot=opt$output_file;

cat("Input Filename: ", InputFile, "\n");
cat("Output Filename Root: ", OutputFileRoot, "\n");
cat("\n");

##############################################################################

load_factors=function(fname){
	factors=read.table(fname, header=TRUE, check.names=FALSE, as.is=T, 
		comment.char="", quote="\"", sep="\t");

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

cnames=colnames(data);
clean_cnames=make.names(cnames, unique=T);
colnames(data)=clean_cnames;

cat("\n");
cat("Writing output with cleaned names.\n");
write_factors(paste(OutputFileRoot, ".rfriendly.tsv", sep=""), data);

map=cbind(cnames, clean_cnames);
colnames(map)=c("Original", "RFriendly"); 

cat("\n");
cat("Writing mapping of original to cleaned names.\n");
write_factors(paste(OutputFileRoot, ".rfriendly.map", sep=""), map);

##############################################################################

cat("\nDone.\n");

print(warnings());
q(status=0);
