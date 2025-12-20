#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);
options(width=120);

params=c(
	"input_file", "i", 1, "character",
	"mapping_file", "m", 1, "character",
	"output_file", "o", 1, "character",
	"concat_orig", "c", 2, "logical",
	"make_rfriendly", "r", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input tsv file>\n",
	"	-m <mapping file>\n",
	"	-o <output file root>\n",
	"	[-c (concatenate/prefix w/ original ID)]\n",
	"	[-r (make final name R friendly with make.names)]\n",
	"\n",
	"This script will read in the TSV file and then rename the column names\n",
	"based on the specified map.\n",
	"\n");

if(
	!length(opt$input_file)||
	!length(opt$mapping_file)||
	!length(opt$output_file)
	){
	cat(usage);
	q(status=-1);
}

InputFile=opt$input_file;
OutputFile=opt$output_file;
MappingFile=opt$mapping_file;
MakeRFriendly=opt$make_rfriendly;

ConcatOrig=F;
if(length(opt$concat_orig)){
	ConcatOrig=T;
}

MakeRFriendly=F;
if(length(opt$make_rfriendly)){
	MakeRFriendly=T;
}


cat("Input Filename: ", InputFile, "\n");
cat("Mapping Filename: ", MappingFile, "\n");
cat("Output Filename: ", OutputFile, "\n");
cat("\n");
cat("Concat Original ID: ", ConcatOrig, "\n");
cat("Make R Friendly: ", MakeRFriendly, "\n");
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

#------------------------------------------------------------------------------

load_mapping=function(fname){
	tab=read.table(fname, sep="\t", header=T, comment.char="#", quote="", row.names=NULL);
	mapping=as.list(setNames(as.vector(tab[,2]), tab[,1]));
	return(mapping);
}

remap=function(inval, map){
	num_names=length(inval);
	mapped=rep("", num_names);
	for(i in 1:num_names){
		mapping=map[[inval[i]]];
		if(is.null(mapping)){
			mapped[i]=inval[i];
		}else{
			if(ConcatOrig){
				mapped[i]=paste(cnames[i], ".", mapping, sep="");
			}else{
				mapped[i]=mapping;
			}

		}

		if(MakeRFriendly){
			mapped[i]=make.names(mapped[i]);
		}
	}

	return(mapped);
}

##############################################################################

data=load_factors(InputFile);
mapping=load_mapping(MappingFile);

cnames=colnames(data);
remapped_names=remap(cnames, mapping);
colnames(data)=remapped_names;

write_factors(table=data, fname=OutputFile);

##############################################################################

cat("\nDone.\n");

print(warnings());
q(status=0);
