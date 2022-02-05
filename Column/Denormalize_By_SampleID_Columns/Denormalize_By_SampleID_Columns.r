#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);
options(width=120);

params=c(
	"input", "i", 1, "character",
	"sampid_cname_list", "s", 1, "character",
	"output", "o", 1, "character",
	"new_sampid_cname", "c", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];


DEFAULT_SAMP_ID_COLNAME="SampleID";

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input tab-separated column data file>\n",
	"	-s <comma separated list of column names>\n",
	"	-o <output tab_separated column data file>\n",
	"	[-c <new sample id column name, default=", DEFAULT_SAMP_ID_COLNAME, "]\n",
	"\n",
	"This script will read in a metadata file\n",
	"and replicate rows for each sample ID in the specified columns.\n",
	"\n",
	"This will allow multiple sample ID to be attributed to the same row.\n",
	"\n",
	"For example, a row will be replicated:\n",
	"\n",
	"   Input:\n",
	"      Age     Sex        OralSampID       StoolSampID\n",
	"       50       M     1234.12345.OW     1234.12345.ST\n",
	"\n",
	"  Output:\n",
	"      SampleID     SampleIDType    Age     Sex        OralSampID       StoolSampID\n",
	" 1234.12345.OW       OralSampID     50       M     1234.12345.OW     1234.12345.ST\n",
	" 1234.12345.ST      StoolSampID     50       M     1234.12345.OW     1234.12345.ST\n",
	"\n");

if(!length(opt$input) || !length(opt$sampid_cname_list) || !length(opt$output)){
	cat(usage);
	q(status=-1);
}

InputFName=opt$input;
SampidColnameList=opt$sampid_cname_list;
OutputFName=opt$output;

NewSampleIDColName=DEFAULT_SAMP_ID_COLNAME;
if(length(opt$new_sampid_cname)){
	NewSampleIDColName=opt$new_sampid_cname;
}

cat("Input Filename: ", InputFName, "\n");
cat("Sample ID Column List: ", SampidColnameList, "\n");
cat("Output Filename: ", OutputFName, "\n");
cat("New Sample ID Column Name: ", NewSampleIDColName, "\n");

targeted_columns=strsplit(SampidColnameList, ",")[[1]];
cat("\nTargeted ID Columns:\n");
print(targeted_columns);
cat("\n");

if(any(NewSampleIDColName==targeted_columns)){
	cat("Error: New Sample ID column name is the same as one of the targets.\n");
	quit(-1);
}

##############################################################################

load_factors=function(fname){
	factors=data.frame(read.delim(fname,  header=TRUE, check.names=FALSE, 
		row.names=NULL,
		as.is=T, comment.char="", quote="", sep="\t"));

	#print(factors);

	dimen=dim(factors);
	cat("Rows Loaded: ", dimen[1], "\n");
	cat("Cols Loaded: ", dimen[2], "\n");

	return(factors);
}

write_factors=function(fname, table){

	dimen=dim(table);
	cat("Rows Exporting: ", dimen[1], "\n");
	cat("Cols Exporting: ", dimen[2], "\n");
	
	write.table(table, fname, quote=F, row.names=F, sep="\t");

}

##############################################################################

# Load factors
factors=load_factors(InputFName);
num_input_rows=nrow(factors);

cat("\nNumber of Input Rows:", num_input_rows, "\n\n");

input_colnames=colnames(factors);

intersect_cnames=intersect(input_colnames, targeted_columns);
if(length(intersect_cnames)!=length(targeted_columns)){
	cat("Error, could not find all targeted columns.\n");
}


outmatrix=character();
outsampid=character();
outsampidtype=character();

for(inrow in 1:num_input_rows){

	for(tcol in targeted_columns){
	
		samp_id=factors[inrow, tcol];
		if(!is.na(samp_id)){
			outmatrix=rbind(outmatrix, factors[inrow,]);
			outsampid=c(outsampid, samp_id);
			outsampidtype=c(outsampidtype, tcol);
		}
	}	
}

outjoined=cbind(outsampid, outsampidtype, outmatrix);
cnames=colnames(outjoined);
cnames[1]=NewSampleIDColName;
cnames[2]=paste(NewSampleIDColName, "Type", sep="");
colnames(outjoined)=cnames;

#print(outjoined);

##############################################################################

write_factors(OutputFName, outjoined);

##############################################################################

cat("\nDone.\n");

print(warnings());
q(status=0);
