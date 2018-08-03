#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);
options(width=120);

params=c(
	"input", "f", 1, "character",
	"condition", "c", 1, "character",
	"output", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <input factor file>\n",
	"	-c <condition, eg \"isFemale==1 & SampleType==\\\"Stool\\\"\">\n",
	"	-o <output factor file with samples meeting condition>\n",
	"\n",
	"This script will subset lines out the factor file based on the\n",
	"R expression specified in the condition option.\n",
	"\n");

if(!length(opt$input) || !length(opt$condition) || !length(opt$output)){
	cat(usage);
	q(status=-1);
}

InputFactorFile=opt$input;
Condition=opt$condition;
OutputFactorFile=opt$output;

cat("Input Factor Filename: ", InputFactorFile, "\n");
cat("Condition: ", Condition, "\n");
cat("Output Factor Filename: ", OutputFactorFile, "\n");

##############################################################################

load_factors=function(fname){
	factors=data.frame(read.table(fname,  header=TRUE, check.names=FALSE, as.is=T, comment.char="", quote="", sep="\t"));

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
factors=load_factors(InputFactorFile);
num_samples=nrow(factors);

keep_ix=eval(parse(text=Condition), envir=factors);

cat("Keep Index:\n");
print(keep_ix);

write_factors(OutputFactorFile, factors[keep_ix,,drop=F]);

##############################################################################

cat("\nDone.\n");

print(warnings());
q(status=0);
