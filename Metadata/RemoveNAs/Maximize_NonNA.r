#!/usr/bin/env Rscript

###############################################################################

library('getopt');
source("~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r");

params=c(
        "input", "i", 1, "character",
        "output", "o", 1, "character",
	"num_iterations", "n", 2, "numeric",
	"required_variables", "r", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

NUM_ITERATIONS=50000;

usage = paste(
        "\nUsage:\n", script_name, "\n",
        "       -i <factor file>\n",
        "       -o <output factor filename root>\n",
	"	[-n <number of iterations, default=", NUM_ITERATIONS, ">\n",
	"	[-r <required variables list>\n",
        "\n",
	"This script will read in a factor file, and then\n",
	"remove samples or factors to maximize the number of\n",
	"non-NAs values remaining.\n",
        "\n");

if(!length(opt$input) || !length(opt$output)){
        cat(usage);
        q(status=-1);
}

InputFname=opt$input;
OutputFname=opt$output;

NumIterations=NUM_ITERATIONS;
if(length(opt$num_iterations)){
	NumIterations=opt$num_iterations;	
}

RequiredVariablesFilename="";
if(length(opt$required_variables)){
	RequiredVariablesFilename=opt$required_variables;	
}

cat("Input Filename: ", InputFname, "\n");
cat("Output Filename: ", OutputFname, "\n");
cat("Num Iterations: ", NumIterations, "\n");
cat("Required Variables Filename: ", RequiredVariablesFilename, "\n");

###############################################################################

load_factor_file=function(fn){
        inmat=read.delim(fn, sep="\t", header=TRUE, row.names=1, check.names=F, 
		comment.char="", quote="");

        # Changes spaces to underscore
        var_names=colnames(inmat);
        var_names=gsub(" ", "_", var_names);
        colnames(inmat)=var_names;

        cat("  Num Factors: ", ncol(inmat), "\n", sep="");
        cat("  Num Samples: ", nrow(inmat), "\n", sep="");
        return(inmat);
}

###############################################################################

factors=load_factor_file(InputFname);

required_arr=NULL;
if(RequiredVariablesFilename!=""){
	required_arr=load_list(RequiredVariablesFilename);
}

if(any(is.na(factors))){
        cat("NAs's found in factors...\n");

        script_path=paste(head(strsplit(script_name, "/")[[1]], -1), collapse="/");
        source(paste(script_path, "/../../Metadata/RemoveNAs/Remove_NAs.r", sep=""));
        factors=remove_sample_or_factors_wNA_parallel(
		factors, 
		required_variables=required_arr,
		num_trials=NumIterations, 
		num_cores=64, 
		outfile=OutputFname);

}else{
	cat("No NA's found in factors.  Exiting with no new output.\n");
}



