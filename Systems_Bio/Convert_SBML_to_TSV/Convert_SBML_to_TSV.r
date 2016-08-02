#!/usr/bin/env Rscript

args=commandArgs(TRUE)

library(sybil)
library(sybilSBML);

arg_count=1;
while(!(is.na(args[arg_count]))){
	
	fname=args[arg_count];
	fname_root=gsub(".xml", "", fname);

	cat("Working on: ", fname, "\n");
	cat("  Building: ", fname_root, "\n");

	model=readSBMLmod(fname);
	modelorg2tsv(model, fname_root, "tsv");
	
	arg_count=arg_count+1;
}

cat("Done.\n");
