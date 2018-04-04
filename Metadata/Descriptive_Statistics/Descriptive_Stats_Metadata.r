#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
library(plotrix);

options(useFancyQuotes=F);
options(digits=5)

params=c(
	"factors", "f", 1, "character",
	"outputroot", "o", 2, "character",
	"exclude", "x", 2, "character",
	"include", "n", 2, "character",
	"response", "y", 2, "character"
);

DEF_MIN_NONNA_PROP=.65;

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <factors>\n",
	"\n",
	"	User specify which subset of variables keep/remove:\n",
	"	[--include <subset variables list filename>]\n",
	"	[--exclude <file with variables to exclude>]\n",
	"\n",
	"This script will read in a metadata file and generate some\n",
	"summary statistics on each column.\n",
	"\n",
	"	1.) Mean, Stdev, 95%CI, Min, Max, Num Samples \n",
	"\n");

if(!length(opt$factors)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputFnameRoot=gsub(".tsv", "", opt$factors);
}else{
	OutputFnameRoot=opt$outputroot;
}
OutputFnameRoot=paste(OutputFnameRoot, ".summary_stats", sep="");

if(!length(opt$include)){
	VariableIncludeListFname="";
}else{
	VariableIncludeListFname=opt$include;
}

if(!length(opt$exclude)){
	VariableExcludeListFname="";
}else{
	VariableExcludeListFname=opt$exclude;
}

if(!length(opt$min_nonna_pro)){
	MinNonNAProp=DEF_MIN_NONNA_PROP;
}else{
	MinNonNAProp=opt$min_nonna_pro;
}

FactorsFname=opt$factors;

cat("Factors Filename/Metadata: ", FactorsFname, "\n", sep="");
cat("Output Filename Root: ", OutputFnameRoot, "\n", sep="");
cat("\n");

if(VariableIncludeListFname!=""){
	cat("Using subset of variables from: ", VariableIncludeListFname, " (Inclusion List)\n");
}

if(VariableExcludeListFname!=""){
	cat("Excluding subset of variables from: ", VariableExcludeListFname, " (Exclusion List)\n");
}

###############################################################################

load_factors=function(fname){
	factors=as.data.frame(read.delim(fname,  header=TRUE, row.names=1, check.names=FALSE, sep="\t"));
	return(factors);
}

load_factors_as_text=function(fname){
	factors=read.delim(fname,  header=TRUE, row.names=1, check.names=FALSE, sep="\t", colClasses="character");
	return(factors);
}

##############################################################################

plot_text=function(strings, max_lines=75){

	plot_page=function(strings){

		num_lines=length(strings);

		top=max_lines;

		plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
			xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
			);
		for(i in 1:num_lines){
			#cat(strings[i], "\n", sep="");
			text(0, top-i*6, strings[i], pos=4, cex=1);
		}

	}

	num_lines=length(strings);
	num_pages=ceiling(num_lines / max_lines);
	#cat("Num Pages: ", num_pages, "\n");
	for(page_ix in 1:num_pages){
		start=(page_ix-1)*max_lines+1;
		end=start+max_lines-1;
		end=min(end, num_lines);
		##print(c(start,end));
		plot_page(strings[start:end]);
	}
}

##############################################################################
# Main Program Starts Here!
##############################################################################

pdf(paste(OutputFnameRoot, ".pdf", sep=""), height=11, width=8.5);

par(family="Courier");

# Load factors
cat("Loading Factors...\n");
factors=load_factors(FactorsFname);
factor_names=colnames(factors);
num_factors=ncol(factors);
cat(num_factors, " Factor(s) Loaded:\n", sep="");
print(factor_names);
cat("\n");

# Subset factors
if(VariableIncludeListFname!=""){
	variable_subset=scan(VariableIncludeListFname, what=character());
	cat("Variable Inclusion List:\n");
	print(variable_subset);
	cat("\n");
	shared_variables=intersect(factor_names, variable_subset);
	cat("Identified:\n");
	print(shared_variables);

	factors=factors[,shared_variables];
	factor_names=colnames(factors);
	num_factors=ncol(factors);
}

if(VariableExcludeListFname!=""){
	variable_subset=scan(VariableExcludeListFname, what=character());
	cat("Variable Exclusion List:\n");
	print(variable_subset);
	cat("\n");
	remaining_variables=setdiff(factor_names, variable_subset);
	cat("Remaining:\n");
	print(remaining_variables);

	factors=factors[,remaining_variables];
	factor_names=colnames(factors);
	num_factors=ncol(factors);
}

factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);
num_samples=num_factor_samples;

###############################################################################

var_dep=list();

par(mfrow=c(3, 2));
bottom_pad=6
par(mar=c(bottom_pad,5,3,1));

for(i in 1:num_factors){
	cur_factor=factors[,i];
	categories=unique(unique(cur_factor));
	num_cat= length(categories);
	numNAs=sum(is.na(cur_factor));
	percNA=round(numNAs/num_samples*100, 2);
	
	cat("\n");
	cat("'", factor_names[i], "' has ", num_cat, " unique values, ", numNAs, " (", percNA, "%) NAs\n", sep="");
	cat("\tUnique: ", paste(head(categories, n=10), collapse=", "), sep="");
	if(num_cat>10){
		cat(", ...");
	}
	cat("\tExample: ", paste(sample(cur_factor, 10), collapse=", "), sep="");
	if(num_cat>10){
		cat(", ...");
	}
	cat("\n");

	isFactor=is.factor(cur_factor);
	cat("Is Factor:", isFactor, "\n");

	nonna_ix=!is.na(cur_factor);
	cur_factor=cur_factor[nonna_ix];
	n=length(cur_factor);

	if(n==0){
		plot(0,0, type="n", xlab="", ylab="", xaxt="n", yaxt="n", main=factor_names[i]);
		text(0,0, "No samples with non-NA values.");
		plot_text(c(
			"All data for variable is NA."
		));

	}else if(isFactor){
		t=table(cur_factor);

		cat_labels=names(t);
		max_lablen=max(nchar(cat_labels));

		if(max_lablen>5){
			par(mar=c(max(bottom_pad, max_lablen/1.7),5,1,1));
		}

		barplot(t, xlab="", ylab="Frequency", main=factor_names[i], col="white", las=2);

		if(max_lablen>5){
			par(mar=c(bottom_pad,5,1,1));
		}

		counts=cbind(t, round(t/n, 4));
		colnames(counts)=c("Counts", "Proportions");
		print(counts);
		num_uniq_cat=nrow(counts);

		
		if(nrow(counts)>10){
			table_summary="More than 10 categories...";
		}else{
			table_summary=capture.output(print(counts));
		}

		plot_text(c(
			paste("N = ", n, " (non-NA)", sep=""),
			paste("Num NAs = ", numNAs, " (", percNA, "%)", sep=""), 
			paste("Num Categories = ", num_uniq_cat, ""),
			"",
			table_summary
		));



	}else{

		mean=mean(cur_factor);
		stdev=sd(cur_factor);
		stderr=stdev/sqrt(n);
		crit=qt(1-(.05/2), n-1);
		intvl=crit*stderr;
		ci95=round(c(mean-intvl, mean+intvl), 4);
		min=min(cur_factor);
		max=max(cur_factor);

		hist(cur_factor, ylab="Frequency", xlab="", main=factor_names[i]);

		plot_text(c(
			paste("N = ", n, " (non-NA)", sep=""),
			paste("Num NAs = ", numNAs, " (", percNA, "%)", sep=""), 
			"",
			paste("mean = ", round(mean, 4), sep=""),
			paste("std dev = ", round(stdev, 4), sep=""),
			paste("95% CI (of mean) = (", ci95[1], ", ", ci95[2], ")", sep=""),
			paste("   i.e. +/- ", round(intvl, 4), sep=""),
			paste("min = ", min, sep=""),
			paste("max = ", max, sep="")
			));

		
	}
}

##############################################################################

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
