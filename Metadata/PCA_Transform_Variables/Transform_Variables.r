#!/usr/bin/env Rscript

###############################################################################

library('getopt');

options(useFancyQuotes=F);
options(digits=5)

params=c(
	"factors", "f", 1, "character",
	"outputroot", "o", 1, "character",
	"targets", "t", 1, "character",
	"subject_id_colname", "S", 2, "character",
	"sample_id_colname", "s", 2, "character",
	"export_original", "O", 2, "logical"
);

NORM_PVAL_CUTOFF=.2;

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

script_path=paste(head(strsplit(script_name, "/")[[1]], -1), collapse="/");
source(paste(script_path, "/Test_and_Apply_Transforms.r", sep=""));

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <factors/metadata file name>\n",
	"	-o <output filename root>\n",
	"	-t <targets variable list>\n",
	"\n",
	"	[-S <Subject ID, column name to export>]\n",
	"	[-s <Sample ID, column name to export>]\n",
	"	[-O (export Original predictor variable list)\n",
	"\n",
	"This script will take in the factors/metadata file and\n",
	"a list of the targets and test and apply the log, sqrt and -1/x\n",
	"transform.\n",
	"\n",
	"The subject,sample ID, and original targeted columns may also be exported\n",	
	"if requested.\n",
	"\n");

if(
	!length(opt$factors) || 
	!length(opt$outputroot) || 
	!length(opt$targets)
){
	cat(usage);
	q(status=-1);
}

FactorsFname=opt$factors;
OutputFnameRoot=opt$outputroot;
TargetList=opt$targets;

ExportOrig=F;
ExportSubjectID="";
ExportSampleID="";

if(length(opt$export_orig)){
	ExportOrig=T;
}
if(length(opt$subject_id_colname)){
	ExportSubjectID=opt$subject_id_colname;
}
if(length(opt$sample_id_colname)){
	ExportSampleID=opt$sample_id_colname;
}

param_text=capture.output({
	cat("\n");
	cat("Factor File Name: ", FactorsFname, "\n");
	cat("Output File Name Root: ", OutputFnameRoot, "\n");
	cat("Target Variable List File Name: ", TargetList, "\n");
	cat("\n");
	cat("Export Sample ID: ", ExportSampleID, "\n");
	cat("Export Subject ID: ", ExportSubjectID, "\n");
	cat("Export Original: ", ExportOrig, "\n");
	cat("\n");
});

print(param_text, quote=F);

###############################################################################

load_factors=function(fname){
	factors=as.data.frame(read.delim(fname,  header=TRUE, row.names=NULL, 
		stringsAsFactors=F,
		check.names=FALSE, sep="\t"));
	return(factors);
}

load_list=function(fname){
	cat("Loading: ", fname, "\n");
	lst=read.delim(fname, header=F, check.names=F, comment.char="#", as.is=T);
	return(lst[,1]);	
}

plot_text=function(strings){

	orig.par=par(no.readonly=T);

	par(mfrow=c(1,1));
	par(family="Courier");
	par(oma=rep(.5,4));
	par(mar=rep(0,4));

	num_lines=length(strings);

	top=max(as.integer(num_lines), 52);

	plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
		xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
		);

	text_size=max(.01, min(.8, .8 - .003*(num_lines-52)));
	#print(text_size);

	for(i in 1:num_lines){
		#cat(strings[i], "\n", sep="");
		strings[i]=gsub("\t", "", strings[i]);
		text(0, top-i, strings[i], pos=4, cex=text_size);
	}

	par(orig.par);
}

##############################################################################
# Main Program Starts Here!
##############################################################################

pdf(paste(OutputFnameRoot, ".transformations.pdf", sep=""), height=11, width=8.5);

plot_text(param_text);

# Load factors
cat("Loading Factors...\n");
loaded_factors=load_factors(FactorsFname);
loaded_factor_names=colnames(loaded_factors);
loaded_sample_names=rownames(loaded_factors);

cat("Loaded factors:\n");
print(loaded_factor_names);
cat("\n");
cat("Loaded sample ids:\n");
print(loaded_sample_names);
cat("\n");

cat("Loading Targets...\n");
targets=load_list(TargetList);
print(targets);

cat("\n");
missing_targets=setdiff(targets, loaded_factor_names);
if(length(missing_targets)){
	cat("WARNING: Missing Targets:\n");
	print(missing_targets);
}else{
	cat("OK:  All Targets Found.\n");
}
target_mat=loaded_factors[,targets];

cat("Testing Targeted Variables for normality.\n");
transformed_mat=test_and_apply_normalizing_transforms(target_mat, NORM_PVAL_CUTOFF);

##############################################################################

out_mat=numeric();
out_colnames=character();

if(ExportSampleID!=""){
	cat("Appending :", ExportSampleID, "\n");
	out_mat=cbind(out_mat, loaded_factors[,ExportSampleID]);
	out_colnames=c(out_colnames, ExportSampleID);
}

if(ExportSubjectID!=""){
	cat("Appending :", ExportSubjectID, "\n");
	out_mat=cbind(out_mat, loaded_factors[,ExportSubjectID]);
	out_colnames=c(out_colnames, ExportSubjectID);
}

if(ExportOrig){
	cat("Appending Original Targeted Values.\n");
	out_mat=cbind(out_mat, target_mat);
	out_colnames=c(out_colnames, colnames(target_mat));
}

cat("Appending Transformed Values.\n");
out_mat=cbind(out_mat, transformed_mat);
out_colnames=c(out_colnames, colnames(transformed_mat));

cat("Writing output file...\n");
colnames(out_mat)=out_colnames;

fname=paste(OutputFnameRoot, ".transformed.tsv", sep="");
write.table(out_mat, file=fname, col.names=T, row.names=F, sep="\t");

##############################################################################

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
