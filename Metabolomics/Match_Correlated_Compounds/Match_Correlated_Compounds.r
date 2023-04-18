#!/usr/bin/env Rscript

###############################################################################

library('getopt');

options(useFancyQuotes=F);
options(digits=5)

params=c(
	"factors", "f", 1, "character",
	"targets", "t", 1, "character",
	"correlation_cutoff", "c", 2, "numeric",
	"outputroot", "o", 1, "character"
);

DEFAULT_CORRELATION=0.8

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-t <targets reference variable list>\n",
	"	-f <factors/metadata 'database' file name>\n",
	"	-c <Correlation cutoff, default = ", DEFAULT_CORRELATION, "\n",
	"	-o <output filename root>\n",
	"\n",
	"This script will take in the factors/metadata file and\n",
	"a list of variable references.  Ideally, the factors/metadata\n",
	"should already be transformed to be close to normally distributed.\n",
	"\n",
	"The variables in the list will be searched against the variables\n",
	"in the reference file to find other variables that are highly positively\n",
	"or negatively correlated to it.\n",
	"\n",
	"For each target, a list of the top 'hits' will be reported.\n",
	"\n",
	"\n");

if(
	!length(opt$factors) || 
	!length(opt$targets) ||
	!length(opt$outputroot)
){
	cat(usage);
	q(status=-1);
}

FactorsFname=opt$factors;
OutputFnameRoot=opt$outputroot;
TargetList=opt$targets;
CorrelationCutoff=DEFAULT_CORRELATION;

if(length(opt$correlation_cutoff)){
	CorrelationCutoff=opt$correlation_cutoff;
}

param_text=capture.output({
	cat("\n");
	cat("Factor File Name: ", FactorsFname, "\n");
	cat("Output File Name Root: ", OutputFnameRoot, "\n");
	cat("Target Variable List File Name: ", TargetList, "\n");
	cat("Correlation Cutoff: ", CorrelationCutoff, "\n");
	cat("\n");
});

print(param_text, quote=F);

###############################################################################

load_factors=function(fname){
	factors=as.data.frame(read.delim(fname,  header=TRUE, row.names=1, 
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

pdf(paste(OutputFnameRoot, ".correlations.pdf", sep=""), height=11, width=8.5);

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

num_targets=length(targets);

##############################################################################

targ_values=loaded_factors[,targets];	

num_var=ncol(loaded_factors);
num_samples=nrow(loaded_factors);
cat("Num Variables: ", num_var, "\n");
cat("Num Samples: ", num_samples, "\n");
cat("Num Targets: ", num_targets, "\n");
cat("\n");
cat("Calculating correlation matrix: \n");
cor_mat=cor(targ_values, loaded_factors);

print(dim(cor_mat));

outfn=paste(OutputFnameRoot, ".correlations.tsv", sep="");
fh=file(outfn, "w");

cat(file=fh, "Target\tCorrelatedName\tCorrelation\n");

par(mfrow=c(3,2));
for(i in 1:num_targets){

	target_name=targets[i];
	cat("Looking at: ", target_name, "\n");

	values=cor_mat[i,,drop=F];
	names(values)=colnames(cor_mat);

	hist(values, breaks=seq(-1,1,length.out=30), main=target_name);

	high_cor_ix=(values>=CorrelationCutoff) | (values<=-CorrelationCutoff);
	high_cor_val=values[high_cor_ix];	

	#print(high_cor_ix);
	#print(high_cor_val);
	sorted_high_cor=sort(high_cor_val, decreasing=T);
	num_high_cor=length(sorted_high_cor);
	sorted_high_cor_names=names(sorted_high_cor);
	print(sorted_high_cor);

	hist(sorted_high_cor, breaks=seq(-1,1,length.out=40),
		main=paste("Histogram of | Extreme Correlations | > ", CorrelationCutoff, sep="")
		);
	title(main=paste("Total Extremes: ", num_high_cor), line=0, cex.main=1);

	cat(file=fh, target_name, ":\n", sep="");

	for(i in 1:num_high_cor){
		cat(file=fh, "\t", sorted_high_cor_names[i], "\t", sorted_high_cor[i], "\n", sep="");
	}

	cat("\n\n");	
}

close(fh);


##############################################################################

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
