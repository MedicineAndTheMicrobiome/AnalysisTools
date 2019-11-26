#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
library('vegan');

options(useFancyQuotes=F);
options(width=80);

DEF_NUM_TOP_CAT=25;

source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');

params=c(
	"summary_file", "s", 1, "character",
	"factor_file", "f", 1, "character",
	"sample_coln", "S", 1, "character",
	"replicate_coln", "R", 1, "character",
	"outputroot", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table>\n",
	"	-f <factors (replicate map)>\n",
	"	-S <Sample Name column name>\n",
	"	-R <Replicate Name column name>\n",
	"	-o <output filename root>\n",
	"\n",
	"	[-R (pay attention to 'remaining' category)]\n",
	"	[-v number of ALR categories to include, default=", DEF_NUM_TOP_CAT, ">]\n",
	"	[-x <shorten category names, with separator in double quotes (default=\"\">)]\n",
	"\n",
	"This script will perform an analyses across a set of\n",
	"of replicates using bootstrapping to look at the effects\n",
	"of sequencing depth.\n",
	"\n",
	"Only a few columns are used in the factor file:\n",
	"	1.) Sample ID (to look up in the summary table)\n",
	"	2.) Sample Name (Underlying sample being replicated)\n",
	"	3.) Replicate Name (e.g. Run ID)\n",
	"\n");

if(
	!length(opt$summary_file) || 
	!length(opt$factor_file) || 
	!length(opt$sample_coln) || 
	!length(opt$replicate_coln) ||
	!length(opt$outputroot)
){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputRoot=gsub(".summary_table.xls", "", opt$summary_file);
	OutputRoot=gsub(".summary_table.tsv", "", OutputRoot);
}else{
	OutputRoot=opt$outputroot;
}

SummaryFile=opt$summary_file;
FactorsFile=opt$factor_file;
OutputRoot=opt$outputroot;
Replicate_ColumnName=opt$replicate_coln;
Sample_ColumnName=opt$sample_coln;


cat("\n");
cat("Summary File : ", SummaryFile, "\n", sep="");
cat("Factors File: ", FactorsFile, "\n", sep="");
cat("Output File: ", OutputRoot, "\n", sep="");
cat("Replicate ID Column Name: ", Replicate_ColumnName, "\n", sep="");
cat("Sample Name Column Name: ", Sample_ColumnName, "\n", sep="");
cat("\n");

pdf(paste(OutputRoot, ".rep_analys.correl.pdf", sep=""), height=14, width=8.5);

##############################################################################

load_factors=function(fname){
	cat("Loading Factors: ", fname, "\n");
	factors=data.frame(read.table(fname, sep="\t",  header=TRUE, row.names=1, 
		check.names=FALSE, comment.char="", quote=""));
	factor_names=colnames(factors);

	ignore_idx=grep("^IGNORE\\.", factor_names);

	if(length(ignore_idx)!=0){
		return(factors[-ignore_idx]);
	}else{
		return(factors);
	}
}

load_summary_file=function(fname){
	cat("Loading Summary Table: ", fname, "\n");
	inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, 
		comment.char="", row.names=1))
	counts_mat=inmat[,2:(ncol(inmat))];
	return(counts_mat);
}

plot_text=function(strings){
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
}

##############################################################################

plot_text(c(
	paste("Summary File : ", SummaryFile, sep=""),
	paste("Factors File: ", FactorsFile, sep=""),
	paste("Output File: ", OutputRoot, sep=""),
	paste("Replicate ID Column Name: ", Replicate_ColumnName, "\n", sep=""),
	paste("Sample Name Column Name: ", Sample_ColumnName, "\n", sep="")
));

# Load matrix
counts=load_summary_file(SummaryFile);
num_taxa=ncol(counts);
num_samples=nrow(counts);

##############################################################################

# Load factors
factors=load_factors(FactorsFile);
factor_names=colnames(factors);
num_factors=ncol(factors);
factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);

cat("\n");
cat(num_factors, " Factor(s) Loaded:\n", sep="");
#print(factor_names);
cat("\n");

##############################################################################

# Reconcile factors with samples
factor_sample_ids=rownames(factors);
counts_sample_ids=rownames(counts);

shared_sample_ids=intersect(factor_sample_ids, counts_sample_ids);
num_shared_sample_ids=length(shared_sample_ids);
num_factor_sample_ids=length(factor_sample_ids);
num_counts_sample_ids=length(counts_sample_ids);

cat("Num counts sample IDs: ", num_counts_sample_ids, "\n");
cat("Num factor sample IDs: ", num_factor_sample_ids, "\n");
cat("Num shared sample IDs: ", num_shared_sample_ids, "\n");
cat("\n");

cat("Samples missing from count information:\n");
print(setdiff(factor_sample_ids, counts_sample_ids));
cat("\n");
cat("Samples missing from factor information:\n");
print(setdiff(counts_sample_ids, factor_sample_ids));
cat("\n");
cat("Total samples shared: ", num_shared_sample_ids, "\n");

shared_sample_ids=sort(shared_sample_ids);

##############################################################################
# For each sample name

replicate_group=as.character(factors[, Replicate_ColumnName]);
sample_names=factors[, Sample_ColumnName];
sample_ids=rownames(factors);

uniq_samp_names=unique(sample_names);
uniq_repl_names=unique(replicate_group);

num_uniq_samp=length(uniq_samp_names);
num_uniq_repl=length(uniq_repl_names);
num_samp_ids=length(sample_ids);

cat("Number of Sample ID (w/ repl info): ", num_samp_ids, "\n");
cat("Number of Unique Sample Names: ", num_uniq_samp, "\n");
cat("\n");
cat("Number of Unique Replicate Group: ", num_uniq_repl, "\n");
print(uniq_repl_names);

##############################################################################

signf_char=function(x){
	if(x<.001){
		return("***");
	}else if(x<.01){
		return("**");
	}else if(x<.05){
		return("*");
	}else if(x<.1){
		return(".");
	}else{
		return("");
	}
}

##############################################################################

#print(uniq_samp_names);
#print(uniq_repl_names);

rep_samp_matrix=matrix(character(), nrow=num_uniq_samp, ncol=num_uniq_repl);
rownames(rep_samp_matrix)=uniq_samp_names;
colnames(rep_samp_matrix)=uniq_repl_names;

rep_samp_depth_matrix=matrix(numeric(), nrow=num_uniq_samp, ncol=num_uniq_repl);
rownames(rep_samp_depth_matrix)=uniq_samp_names;
colnames(rep_samp_depth_matrix)=uniq_repl_names;

samp_ids=rownames(factor);

for(i in 1:num_samp_ids){
	samp_id=sample_ids[i];
	run_id=as.character(factors[i,Replicate_ColumnName]);
	samp_name=as.character(factors[i,Sample_ColumnName]);
	rep_samp_matrix[samp_name, run_id]=samp_id;
	rep_samp_depth_matrix[samp_name, run_id]=sum(counts[samp_id,]);	
	
}

print(rep_samp_matrix);
print(rep_samp_depth_matrix);

correl_val_matrix=matrix(NA, nrow=num_uniq_repl, ncol=num_uniq_repl);
colnames(correl_val_matrix)=uniq_repl_names;
rownames(correl_val_matrix)=uniq_repl_names;

correl_pval_matrix=matrix(NA, nrow=num_uniq_repl, ncol=num_uniq_repl);
colnames(correl_pval_matrix)=uniq_repl_names;
rownames(correl_pval_matrix)=uniq_repl_names;


log_correl_val_matrix=matrix(NA, nrow=num_uniq_repl, ncol=num_uniq_repl);
colnames(log_correl_val_matrix)=uniq_repl_names;
rownames(log_correl_val_matrix)=uniq_repl_names;

log_correl_pval_matrix=matrix(NA, nrow=num_uniq_repl, ncol=num_uniq_repl);
colnames(log_correl_pval_matrix)=uniq_repl_names;
rownames(log_correl_pval_matrix)=uniq_repl_names;

for(repA in uniq_repl_names){
	for(repB in uniq_repl_names){	
		v1=rep_samp_depth_matrix[,repA];
		v2=rep_samp_depth_matrix[,repB];
	
		cor.res=cor.test(v1, v2);	
		correl_val_matrix[repA, repB]=cor.res$estimate;
		correl_pval_matrix[repA, repB]=cor.res$p.value;

		cor.res=cor.test(log10(v1+1), log10(v2+1));	
		log_correl_val_matrix[repA, repB]=cor.res$estimate;
		log_correl_pval_matrix[repA, repB]=cor.res$p.value;
	}
}

plot_text(c(
	"Depth:",
	"Correlations:",
	capture.output(print(signif(correl_val_matrix, 4), quote=F)),
	"",
	"P-Values:", 
	capture.output(print(signif(correl_pval_matrix, 4), quote=F)),
	"",
	"",
	"Log 10 (Depth):",
	"Correlations:",
	capture.output(print(signif(log_correl_val_matrix, 4), quote=F)),
	"",
	"P-Values:", 
	capture.output(print(signif(log_correl_pval_matrix, 4), quote=F))
));

repcol=rainbow(num_uniq_repl, start=0, end=.8);
names(repcol)=uniq_repl_names;

par(mar=c(4,4,5,1));
par(mfrow=c(4,1));
for(repA in uniq_repl_names){
	v1=rep_samp_depth_matrix[,repA];
	notna_ix=!is.na(v1);

	byA_mat=rep_samp_depth_matrix[notna_ix,];

	v1_ord=order(byA_mat[,repA]);
	byA_mat=byA_mat[v1_ord,];

	numrows=nrow(byA_mat);
	logval=log10(byA_mat+1);
	maxval=max(logval[!is.na(logval)]);
	plot(0, 0, type="n", ylim=c(0, maxval), xlim=c(0, numrows),
		main=paste("Ordered by ", repA, sep=""),
		ylab="Log10(Depth)", xlab="Samples (Ordered)");
	
	for(repB in uniq_repl_names){
		if(repB==repA){
			next;
		}
	
		points(1:numrows, logval[,repB], type="l", col=repcol[repB]);
		points(1:numrows, logval[,repB], type="l", col="grey", lwd=.5);

	}
	points(1:numrows, logval[,repA], type="l", col=repcol[repA], lwd=3);
	points(1:numrows, logval[,repA], type="l", col="black", lwd=.5);
}



##############################################################################

cat("Done.\n");
print(warnings());
q(status=0);
