#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
options(useFancyQuotes=F);

source("~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r");
source("~/git/AnalysisTools/Metadata/OutputFileLibrary/OutputFileLibrary.r");

NUM_VAR=20;

params=c(
	"input_fn", "i", 1, "character",
	"reference_id_fn", "r", 1, "character",
	"output_root", "o", 1, "character",
	"start_col", "b", 2, "numeric",
	"end_col", "e", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input alr/variable data filename>\n",
	"	-r <reference sample/subject IDs filename>\n",
	"	-o <output root>\n",
	"\n",
	"    Column Names in Factor File:\n",
	"	[-b <begin/first column number, default=first>]\n",
	"	[-e <end/last column number, default=last>\n",
	"\n",
	"This script will compute the difference and pvalue\n",
	"between the reference and the other samples\n",
	"\n",
	"The first column is the sample/ID, so the first data\n",
	"column is 2.  The last column is the number of columns\n",
	"in the file.\n",
	"\n");

if(	!length(opt$input_fn) || 
	!length(opt$reference_id_fn) || 
	!length(opt$output_root)){

	cat(usage);
	q(status=-1);
}

StartCol=2;
EndCol=Inf;

if(length(opt$start_col)){
	StartCol=opt$start_col;
}

if(length(opt$end_col)){
	EndCol=opt$end_col;
}

InputDataFilename=opt$input_fn;
ReferenceListFilename=opt$reference_id_fn;
OutputFilenameRoot=opt$output_root;


pdf(paste(OutputFilenameRoot, ".top_diff.pdf", sep=""), height=11, width=8.5);

plot_text(c(
	paste("Data File : ", InputDataFilename, sep=""),
	paste("Reference File: ", ReferenceListFilename, sep=""),
	paste("Output File Root: ", OutputFilenameRoot, sep=""),
	"",
	paste("Start Col: ", StartCol, sep=""),
	paste("End Col: ", EndCol, sep="")
), echo=T);


##############################################################################

load_factors=function(fname){
	factors=data.frame(read.table(fname,  header=TRUE, row.names=1, check.names=FALSE));
	factor_names=colnames(factors);

	ignore_idx=grep("^IGNORE\\.", factor_names);

	if(length(ignore_idx)!=0){
		return(factors[-ignore_idx]);
	}else{
		return(factors);
	}
}

##############################################################################

data=load_factors(InputDataFilename);

num_samples=nrow(data);
num_variables=ncol(data);

cat("Num Samples: ", num_samples, "\n");
cat("Num Variables: ", num_variables, "\n");

EndCol=min(EndCol, num_variables+1);

cat("Sample/Subject IDs:\n");
sample_ids=rownames(data);
print(sample_ids);

cat("\n");
cat("Variable Names in File:\n");
variable_names=colnames(data);
cat("  Beginning:\n");
print(head(variable_names));
cat("...\n  End:\n");
print(tail(variable_names));
cat("\n");

#------------------------------------------------------------------------------
# Adjust the targeted column indices so they line up with what was read in

var_start_ix=StartCol-1;
var_end_ix=EndCol-1;

cat("Actual columns identified:\n");
cat("  Start: ", var_start_ix, "\n");
cat("    End: ", var_end_ix, "\n");

selected_data=data[,var_start_ix:var_end_ix,drop=F];

cat("\n");
cat("Variable Names in File:\n");
variable_names=colnames(selected_data);
cat("  Beginning:\n");
print(head(variable_names));
cat("...\n  End:\n");
print(tail(variable_names));
cat("\n");

#############################################################################	

reference_ids=load_list(ReferenceListFilename, "Reference ID list");
experimental_ids=setdiff(rownames(selected_data), reference_ids);

#############################################################################	

num_var=ncol(selected_data);

reference=selected_data[reference_ids,,drop=F];
experimental=selected_data[experimental_ids,,drop=F];

#------------------------------------------------------------------------------

hdr=c("CombMeanALR", "MedianRef", "MedianExp", "DiffMedian", "WilcoxPval",
	"signf_0.100", "signf_0.050", "signf_0.010", "signf_0.001"	
	);
hdr_len=length(hdr);
diff_tab=matrix(0, nrow=num_var, ncol=hdr_len);
colnames(diff_tab)=hdr;
rownames(diff_tab)=variable_names;
 
for(var_ix in variable_names){
	
	all_val=selected_data[,var_ix];
	ref_val=reference[,var_ix];
	exp_val=experimental[,var_ix];

	median_ref=median(ref_val);
	median_exp=median(exp_val);

	wcx_res=wilcox.test(ref_val, exp_val);
	#print(wcx_res);
	#print(wcx_res$p.value);
	pval=wcx_res$p.value;

	diff_tab[var_ix,"CombMeanALR"]=mean(all_val);
	diff_tab[var_ix,"MedianRef"]=median_ref;
	diff_tab[var_ix,"MedianExp"]=median_exp;

	diff_tab[var_ix,"DiffMedian"]=median_exp-median_ref;
	diff_tab[var_ix,"WilcoxPval"]=pval;
	diff_tab[var_ix,"signf_0.100"]=pval<0.100;
	diff_tab[var_ix,"signf_0.050"]=pval<0.050;
	diff_tab[var_ix,"signf_0.010"]=pval<0.010;
	diff_tab[var_ix,"signf_0.001"]=pval<0.001;

}

print(diff_tab);
sigf_ix=order(diff_tab[,"WilcoxPval"], decreasing=F);

ordered_diff_tab=diff_tab[sigf_ix,];

#############################################################################	
# Generate histograms ordered by decreasing abundance

par(mfrow=c(4,2));
for(var_ix in rownames(ordered_diff_tab)){
	all_val=selected_data[,var_ix];
	ref_val=reference[,var_ix];
	exp_val=experimental[,var_ix];

	median_ref=ordered_diff_tab[var_ix,"MedianRef"];
	median_exp=ordered_diff_tab[var_ix,"MedianExp"];
	pval=ordered_diff_tab[var_ix,"WilcoxPval"];

	all_hist=hist(all_val, breaks=20, plot=F);

	hist(ref_val, main=sprintf("%s: Reference", var_ix), breaks=all_hist$breaks);
	title(main=sprintf("Median: %6.4f ", median_ref), cex.main=.8, font.main=3, line=1);
	abline(v=median_ref, col="black", lwd=2.2);
	abline(v=median_ref, col="grey", lwd=1.6, lty="dashed");

	histcol="grey";
	if(pval<0.1){
		histcol=ifelse(median_exp>median_ref, "red", "blue");
	}

	hist(exp_val, main=sprintf("%s: Experimental", var_ix), 
		breaks=all_hist$breaks,
		col=histcol
		);
	title(main=sprintf("Median: %6.4f", median_exp), cex.main=.8, font.main=3, line=1);
	title(main=sprintf("p-val: %6.4f", pval), cex.main=.8, font.main=3, line=0.5);
	abline(v=median_ref, col="black", lwd=2.2);
	abline(v=median_ref, col="grey", lwd=1.6, lty="dashed");
	abline(v=median_exp, col="black", lwd=2.2);
	abline(v=median_exp, col=histcol, lwd=1.6, lty="dashed");
}

#############################################################################	
# Generate Volcano Plot

par(mfrow=c(1,1));
par(mar=c(5,5,5,5));
plot(diff_tab[,"DiffMedian"], -log10(diff_tab[,"WilcoxPval"]),
	xlab="ALR Difference From Reference (Exp-Ref)", ylab="Wilcoxon -log10(P-value)",
	main="Volcano Plot: Difference vs Significance");
abline(v=0, lwd=1.5, col="black", lty="dashed");
sig_levels=c(0.05,0.1,0.01,0.001);
abline(h=-log10(sig_levels), col="blue", lty="dotted");
axis(side=4, at=-log10(sig_levels), label=sig_levels);

#############################################################################	

sum_signf_0.100=sum(diff_tab[,"signf_0.100"]);
sum_signf_0.050=sum(diff_tab[,"signf_0.050"]);
sum_signf_0.010=sum(diff_tab[,"signf_0.010"]);
sum_signf_0.001=sum(diff_tab[,"signf_0.001"]);

plot_text(c(
	paste("Num Significant:"),
	sprintf("   < 0.100 : %i", sum_signf_0.100), 
	sprintf("   < 0.050 : %i", sum_signf_0.050), 
	sprintf("   < 0.010 : %i", sum_signf_0.010), 
	sprintf("   < 0.001 : %i", sum_signf_0.001),
	""
), echo=T);

#############################################################################	

id_type=strsplit(readLines(InputDataFilename,n=1), "\t")[[1]][1];

cat("ID Type used as column name of rownames:\n");
print(id_type);

for(signf in c("signf_0.100", "signf_0.050", "signf_0.010", "signf_0.001")){

	signf_ix = diff_tab[,signf]==1;
	chosen_var=rownames(diff_tab[signf_ix,,drop=F]);
	out_tab=cbind(rownames(data), data[, chosen_var, drop=F]);
	colnames(out_tab)=c(id_type, chosen_var);
	write.table(
		out_tab, 
		paste(OutputFilenameRoot, ".top_diff.", signf, ".tsv", sep=""),
		row.names=F,
		sep="\t", quote=F
		);

}

#############################################################################	

dev.off();

cat("Done.\n");
print(warnings());
q(status=0);
