#!/usr/bin/env Rscript

library('getopt');

source("~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r");

params=c(
	"summary_table", "s", 1, "character",
	"output_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-s <input summary_table.tsv>\n",
	"	-o <output file root>\n",
	"\n",	
	"\n",
	"\n");

if(!length(opt$summary_table)){
	cat(usage);
	q(status=-1);
}

options(width=200);

###############################################################################

InputSummaryTable=opt$summary_table;
OutputRoot=opt$output_root;
cnt_adjmt=0.01;

cat("\n");
cat("Input Summary Table: ", InputSummaryTable, "\n");
cat("Output Root: ", OutputRoot, "\n");
cat("Count Adjustment: ", cnt_adjmt, "\n");
cat("\n");

###############################################################################

orig_counts=load_summary_file(InputSummaryTable);

adj_counts=orig_counts+cnt_adjmt;

normalized=normalize(adj_counts);

num_samples=nrow(normalized);
num_categories=ncol(normalized);

#------------------------------------------------------------------------------

mean_abund=apply(normalized, 2, mean);
sort_ix=order(mean_abund, decreasing=T);
normalized=normalized[,sort_ix,drop=F];
adj_counts=adj_counts[,sort_ix,drop=F];
orig_counts=orig_counts[,sort_ix,drop=F];

cat_names=colnames(normalized);

#------------------------------------------------------------------------------

hasRemaining=any(cat_names=="Remaining");

if(hasRemaining){
	cat("Remaining Identified.\n");
}else{
	cat("No Remaining Found.\n");
}

cat("Example Categories:\n");
print(head(cat_names));
cat("\n");

###############################################################################

analyze_points=function(cnt, abd, trans, clr_trans){

	cat_name=colnames(trans);
	cat("Working on: ", cat_name, "\n");

	trans_type=ifelse(clr_trans, "CLR", "ALR");

	zero_ct_smp_ix=(cnt==0);
	num_zeros=sum(zero_ct_smp_ix);
	num_samp=nrow(cnt);
	perc_zeros=num_zeros/num_samp*100.0;
	cat("Num Zeros: ", num_zeros, "\n");

	hrec=hist(log10(cnt+.5), main="", xlab="Log10(Counts)", breaks=20, border=NULL);
	hrec=hist(log10(cnt[zero_ct_smp_ix]+.5), col="red", add=T, breaks=hrec$breaks, border=NULL);

	hrec=hist(log10(abd), main=cat_name, xlab="Log10(Abundance)", breaks=20, lwd=.1, border=NULL);
	hrec=hist(log10(abd[zero_ct_smp_ix]), col="red", add=T, breaks=hrec$breaks, lwd=.1, border=NULL);
	title(main=cat_name);
	title(main=paste("Zero Count/Abundance Samples: ", num_zeros, 
			" ( ", sprintf("%3.2f", perc_zeros), "% )", sep=""),
		cex.main=.8, col.main="red", font.main=2, line=.90
		);

	hrec=hist(trans, main="", xlab=trans_type, breaks=20, lwd=.1, border=NULL);
	hrec=hist(trans[zero_ct_smp_ix], col="red", add=T, breaks=hrec$breaks, lwd=.1, border=NULL);

	cat("-----------------------------------------------------------\n");

}

center_log_transform=function(norm_mat){
	
	calc_geom_mean=function(x){exp(mean(log(x)))};	

	geom_mean=apply(norm_mat, 1, calc_geom_mean);
	mult=sweep(norm_mat, 1, geom_mean, "/");
	lr=log(mult);
	dimnames(lr)=dimnames(norm_mat);
	return(lr);

}

additive_log_transform=function(norm_mat, rem_catnm="Remaining"){

	cat_names=colnames(norm_mat);
	rem_ix=cat_names==rem_catnm;

	rem=norm_mat[, rem_catnm];

	target_mat=norm_mat[,!rem_ix, drop=F];
	mult=sweep(target_mat, 1, rem, "/");
	lr=log(mult);
	dimnames(lr)=dimnames(target_mat);
	return(lr);

}

###############################################################################

if(hasRemaining){
	cat("Performing Additive Log Ratio.\n");
	lr_norm=additive_log_transform(normalized);

	rmn_ix=(cat_names=="Remaining");
	normalized=normalized[,!rmn_ix];
	orig_counts=orig_counts[,!rmn_ix];
	adj_counts=adj_counts[,!rmn_ix];
}else{
	cat("Performing Centered Log Ratio.\n");
	lr_norm=center_log_transform(normalized);
}


out_fn=paste(OutputRoot, ".screen.pdf", sep="");
pdf(out_fn, height=11, width=8.5);


par(mfrow=c(4,3));


num_analys_cat=ncol(normalized);
for(i in 1:num_analys_cat){
	ct=orig_counts[,i,drop=F];
	dp=normalized[,i,drop=F];
	lr=lr_norm[,i,drop=F];
	analyze_points(ct, dp, lr, clr=!hasRemaining);
}


###############################################################################



cat("Done.\n")
print(warnings());

q(status=0)
