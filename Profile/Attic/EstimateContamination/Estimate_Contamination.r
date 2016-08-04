#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"exp_ids_fn", "e", 1, "character",
	"con_ids_fn", "c", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

AbundanceCutoff=1;

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input samples summary table.xls>\n",
	"	-e <experimental sample IDs list>\n",
	"	-c <contamiment/negative control sample IDs list>\n",
	"\n",	
	"This script will use the contaminant profiles/negative controls\n",
	"to adjust the taxa proportions from the experimental samples.\n",
	"\n",
	"\n", sep="");

if(!length(opt$input_file) || !length(opt$exp_ids_fn) || !length(opt$con_ids_fn)){
	cat(usage);
	q(status=-1);
}

###############################################################################

InputSummaryFile=opt$input_file;
ExperimentalIDsFn=opt$exp_ids_fn;
ContaminantIDsFn=opt$con_ids_fn;
OutputFileRoot=gsub("\\.summary_table\\.xls", "", opt$input_file);

cat("\n")
cat("Input Summary File Table: ", InputSummaryFile, "\n", sep="");
cat("Experimental IDs: ", ExperimentalIDsFn, "\n", sep="");
cat("Contaminant IDs: ", ContaminantIDsFn, "\n", sep="");
cat("Output Filename Root: ", OutputFileRoot, "\n", sep="");
cat("\n");


###############################################################################
###############################################################################

load_summary_table=function(summary_table_fn){
	# Load data
	cat("Loading Matrix (", summary_table_fn, ") ...\n", sep="");
	inmat=as.matrix(read.table(summary_table_fn, sep="\t", header=TRUE, check.names=FALSE, row.names=1))

	#cat("\nOriginal Matrix:\n")
	#print(inmat);

	# Grab columns we need into a vector, ignore totals, we won't trust it.
	counts_mat=inmat[,2:(ncol(inmat))];
	#cat("\nCounts Matrix:\n");
	#print(counts_mat);

	num_samples=nrow(counts_mat);
	num_categories=ncol(counts_mat);
	sample_names=rownames(counts_mat);

	cat("\n");
	cat("Num Samples: ", num_samples, "\n");
	cat("Num Categories: ", num_categories, "\n");
	cat("\n");
	return(counts_mat);
}

###############################################################################

load_ids=function(list_fn){
	cat("Loading List (", list_fn, ") ...\n", sep="");
        list=as.vector(read.table(list_fn)[,1]);
	#print(list);
	return(list);
}

###############################################################################

normalize=function(counts){
	# Sum sample totals
	sample_totals=numeric();
	num_samples=nrow(counts);
	num_categories=ncol(counts);

	for(i in 1:num_samples){
		sample_totals[i]=sum(counts_mat[i,]);
	}
	#print(sample_totals);

	# normalize, to compute probabilities
	normalized=matrix(0, nrow=num_samples, ncol=num_categories);
	for(i in 1:num_samples){
		normalized[i,]=counts_mat[i,]/sample_totals[i];
	}

	# Preserve the names
	colnames(normalized)=colnames(counts);
	rownames(normalized)=rownames(counts);

	return(normalized);
}

###############################################################################

plot_RAcurve=function(ordered_composition, title="", top=0, max=1){
	#cat("Plotting: ", title, "\n", sep="");
	names=names(ordered_composition);
	num_cat=length(names);
	if(top<=0){
		top=num_cat;
	}
	barplot(ordered_composition[1:top], names.arg=names[1:top], las=2, cex.names=.5, main=paste("\n\n",title, sep=""), ylim=c(0,max));
}

###############################################################################

fit_contaminant_mixture_model=function(contam_comp, exper_comp){

	num_cat=length(contam_comp);
	cont_idx=contam_comp>0;

	objective=function(p){
		sqr_diff=(exper_comp-p*contam_comp)^2;
		score=sum(sqr_diff[cont_idx]);
		return(score);
	}

	optim_res=optim(par=0, objective, method="Brent", lower=0, upper=3);
	contam_mult=optim_res$par;
	cat("Contaminant Multipler: ", contam_mult, "\n", sep="");
	
	target_removal=contam_mult*contam_comp;
	adjusted_comp=exper_comp-target_removal;
	negative_comp=adjusted_comp<0;

	# Compute normalized filtered composition
	positive_filt_comp=adjusted_comp;
	positive_filt_comp[negative_comp]=0;
	normalized_filt_comp=positive_filt_comp/sum(positive_filt_comp);

	# Compute percent removed from each category
	percent_removed=exper_comp-positive_filt_comp;

	# Output results in structure/list
	fit=list();
	fit$contaminant_multiplier=contam_mult;
	fit$normalized_filtered_composition=normalized_filt_comp;
	fit$percent_removed=percent_removed;
	return(fit);
}

###############################################################################

counts_mat=load_summary_table(InputSummaryFile);
counts_total=apply(counts_mat, 1, sum);

if(T){
	# Clean names
	long_names=colnames(counts_mat);
	split_names=strsplit(long_names, " ");
	num_names=length(long_names);
	short_names=character();
	for(i in 1:num_names){
		short_names[i]=tail(split_names[[i]],1);
	}
	colnames(counts_mat)=short_names;
}

normalized_mat=normalize(counts_mat);

#cat("Normalized.\n");
#print(normalized_mat);
#cat("\n");

exp_ids=load_ids(ExperimentalIDsFn);
cat("Experimental IDs:\n");
print(exp_ids);
cat("\n");

cont_ids=load_ids(ContaminantIDsFn);
cat("Contaminant IDs:\n");
print(cont_ids);
cat("\n");

###############################################################################

cont_norm_mat=normalized_mat[cont_ids,, drop=F];
exp_norm_mat=normalized_mat[exp_ids,, drop=F];

cont_pool=apply(cont_norm_mat, 2, mean);
exp_pool=apply(exp_norm_mat, 2, mean);

order=order(exp_pool, decreasing=T);
#order=order(cont_pool, decreasing=T);

###############################################################################

pdf(paste(OutputFileRoot, ".contaminant_analysis.pdf", sep=""), height=11, width=8.5);

par(mfrow=c(2,1));
par(oma=c(3, 1, .5, 1));

max=max(c(exp_pool, cont_pool))*1.2;
top_disp=100;
plot_RAcurve(exp_pool[order], title="Pooled Experimental", top=top_disp, max=max);
plot_RAcurve(cont_pool[order], title="Pooled Contaminant", top=top_disp, max=max);

###############################################################################

par(mfrow=c(4,1));
num_exp_samples=length(exp_ids);
fits=list();
for(i in 1:num_exp_samples){

	exp_samp=exp_norm_mat[i,];

	fits[[i]]=fit_contaminant_mixture_model(cont_pool, exp_samp);
	fit=fits[[i]];
	cont_mult=fit$contaminant_multiplier;
	filtered=fit$normalized_filtered_composition;
	perc_rem=fit$percent_removed;

	# Plot graphs in PDF file
	max=max(c(exp_samp, cont_pool, filtered))*1.2;
	plot_RAcurve(exp_samp[order], title=sprintf("Input Experimental Sample: %s", exp_ids[i]), top=top_disp, max=max);
	plot_RAcurve(cont_pool[order], title="Pooled Contaminant", top=top_disp, max=max);
	plot_RAcurve(perc_rem[order], title=sprintf("Contaminants Removed\n(Multiplier: %3.2f)", cont_mult), top=top_disp, max=max);
	plot_RAcurve(filtered[order], title="Filtered Experimental Sample", top=top_disp, max=max);

}

###############################################################################
# Output summary file table

filtered_counts_fh=file(paste(OutputFileRoot, ".contam_filt.summary_table.xls", sep=""), "w");

# Output header info
cat(file=filtered_counts_fh, "sample_id", "total", paste(colnames(counts_mat), collapse="\t"), sep="\t");
cat(file=filtered_counts_fh, "\n");

# Output sample info
for(i in 1:num_exp_samples){
	
	sample_id=exp_ids[i];
	fit=fits[[i]];

	filtered_counts=round(counts_mat[sample_id,]-fit$percent_removed*counts_total[sample_id]);

	total=sum(filtered_counts);
	cat(file=filtered_counts_fh, sample_id, total, paste(filtered_counts, collapse="\t"), sep="\t");
	cat(file=filtered_counts_fh, "\n");
}
close(filtered_counts_fh);

###############################################################################
# Output statistics on what was filtered

statistics_fh=file(paste(OutputFileRoot, ".contam_filt.stats.csv", sep=""), "w");
cat(file=statistics_fh, "Sample_ID", "Input_Counts", "Filtered_Counts", "Percent_Filtered", sep=",");
cat(file=statistics_fh, "\n");

for(i in 1:num_exp_samples){
	
	sample_id=exp_ids[i];
	fit=fits[[i]];

	filtered_counts=round(counts_mat[sample_id,]-fit$percent_removed*counts_total[sample_id]);
	total_filtered=sum(filtered_counts);
	perc_filtered=100*(1-total_filtered/counts_total[sample_id]);

	cat(file=statistics_fh, sample_id, counts_total[sample_id], total_filtered, sprintf("%3.2f", perc_filtered), sep=",");
	cat(file=statistics_fh, "\n");

	
}

close(statistics_fh);

###############################################################################

cat("Done.\n\n")
if(!is.null(warnings())){
	print(warnings());
}

q(status=0)
