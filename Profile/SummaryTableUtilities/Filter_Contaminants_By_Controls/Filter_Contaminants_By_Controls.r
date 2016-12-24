#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_summary_table", "i", 1, "character",
	"avg_contam", "a", 2, "character",
	"paired_contam", "p", 2, "character",
	"plevel", "l", 2, "numeric",
	"output_fname_root", "o", 2, "character",
	"short_name_delim", "d", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

AbundanceCutoff=1;
DEFAULT_PLEVEL=.5;
DEFAULT_DELIM=";";

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input samples summary table.tsv>\n",
	"\n",
	"	One of:\n",
	"	[-a <contamiment/negative control sample IDs list>]\n",
	"	[-p <list of experimental/contaminant paired sample IDs>]\n",
	"\n",
	"	[-l <p-norm used in the objective function, default=", DEFAULT_PLEVEL, ">]\n",
	"	[-o <output filename root>]\n",
	"	[-d <name shortening delimitor, default=", DEFAULT_DELIM, ">]\n",
	"\n",	
	"This script will use the contaminant profiles/negative controls\n",
	"to adjust the taxa proportions from the experimental samples.\n",
	"The goal is essentially to try to find the distribution of the\n",
	"contaminant sample inside the experimental sample and then substract\n",
	"the contaminant distribution from the experimental distribution.\n",
	"\n",
	"The output is a summary table consisting only of the filtered\n",
	"experimental samples, and some diagnostic pdf file.\n",
	"The counts in the summary table is be adjusted to reflected\n",
	"the reduced proportion of reads per sample.\n",
	"\n",
	"The algorithm works by using a mixture model consisting of\n",
	"an unknown ratio of experimental to contaminant.\n",
	"\n",
	"	PureExperimental+ c*Contaminant = ObservedExperimental\n",
	"\n",
	"The goal is to try to calculate the coefficient c.\n",
	"We will calculate c by minimizing the following objective function:\n",
	"\n",
	"	Sum (abs(Exp[i]-c*Cont[i]))^p\n",
	"	 where i is only the categories in the Cont\n",
	"\n",
	"	With the constraint that c is in [0,1].\n",
	"\n",
	"	p is the p-norm degree:\n",
	"	When p=2, the fit is a least squares.\n",
	"	When p=1, the fit is absolute value.\n",
	"\n",
	"	The smaller the p, e.g. .5, the less the effect of outliers,\n",
	"	but the fit can become less robust.\n",
	"	For contaminant removal, outliers may be actual levels of\n",
	"	non-contaminant, so we want the fit to be spread out across\n",
	"	all categories, and not dominated by a few categories\n",
	"	that fit poorly.\n",
	"\n",
	"\n",
	"There are two options to apply the contaminant filter:\n",
	"	1.) Averaging (-a): Take an average of contaminants and apply a single\n",
	"		profile as the contaminent to each experimental sample.\n",
	"	2.) Paired (-p): If each sample is paired with a negative control,\n",
	"		then apply a single negative control with a single\n",
	"		experimental sample, e.g. BAL / BAL Control.\n",
	"\n",
	"The format for the -a negative control file:\n",
	"	<contam1 sample id>\\n\n",
	"	<contam2 sample id>\\n\n",
	"	<contam3 sample id>\\n\n",
	"	...\n",
	"	<contamn sample id>\\n\n",
	"\n",
	"The format for the -p paired exp/cont file:\n",
	"	<exper1 sample id>\\t<contam1 sample id>\\n\n",	
	"	<exper2 sample id>\\t<contam2 sample id>\\n\n",	
	"	<exper3 sample id>\\t<contam3 sample id>\\n\n",	
	"	...\n",
	"	<expern sample id>\\t<contamn sample id>\\n\n",	
	"\n",
	"\n", sep="");

if(!length(opt$input_summary_table)){
	cat(usage);
	q(status=-1);
}

###############################################################################

InputSummaryTable=opt$input_summary_table;
AvgContamFName=opt$avg_contam;
PairedContamFName=opt$paired_contam;
PLevel=opt$plevel;
OutputFileRoot=opt$output_fname_root;
ShortNameDelim=opt$short_name_delim;

if(!length(OutputFileRoot)){
	OutputFileRoot=gsub("\\.summary_table\\.xls", "", InputSummaryTable);
	OutputFileRoot=gsub("\\.summary_table\\.tsv", "", OutputFileRoot);
}

if(!length(PLevel)){
	PLevel=DEFAULT_PLEVEL;
}

if(!length(ShortNameDelim)){
	ShortNameDelim=DEFAULT_DELIM;
}

cat("\n")
cat("Input Summary File Table: ", InputSummaryTable, "\n", sep="");
cat("Output Filename Root: ", OutputFileRoot, "\n", sep="");
cat("P-norm level: ", PLevel, "\n", sep="");
cat("Short Name Delim: ", ShortNameDelim, "\n", sep="");
cat("\n");

doPaired=NULL;
if(length(PairedContamFName)){
	cat("Performing paired contaminant removal: ", PairedContamFName, "\n");
	doPaired=T;
}else{
	cat("Performing averaged contaminant removal: ", AvgContamFName, "\n");
	doPaired=F;
}


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

load_pairs=function(pairs_fn){
	cat("Loading Pairs (", pairs_fn, ") ...\n", sep="");
	pairs=read.delim(pairs_fn, sep="\t", row.names=1, header=F, as.is=T);
	src_names=rownames(pairs);
	dst_names=pairs[,1];
	names(dst_names)=src_names;
	return(dst_names);
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
	barplot(ordered_composition[1:top], names.arg=names[1:top], las=2, 
		cex.names=.75, cex.axis=.75,
		main=paste("\n\n",title, sep=""), ylim=c(0,max));
}

###############################################################################

fit_contaminant_mixture_model=function(contam_comp, exper_comp, plevel){

	num_cat=length(contam_comp);
	cont_idx=contam_comp>0;

	# Extract the categories that are non-zero in contaminant
	expr_sub=exper_comp[cont_idx];
	cont_sub=contam_comp[cont_idx];

	# Define the objective function
	objective=function(p){
		diff=(abs((1-p)*expr_sub-p*cont_sub))^plevel;
		score=sum(diff);
		return(score);
	}

	# Perform the optimization
	optim_res=optim(par=0, objective, method="Brent", lower=0, upper=1);
	contam_prop=optim_res$par;
	cat("Contaminant Proportion: ", contam_prop, "\n", sep="");
	contam_scale=contam_prop/(1-contam_prop);
	cat("Contaminant Scale: ", contam_scale, "\n", sep="");

	# Remove contaminant from experimentals
	adjusted_comp=exper_comp - contam_scale*contam_comp;

	# Remove/zero out negative composition
	negative_cat_ix=adjusted_comp<0;
	adjusted_comp[negative_cat_ix]=0;

	# Compute normalized filtered composition (so distribution adds up to 1 again)
	normalized_filt_comp=adjusted_comp/sum(adjusted_comp);

	# Estimate Proportion removed
	proportion_removed=sum(exper_comp-adjusted_comp);

	# Output results in structure/list
	fit=list();
	fit$contaminant_multiplier=contam_scale;
	fit$contaminant_proportion=contam_prop;
	fit$normalized_filtered_composition=normalized_filt_comp;
	fit$proportion_removed=proportion_removed;
	return(fit);
}

###############################################################################

counts_mat=load_summary_table(InputSummaryTable);
counts_total=apply(counts_mat, 1, sum);
summary_table_sample_ids=rownames(counts_mat);

# Clean names
long_names=colnames(counts_mat);
split_names=strsplit(long_names, ShortNameDelim);
num_names=length(long_names);
short_names=character();
for(i in 1:num_names){
	short_names[i]=tail(split_names[[i]],1);
}

normalized_mat=normalize(counts_mat);


# Load IDs
AvgContamFName=opt$avg_contam;
PairedContamFName=opt$paired_contam;

avg_cont_dist=NULL;
experm_samples=character();
contam_samples=character();

if(doPaired){
	pairs=load_pairs(PairedContamFName);
	experm_samples=intersect(names(pairs), summary_table_sample_ids);
	contam_samples=intersect(pairs, summary_table_sample_ids);
}else{
	ids_array=load_ids(AvgContamFName);
	overlapping_ids=intersect(ids_array, summary_table_sample_ids);

	contam_mat=normalized_mat[overlapping_ids,, drop=F];
	avg_cont_dist=apply(normalized_mat, 2, mean);

	experm_samples=setdiff(summary_table_sample_ids, ids_array);
	contam_samples=overlapping_ids;
}

cat("Experimental Samples:\n");
print(experm_samples);

cat("Contaminant Samples:\n");
print(contam_samples);

# Reorder the normalized matrix so that we'll see both exp and ctl 
# categories in the rank abundance plot
print(experm_samples);
print(rownames(normalized_mat));
mean_exp_abund=apply(normalized_mat[experm_samples,,drop=F], 2, mean);
mean_con_abund=apply(normalized_mat[contam_samples,,drop=F], 2, mean);
mean_both_abund=(mean_exp_abund+mean_con_abund)/2;
mean_order=order(mean_both_abund, decreasing=T);

normalized_mat=normalized_mat[,mean_order];
short_names=short_names[mean_order];
original_names=colnames(normalized_mat);
colnames(normalized_mat)=short_names;

if(!doPaired){
	avg_cont_dist=avg_cont_dist[mean_order];
	names(avg_cont_dist)=short_names;
}


###############################################################################

pdf(paste(OutputFileRoot, ".dist_contam.pdf", sep=""), height=8.5, width=14);

top_cat_to_plot=35;

# Plot 3 columns of images:
# a.) Experimental
# b.) Control
# c.) Subracted

###############################################################################

par(mfrow=c(4,3));
mar=par()$mar;
mar=c(8.1, 3.1, 2.1, 1);
#mar=c(6.1, 3.1, 3.1, 1);
par(mar=mar);

#num_exp_samples=length(exp_ids);
fits=list();
for(exp_samp_id in experm_samples){

	exp_dist=normalized_mat[exp_samp_id,];

	if(doPaired){
		ctl_name=pairs[exp_samp_id];
		ctl_dist=normalized_mat[ctl_name,];
	}else{
		ctl_name="average control";
		ctl_dist=avg_cont_dist;
	}

	fit=fit_contaminant_mixture_model(ctl_dist, exp_dist, PLevel);
	mult=fit$contaminant_multiplier;
	prop=fit$contaminant_proportion;
	filt_dist=fit$normalized_filtered_composition;
	removed=fit$proportion_removed;

	max_abund=max(exp_dist, ctl_dist, filt_dist);
	max_y=max_abund*1.1;
	
	plot_RAcurve(exp_dist, title=paste("Experimental:", exp_samp_id), top=top_cat_to_plot, max=max_y);
	plot_RAcurve(ctl_dist, title=paste("Control:", ctl_name), top=top_cat_to_plot, max=max_y);
	plot_RAcurve(filt_dist, title="Corrected", top=top_cat_to_plot, max=max_y);
	mtext(paste("Proportion Removed:", round(removed, 3)), line=-1.75, outer=F, cex=.5);
	mtext(paste("Multiplier:", round(mult, 3)), line=-2.5, outer=F, cex=.5);

}

###############################################################################
# Output summary file table

filtered_counts_fh=file(paste(OutputFileRoot, ".contam_filt.summary_table.tsv", sep=""), "w");

# Output header info
cat(file=filtered_counts_fh, "sample_id", "total", paste(colnames(counts_mat), collapse="\t"), sep="\t");
cat(file=filtered_counts_fh, "\n");

# Output sample info
#for(i in 1:num_exp_samples){
#	
#	sample_id=exp_ids[i];
#	fit=fits[[i]];
#
#	filtered_counts=round(counts_mat[sample_id,]-fit$percent_removed*counts_total[sample_id]);
#
#	total=sum(filtered_counts);
#	cat(file=filtered_counts_fh, sample_id, total, paste(filtered_counts, collapse="\t"), sep="\t");
#	cat(file=filtered_counts_fh, "\n");
#}
#close(filtered_counts_fh);
#
################################################################################
## Output statistics on what was filtered
#
#statistics_fh=file(paste(OutputFileRoot, ".contam_filt.stats.csv", sep=""), "w");
#cat(file=statistics_fh, "Sample_ID", "Input_Counts", "Filtered_Counts", "Percent_Filtered", sep=",");
#cat(file=statistics_fh, "\n");
#
#for(i in 1:num_exp_samples){
#	
#	sample_id=exp_ids[i];
#	fit=fits[[i]];
#
#	filtered_counts=round(counts_mat[sample_id,]-fit$percent_removed*counts_total[sample_id]);
#	total_filtered=sum(filtered_counts);
#	perc_filtered=100*(1-total_filtered/counts_total[sample_id]);
#
#	cat(file=statistics_fh, sample_id, counts_total[sample_id], total_filtered, sprintf("%3.2f", perc_filtered), sep=",");
#	cat(file=statistics_fh, "\n");
#
#	
#}
#
#close(statistics_fh);

###############################################################################

cat("Done.\n\n")
if(!is.null(warnings())){
	print(warnings());
}

q(status=0)
