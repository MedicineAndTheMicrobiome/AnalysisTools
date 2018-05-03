#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"top_categories", "t", 2, "character",
	"output_file", "o", 2, "character",
	"shorten_category_names", "s", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

DEF_DIVERSITY="tail";
TOP_CATEGORIES=4;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	[-t <top categories to display, default=", TOP_CATEGORIES, ">]\n",
	"	[-o <output file root name>]\n",
	"	[-s <shorten category names, with separator in double quotes (default=\"\")>]\n",
	"\n",
	"	This script will readin the summary table and generate pair-wise\n",
	"	comparisons across all samples.\n",
	"\n",
	"	The p-values reflect the variation from resampling the same\n",
	"	subject and should be use to infer the characteristics of the\n",
	"	the cohort the subject was selected from.\n",	
	"\n",
	"\n", sep="");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
FactorFileName=opt$factor_file;

NumTopCategories=TOP_CATEGORIES;
if(length(opt$top_categories)){
	NumTopCategories=opt$top_categories;
}

if(length(opt$output_file)>0){
	OutputFileRoot=opt$output_file;
}else{
	OutputFileRoot=InputFileName;
	OutputFileRoot=gsub("\\.summary_table\\.tsv$", "", OutputFileRoot);
	OutputFileRoot=gsub("\\.summary_table\\.xls$", "", OutputFileRoot);
	cat("No output file root specified.  Using input file name as root.\n");
}

if(length(opt$shorten_category_names)){
        ShortenCategoryNames=opt$shorten_category_names;
}else{
        ShortenCategoryNames="";
}

if(ShortenCategoryNames==TRUE){
        cat("Error:  You need to specify a delimitor to split the category names.\n");
        cat("        i.e., this -x option is not a flag, it requires a parameter.\n");
        quit(status=-1);
}

###############################################################################

OutputPDF = paste(OutputFileRoot, ".pw_cmp.pdf", sep="");
pdf(OutputPDF,width=8.5,height=11)

param_txt=capture.output(
	{
	cat("Parameters:\n");
	cat("  Input File Name: ", InputFileName, "\n", sep="");
	cat("  Output File Root Name: ", OutputFileRoot, "\n", sep="");
	cat("  Shorten Categories: ", ShortenCategoryNames, "\n", sep="");
	cat("  Output PDF file name: ", OutputPDF, "\n", sep="");
	cat("\n");
	});
cat(paste(param_txt, collapse="\n"));

###############################################################################

load_summary_file=function(fname){
        cat("Loading Summary Table: ", fname, "\n");
        inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", row.names=1))
        counts_mat=inmat[,2:(ncol(inmat))];
        return(counts_mat);
}

normalize=function(counts){
        totals=apply(counts, 1, sum);
        num_samples=nrow(counts);
        normalized=matrix(0, nrow=nrow(counts), ncol=ncol(counts));

        for(i in 1:num_samples){
                normalized[i,]=counts[i,]/totals[i];
        }

        colnames(normalized)=colnames(counts);
        rownames(normalized)=rownames(counts);
        return(normalized);
}

plot_text=function(strings){

	orig.par=par(no.readonly=T);
        par(family="Courier");
        par(oma=rep(.1,4));
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

###############################################################################

# Get color assignments
get_colors=function(num_col, alpha=1){
	colors=hsv(seq(0,1,length.out=num_col+1), c(1,.5), c(1,.75,.5), alpha=alpha);
	color_mat_dim=ceiling(sqrt(num_col));
	color_pad=rep("grey", color_mat_dim^2);
	color_pad[1:num_col]=colors[1:num_col];
	color_mat=matrix(color_pad, nrow=color_mat_dim, ncol=color_mat_dim);
	colors=as.vector(t(color_mat));
	colors=colors[colors!="grey"];
}

# Generate bootstrap samples
make_bootstrap_samples=function(counts_summary_table, num_bs, normalize=F){

	num_samples=nrow(counts_summary_table);
	num_categoires=ncol(counts_summary_table);
	sample_names=rownames(counts_summary_table);

	totals=apply(counts_summary_table, 1, sum);

	# Normalize counts
        normalized=matrix(0, nrow=nrow(counts_summary_table), ncol=ncol(counts_summary_table));
	colnames(normalized)=colnames(counts_summary_table);
        for(i in 1:num_samples){
                normalized[i,]=counts_summary_table[i,]/totals[i];
	}

	bs_samp_list=list();
		
	for(i in 1:num_samples){
		cur_samp_id=sample_names[i];

		resamp_mat=matrix(0, nrow=num_bs, ncol=num_categoires);
		colnames(resamp_mat)=colnames(counts_summary_table);
		rownames(resamp_mat)=paste(cur_samp_id, sprintf("%05i", 1:num_bs), sep="");

		for(bsix in 1:num_bs){
			counts=as.vector(rmultinom(1, size=totals[i], prob=normalized[i,]));

			if(normalize==T){
				resamp_mat[bsix,]=counts/totals[i];
			}else{
				resamp_mat[bsix,]=counts;
			}
		}

		bs_samp_list[[cur_samp_id]]=resamp_mat;

		#print(counts_summary_table[i,]);
		#print(bs_samp_list[[cur_samp_id]]);
		#cat("---------------------------------------------\n");

	}
	
	return(bs_samp_list);

}

calc_sorted_pair_ratios=function(counts_a, counts_b){
	# calculate median and 95% ci of log ratios of (a+1)/(b+1)
	num_samples=nrow(counts_a);
	num_categories=ncol(counts_a);

	if(num_samples!=nrow(counts_b)){
		cat("Error: Must have identical samples sizes.\n");
		quit(-1);
	}
	if(num_categories!=ncol(counts_b)){
		cat("Error: Must have identical number of categories.\n");
		quit(-1);
	}

	medians=rep(0, num_categories);
	lb95=rep(0, num_categories);
	ub95=rep(0, num_categories);
	
	# Calculate median and bounds across categories
	for(cat_ix in 1:num_categories){

		# Calculate logratios across samples
		lograt=rep(0, num_samples);
		for(samp_ix in 1:num_samples){

			counts_a_p1=counts_a[samp_ix,cat_ix]+1;
			counts_b_p1=counts_b[samp_ix,cat_ix]+1;

			lograt[samp_ix]=log2(counts_a_p1/counts_b_p1);

		}
		
		# Calculate mediand and 95% CI
		stats=quantile(lograt, c(.025, .5, .975));
		lb95[cat_ix]=stats[1];
		medians[cat_ix]=stats[2];
		ub95[cat_ix]=stats[3];
	}

	# Sort
	ix=order(medians, decreasing=T);

	# Prep return values
	return_list=list();
	return_list[["categories"]]=colnames(counts_a)[ix];
	return_list[["medians"]]=medians[ix];
	return_list[["LB95"]]=lb95[ix];
	return_list[["UB95"]]=ub95[ix];

	return(return_list);

}

plot_ratios=function(ratio_rec, samp_a_name, samp_b_name, title, which){

	cat("Plotting ratios for: ", samp_a_name, " vs ", samp_b_name, ": ", title,"\n");

	if(all(which==F)){
		cat("Nothing to plot...\n");
		return();
	}
	
	category_names=ratio_rec[["categories"]][which];
	medians=ratio_rec[["medians"]][which];
	LB95=ratio_rec[["LB95"]][which];
	UB95=ratio_rec[["UB95"]][which];
	num_cat_to_plot=sum(which);

	# Compute plot ranges
	ymin=min(LB95, UB95, medians);
	ymax=max(LB95, UB95, medians);
	yrange=ymax-ymin

	# If we don't go above or below zero at all, then don't pad with 20%
	ylimits=c(ymin-yrange*.05, ymax+yrange*.05);
	if(ymax==0){
		ylimits[2]=.01;
	}
	if(ymin==0){
		ylimits[1]=-.01;
	}

	label_size=40/num_cat_to_plot
	label_size=min(1, label_size);

	plot_name=paste(samp_a_name, " vs. ", samp_b_name, ": ", title, sep="");
	mids=barplot(medians, names.arg=category_names, 
		main=plot_name, ylim=ylimits, cex.names=label_size,
		las=2);

	
	for(i in 1:num_cat_to_plot){
		points(c(mids[i], mids[i]), c(LB95[i], UB95[i]), 
		cex=label_size, type="b", lwd=.5, col="grey");
	}
	
}

###############################################################################

counts_mat=load_summary_file(InputFileName);
num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);

st_info_txt=capture.output(
	{
	cat("Summary Table:\n");
	cat("  Num Samples: ", num_samples, "\n", sep="");
	cat("  Num Categories: ", num_categories, "\n", sep="");
	cat("\n");
	}
);

cat(paste(st_info_txt, collapse="\n"));

plot_text(c(
	param_txt,
	"",
	st_info_txt
	));

###############################################################################

sample_names=rownames(counts_mat);

if(ShortenCategoryNames!=""){
        full_names=colnames(counts_mat);
        splits=strsplit(full_names, ShortenCategoryNames);
        short_names=character();
        for(i in 1:length(full_names)){
                short_names[i]=tail(splits[[i]], 1);

                short_names[i]=gsub("_unclassified$", "_uncl", short_names[i]);
                short_names[i]=gsub("_group", "_grp", short_names[i]);
        }
        colnames(counts_mat)=short_names;

        cat("Names have been shortened.\n");
}

normalized_mat=normalize(counts_mat);

###############################################################################

idx=1:num_categories;

cat("Generating bootstrapped samples:\n");
bs_samp_list=make_bootstrap_samples(counts_mat, num_bs=2000, normalize=F);

par(mar=c(40,4,4,1));

cat("Starting pairwise comparisons:\n");
for(samp_a_ix in 1:num_samples){

	samp_a_name=sample_names[samp_a_ix];

	for(samp_b_ix in 1:num_samples){

		if(samp_b_ix==samp_a_ix){
			next;
		}

		samp_b_name=sample_names[samp_b_ix];

		cat(samp_a_name, " vs. ", samp_b_name, "\n", sep="");

		ratios=calc_sorted_pair_ratios(bs_samp_list[[samp_a_name]], bs_samp_list[[samp_b_name]]);


		plot(0,0, type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="");
		text(0,0, paste(samp_a_name, "\nvs.\n", samp_b_name, sep=""), cex=4);

		cat("Plotting histogram...\n");
		hist(ratios[["medians"]], 
			xlab="Log2 Ratios",
			main=paste(samp_a_name, " vs. ", samp_b_name, 
				": Histogram of all Log Ratios", sep=""));

		cat("Plotting top...\n");
		for(top in c(10, 15, 25, 40, 100, 101)){

			plot_ratios(
				ratios, samp_a_name, samp_b_name,
				title=sprintf("Top %i", top),
				which=(idx<=top)
			);

		}

		cat("Plotting subsets...\n");
		plot_ratios(
			ratios, samp_a_name, samp_b_name,
			title="Log Ratio Medians > 0",
			which=ratios[["medians"]]>0
		);

		plot_ratios(
			ratios, samp_a_name, samp_b_name,
			title="Log Ratio Medians < 0",
			which=ratios[["medians"]]<0
		);

		plot_ratios(
			ratios, samp_a_name, samp_b_name,
			title="Log Ratio LB 95 > 0",
			which=ratios[["LB95"]]>0
		);

		plot_ratios(
			ratios, samp_a_name, samp_b_name,
			title="Log Ratio LB < 0 and UB > 0",
			which=(ratios[["LB95"]]<0) & (ratios[["UB95"]]>0)
		);

		plot_ratios(
			ratios, samp_a_name, samp_b_name,
			title="Log Ratio UB 95 < 0",
			which=ratios[["UB95"]]<0
		);

		
		
	}
}


###############################################################################

dev.off();

cat("Done.\n")
warn=warnings();
if(length(warn)){
	print(warn);
}
q(status=0);

