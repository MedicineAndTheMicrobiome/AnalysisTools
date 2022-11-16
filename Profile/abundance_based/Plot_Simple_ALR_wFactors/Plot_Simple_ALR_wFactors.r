#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"factor_file", "f", 1, "character",
	"factor_subset_fn", "M", 2, "character",
	"top_categories", "T", 2, "numeric",
	"output_file", "o", 2, "character",
	"shorten_category_names", "x", 2, "character",
	"num_rows", "r", 2, "numeric",
	"num_cols", "c", 2, "numeric",
	"tag_name", "t", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

TOP_CATEGORIES=20;
NUM_ROWS=6;
NUM_COLS=5;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	-f <factor file>\n",
	"	[-M <filename for subset of factors/variables to analyze>]\n",
	"	[-T <top categories to analyze, default=", TOP_CATEGORIES, ">]\n",
	"	[-o <output file root name>]\n",
	"	[-x <shorten category names, with separator in double quotes (default=\"\")>]\n",
	"\n",
	"	[-r <number of rows per page, default=", NUM_ROWS, "\n",
	"	[-c <number of rows per page, default=", NUM_COLS, "\n",
	"\n",
	"	[-t <tag name>]\n",
	"\n",
	"	This script will read in the summary table and the factor file.\n",
	"\n",
	"\n", sep="");

if(!length(opt$input_file) || !length(opt$factor_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
FactorFileName=opt$factor_file;

FactorSubset="";
if(length(opt$factor_subset)){
	FactorSubset=opt$factor_subset;
}

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

NumTopCategories=TOP_CATEGORIES;
if(length(opt$top_categories)){
	NumTopCategories=opt$top_categories;
}

NumPlotRows=NUM_ROWS;
if(length(opt$num_rows)){
	NumPlotRows=opt$num_rows;
}

NumPlotCols=NUM_COLS;
if(length(opt$num_cols)){
	NumPlotCols=opt$num_cols;
}

if(length(opt$tag_name)){
        TagName=opt$tag_name;
        cat("Setting TagName Hook: ", TagName, "\n");
        setHook("plot.new",
                function(){
                        cat("Hook called.\n");
                        if(par()$page==T){
                                oma_orig=par()$oma;
                                exp_oma=oma_orig;
                                exp_oma[1]=max(exp_oma[1], 1);
                                par(oma=exp_oma);
                                mtext(paste("[", TagName, "]", sep=""), side=1, line=exp_oma[1]-1,
                                        outer=T, col="steelblue4", font=2, cex=.8, adj=.97);
                                par(oma=oma_orig);
                        }
                }, "append");

}else{
        TagName="";
}

###############################################################################

OutputFileRoot=paste(OutputFileRoot, ".top", NumTopCategories,  sep="");
OutputPDF = paste(OutputFileRoot, ".alr_plots.pdf", sep="");
cat("Output PDF file name: ", OutputPDF, "\n", sep="");

pdf(OutputPDF,width=8.5, height=11);

###############################################################################

load_factors=function(fname){
        factors=data.frame(read.table(fname,  header=TRUE, check.names=FALSE, 
		row.names=1, comment.char="", quote="", sep="\t"));

        dimen=dim(factors);
        cat("Rows Loaded: ", dimen[1], "\n");
        cat("Cols Loaded: ", dimen[2], "\n");
        return(factors);
}

load_summary_file=function(fname){
        cat("Loading Summary Table: ", fname, "\n");
        inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, 
		comment.char="", row.names=1))

        counts_mat=inmat[,2:(ncol(inmat))];
        return(counts_mat);
}

load_list=function(filename){
        val=scan(filename, what=character(), comment.char="#");
        return(val);
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

orig_factors_mat=load_factors(FactorFileName);
#print(factors_mat);
if(FactorSubset!=""){
	subset_list=load_list(FactorSubset);	
	orig_factors_mat=orig_factors_mat[,subset_list,drop=F];
}

###############################################################################

orig_counts_mat=load_summary_file(InputFileName);

###############################################################################
# Reconcile

orig_factors_samples=rownames(orig_factors_mat);
orig_counts_samples=rownames(orig_counts_mat);
shared=intersect(orig_factors_samples, orig_counts_samples);

cat("\n\n");
cat("Samples not represented in summary table file:\n");
excl_to_st=setdiff(orig_counts_samples, shared);
print(excl_to_st);
cat("Samples not represented in offsets file:\n");
excl_to_fct=setdiff(orig_factors_samples, shared);
print(excl_to_fct);
cat("\n\n");

num_shared=length(shared);
cat("Number of Shared Samples: ", num_shared, "\n");

factors_mat=orig_factors_mat[shared,, drop=F];
counts_mat=orig_counts_mat[shared,, drop=F];

###############################################################################

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

###############################################################################

counts_mat=counts_mat+.5;
normalized_mat=normalize(counts_mat);

###############################################################################

max_excl_show=50;

if(length(excl_to_st)>max_excl_show){
       excl_to_st_andmore="...(more)";
}else{
       excl_to_st_andmore="";
}

if(length(excl_to_fct)>max_excl_show){
       excl_to_fct_andmore="...(more)";
}else{
       excl_to_fct_andmore="";
}

plot_text(c(
	paste("Summary Table File: ", InputFileName),
	paste("Factor File: ", FactorFileName),
	paste("Include List: ", FactorSubset),
	paste("Output File Root: ", OutputFileRoot),
	"",
	"Summary Table:",
	paste("      Num Samples:", nrow(orig_counts_mat)),
	paste("   Num Categories:", ncol(orig_counts_mat)),
	"",
	"Factor Table:",
	paste("      Num Samples:", nrow(orig_factors_mat)),
	paste("      Num Factors:", ncol(orig_factors_mat)),
	"",
	paste("Shared Samples:", num_shared),
	"",
	paste("Number of Top Categories from each sample to summarize:", NumTopCategories),
	"",
	"Samples exclusive to Summary Table:",
	capture.output(print(head(excl_to_st, max_excl_show))),
	excl_to_st_andmore,
	"",
	"Samples exclusive to Factor Table:",
	capture.output(print(head(excl_to_fct, max_excl_show))),
	excl_to_fct_andmore
));

###############################################################################

num_categories=ncol(normalized_mat);
factor_names=colnames(factors_mat);
num_factors=ncol(factors_mat);

cat("Target Factors:\n");
print(factor_names);

cat("Num Factors: ", num_factors, "\n", sep="");

#------------------------------------------------------------------------------

compute_alr=function(st, top){
	
	#cat("Num Top: ", top, "\n");
	
	mean_abd=apply(st, 2, mean);
	order_ix=order(mean_abd, decreasing=T);
	st_ordered=st[,order_ix, drop=F];
	Remaining=apply(st_ordered, 1, function(x){ 1-sum(x[1:top]);});
	
	alr=matrix(NA, nrow=nrow(st), ncol=top);
	rownames(alr)=rownames(st);
	colnames(alr)=colnames(st_ordered)[1:top];
	
	for(i in 1:nrow(st)){
		alr[i,]=log(st_ordered[i, 1:top]/Remaining[i]);
		#print(alr[i,]);
	}
	
	#print(mean_abd);
	#print(st_ordered);
	#print(Remaining);
	#print(alr);

	return(alr);

}

get_color=function(pval){
	slp_col="white";
	if(pval<0.001){
		slp_col="red";
	}else if(pval<0.01){
		slp_col="purple";
	}else if(pval<0.05){
		slp_col="blue";
	}else if(pval<0.1){
		slp_col="green";
	}
	return(slp_col);
}


for(i in 1:num_factors){
		
	factor_name=factor_names[i];
	factor_data=factors_mat[,i,drop=F];

	nona_ix=!is.na(factor_data);
	factor_data=factor_data[nona_ix,,drop=F];
	sample_ids=rownames(factor_data);

	cat("Plotting: ", factor_name, "\n");
	#print(factor_data);
	#print(sample_ids);
	
	stn=normalized_mat[sample_ids,,drop=F];

	alr_val=compute_alr(stn, top=NumTopCategories);
	alr_names=colnames(alr_val);

	par(mfrow=c(6,5));
	par(oma=c(2,2,3,2));
	par(mar=c(2,3,3,.1));
	for(p in 1:NumTopCategories){
		hist(alr_val[,p], main=alr_names[p]);
	}
	mtext(
		paste("Histograms for ", factor_name, " w/o NAs", sep=""),
		side=3, outer=T,
		cex=2, font=2
	);

	slp_lwd=2;

	#######################################################################
	# Plot individual

	par(mfrow=c(NumPlotRows, NumPlotCols));
	for(p in 1:NumTopCategories){
		#print(factor_data);
		#print(alr_val[,p]);

		x=factor_data[,1];
		y=alr_val[,p];

		data_range=range(x);
		margin=diff(data_range)*.15;
	
		lmfit=lm(y~x);
		sumfit=summary(lmfit);
		slope=sumfit$coefficients["x", "Estimate"];
		pval=sumfit$coefficients["x", "Pr(>|t|)"];
	
		plot(x, y, 
			xlim=c(data_range[1]-margin, data_range[2]+margin),
			main=alr_names[p]
		);
		title(main=paste("slope=", round(slope,3) , ", p-val=", round(pval,3), sep=""),
			line=.35, cex.main=.75, font=1
		);

		slp_col=get_color(pval);

		abline(lmfit, col="black", lwd=slp_lwd+1);
		abline(lmfit, col=slp_col, lty="longdash", lwd=slp_lwd);
	}
	mtext(
		paste("Scatter Plots for ", factor_name, sep=""),
		side=3, outer=T,
		cex=2, font=2
	);


	#######################################################################
	# Plot with a common Y/ALR axis

	par(mfrow=c(NumPlotRows, NumPlotCols));
	par(mar=c(2,1.5,3,.1));
	all_alr_range=range(alr_val);
	armar=diff(all_alr_range)*.15;
	for(p in 1:NumTopCategories){
		#print(factor_data);
		#print(alr_val[,p]);

		x=factor_data[,1];
		y=alr_val[,p];

		data_range=range(x);
		margin=diff(data_range)*.15;
	
		lmfit=lm(y~x);
		sumfit=summary(lmfit);
		slope=sumfit$coefficients["x", "Estimate"];
		pval=sumfit$coefficients["x", "Pr(>|t|)"];

		if((p-1)%%NumPlotCols==0){
			first_col=T;	
		}else{
			first_col=F;
		}
	
		plot(x, y, 
			xlim=c(data_range[1]-margin, data_range[2]+margin),
			ylim=c(all_alr_range[1]-armar, all_alr_range[2]+armar),
			main=alr_names[p],
			bty="n",
			yaxt=ifelse(first_col, "s", "n"),
			ylab=ifelse(first_col, "ALR", "")
		);


		title(main=paste("slope=", round(slope,3) , ", p-val=", round(pval,3), sep=""),
			line=.35, cex.main=.75, font=1
		);

		slp_col=get_color(pval);

		abline(lmfit, col="black", lwd=slp_lwd+1);
		abline(lmfit, col=slp_col, lty="longdash", lwd=slp_lwd);
	}
	mtext(
		paste("Scatter Plots for ", factor_name, " Common ALR Axis", sep=""),
		side=3, outer=T,
		cex=2, font=2
	);

	#######################################################################
	# Plot with rank abundance style
	cat("Plotting in rank abundance style.\n");

	par(mfrow=c(NumPlotRows/2, 1));
	par(mar=c(10,4,4,2.5));
	all_alr_range=range(alr_val);
	alrmar=diff(all_alr_range)*.15;

	plot(0, type="n", 
		xlim=c(.25, NumTopCategories+.25), 
		ylim=c(all_alr_range[1]-armar, all_alr_range[2]+armar),
		bty="l",
		xaxt="n",
		xlab="", ylab="ALR"
	);

	x_range=range(factor_data);
	x_mid=(x_range[2]-x_range[1])/2;
	x_span=x_range[2]-x_range[1];

	transform_pts=function(x, mid, span, scale, offset){
		(x-mid)/span*scale+offset
	}

	cat("x-mid = ", x_mid, " / x-span", x_span, "\n");

	scale=.75;

	for(p in 1:NumTopCategories){

		cat(alr_names[p], "\n");
		x=factor_data[,1];
		y=alr_val[,p];

		min_y=min(y);

		lmfit=lm(y~x);
		sumfit=summary(lmfit);

		slope=sumfit$coefficients["x", "Estimate"];
		intercept=sumfit$coefficients["(Intercept)", "Estimate"];
		pval=sumfit$coefficients["x", "Pr(>|t|)"];

		slp_col=get_color(pval);

		linex_start=transform_pts(x_range[1], x_mid, x_span, scale, p);
		linex_end=transform_pts(x_range[2], x_mid, x_span, scale, p);

		liney_start=slope*x_range[1]+intercept;
		liney_end=slope*x_range[2]+intercept;

		# Fill
		polygon(
			c(linex_start, linex_end, linex_end, linex_start), 
			c(liney_start, liney_end, all_alr_range[1], all_alr_range[1]), 
			col="blue"
		);

		# slope
		points(c(linex_start, linex_end), c(liney_start, liney_end), 
			col="white", lwd=slp_lwd+2, type="l");
		points(c(linex_start, linex_end), c(liney_start, liney_end), 
			col="black", lwd=slp_lwd+1, type="l");
		points(c(linex_start, linex_end), c(liney_start, liney_end), 
			col=slp_col, lwd=slp_lwd, lty="dashed", type="l");

		# Sample Glyphs
		trans=transform_pts(x, x_mid, x_span, scale, p);
		points(trans, y, cex=.7);

		#abline(lmfit, col="black", lwd=slp_lwd+1);
		#abline(lmfit, col=slp_col, lty="longdash", lwd=slp_lwd);
		#mtext(alr_names[p], 
	}

	bar_width=scale;
	label_size=min(c(1,.7*bar_width/par()$cxy[1]));
	text(
		(1:NumTopCategories)-par()$cxy[1]/2, 
		rep(-par()$cxy[2]/2, NumTopCategories)+(all_alr_range[1]-armar*1.05),
		alr_names, srt=-45, xpd=T, pos=4, cex=label_size);

	mtext(
		paste("Scatter Plots for ", factor_name, " Rank ALR Style", sep=""),
		side=3, outer=T,
		cex=2, font=2
	);

	cat("\n");
	
}

dev.off();


###############################################################################

cat("Done.\n")
warn=warnings();
if(length(warn)){
	print(warn);
}
q(status=0)
