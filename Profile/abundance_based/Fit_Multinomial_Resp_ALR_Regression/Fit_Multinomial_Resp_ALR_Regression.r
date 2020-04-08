#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
library(car);

source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');

options(useFancyQuotes=F);


params=c(
	"summary_file", "s", 1, "character",
	"num_variables", "p", 2, "numeric",
	"additional_categories", "a", 2, "character",
	"contains_remaining", "R", 2, "logical",
	"shorten_category_names", "x", 2, "character",

	"factors", "f", 1, "character",
	"covariates", "c", 1, "character",
	"responses", "y", 1, "character",
	"required", "q", 2, "character",
	"reference_levels", "r", 2, "character",

	"time_column", "t", 2, "character",
	"subject_column", "j", 2, "character",

	"outputroot", "o", 2, "character"

);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

NUM_TOP_CAT=30;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"  Summary Table Parameters:\n",
	"	-s <summary file table for taxa/function (used as Predictor/X's)>\n",
	"	[-p <number of top taxonomic/categorical variables, default=", NUM_TOP_CAT, ">]\n",
	"	[-a <list of additional categories (from summary table) to include (filename)>]\n",
	"	[-R (pay attention to 'remaining' category)]\n",
	"	[-x <shorten category names, with separator in double quotes (default=\"\")>]\n",
	"\n",
	"  Factor/Metadata Parameters:\n",
	"	-f <factors file, contains covariates and multinomial Y>\n",
	"	-c <list of covariate X's names to select from factor file (filename)>\n",
	"	-y <Categorical (multinomial) response Y's to select from factor file (filename)>\n",
	"	[-q <required list of variables to include after NA removal>]\n",
	"	[-r <reference levels file for variables in factor file>]\n",
	"	[-t <time variable>]\n",
	"  Time-related Variables:\n",
	"	[-t <name of time column, default=NA>]\n",
	"	[-j <name of subject column, default=NA>]\n",
	"\n",
	"  Output Parameters:\n",
	"	-o <output filename root>\n",
	"\n",
	"This script will fit the following model:\n",
	"\n",
	" Multinomial Response = covariates + MALR(top-p taxonomy)\n",
	"\n",
	"If the -R flag is set, a 'remaining' category will be included in the denominator\n",
	"	independent of how large it is.  I.e., do not use it as one of the predictor variables.\n",
	"\n",
	"\n",
	"If the time options, -t and -j, are not specified, it assumed just one time will be considered\n",
	"  i.e., cross-sectional.\n",
	"\n", sep="");

if(!length(opt$summary_file) || !length(opt$factors) || !length(opt$covariates) || !length(opt$responses)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputRoot=gsub(".summary_table.xls", "", opt$summary_file);
	OutputRoot=gsub(".summary_table.tsv", "", OutputRoot);
}else{
	OutputRoot=opt$outputroot;
}

if(!length(opt$num_variables)){
	NumVariables=NUM_TOP_CAT;
}else{
	NumVariables=opt$num_variables;
}

if(!length(opt$reference_levels)){
        ReferenceLevelsFile="";
}else{
        ReferenceLevelsFile=opt$reference_levels;
}

if(length(opt$contains_remaining)){
	UseRemaining=T;
}else{
	UseRemaining=F;
}

if(length(opt$shorten_category_names)){
	ShortenCategoryNames=opt$shorten_category_names;
}else{
	ShortenCategoryNames="";
}

if(length(opt$required)){
	RequiredFile=opt$required;
}else{
	RequiredFile="";
}

if(length(opt$additional_categories)){
	AdditionalCatFile=opt$additional_categories;
}else{
	AdditionalCatFile="";
}

if(length(opt$time_column)){
	TimeColumn=opt$time_column;
}else{
	TimeColumn=character(0);
}

if(length(opt$subject_column)){
	SubjectColumn=opt$subject_column;
}else{
	SubjectColumn=character(0);
}

SummaryFile=opt$summary_file;
FactorsFile=opt$factors;
CovariatesFile=opt$covariates;
ResponseColname=opt$responses;

cat("\n");
cat("      Summary File: ", SummaryFile, "\n", sep="");
cat("      Factors File: ", FactorsFile, "\n", sep="");
cat("   Covariates File: ", CovariatesFile, "\n", sep="");
cat("  Response Colname: ", ResponseColname, "\n", sep="");
cat("     Required File: ", RequiredFile, "\n", sep="");
cat("   Additional File: ", AdditionalCatFile, "\n", sep="");
cat("\n");
cat("      Time Colname: ", TimeColumn, "\n", sep="");
cat("   Subject Colname: ", SubjectColumn, "\n", sep="");
cat("\n");
cat("    Output File: ", OutputRoot, "\n", sep="");
cat("\n");
cat("Number of MALR Variables: ", NumVariables, "\n", sep="");
cat("Reference Levels File: ", ReferenceLevelsFile, "\n", sep="");
cat("Use Remaining? ", UseRemaining, "\n");
cat("Shorten Category Names: '", ShortenCategoryNames, "'\n", sep="");
cat("\n");

if(ShortenCategoryNames==TRUE){
	cat("Error:  You need to specify a delimitor to split the category names.\n");
	cat("        i.e., this -x option is not a flag, it requires a parameter.\n");
	quit(status=-1);
}

options(width=100);
cat("Text Line Width: ", options()$width, "\n", sep="");

##############################################################################

load_factors=function(fname){
	cat("Loading Factors/Metadata: ", fname, "\n", sep="");
	factors=data.frame(read.table(fname,  sep="\t", header=TRUE, row.names=1, check.names=FALSE));
	factor_names=colnames(factors);

	ignore_idx=grep("^IGNORE\\.", factor_names);

	if(length(ignore_idx)!=0){
		return(factors[-ignore_idx]);
	}else{
		return(factors);
	}
}

load_summary_file=function(fname){
	inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", 
		quote="", row.names=1))
	counts_mat=inmat[,2:(ncol(inmat))];

	# Clean category names a little
	cat_names=colnames(counts_mat);
	cat_names=gsub("-", "_", cat_names);
	cat_names=gsub("\\[", "", cat_names);
	cat_names=gsub("\\]", "", cat_names);
	colnames(counts_mat)=cat_names;
	
	cat("Num Categories in Summary Table: ", ncol(counts_mat), "\n", sep="");
	return(counts_mat);
}

load_reference_levels_file=function(fname){
        inmat=as.matrix(read.table(fname, sep="\t", header=F, check.names=FALSE, comment.char="#", row.names=1))
        colnames(inmat)=c("ReferenceLevel");
        print(inmat);
        cat("\n");
        if(ncol(inmat)!=1){
                cat("Error reading in reference level file: ", fname, "\n");
                quit(status=-1);
        }
        return(inmat);
}

relevel_factors=function(factors, ref_lev_mat){
        num_factors_to_relevel=nrow(ref_lev_mat);
        relevel_names=rownames(ref_lev_mat);
	factor_names=colnames(factors);
        for(i in 1:num_factors_to_relevel){
		
		target_relev_name=relevel_names[i];
		if(any(target_relev_name==factor_names)){
			tmp=factors[,target_relev_name];
			#print(tmp);
			tmp=relevel(tmp, ref_lev_mat[i, 1]);
			#print(tmp);
			factors[,target_relev_name]=tmp;
		}else{
			cat("Note: ", target_relev_name, 
				" not in model.  Ignoring reference releveling.\n\n", sep="");
		}
        }
        return(factors);
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

load_list=function(filename){
	cat("Loading List: ", filename, "\n", sep="");
	val=scan(filename, what=character(), comment.char="#");
	return(val);
}

extract_top_categories=function(ordered_normalized, top, additional_cat=c()){

        num_samples=nrow(ordered_normalized);
        num_categories=ncol(ordered_normalized);

        cat("Samples: ", num_samples, "\n");
        cat("Categories: ", num_categories, "\n");

        num_top_to_extract=min(num_categories-1, top);

        cat("Top Requested to Extract: ", top, "\n");
        cat("Columns to Extract: ", num_top_to_extract, "\n");

        # Extract top categories requested
        top_cat=ordered_normalized[,1:num_top_to_extract];

        if(length(additional_cat)){
                cat("Additional Categories to Include:\n");
                print(additional_cat);
        }else{
                cat("No Additional Categories to Extract.\n");
        }

        # Extract additional categories
        # :: Make sure we can find the categories
        available_cat=colnames(ordered_normalized);
        missing_cat=setdiff(additional_cat, available_cat);
        if(length(missing_cat)){
                cat("Error: Could not find categories: \n");
                print(missing_cat);
                quit(status=-1);
        }

        # :: Remove categories we have already extracted in the top N
        already_extracted_cat=colnames(top_cat);
        extra_cat=setdiff(additional_cat, already_extracted_cat);

        num_extra_to_extract=length(extra_cat);
        cat("Num Extra Categories to Extract: ", num_extra_to_extract, "\n");

        # Allocate/Prepare output matrix
        num_out_mat_cols=num_top_to_extract+num_extra_to_extract+1;
        out_mat=matrix(0, nrow=num_samples, ncol=num_out_mat_cols);
        rownames(out_mat)=rownames(ordered_normalized);
        colnames(out_mat)=c(already_extracted_cat, extra_cat, "Remaining");

        # Copy over top and additional categories, and compute remainding
        all_cat_names=c(already_extracted_cat, extra_cat);
        out_mat[,all_cat_names]=ordered_normalized[,all_cat_names];
        out_mat[,"Remaining"]=apply(out_mat, 1, function(x){1-sum(x)});

        return(out_mat);

}

additive_log_ratio=function(ordered_matrix){
# Assumes last column will be the denominator

	num_cat=ncol(ordered_matrix);
	num_samp=nrow(ordered_matrix);

	denominator=ordered_matrix[,num_cat];
	alr_mat=matrix(0, nrow=num_samp, ncol=(num_cat-1));
	
	for(i in 1:num_samp){
		alr_mat[i,]=log(ordered_matrix[i,1:(num_cat-1)]/denominator[i]);
	}

	rownames(alr_mat)=rownames(ordered_matrix)
	colnames(alr_mat)=head(colnames(ordered_matrix), num_cat-1);

	alr_struct=list();
	alr_struct[["transformed"]]=alr_mat;
	alr_struct[["denominator"]]=denominator;

	return(alr_struct);
}

plot_text=function(strings){
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
}

title_page=function(title, subtitle=""){

	par(mfrow=c(1,1));
	plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n",  xaxt="n", yaxt="n",
		xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
		);

	if(subtitle!=""){
		text(.5, .5, title, adj=c(.5,-1), cex=4, font=2); 
		text(.5, .5, subtitle, adj=c(.5, 1), cex=2); 
	}else{
		text(.5, .5, title, adj=c(.5,.5), cex=4, font=2); 
	}
}

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, deci_pts=4, 
	label_zeros=T, counts=F, value.cex=2, 
	plot_col_dendr=F,
	plot_row_dendr=F
){

        num_row=nrow(mat);
        num_col=ncol(mat);

	row_names=rownames(mat);
	col_names=colnames(mat);

	orig.par=par(no.readonly=T);

	cat("Painting Matrix: ", title, "\n");
        cat("Num Rows: ", num_row, "\n");
        cat("Num Cols: ", num_col, "\n");


	if(num_row==0 || num_col==0){
		plot(0, type="n", xlim=c(-1,1), ylim=c(-1,1), xaxt="n", yaxt="n", bty="n", xlab="", ylab="",
			main=title);
		text(0,0, "No data to plot...");
		return();
	}

	# Flips the rows, so becuase origin is bottom left
        mat=mat[rev(1:num_row),, drop=F];

	# Generate a column scheme
        num_colors=50;
        color_arr=rainbow(num_colors, start=0, end=4/6);
        if(high_is_hot){
                color_arr=rev(color_arr);
        }

	# Provide a means to map values to an (color) index 
        remap=function(in_val, in_range, out_range){
                in_prop=(in_val-in_range[1])/(in_range[2]-in_range[1])
                out_val=in_prop*(out_range[2]-out_range[1])+out_range[1];
                return(out_val);
        }

	# If range is not specified, find it based on the data
        if(is.na(plot_min)){
                plot_min=min(mat);
        }
        if(is.na(plot_max)){
                plot_max=max(mat);
        }
        cat("Plot min/max: ", plot_min, "/", plot_max, "\n");

	# Get Label lengths
	row_max_nchar=max(nchar(row_names));
	col_max_nchar=max(nchar(col_names));
	cat("Max Row Names Length: ", row_max_nchar, "\n");
	cat("Max Col Names Length: ", col_max_nchar, "\n");

	##################################################################################################
	
	get_dendrogram=function(in_mat, type){
		if(type=="row"){
			dendist=dist(in_mat);
		}else{
			dendist=dist(t(in_mat));
		}
		
		get_clstrd_leaf_names=function(den){
		# Get a list of the leaf names, from left to right
			den_info=attributes(den);
			if(!is.null(den_info$leaf) && den_info$leaf==T){
				return(den_info$label);
			}else{
				lf_names=character();
				for(i in 1:2){
					lf_names=c(lf_names, get_clstrd_leaf_names(den[[i]]));
				}
				return(lf_names);
			}
		}

		hcl=hclust(dendist, method="ward.D2");
		dend=list();
		dend[["tree"]]=as.dendrogram(hcl);
		dend[["names"]]=get_clstrd_leaf_names(dend[["tree"]]);
		return(dend);
	}


	##################################################################################################
	# Comput Layouts
	col_dend_height=ceiling(num_row*.1);
	row_dend_width=ceiling(num_col*.2);
	
	heatmap_height=num_row;
	heatmap_width=num_col;

	if(num_row==1){
		plot_row_dendr=F;
	}
	if(num_col==1){
		plot_col_dendr=F;
	}

	# Don't plot dendrogram if there are any NAs in the matrix
	if(any(is.na(mat))){
		plot_col_dendr=F;
		plot_row_dendr=F;
	}

	if(plot_col_dendr && plot_row_dendr){
		layoutmat=matrix(
			c(
			rep(c(rep(4, row_dend_width), rep(3, heatmap_width)), col_dend_height),
			rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height)
			), byrow=T, ncol=row_dend_width+heatmap_width);

		col_dendr=get_dendrogram(mat, type="col");
		row_dendr=get_dendrogram(mat, type="row");

		mat=mat[row_dendr[["names"]], col_dendr[["names"]], drop=F];
		
	}else if(plot_col_dendr){
		layoutmat=matrix(
			c(
			rep(rep(2, heatmap_width), col_dend_height),
			rep(rep(1, heatmap_width), heatmap_height)
			), byrow=T, ncol=heatmap_width); 

		col_dendr=get_dendrogram(mat, type="col");
		mat=mat[, col_dendr[["names"]], drop=F];
		
	}else if(plot_row_dendr){
		layoutmat=matrix(
			rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height),
			byrow=T, ncol=row_dend_width+heatmap_width);

		row_dendr=get_dendrogram(mat, type="row");
		mat=mat[row_dendr[["names"]],,drop=F];
	}else{
		layoutmat=matrix(
			rep(1, heatmap_height*heatmap_width), 
			byrow=T, ncol=heatmap_width);
	}

	#print(layoutmat);
	layout(layoutmat);

	##################################################################################################
	
	par(oma=c(col_max_nchar*.60, 0, 3, row_max_nchar*.60));
	par(mar=c(0,0,0,0));
        plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), xaxt="n", yaxt="n", bty="n", xlab="", ylab="");
	mtext(title, side=3, line=0, outer=T, font=2);

        # x-axis
        axis(side=1, at=seq(.5, num_col-.5, 1), labels=colnames(mat), las=2, line=-1.75);
        axis(side=4, at=seq(.5, num_row-.5, 1), labels=rownames(mat), las=2, line=-1.75);

        if(log_col){
                plot_min=log10(plot_min+.0125);
                plot_max=log10(plot_max+.0125);
        }

        for(x in 1:num_col){
                for(y in 1:num_row){

                        if(log_col){
                                col_val=log10(mat[y,x]+.0125);
                        }else{
                                col_val=mat[y,x];
                        }

                        remap_val=remap(col_val, c(plot_min, plot_max), c(1, num_colors));
                        col_ix=ceiling(remap_val);

                        rect(x-1, y-1, (x-1)+1, (y-1)+1, border=NA, col=color_arr[col_ix]);

                        if(mat[y,x]!=0 || label_zeros){
                                if(counts){
                                        text_lab=sprintf("%i", mat[y,x]);
                                }else{
                                        text_lab=sprintf(paste("%0.", deci_pts, "f", sep=""), mat[y,x]);
                                }
                                text(x-.5, y-.5, text_lab, srt=atan(num_col/num_row)/pi*180, cex=value.cex, font=2);
                        }
                }
        }

	##################################################################################################

	par(mar=c(0, 0, 0, 0));

	if(plot_row_dendr && plot_col_dendr){
		rdh=attributes(row_dendr[["tree"]])$height;
		cdh=attributes(col_dendr[["tree"]])$height;
		plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", bty="n", xlim=c(rdh, 0));
		plot(col_dendr[["tree"]], leaflab="none",xaxt="n", yaxt="n", bty="n", ylim=c(0, cdh));
		plot(0,0, type="n", bty="n", xaxt="n", yaxt="n");
		#text(0,0, "Placeholder");
	}else if(plot_row_dendr){
		rdh=attributes(row_dendr[["tree"]])$height;
		plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", bty="n", xlim=c(rdh, 0));
		#text(0,0, "Row Dendrogram");
	}else if(plot_col_dendr){
		cdh=attributes(col_dendr[["tree"]])$height;
		plot(col_dendr[["tree"]], leaflab="none", xaxt="n", yaxt="n", bty="n", ylim=c(0, cdh));
		#text(0,0, "Col Dendrogram");
	}

	par(orig.par);

}


plot_histograms=function(table){
	num_cols=ncol(table);	
	orig.par=par(no.readonly=T);

	par(mfrow=c(5,3));
	par(mar=c(2,2,2,2));
	par(oma=c(2,2,2,2));
	colname=colnames(table);
	for(i in 1:num_cols){
		vals=table[,i];
		if(mode(vals)!="numeric" || is.factor(vals)){
			vals=as.factor(vals);
			barplot(prop.table(table(vals)), main=colname[i], col="white");
		}else{
			hist(vals, main=colname[i], xlab="values", ylab="");
		}
	}

	par(orig.par);
}

add_sign_col=function(coeff){

	if(nrow(coeff)==0){
		return(coeff);
	}else{

		cnames=colnames(coeff);
		pval_ix=which(cnames=="Pr(>|t|)");
		pval=coeff[,pval_ix];

		sig_char=function(val){
			if(!is.null(val) && !is.nan(val) && !is.na(val)){
				if(val <= .0001){ return("***");}
				if(val <= .001 ){ return("** ");}
				if(val <= .01  ){ return("*  ");}
				if(val <= .05  ){ return(":  ");}
				if(val <= .1   ){ return(".  ");}
				return(" ");
			}else{
				return(" ");
			}
		}

		sig_arr=sapply(pval, sig_char);
		fdr=round(p.adjust(pval, method="fdr"), 4);
		fdr_sig_char=sapply(fdr, sig_char);
		#out_mat=cbind(coeff, fdr);

		fmt=function(x){
			return(formatC(x, format="f", digits=4, width=7));
		}

		out_mat=cbind(
			fmt(coeff), sig_arr, fmt(fdr), fdr_sig_char);
		colnames(out_mat)=c(cnames, "Signf", "FDR", "Signf");
		
		return(out_mat);
	}	
}

sig_char=function(val){
	if(!is.null(val) && !is.nan(val) && !is.na(val)){
		if(val <= .0001){ return("***");}
		if(val <= .001 ){ return("** ");}
		if(val <= .01  ){ return("*  ");}
		if(val <= .05  ){ return(":  ");}
		if(val <= .1   ){ return(".  ");}
	}
	return(" ");
}

add_signf_char_to_matrix=function(mat){
	numrows=nrow(mat);
	numcols=ncol(mat);

	outmat=matrix(sprintf("%0.4f",mat), ncol=numcols);
	colnames(outmat)=colnames(mat);
	rownames(outmat)=rownames(mat);
	outmat=cbind(rbind(outmat, "", ""),"", "");

	sigfun=function(x){
		mx=min(x, na.rm=T);
		return(sig_char(mx));
	}

	sigrow=apply(mat, 1, sigfun);
	sigcol=apply(mat, 2, sigfun);

	outmat[1:numrows, numcols+2]=sigrow;
	outmat[numrows+2, 1:numcols]=sigcol;

	return(outmat);

}

plot_fit=function(fit, sumfit, i=1){
	par.orig=par(no.readonly=T);

	observed=as.numeric(as.matrix(fit$model)[,i, drop=F]);
	predicted=as.numeric(as.matrix(fit$fitted.values)[,i, drop=F]);

	adjsqrd=signif(sumfit[[i]]$adj.r.squared, 4);
	fstat=sumfit[[i]]$fstatistic;
	pval=signif(1-pf(fstat["value"], fstat["numdf"], fstat["dendf"]), 4);

	name=colnames(predicted);

	par(mar=c(15, 15, 15, 15));

	obs_pred_range=range(c(observed, predicted));
	span=diff(obs_pred_range);
	pad=span*.07;

	plot(observed, predicted, 
		xlim=c(obs_pred_range[1]-pad, obs_pred_range[2]+pad),
		ylim=c(obs_pred_range[1]-pad, obs_pred_range[2]+2*pad),
		main="", xlab="Observed", ylab="Predicted");

	mtext(name, line=5.5, font=2, cex=3);
	mtext("Predicted vs. Observed", line=4, font=2, cex=1.1);
	mtext(paste("Adjusted R^2 = ", adjsqrd, sep=""),line=2, cex=.9);
	mtext(paste("Model F-stat p-value = ", pval, sep=""), line=1, cex=.9);

	abline(a=0, b=1, col="blue");
	points(lowess(observed, predicted), type="l", col="red");

	
	legend(obs_pred_range[1]+pad, obs_pred_range[2]+pad, legend=c("Model", "Ideal"), 
		fill=c("red", "blue"), bty="n");

	par(par.orig);
}

plot_ts_stat_table=function(stat_mat, 
	title="", subtitle="",
	grp_colors=NULL, 
	plot_tmax=NULL, plot_tmin=NULL,
	plot_ymax=NULL, plot_ymin=NULL,
	nlog10_reflines=F,
	zero_refline=F
	){

	#print(stat_mat);
	num_times=ncol(stat_mat);
	num_rows=nrow(stat_mat);

	grp_names=rownames(stat_mat);
	time_names=colnames(stat_mat)
	time_values=as.numeric(time_names);

	cat("Times:\n");
	print(time_values);

	stat_range=range(stat_mat, na.rm=T);
	time_range=range(time_values, na.rm=T);


	if(is.null(plot_ymax)){
		plot_ymax=stat_range[2];
	}
	if(is.null(plot_ymin)){
		plot_ymin=stat_range[1];
	}

	if(num_rows==0){
		plot_ymax=0;
		plot_ymin=0;
	}

	if(is.null(plot_tmax)){
		plot_xmax=time_range[2];
	}
	if(is.null(plot_tmin)){
		plot_xmin=time_range[1];
	}

	if(nlog10_reflines){
		ref_lines=-log10(c(.1,.05,.01));
		plot_ymax=max(plot_ymax, ref_lines);
	}


	if(is.null(grp_colors)){
		grp_colors=rainbow(num_rows, start=0, end=5/6, alpha=0.33);
		names(grp_colors)=grp_names;
		print(grp_colors);
	}

	if(nlog10_reflines){
		par(mar=c(4, 4, 4, 4));
	}else{
		par(mar=c(4, 4, 4, 1));
	}

	plot(0, type="n", xlab="Time", ylab="", 
		xlim=c(plot_xmin, plot_xmax), ylim=c(plot_ymin, plot_ymax));

	title(main=title);
	title(main=subtitle, line=.8, cex.main=.85);

	if(nlog10_reflines){
		abline(h=ref_lines, lty=c(3,2,5), col="grey");
		axis(side=4, at=ref_lines, labels=c("0.10", "0.05", "0.01"), las=2);
	}

	if(zero_refline){
		abline(h=0, col="grey");
	}


	if(num_rows>0){
		for(i in 1:num_rows){
			points(time_values, stat_mat[i,], type="b", col=grp_colors[i], lwd=4);
		}
	}

}

plot_group_legend=function(color_map){
	par(mar=c(0,0,0,0));
	plot(0, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), bty="n", xaxt="n", yaxt="n");

	num_entries=length(color_map);
	leg_nam=names(color_map);

	if(num_entries>12){
		cutoff=floor(seq(1,num_entries, length.out=4));

		c1r=cutoff[1]:cutoff[2];
		c2r=(cutoff[2]+1):cutoff[3];
		c3r=(cutoff[3]+1):cutoff[4];
		
		legend(0,1, legend=leg_nam[c1r], fill=color_map[c1r], bty="n");
		legend(.33,1, legend=leg_nam[c2r], fill=color_map[c2r], bty="n");
		legend(.66,1, legend=leg_nam[c3r], fill=color_map[c3r], bty="n");
	}else{
		if(num_entries>0){
			legend(0,1, legend=leg_nam, fill=color_map, bty="n");
		}
	}
}

##############################################################################
##############################################################################

pdf(paste(OutputRoot, ".multn_resp_alr_pred.pdf", sep=""), height=11, width=9.5);

# Load summary file table counts 
cat("Loading summary table...\n");
counts=load_summary_file(SummaryFile);

# Remove zero count samples
tot=apply(counts, 1, sum);
nonzero=tot>0;
if(!(all(nonzero))){
	cat("WARNING: Zero count samples found:\n");
	samp_names=rownames(counts);
	print(samp_names[!nonzero]);
	cat("\n");
	counts=counts[nonzero,,drop=F];
}

num_taxa=ncol(counts);
num_samples=nrow(counts);
#print(counts);

# Shorten cateogry names
if(ShortenCategoryNames!=""){
	full_names=colnames(counts);
	splits=strsplit(full_names, ShortenCategoryNames);
	short_names=character();
	for(i in 1:length(full_names)){
		short_names[i]=tail(splits[[i]], 1);

		short_names[i]=gsub("_unclassified$", "_uncl", short_names[i]);
		short_names[i]=gsub("_group", "_grp", short_names[i]);
                short_names[i]=gsub("\\[", "", short_names[i]);
                short_names[i]=gsub("\\]", "", short_names[i]);
                short_names[i]=gsub("\\(", "", short_names[i]);
                short_names[i]=gsub("\\)", "", short_names[i]);
	}
	colnames(counts)=short_names;
	cat("Names have been shortened.\n");
}

# Normalize
counts=counts+.5;
normalized=normalize(counts);
#print(normalized);

if(UseRemaining){
	category_names=colnames(counts);	
	uc_cat_names=toupper(category_names);
	remaining_ix=which(uc_cat_names=="REMAINDER" | uc_cat_names=="REMAINING");
	if(length(remaining_ix)!=1){
		cat("*******************************************************\n");
		cat("*  WARNING:  Could not identify remaining column.     *\n");
		cat("*******************************************************\n");
		UseRemaining=F;
	}else{
		cat("Remaining original column: ", remaining_ix, "\n");
		# Take out "remaining" column so it doesn't end up as a top column
		normalized_remaining_col_dat=normalized[,remaining_ix, drop=F];
		normalized=normalized[,-remaining_ix];
	}
}

# Reorder by abundance
cat("Reordering summary table categories by abundance...\n");
mean_abund=apply(normalized, 2, mean);
ix=order(mean_abund, decreasing=TRUE);
normalized=normalized[,ix];
mean_abund=mean_abund[ix];

if(UseRemaining){
	normalized=cbind(normalized, normalized_remaining_col_dat);
	mean_abund=c(mean_abund, mean(normalized_remaining_col_dat));
}

sorted_taxa_names=colnames(normalized);

num_top_taxa=NumVariables;
num_top_taxa=min(c(num_top_taxa, num_taxa));
prop_abundance_represented=sum(mean_abund[1:num_top_taxa]);

cat("\nThe top ", num_top_taxa, " taxa are:\n", sep="");
for(i in 1:num_top_taxa){
	cat("\t", sorted_taxa_names[i], "\t[", mean_abund[i], "]\n", sep="");
}
cat("\n");

cat("Accounting for ", prop_abundance_represented, " of taxa.\n", sep="");
cat("\n");

##############################################################################

# Load factors
cat("Loading Factors...\n");
factors=load_factors(FactorsFile);
factor_names=colnames(factors);
num_factors=ncol(factors);
factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);

cat("\n");
cat(num_factors, " Factor(s) Loaded:\n", sep="");
print(factor_names);
cat("\n");

# Load covariate/response list
covariates_arr=load_list(CovariatesFile);
cat("Covariates Variables:\n");
print(covariates_arr);
cat("\n");

responses_arr=ResponseColname;

# Confirm we can find the response column
if(length(intersect(ResponseColname, factor_names))){
	responses=factors[,ResponseColname];
	resp_tab=table(responses);
	cat("Response Column Name: ", ResponseColname, " found in factor file.\n", sep="");
	print(resp_tab);
}else{
	cat("Error: ", ResponseColname, " not found in factor file.\n", sep="");
	quit(-1);
}


required_arr=NULL;
if(""!=RequiredFile){
	required_arr=load_list(RequiredFile);
	cat("Required Variables:\n");
	print(required_arr);
	cat("\n");
}else{
	cat("No Required Variables specified...\n");
}

plot_text(c(
	paste("     Summary File: ", SummaryFile, sep=""),
	paste("     Factors File: ", FactorsFile, sep=""),
	paste("  Covariates File: ", CovariatesFile, sep=""),
	paste(" Response Colname: ", ResponseColname, sep=""),
	paste("    Required File: ", RequiredFile, sep=""),
	paste("  Additional File: ", AdditionalCatFile, sep=""),
	"\n",
	paste("Output File: ", OutputRoot, sep=""),
	"\n",
	paste("Number of MALR Variables: ", NumVariables, sep=""),
	paste("Reference Levels File: ", ReferenceLevelsFile, sep=""),
	paste("Use Remaining? ", UseRemaining, "\n"),
	paste("Shorten Category Names: ", ShortenCategoryNames, sep="")
));
	
plot_text(c(
	"Variables Targeted:",
	"",
	"Responses:",
	capture.output(print(resp_tab)),
	"",
	"Covariates:",
	capture.output(print(covariates_arr)),
	"",
	"Required Variables:",
	 capture.output(print(required_arr))
));


overlapping_variables=intersect(responses_arr, covariates_arr);
if(length(overlapping_variables)!=0){
	cat("Error:  You have the same variable names in response and covariate:\n");
	print(overlapping_variables);	
	quit(-1);
}

kept_variables=unique(c(responses_arr, covariates_arr, TimeColumn, SubjectColumn));
kept_factors=factors[,kept_variables, drop=F];

cat("Kept Variables out of factors:\n");
print(kept_variables);

# Relevel factor levels
if(ReferenceLevelsFile!=""){
        ref_lev_mat=load_reference_levels_file(ReferenceLevelsFile)
	print(ref_lev_mat);
        kept_factors=relevel_factors(kept_factors, ref_lev_mat);
	cat("Releveling factors levels to specified reference.\n");
	
}else{
        cat("No Reference Levels File specified.\n");
}

response_factors=kept_factors[,ResponseColname, drop=F];
summary(response_factors);
covariate_factors=kept_factors[,covariates_arr, drop=F];
subject_ids=kept_factors[,SubjectColumn, drop=F];
time_ids=kept_factors[,TimeColumn, drop=F];

resp_cov_factors=cbind(subject_ids, time_ids, response_factors, covariate_factors);


##############################################################################
# Reconcile factors with samples
factor_sample_ids=rownames(resp_cov_factors);
counts_sample_ids=rownames(counts);

#print(factor_sample_id);
#print(counts_sample_id);

shared_sample_ids=intersect(factor_sample_ids, counts_sample_ids);
num_shared_sample_ids=length(shared_sample_ids);
num_factor_sample_ids=length(factor_sample_ids);
num_counts_sample_ids=length(counts_sample_ids);

cat("Num counts sample IDs: ", num_counts_sample_ids, "\n");
cat("Num factor sample IDs: ", num_factor_sample_ids, "\n");
cat("Num shared sample IDs: ", num_shared_sample_ids, "\n");
cat("\n");

cat("Samples missing from summary table/count information:\n");
print(setdiff(factor_sample_ids, counts_sample_ids));
cat("\n");
cat("Samples missing from factor information:\n");
print(setdiff(counts_sample_ids, factor_sample_ids));
cat("\n");
cat("Total samples shared: ", num_shared_sample_ids, "\n");

shared_sample_ids=sort(shared_sample_ids);

# Reorder data by sample id
normalized=normalized[shared_sample_ids,];
num_samples=nrow(normalized);
recon_factors=resp_cov_factors[shared_sample_ids,,drop=F];

#factors_wo_nas=remove_sample_or_factors_wNA(recon_factors);
num_samples_recon=nrow(recon_factors);
num_factors_recon=ncol(recon_factors);

factors_wo_nas_res=remove_sample_or_factors_wNA_parallel(recon_factors, required=required_arr, num_trials=640000, num_cores=64, outfile=paste(OutputRoot, ".noNAs", sep=""));

factors_wo_nas=factors_wo_nas_res$factors;
factor_names_wo_nas=colnames(factors_wo_nas);
factor_sample_ids=rownames(factors_wo_nas);
normalized=normalized[factor_sample_ids,];
num_samples_wo_nas=nrow(factors_wo_nas);
num_factors_wo_nas=ncol(factors_wo_nas);

responses_arr=intersect(responses_arr, factor_names_wo_nas);
covariates_arr=intersect(covariates_arr, factor_names_wo_nas);

if(length(responses_arr)==0){
	cat("Error: None of your response variables survived NA remove.  You need to make some of them required.\n");
	quit(status=-1);
}
if(length(covariates_arr)==0){
	cat("Warning: None of your covariates survived NA remove.  You need to make some of them required.\n");
}

response_factors=factors_wo_nas[,responses_arr, drop=F];
covariate_factors=factors_wo_nas[,covariates_arr, drop=F];

##############################################################################

if(num_top_taxa>= num_taxa){
	num_top_taxa = (num_taxa-1);
	cat("Number of taxa to work on was changed to: ", num_top_taxa, "\n");
}

##############################################################################

if(AdditionalCatFile!=""){
	additional_categories=load_list(AdditionalCatFile);
}else{
	additional_categories=c();
}

cat_abundances=extract_top_categories(normalized, num_top_taxa, additional_cat=additional_categories);
resp_alr_struct=additive_log_ratio(cat_abundances);
alr_categories_val=resp_alr_struct$transformed;
alr_cat_names=colnames(alr_categories_val);

plot_text(c(
	"Acceptable Variables after NA Removal:",
	"",
	"Responses:",
	capture.output(print(responses_arr)),
	"",
	"Covariates:",
	capture.output(print(covariates_arr)),
	"",
	"",
	paste("Num Reconciled Samples: ", num_samples_recon, sep=""),
	paste("Num Reconciled Factors: ", num_factors_recon, sep=""),
	"",
	paste("Num Samples w/o NAs: ", num_samples_wo_nas, sep=""),
	paste("Num Factors w/o NAs: ", num_factors_wo_nas, sep="")
));

plot_text(c(
	paste("ALR Categories (Top ", num_top_taxa, ")", sep=""),
	capture.output(print(alr_cat_names))
));

#response_factors=as.matrix(response_factors);
dafr_predictors_factors=as.data.frame(cbind(covariate_factors, alr_categories_val));

cat("Response Summary:\n");
s=summary(response_factors);
plot_text(c(
	"Response Summary:",
	"\n",
	capture.output(print(s))
));

cat("Plotting Response Histograms:\n");
#print(response_factors);
plot_histograms(response_factors);

cat("\n");
cat("Covariate Summary:\n");
s=summary(covariate_factors);
plot_text(c(
	"Covariates Summary:",
	"\n",
	capture.output(print(s))
));
if(length(covariates_arr)>0){
	cat("Plotting Covariate Histograms:\n");
	print(covariate_factors);
	plot_histograms(covariate_factors);
}

cat("\n");
cat("ALR Category Summary:\n");
s=summary(alr_categories_val);
plot_text(c(
	"ALR Categories Summary:",
	"\n",
	capture.output(print(s))
));
cat("Plotting ALR Category Histograms:\n");
#print(alr_categories_val);
plot_histograms(alr_categories_val);

###############################################################################

# Remove variables with no information
subject_ids=character();
time_ids=numeric();


cat("Setting up subject id variable...\n");
if(length(SubjectColumn)){
	subject_ids=factors_wo_nas[, SubjectColumn];
}else{
	cat("Subject IDs not specified, using sample IDs.\n");
	subject_ids=rownames(factors_wo_nas);
}	

cat("Setting up times series variable...\n");
if(length(TimeColumn)){
	time_ids=factors_wo_nas[, TimeColumn];
}else{
	cat("Time IDs not specified, assuming single cross-section.\n");
	time_ids=rep(0, length(subject_ids));
}

cat("\n");
cat("All Time IDs:\n");
print(time_ids);
cat("All Subject IDs:\n");
print(subject_ids);

unique_time_ids=sort(unique(time_ids));
unique_subject_ids=sort(unique(subject_ids));
unique_responses=sort(unique(factors_wo_nas[,ResponseColname]));

cat("\n");
cat("Unique Time Points:\n");
print(unique_time_ids);
cat("Unique Subject IDs:\n");
print(unique_subject_ids);
cat("\n");
cat("Unique Response Groups:\n");
print(unique_responses);

###############################################################################

require(nnet);
require(generalhoslem);

process_model=function(fit, resp){
	
	res=list();
	res[["fit"]]=fit;
	res[["summary"]]=summary(fit);

	coeff=res[["summary"]]$coefficients
	stderr=res[["summary"]]$standard.errors;
	z=coeff/stderr;
	res[["pvalues"]]=(1-pnorm(abs(z), 0, 1))*2;

	# See M. W. Fagerland and D. W. Hosmer,
	# "A generalized Hosmer–Lemeshowgoodness-of-fit test for multinomial logisticregression models"
	# Null distribution is if there was a good fit, so low p-values are bad fits.
	#res[["gof"]]=logitgof(resp, fitted(fit));
	res[["NegLogLikelihood"]]=res[["summary"]]$value;
	
	return(res);

}

fit_info=list();

for(cur_time_id in unique_time_ids){

	cat("\n*******************************************************************\n");
	cat("Working on time: ", cur_time_id, "\n");

	cur_time_ix=(cur_time_id==time_ids);
	cur_time_str=sprintf("%02g", cur_time_id);

	cur_factors=factors_wo_nas[cur_time_ix,,drop=F];
	cur_alr=as.data.frame(alr_categories_val[cur_time_ix,,drop=F]);

	num_samples_at_curtime=nrow(cur_factors);
	cat("Num Samples: ", num_samples_at_curtime, "\n", sep="");
	cat("\n");

	cur_responses=cur_factors[,ResponseColname];
	cur_predictors=cur_factors[,covariates_arr];

	fit_info[[cur_time_str]]=list();

	# Fit intercept only
	cat("Intercept Only:\n");
	null_model_str=paste("cur_responses ~ 1");
	null_mlr_fit=multinom(as.formula(null_model_str), data=cur_factors);
	fit_info[[cur_time_str]][["intercept_only"]]=process_model(null_mlr_fit, cur_responses);
	cat("\n");

	# Fit covariates
	cat("Covariates Only:\n");
	cov_model_str=paste("cur_responses ~ ", paste(covariates_arr, collapse=" + ", sep=""));
	cov_mlr_fit=multinom(as.formula(cov_model_str), data=cur_factors);
	fit_info[[cur_time_str]][["covariates_only"]]=process_model(cov_mlr_fit, cur_responses);
	cat("\n");

	# Fit alr categories
	cat("ALR Categories Only:\n");
	alr_model_str=paste("cur_responses ~ ", paste(alr_cat_names, collapse=" + ", sep=""));
	alr_mlr_fit=multinom(as.formula(alr_model_str), data=cur_alr);
	fit_info[[cur_time_str]][["alr_only"]]=process_model(alr_mlr_fit, cur_responses);
	cat("\n");

	# Fit combined 
	cat("Full Combined:\n");
	comb_model_str=paste("cur_responses ~ ", paste(c(covariates_arr, alr_cat_names), collapse=" + ", sep=""));
	comb_mlr_fit=multinom(as.formula(comb_model_str), data=cbind(cur_factors, cur_alr));
	fit_info[[cur_time_str]][["alr_and_covariates"]]=process_model(comb_mlr_fit, cur_responses);
	cat("\n");

	# Overall statistics
	fit_info[[cur_time_str]][["responses"]]=table(as.character(cur_responses));
	fit_info[[cur_time_str]][["num_responses"]]=length(cur_responses);

	#print(fit_info[[cur_time_str]][["combined"]][["pvalues"]]);

}

# Fit variables:
# [1] "n"             "nunits"        "nconn"         "conn"          "nsunits"       "decay"
# [7] "entropy"       "softmax"       "censored"      "value"         "wts"           "convergence"
#[13] "fitted.values" "residuals"     "lev"           "call"          "terms"         "weights"
#[19] "deviance"      "rank"          "lab"           "coefnames"     "vcoefnames"    "xlevels"
#[25] "edf"           "AIC"


# Extract out key statistics into simpler data structures
model_types=c("covariates_only", "alr_only", "alr_and_covariates");
#model_types=c("intercept_only", "covariates_only", "alr_only", "alr_and_covariates");
time_str_ids=sprintf("%02g", unique_time_ids);

num_time_pts=length(unique_time_ids);
num_modeltypes=length(model_types);
num_unique_responses=length(unique_responses);

# Assign group colors
grp_colors=rainbow(num_unique_responses, start=0, end=5/6, alpha=3/8);
names(grp_colors)=unique_responses;

# Assign model colors
model_colors=rainbow(num_modeltypes, start=0+1/12, end=5/6+1/12, alpha=3/8);
names(model_colors)=model_types;


# Define generic matrix, so we can copy it empty for different variables
stat_matrix=matrix(NA, nrow=num_modeltypes, ncol=num_time_pts);
stat_arr=rep(NA, num_time_pts);

colnames(stat_matrix)=time_str_ids;
rownames(stat_matrix)=model_types;
names(stat_arr)=time_str_ids;

aic_matrix=stat_matrix;
sampsize_arr=stat_arr;

resp_grps=matrix(NA, ncol=num_time_pts, nrow=num_unique_responses);
colnames(resp_grps)=time_str_ids;
rownames(resp_grps)=unique_responses;

pvalues=list();
coefficients=list();

for(cur_time_str  in time_str_ids){

	cat("Extracting Time: ", cur_time_str, "\n");

	sampsize_arr[cur_time_str]=fit_info[[cur_time_str]][["num_responses"]];

	resp_grps[unique_responses, cur_time_str]=fit_info[[cur_time_str]][["responses"]][unique_responses];

	pvalues[[cur_time_str]]=list();
	coefficients[[cur_time_str]]=list();
	for(modix in model_types){
		cat("Extracting: ", modix, "\n");

		aic_matrix[modix, cur_time_str]=fit_info[[cur_time_str]][[modix]][["fit"]][["AIC"]];

		pvalues[[cur_time_str]][[modix]]=fit_info[[cur_time_str]][[modix]][["pvalues"]];
		coefficients[[cur_time_str]][[modix]]=fit_info[[cur_time_str]][[modix]][["Coefficients"]]=
			summary(fit_info[[cur_time_str]][[modix]][["fit"]])$coefficients;

	}

}

cat("\n");
cat("Sample Sizes:\n");
print(sampsize_arr);

cat("\n");
cat("Response Group Sizes:\n");
print(resp_grps);

layout_m=matrix(c(1,1,1,1,1,2), nrow=6, ncol=1);
layout(layout_m);
plot_ts_stat_table(resp_grps, title="Response Group Sizes", subtitle="Number of Samples per Group Over Time", grp_colors=grp_colors);
plot_group_legend(grp_colors);

# Plot model fits

cat("AIC:\n");
print(aic_matrix);
layout_m=matrix(c(1,1,1,1,1,2), nrow=6, ncol=1);
layout(layout_m);
plot_ts_stat_table(aic_matrix, title="AIC", subtitle="Lower Values, Better Fit", grp_colors=model_colors);
plot_group_legend(model_colors);


# Plot all coefficients 
# Plot all pvalues
# Plot coefficients at .1, .05, and .01

cat("----------------------------------------------------------------------------\n");

#print(coefficients);

first_time_pt=time_str_ids[1];
cat("First Time Point:", first_time_pt, "\n");

cat("Excluding Reference: ", as.character(unique_responses)[1], "\n\n");
response_no_reference=unique_responses[-1];

coef_bytime_list=list();
pval_bytime_list=list();


# for each response category
for(resp_ix in response_no_reference){
	cat("Extracting Response: ", resp_ix, "\n");
	# for each 4 model types

	coef_bytime_list[[resp_ix]]=list();
	pval_bytime_list[[resp_ix]]=list();

	for(model_ix in model_types){

		cat("Extracting Model: ", model_ix, "\n");

		cur_tab=coefficients[[first_time_pt]][[model_ix]];
		
		inc_var_names=colnames(cur_tab);
		num_ivn=length(inc_var_names);

		coef_by_time_matrix=matrix(NA, nrow=num_ivn, ncol=num_time_pts);
		rownames(coef_by_time_matrix)=inc_var_names;
		colnames(coef_by_time_matrix)=time_str_ids;

		pval_by_time_matrix=matrix(NA, nrow=num_ivn, ncol=num_time_pts);
		rownames(pval_by_time_matrix)=inc_var_names;
		colnames(pval_by_time_matrix)=time_str_ids;


		for(time_ix in time_str_ids){
			for(ivn_ix in inc_var_names){

				#cat("Transfering: ", time_ix, "/", ivn_ix, "/", resp_ix, "\n");

				cur_coef_tab=coefficients[[time_ix]][[model_ix]]
				responses_available=rownames(cur_coef_tab);

				if(any(resp_ix==responses_available)){
					coef_by_time_matrix[ivn_ix, time_ix]=
						coefficients[[time_ix]][[model_ix]][resp_ix, ivn_ix];

					pval_by_time_matrix[ivn_ix, time_ix]=
						pvalues[[time_ix]][[model_ix]][resp_ix, ivn_ix];
				}else{
					coef_by_time_matrix[ivn_ix, time_ix]=
						NA;
					pval_by_time_matrix[ivn_ix, time_ix]=
						NA;
				}
			}
		}

		coef_bytime_list[[resp_ix]][[model_ix]]=coef_by_time_matrix;
		pval_bytime_list[[resp_ix]][[model_ix]]=pval_by_time_matrix;
	}
}

neglog10_trans_mat=function(mat){
	minnonzero=min(mat[mat!=0], na.rm=T);		
	cat("Min Nonzero Pvalue: ", minnonzero, "\n");
	mat[mat==0]=minnonzero/10;
	nlog10pval=-log10(mat);
	return(nlog10pval);
}

for(resp_ix in response_no_reference){
		
	title_page(paste("Response:\n", resp_ix));

	for(model_ix in model_types){

		title_page(resp_ix, subtitle=paste("Model Type:", model_ix));

		cur_coef_by_time_matrix=coef_bytime_list[[resp_ix]][[model_ix]];
		cur_pval_by_time_matrix=pval_bytime_list[[resp_ix]][[model_ix]];

		num_model_coeff=nrow(cur_coef_by_time_matrix);

		var_colors=rainbow(num_model_coeff, start=0+1/12, end=5/6+1/12, alpha=3/8);
		names(var_colors)=rownames(cur_coef_by_time_matrix);

		#----------------------------------------------------------------------------
		# Plot coefficients
		par(mfrow=c(1,1));
		coef_tab=capture.output(print(cur_coef_by_time_matrix, digits=4));
		plot_text(c(
			paste("Reponse: ", resp_ix, "  Model: ", model_ix),
			"Coefficients:",
			"",
			coef_tab
		));

		layout_m=matrix(c(1,1,1,2), nrow=4, ncol=1);
		layout(layout_m);
		plot_ts_stat_table(cur_coef_by_time_matrix, title=paste(resp_ix, ", Coefficients: ", model_ix, "\n"),
			grp_colors=var_colors, zero_refline=T);
		plot_group_legend(var_colors);

		#----------------------------------------------------------------------------
		# Plot p-values
		par(mfrow=c(1,1));
		pval_signf_mat=add_signf_char_to_matrix(cur_pval_by_time_matrix);
		pval_tab=capture.output(print(pval_signf_mat, quote=F));
		plot_text(c(
			paste("Reponse: ", resp_ix, "  Model: ", model_ix),
			"P-Values:",
			"",
			pval_tab
		));

		nlog10pval=neglog10_trans_mat(cur_pval_by_time_matrix);

		layout_m=matrix(c(1,1,1,2), nrow=4, ncol=1);
		layout(layout_m);
		plot_ts_stat_table(nlog10pval, title=paste(resp_ix, ", -log10(P-Values): ", model_ix, "\n"),
			grp_colors=var_colors, plot_ymin=0, nlog10_reflines=T);
		plot_group_legend(var_colors);

		#----------------------------------------------------------------------------
		# Pull out predictors that were significant p<.1

		signf_coef_ix=apply(cur_pval_by_time_matrix, 1, function(x){ min(x, na.rm=T)<0.1 });
		print(signf_coef_ix);
		
		signf_coef=cur_coef_by_time_matrix[signf_coef_ix,,drop=F];
		signf_pval=cur_pval_by_time_matrix[signf_coef_ix,,drop=F];
		signf_colr=var_colors[signf_coef_ix];

		
		layout_m=matrix(c(1,1,1,2), nrow=4, ncol=1);
		layout(layout_m);
		plot_ts_stat_table(signf_coef, 
			title=paste(resp_ix, ", Significant Coefficients: ", model_ix, "\n"),
			grp_colors=signf_colr, zero_refline=T);
		plot_group_legend(signf_colr);
		
		layout_m=matrix(c(1,1,1,2), nrow=4, ncol=1);
		layout(layout_m);
		nlog10pval=neglog10_trans_mat(signf_pval);
		plot_ts_stat_table(nlog10pval, 
			title=paste(resp_ix, ", Significant -log10(P-Values): ", model_ix, "\n"),
			grp_colors=signf_colr, plot_ymin=0, nlog10_reflines=T);
		plot_group_legend(signf_colr);


	}
}

quit();

hereend








































covariate_formula_str=paste(covariates_arr, collapse=" + ");
alr_category_formula_str=paste(alr_cat_names, collapse=" + ");
cat("\n");
cat("Covariates: \n", covariate_formula_str, "\n", sep="");
cat("\n");
cat("ALR Predictors: \n", alr_category_formula_str, "\n", sep="");
cat("\n");

if(nchar(covariate_formula_str)){
	formula_str=paste("response_factors_dfr ~ ", covariate_formula_str, " + ", alr_category_formula_str, sep="");
	reduced_formula_str=paste("response_factors_dfr ~ ", covariate_formula_str, sep="");
}else{
	formula_str=paste("response_factors_dfr ~ ", alr_category_formula_str, sep="");
	reduced_formula_str=paste("response_factors_dfr ~ 1", sep="");
}

response_factors_dfr=as.matrix(response_factors);

cat("Fitting Full Model...\n");
lmfit=lm(as.formula(formula_str), data=dafr_predictors_factors);

cat("Fitting Reduced Model..\n");
reduced_lmfit=lm(as.formula(reduced_formula_str), data=dafr_predictors_factors);

lm_summaries=summary(lmfit);
reduced_lm_summaries=summary(reduced_lmfit);

num_responses=length(responses_arr);

manova_res=c();
if(num_responses > 1 && (num_responses < lmfit$df.residual)){
	tryCatch({
		manova_res=anova(lmfit);
	}, error=function(e){
		cat("ERROR: MANOVA Could not be perfomred.\n");
		print(e);
		cat("Things to do:\n");
		cat(" 1.) Check for high correlations amongst response variables.\n");
	});
		
}else{
	cat("WARNING: MANOVA Could not be performed because too many responses for number of samples...\n");
}

if(num_responses==1){

	#tmp=lmfit;
	#lmfit=list();
	#lmfit[[1]]=tmp;

	tmp=lm_summaries;
	lm_summaries=list();
	lm_summaries[[1]]=tmp;
	names(lm_summaries)=paste("Response", responses_arr[1]);

	tmp=reduced_lm_summaries;
	reduced_lm_summaries=list();
	reduced_lm_summaries[[1]]=tmp;
	names(reduced_lm_summaries)=paste("Response", responses_arr[1]);
	
}

#print(lm_summaries[[1]]$adj.r.squared);
#quit();
		
responses=names(lm_summaries);
reduced_responses=names(reduced_lm_summaries);

print(names(lm_summaries));
response_fit_names=gsub("^Response ", "", names(lm_summaries));

# Calculate reduced/full model pvalues
calc_model_imprv=function(mv_resp_name, reduced_form_str, full_form_str, pred_datafr, resp_datafr){
	#print(reduced_form_str);
	#print(full_form_str);
	#print(pred_datafr);
	#print(resp_datafr);

	reduced_form_str=gsub(mv_resp_name, "tmp_resp", reduced_form_str);
	full_form_str=gsub(mv_resp_name, "tmp_resp", full_form_str);

	num_resp=ncol(resp_datafr);
	resp_names=colnames(resp_datafr);

	improv_mat=matrix(NA, nrow=num_resp, ncol=7);
	rownames(improv_mat)=resp_names;
	colnames(improv_mat)=c("Full R^2", "Reduced R^2", "Full Adj R^2", "Reduced Adj R^2",
		"Diff Adj. R^2", "Perc Improvement", "Diff ANOVA P-Value"); 

	for(i in 1:num_resp){
		tmp_resp=resp_datafr[,i];
		all_data=cbind(tmp_resp, pred_datafr);
		full_fit=lm(as.formula(full_form_str), data=all_data);
		full_sum=summary(full_fit);
		reduced_fit=lm(as.formula(reduced_form_str), data=all_data);
		reduced_sum=summary(reduced_fit);

		anova_res=anova(reduced_fit, full_fit);
		anova_pval=anova_res["2", "Pr(>F)"];

		improv_mat[i,]=c(
			full_sum$r.squared,
			reduced_sum$r.squared,
			full_sum$adj.r.squared,
			reduced_sum$adj.r.squared,
			full_sum$adj.r.squared-reduced_sum$adj.r.squared,
			100*(full_sum$adj.r.squared-reduced_sum$adj.r.squared)/reduced_sum$adj.r.squared,
			anova_pval);

	}

	return(improv_mat);

}

improv_mat=calc_model_imprv("response_factors_dfr", reduced_formula_str, formula_str, 
	dafr_predictors_factors, response_factors_dfr);

###############################################################################
# Accumulate univariate results into matrix and generate heat map

num_coefficients=nrow(lm_summaries[[1]]$coefficients);
coefficient_names=rownames(lm_summaries[[1]]$coefficients);
summary_var_names=rownames(lm_summaries[[1]]$coefficients);
covariate_coefficients=setdiff(coefficient_names, c("(Intercept)", alr_cat_names));

# In case some alr categories are not calculable...
shrd_alr_names=intersect(alr_cat_names, coefficient_names);

summary_res_coeff=matrix(0, nrow=num_coefficients, ncol=num_responses, dimnames=list(coefficient_names, response_fit_names));
summary_res_pval=matrix(0, nrow=num_coefficients, ncol=num_responses, dimnames=list(coefficient_names, response_fit_names));
summary_res_rsqrd=matrix(0, nrow=num_responses, ncol=3, dimnames=list(response_fit_names, 
	c("Full Model", "Reduced Model", "Difference")));

for(i in 1:num_responses){
	cat("\n\nWorking on: ", responses[i], "\n");
	univar_summary=lm_summaries[[i]]$coefficients;

	rsquared=signif(c(lm_summaries[[i]]$r.squared, lm_summaries[[i]]$adj.r.squared),4);
	reduced_rsquared=signif(c(reduced_lm_summaries[[i]]$r.squared, reduced_lm_summaries[[i]]$adj.r.squared),4);
	
	fstat=lm_summaries[[i]]$fstatistic;
	pval=1-pf(fstat[["value"]], fstat[["numdf"]], fstat[["dendf"]]);

	rsqrd_diff=rsquared[1]-reduced_rsquared[1];
	adj_rsqrd_diff=signif(rsquared[2]-reduced_rsquared[2], 4);


	univar_summary_wsignf=add_sign_col(univar_summary[shrd_alr_names,,drop=F]);

	plot_text(c(
		paste("Univariate: ", responses[i], sep=""),
		"Covariates portion of Full Model (Covariates + ALR):",
		"",
		capture.output(
			print(univar_summary[covariate_coefficients,,drop=F], digits=4)
		)
	));

	plot_text(c(
		paste("Univariate: ", responses[i], sep=""),
		"ALR Predictors portion of Full Model (Covariates + ALR):",
		"\n",
		capture.output(print(univar_summary_wsignf, quote=F))
	));

	plot_text(c(
		paste("Univariate: ", responses[i], sep=""),
		"Full Model (covariates + ALR) Summary:",
		paste("Multiple R-squared: ", rsquared[1], ", Adjusted R-squared: ", rsquared[2], sep=""),
		paste("F-statistic: ", fstat[["value"]], " on ", 
					fstat[["numdf"]], " and ", 
					fstat[["dendf"]], " DF, p-value: ", pval, sep=""),
		"",
		"",
		"Reduced Model (covariates only) Summary:",
		paste("Multiple R-squared: ", reduced_rsquared[1], ", Adjusted R-squared: ", reduced_rsquared[2], sep=""),
		"",
		"Difference (contribution of ALR):",
		paste("Multiple R-squared: ", rsqrd_diff, 
			", Adjusted R-squared: ", adj_rsqrd_diff, sep=""),
		"",
		"(Positive Adjusted R-squared differences suggests that including ALR predictors improved the model)"
	));

	plot_fit(lmfit, lm_summaries, i);

	summary_res_coeff[,i]=univar_summary[,"Estimate"];
	summary_res_pval[,i]=univar_summary[,"Pr(>|t|)"];
	summary_res_rsqrd[i,]=c(rsquared[2], reduced_rsquared[2], adj_rsqrd_diff);
}

#summary_res_coeff=round(summary_res_coeff,2);
summary_res_pval_rnd=round(summary_res_pval,2);

cat("Coefficients Matrix:\n");
print(summary_res_coeff);
#par(oma=c(10,14,5,1));
if(length(covariate_coefficients)>0){
	paint_matrix(summary_res_coeff[covariate_coefficients,,drop=F], title="Covariate Coefficients",value.cex=2, deci_pts=2);
}

# Variations of ALR Predictor Coefficients
paint_matrix(summary_res_coeff[shrd_alr_names,,drop=F], title="ALR Predictors Coefficients (By Decreasing Abundance)", 
	value.cex=1, deci_pts=2, plot_row_dendr=F, plot_col_dendr=F);
paint_matrix(summary_res_coeff[shrd_alr_names,,drop=F], title="ALR Predictors Coefficients (ALR Clusters)", 
	value.cex=1, deci_pts=2, plot_row_dendr=T, plot_col_dendr=F);
paint_matrix(summary_res_coeff[shrd_alr_names,,drop=F], title="ALR Predictors Coefficients (Response Clusters)", 
	value.cex=1, deci_pts=2, plot_row_dendr=F, plot_col_dendr=T);
paint_matrix(summary_res_coeff[shrd_alr_names,,drop=F], title="ALR Predictors Coefficients (ALR and Response Clusters)", 
	value.cex=1, deci_pts=2, plot_row_dendr=T, plot_col_dendr=T);

cat("\nP-values:\n");
#print(summary_res_pval);
#par(oma=c(10,14,5,1));

if(length(covariate_coefficients)>0){
	paint_matrix(summary_res_pval_rnd[covariate_coefficients,,drop=F], title="Covariate P-values", 
		plot_min=0, plot_max=1, high_is_hot=F, value.cex=2, deci_pts=2);
}

# Variations of ALR Predictor P-values
paint_matrix(summary_res_pval_rnd[shrd_alr_names,,drop=F], title="ALR Predictors P-values (By Decreasing Abundance)", 
	plot_min=0, plot_max=1, high_is_hot=F, value.cex=1, deci_pts=2);


mask_matrix=function(val_mat, mask_mat, mask_thres, mask_val){
        masked_matrix=val_mat;
        masked_matrix[mask_mat>mask_thres]=mask_val;
        return(masked_matrix);
}


# Mask coefficients at various pvalue
signf_coef=mask_matrix(summary_res_coeff[shrd_alr_names,,drop=F], summary_res_pval_rnd[shrd_alr_names,,drop=F], .10, 0);
paint_matrix(signf_coef, title="Significant ALR Predictors Coefficients (p-value < .10)", 
	value.cex=1, deci_pts=2, plot_row_dendr=F, plot_col_dendr=F, label_zeros=F);

signf_coef=mask_matrix(summary_res_coeff[shrd_alr_names,,drop=F], summary_res_pval_rnd[shrd_alr_names,,drop=F], .05, 0);
paint_matrix(signf_coef, title="Significant ALR Predictors Coefficients (p-value < .05)", 
	value.cex=1, deci_pts=2, plot_row_dendr=F, plot_col_dendr=F, label_zeros=F);

signf_coef=mask_matrix(summary_res_coeff[shrd_alr_names,,drop=F], summary_res_pval_rnd[shrd_alr_names,,drop=F], .01, 0);
paint_matrix(signf_coef, title="Significant ALR Predictors Coefficients (p-value < .01)", 
	value.cex=1, deci_pts=2, plot_row_dendr=F, plot_col_dendr=F, label_zeros=F);


paint_matrix(summary_res_pval_rnd[shrd_alr_names,,drop=F], title="ALR Predictors P-values (ALR Clusters)", 
	plot_min=0, plot_max=1, high_is_hot=F, value.cex=1, deci_pts=2,
	plot_row_dendr=T
);
paint_matrix(summary_res_pval_rnd[shrd_alr_names,,drop=F], title="ALR Predictors P-values (Response Clusters)", 
	plot_min=0, plot_max=1, high_is_hot=F, value.cex=1, deci_pts=2,
	plot_col_dendr=T
);
paint_matrix(summary_res_pval_rnd[shrd_alr_names,,drop=F], title="ALR Predictors P-values (ALR and Response Clusters)", 
	plot_min=0, plot_max=1, high_is_hot=F, value.cex=1, deci_pts=2,
	plot_row_dendr=T, plot_col_dendr=T
);

###############################################################################

paint_matrix(summary_res_rsqrd, title="Univariate Adjusted R-Squared");

# Report Full/Reduced Improvement ANOVA P-values
Signf=sapply(improv_mat[,"Diff ANOVA P-Value"], sig_char);
plot_text(c(
	"ANOVA Comparing Full (Cov + ALR) and Reduced (Cov only) Models for Improvement",
	"",
	capture.output(print(round(improv_mat[,c(1,2)], 3))),
	"",
	capture.output(print(round(improv_mat[,c(3,4)], 3))),
	"",
	capture.output(print(cbind(
		apply(improv_mat[,c(5,6,7), drop=F], 1:2, function(x){sprintf("%8.4f",x)}), 
		Signf), quote=F))
));

# Report MANOVA
manova_pval_mat=matrix(NA, nrow=length(alr_cat_names), ncol=1, dimnames=list(alr_cat_names, "Pr(>F)"));
if(length(manova_res)>0){
	plot_text(c(
		"MANOVA",
		capture.output(print(manova_res))
	));

	manova_pval=manova_res[,"Pr(>F)", drop=F]
	manova_pval_mat[alr_cat_names,]=manova_pval[alr_cat_names,];
	#par(oma=c(10,14,5,1));
	paint_matrix(manova_pval_mat[alr_cat_names,,drop=F], title="ALR Predictors MANOVA", 
		plot_min=0, plot_max=1, high_is_hot=F, value.cex=1, deci_pts=3);
}

###############################################################################
# Write summary to file
fh=file(paste(OutputRoot, ".alr_as_pred.review.tsv", sep=""), "w");
cat(file=fh, c("Name:", OutputRoot, ""), sep="\t");
cat(file=fh, "\n");
cat(file=fh, "\n");
cat(file=fh, c("Predictors:"), sep="\t");
cat(file=fh, "\n");

category_names=rownames(manova_pval_mat);

# MANOVA Covariates pvalues
cat(file=fh, paste("[Covariates]", "p-value:", "signif:", sep="\t"), "\n", sep="");
for(cov in covariates_arr){
	if(length(manova_res>0)){
		pv=signif(manova_res[cov,"Pr(>F)"], 5);
		cat(file=fh, c(cov, pv, sig_char(pv)), sep="\t");
	}else{
		cat(file=fh, c(cov, NA, ""), sep="\t");
	}
	cat(file=fh, "\n");
}

# MANOVA ALR pvalues
cat(file=fh, c("[ALR Categories]", "p-value:", "signif:"), sep="\t");
cat(file=fh, "\n");
for(cat_name in alr_cat_names){
	if(length(manova_res>0)){
		pv=signif(manova_res[cat_name,"Pr(>F)"], 5);
		cat(file=fh, c(cat_name, pv, sig_char(pv)), sep="\t");
	}else{
		cat(file=fh, c(cat_name, NA, ""), sep="\t");
	}
	cat(file=fh, "\n");
}
cat(file=fh, "\n");

# Response pvalues, R^2,  
num_resp=nrow(summary_res_rsqrd);
resp_name=rownames(summary_res_rsqrd);
cat(file=fh, c("Univariate Adj R-Sqrd:"),"\n");
cat(file=fh, paste("Resp Name:", "Full Model:", "ALR Contrib:", sep="\t"), "\n", sep="");
for(i in 1:num_resp){
	cat(file=fh, paste(
		resp_name[i], 
		summary_res_rsqrd[i, "Full Model"],
		summary_res_rsqrd[i, "Difference"],
	sep="\t"), "\n", sep="");
}
close(fh);

###############################################################################
# Write ALR Predictor Coefficients to file

fh=file(paste(OutputRoot, ".alr_as_pred.alr.covr.coefficients.tsv", sep=""), "w");
num_out_col=ncol(summary_res_coeff);
response_names=colnames(summary_res_coeff);

cat(file=fh, "Grouping:\t", paste(rep(OutputRoot, num_out_col), collapse="\t"), sep="");
cat(file=fh, "\n");
cat(file=fh, "Responses:\t", paste(response_names, collapse="\t"), sep="");
cat(file=fh, "\n\n");

for(var in covariate_coefficients){
	cat(file=fh, var, paste(
		round(summary_res_coeff[var, response_names], 4), collapse="\t"), sep="\t");
	cat(file=fh, "\n");
}

cat(file=fh, "\n");

for(var in shrd_alr_names){
	cat(file=fh, var, paste(
		round(summary_res_coeff[var, response_names], 4), collapse="\t"), sep="\t");
	cat(file=fh, "\n");
}

close(fh);


###############################################################################
# Write ALR p-values to file
# Format: Factors as columns, taxa (predictor) as rows

exp_tab=summary_res_pval[shrd_alr_names,,drop=F];
write.table(exp_tab,  file=paste(OutputRoot, ".alr_as_pred.pvals.tsv", sep=""), sep="\t", quote=F, col.names=NA, row.names=T);

exp_tab=summary_res_coeff[shrd_alr_names,,drop=F];
write.table(exp_tab,  file=paste(OutputRoot, ".alr_as_pred.coefs.tsv", sep=""), sep="\t", quote=F, col.names=NA, row.names=T);


##############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
