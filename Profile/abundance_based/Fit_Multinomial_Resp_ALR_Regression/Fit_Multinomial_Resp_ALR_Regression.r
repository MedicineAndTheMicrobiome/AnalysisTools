#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
library(car);

source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');

options(useFancyQuotes=F);
options(width=200);


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
	"	-y <Categorical (multinomial) response Y's to select from factor file (column name)>\n",
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
	TimeColumn=character();
}

if(length(opt$subject_column)){
	SubjectColumn=opt$subject_column;
}else{
	SubjectColumn=character();
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

	if(ncol(inmat)!=1){
		cat("Reference levels requires 2 columns, Variable Name and Reference Level.\n");
		quit(status=-1);
	}

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
		relevel_target=relevel_names[i];

		if(length(intersect(relevel_target, factor_names))){
			target_level=ref_lev_mat[i, 1];
			tmp=factors[,relevel_target];
			if(length(intersect(target_level, tmp))){
				tmp=relevel(tmp, target_level);
    				factors[,relevel_target]=tmp;
			}else{
				cat("WARNING: Target level '", target_level,
					"' not found in '", relevel_target, "'!!!\n", sep="");
			}
		}else{
			cat("WARNING: Relevel Target Not Found: '", relevel_target, "'!!!\n", sep="");
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

title_page=function(title, subtitle="", notes=""){

	par(mfrow=c(1,1));
	plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n",  xaxt="n", yaxt="n",
		xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
		);

	text(.5, .8, title, adj=c(.5, -1), cex=4, font=2); 
	text(.5, .5, subtitle, adj=c(.5, 0), cex=2, font=1); 
	text(.5, .4, notes, adj=c(.5, 1), cex=1, font=3); 
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
                plot_min=min(c(mat,Inf), na.rm=T);
        }
        if(is.na(plot_max)){
                plot_max=max(c(mat,-Inf), na.rm=T);
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
	#if(any(is.na(mat))){
	#	plot_col_dendr=F;
	#	plot_row_dendr=F;
	#}

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
			barplot(prop.table(table(vals)), main=colname[i], col="white", las=2);
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

add_signf_char_to_matrix=function(pval_mat, coef_mat){
	numrows=nrow(pval_mat);
	numcols=ncol(pval_mat);

	out_list=list();

	if(numrows==0){
		out_list[["pvalue"]]=pval_mat;
		out_list[["coeff"]]=coef_mat;
		return(out_list);
	}else{

		pval_outmat=matrix(sprintf("%6.3f",pval_mat), ncol=numcols);
		colnames(pval_outmat)=colnames(pval_mat);
		rownames(pval_outmat)=rownames(pval_mat);
		pval_outmat=cbind(rbind(pval_outmat, "", ""),"", "");

		coef_outmat=matrix(sprintf("%6.3f",coef_mat), ncol=numcols);
		colnames(coef_outmat)=colnames(coef_mat);
		rownames(coef_outmat)=rownames(coef_mat);
		coef_outmat=cbind(rbind(coef_outmat, "", ""),"", "");

		sigfun=function(x){
			mx=min(c(x, Inf), na.rm=T);
			return(sig_char(mx));
		}

		sigrow=apply(pval_mat, 1, sigfun);
		sigcol=apply(pval_mat, 2, sigfun);

		pval_outmat[1:numrows, numcols+2]=sigrow;
		pval_outmat[numrows+2, 1:numcols]=sigcol;

		coef_outmat[1:numrows, numcols+2]=sigrow;
		coef_outmat[numrows+2, 1:numcols]=sigcol;

		out_list[["pvalue"]]=pval_outmat;
		out_list[["coeff"]]=coef_outmat;

		return(out_list);
	}

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
	zero_refline=F,
	one_refline=F,
	label=F
	){

	num_times=ncol(stat_mat);
	num_rows=nrow(stat_mat);

	grp_names=rownames(stat_mat);
	time_names=colnames(stat_mat)
	time_values=as.numeric(time_names);

	cat("Times:\n");
	print(time_values);

	if(nrow(stat_mat)==0){
		stat_range=c(Inf, -Inf);
	}else{
		stat_range=range(stat_mat, na.rm=T);
	}
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

	if(zero_refline){
		plot_ymax=max(plot_ymax, 0);
		plot_ymin=min(plot_ymin, 0);
	}
	if(one_refline){
		plot_ymax=max(plot_ymax, 1);
		plot_ymin=min(plot_ymin, 1);
	}


	if(is.null(grp_colors)){
		grp_colors=rainbow(num_rows, start=0, end=5/6, alpha=0.33);
		names(grp_colors)=grp_names;
		print(grp_colors);
	}

	if(nlog10_reflines){
		par(mar=c(4, 4, 5, 4));
	}else{
		par(mar=c(4, 4, 5, 1));
	}

	plot(0, type="n", xlab="Time", ylab="", 
		xlim=c(plot_xmin, plot_xmax), ylim=c(plot_ymin, plot_ymax));

	title(main=title, cex.main=3, font.main=2);
	title(main=subtitle, line=.8, cex.main=.85);

	if(nlog10_reflines){
		abline(h=ref_lines, lty=c(3,2,5), col="grey");
		axis(side=4, at=ref_lines, labels=c("0.10", "0.05", "0.01"), las=2);
	}

	if(zero_refline){
		abline(h=0, col="grey");
	}

	if(one_refline){
		abline(h=1, col="grey");
	}


	if(num_rows>0){
		for(i in 1:num_rows){
			points(time_values, stat_mat[i,], type="b", col=grp_colors[i], lwd=4);
		}
	}else{
		text((plot_xmin+plot_xmax)/2, (plot_ymin+plot_ymax)/2, "[This plot intentionally left blank.]", cex=2);
	}

	if(label==T){
		stat_mag_mat=abs(stat_mat);
		max_time=apply(stat_mag_mat, 1, function(x){ min(which(max(abs(x), na.rm=T)==x))});

		if(num_rows>0){
			for(i in 1:num_rows){
				varname=grp_names[i];
				cat(time_values[max_time[i]], " / ", 
					stat_mat[i, max_time[i]], "/", varname, "\n");


				xpos=time_values[max_time[i]];
				ypos=stat_mat[i, max_time[i]];

				relpos=c(xpos/plot_xmax, ypos/plot_ymax);

				text(xpos, ypos, paste(varname, " [", signif(ypos, 2), "]", sep=""), adj=relpos); 
			}
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

plot_predictions_over_time=function(pred_time_mat, subj_resp_map, resp_name, model_name, resp_colors, grp_tick=F){

	cat("Plotting Predictions Over Time:\n");

	print(pred_time_mat);
	print(subj_resp_map);
	print(resp_name);
	print(model_name);
	print(resp_colors);

	num_time_pts=ncol(pred_time_mat);
	num_samples=nrow(pred_time_mat);

	time_pts=colnames(pred_time_mat);
	sample_ids=rownames(pred_time_mat);

	uniq_resp_names=unique(subj_resp_map);

	time_vals=as.numeric(time_pts);
	time_min=min(time_vals);
	time_max=max(time_vals);
	time_diff=(time_max-time_min);
	xpad=time_diff*.05;

	if(num_time_pts==1){
		min_time_sep=1;
	}else{
		min_time_sep=min(diff(time_vals));
	}

	layout_mat=matrix(c(1,1,2,2,3), ncol=1);
	layout(layout_mat);
	par(mar=c(4,4,4,1));
	
	jitter=rnorm(num_samples, 0, min_time_sep/6); 

	#############################################################################
	plot(0,0, type="n", 
		xlim=c(time_min-xpad, time_max+xpad),
		ylim=c(0,1),
		xlab="Time",
		ylab=paste("Probability of Being Classified as ", resp_name, sep="")
	);

	abline(h=c(0,.5,1), col="grey", lty=2);

	for(t in time_pts){
		tval=as.numeric(t);

		abline(v=tval, col="grey");

		for(s in 1:num_samples){	
			color=resp_colors[subj_resp_map[sample_ids[s]]];
			points(tval+jitter[s], pred_time_mat[s,t], col=color);
		}
	}
	title(main=paste("Response: ", resp_name, " / Model: ", model_name, " Prediction in Probabilities", sep=""));

	#############################################################################

	logit_matrix=log(pred_time_mat/(1-pred_time_mat));
	numeric_ix=is.finite(logit_matrix);	
	numeric_val=logit_matrix[numeric_ix];

	logit_mag=max(abs(numeric_val));
	if(logit_mag==Inf){
		logit_mag=1e300;
	}else if(logit_mag==-Inf){
		logit_mag=1e300;
	}
	logit_matrix[logit_matrix==Inf]=logit_mag+1;
	logit_matrix[logit_matrix==-Inf]=-logit_mag-1;
	logit_max=logit_mag+1;


	plot(0,0, type="n", 
		xlim=c(time_min-xpad, time_max+xpad),
		ylim=c(-logit_max,logit_max),
		xlab="Time",
		ylab=paste("Logit(Probability of ", resp_name, " Classification)", sep="")
	);

	abline(h=0, col="grey", lty=2);

	for(t in time_pts){

		tval=as.numeric(t);

		if(grp_tick){
			cur_time_val=logit_matrix[,t,drop=F];
			
			for(rpix in uniq_resp_names){
				rp_sbj=names(subj_resp_map[subj_resp_map==rpix]);
				mn=mean(cur_time_val[rp_sbj,], na.rm=T);

				points(c(tval-min_time_sep/8, tval+min_time_sep/8), 
					c(mn, mn), col=resp_colors[rpix], lwd=2, type="l");
			}

		}


		abline(v=tval, col="grey");

		for(s in 1:num_samples){	
			color=resp_colors[subj_resp_map[sample_ids[s]]];
			lo=logit_matrix[s,t];
			points(tval+jitter[s], lo, col=color);
		}
	}
	title(main=paste("Response: ", resp_name, " / Model: ", model_name, 
		" Predictions in logit(Probability)s", sep=""));

	#############################################################################

	plot_group_legend(resp_colors)

}

plot_signif_variables_over_time=function(factors, resp_cnm, time_cnm, sbj_cnm, 
	resp, targ_var, times, sbj_to_resp_map, resp_to_color_map, model_ix){

	if(length(targ_var)==0){
		plot(0,0, type="n",
			xlim=c(-1,1), ylim=c(-1,1),
			xlab="", ylab="", xaxt="n", yaxt="n",
			main="", bty="n");
		text(0,0, 
			paste(resp, " / ", model_ix, ":\n\nNo significant Predictors\nto Plot.", sep=""), 
			cex=3);
		return();
	}

	targ_var=intersect(targ_var, colnames(factors));
	signf_factors=factors[,targ_var,drop=F];

	responses_val=as.character(factors[,resp_cnm]);

	time_val=as.numeric(as.character(factors[,time_cnm]));
	uniq_times=sort(unique(time_val));
	num_uniq_times=length(uniq_times);

	if(num_uniq_times>1){
		min_time_span=min(diff(uniq_times));
	}else{
		min_time_span=1;
	}

	sbj_val=as.character(factors[,sbj_cnm]);
	uniq_sbj=unique(sbj_val);
	num_sbj=length(unique(sbj_val));

	resp_group_val=sbj_to_resp_map[sbj_val];
	uniq_resp_grps=sort(unique(resp_group_val));
	num_uniq_resp_grps=length(uniq_resp_grps);

	#print(signf_factors);
	#print(responses_val);
	#print(time_val);
	#print(sbj_val);

	time_range=range(times);
	xpad=(time_range[2]-time_range[1])*.05;

	num_plots=length(targ_var);

	num_rows_per_page=5;

	par(mfrow=c(num_rows_per_page,1));

	jitter=rnorm(num_sbj, 0, min_time_span/8);
	names(jitter)=uniq_sbj;
	cat("Jitter:\n");
	print(jitter);

	plot_ix=0;

	if(max(nchar(uniq_sbj))<=3){
		label_pts=T;
	}else{
		label_pts=F;
	}

	for(var_ix in targ_var){

		vals=signf_factors[,var_ix];

		grp_by_time_means=matrix(NA, nrow=num_uniq_resp_grps, ncol=num_uniq_times);
		colnames(grp_by_time_means)=as.character(uniq_times);
		rownames(grp_by_time_means)=uniq_resp_grps;
		
		# Calculate means for groups over time
		for(rgrp in uniq_resp_grps){
			for(t in uniq_times){
				cur_time_ix=(time_val==t);
				cur_rgrp_ix=(resp_group_val==rgrp);
				grp_by_time_means[rgrp, as.character(t)]=
					mean(vals[cur_time_ix & cur_rgrp_ix], na.rm=T);
			}
		}
		
		val_rng=range(vals);
		rng_spn=val_rng[2]-val_rng[1];
		ypad=rng_spn*.1;

		par(mar=c(2,2,4,1));
		plot(0,0, type="n", 
			xlim=c(time_range[1]-xpad, time_range[2]+xpad),
			ylim=c(val_rng[1]-ypad, val_rng[2]+ypad),
			main=paste(resp, "/", model_ix, ":  ", var_ix, " Values", sep=""), 
			xlab="time", ylab="");

		abline(v=uniq_times, col="grey");

		# Plot mean lines
		for(rgrp in uniq_resp_grps){
			if(num_uniq_times>1){
				points(uniq_times, grp_by_time_means[rgrp,], type="l", col=resp_to_color_map[rgrp]);
			}else{
				points(
					c(uniq_times-1/8, uniq_times+1/8),
					rep(grp_by_time_means[rgrp,],2),
					 type="l", col=resp_to_color_map[rgrp], lwd=2);
			}
		}

		# Plot scattered points
		for(sbj_ix in uniq_sbj){
			cur_sbj=(sbj_val==sbj_ix);			

			cur_vals=vals[cur_sbj];
			cur_time=time_val[cur_sbj];

			time_order=order(cur_time);
			cur_vals=cur_vals[time_order];
			cur_time=cur_time[time_order];

			resp_grp=sbj_to_resp_map[sbj_ix];
			sbj_col=resp_to_color_map[resp_grp];

			points(cur_time+jitter[sbj_ix], cur_vals, type="p", col=sbj_col);

			if(label_pts){
				text(cur_time+jitter[sbj_ix], cur_vals, sbj_ix, cex=.3);
			}
		}

		plot_ix=plot_ix+1;
		if(plot_ix == (num_rows_per_page-1)){
			plot_group_legend(resp_to_color_map);
			plot_ix=0;	
		}
	}

	if(plot_ix<num_rows_per_page && plot_ix!=0){
		plot_group_legend(resp_to_color_map);
	}

}

output_predictions_over_time=function(
	output_fn_root,
	time_mat, 
	subj_id_to_resp_map, resp, model){

	cat("\n\n");
	cat("Inside: output_predictions_over_time\n");

	cat("Output File: ", output_fn_root, "\n");
	cat("Response: ", resp, "\n");
	cat("Model: ", model, "\n");

	#cat("Time Mat:\n");
	#print(time_mat);

	#cat("Sbj->Resp Map:\n");
	#print(subj_id_to_resp_map);

	# Make directories for all the different combinations of output
	predictions_dir=paste(output_fn_root, ".multn_predictions", sep="");
	if(!dir.exists(predictions_dir)){
		dir.create(predictions_dir);
	}

	model_dir=paste(predictions_dir, "/model_[", model, "]", sep="");
	if(!dir.exists(model_dir)){
		dir.create(model_dir);
	}

	responses_dir=paste(model_dir, "/pred_as_[", resp, "]", sep="");
	if(!dir.exists(responses_dir)){
		dir.create(responses_dir);
	}

	# Output separate lists for each sample/resp category type	
	resp_types=unique(subj_id_to_resp_map);
	for(rtype in resp_types){
		cat("Obs Resp Types: ", rtype, "\n");
		samp_ids=(rtype==subj_id_to_resp_map);
		probs_mat=time_mat[samp_ids,,drop=F];
		print(probs_mat);

		full_path=paste(responses_dir, "/obsr_as_[", rtype, "].tsv", sep="");
		write.table(probs_mat, full_path, sep="\t", col.names=NA, row.names=T, quote=F);
	}
	
	cat("Leaving: output_predictions_over_time\n");

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

required_arr=ResponseColname;

if(""!=RequiredFile){
	list_required=load_list(RequiredFile);
	cat("List Required Variables:\n");
	print(list_required);
	cat("\n");
	required_arr=c(required_arr, list_required);
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

if(length(TimeColumn)){
	resp_cov_factors=cbind(subject_ids, kept_factors[,TimeColumn, drop=F], response_factors, covariate_factors);
}else{
	resp_cov_factors=cbind(subject_ids, response_factors, covariate_factors);
}

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

response_reference_level=as.character(levels(response_factors[,1])[1]);
cat("Response Summary:\n");
s=summary(response_factors);
plot_text(c(
	"Response Summary:",
	"\n",
	capture.output(print(s)),
	"",
	paste("Reference Level: ", response_reference_level, "\n", sep="")
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

cat("Setting up time id variable...\n");
if(length(TimeColumn)==0){
	cat("No Time Column Specified... Assuming Cross-sectional Analysis.\n");
	time_ids=rep(0, nrow(factors_wo_nas));
	factors_wo_nas=cbind(factors_wo_nas, time_ids);
	TimeColumn="time_ids";
}else{
	time_ids=factors_wo_nas[,TimeColumn];
}

subject_ids=character();
cat("Setting up subject id variable...\n");
if(length(SubjectColumn)){
	subject_ids=factors_wo_nas[, SubjectColumn];
}else{
	cat("Subject IDs not specified, using sample IDs.\n");
	subject_ids=rownames(factors_wo_nas);
	factors_wo_nas=cbind(factors_wo_nas, subject_ids);
	SubjectColumn="subject_ids";
}	

cat("\n");
cat("All Time IDs:\n");
print(time_ids);
cat("All Subject IDs:\n");
print(subject_ids);

responses=factors_wo_nas[,ResponseColname];

unique_time_ids=sort(unique(time_ids));
unique_subject_ids=sort(unique(subject_ids));
unique_responses=sort(unique(responses));

num_unique_time_ids=length(unique_time_ids);
num_unique_subject_ids=length(unique_subject_ids);
num_unique_responses=length(unique_responses);

cat("\n");
cat("Unique Time Points:\n");
print(unique_time_ids);
cat("Unique Subject IDs:\n");
print(unique_subject_ids);
cat("\n");
cat("Unique Response Groups:\n");
print(unique_responses);


# Setup sample to subject map
sample_id_to_subject_map=factors_wo_nas[, c(SubjectColumn, ResponseColname), drop=F];
cat("Sample to Subject Map:\n");
print(sample_id_to_subject_map);
plot_text(c(
	paste("Sample ID to Subject (", SubjectColumn, ") Map:", sep=""),
	"",
	capture.output(print(sample_id_to_subject_map, quotes=F))
));

# Setup subject response map
subject_id_to_response_map=rep(NA, num_unique_subject_ids);
names(subject_id_to_response_map)=unique_subject_ids;

for(sbj_ix in unique_subject_ids){
	ix=min(which(sbj_ix==factors_wo_nas[, SubjectColumn]));
	subject_id_to_response_map[sbj_ix]=as.character(factors_wo_nas[ix,ResponseColname]);
}

cat("\n");
cat("Subject to Response Map:\n");
print(subject_id_to_response_map);
plot_text(c(
	paste("Subject (", SubjectColumn, ") to Response (", ResponseColname, ") Map:", sep=""),
	"",
	capture.output(print(subject_id_to_response_map, quotes=F))
));

###############################################################################

sample_sizes_resp_time=matrix(NA, nrow=num_unique_responses, ncol=num_unique_time_ids);
rownames(sample_sizes_resp_time)=unique_responses;
colnames(sample_sizes_resp_time)=as.character(unique_time_ids);

for(rix in unique_responses){
	for(tix in as.character(unique_time_ids)){
		cur_resp_ix=(rix==responses);
		cur_time_ix=(as.numeric(tix)==time_ids);
		tm_rsp_ix=(cur_time_ix & cur_resp_ix);

		#cat(rix, " ", tix, " ", tm_rsp_ix, "\n");
		sample_sizes_resp_time[rix,tix]=sum(tm_rsp_ix);
	}
}

###############################################################################

require(nnet);
require(generalhoslem);

process_model=function(fit, null_fit=NULL, num_samples=NULL){
	
	res=list();
	res[["fit"]]=fit;
	res[["summary"]]=summary(fit);

	coeff=res[["summary"]]$coefficients
	stderr=res[["summary"]]$standard.errors;
	z=coeff/stderr;
	res[["pvalues"]]=(1-pnorm(abs(z), 0, 1))*2;

	res[["num_coeff"]]=nrow(coeff);
	res[["num_respGrps"]]=ncol(coeff);
	res[["num_predictors_exclIntercept"]]=(res[["num_coeff"]]-1)*res[["num_respGrps"]];
	res[["num_predictors"]]=(res[["num_coeff"]])*res[["num_respGrps"]];
	res[["num_samples"]]=num_samples;

	pred_names=colnames(res[["pvalues"]]);
	pred_names=setdiff(pred_names, "(Intercept)");

	signf_counts=res[["pvalues"]][,pred_names, drop=F]<0.1;
	res[["num_signif"]]=apply(signf_counts, 1, function(x){sum(x, na.rm=T);});

	# See M. W. Fagerland and D. W. Hosmer,
	# "A generalized Hosmer–Lemeshowgoodness-of-fit test for multinomial logisticregression models"
	# Null distribution is if there was a good fit, so low p-values are bad fits.
	#res[["gof"]]=logitgof(resp, fitted(fit));
	res[["NegLogLikelihood"]]=res[["summary"]]$value;
	res[["PackageAIC"]]=fit["AIC"];

	# McFadden's Adj pseudo-R^2 is 1-((LLfull-K)/LLintercep), but out of the computational
	# fit we get -LL, so the formula below is modified.

	if(!is.null(null_fit)){
	
		full_negLL=res[["NegLogLikelihood"]];
		rcdc_negLL=null_fit[["NegLogLikelihood"]];

		res[["R2_McFadden"]]=1-(full_negLL/rcdc_negLL);

		K=res[["num_predictors_exclIntercept"]];

		res[["R2_McFadden_Adj"]]=1-(full_negLL+K)/rcdc_negLL;
	
		res[["R2_CoxSnell"]]=1-(exp(full_negLL-rcdc_negLL))^(2/num_samples);

	}else{
		cat("Pseudo-R^2's not calculated for Intercept-Only (NULL) Model.\n");
	}

	return(res);

}

fit_info=list();

standardize_columns=function(datamat, target_var){
	out_mat=datamat;

	for(v in target_var){
		val=datamat[,v];
		if(is.numeric(val)){
			mean=mean(val);
			sd=sd(val);
			if(sd==0){
				norm=val-mean;
			}else{
				norm=(val-mean)/sd;
			}
			out_mat[,v]=norm;
		}else{
			out_mat[,v]=val;
		}
	}
	
	return(out_mat);
}

standardize=T;

for(cur_time_id in unique_time_ids){

	cat("\n*******************************************************************\n");
	cat("Working on time: ", cur_time_id, "\n");

	cur_time_ix=(cur_time_id==time_ids);
	cur_time_str=sprintf("%02g", cur_time_id);

	cur_factors=factors_wo_nas[cur_time_ix,,drop=F];
	cur_alr=as.data.frame(alr_categories_val[cur_time_ix,,drop=F]);

	cur_responses=cur_factors[,ResponseColname];
	cur_predictors=cur_factors[,covariates_arr];

	if(standardize){
		cur_factors=standardize_columns(cur_factors, target=covariates_arr);
		cur_alr=standardize_columns(cur_alr, target=alr_cat_names);
	}

	num_samples_at_curtime=nrow(cur_factors);
	cat("Num Samples: ", num_samples_at_curtime, "\n", sep="");
	cat("\n");

	fit_info[[cur_time_str]]=list();

	# Fit intercept only
	cat("Intercept Only:\n");
	null_model_str=paste("cur_responses ~ 1");
	null_mlr_fit=multinom(as.formula(null_model_str), data=cur_factors);
	fit_info[[cur_time_str]][["intercept_only"]]=process_model(null_mlr_fit, NULL, num_samples_at_curtime);
	cat("\n");

	# Fit covariates
	cat("Covariates Only:\n");
	cov_model_str=paste("cur_responses ~ ", paste(covariates_arr, collapse=" + ", sep=""));
	cov_mlr_fit=multinom(as.formula(cov_model_str), data=cur_factors);
	fit_info[[cur_time_str]][["covariates_only"]]=
		process_model(cov_mlr_fit, fit_info[[cur_time_str]][["intercept_only"]], num_samples_at_curtime);
	cat("\n");

	# Fit alr categories
	cat("ALR Categories Only:\n");
	alr_model_str=paste("cur_responses ~ ", paste(alr_cat_names, collapse=" + ", sep=""));
	alr_mlr_fit=multinom(as.formula(alr_model_str), data=cur_alr);
	fit_info[[cur_time_str]][["alr_only"]]=
		process_model(alr_mlr_fit, fit_info[[cur_time_str]][["intercept_only"]], num_samples_at_curtime);
	cat("\n");

	# Fit combined 
	cat("Full Combined:\n");
	comb_model_str=paste("cur_responses ~ ", paste(c(covariates_arr, alr_cat_names), collapse=" + ", sep=""));
	comb_mlr_fit=multinom(as.formula(comb_model_str), data=cbind(cur_factors, cur_alr));
	fit_info[[cur_time_str]][["alr_and_covariates"]]=
		process_model(comb_mlr_fit, fit_info[[cur_time_str]][["intercept_only"]], num_samples_at_curtime);
	cat("\n");

	# Overall statistics
	fit_info[[cur_time_str]][["responses"]]=table(as.character(cur_responses));
	fit_info[[cur_time_str]][["num_responses"]]=length(cur_responses);

	#print(fit_info[[cur_time_str]][["combined"]][["pvalues"]]);

}

print(fit_info);

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
mcfadden_r2_mat=stat_matrix;
mcfadden_r2adj_mat=stat_matrix;
coxsnell_r2_mat=stat_matrix;

sampsize_arr=stat_arr;

resp_grps=matrix(NA, ncol=num_time_pts, nrow=num_unique_responses);
colnames(resp_grps)=time_str_ids;
rownames(resp_grps)=unique_responses;

pvalues=list();
coefficients=list();
predicted=list();


for(cur_time_str  in time_str_ids){

	cat("Extracting Time: ", cur_time_str, "\n");

	sampsize_arr[cur_time_str]=fit_info[[cur_time_str]][["num_responses"]];

	resp_grps[unique_responses, cur_time_str]=fit_info[[cur_time_str]][["responses"]][unique_responses];

	pvalues[[cur_time_str]]=list();
	coefficients[[cur_time_str]]=list();
	predicted[[cur_time_str]]=list();

	for(modix in model_types){
		cat("Extracting: ", modix, "\n");

		aic_matrix[modix, cur_time_str]=fit_info[[cur_time_str]][[modix]][["fit"]][["AIC"]];

		pvalues[[cur_time_str]][[modix]]=fit_info[[cur_time_str]][[modix]][["pvalues"]];

		coefficients[[cur_time_str]][[modix]]=fit_info[[cur_time_str]][[modix]][["Coefficients"]]=
			summary(fit_info[[cur_time_str]][[modix]][["fit"]])$coefficients;

		predicted[[cur_time_str]][[modix]]=
			fit_info[[cur_time_str]][[modix]][["fit"]][["fitted.values"]];

		mcfadden_r2_mat[modix, cur_time_str]=
			fit_info[[cur_time_str]][[modix]][["R2_McFadden"]];
		mcfadden_r2adj_mat[modix, cur_time_str]=
			fit_info[[cur_time_str]][[modix]][["R2_McFadden_Adj"]];
		coxsnell_r2_mat[modix, cur_time_str]=
			fit_info[[cur_time_str]][[modix]][["R2_CoxSnell"]];

	}

}

cat("\n");
cat("Sample Sizes:\n");
print(sampsize_arr);

cat("\n");
cat("Response Group Sizes:\n");
print(sample_sizes_resp_time);

layout_m=matrix(c(1,1,1,1,1,2), nrow=6, ncol=1);
layout(layout_m);

plot_ts_stat_table(sample_sizes_resp_time, title="Response Group Sizes", 
	subtitle="Number of Samples per Group Over Time", grp_colors=grp_colors, plot_ymin=0, label=T);
plot_group_legend(grp_colors);

# Plot model fits

cat("AIC:\n");
print(aic_matrix);
layout_m=matrix(c(1,1,1,1,1,2), nrow=6, ncol=1);
layout(layout_m);

plot_ts_stat_table(aic_matrix, 
	title="AIC", subtitle="Lower Values, Better Fit", 
	grp_colors=model_colors, label=T);
plot_group_legend(model_colors);

plot_ts_stat_table(mcfadden_r2_mat, 
	zero_refline=T,one_refline=T,
	title="McFadden's Pseudo-R^2", subtitle="Closer to 1 the better", 
	grp_colors=model_colors, label=T);
plot_group_legend(model_colors);

plot_ts_stat_table(mcfadden_r2adj_mat, 
	zero_refline=T,one_refline=T,
	title="McFadden's (Adjusted) Pseudo-R^2", subtitle="(Adjusted for Num Predictors)", 
	grp_colors=model_colors, label=T);
plot_group_legend(model_colors);

plot_ts_stat_table(coxsnell_r2_mat, 
	zero_refline=T,one_refline=T,
	title="Cox & Snell's Pseudo-R^2", subtitle="(Adjusted for Sample Size)", 
	grp_colors=model_colors, label=T);
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
pred_bytime_list=list();

# for each response category
for(resp_ix in response_no_reference){
	cat("Extracting Response: ", resp_ix, "\n");
	# for each 4 model types

	coef_bytime_list[[resp_ix]]=list();
	pval_bytime_list[[resp_ix]]=list();
	pred_bytime_list[[resp_ix]]=list();

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

		pred_by_time_matrix=matrix(NA, nrow=num_unique_subject_ids, ncol=num_time_pts);
		rownames(pred_by_time_matrix)=unique_subject_ids;
		colnames(pred_by_time_matrix)=time_str_ids;

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

			cur_pred_mat=predicted[[time_ix]][[model_ix]];
			cur_samp_ids=rownames(cur_pred_mat);

			avail_responses=colnames(cur_pred_mat);
			if(any(resp_ix==avail_responses)){
				for(smp_ix in cur_samp_ids){
					sbj_id=sample_id_to_subject_map[smp_ix, SubjectColumn];
					pred_by_time_matrix[sbj_id, time_ix]=cur_pred_mat[smp_ix, resp_ix];
				}
			}

		}

		# remove intercept
		itcp_ix=which("(Intercept)"==rownames(coef_by_time_matrix));
		coef_by_time_matrix=coef_by_time_matrix[-itcp_ix,,drop=F];
		pval_by_time_matrix=pval_by_time_matrix[-itcp_ix,,drop=F];

		coef_bytime_list[[resp_ix]][[model_ix]]=coef_by_time_matrix;
		pval_bytime_list[[resp_ix]][[model_ix]]=pval_by_time_matrix;
		pred_bytime_list[[resp_ix]][[model_ix]]=pred_by_time_matrix;
	}
}

neglog10_trans_mat=function(mat){
	minnonzero=min(c(mat[mat!=0], Inf), na.rm=T);		
	cat("Min Nonzero Pvalue: ", minnonzero, "\n");
	mat[mat==0]=minnonzero/10;
	nlog10pval=-log10(mat);
	return(nlog10pval);
}

outwidth=2000;


#factors_wo_nas,alr_categories_val
if(any(rownames(factors_wo_nas)!=rownames(alr_categories_val))){
	cat("Error: Factor sample IDs don't match ALR sample IDs\n");
	quit(-1)
}else{
	alr_and_covariates_matrix=cbind(factors_wo_nas,alr_categories_val);	
}

list_to_text=function(lis){
	groups=names(lis);
	
	if(length(groups)==0){
		return("No Values.");
	}

	out_line=c();
	for(g in groups){
		out_line=c(out_line, paste("[", g, "]"));
		arr=lis[[g]];
		for(a in arr){
			out_line=c(out_line, paste("   ", a));
		}	
		out_line=c(out_line, "");
	}
	return(out_line);
}

for(resp_ix in response_no_reference){
		
	title_page(
		resp_ix, 
		subtitle="as a response", 
		notes=paste("Model Diagnostics:", "", paste(model_types, collapse="\n"), "", 
			paste("Comparison of Models for ", resp_ix, sep=""), sep="\n")
	);

	signf_coef_by_model=list();
	signf_pval_by_model=list();

	for(model_ix in model_types){

		title_page(resp_ix, 
			subtitle=paste("Model:", model_ix),
			note=paste(c(
				"All Coefficients Table", 
				"All Coefficients Dendrogram/Heatmap",
				"All Coefficients Plot",
				"All P-Values Table", 
				"All P-Values Dendrogram/Heatmap",
				"All -log10(P-Values) Plot", 
				"",
				"Only Significant Coefficients Table", 
				"Only Significant Coefficients Plot", 
				"Only Significant P-Values Table", 
				"Only Significant P-Values Plot", 
				"",
				"Predictions (Probability & Logit)", 
				"Significant Predictors Plots"
			), collapse="\n")
		);

		cur_coef_by_time_matrix=coef_bytime_list[[resp_ix]][[model_ix]];
		cur_pval_by_time_matrix=pval_bytime_list[[resp_ix]][[model_ix]];
		cur_pred_by_time_matrix=pred_bytime_list[[resp_ix]][[model_ix]];

		annotated_signf_mat_list=
			add_signf_char_to_matrix(cur_pval_by_time_matrix, cur_coef_by_time_matrix);
	
		num_model_coeff=nrow(cur_coef_by_time_matrix);

		var_colors=rainbow(num_model_coeff, start=0+1/12, end=5/6+1/12, alpha=3/8);
		names(var_colors)=rownames(cur_coef_by_time_matrix);

		#----------------------------------------------------------------------------
		# Plot coefficients
		par(mfrow=c(1,1));
		coef_tab=capture.output(print(annotated_signf_mat_list[["coeff"]], quote=F));
		plot_text(c(
			paste("Response: ", resp_ix, "  Model: ", model_ix),
			"",
			paste("Reference Level: ", response_reference_level, sep=""),
			"Coefficients:",
			"",
			coef_tab
		));

		paint_matrix(
			mat=cur_coef_by_time_matrix, 
			title=paste("Response: ", resp_ix, " / Model:  ", model_ix, " :  Coefficients", sep=""), 
			value.cex=1.5, plot_row_dendr=T);

		layout_m=matrix(c(1,1,1,2), nrow=4, ncol=1);
		layout(layout_m);
		plot_ts_stat_table(cur_coef_by_time_matrix, 
			title=paste("All Coefficients"),
			subtitle=paste("Response: ", resp_ix, " / Model:  ", model_ix, sep=""),
			grp_colors=var_colors, zero_refline=T);
		plot_group_legend(var_colors);

		#----------------------------------------------------------------------------
		# Plot p-values
		par(mfrow=c(1,1));
		pval_tab=capture.output(print(annotated_signf_mat_list[["pvalue"]], quote=F));

		plot_text(c(
			paste("Response: ", resp_ix, " / Model:  ", model_ix, sep=""),
			"",
			"P-Values:",
			"",
			pval_tab
		));

		paint_matrix(
			mat=cur_pval_by_time_matrix, 
			title=paste("Response: ", resp_ix, " / Model:  ", model_ix, "  P-Values", sep=""), 
			value.cex=1.5, plot_row_dendr=T, high_is_hot=F);

		nlog10pval=neglog10_trans_mat(cur_pval_by_time_matrix);

		layout_m=matrix(c(1,1,1,2), nrow=4, ncol=1);
		layout(layout_m);
		plot_ts_stat_table(nlog10pval, 
			title=paste("All -log10(P-Values)"),
			subtitle=paste("Response: ", resp_ix, " / Model:  ", model_ix, sep=""),
			grp_colors=var_colors, plot_ymin=0, nlog10_reflines=T);
		plot_group_legend(var_colors);

		#----------------------------------------------------------------------------
		# Pull out predictors that were significant p<.1

		signf_coef_ix=apply(cur_pval_by_time_matrix, 1, function(x){ min(c(x, Inf), na.rm=T)<0.1 });
		
		signf_coef=cur_coef_by_time_matrix[signf_coef_ix,,drop=F];
		signf_pval=cur_pval_by_time_matrix[signf_coef_ix,,drop=F];
		signf_colr=var_colors[signf_coef_ix];
		num_signf=sum(signf_coef_ix);

		annotated_signf_mat_list=
			add_signf_char_to_matrix(signf_pval, signf_coef);

		par(mfrow=c(1,1));
		coef_tab=capture.output(
			print(annotated_signf_mat_list[["coeff"]], digits=4, quote=F, width=outwidth));

		plot_text(c(
			paste("Response: ", resp_ix, " / Model: ", model_ix),
			"",
			"Significant Coefficients:",
			"",
			coef_tab
		));
		
		layout_m=matrix(c(1,1,1,2), nrow=4, ncol=1);
		layout(layout_m);
		plot_ts_stat_table(signf_coef, 
			title=paste("Only Significant Coefficients:", sep=""),
			subtitle=paste("Response: ", resp_ix, " / Model:  ", model_ix, sep=""),
			grp_colors=signf_colr, zero_refline=T,
			label=T);
		plot_group_legend(signf_colr);


		par(mfrow=c(1,1));
		pval_tab=capture.output(
			print(annotated_signf_mat_list[["pvalue"]], quote=F, width=outwidth, row.names=F));

		plot_text(c(
			paste("Response: ", resp_ix, " / Model: ", model_ix, sep=""),
			"",
			"Significant P-Values:",
			"",
			pval_tab
		));
		
		layout_m=matrix(c(1,1,1,2), nrow=4, ncol=1);
		layout(layout_m);
		nlog10pval=neglog10_trans_mat(signf_pval);
		plot_ts_stat_table(nlog10pval, 
			title=paste("Only Significant -log10(P-Values)", sep=""),
			subtitle=paste("Response: ", resp_ix, " / Model:  ", model_ix, sep=""),
			grp_colors=signf_colr, plot_ymin=0, nlog10_reflines=T,
			label=T);
		plot_group_legend(signf_colr);

		#----------------------------------------------------------------------------

		if(num_signf>0){
			signf_coef_by_model[[model_ix]]=coef_tab;
			signf_pval_by_model[[model_ix]]=pval_tab;
		}else{
			signf_coef_by_model[[model_ix]]="No significant variables";
			signf_pval_by_model[[model_ix]]="No significant variables";
		}

		plot_predictions_over_time(
			cur_pred_by_time_matrix, 
			subject_id_to_response_map, resp_ix, model_ix, grp_colors, grp_tick=T);

		plot_signif_variables_over_time(
			alr_and_covariates_matrix,
			ResponseColname,
			TimeColumn,
			SubjectColumn,
			resp_ix,
			rownames(signf_coef),
			unique_time_ids,
			subject_id_to_response_map,
			grp_colors,
			model_ix
		);

		output_predictions_over_time(
			OutputRoot,
			cur_pred_by_time_matrix, 
			subject_id_to_response_map, resp_ix, model_ix);
	
	}

	# Output 


	par(mfrow=c(1,1));

	signf_coef_as_text=list_to_text(signf_coef_by_model);
	signf_pval_as_text=list_to_text(signf_pval_by_model);
	plot_text(c(
		paste("Response: ", resp_ix, "    Reference: ", response_reference_level),
		"",
		"Comparison of Models:",
		"",
		"-------------------------------------------------------------------------------------",
		"",
		"Significant Coefficients:",
		"",
		signf_coef_as_text,
		"-------------------------------------------------------------------------------------",
		"",
		"Significant P-Values:",
		"",
		signf_pval_as_text
	));
	

}

plot_variables_by_response_venn=function(coef_bytime, pval_bytime){

	responses=names(coef_bytime);
	models=names(coef_bytime[[responses[1]]]);

	print(responses);
	print(models);

	#print(names(coef_bytime));
	#print(coef_bytime[["AlloTx_noTAC"]][["covariates_only");


	# Get the variable list from each model 
	model_variable_list=list();
	for(model_ix in models){
		model_variable_list[[model_ix]]=rownames(coef_bytime[[responses[1]]][[model_ix]]);
	}


	var_by_resp_list=list();
	for(model_ix in models){

		resp_var_mat=
			matrix(NA, nrow=length(model_variable_list[[model_ix]]), ncol=length(responses));
		rownames(resp_var_mat)=model_variable_list[[model_ix]];
		colnames(resp_var_mat)=responses;

		for(resp_ix in responses){

			coef_tab=coef_bytime[[resp_ix]][[model_ix]];
			pval_tab=pval_bytime[[resp_ix]][[model_ix]];
			
			variables=rownames(coef_tab);

			signf=apply(pval_tab, 1,  
				function(x){
					min(c(x,Inf),na.rm=T);
				});

			resp_var_mat[variables, resp_ix]=signf[variables];

		}

		var_by_resp_list[[model_ix]]=resp_var_mat;
	}

	title_page(
		"Heatmap Comparisons", 
		paste("\nAll P-Values\n\nBetween Response Groups:", "", 
		paste(response_no_reference, collapse="\n"), sep="\n"),
		paste("Across Models:", "", paste(models, collapse="\n"), sep="\n")
	);

	for(model_ix in models){
		paint_matrix(
			mat=var_by_resp_list[[model_ix]], 
			title=paste(model_ix, ": p-values", sep=""),
			plot_min=0, plot_max=1,
			high_is_hot=F
		);
	}

	title_page(
		"Exclusivity Heatmaps",
		paste("\nSignificant P-values (p<0.1)\n\nBetween Response Groups:",
			paste(response_no_reference, collapse="\n"), sep="\n"),
		paste("Across Models:", paste(models, collapse="\n"), sep="\n")
	);

	for(model_ix in models){

		signf_mat=var_by_resp_list[[model_ix]] < .1;
		print(signf_mat);

		var_names=rownames(signf_mat);
		excl_to_resp_list=list();
		incl_to_all_resp=character();
	
		# Inclusive to all responses
		incl_all_ix=apply(signf_mat, 1, all);
		incl_to_all_resp=var_names[incl_all_ix];
	
		# Exclusive to one response
		tot_incl_resp=apply(signf_mat, 1, function(x){sum(x)==1});
		for(resp_ix in response_no_reference){
			excl_to_resp_list[[resp_ix]]=var_names[tot_incl_resp & signf_mat[,resp_ix]];
		}

		paint_matrix(
			mat=signf_mat, 
			title=paste(model_ix, ": p-values < 0.1", sep=""),
			counts=T, label_zeros=F,
			plot_min=0, plot_max=1,
			high_is_hot=T
		);

		plot_text(c(
			paste("Model: ", model_ix, sep=""),
			"Exclusivity of Predictors:",
			"",
			list_to_text(excl_to_resp_list),
			"",
			"",
			"Shared Across All Response Groups:",
			"",
			list_to_text(incl_to_all_resp)
		));
	}
	
}

plot_variables_by_response_venn(coef_bytime_list, pval_bytime_list);

###############################################################################

title_page("Association Comparisons", "Between Responses", "Over Time");

mask_matrix=function(val_mat, mask_mat, mask_thres, mask_val){
        masked_matrix=val_mat;
        masked_matrix[mask_mat>mask_thres]=mask_val;
        return(masked_matrix);
}

for(model_ix in model_types){

	title_page(
		model_ix,
		paste("Responses: ", "", paste(response_no_reference, collapse="\n"), sep="\n"),
		paste(c("All Coefficients / P-Values Table", "All Coefficient Heatmap", "Significant Coefficients"),
			collapse="\n")
	)

	for(time_ix in time_str_ids){

		coef_mat=t(fit_info[[time_ix]][[model_ix]][["Coefficients"]]);
		pval_mat=t(fit_info[[time_ix]][[model_ix]][["pvalues"]]);

		pred_names_noIntc=setdiff(rownames(coef_mat), "(Intercept)");
		coef_mat=coef_mat[pred_names_noIntc,, drop=F];
		pval_mat=pval_mat[pred_names_noIntc,, drop=F];
		
		title=paste("Model: ", model_ix, "  Time: ", time_ix);

		plot_text(c(
			title,
			"",
			"Coefficients",
			capture.output(print(coef_mat)),
			"",
			"P-Values",
			capture.output(print(pval_mat))
		));

		paint_matrix(coef_mat, title=paste(title, " Coefficients (All)"));


		filt_coef_mat=mask_matrix(coef_mat, pval_mat, .1, 0);
		paint_matrix(filt_coef_mat, title=paste(title, " Coefficients (P-values < 0.1)"), label_zeros=F);

	}
}

###############################################################################


signf_mat=matrix(0, nrow=num_modeltypes, ncol=length(response_no_reference));
rownames(signf_mat)=model_types;
colnames(signf_mat)=response_no_reference;

signf_mat_wtime=matrix("", nrow=num_modeltypes, ncol=length(response_no_reference));
rownames(signf_mat_wtime)=model_types;
colnames(signf_mat_wtime)=response_no_reference;

for(model_ix in model_types){
	for(time_ix in time_str_ids){
		for(rsp_ix in response_no_reference){

			# cumulative
			num_sig=fit_info[[time_ix]][[model_ix]][["num_signif"]][rsp_ix];
			if(is.finite(num_sig)){
				signf_mat[model_ix, rsp_ix]=num_sig+signf_mat[model_ix, rsp_ix];
			}

			# by time
			if(signf_mat_wtime[model_ix, rsp_ix]==""){
				signf_str=num_sig;
			}else{
				signf_str=paste(signf_mat_wtime[model_ix, rsp_ix], ",", num_sig, sep="");
			}
			signf_mat_wtime[model_ix, rsp_ix]=signf_str;
		}
	}
}

total_signf=sum(signf_mat);

if(num_unique_time_ids>1){
	out_arr=c(
		paste("Total Significant Predictors: ", total_signf, sep=""),
		"",
		"Number of Significant Predictors:",
		capture.output(print(signf_mat)),
		"",
		"Number of Significant Predictors by Time (comma-separated):",
		capture.output(print(signf_mat_wtime, quote=F))
	);

	

}else{
	out_arr=c(
		paste("Total Significant Predictors: ", total_signf, sep=""),
		"",
		"Number of Significant Predictors:",
		capture.output(print(signf_mat))
	);
}

plot_text(out_arr);

# Output summary in text
fh=file(paste(OutputRoot, ".signf_summary.txt", sep=""), "w");
cat(file=fh, "\nRun Name: ", OutputRoot, "\n", sep="");	
cat(file=fh, "\n");	
cat(file=fh, paste(out_arr, collapse="\n"));
cat(file=fh, "\n\n");	
close(fh);


###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
