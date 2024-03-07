#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
library(car);
library(pscl);

source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');

options(useFancyQuotes=F);


params=c(
	"summary_file", "s", 1, "character",
	"factors", "f", 1, "character",
	"covariates", "c", 1, "character",
	"responses", "y", 1, "character",
	"required", "q", 2, "character",

	"num_variables", "p", 2, "numeric",
	"additional_categories", "a", 2, "character",
	"reference_levels", "r", 2, "character",
	"outputroot", "o", 2, "character",

	"contains_remaining", "R", 2, "logical",
	"shorten_category_names", "x", 2, "character",

	"tag_name", "t", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

NUM_TOP_CAT=30;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table for taxa/function (used as Predictor/X's)>\n",
	"	-f <factors file, contains covariates and multivariate Y>\n",
	"	-c <list of covariate X's names to select from factor file (filename)>\n",
	"	-y <list of response Y's names to select from factor file (filename)>\n",
	"	[-q <required list of variables to include after NA removal>]\n",
	"\n",
	"	[-p <number of top taxonomic/categorical variables, default=", NUM_TOP_CAT, ">]\n",
	"	[-a <list of additional categories (from summary table) to include (filename)>]\n",
	"	[-r <reference levels file for Y's in factor file>]\n",
	"	[-o <output filename root>]\n",
	"\n",
	"	[-R (pay attention to 'remaining' category)]\n",
	"	[-x <shorten category names, with separator in double quotes (default=\"\")>]\n",
	"\n",
	"	[-t <tag name>]\n",
	"\n",
	"This script will fit the following model:\n",
	"\n",
	" Multivariate Response = covariates + MALR(top-p taxonomy)\n",
	"\n",
	"If the -R flag is set, a 'remaining' category will be be included in the denominator\n",
	"	independent of how large it is.  I.e., do not use it as one of the response variables.\n",
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

if(length(opt$tag_name)){
	TagName=opt$tag_name;
	cat("Setting TagName Hook: ", TagName, "\n");
	setHook("plot.new", 
		function(){
			#cat("Hook called.\n");
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


SummaryFile=opt$summary_file;
FactorsFile=opt$factors;
CovariatesFile=opt$covariates;
ResponseFile=opt$responses;

cat("\n");
cat("   Summary File: ", SummaryFile, "\n", sep="");
cat("   Factors File: ", FactorsFile, "\n", sep="");
cat("Covariates File: ", CovariatesFile, "\n", sep="");
cat("  Response File: ", ResponseFile, "\n", sep="");
cat("  Required File: ", RequiredFile, "\n", sep="");
cat("Additional File: ", AdditionalCatFile, "\n", sep="");
cat("\n");
cat("Output File: ", OutputRoot, "\n", sep="");
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
	factors=data.frame(read.table(fname,  sep="\t", header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=TRUE));
	factor_names=colnames(factors);

	ignore_idx=grep("^IGNORE\\.", factor_names);

	if(length(ignore_idx)!=0){
		return(factors[-ignore_idx]);
	}else{
		return(factors);
	}
}

load_summary_file=function(fname){
	inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", quote="", row.names=1))
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
	par(oma=rep(1,4));
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

plot_fit=function(fit, sumfit, mod_stats, mod_anovas, resp_name){

	#print(names(fit));
	#print(names(sumfit));

	par.orig=par(no.readonly=T);

	observed=fit[["model"]][,resp_name];
	predicted=fit[["fitted.values"]];
	
	model_fam=fit$family[[1]];
	mcfadden_r2=mod_stats[["full"]][["mcfadden"]];

	#fstat=sumfit[[i]]$fstatistic;
	#model_pval=1-pf(fstat["value"], fstat["numdf"], fstat["dendf"]);

	#name=colnames(predicted);

	par(mar=c(15, 15, 15, 15));

	obs_pred_range=range(c(observed, predicted));
	span=diff(obs_pred_range);
	pad=span*.07;

	plot(observed, predicted, 
		xlim=c(obs_pred_range[1]-pad, obs_pred_range[2]+pad),
		ylim=c(obs_pred_range[1]-pad, obs_pred_range[2]+2*pad),
		main="", xlab="Observed", ylab="Predicted");

	abline(a=0, b=1, col="blue");
	points(lowess(observed, predicted), type="l", col="red");
	legend(obs_pred_range[1]+pad, obs_pred_range[2]+pad, legend=c("Model", "Ideal"), 
		fill=c("red", "blue"), bty="n");

	mtext(resp_name, line=8.5, font=2, cex=3);
	mtext("Predicted vs. Observed", line=7, font=2, cex=1.1);

	# Annotate with some stats
	mtext(sprintf("Model Family: %s", model_fam), line=6, cex=.8);
	mtext(sprintf("McFadden's R^2 = %4.3f", mcfadden_r2), line=5, cex=.9);
	mtext(sprintf("Full Model (vs. Null Model) P-val: %6.3g",
		mod_anovas[["null_full"]]), line=3, cex=.9);
	mtext(sprintf("ALR's Contrib to Model Improvement (Full vs. Reduced Model): P-val = %6.3g", 
		mod_anovas[["reduced_full"]]), line=2, cex=.9);
	mtext(sprintf("Covariates-Only Model (Reduced vs. Null Model): P-val = %6.3g", 
		mod_anovas[["null_reduced"]]), line=1, cex=.9);

	par(par.orig);
}

plot_predictor_barplot=function(coef_tab, resp_name){
	
	par.orig=par(no.readonly=T);
	print(coef_tab);

	# Remove intercept from table
	predictors=setdiff(rownames(coef_tab), "(Intercept)");
	coef_tab=coef_tab[predictors,,drop=F];

	# Keep only pred w/ pval<.1
	pval=coef_tab[,4,drop=F];
	abv_cutoff_ix=pval<0.1;
	coef_tab=coef_tab[abv_cutoff_ix,,drop=F];

	# Sort table by pval
	pval=coef_tab[,4,drop=F];
	order_inc_pval_ix=order(pval, decreasing=F);
	coef_tab=coef_tab[order_inc_pval_ix,,drop=F];

	# Pull out coef/est
	coef=coef_tab[,"Estimate"];
	pval=coef_tab[,4];
	pred=rownames(coef_tab);	
	cols=sapply(coef, function(x){ ifelse(x>0, "green", "red")});

	print(coef);
	print(pval);
	print(pred);



	ref_cutoffs=c(1, 0.1, 0.05, 0.01, 0.001);
	min_pval=min(pval);

	nlp_ref_cutoffs=-log10(ref_cutoffs);
	nlp_min_pval=-log10(min(min_pval, 0.001));
	nlp_pvals=-log10(pval);

	par(mfrow=c(1,1));
	par(mar=c(15,5,10,5));
	mids=barplot(nlp_pvals, names.arg="", main="", ylab="-Log10(p-value)",
		col=cols,
		ylim=c(0, nlp_min_pval));
	abline(h=nlp_ref_cutoffs, lty="dashed", col="blue");
	axis(1, at=mids, labels=pred, las=2);
	axis(4, at=nlp_ref_cutoffs, labels=ref_cutoffs, las=2);
	mtext("P-values", side=4, line=3.5, cex=.8);	

	title(main="most significant predictors", line=1.5, cex.main=1.5, font.main=1);
	title(main=resp_name, line=3, cex.main=3, font.main=2);

	par(par.orig);

}

##############################################################################
##############################################################################

pdf(paste(OutputRoot, ".alr_as_pred.pdf", sep=""), height=11, width=9.5);

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
responses_arr=load_list(ResponseFile);
cat("Multivariate Response Variables:\n");
print(responses_arr);
cat("\n");

covariates_arr=load_list(CovariatesFile);
cat("Covariates Variables:\n");
print(covariates_arr);
cat("\n");

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
	paste("   Summary File: ", SummaryFile, sep=""),
	paste("   Factors File: ", FactorsFile, sep=""),
	paste("Covariates File: ", CovariatesFile, sep=""),
	paste("  Response File: ", ResponseFile, sep=""),
	paste("  Required File: ", RequiredFile, sep=""),
	paste("Additional File: ", AdditionalCatFile, sep=""),
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
	capture.output(print(responses_arr)),
	"",
	"Covariates:",
	capture.output(print(covariates_arr)),
	"",
	"Required Variables:",
	 capture.output(print(required_arr))
));

if(length(responses_arr)==0){
	msg=c(
		"For ALR as a predictor, a list of response (y) variables needs to be provided.",
		"The Response file was provided, but was empty.",
		""
	);
	plot_text(msg);
	cat(msg, sep="\n");
	quit(status=0);
}


overlapping_variables=intersect(responses_arr, covariates_arr);
if(length(overlapping_variables)!=0){
	cat("Error:  You have the same variable names in response and covariate:\n");
	print(overlapping_variables);	
	quit(status=-1);
}

kept_variables=union(responses_arr, covariates_arr);
kept_factors=factors[,kept_variables, drop=F];

# Relevel factor levels
if(ReferenceLevelsFile!=""){
        ref_lev_mat=load_reference_levels_file(ReferenceLevelsFile)
        kept_factors=relevel_factors(kept_factors, ref_lev_mat);
}else{
        cat("No Reference Levels File specified.\n");
}

response_factors=kept_factors[,responses_arr, drop=F];
summary(response_factors);
covariate_factors=kept_factors[,covariates_arr, drop=F];
resp_cov_factors=cbind(response_factors, covariate_factors);

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

num_response_variables=length(responses_arr);
num_covariates_variables=length(covariates_arr);

##############################################################################

# Assign 0's to values smaller than smallest abundance across entire dataset
#min_assay=min(normalized[normalized!=0]);
#cat("Lowest non-zero value: ", min_assay, "\n", sep="");
#zero_replacment=min_assay/10;
#cat("Substituting 0's with: ", zero_replacment, "\n", sep="");
#normalized[normalized==0]=zero_replacment;

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
predictors_factors_mat=cbind(covariate_factors, alr_categories_val);

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

is_numeric_response_factors=logical();
for(i in 1:ncol(response_factors)){
	is_numeric_response_factors[i]=is.numeric(response_factors[,i]);
}
cor_mat=cor(response_factors[, is_numeric_response_factors, drop=F]);
par(oma=c(1,1,1,1));
paint_matrix(cor_mat, title="Response Correlations", value.cex=.7, plot_min=-1, plot_max=1 );

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

# Determine whether response is normal or binomial, or other?
cat("Determining whether response variable is continuous or boolean...\n");

assign_model_families=function(resp_var_mat){

	num_resp=ncol(resp_var_mat);

	mod_fam=character(num_resp);
	names(mod_fam)=colnames(resp_var_mat);

	for(i in 1:num_resp){
		cur_resp=resp_var_mat[,i];

		uniq_resp=unique(cur_resp);
		num_uniq_resp=length(uniq_resp);
		min_resp=min(uniq_resp);
		max_resp=max(uniq_resp);

		if(num_uniq_resp==2 && min_resp==0 && max_resp==1){
			mod_fam[i]="binomial";
		}else{
			mod_fam[i]="gaussian";
		}
	}

	return(mod_fam);
}

model_families=assign_model_families(response_factors);
num_bool_resp=sum(model_families=="binomial");
all_gaussian_resp=all(model_families=="gaussian");

cat("Model Families:\n");
print(model_families);

run_glm=F;
if(num_bool_resp){
	cat("Boolean responses identified.  Running GLM's for logistic regression.\n");
	run_glm=T;
}

###############################################################################
# Set up and run univariate and manova

covariate_formula_str=paste(covariates_arr, collapse=" + ");
alr_category_formula_str=paste(alr_cat_names, collapse=" + ");
cat("\n");
cat("Covariates: \n", covariate_formula_str, "\n", sep="");
cat("\n");
cat("ALR Predictors: \n", alr_category_formula_str, "\n", sep="");
cat("\n");

#if(nchar(covariate_formula_str)){
#	formula_str=paste("response_factors_dfr ~ ", covariate_formula_str, " + ", 
#		alr_category_formula_str, sep="");
#	reduced_formula_str=paste("response_factors_dfr ~ ", covariate_formula_str, sep="");
#}else{
#	formula_str=paste("response_factors_dfr ~ ", alr_category_formula_str, sep="");
#	reduced_formula_str=paste("response_factors_dfr ~ 1", sep="");
#}

if(length(covariates_arr)){
	# If there are covariates to include, add them to the predictors
	full_pred_str=paste(covariate_formula_str, " + ", alr_category_formula_str, sep="");
	reduced_pred_str=covariate_formula_str;
}else{
	# If there are no covariates, just have the alr categories
	full_pred_str=paste(alr_category_formula_str, sep="");
	reduced_pred_str="1";
}

cat("Full Model Predictors:\n");
print(full_pred_str);
cat("\n");
cat("Reduced Model Predictors:\n");
print(reduced_pred_str);
cat("\n");

response_factors_dfr=as.matrix(response_factors);

#--------------------------------------------------

num_pred_var=ncol(predictors_factors_mat);
num_samples=nrow(predictors_factors_mat);
cat("Number of Predictors (covariates + ALR): ", num_pred_var, "\n");
cat("Number of Samples: ", num_samples, "\n");

#------------------------------------------------------------------------------

compute_mod_stats=function(fit, sumfit){

	ms=list();

	ms[["mcfadden"]]=1-sumfit$deviance/sumfit$null.deviance;
	ms[["aic"]]=fit$aic;

	model_family=fit$family$family;
	if(model_family=="gaussian"){

		ss_res=sum(fit$residuals^2);
		obs_resp=fit$y;
		ss_tot=sum((obs_resp - mean(obs_resp))^2);
		unadj_r2=1-(ss_res/ss_tot);

		ms[["unadjusted"]]=unadj_r2;

		num_obs=length(obs_resp);
		num_pred_var=length(fit$coefficients)-1; # exclude intercept
		adj_r2=1-((1-unadj_r2)*(num_obs-1)/(num_obs-num_pred_var-1));

		ms[["adjusted"]]=adj_r2;

	}else{
		
		ms[["unadjusted"]]=NA;
		ms[["adjusted"]]=NA;
	}

	return(ms);

}

#------------------------------------------------------------------------------

if(run_glm){
	cat("Running glm...\n");

	resp_pred_mat=cbind(response_factors, predictors_factors_mat);

	glm_full_fit_list=list();
	glm_full_sumfit_list=list();

	glm_reduced_fit_list=list();
	glm_reduced_sumfit_list=list();

	glm_null_fit_list=list();
	glm_null_sumfit_list=list();

	mod_stats_list=list();
	model_anova_list=list();

	num_resp_var=ncol(response_factors);
	
	resp_names=colnames(response_factors);
	for(rsp_ix in 1:num_resp_var){
	
		cur_resp_name=resp_names[rsp_ix];
		cat("Fitting ", cur_resp_name, " Models:\n");


		full_form_str=paste(cur_resp_name, "~", full_pred_str);
		reduced_form_str=paste(cur_resp_name, "~", reduced_pred_str);
		null_form_str=paste(cur_resp_name, "~1");

		#---------------------------------------------------------------
		# Full Model

		glm_full_fit=glm(as.formula(full_form_str), data=resp_pred_mat, 
			family=model_families[cur_resp_name]);

		glm_full_sumfit=summary(glm_full_fit);

		#print(glm_full_fit);
		#print(glm_full_sumfit);

		#---------------------------------------------------------------
		# Reduced Model

		glm_reduced_fit=glm(as.formula(reduced_form_str), data=resp_pred_mat, 
			family=model_families[cur_resp_name]);

		glm_reduced_sumfit=summary(glm_reduced_fit);

		#print(glm_reduced_fit);
		#print(glm_reduced_sumfit);

		#---------------------------------------------------------------
		# Null Model

		glm_null_fit=glm(as.formula(null_form_str), data=resp_pred_mat, 
			family=model_families[cur_resp_name]);

		glm_null_sumfit=summary(glm_null_fit);

		#print(glm_null_fit);
		#print(glm_null_sumfit);

		#---------------------------------------------------------------
		# Calculate R^2
		full_r2s=compute_mod_stats(glm_full_fit, glm_full_sumfit);
		reduced_r2s=compute_mod_stats(glm_reduced_fit, glm_reduced_sumfit);

		mod_stat=list();
		mod_stat[["full"]]=full_r2s;
		mod_stat[["reduced"]]=reduced_r2s;

		#---------------------------------------------------------------
		# Calculate ANOVA
		anova_null_reduced_res=anova(glm_null_fit, glm_reduced_fit, test="Chi");
                anova_null_reduced_pval=anova_null_reduced_res[["Pr(>Chi)"]][2];

		anova_null_full_res=anova(glm_null_fit, glm_full_fit, test="Chi");
                anova_null_full_pval=anova_null_full_res[["Pr(>Chi)"]][2];

		anova_reduced_full_res=anova(glm_reduced_fit, glm_full_fit, test="Chi");
                anova_reduced_full_pval=anova_reduced_full_res[["Pr(>Chi)"]][2];

		#---------------------------------------------------------------
		# Save models to list

		glm_full_fit_list[[cur_resp_name]]=glm_full_fit;
		glm_full_sumfit_list[[cur_resp_name]]=glm_full_sumfit;
		glm_reduced_fit_list[[cur_resp_name]]=glm_reduced_fit;
		glm_reduced_sumfit_list[[cur_resp_name]]=glm_reduced_sumfit;
		glm_null_fit_list[[cur_resp_name]]=glm_null_fit;
		glm_null_sumfit_list[[cur_resp_name]]=glm_null_sumfit;

		mod_stats_list[[cur_resp_name]]=mod_stat;

		mod_anovas=list();
		mod_anovas[["null_reduced"]]=anova_null_reduced_pval;
		mod_anovas[["null_full"]]=anova_null_full_pval;
		mod_anovas[["reduced_full"]]=anova_reduced_full_pval;
		model_anova_list[[cur_resp_name]]=mod_anovas;


		cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	
	}

	# Copy all the residuals from the glm fit to a single matrix
	num_samples=nrow(response_factors);
	residuals_mat=matrix(NA, nrow=num_samples, ncol=num_resp_var);
	colnames(residuals_mat)=resp_names;
	for(rsp_ix in 1:num_resp_var){
		cur_resp_name=resp_names[rsp_ix];
		residuals_mat[,rsp_ix]=glm_full_fit_list[[cur_resp_name]]$residuals;
	}
	#print(residuals_mat);

	# Run manova
	residual_df=(glm_full_fit_list[[1]]$df.residual);
	cat("DF of Residuals: ", residual_df, "\n", sep="");
	if(all_gaussian_resp && num_resp_var > 1 && (num_resp_var < residual_df)){

		cat("Conditions met to running MANOVA...\n");
		tryCatch({

			manova_res=manova(
				as.formula(paste("residuals_mat~", full_pred_str)), 
				data=resp_pred_mat);
			manova_sumres=summary(manova_res);

		}, error=function(e){
			cat("Error: MANOVA was attempted but not successful.\n");
			print(e);
		});

		print(manova_res);
		print(summary(manova_res));

	}else{

		cat("Conditions unmet to running MANOVA...\n");
		manova_res=NA;
		manova_sumres=NA;
	}

}

# glm_full_fit_list
# glm_full_sumfit_list
# glm_reduced_fit_list
# glm_reduced_sumfit_list
# mod_stat_list
# null_anova_list
# manova_res
# manova_sumres 

###############################################################################

accumulate_model_improvement_stats=function(mod_stats_lst, model_anova_lst){

	# "mcfadden" "aic" "unadjusted" "adjusted"

	resp_vars=names(model_anova_lst);
	num_resp_vars=length(resp_vars);
	
	headers=c(
		"Reduced R2", "Full R2", "Diff R2",
		"Reduced Adj R2", "Full Adj R2", "Diff Adj R2",
		"Reduced McFadden", "Full McFadden", "Diff McFadden",
		"Reduced AIC", "Full AIC", "Diff AIC",
		"ANOVA Full-Reduced P-Val",
		"ANOVA Reduced-Null P-Val",
		"ANOVA Full-Null P-Val"
	);

	summary_stat_mat=matrix(NA, nrow=num_resp_vars, ncol=length(headers));
	colnames(summary_stat_mat)=headers;
	rownames(summary_stat_mat)=resp_vars;

	for(var in resp_vars){

		reduced_stats=mod_stats_lst[[var]][["reduced"]];
		full_stats=mod_stats_lst[[var]][["full"]];
		model_anova=model_anova_lst[[var]];

		summary_stat_mat[var,]=c(
			reduced_stats[["unadjusted"]], full_stats[["unadjusted"]], 			
				full_stats[["unadjusted"]]-reduced_stats[["unadjusted"]],

			reduced_stats[["adjusted"]], full_stats[["adjusted"]], 			
				full_stats[["adjusted"]]-reduced_stats[["adjusted"]],

			reduced_stats[["mcfadden"]], full_stats[["mcfadden"]], 			
				full_stats[["mcfadden"]]-reduced_stats[["mcfadden"]],

			reduced_stats[["aic"]], full_stats[["aic"]], 			
				full_stats[["aic"]]-reduced_stats[["aic"]],

			model_anova[["reduced_full"]],
			model_anova[["null_reduced"]],
			model_anova[["null_full"]]
		);
		
	}

	# Generate a formatted table for printing
	summary_stat_char_mat=matrix("", nrow=num_resp_vars, ncol=length(headers)+4);
	colnames(summary_stat_char_mat)=c(headers, 
		"Substantial AIC Improvement", 
		"Signf FR", 
		"Signf RN", 
		"Signf FN"
		);
	rownames(summary_stat_char_mat)=resp_vars;

	for(var in 1:9){
		summary_stat_char_mat[,var]=sprintf("%5.3f", summary_stat_mat[,var]);
	}
	for(var in 10:12){
		summary_stat_char_mat[,var]=sprintf("%8.1f", summary_stat_mat[,var]);
	}

	format_pval=function(x){
		res=ifelse(x>.0001,
			sprintf("%11.3f", x),
			sprintf("%11.03g", x));
		return(res);
	}

	summary_stat_char_mat[,"ANOVA Full-Reduced P-Val"]=
		format_pval(summary_stat_mat[,"ANOVA Full-Reduced P-Val"]);
	summary_stat_char_mat[,"ANOVA Reduced-Null P-Val"]=
		format_pval(summary_stat_mat[,"ANOVA Reduced-Null P-Val"]);
	summary_stat_char_mat[,"ANOVA Full-Null P-Val"]=
		format_pval(summary_stat_mat[,"ANOVA Full-Null P-Val"]);
	
	summary_stat_char_mat[,"Substantial AIC Improvement"]=
		ifelse(summary_stat_mat[,"Diff AIC"]<(-2), "*", "");

	summary_stat_char_mat[,"Signf FR"]=
		sapply(summary_stat_mat[,"ANOVA Full-Reduced P-Val"], sig_char);
	summary_stat_char_mat[,"Signf RN"]=
		sapply(summary_stat_mat[,"ANOVA Reduced-Null P-Val"], sig_char);
	summary_stat_char_mat[,"Signf FN"]=
		sapply(summary_stat_mat[,"ANOVA Full-Null P-Val"], sig_char);
		

	print(summary_stat_mat);
	print(summary_stat_char_mat, quote=F);

	# Save and return results
	results=list();
	results[["numeric"]]=summary_stat_mat;
	results[["character"]]=summary_stat_char_mat;

	return(results);


}

mod_imprv_mat=accumulate_model_improvement_stats(mod_stats_list, model_anova_list);

reformat_coef_tab=function(coef_tab){

	num_col=ncol(coef_tab);
	num_coef=nrow(coef_tab);

	out_tab=matrix(character(), nrow=num_coef, ncol=num_col+1);
	rownames(out_tab)=rownames(coef_tab);
	colnames(out_tab)=c(colnames(coef_tab), "Signf");

	for(i in 1:3){
		out_tab[,i]=sprintf("%8.04f", coef_tab[,i]);
	}

	out_tab[,4]=ifelse(coef_tab[,4]>0.0001,
		sprintf("%11.3f", coef_tab[,4]),
		sprintf("%11.03g", coef_tab[,4]));

	out_tab[,"Signf"]=sapply(coef_tab[,4], sig_char);

	return(out_tab);
}

###############################################################################
# Accumulate univariate results into matrix and generate heat map

# glm_full_fit_list
# glm_full_sumfit_list
# glm_reduced_fit_list
# glm_reduced_sumfit_list
# mod_stat_list
# null_anova_list
# reduced_anova_list
# manova_res
# manova_sumres

num_responses=length(glm_full_fit_list);
response_names=names(glm_full_fit_list);
example_coef_mat=glm_full_sumfit_list[[1]]$coefficients;

# Report individual fits
for(i in 1:num_responses){

	cur_response=response_names[i];
	cat("\n\nWorking on: ", cur_response, "\n");

	full_coef_tab=glm_full_sumfit_list[[i]]$coefficients;
	reduced_coef_tab=glm_reduced_sumfit_list[[i]]$coefficients;

	cat("Full Model Coefficients:\n");
	print(full_coef_tab);

	cat("Reduced Model Coefficients:\n");
	print(reduced_coef_tab);

	reduced_predictors=rownames(reduced_coef_tab);
	cat("Reduced Predictors:\n");
	print(reduced_predictors);

	full_coef_tab_char=reformat_coef_tab(full_coef_tab);
	reduced_coef_tab_char=reformat_coef_tab(reduced_coef_tab);

	alr_comp=setdiff(rownames(full_coef_tab), reduced_predictors);
	
	#print(full_coef_tab_char, quote=F);
	#print(reduced_coef_tab_char, quote=F);

	mod_stats=mod_stats_list[[cur_response]];
	print(mod_stats);

	mod_anovas=model_anova_list[[cur_response]];
	cat("Reduced ANOVA: ", mod_anovas[["reduced_full"]], "\n");
	cat("Null ANOVA: ", mod_anovas[["null_full"]], "\n");
	cat("Cov ANOVA: ", mod_anovas[["null_reduced"]], "\n");

	full_mcfadden=mod_stats[["full"]][["mcfadden"]];
	reduced_mcfadden=mod_stats[["reduced"]][["mcfadden"]];

	full_aic=mod_stats[["full"]][["aic"]];
	reduced_aic=mod_stats[["reduced"]][["aic"]];
	aic_diff=reduced_aic-full_aic;

	#----------------------------------------------------------------------

	full_reduced_coef_comp=c(
		paste("Univariate Response: ", cur_response, sep=""),
		"",
		"Comparison of the Covariates portion of w/ and w/o ALR:",
		"",
		"---------------------------------------------------------------------------------------",
		"",
		"Predictors of Reduced Model (Covariates-only Model):",
		sprintf("McFadden's R^2: %4.3f, AIC: %.1f", reduced_mcfadden, reduced_aic),
		"",
		capture.output(print(reduced_coef_tab_char, quote=F)),
		"",
		"---------------------------------------------------------------------------------------",
		"",

		"Predictors Portion of Full Model (Covariates + ALR Model):",
		"  (Contribution of Covariates when controlling for ALR Components.)",
		sprintf("McFadden's R^2: %4.3f, AIC: %.1f", full_mcfadden, full_aic),
		"",
		capture.output(print(full_coef_tab_char[reduced_predictors,,drop=F], quote=F)),
		"",
		"(See next page for ALR components of Full Model.)",
		"---------------------------------------------------------------------------------------",
		"",
		paste("Model Improvement ANOVA (Full vs. Reduced): ", 
			sprintf("%6.02g", mod_anovas[["reduced_full"]]), 
			" ", sig_char(mod_anovas[["reduced_full"]]), sep=""),
		"",
		paste("Diff. in AIC: ", round(reduced_aic,2), " [reduced] - ", round(full_aic,2), 
			" [full] = ", round(aic_diff,2), " [dAIC]", sep=""),
		"If dAIC is > 2, then full model is substantially better than reduced model.",
		paste("Full model is ", 
			ifelse(aic_diff>2, "", "NOT "), 
			"substantially better then reduced model.", sep=""),
		"",
		sprintf("McFadden's R^2: %4.3f [full] - %4.3f [reduced] = %4.3f [delta]",
			full_mcfadden, reduced_mcfadden, (full_mcfadden-reduced_mcfadden), sep="")	
		
	);

	print(full_reduced_coef_comp, quote=F);
	cat("\n");
	plot_text(full_reduced_coef_comp);

	#----------------------------------------------------------------------

	alr_components_coef_comp=c(
		paste("Univariate Response: ", cur_response, sep=""),
		"",
		"ALR Predictors portion of Full Model (Covariates + ALR Model):",
		"",
		capture.output(print(full_coef_tab_char[alr_comp,,drop=F], quote=F))
	);

	print(alr_components_coef_comp, quote=F);
	plot_text(alr_components_coef_comp);

	#----------------------------------------------------------------------

	plot_predictor_barplot(glm_full_sumfit_list[[i]]$coefficients, cur_response);

	plot_fit(glm_full_fit_list[[i]], glm_full_sumfit_list[[i]], 
		mod_stats, mod_anovas, cur_response);

}

###############################################################################

cat("Accumulating coefficients and p-values into matrices.\n");

# Some times coefficients have been recoded by R into dummy variables, so we can't
#   just use the predictor variables names.

num_coefficients=nrow(example_coef_mat);
coefficient_names=rownames(example_coef_mat);

covariate_coefficients=setdiff(coefficient_names, c("(Intercept)", alr_cat_names));
num_covariate_coefficients=length(covariate_coefficients);

# In case, some alr categories are not calculable...
shrd_alr_names=intersect(alr_cat_names, coefficient_names);

cat("Num Coefficients: ", num_coefficients, "\n");
cat("Num Response Variables: ", num_responses, "\n");
cat("Num Covariate Predictors: ", num_covariate_coefficients, "\n");

# Coef and Pval for full model
summary_res_coef=matrix(0, nrow=num_coefficients, ncol=num_responses, 
	dimnames=list(coefficient_names, response_names));
summary_res_pval=matrix(0, nrow=num_coefficients, ncol=num_responses, 
	dimnames=list(coefficient_names, response_names));

# Coef and Pval for reduced model
summary_res_reduced_coef=matrix(0, nrow=num_covariate_coefficients, ncol=num_responses,
	dimnames=list(covariate_coefficients, response_names));
summary_res_reduced_pval=matrix(0, nrow=num_covariate_coefficients, ncol=num_responses,
	dimnames=list(covariate_coefficients, response_names));

# R^2 Table
summary_res_rsqrd=matrix(0, nrow=num_responses, ncol=3,
	dimnames=list(response_names, c("Reduced R^2", "Full R^2", "Difference")));

for(i in 1:num_responses){

	cur_response=response_names[i];

	# Full model coef/pvals
	fit_coef=glm_full_sumfit_list[[i]][["coefficients"]];
	pred=rownames(fit_coef);
	summary_res_coef[pred,cur_response]=fit_coef[pred,"Estimate"];
	summary_res_pval[pred,cur_response]=fit_coef[pred, 4]; # Pr(>|z|) or Pr(>|t|)

	# Reduce/Full Model Mcfaddens's
	mod_stats=mod_stats_list[[cur_response]];
	summary_res_rsqrd[cur_response,"Reduced R^2"]=mod_stats[["reduced"]][["mcfadden"]];
	summary_res_rsqrd[cur_response,"Full R^2"]=mod_stats[["full"]][["mcfadden"]];
	summary_res_rsqrd[cur_response,"Difference"]=
		mod_stats[["full"]][["mcfadden"]]-mod_stats[["reduced"]][["mcfadden"]];

}

# Remove (Intercept)
keep_coef=setdiff(rownames(summary_res_coef), "(Intercept)");
summary_res_coef=summary_res_coef[keep_coef,,drop=F];
summary_res_pval=summary_res_pval[keep_coef,,drop=F];

cat("Coefficients:\n");
print(summary_res_coef);
cat("\n");
cat("P-values:\n");
print(summary_res_pval);
cat("\n");
cat("McF's R^2:\n");
print(summary_res_rsqrd);

###############################################################################

mask_matrix=function(val_mat, mask_mat, mask_thres, mask_val){
        masked_matrix=val_mat;
        masked_matrix[mask_mat>mask_thres]=mask_val;
        return(masked_matrix);
}

#------------------------------------------------------------------------------

#summary_res_coeff=round(summary_res_coeff,2);
summary_res_pval_rnd=round(summary_res_pval,2);

# Plot covariates from full model
if(num_covariate_coefficients>0){

	cov_coef_mat=summary_res_coef[covariate_coefficients,,drop=F];
	cov_pval_mat=summary_res_pval[covariate_coefficients,,drop=F];

	paint_matrix(cov_coef_mat,
		title="Covariate Coefficients",value.cex=2, deci_pts=2);
	paint_matrix(cov_pval_mat,
		title="Covariate P-values",value.cex=2, deci_pts=2,
		plot_min=0, plot_max=1, high_is_hot=F);

	cov_coef_mat_100=mask_matrix(cov_coef_mat, cov_pval_mat, 0.10, 0);
	cov_coef_mat_050=mask_matrix(cov_coef_mat, cov_pval_mat, 0.05, 0);
	cov_coef_mat_010=mask_matrix(cov_coef_mat, cov_pval_mat, 0.01, 0);
	cov_coef_mat_001=mask_matrix(cov_coef_mat, cov_pval_mat, 0.001, 0);

	paint_matrix(cov_coef_mat_100, title="Covariate Coefficients: p-val < 0.10",
		value.cex=2, deci_pts=2, label_zeros=F);
	paint_matrix(cov_coef_mat_050, title="Covariate Coefficients: p-val < 0.05",
		value.cex=2, deci_pts=2, label_zeros=F);
	paint_matrix(cov_coef_mat_010, title="Covariate Coefficients: p-val < 0.01",
		value.cex=2, deci_pts=2, label_zeros=F);
	paint_matrix(cov_coef_mat_001, title="Covariate Coefficients: p-val < 0.001",
		value.cex=2, deci_pts=2, label_zeros=F);
}

# Variations of ALR Predictor Coefficients
paint_matrix(summary_res_coef[shrd_alr_names,,drop=F], 
	title="ALR Predictors Coefficients (By Decreasing Abundance)", 
	value.cex=1, deci_pts=2, plot_row_dendr=F, plot_col_dendr=F);
paint_matrix(summary_res_coef[shrd_alr_names,,drop=F], 
	title="ALR Predictors Coefficients (ALR Clusters)", 
	value.cex=1, deci_pts=2, plot_row_dendr=T, plot_col_dendr=F);
paint_matrix(summary_res_coef[shrd_alr_names,,drop=F], 
	title="ALR Predictors Coefficients (Response Clusters)", 
	value.cex=1, deci_pts=2, plot_row_dendr=F, plot_col_dendr=T);
paint_matrix(summary_res_coef[shrd_alr_names,,drop=F], 
	title="ALR Predictors Coefficients (ALR and Response Clusters)", 
	value.cex=1, deci_pts=2, plot_row_dendr=T, plot_col_dendr=T);


# Variations of ALR Predictor P-values
paint_matrix(summary_res_pval_rnd[shrd_alr_names,,drop=F], 
	title="ALR Predictors P-values (By Decreasing Abundance)", 
	plot_min=0, plot_max=1, high_is_hot=F, value.cex=1, deci_pts=2);


###############################################################################
# Mask coefficients at various pvalue
signf_coef=mask_matrix(summary_res_coef[shrd_alr_names,,drop=F], 
	summary_res_pval_rnd[shrd_alr_names,,drop=F], .10, 0);
paint_matrix(signf_coef, title="Significant ALR Predictors Coefficients (p-value < .10)", 
	value.cex=1, deci_pts=2, plot_row_dendr=F, plot_col_dendr=F, label_zeros=F);

signf_coef=mask_matrix(summary_res_coef[shrd_alr_names,,drop=F], 
	summary_res_pval_rnd[shrd_alr_names,,drop=F], .05, 0);
paint_matrix(signf_coef, title="Significant ALR Predictors Coefficients (p-value < .05)", 
	value.cex=1, deci_pts=2, plot_row_dendr=F, plot_col_dendr=F, label_zeros=F);

signf_coef=mask_matrix(summary_res_coef[shrd_alr_names,,drop=F], 
	summary_res_pval_rnd[shrd_alr_names,,drop=F], .01, 0);
paint_matrix(signf_coef, title="Significant ALR Predictors Coefficients (p-value < .01)", 
	value.cex=1, deci_pts=2, plot_row_dendr=F, plot_col_dendr=F, label_zeros=F);


paint_matrix(summary_res_pval_rnd[shrd_alr_names,,drop=F], 
	title="ALR Predictors P-values (ALR Clusters)", 
	plot_min=0, plot_max=1, high_is_hot=F, value.cex=1, deci_pts=2,
	plot_row_dendr=T
);
paint_matrix(summary_res_pval_rnd[shrd_alr_names,,drop=F], 
	title="ALR Predictors P-values (Response Clusters)", 
	plot_min=0, plot_max=1, high_is_hot=F, value.cex=1, deci_pts=2,
	plot_col_dendr=T
);
paint_matrix(summary_res_pval_rnd[shrd_alr_names,,drop=F], 
	title="ALR Predictors P-values (ALR and Response Clusters)", 
	plot_min=0, plot_max=1, high_is_hot=F, value.cex=1, deci_pts=2,
	plot_row_dendr=T, plot_col_dendr=T
);

###############################################################################

paint_matrix(summary_res_rsqrd, title="Univariate McFadden's R-Squared");

###############################################################################
# Report model comparison statistics

#print(summary_res_rsqrd);

brkt=function(mat){
	cnames=colnames(mat);
	colnames(mat)=sapply(cnames, function(x){ paste("[", x, "]", sep="");});
	return(mat);	
}

print(mod_imprv_mat);

r2mcf_mat=mod_imprv_mat[["character"]][,
	c("Reduced McFadden","Full McFadden","Diff McFadden",
	"ANOVA Full-Reduced P-Val","Signf FR"),drop=F];
r2_mat=mod_imprv_mat[["character"]][,
	c("Reduced R2","Full R2","Diff R2",
	"ANOVA Full-Reduced P-Val","Signf FR"),drop=F];
r2adj_mat=mod_imprv_mat[["character"]][,
	c("Reduced Adj R2","Full Adj R2","Diff Adj R2",
	"ANOVA Full-Reduced P-Val","Signf FR"),drop=F];
aic_mat=mod_imprv_mat[["character"]][,
	c("Reduced AIC","Full AIC","Diff AIC","Substantial AIC Improvement"),drop=F];

allmod_mat=mod_imprv_mat[["character"]][,
	c("Reduced McFadden","Full McFadden", "Diff McFadden",
	"ANOVA Reduced-Null P-Val","Signf RN",
	"ANOVA Full-Reduced P-Val","Signf FR",
	"ANOVA Full-Null P-Val","Signf FN"),
	drop=F];



fn_signf_ix=order(mod_imprv_mat[["numeric"]][,"ANOVA Full-Null P-Val"], decreasing=F);
rn_signf_ix=order(mod_imprv_mat[["numeric"]][,"ANOVA Reduced-Null P-Val"], decreasing=F);
fr_signf_ix=order(mod_imprv_mat[["numeric"]][,"ANOVA Full-Reduced P-Val"], decreasing=F);

plot_text(c(
	"Model Comparisons Using McFadden's R^2",
	"  (Use for comparing GLMs.)",
	"",
	capture.output(print(brkt(r2mcf_mat), quote=F))
	));

plot_text(c(
	"Model Comparisons Using Standard Unadjusted R^2",
	"  (Use for comparing LMs, NAs for binary responses.)",
	"",
	capture.output(print(brkt(r2_mat), quote=F))
	));

plot_text(c(
	"Model Comparisons Using Standard Adjusted R^2",
	"  (Use for comparing LMs, NAs for binary responses.)",
	"",
	capture.output(print(brkt(r2adj_mat), quote=F))
	));

plot_text(c(
	"Model Comparisons Using AIC",
	"  (Use for comparing GLMs and LMs.)",
	"  (Model improvements are consider significant when diff >2)",
	"",
	capture.output(print(brkt(aic_mat), quote=F))
	));

#------------------------------------------------------------------------------

FN_mat=allmod_mat[fn_signf_ix,c("Full McFadden","ANOVA Full-Null P-Val","Signf FN"),drop=F];
RN_mat=allmod_mat[rn_signf_ix,c("Reduced McFadden","ANOVA Reduced-Null P-Val","Signf RN"),drop=F];
FR_mat=allmod_mat[fr_signf_ix,c("Diff McFadden","ANOVA Full-Reduced P-Val","Signf FR"),drop=F];

plot_text(c(
	"Comparison of Model Strengths:",
	"  (Sorted by Strongest Full Model: Full vs. Null)",
	"",
	capture.output(print(brkt(FN_mat), quote=F))
	));

plot_text(c(
	"Comparison of Model Strengths",
	"  (Sorted by Strongest Covariates-Only Model: Reduced vs. Null)",
	"",
	capture.output(print(brkt(RN_mat), quote=F))
	));

plot_text(c(
	"Comparison of Model Strengths",
	"  (Sorted by Strongest ALR Improvement:  Full vs. Reduced)",
	"",
	capture.output(print(brkt(FR_mat), quote=F))
	));

#------------------------------------------------------------------------------

pull_abbrev_tab=function(full_coef_tab){

	#cat("Original Coef Tab:\n");
	#print(full_coef_tab);

	# Remove intercept
	predictors=setdiff(rownames(full_coef_tab), "(Intercept)");
	cur_coef_tab=full_coef_tab[predictors,,drop=F];

	# Remove if pval>.1
	signf_ix=(cur_coef_tab[,4] < 0.1);
	cur_coef_tab=cur_coef_tab[signf_ix, c(1,4), drop=F];

	# Sort
	signf_ord_ix=order(cur_coef_tab[,2], decreasing=F);
	cur_coef_tab=cur_coef_tab[signf_ord_ix,,drop=F];

	return(cur_coef_tab);

}

draw_list=function(abbr_tab, coloffset){

	chw=par()$cxy[1]/1.1; # width is over estimated
	chh=par()$cxy[2];	
	
	ch_rows=1/chh;
	ch_cols=1/chw;

	#cat("Rows: ", ch_rows, " / Cols: ", ch_cols, "\n");

	# Include cutoffs separators into table
	cutoffs=c(-1, 0.001, 0.01, 0.05);
	cutoff_entries=cbind(0, cutoffs);
	rownames(cutoff_entries)=c(" < 0.001 :", " < 0.01 :", " < 0.05 :", " < 0.1 :");
	abbr_tab_wcutoffs=rbind(abbr_tab, cutoff_entries);
	order_ix=order(abbr_tab_wcutoffs[,2]);
	abbr_tab_wcutoffs=abbr_tab_wcutoffs[order_ix,,drop=F];

	print(abbr_tab_wcutoffs);

	# Calculate font size adjustments based on number of items in list
	num_items=nrow(abbr_tab_wcutoffs);
	vadj=1;
	# If too many items, reduce size
	if(num_items>ch_rows){
		ch_rows=num_items;
		vadj=num_items/ch_rows;
	}

	item_names=rownames(abbr_tab_wcutoffs);

	row_id=0;
	for(i in 1:num_items){

		item_len=nchar(item_names[i]);
		hadj=1;
		# If size was adjust and the item is still to long, apply another adjustment
		if(vadj*item_len*chw > 1){
			hadj=1/(vadj*item_len*chw);
		}

		estimate=abbr_tab_wcutoffs[i,"Estimate"];
		if(estimate>0){
			draw_bullets=T;
			est_bgcol="green";
			est_sycol="black";
			est_ch="+";
			est_labfont=1;
			est_labcol="black";
			co_adj=1;
		}else if(estimate<0){
			draw_bullets=T;
			est_bgcol="red";
			est_sycol="white";
			est_ch="-";
			est_labfont=1;
			est_labcol="black";
			co_adj=1;
		}else{
			draw_bullets=F;
			est_ch="";
			est_labfont=3;
			est_labcol="blue";
			co_adj=.7;
		}


		# Mark +/- bullet glyphs
		if(draw_bullets){
			points(coloffset+chw*.6, 1-chh*row_id, pch=21, cex=1.3*vadj, 
				col="black", bg=est_bgcol);
			points(coloffset+chw*.6, 1-chh*row_id, pch=est_ch, cex=1*vadj, 
				col=est_sycol);
		}

		# Draw response name
		text(x=coloffset+chw*.6, y=1-chh*row_id, item_names[i], cex=hadj*vadj*co_adj, 
			pos=4, font=est_labfont, col=est_labcol);

		row_id=row_id+1;
	}

	
	
}

plot_summary_lists=function(sumfit_list, columns_pp=4){

	par.orig=par(no.readonly=T);

	resp_names=names(sumfit_list);
	
	cat("Columns per page: ", columns_pp, "\n");


	par(mar=c(.5, .5, 3, .5));
	par(mfrow=c(2,1));

	
	resp_count=0;
	for(cur_res in resp_names){

		resp_col=resp_count %% columns_pp;

		if(resp_col==0){
			plot(0, type="n", xlim=c(0,columns_pp), ylim=c(0,1), bty="n",
				xaxt="n", yaxt="n", xlab="", ylab="");
			abline(v=seq(0,columns_pp));
		}

		cat("Current Response: ", cur_res, "\n");

		# Label response name over column
		mtext(cur_res, side=3, line=.1, at=resp_col+.5, cex=1.1, font=2);

		# Extract abbreviated tables from coefficients in sumfit record
		abtab=pull_abbrev_tab(sumfit_list[[cur_res]]$coefficients);

		# Print the response's predictor in the specified column
		draw_list(abtab, resp_col);
		
		resp_count=resp_count+1;
	}
	
	par(par.orig);

}

plot_summary_lists(glm_full_sumfit_list, 5);

quit();

###############################################################################


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

fh=file(paste(OutputRoot, ".alr_as_pred.all_pred.coef_and_pval.tsv", sep=""), "w");
num_out_col=ncol(summary_res_coef);
response_names=colnames(summary_res_coef);

cat(file=fh, "Grouping:\t", paste(rep(OutputRoot, num_out_col), collapse="\t"), sep="");
cat(file=fh, "\n");
cat(file=fh, "Responses:\t", paste(response_names, collapse="\t"), sep="");
cat(file=fh, "\n\n");

# Coefficients
cat(file=fh, "[Coefficients]\n");
for(var in covariate_coefficients){
	cat(file=fh, var, paste(
		round(summary_res_coef[var, response_names], 4), collapse="\t"), sep="\t");
	cat(file=fh, "\n");
}

cat(file=fh, "\n");

for(var in shrd_alr_names){
	cat(file=fh, var, paste(
		round(summary_res_coef[var, response_names], 4), collapse="\t"), sep="\t");
	cat(file=fh, "\n");
}

cat(file=fh, "\n\n");

# P-values
cat(file=fh, "[P-values]\n");
for(var in covariate_coefficients){
	cat(file=fh, var, paste(
		round(summary_res_pval[var, response_names], 4), collapse="\t"), sep="\t");
	cat(file=fh, "\n");
}

cat(file=fh, "\n");

for(var in shrd_alr_names){
	cat(file=fh, var, paste(
		round(summary_res_pval[var, response_names], 4), collapse="\t"), sep="\t");
	cat(file=fh, "\n");
}

close(fh);

###############################################################################
# Write Reduced Model (covariates only) to file

fn=paste(OutputRoot, ".alr_as_pred.reduced.cov_only.coef_and_pval.tsv", sep="");

fh=file(fn, "w");
cat(file=fh, "Reduced Model (Covariates-only) Regression Results:\n\n");
cat(file=fh, "[Coefficients]\n\nCovariates\t");
close(fh);

write.table(file=fn, x=summary_res_reduced_coef, append=T, quote=F, sep="\t");

fh=file(fn, "a");
cat(file=fh, "\n\n[P-values]\n\nCovariates\t");
close(fh);

write.table(file=fn, x=summary_res_reduced_pval, append=T, quote=F, sep="\t");


###############################################################################
# Write ALR p-values to file
# Format: Factors as columns, taxa (predictor) as rows

exp_tab=summary_res_pval[shrd_alr_names,,drop=F];
write.table(exp_tab,  file=paste(OutputRoot, ".alr_as_pred.pvals.tsv", sep=""), sep="\t", quote=F, col.names=NA, row.names=T);

exp_tab=summary_res_coef[shrd_alr_names,,drop=F];
write.table(exp_tab,  file=paste(OutputRoot, ".alr_as_pred.coefs.tsv", sep=""), sep="\t", quote=F, col.names=NA, row.names=T);

###############################################################################
# Write R^2 Reduced/Full Model stats

write.table(improv_mat, file=paste(OutputRoot, ".alr_as_pred.rsqrd.tsv", sep=""), sep="\t", quote=F, col.names=NA, row.names=T);

###############################################################################
# Write MANOVA to file

if(length(manova_res)){
	num_variables=nrow(manova_res)-1;
	outmat=matrix("", nrow=num_variables, ncol=3);
	colnames(outmat)=c(TagName, "Pr(>F)", "Signf");
	varnames=unlist(rownames(manova_res));
	pvals=unlist(manova_res["Pr(>F)"]);
	outmat[,TagName]=varnames[1:num_variables];
	outmat[,"Pr(>F)"]=sprintf("%4.4f", pvals[1:num_variables]);
	outmat[,"Signf"]=sapply(pvals[1:num_variables], sig_char); 
}else{
	outmat=matrix("-", nrow=1, ncol=2);
	colnames(outmat)=c(TagName, "Pr(>F)");
}
write.table(outmat, file=paste(OutputRoot, ".alr_as_pred.anova.summary.tsv", sep=""), sep="\t", quote=F, col.names=T, row.names=F);

###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
