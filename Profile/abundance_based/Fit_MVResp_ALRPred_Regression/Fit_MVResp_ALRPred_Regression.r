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
	"factors", "f", 1, "character",
	"covariates", "c", 1, "character",
	"responses", "y", 1, "character",
	"required", "q", 2, "character",

	"num_variables", "p", 2, "numeric",
	"additional_categories", "a", 2, "character",
	"reference_levels", "r", 2, "character",
	"outputroot", "o", 2, "character",

	"contains_remaining", "R", 2, "logical",
	"shorten_category_names", "x", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

NUM_TOP_CAT=35;

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
	inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", quote="", row.names=1))
	counts_mat=inmat[,2:(ncol(inmat))];

	# Clean category names a little
	cat_names=colnames(counts_mat);
	cat_names=gsub("-", "_", cat_names);
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
			cat("Note: ", target_relev_name, " not in model.  Ignoring reference releveling.\n\n", sep="");
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

additive_log_rato=function(ordered_matrix){
# Assumes last column will be the denominator

	num_cat=ncol(ordered_matrix);
	num_samp=nrow(ordered_matrix);

	denominator=ordered_matrix[,num_cat];
	alr_mat=matrix(0, nrow=num_samp, ncol=(num_cat-1));
	
	for(i in 1:num_samp){
		alr_mat[i,]=log(ordered_matrix[i,1:(num_cat-1)]/denominator[i]);
		#print(alr_mat[i,]);
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
	cnames=colnames(coeff);
	pval_ix=which(cnames=="Pr(>|t|)");
	pval=coeff[,pval_ix];

	sig_char=function(val){
		if(val <= .0001){ return("***");}
		if(val <= .001 ){ return("** ");}
		if(val <= .01  ){ return("*  ");}
		if(val <= .05  ){ return(":  ");}
		if(val <= .1   ){ return(".  ");}
		return(" ");
	}

	sig_arr=sapply(pval, sig_char);
	fdr=round(p.adjust(pval, method="fdr"), 4);
	fdr_sig_char=sapply(fdr, sig_char);
	#out_mat=cbind(coeff, fdr);
	out_mat=cbind(coeff, sig_arr, fdr, fdr_sig_char);
	colnames(out_mat)=c(cnames, "Signf", "FDR", "Signf");
	
	return(formatC(out_mat, format="f", digits=5,));
	
}


plot_fit=function(fit, sumfit, i=1){
	par.orig=par(no.readonly=T);

	observed=as.matrix(fit$model)[,i, drop=F];
	predicted=as.matrix(fit$fitted.values)[,i, drop=F];

	adjsqrd=signif(sumfit[[i]]$adj.r.squared, 4);
	fstat=sumfit[[i]]$fstatistic;
	pval=signif(1-pf(fstat["value"], fstat["numdf"], fstat["dendf"]), 4);

	name=colnames(predicted);

	par(mar=c(15, 15, 15, 15));
	plot(observed, predicted, main="", xlab="Observed", ylab="Predicted");

	mtext(name, line=5.5, font=2, cex=3);
	mtext("Predicted vs. Observed", line=4, font=2, cex=1.1);
	mtext(paste("Adjusted R^2 = ", adjsqrd, sep=""),line=2, cex=.9);
	mtext(paste("Model F-stat p-value = ", pval, sep=""), line=1, cex=.9);

	abline(a=0, b=1, col="blue");

	par(par.orig);
}

##############################################################################
##############################################################################

pdf(paste(OutputRoot, ".mvr_alrp.pdf", sep=""), height=11, width=9.5);

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
	}
	colnames(counts)=short_names;
	cat("Names have been shortened.\n");
}

# Normalize
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

print(setdiff(covariates_arr, factor_names));
print(setdiff(responses_arr, factor_names));

overlapping_variables=intersect(responses_arr, covariates_arr);
if(length(overlapping_variables)!=0){
	cat("Error:  You have the same variable names in response and covariate:\n");
	print(overlapping_variables);	
	quit(-1);
}

kept_variables=union(responses_arr, covariates_arr);
kept_factors=factors[,kept_variables, drop=F];

# Relevel factor levels
if(ReferenceLevelsFile!=""){
        ref_lev_mat=load_reference_levels_file(ReferenceLevelsFile)
        factors=relevel_factors(kept_factors, ref_lev_mat);
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

##############################################################################

# Assign 0's to values smaller than smallest abundance across entire dataset
min_assay=min(normalized[normalized!=0]);
cat("Lowest non-zero value: ", min_assay, "\n", sep="");
zero_replacment=min_assay/10;
cat("Substituting 0's with: ", zero_replacment, "\n", sep="");
normalized[normalized==0]=zero_replacment;

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
resp_alr_struct=additive_log_rato(cat_abundances);
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

#print(response_factors);
response_factors=as.matrix(response_factors);
dafr_predictors_factors=as.data.frame(cbind(covariate_factors, alr_categories_val));

cat("Response Summary:\n");
s=summary(response_factors);
plot_text(c(
	"Response Summary:",
	"\n",
	capture.output(print(s))
));
cat("Plotting Response Histograms:\n");
print(response_factors);
plot_histograms(response_factors);

cor_mat=cor(response_factors);
par(oma=c(1,1,1,1));
paint_matrix(cor_mat, title="Response Correlations");

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
print(alr_categories_val);
plot_histograms(alr_categories_val);

###############################################################################
# Set up and run univariate and manova

covariate_formula_str=paste(covariates_arr, collapse=" + ");
alr_category_formula_str=paste(alr_cat_names, collapse=" + ");
cat("\n");
cat("Covariates: \n", covariate_formula_str, "\n", sep="");
cat("\n");
cat("ALR Predictors: \n", alr_category_formula_str, "\n", sep="");
cat("\n");

if(nchar(covariate_formula_str)){
	formula_str=paste("response_factors ~ ", covariate_formula_str, " + ", alr_category_formula_str, sep="");
	reduced_formula_str=paste("response_factors ~ ", covariate_formula_str, sep="");
}else{
	formula_str=paste("response_factors ~ ", alr_category_formula_str, sep="");
	reduced_formula_str=paste("response_factors ~ 1", sep="");
}

lmfit=lm(as.formula(formula_str), data=dafr_predictors_factors);
reduced_lmfit=lm(as.formula(reduced_formula_str), data=dafr_predictors_factors);

lm_summaries=summary(lmfit);
reduced_lm_summaries=summary(reduced_lmfit);

print(lmfit);
print(lm_summaries);

num_responses=length(responses_arr);

if(num_responses > 1 && (num_responses < lmfit$df.residual)){
	manova_res=anova(lmfit);
}else{
	manova_res=c();
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

###############################################################################
# Accumulate univariate results into matrix and generate heat map

num_coefficients=nrow(lm_summaries[[1]]$coefficients);
coefficient_names=rownames(lm_summaries[[1]]$coefficients);
summary_var_names=rownames(lm_summaries[[1]]$coefficients);
covariate_coefficients=setdiff(coefficient_names, c("(Intercept)", alr_cat_names));

summary_res_coeff=matrix(0, nrow=num_coefficients, ncol=num_responses, dimnames=list(coefficient_names, response_fit_names));
summary_res_pval=matrix(0, nrow=num_coefficients, ncol=num_responses, dimnames=list(coefficient_names, response_fit_names));
summary_res_rsqrd=matrix(0, nrow=num_responses, ncol=3, dimnames=list(response_fit_names, 
	c("Full Model", "Reduced Model", "Difference")));

for(i in 1:num_responses){
	cat("\n\nWorking on: ", responses[i], "\n");
	univar_summary=round(lm_summaries[[i]]$coefficients, 4);

	rsquared=signif(c(lm_summaries[[i]]$r.squared, lm_summaries[[i]]$adj.r.squared),4);
	reduced_rsquared=signif(c(reduced_lm_summaries[[i]]$r.squared, reduced_lm_summaries[[i]]$adj.r.squared),4);
	
	fstat=lm_summaries[[i]]$fstatistic;
	pval=1-pf(fstat[["value"]], fstat[["numdf"]], fstat[["dendf"]]);

	rsqrd_diff=rsquared[1]-reduced_rsquared[1];
	adj_rsqrd_diff=signif(rsquared[2]-reduced_rsquared[2], 4);

	plot_text(c(
		paste("Univariate: ", responses[i], sep=""),
		"Covariates portion of Full Model (Covariates + ALR):",
		"",
		capture.output(print(univar_summary[covariate_coefficients,,drop=F]))
	));

	plot_text(c(
		paste("Univariate: ", responses[i], sep=""),
		"ALR Predictors portion of Full Model (Covariates + ALR):",
		"\n",
		capture.output(print(add_sign_col(univar_summary[alr_cat_names,,drop=F]), quote=F))
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

summary_res_coeff=round(summary_res_coeff,2);
summary_res_pval=round(summary_res_pval,2);

cat("Coefficients Matrix:\n");
print(summary_res_coeff);
#par(oma=c(10,14,5,1));
if(length(covariate_coefficients)>0){
	paint_matrix(summary_res_coeff[covariate_coefficients,,drop=F], title="Covariate Coefficients",value.cex=2, deci_pts=2);
}

# Variations of ALR Predictor Coefficients
paint_matrix(summary_res_coeff[alr_cat_names,,drop=F], title="ALR Predictors Coefficients (By Decreasing Abundance)", 
	value.cex=2, deci_pts=2, plot_row_dendr=F, plot_col_dendr=F);
paint_matrix(summary_res_coeff[alr_cat_names,,drop=F], title="ALR Predictors Coefficients (ALR Clusters)", 
	value.cex=2, deci_pts=2, plot_row_dendr=T, plot_col_dendr=F);
paint_matrix(summary_res_coeff[alr_cat_names,,drop=F], title="ALR Predictors Coefficients (Response Clusters)", 
	value.cex=2, deci_pts=2, plot_row_dendr=F, plot_col_dendr=T);
paint_matrix(summary_res_coeff[alr_cat_names,,drop=F], title="ALR Predictors Coefficients (ALR and Response Clusters)", 
	value.cex=2, deci_pts=2, plot_row_dendr=T, plot_col_dendr=T);

cat("\nP-values:\n");
print(summary_res_pval);
#par(oma=c(10,14,5,1));

if(length(covariate_coefficients)>0){
	paint_matrix(summary_res_pval[covariate_coefficients,,drop=F], title="Covariate P-values", 
		plot_min=0, plot_max=1, high_is_hot=F, value.cex=2, deci_pts=2);
}

# Variations of ALR Predictor P-values
paint_matrix(summary_res_pval[alr_cat_names,,drop=F], title="ALR Predictors P-values (By Decreasing Abundance)", 
	plot_min=0, plot_max=1, high_is_hot=F, value.cex=2, deci_pts=2);

paint_matrix(summary_res_pval[alr_cat_names,,drop=F], title="ALR Predictors P-values (ALR Clusters)", 
	plot_min=0, plot_max=1, high_is_hot=F, value.cex=2, deci_pts=2,
	plot_row_dendr=T
);
paint_matrix(summary_res_pval[alr_cat_names,,drop=F], title="ALR Predictors P-values (Response Clusters)", 
	plot_min=0, plot_max=1, high_is_hot=F, value.cex=2, deci_pts=2,
	plot_col_dendr=T
);
paint_matrix(summary_res_pval[alr_cat_names,,drop=F], title="ALR Predictors P-values (ALR and Response Clusters)", 
	plot_min=0, plot_max=1, high_is_hot=F, value.cex=2, deci_pts=2,
	plot_row_dendr=T, plot_col_dendr=T
);

###############################################################################

paint_matrix(summary_res_rsqrd, title="Univariate Adjusted R-Squared");

if(length(manova_res)>0){
	print(manova_res);
	plot_text(c(
		"MANOVA",
		capture.output(print(manova_res))
	));

	manova_pval_mat=matrix(0, nrow=length(alr_cat_names), ncol=1, dimnames=list(alr_cat_names, "Pr(>F)"));
	manova_pval=manova_res[,"Pr(>F)", drop=F]
	manova_pval_mat[alr_cat_names,]=manova_pval[alr_cat_names,];
	print(manova_pval_mat);
	#par(oma=c(10,14,5,1));
	paint_matrix(manova_pval_mat[alr_cat_names,,drop=F], title="ALR Predictors MANOVA", 
		plot_min=0, plot_max=1, high_is_hot=F, value.cex=1, deci_pts=3);

}
###############################################################################

#print(cat_abundances);
#print(alr_cat_names);
#print(manova_pval_mat);
#print(mean_abund);

plot_rank_abund=function(abundances, pvals, range=c(1,10), title="", ylim=NULL){

	num_colors=50;
	color_arr=rainbow(num_colors, start=0, end=4/6);

        remap=function(in_val, in_range, out_range){
                in_prop=(in_val-in_range[1])/(in_range[2]-in_range[1])
                out_val=in_prop*(out_range[2]-out_range[1])+out_range[1];
                return(out_val);
        }

	if(is.null(ylim)){
		ylim=c(0, max(abundances)*1.1);
	}
	
	category_colors=remap(pvals, c(0,1), c(1, num_colors));
	category_names=names(abundances);

	range_arr=range[1]:range[2];

	proportion_rep=sum(abundances[range_arr]);
	num_cat=range[2]-range[1];
	mids=barplot(abundances[range_arr], col=color_arr[category_colors[range_arr]], names.arg="", 
		yaxt="n",
		main=title, ylim=ylim);
	mtext(paste("Percent Represented: ", round(proportion_rep*100,2), "%", sep=""), side=3, line=0);

	# Calculate labels sizes/positions
	spacing=mids[2]-mids[1];

	gpar=par();
	dims=gpar$usr; # xlim/xlim of plot
	cxy=gpar$cxy;  # character size in plot dimensions

	# Labels
	resize=min(1,30/num_cat);
	text(mids-cxy[1]*.75, -cxy[2]/2, category_names[range_arr], srt=-45, xpd=T, cex=resize, pos=4);

	# y axis
	axis(2, at=seq(0,1,.05), las=2)

	# Legend
	legend_colors=remap(c(.05, .1, .5, 1), c(0,1), c(1, num_colors));
	legend(dims[2]-(dims[2]-dims[1])/5, dims[4]-(dims[4]-dims[3])/5,
		legend=c("0.05", "0.10", "0.50", "1.00"),
		fill=color_arr[legend_colors],
		cex=.75,
		title="p-value colors:",
		bty="n"
	);
	
}

if(length(manova_res)>0){

	par(oma=c(1,1,1,1));
	par(mar=c(15,3,4,15));
	par(mfrow=c(2,1));

	max_mean_abund=max(mean_abund)*1.2;
	plot_rank_abund(mean_abund, manova_pval_mat, range=c(1,num_top_taxa), "All Top Categories", ylim=c(0, max_mean_abund));

	if(num_top_taxa>=10){
		plot_rank_abund(mean_abund, manova_pval_mat, range=c(1,10), "Top 10 Categories", ylim=c(0, max_mean_abund));
	}

	if(num_top_taxa>=20){
		plot_rank_abund(mean_abund, manova_pval_mat, range=c(1,20), "Top 20 Categories", ylim=c(0, max_mean_abund));
	}

	if(num_top_taxa>=40){
		plot_rank_abund(mean_abund, manova_pval_mat, range=c(1,40), "Top 40 Categories", ylim=c(0, max_mean_abund));
	}

	plot_rank_abund(mean_abund, manova_pval_mat, range=c(1,floor(num_top_taxa/2)), 
		"Top Half of Categories", ylim=c(0, max_mean_abund));
	plot_rank_abund(mean_abund, manova_pval_mat, range=c(floor(num_top_taxa/2)+1,num_top_taxa), 
		"Bottom Half of Categories", ylim=c(0, max_mean_abund));
}


##############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
