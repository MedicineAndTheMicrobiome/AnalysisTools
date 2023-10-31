#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
library(car);

source('~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r');

options(useFancyQuotes=F);


params=c(
	"summary_file", "s", 1, "character",

	"pairings", "p", 1, "character",
	"factors", "f", 1, "character",
	"factor_samp_id_name", "F", 1, "character",
	"model_var", "M", 1, "character",
	"required", "q", 2, "character",
	"response", "e", 1, "character",
	"predictor", "g", 1, "character",

	"num_resp_var", "u", 2, "numeric",
	"num_pred_var", "v", 2, "numeric",
	"alr_list_file", "a", 2, "character",

	"reference_levels", "c", 2, "character",
	"outputroot", "o", 2, "character",

	"contains_remaining", "R", 2, "logical",
	"shorten_category_names", "x", 2, "character",
	
	"tag_name", "t", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

NUM_TOP_RESP_CAT=35;
NUM_TOP_PRED_CAT=10;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table>\n",
	"\n",
	"	-p <pairings map, pairing Resp and Pred sample IDs. Must have header/column names>\n",
	"	-f <factors file, contains covariates and factors>\n",
	"       -F <column name of sample ids in factor file>\n",
	"	-M <list of covariate X's names to include in the model from the factor file>\n",
	"	-e <response ALR name, (column name in pairings file)\n",
	"	-g <predictor ALR name, (column name in pairings file)\n",
	"	[-q <required list of variables to include after NA removal>]\n",
	"\n",
	"	[-u <number of top response categories to analyze, default=", NUM_TOP_RESP_CAT, ">]\n",
	"	[-v <number of top predictor (as ALR) categories to include, default=", NUM_TOP_PRED_CAT, ">]\n",
	"	[-a <list of ALR categories to use in additon to top>]\n",
	"\n",
	"	[-o <output filename root>]\n",
	"\n",
	"	[-c <reference levels file for Y's in factor file>]\n",
	"	[-R (pay attention to 'remaining' category)]\n",
	"	[-x <shorten category names, with separator in double quotes (default=\"\")>]\n",
	"\n",
	"	[-t <tag name>]\n",
	"\n",
	"This script will fit the following 2 models for each category, i = 1 to ", NUM_TOP_RESP_CAT, ":\n",
	"  where the top P=", NUM_TOP_PRED_CAT, " categories are included.\n",
	"\n",
	"	1.) ALR_Resp[i] = ALR_Pred[i] + ALR_Pred[top-p] + covariates\n",
	"\n",
	"If the -R flag is set, a 'remaining' category will be be included in the denominator\n",
	"	independent of how large it is.  I.e., do not use it as one of the response variables.\n",
	"\n", sep="");

if(!length(opt$summary_file) || !length(opt$factors) || !length(opt$model_var) || 
	 !length(opt$response) || !length(opt$pairings) ){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputRoot=gsub(".summary_table.xls", "", opt$summary_file);
	OutputRoot=gsub(".summary_table.tsv", "", OutputRoot);
}else{
	OutputRoot=opt$outputroot;
}

if(!length(opt$num_resp_var)){
	NumRespVariables=NUM_TOP_RESP_CAT;
}else{
	NumRespVariables=opt$num_resp_var;
}

if(!length(opt$num_pred_var)){
	NumPredVariables=NUM_TOP_PRED_CAT;
}else{
	NumPredVariables=opt$num_pred_var;
}

NumMaxALRVariables=max(NumRespVariables, NumPredVariables);

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

if(length(opt$alr_list_file)){
	ALRCategListFile=opt$alr_list_file;
}else{
	ALRCategListFile="";
}

if(length(opt$factor_samp_id_name)){
        FactorSampleIDName=opt$factor_samp_id_name;
}else{
        FactorSampleIDName=1;
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
ModelVarFile=opt$model_var;
PairingsFile=opt$pairings;
ResponseName=opt$response;
PredictorName=opt$predictor;



OutputRoot=paste(OutputRoot, ".p_", PredictorName, ".r_", ResponseName, sep="");

cat("\n");
cat("         Summary File: ", SummaryFile, "\n", sep="");
cat("         Factors File: ", FactorsFile, "\n", sep="");
cat(" Model Variables File: ", ModelVarFile, "\n", sep="");
cat("        Pairings File: ", PairingsFile, "\n", sep="");
cat("        Response Name: ", ResponseName, "\n", sep="");
cat("       Predictor Name: ", PredictorName, "\n", sep="");
cat("          Output File: ", OutputRoot, "\n", sep="");
cat("\n");
cat("Number of Predictor Variables: ", NumPredVariables, "\n", sep="");
cat(" Number of Response Variables: ", NumRespVariables, "\n", sep="");
cat("           Max ALR Var to Fit: ", NumMaxALRVariables, "\n", sep="");
cat("\n");
cat("List of ALR Categories to Include (instead of using Top):", ALRCategListFile, "\n", sep="");
cat("\n");
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
##############################################################################

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

mask_matrix=function(val_mat, mask_mat, mask_thres, mask_val){
	masked_matrix=val_mat;
	masked_matrix[mask_mat>mask_thres]=mask_val;
	return(masked_matrix);
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
	label_zeros=T, counts=F, value.cex=1, 
	plot_col_dendr=F,
	plot_row_dendr=F
){

	cat("Working on: ", title, "\n");

        num_row=nrow(mat);
        num_col=ncol(mat);

	if(num_row==0 || num_col==0){
		cat("Nothing to plot.\n");
		return();
	}

	any_nas=any(is.na(mat));

	if(num_row==1 || any_nas){
		plot_row_dendr=F;
	}
	if(num_col==1 || any_nas){
		plot_col_dendr=F;
	}

	row_names=rownames(mat);
	col_names=colnames(mat);

	orig.par=par(no.readonly=T);

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
                plot_min=min(mat, na.rm=T);
        }
        if(is.na(plot_max)){
                plot_max=max(mat, na.rm=T);
        }

	if(plot_min>=-1 && plot_max<=1){
		fractions_only=T;	
	}else{
		fractions_only=F;
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
		mat=mat[row_dendr[["names"]],, drop=F];
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

                        if(is.na(mat[y,x]) || mat[y,x]!=0 || label_zeros){
                                if(counts){
                                        text_lab=sprintf("%i", mat[y,x]);
                                }else{
                                        text_lab=sprintf(paste("%0.", deci_pts, "f", sep=""), mat[y,x]);
					if(fractions_only){
						if(!is.na(mat[y,x])){
							if(mat[y,x]==-1 || mat[y,x]==1){
								text_lab=as.integer(mat[y,x]);	
							}else{
								text_lab=gsub("0\\.","\\.", text_lab);
							}
						}
					}
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
		if(val == 0){ return("***");}
		if(val <= .001){ return("** ");}
		if(val <= .01){ return("*  ");}
		if(val <= .05){ return(".  ");}
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

##############################################################################
##############################################################################

# Open main output file
pdf(paste(OutputRoot, ".2smp_alr.pdf", sep=""), height=11, width=9.5);

NumMaxALRVariables=max(NumPredVariables, NumRespVariables);
input_files=load_and_reconcile_files(
		sumtab=list(fn=SummaryFile, shorten_cat_names_char=ShortenCategoryNames, 
			return_top=NumMaxALRVariables, specific_cat_fn=ALRCategListFile),
		factors=list(fn=FactorsFile),
		pairs=list(fn=PairingsFile, a_cname=ResponseName, b_cname=PredictorName),
		covariates=list(fn=ModelVarFile),
		grpvar=list(fn=""),
		reqvar=list(fn=RequiredFile)
	);

counts=input_files[["SummaryTable_counts"]];
factors=input_files[["Factors"]];
good_pairs_map=input_files[["PairsMap"]];
model_var_arr=input_files[["Covariates"]];
required_arr=input_files[["RequiredVariables"]];

##############################################################################
##############################################################################

# Normalize
cat("Normalizing counts...\n");
counts=counts+.5;
normalized=normalize(counts);

##############################################################################

resp_alr_struct=additive_log_rato(normalized);
alr_categories_val=resp_alr_struct$transformed;
alr_cat_names=colnames(alr_categories_val);
NumRespVariables=ncol(alr_categories_val);

cat("\n");
cat("ALR Category Summary:\n");
s=summary(alr_categories_val);
print(s);
plot_text(c(
	"ALR Categories Summary:",
	"\n",
	capture.output(print(s))
));
plot_histograms(alr_categories_val);

##############################################################################

# Order/Split the data....
# Order the pairings map by response sample IDs
response_sample_ids=good_pairs_map[,ResponseName];
predictor_sample_ids=good_pairs_map[,PredictorName];

# Extract the predictor ALR and factors values in the right order
response_alr=alr_categories_val[response_sample_ids,,drop=F];
predictor_alr=alr_categories_val[predictor_sample_ids,,drop=F];

##############################################################################

# Plot relationship between predictors and response from same taxa
alr_names=colnames(alr_categories_val);
plots_per_page=6
par(mfrow=c(plots_per_page,2));
par(oma=c(0,0,0,0));
par(mar=c(4,4,3,2));
median_delta=numeric(NumRespVariables);
names(median_delta)=alr_names;
for(cat_ix in 1:NumRespVariables){

	cur_alr_resp=response_alr[,cat_ix];
	cur_alr_pred=predictor_alr[,cat_ix];


	if((cat_ix %% plots_per_page)==0 || cat_ix==NumRespVariables){
		bottom_label=PredictorName;
	}else{
		bottom_label="";
	}

	# Plot matching pred/resp off reference line
	plot(cur_alr_pred, cur_alr_resp, type="n", main=alr_names[cat_ix],
		xlab=bottom_label, ylab=ResponseName);
	abline(0, 1, col="blue");
	points(cur_alr_pred, cur_alr_resp, main=alr_names[cat_ix]);

	deltas=cur_alr_resp-cur_alr_pred;
	median_delta[cat_ix]=median(deltas);
	hist(deltas, main=paste(alr_names[cat_ix], " ALR: ", ResponseName, " - ", PredictorName, sep=""), 
		xlab="");
	
}

# Generate bar plots for differences between response and predictors
par(mfrow=c(2,1));
par(mar=c(20,4,3,1));
color_arr=rainbow(NumRespVariables, start=0, end=4/6);
barplot(median_delta, horiz=F, las=2, main="", col=color_arr);
title(main=paste("Median Difference between '", ResponseName, "' and '", PredictorName, "'", sep=""), line=2);
title(main="(Ordered by decreasing abundance)", line=1, cex=.7);


dec_del_ix=order(median_delta, decreasing=T);
barplot(median_delta[dec_del_ix], horiz=F, las=2, main="", col=color_arr[dec_del_ix]);
title(main=paste("Median Difference between '", ResponseName, "' and '", PredictorName, "'", sep=""), line=2);
title(main="(Ordered by decreasing median difference)", line=1, cex=.7);

par(mfrow=c(1,1));
#print(response_alr);
#print(predictor_alr);
#print(factors);


# Store the results
num_model_pred=length(model_var_arr);

cov_model_string=paste(model_var_arr, collapse="+");
mmat=model.matrix(as.formula(paste("~", cov_model_string, "-1")), data=as.data.frame(factors));
cov_coeff_names=colnames(mmat);
num_cov_coeff_names=length(cov_coeff_names);

cat("Anticipated Coefficient Names:\n");
print(cov_coeff_names);
cat("\n");

category_alr_coef_mat=matrix(NA, nrow=NumRespVariables, ncol=NumRespVariables,
	dimnames=list(alr_cat_names, alr_cat_names));
category_alr_pval_mat=matrix(NA, nrow=NumRespVariables, ncol=NumRespVariables,
	dimnames=list(alr_cat_names, alr_cat_names));

covariates_coef_mat=matrix(NA, nrow=NumRespVariables, ncol=num_cov_coeff_names, 
	dimnames=list(alr_cat_names, cov_coeff_names));
covariates_pval_mat=matrix(NA, nrow=NumRespVariables, ncol=num_cov_coeff_names, 
	dimnames=list(alr_cat_names, cov_coeff_names));

rsqrd_mat=matrix(NA, nrow=NumRespVariables, ncol=4, 
	dimnames=list(alr_cat_names[1:NumRespVariables], c("R^2", "Adj-R^2", "Reduced Adj-R^2", "ALR Contrib")));

model_pval_mat=matrix(NA, nrow=NumRespVariables, ncol=1, 
	dimnames=list(alr_cat_names[1:NumRespVariables], c("F-statistic P-value")));

# Fit the regression model

cat("Num Response ALR Variables: ", NumRespVariables, "\n");
cat("Num Predictor ALR Variables: ", NumPredVariables, "\n");

all_pred_alr_names=colnames(predictor_alr);

if(NumPredVariables>0){
	top_alr_pred_names=all_pred_alr_names[1:NumPredVariables];
}else{
	top_alr_pred_names=c();
}
cat("Top ALR Predictors Available: \n");
print(top_alr_pred_names);

cat("*************************************************\n\n");

for(resp_ix in 1:NumRespVariables){
	alr_resp=response_alr[,resp_ix,drop=F];
	resp_cat_name=colnames(alr_resp);
	alr_resp=as.vector(alr_resp);
	
	cat("Fitting: ", resp_ix, ".) ", resp_cat_name, "\n");

	# Always include the response ALR category in the predictor

	if(any(resp_cat_name == top_alr_pred_names)){
		cat("Current response variable is in list of top predictors to include.\n");
		alr_predictors=top_alr_pred_names;
	}else{
		cat("Current response variable is NOT in list of top predictors to include,\n");
		cat("  so including Top N-1 categories as well.\n");
		if(NumPredVariables>1){
			alr_predictors=c(resp_cat_name, top_alr_pred_names[1:(NumPredVariables-1)]);
		}else{
			alr_predictors=resp_cat_name;
		}
	}
	cat("ALR Predictors to include: \n");
	print(alr_predictors);

	alr_pred=predictor_alr[,alr_predictors, drop=F];
	alr_pred_names=colnames(alr_pred);
	model_pred_df=as.data.frame(cbind(alr_pred, factors));
	num_samples_in_df=nrow(model_pred_df);
	cat("Num samples: ", num_samples_in_df, "\n", sep="");

	# Build formula string for full and reduced model

	full_variables=c(alr_pred_names,model_var_arr);
	num_full_variables=length(full_variables);

	model_str=paste("alr_resp ~ ", paste(full_variables, collapse=" + "), sep="");

	num_reduced_variables=length(model_var_arr);
	if(num_reduced_variables==0){
		model_reduced_str=paste("alr_resp ~ 1", sep="");
	}else{
		model_reduced_str=paste("alr_resp ~ ", paste(model_var_arr, collapse=" + "), sep="");
	}

	cat("Model String: \n");
	cat("Full (", num_full_variables, " variables): \n", sep="");
	print(model_str);
	cat("Reduced (", num_reduced_variables, " variables):\n", sep="");
	print(model_reduced_str);

	# Fit full and reduced model
	cat("Fitting Full...\n");
	lm_fit=lm(as.formula(model_str), data=model_pred_df);

	

	cat("Fitting Reduced...\n");
	lm_reduced_fit=lm(as.formula(model_reduced_str), data=model_pred_df);

	# Summarize full and reduced model
	sum_fit=summary(lm_fit);
	sum_reduced_fit=summary(lm_reduced_fit);

	# ANOVA on full model
	anova_res=anova(lm_fit);

	plot_text(c(
		paste(resp_ix, ".) ", resp_cat_name, ":", sep=""),
		"",
		capture.output(print(sum_fit))
		)
	);

	plot_text(c(
		paste(resp_ix, ".) ", resp_cat_name, ":", sep=""),
                "",
		capture.output(print(anova_res))
		)
	);

	tryCatch({
		mmps(lm_fit);
	}, error=function(e){
		print(e);
		plot_text(capture.output(print(e)));
	});

	model_coef_names=setdiff(rownames(sum_fit$coefficients), "(Intercept)");

	# Save the ALR result stats
	cat_names=intersect(model_coef_names, alr_pred_names);
	category_alr_coef_mat[resp_cat_name, cat_names]=sum_fit$coefficients[cat_names,"Estimate"];
	category_alr_pval_mat[resp_cat_name, cat_names]=sum_fit$coefficients[cat_names,"Pr(>|t|)"];

	# Save the covariate result stats
	cat_names=intersect(model_coef_names, cov_coeff_names);
	covariates_coef_mat[resp_cat_name, cat_names]=sum_fit$coefficients[cat_names,"Estimate"];
	covariates_pval_mat[resp_cat_name, cat_names]=sum_fit$coefficients[cat_names,"Pr(>|t|)"];

	rsqrd_mat[resp_cat_name, "R^2"]=sum_fit$r.squared;
	rsqrd_mat[resp_cat_name, "Adj-R^2"]=sum_fit$adj.r.squared;
	rsqrd_mat[resp_cat_name, "Reduced Adj-R^2"]=sum_reduced_fit$adj.r.squared;
	rsqrd_mat[resp_cat_name, "ALR Contrib"]=sum_fit$adj.r.squared-sum_reduced_fit$adj.r.squared;

	model_pval_mat[resp_cat_name, "F-statistic P-value"]=
		1-pf(sum_fit$fstatistic[1], sum_fit$fstatistic[2], sum_fit$fstatistic[3]);

	cat("\n\n*************************************************\n\n");

}

###############################################################################

all.nas=apply(covariates_coef_mat, 2, function(x){all(is.na(x))});
covariates_coef_mat=covariates_coef_mat[,!all.nas,drop=F];

all.nas=apply(covariates_pval_mat, 2, function(x){all(is.na(x))});
covariates_pval_mat=covariates_pval_mat[,!all.nas,drop=F];

print(rsqrd_mat);

par(oma=c(2,1,5,2));

# ALR Coefficients
paint_matrix(category_alr_coef_mat, 
	title=paste("Top ", NumPredVariables, " '", PredictorName,"' Predictor ALR Coefficients for Top ", 
		NumRespVariables, " '", ResponseName, "' Responses ALR", sep=""), 
	deci_pts=2, value.cex=.8);
mtext(PredictorName, side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# ALR Coefficients Clustered
paint_matrix(category_alr_coef_mat[,1:NumPredVariables, drop=F], 
	title=paste("Top ", NumPredVariables, " '", PredictorName,"' Predictor ALR Coefficients for Top ",
		NumRespVariables, " '", ResponseName, "' Responses ALR", sep=""), 
	 plot_col_dendr=T, plot_row_dendr=T,
	deci_pts=2, value.cex=.8);
mtext(PredictorName, side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# ALR P-Values
paint_matrix(category_alr_pval_mat, plot_min=0, plot_max=1,
	title=paste("Top ", NumPredVariables, " '", PredictorName,"' Predictor ALR P-Values for Top ", 
		NumRespVariables, " '", ResponseName, "' Responses ALR", sep=""), 
	high_is_hot=F, deci_pts=2, value.cex=.8);
mtext(PredictorName, side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# ALR Coefficients w/ P-value masked
limit=max(abs(range(category_alr_coef_mat)));
category_alr_coef_masked_mat=mask_matrix(
	val_mat=category_alr_coef_mat, 
	mask_mat=category_alr_pval_mat, 
	mask_thres=0.05, 
	mask_val=0.0);
paint_matrix(category_alr_coef_masked_mat,
	title=paste("Top ", NumPredVariables, " '", 
		PredictorName,"' Predictor ALR Coeff for Top ", NumRespVariables, " '", ResponseName,
		"' Responses ALR P-Val(<.05) Maskd", sep=""), 
	label_zeros=F, high_is_hot=F, deci_pts=2, value.cex=.8);
mtext(PredictorName, side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# ALR P-Values Clustered
paint_matrix(category_alr_pval_mat[,1:NumPredVariables, drop=F], plot_min=0, plot_max=1,
	title=paste("Top ", NumPredVariables, " '", PredictorName,"' Predictor ALR P-Values for Top ", 
		NumRespVariables, " '", ResponseName, "' Responses ALR", sep=""), 
	plot_col_dendr=T, plot_row_dendr=T,
	high_is_hot=F, deci_pts=2, value.cex=.8);
mtext(PredictorName, side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# Covariate Coefficients
paint_matrix(covariates_coef_mat, 
	title=paste("Predictor Coefficients for Top ", NumRespVariables, " '", ResponseName, "' Categories", sep=""),
	deci_pts=2, value.cex=.8);
mtext("Covariates & Predictors", side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# Covariate Coefficients w/ Clustering
paint_matrix(covariates_coef_mat, 
	title=paste("Covariates Coefficients for Top ", NumRespVariables, " '", ResponseName, "' Categories", sep=""),
	plot_col_dendr=T, plot_row_dendr=T,
	deci_pts=2, value.cex=.8);
mtext("Covariates & Predictors", side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# Covariate P-Values
paint_matrix(covariates_pval_mat, plot_min=0, plot_max=1, 
	title=paste("Covariates P-Values for Top ", NumRespVariables, " '", ResponseName, "' Categories", sep=""),
	high_is_hot=F, deci_pts=2, value.cex=.8);
mtext("Covariates & Predictors", side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# Covariate Coefficients w/ P-value masked
covariates_coef_masked_mat=mask_matrix(
        val_mat=covariates_coef_mat, mask_mat=covariates_pval_mat,
        mask_thres=0.1, mask_val=0.0);
paint_matrix(covariates_coef_masked_mat,,
        title=paste("Predictor Coeff for Top ", NumRespVariables, " '", ResponseName,
                "' Responses ALR P-Val(<.10) Maskd", sep=""),
        label_zeros=F, high_is_hot=F, deci_pts=2, value.cex=.8);
mtext("Covariates & Predictors", side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

covariates_coef_masked_mat=mask_matrix(
        val_mat=covariates_coef_mat, mask_mat=covariates_pval_mat,
        mask_thres=0.05, mask_val=0.0);
paint_matrix(covariates_coef_masked_mat,,
        title=paste("Predictor Coeff for Top ", NumRespVariables, " '", ResponseName,
                "' Responses ALR P-Val(<.05) Maskd", sep=""),
        label_zeros=F, high_is_hot=F, deci_pts=2, value.cex=.8);
mtext("Covariates & Predictors", side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

covariates_coef_masked_mat=mask_matrix(
        val_mat=covariates_coef_mat, mask_mat=covariates_pval_mat,
        mask_thres=0.01, mask_val=0.0);
paint_matrix(covariates_coef_masked_mat,,
        title=paste("Predictor Coeff for Top ", NumRespVariables, " '", ResponseName,
                "' Responses ALR P-Val(<.01) Maskd", sep=""),
        label_zeros=F, high_is_hot=F, deci_pts=2, value.cex=.8);
mtext("Covariates & Predictors", side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# Covariate P-Values Clustered
paint_matrix(covariates_pval_mat, plot_min=0, plot_max=1, 
	title=paste("Covariates P-Values for Top ", NumRespVariables, " '", ResponseName, "' Categories", sep=""),
	plot_col_dendr=T, plot_row_dendr=T,
	high_is_hot=F, deci_pts=2, value.cex=.8);
mtext("Covariates & Predictors", side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# R^2
paint_matrix(rsqrd_mat, plot_min=0, plot_max=1, 
	title=paste("Explained Variation for Top ", NumRespVariables, " Responses: R^2", sep=""));
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# Model pval
paint_matrix(model_pval_mat, plot_min=0, plot_max=1, high_is_hot=F, 
	title=paste("Model F-Statistic P-Values for Top ", NumRespVariables, sep=""));
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# Plot Pred->Resp for ALR
plot_pred_resp_bar=function(coef_mat, pval_mat, title=""){
	val=diag(coef_mat);
	pval=diag(pval_mat);
	cat_names=rownames(coef_mat);

	ord=order(val, decreasing=T);

	val=val[ord];
	pval=pval[ord];
	cat_names=cat_names[ord];

        num_colors=50;
        color_arr=rainbow(num_colors, start=0, end=4/6);

        # Provide a means to map values to an (color) index
        remap=function(in_val, in_range, out_range){
                in_prop=(in_val-in_range[1])/(in_range[2]-in_range[1])
                out_val=in_prop*(out_range[2]-out_range[1])+out_range[1];
                return(out_val);
        }

	colors=ceiling(remap(pval, c(0,1), c(1,num_colors)));

	barplot(val, names.arg=cat_names, las=2, col=color_arr[colors],
		main=title
		);
}

par(oma=c(20, 5, 5, 1));
plot_pred_resp_bar(category_alr_coef_mat, category_alr_pval_mat);
mtext(paste("Predictability of '", ResponseName, "' ALR based on '", PredictorName, "' ALR", sep=""),
	side=3, font=2, line=2
	);
mtext("After Controlling for Covariates & Predictors", side=3, font=2, line=1);
mtext("Regression Coefficient", side=2, font=1, line=3.75);

###############################################################################

write.table(category_alr_pval_mat[,1:NumPredVariables, drop=F],
	file=paste(OutputRoot, ".alr_as_pred.nolabels.pvals.tsv", sep=""), 
	sep="\t", quote=F, col.names=NA, row.names=T);

write.table(category_alr_coef_mat[,1:NumPredVariables, drop=F],
	file=paste(OutputRoot, ".alr_as_pred.nolabels.coefs.tsv", sep=""), 
	sep="\t", quote=F, col.names=NA, row.names=T);

# Write out predictor ALR predictors
rownames(category_alr_pval_mat)=paste(ResponseName,".",rownames(category_alr_pval_mat), sep="");
colnames(category_alr_pval_mat)=paste(PredictorName,".",colnames(category_alr_pval_mat), sep="");

rownames(category_alr_coef_mat)=paste(ResponseName,".",rownames(category_alr_coef_mat), sep="");
colnames(category_alr_coef_mat)=paste(PredictorName,".",colnames(category_alr_coef_mat), sep="");

write.table(category_alr_pval_mat[,1:NumPredVariables, drop=F],
	file=paste(OutputRoot, ".alr_as_pred.pvals.tsv", sep=""), 
	sep="\t", quote=F, col.names=NA, row.names=T);

write.table(category_alr_coef_mat[,1:NumPredVariables, drop=F],
	file=paste(OutputRoot, ".alr_as_pred.coefs.tsv", sep=""), 
	sep="\t", quote=F, col.names=NA, row.names=T);

write.table(t(category_alr_pval_mat[,1:NumPredVariables, drop=F]),
	file=paste(OutputRoot, ".alr_as_pred.tp.pvals.tsv", sep=""), 
	sep="\t", quote=F, col.names=NA, row.names=T);

write.table(t(category_alr_coef_mat[,1:NumPredVariables, drop=F]),
	file=paste(OutputRoot, ".alr_as_pred.tp.coefs.tsv", sep=""), 
	sep="\t", quote=F, col.names=NA, row.names=T);


###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
