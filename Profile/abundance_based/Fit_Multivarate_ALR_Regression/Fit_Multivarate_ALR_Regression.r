#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
library(car);
options(useFancyQuotes=F);

params=c(
	"summary_file", "s", 1, "character",
	"factors", "f", 1, "character",
	"num_variables", "p", 2, "numeric",
	"reference_levels", "r", 2, "character",
	"outputroot", "o", 2, "character",
	"model", "m", 2, "character",
	"contains_remaining", "R", 2, "logical",
	"shorten_category_names", "x", 2, "logical",
	"test_run", "t", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table>\n",
	"	-f <factors>\n",
	"	[-r <reference levels file>]\n",
	"	[-p <number of variables>]\n",
	"	[-o <output filename root>]\n",
	"	[-m <model formula string>]\n",
	"	[-R (pay attention to 'remaining' category)]\n",
	"\n",
	"	[-x (shorten category names)]\n",
	"	[-t (test run flag)]\n",
	"\n",
	"This script will read in the summary file table, then perform\n",
	"a multivariate logistic regression on the the top categories\n",
	"using the factors/predictors in the factor file.\n",
	"\n",
	"For testing purposes, you can append the factor name with 'IGNORE.' and\n",
	"the factor will be ignored.\n",
	"\n",
	"The -m option will allow you to specify a more sophisticated regression\n",
	"model.\n",
	"	For example: \n",
	"		-m \"F1 + F2 + F1*F2\"\n",
	"\n",
	"	Will Fit: \n",
	"		ALR = b0 + b1*F1 + b2*F2 + b3*F1*F2\n",
	"\n",
	"Without the -m option, all factors in the factors file will be fit as main effects.\n",
	"\n",
	"If the -R flag is set, a 'remaining' category will be be included in the denominator\n",
	"	independent of how large it is.  I.e., do not use it as one of the response variables.\n",
	"\n");

if(!length(opt$summary_file) || !length(opt$factors)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputRoot=gsub(".summary_table.xls", "", opt$summary_file);
}else{
	OutputRoot=opt$outputroot;
}

if(!length(opt$num_variables)){
	NumVariables=20;
}else{
	NumVariables=opt$num_variables;
}

if(!length(opt$reference_levels)){
        ReferenceLevelsFile="";
}else{
        ReferenceLevelsFile=opt$reference_levels;
}

if(length(opt$model)){
	Model=opt$model;
}else{
	Model="All Factors";
}

if(length(opt$contains_remaining)){
	UseRemaining=T;
}else{
	UseRemaining=F;
}

if(length(opt$shorten_category_names)){
	ShortenCategoryNames=T;
}else{
	ShortenCategoryNames=F;
}

if(length(opt$test_run)){
	TestRun=T;
}else{
	TestRun=F;
}

SummaryFile=opt$summary_file;
FactorsFile=opt$factors;

cat("\n");
cat("Summary File : ", SummaryFile, "\n", sep="");
cat("Factors File: ", FactorsFile, "\n", sep="");
cat("Number of Response Variables: ", NumVariables, "\n", sep="");
cat("Reference Levels File: ", ReferenceLevelsFile, "\n", sep="");
cat("Output File: ", OutputRoot, "\n", sep="");
cat("Model: ", Model, "\n", sep="");
cat("Use Remaining? ", UseRemaining, "\n");
cat("\n");

if(TestRun){
	cat("***********************************************\n");
	cat("*  Test Run                                   *\n");
	cat("***********************************************\n");
	rnd=paste(".", sprintf("%i",sample(1000, 1)), sep="");
}else{
	rnd="";
}

options(width=120);
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
        for(i in 1:num_factors_to_relevel){
                tmp=factors[,relevel_names[i]];
                #print(tmp);
                tmp=relevel(tmp, ref_lev_mat[i, 1]);
                #print(tmp);
                factors[,relevel_names[i]]=tmp;
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

plot_correl_heatmap=function(mat, title="", noPrintZeros=F, guideLines=F){

	if(is.null(dim(mat))){
		cat(title, " Matrix is NULL.  No heatmap generated.\n");
		return;
	}

	cat("Plotting: ", title, "\n");

	par(family="Courier");
	par(oma=c(10, 10, 1, .5));
	par(mar=c(6.1, 5.1, .5, .5));

	# Generate colors from red to blue
        colors=(rainbow(2^16, start=0, end=0.65));

	# Remember that rows and columsn are reversed in the image
	image(1:nrow(mat),1:ncol(mat), mat,
                xaxt="n", yaxt="n",
                xlab="", ylab="",
                col=colors
        );

	# Pad strings
	cnames=paste(colnames(mat), " ", sep="");
	rnames=paste(rownames(mat), " ", sep="");
	
	# Get longest length of each column or row label
	cname_max_len=max(nchar(cnames));
	rname_max_len=max(nchar(rnames));

	# Get the number of rows and columns
	ncols=ncol(mat);
	nrows=nrow(mat);

	cscale=min(c(45/cname_max_len, 55/ncols));
	rscale=min(c(45/rname_max_len, 55/nrows));

	cscale=min(3, cscale);
	rscale=min(3, rscale);

        max_width=max(nchar(sprintf("%.2f",mat)));
	#cell_cex=sqrt(min(c(cscale, rscale))^2);
	cell_cex=(1/max_width)*sqrt(min(c(cscale, rscale))^2);

        for(i in 1:nrow(mat)){
                for(j in 1:ncol(mat)){

			if(!is.na(mat[i,j]) && (noPrintZeros && mat[i,j]==0)){
				# Skip	
			}else{
				str=sprintf("%.2f",mat[i,j]);
				str=gsub("^0\\.",".", str);
				text(i,j,labels=str, cex=cell_cex, srt=45);
			}
                }
        }

	# Plot guidelines
	if(guideLines){
		
		splits=c(2,3,4,5);
		h_remainder=ncols %% splits;
		best_h_split=splits[max(which(h_remainder==min(h_remainder)))];
		if(ncols>best_h_split){
			h_line_pos=seq(best_h_split, ncols, best_h_split)+.5;
			abline(h=h_line_pos, col="black", lty="dashed");
			abline(h=h_line_pos, col="white", lty="dotted");
		}

		v_remainder=nrows %% splits;
		best_v_split=splits[max(which(v_remainder==min(v_remainder)))];
		if(is.na(best_v_split)){
			best_v_split=max(splits);
		}
		if(nrows>best_v_split){
			v_line_pos=seq(best_v_split, nrows, best_v_split)+.5;
			abline(v=v_line_pos, col="black", lty="dashed");
			abline(v=v_line_pos, col="white", lty="dotted");
		}

	}

        # Plot the labels
        mtext(cnames, at=1:ncols, side=2, las=2, cex=cscale);
        mtext(rnames, at=1:nrows, side=1, las=2, cex=rscale);

	# Plot the title
	mtext(title, line=0, at=nrows*.5, side=3, font=2);

}

write_top_categorical_effects_by_factor=function(output_fn, coeff_mat, pval_mat, top_n=10){

	cat("Writing top category effects by factor...\n");

	pval_cat=colnames(pval_mat);
	pval_fac=rownames(pval_mat);
	
	coeff_cat=colnames(coeff_mat);
	coeff_fac=rownames(coeff_mat);

	if(!all(pval_fac==coeff_fac) && !all(pval_cat==coeff_cat)){
		cat("Error: categories or factors don't match up.\n");
		quit(status=-1);
	}else{
		factors=pval_fac;
		categories=pval_cat;
	}

	
	mat_buf=matrix(0, ncol=2, nrow=length(categories));
	colnames(mat_buf)=c("ALR", "p-value");
	rownames(mat_buf)=categories;


	sig_fun_str=function(x){
		if(x <= 0.0001){
			return("****");
		}else if(x <= 0.001){
			return("***");
		}else if(x <= 0.01){
			return("**");
		}else if(x <= 0.05){
			return("*");
		}else if(x <= 0.1){
			return(".")
		}else{
			return("");
		}
	}

	fh=file(output_fn, "w");
	

	for(cur_factor in factors){
		cat("Working on: ", cur_factor, "\n");

		mat_buf[categories,"ALR"]=coeff_mat[cur_factor,categories];	
		mat_buf[categories,"p-value"]=pval_mat[cur_factor,categories];	
		
		# Sort
		sort_ix=order(mat_buf[,"ALR"], decreasing=T);
		mat_buf=mat_buf[sort_ix,];
		signif=sapply(mat_buf[,2], sig_fun_str);
		sort_cat=rownames(mat_buf);

		cat(file=fh, cur_factor, ":,Category,ALR,p-value,Signif\n");

		# Output Top N
		ix=1;
		while(mat_buf[ix,"ALR"]>0 && ix<=top_n){
			vals=c(paste(ix,"+", sep=""), sort_cat[ix], mat_buf[ix,"ALR"], mat_buf[ix, "p-value"], signif[ix]);
			cat(file=fh, paste(vals, collapse=","), "\n");	
			ix=ix+1;
		}

		# Separator
		cat(file=fh, "...\n");

		# Output Bottom N
		num_cats=nrow(mat_buf);
		ix=0;
		while(mat_buf[num_cats-ix,"ALR"]<0 && (ix<top_n)){
			vals=c(paste(ix+1, "-", sep=""), sort_cat[num_cats-ix], mat_buf[num_cats-ix,"ALR"], 
				mat_buf[num_cats-ix, "p-value"], signif[num_cats-ix]);
			cat(file=fh, paste(vals, collapse=","), "\n");	
			ix=ix+1;
		}

		cat(file=fh, "\n\n");
	}


	close(fh);

}


##############################################################################

# Load summary file table counts 
counts=load_summary_file(SummaryFile);
num_taxa=ncol(counts);
num_samples=nrow(counts);
#print(counts);

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
		normalized_remaining_col_dat=normalized[,remaining_ix, drop=F];
		normalized=normalized[,-remaining_ix];
	}
}

# Shorten cateogry names
if(ShortenCategoryNames){
	full_names=colnames(counts);
	splits=strsplit(full_names, " ");
	short_names=character();
	for(i in 1:length(full_names)){
		short_names[i]=tail(splits[[i]], 1);
	}
	colnames(counts)=short_names;
	cat("Names have been shortened.\n");
}


# Reorder by abundance
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
factors=load_factors(FactorsFile);
factor_names=colnames(factors);
num_factors=ncol(factors);
factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);

cat("\n");
cat(num_factors, " Factor(s) Loaded:\n", sep="");
print(factor_names);
cat("\n");

# Remove factors not in model
if(Model!="All Factors"){
	model_pred_vars=gsub("\\+", " ", Model);
	model_pred_vars=gsub("\\*", " ", model_pred_vars);
	model_pred_vars=gsub("\\:", " ", model_pred_vars);
	model_pred_vars=unique(strsplit(model_pred_vars, " ")[[1]]);
	used_factors=intersect(model_pred_vars, factor_names);

	factors=factors[, used_factors, drop=F];
	factor_names=colnames(factors);
	num_factors=ncol(factors);
}

# Relevel factor levels
if(ReferenceLevelsFile!=""){
        ref_lev_mat=load_reference_levels_file(ReferenceLevelsFile)
        factors=relevel_factors(factors, ref_lev_mat);
}else{
        cat("No Reference Levels File specified.\n");
}

continuous_factors=factors;
is_continous_factor=logical(num_factors);

for(f in 1:num_factors){
	level_info=levels(factors[,f]);
	is_continous_factor[f]=is.null(level_info);

	if(is_continous_factor[f]){
		# do nothing
	}else if(length(level_info)==2){
		# Convert two level factors to numeric
		is_continous_factor[f]=TRUE;
		continuous_factors[,f]=as.integer(continuous_factors[,f]);
	}else{
		is_continous_factor[f]=FALSE;
	}
}

continuous_factors=continuous_factors[,is_continous_factor, drop=F];
print(continuous_factors);
factor_correlations=cor(continuous_factors);

#pdf(paste(OutputRoot, ".factor_cor.pdf", sep=""), height=11, width=8.5);
#plot_correl_heatmap(factor_correlations, title="Factor Correlations");
#dev.off();

##############################################################################

# Reconcile factors with samples
factor_sample_ids=rownames(factors);
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

cat("Samples missing from count information:\n");
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
factors=factors[shared_sample_ids,,drop=F];


##############################################################################

extract_top_categories=function(ordered_normalized, top){

	num_samples=nrow(ordered_normalized);
	num_categories=ncol(ordered_normalized);

	cat("Samples: ", num_samples, "\n");
	cat("Categories: ", num_categories, "\n");
	
	num_saved=min(c(num_categories, top+1));

	cat("Top Requested to Extract: ", top, "\n");
	cat("Columns to Extract: ", num_saved, "\n");

	top_cat=matrix(0, nrow=num_samples, ncol=num_saved);
	top=num_saved-1;

	# Extract top categories requested
	top_cat[,1:top]=ordered_normalized[,1:top];

	# Included remaineder as sum of remaining categories
	top_cat[,(top+1)]=apply(
		ordered_normalized[,(top+1):num_categories, drop=F],
		1, sum);

	rownames(top_cat)=rownames(ordered_normalized);
	colnames(top_cat)=c(colnames(ordered_normalized)[1:top], "Remaining");

	return(top_cat);
			
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

##############################################################################

# Assign 0's to values smaller than smallest abundance across entire dataset
min_assay=min(normalized[normalized!=0]);
cat("Lowest non-zero value: ", min_assay, "\n\n", sep="");
zero_replacment=min_assay/10;
normalized[normalized==0]=zero_replacment;

##############################################################################

if(num_top_taxa>= num_taxa){
	num_top_taxa = (num_taxa-1);
	cat("Number of taxa to work on was changed to: ", num_top_taxa, "\n");
}

# Perform ALR transform
responses=extract_top_categories(normalized, num_top_taxa);
resp_alr=additive_log_rato(responses)$transformed;

##############################################################################

plot_text=function(strings){
	par(family="Courier");
	par(oma=rep(.5,4));
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

##############################################################################
pdf(paste(OutputRoot, rnd, ".mlr.mmp.pdf", sep=""), height=11, width=8.5);

# Output the factor correlations
if(ncol(factor_correlations)>0){
	plot_correl_heatmap(factor_correlations, title="Factor Correlations");
}else{
	plot_text("Number of ordinal factors in model is zero, so no factor correlation heatmap was generated.");
}

##############################################################################
# Output reference factor levels

text=character();
text[1]="Reference factor levels:";
text[2]="";

for(i in 1:num_factors){
        fact_levels=levels(factors[,i]);
        if(!is.null(fact_levels)){
                fact_info=paste(factor_names[i], ": ", fact_levels[1], sep="");
        }else{
                fact_info=paste(factor_names[i], ": None (ordered factor)", sep="");
        }
        text=c(text, fact_info);

}

text=c(text, "");
text=c(text, paste("Number of Samples: ", num_samples, sep=""));
text=c(text, "");
text=c(text, "Description of Factor Levels and Samples:");

width_orig=options()$width;
options(width=80);
text=c(text, capture.output(summary(factors)));
options(width=width_orig);

plot_text(text);

##############################################################################

cat("Performing regression.\n");
cat("Extracting: ", num_top_taxa, " + 1 (remaining) categories.\n", sep="");

responses=extract_top_categories(normalized, num_top_taxa);
resp_alr_struct=additive_log_rato(responses);
transformed=resp_alr_struct$transformed;


if(Model!="All Factors"){
	model_pred_str=Model;
}else{
	model_pred_str=paste(factor_names, collapse=" + ");
}


# Try to perform MANOVA
model_string= paste("transformed ~", model_pred_str);
cat("\nFitting this multivariate model: ", model_string, "\n");
text=character();

mv_fit=tryCatch({
	mv_fit=lm(as.formula(model_string), data=factors);
}, error = function(e){
	print(e);
	text[1]="Could not perform multivariate anova on your data because the formula";
	text[2]="did not appear to be a fixed effect only model.";
	plot_text(text);
	cat("(There will be no mv_fit data structure to compute on.)\n");
});

manova_trial=tryCatch({
	manova_res=anova(mv_fit);
	res=list();
	res[["manova"]]=manova_res;
	res[["error"]]=NULL;
	res;
}, error = function(e){
	res=list();
	res[["manova"]]=NULL;
	res[["error"]]=e;
	res;
});

manova_success=ifelse(is.null(manova_trial[["manova"]]), F, T);

if(manova_success){
	manova_res=manova_trial[["manova"]];
	print(manova_res);
	cat("\n");
	manova_txt=capture.output(manova_res);
}else{
	manova_txt=paste("Error performing MANOVA: ", manova_trial[["error"]], sep="");
}

text[1]=paste("Multivariate Regression with ", num_top_taxa, " top taxa", sep="");
text[2]=paste("Proportion of overall mean abundance represented: ", prop_abundance_represented, sep="");
text[3]="";
text=c(text, strsplit(model_string, "(?<=.{80})", perl=T)[[1]]);
text=c(text, "");
text=c(text, manova_txt);
plot_text(text);

###############################################################################
# Compute taxonomic correlations and pvalues

# Compute and Plot Taxonomic correlations
cor_mat=cor(transformed);
plot_correl_heatmap(cor_mat, title="Category Correlations");

# Compute pvalues for correlations, Null Hypothesis is cor=0
cat("Computing P-values for Category Correlations...\n");
pval_matrix=matrix(NA, ncol=num_top_taxa, nrow=num_top_taxa);
colnames(pval_matrix)=colnames(transformed);
rownames(pval_matrix)=colnames(transformed);

pval_vect=numeric(num_top_taxa*(num_top_taxa-1)/2);

num_corr_to_test=0;
for(i in 2:num_top_taxa){
	for(j in 1:(i-1)){
		pval=cor.test(transformed[,i],transformed[,j])$p.value;
		pval_matrix[i,j]=pval;
		pval_matrix[j,i]=pval;
		num_corr_to_test=num_corr_to_test+1;	
		pval_vect[num_corr_to_test]=pval;
		
	}
}

# Plot pvalues and log10(pvalues);
plot_correl_heatmap(pval_matrix, title="Unadjusted Correlation P-values");
plot_correl_heatmap(log10(pval_matrix), title="Unadjusted Correlation Log10(P-values)");

# FDR adjust pvalues
cat("Adjusting P-values for Multiple Testing using Holm.\n");
adjust_pval_vect=p.adjust(pval_vect, method="holm");
fdr_pval_matrix=matrix(NA, ncol=num_top_taxa, nrow=num_top_taxa);
colnames(fdr_pval_matrix)=colnames(transformed);
rownames(fdr_pval_matrix)=colnames(transformed);

num_corr_to_test=0;
for(i in 2:num_top_taxa){
	for(j in 1:(i-1)){
		num_corr_to_test=num_corr_to_test+1;	
		fdr_pval_matrix[i,j]=adjust_pval_vect[num_corr_to_test];
		fdr_pval_matrix[j,i]=fdr_pval_matrix[i,j];
	}
}

# Plot Adjust p-values
plot_correl_heatmap(fdr_pval_matrix, title="Holm Adjusted Correlation P-values");
plot_correl_heatmap((fdr_pval_matrix<0.05)*cor_mat, title="Significant (<0.05) Correlation Holm Adjusted P-values",
	noPrintZeros=T, guideLines=T);

##############################################################################

uv_fit=list();

#model_string=paste("transformed[,1] ~", paste(factor_names, collapse=" + "));
#model_matrix=model.matrix(as.formula(model_string), data=factors);

#model_variables=attributes(model_matrix)$dimnames[[2]];
#num_model_variables=length(model_variables);

#print(model_variables);
#print(num_model_variables);

# Determine many many ANOVA coefficients were analyzed, and store p-values in matrix
model_matrix=(model.matrix(as.formula(model_string), data=factors));
num_coeff=ncol(model_matrix);
cat("Number of Coefficients Expected: ", num_coeff, "\n");
coeff_names=colnames(model_matrix);

uv_pval_mat=matrix(NA, nrow=num_coeff, ncol=num_top_taxa,
		dimnames=list(coeff_names, sorted_taxa_names[1:num_top_taxa]));

uv_coeff_mat=matrix(NA, nrow=num_coeff, ncol=num_top_taxa,
		dimnames=list(coeff_names, sorted_taxa_names[1:num_top_taxa]));


tmp_model_string= paste("transformed[,1] ~", model_pred_str);
test_uv_fit=lm(as.formula(tmp_model_string), data=factors);
anova_factor_names=setdiff(rownames(anova(test_uv_fit)), c("Residuals", "(Intercept)"));
print(anova_factor_names);

uv_anova_pval_mat=matrix(NA, nrow=length(anova_factor_names), ncol=num_top_taxa,
		dimnames=list(anova_factor_names, sorted_taxa_names[1:num_top_taxa]));


# Store R^2 for each taxa
rsqrd=numeric(num_top_taxa);
adj_rsqrd=numeric(num_top_taxa);

for(var_ix in 1:num_top_taxa){
	summary_txt=character();

	cat("\n");
	cat("##########################################################################\n");
	cat("#                                                                        #\n");
	cat("# ", sorted_taxa_names[var_ix], "\n");
	cat("#                                                                        #\n");
	cat("##########################################################################\n");

	ALR_Abundance=transformed[,var_ix];
	model_string= paste("ALR_Abundance ~", model_pred_str);
	cat("Fitting: ", model_string, "\n");

	# Compute Univariate fit
	uv_fit[[var_ix]]=lm(as.formula(model_string), data=factors);
	#print(uv_fit[[var_ix]]);

	# Analyze fit
	uv_summ=summary(uv_fit[[var_ix]]);

	# Save overall fit
	rsqrd[var_ix]=uv_summ$r.squared;
	adj_rsqrd[var_ix]=uv_summ$adj.r.squared;
	
	# Perform univariate ANOVA
	uv_anova=anova(uv_fit[[var_ix]]);
	avail_coef_names=setdiff(rownames(uv_anova), "Residuals");
	uv_anova_pval_mat[avail_coef_names, var_ix]=uv_anova[avail_coef_names, "Pr(>F)"];

	# Identify which coefficients could be used
	regress_table=uv_summ$coefficients;
	avail_coef_names=rownames(regress_table);
	uv_coeff_mat[avail_coef_names ,var_ix]=regress_table[avail_coef_names, "Estimate"];
	uv_pval_mat[avail_coef_names ,var_ix]=regress_table[avail_coef_names, "Pr(>|t|)"];

}

cat("\nUnivariate Regression Coefficients:\n");
print(uv_coeff_mat)
cat("\nUnivariate Regression P-values:\n");
print(uv_pval_mat)
cat("\nUnivariate ANOVA P-values:\n");
print(uv_anova_pval_mat);
cat("\n");

# Plot univariate ANOVA F tests
anova_factors=rownames(uv_anova_pval_mat);
plot_correl_heatmap(uv_anova_pval_mat, title="Univariate F-Tests Pr(>F)", guideLines=T);

# Remove NAs
not_estimable=apply(uv_coeff_mat, 1, function(x){ any(is.na(x))});
coef_names_not_estimable=rownames(uv_coeff_mat)[not_estimable];
uv_coeff_mat=uv_coeff_mat[!not_estimable,];
uv_pval_mat=uv_pval_mat[!not_estimable,];

# remove intercept
uv_coeff_mat=uv_coeff_mat[-1,]; 
uv_pval_mat=uv_pval_mat[-1,] 

# Plot pvalues
plot_correl_heatmap(uv_pval_mat, title="Univariate Coefficients Pr(>|t|)");

# Plot log(pvalues, 10)
log_uv_pval_mat=log(uv_pval_mat,10);
plot_correl_heatmap(log_uv_pval_mat, title="Univariate Coefficients Log10[Pr(>|t|)]");

# Plot Heatmap
if(ncol(log_uv_pval_mat)>=2 && nrow(log_uv_pval_mat)>=2){

	# Get current graphic settings so we can restore them.
	par_oma_before=par()$oma;
	par_mar_before=par()$mar;

	par(oma=c(3, 0, 0, 10));
	par(mar=c(0, 0, 0, 0));
	num_hm_colors=20;
	cl_hm_colors=(rainbow(num_hm_colors, start=0, end=0.65));
	heatmap(log_uv_pval_mat, col=cl_hm_colors, margins=c(27, 7));

	# Plot Legend
	par(mar=c(10,1,10,1), oma=c(0,0,0,0),  mfrow=c(2,1));
	log_uv_pval_range=range(log_uv_pval_mat);
	color_spans=seq(log_uv_pval_range[1], log_uv_pval_range[2], length.out=num_hm_colors);

	barplot(rep(1, num_hm_colors), col=cl_hm_colors, yaxt="n", 
		space=0, names.arg=sprintf("%3.4f", color_spans), las=2, main="HeatMap Legend (Log10[p-values])");
	barplot(rep(1, num_hm_colors), col=cl_hm_colors, yaxt="n", 
		space=0, names.arg=sprintf("%3.4f", 10^color_spans), las=2, main="HeatMap Legend (p-values)");

	# Restore prior graphics settings
	par(oma=par_oma_before, mar=par_mar_before, mfrow=c(1,1));

}else{
	cat("No heatmap generated because p-value matrix is not multi-dimensional.\n");
}

# Plot log pvalues sorted by most signficiant predictor and taxa
pred_ix=order(apply(log_uv_pval_mat, 1, mean));
taxa_ix=order(apply(log_uv_pval_mat, 2, mean));
plot_correl_heatmap(log_uv_pval_mat[pred_ix, taxa_ix, drop=F], title="Sorted Univariate Coefficients Log10[Pr(>|t|)]", guideLines=T);

# Plot R^2
rsqrd_mat=rbind(rsqrd, adj_rsqrd);
rownames(rsqrd_mat)=c("R^2", "Adjusted R^2");
colnames(rsqrd_mat)=sorted_taxa_names[1:num_top_taxa];
plot_correl_heatmap(rsqrd_mat, title="Univariate R Squared");

# Plot univariate coefficients
not_na_coeff=apply(uv_coeff_mat, 1, function(x){!any(is.na(x))});
plot_correl_heatmap(uv_coeff_mat[not_na_coeff,], title="Univariate Coefficients", guideLines=T);

# Significant Univariate Coefficients
uv_pval_vect=as.vector(uv_pval_mat);
adj_uv_pval_vect=p.adjust(uv_pval_vect, method="fdr");
adj_uv_pval_mat=matrix(adj_uv_pval_vect, ncol=ncol(uv_pval_mat));
sig_coeff=(adj_uv_pval_mat<0.05)*uv_coeff_mat[not_na_coeff,, drop=F];
plot_correl_heatmap(sig_coeff, "Significant (FDR) Coefficients", noPrintZeros=T, guideLines=T);

if(length(coef_names_not_estimable)){
	# Output not estimable coefficients
	out_lines=character();
	out_lines[1]="Coefficients not calculatable:";
	out_lines[2]="";
	out_lines=c(out_lines, coef_names_not_estimable);
	plot_text(out_lines);
}


# Write Top categories that have changed to file
write_top_categorical_effects_by_factor(paste(OutputRoot,".top_effects.csv", sep=""), uv_coeff_mat, uv_pval_mat, top_n=20);

##############################################################################

reg_coef_power=function(uv_reg_fit, factor=10, alpha=.05, power=.8){
	
	cat("---------------------------------------------\n");
	cat("Factor (effect size): ", factor, "\n");
	cat("alpha = ", alpha, " / power = ", power, "\n");
	summ_fit=summary(uv_reg_fit);
	coef_names=setdiff(rownames(summ_fit$coefficients), "(Intercept)");
	num_coeff=length(coef_names);
	model_matrix=uv_reg_fit$model;
	n=nrow(model_matrix);

	sigma_reg=summ_fit$sigma;	# Stderr of Residuals

	SSx=numeric(num_coeff);
	sigma_x=numeric(num_coeff);
	sigma_b=numeric(num_coeff);
	needed_sample_size=numeric(num_coeff);
	
	for(i in 1:num_coeff){
		cur_coef=coef_names[i];

		#cat("Current Coefficient: ", cur_coef, "\n");

		components=strsplit(cur_coef, ":")[[1]];
		x_values=apply(model_matrix[,components, drop=F], 1, prod);
		range_x=range(x_values);
		mean_x=mean(x_values);
		SSx[i]=sum((x_values-mean_x)^2);

		sigma_x[i]=sqrt(SSx[i]/n);
		sigma_b[i]=sigma_reg/sqrt(SSx[i]);

		b=log(factor)*diff(range_x);
		
		# Iterate over N, since the critical t value depends on N
		curN=0;
		N=200;
		iterations=0;
		while(curN!=N && iterations<10){
			curN=N;
			df=(curN-(num_coeff+1));
			df=ifelse(df>0, df, 1);
			t_alpha_2=qt(1-(alpha/2), df);
			#t_beta=(abs(b)/sigma_b[i]) - t_alpha_2;
			t_beta=qt(power, df);

			N=ceiling(((t_alpha_2 + t_beta)^2) * ((sigma_reg/(b*sigma_x[i]))^2));	

			if(0){
				cat("df: ", df, "\n");
				cat("t_alpha_2: ", t_alpha_2, "\n");
				cat("t_beta: ", t_beta, "\n");
				cat("sigma_reg: ", sigma_reg, "\n");
				cat("sigma_x: ", sigma_x[i], "\n");
				cat("b: ", b, "\n");
				cat("N[", iterations, "]:", N, "\n\n");
			}

			iterations=iterations+1;
		}


		needed_sample_size[i]=N;
		cat("\tsigma_x: ", sigma_x[i], "\n");
		cat("\tsigma_b: ", sigma_b[i], "\n");
		cat("\tEffect: ", b, "\n");
		cat("\t", cur_coef, " N: ", needed_sample_size[i], "\n");
	}
	
	cat("---------------------------------------------\n");

}

##############################################################################
# Plot univariate analyses

par(oma=c(0,0,3,0));
for(var_ix in 1:num_top_taxa){

	setHook("plot.new", function(){mtext(sorted_taxa_names[var_ix], outer=T, line=-.5);});

	# Output univariate ANOVA results
	summary_txt=c();
	summary_txt[1]="Univariate Regression:";
	summary_txt[2]="";
	summary_txt[3]=paste(var_ix, ".) ", sorted_taxa_names[var_ix], sep="");
	summary_txt[4]="";
	summary_txt[5]=paste("Mean abundance: ",  sprintf("%3.1f%%",mean_abund[var_ix]*100), sep="");
	summary_txt[6]="";
	summary_txt[7]=paste("R^2: ", sprintf("%3.4f", rsqrd_mat[1,var_ix]), sep="");
	summary_txt[8]=paste("Adjusted R^2: ", sprintf("%3.4f", rsqrd_mat[2, var_ix]), sep="");
	summary_txt[9]="";
	summary_txt=c(summary_txt, capture.output(anova(uv_fit[[var_ix]])));
	plot_text(summary_txt);	

	
	# Regenerate summary table after removing NAs
	uv_summary=summary(uv_fit[[var_ix]]);
	tmp_dtfr=as.data.frame(uv_summary$coefficients);
	not_estimable=is.na(tmp_dtfr[,"Estimate"]);
	tmp_dtfr=tmp_dtfr[!not_estimable,];

	# Mark the p-values that are significant
	sig_char=function(x){
		s=character(length(x));
		s[x <= .1]   = ".   ";
		s[x <= .05]  = "*   ";
		s[x <= .01]  = "**  ";
		s[x <= .001] = "*** ";
		s[x <= .0001]= "****";
		return(s);
	}

	signif=sig_char(tmp_dtfr[,4]);
	tmp_dtfr=cbind(tmp_dtfr, signif);
	coeff_txt=capture.output(tmp_dtfr);

	# Sigma
	rse_txt=sprintf(
		"Residual standard error: %5.4f on %i degrees of freedom",
		uv_summary$sigma,
		uv_summary$df[2]
	);

	# F-statistics
	fstat_txt=sprintf(
		"F-statistic: %5.4f on %i and %i DF,  p-value: %5.4f",
		uv_summary$fstatistic[1],
		uv_summary$fstatistic[2],
		uv_summary$fstatistic[3],
		1-pf(uv_summary$fstatistic[1], uv_summary$fstatistic[2], uv_summary$fstatistic[3])
	);

	PowerCalc=F;
	if(PowerCalc){
		print( sorted_taxa_names[var_ix]);
		reg_coef_power(uv_fit[[var_ix]]);
	}
	 
	# Build page contents
	summary_txt=c(
		"Univariate Regression Coefficients for: ", 
		paste("     ", sorted_taxa_names[var_ix], sep=""),
		"",
		coeff_txt,
		"",
		rse_txt,
		fstat_txt
	)
	plot_text(summary_txt);

	# Generate MMPs
	mmps(uv_fit[[var_ix]], main="")

	# Generate sideways histogram
	par(mar=c(5.1,6.1,1,1));
	h=hist(resp_alr[,var_ix], breaks=20, plot=F);
	barplot(h$counts, horiz=T, names.arg=h$mids, space=0, las=2, cex.names=.75,
		ylab="ALR Transformed Abundance", xlab="Counts", main="");
	
	#, main=paste(var_ix, ".) ", sorted_taxa_names[var_ix], 
	#	sprintf(" [%3.1f%%]",mean_abund[var_ix]*100), sep=""));

	setHook("plot.new", NULL, "replace");
}

#############################################################################	
dev.off();

##############################################################################
# Output factor information and MANOVA table

sink(paste(OutputRoot, ".mlr.log.txt", sep=""));

cat("\nFactor information:\n\n");
summary(factors);

cat("\n");
cat(manova_txt, "\n");

sink();

##############################################################################
# Output pvalues 

if(manova_success){
	manova_col=ncol(manova_res);
	manova_row=nrow(manova_res);

	factor_names=attributes(manova_res)$row.names[2:(manova_row-1)];
	pvalues=manova_res[factor_names, "Pr(>F)"];

	pval_fname=paste(OutputRoot, ".pval.tsv", sep="");

	fh=file(pval_fname, "w");

	cat(file=fh, "# Filename", "\t", paste(factor_names, collapse="\t"), "\n", sep="");
	cat(file=fh, OutputRoot, "\t", paste(pvalues, collapse="\t"), "\n", sep="");
	close(fh);
}

##############################################################################
# Output univariate regression coefficient

uv_coeff_fname=paste(OutputRoot, ".uv_coeff.tsv", sep="");
write.table(t(uv_coeff_mat), file=uv_coeff_fname, quote=F, sep="\t", col.names=NA, row.names=T);

##############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
