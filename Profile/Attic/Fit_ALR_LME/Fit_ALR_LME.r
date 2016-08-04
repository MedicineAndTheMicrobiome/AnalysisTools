#!/usr/bin/env Rscript

###############################################################################

NUM_VAR=20;

library(MASS);
library('getopt');
options(useFancyQuotes=F);

params=c(
	"summary_file", "s", 1, "character",
	"factors", "f", 1, "character",
	"model", "m", 2, "character",
	"num_variables", "p", 2, "numeric",
	"reference_levels", "r", 2, "character",
	"outputroot", "o", 2, "character",
	"contains_remaining", "R", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table>\n",
	"	-f <factors>\n",
	"	-m <model formula string>\n",
	"\n",
	"	[-r <reference levels file>]\n",
	"	[-p <number of variables, default = ", NUM_VAR, ">]\n",
	"	[-o <output filename root>]\n",
	"	[-R (pay attention to 'remaining' category)]\n",
	"\n",
	"This script will read in the summary file table, then perform\n",
	"then fit a linear mixed effect model.\n",
	"\n",
	"The -m parameter is required and it specifies the mixed effect\n",
	"model.\n",
	"	For example: \n",
	"		-m \"F1 + F2 + F1*F2\"\n",
	"\n",
	"If the -R flag is set, a 'remaining' category will be be included in the denominator\n",
	"	independent of how large it is.  I.e., do not use it as one of the response variables.\n",
	"\n",
	"The reference levels file specifies which level for the factor should\n",
	"be considered the reference (i.e. 0)\n",
	"\n");

if(!length(opt$summary_file) || !length(opt$factors) || !length(opt$model)){
	cat(usage);
	q(status=-1);
}else{
	SummaryFile=opt$summary_file;
	FactorsFile=opt$factors;
	Model=opt$model;
}

if(!length(opt$outputroot)){
	OutputRoot=gsub(".summary_table.xls$", "", opt$summary_file);
	OutputRoot=gsub(".summary_table.tsv$", "", OutputRoot);
}else{
	OutputRoot=opt$outputroot;
}

if(!length(opt$num_variables)){
	NumVariables=NUM_VAR;
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


cat("\n");
cat("Summary File : ", SummaryFile, "\n", sep="");
cat("Factors File: ", FactorsFile, "\n", sep="");
cat("Model: ", Model, "\n", sep="");
cat("Output File: ", OutputRoot, "\n", sep="");

cat("Number of Response Variables: ", NumVariables, "\n", sep="");
cat("Reference Levels File: ", ReferenceLevelsFile, "\n", sep="");
cat("Use Remaining? ", UseRemaining, "\n");
cat("\n");

options(width=120);
cat("Text Line Width: ", options()$width, "\n", sep="");

##############################################################################

load_factors=function(fname){
	factors=data.frame(read.table(fname, sep="\t", header=TRUE, row.names=1, check.names=FALSE));
	factor_names=colnames(factors);
}

load_summary_file=function(fname){
	inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", quote="", row.names=1))
	counts_mat=inmat[,2:(ncol(inmat))];
	cat("Num    Samples in Summary Table: ", nrow(counts_mat), "\n", sep="");
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

##############################################################################

reconcile_samples=function(factor_table, summary_table){

	cat("Reconciling sample IDs between Factor Table and Summary Table...\n");

	factor_sample_ids=rownames(factor_table);
	sumtab_sample_ids=rownames(summary_table);

	num_factor_sample_ids=length(factor_sample_ids);
	num_sumtab_sample_ids=length(counts_sample_ids);

	#print(factor_sample_id);
	#print(counts_sample_id);

	shared_sample_ids=intersect(factor_sample_ids, sumtab_sample_ids);
	num_shared_sample_ids=length(shared_sample_ids);

	cat("Num summary table sample IDs: ", num_sumtab_sample_ids, "\n");
	cat("       Num factor sample IDs: ", num_factor_sample_ids, "\n");
	cat("       Num shared sample IDs: ", num_shared_sample_ids, "\n");
	cat("\n");

	cat("Samples missing from Summary Table:\n");
	print(setdiff(factor_sample_ids, sumtab_sample_ids));
	cat("\n");

	cat("Samples missing from Factor Table:\n");
	print(setdiff(sumtab_sample_ids, factor_sample_ids));
	cat("\n");

	cat("Total samples shared: ", num_shared_sample_ids, "\n");

	shared_sample_ids=sort(unique(shared_sample_ids));

	reconciled=list():
	reconciled$counts=summary_table[shared_sample_ids,,drop=F];
	reconciled$factors=factor_table[shared_sample_ids,,drop=F];

	return(reconciled);
}

##############################################################################

# Assign 0's to values smaller than smallest abundance across entire dataset
reset_zeros=function(normalized_st){
	min_assay=min(normalized_st[normalized_st!=0]);
	zero_replacment=min_assay/10;

	cat("Lowest non-zero value: ", min_assay, "\n\n", sep="");
	cat("Replacing 0's with: ", zero_replacment, "\n", sep="");

	nozero_mat=normalized_st;
	nozero_mat[normalized_st==0]=zero_replacment;
	return(nozero_mat);
}

##############################################################################

pdf(paste(OutputRoot, ".alr_lme.pdf", sep=""), height=11, width=8.5);

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

if(num_top_taxa>= num_taxa){
	num_top_taxa = (num_taxa-1);
	cat("Number of taxa to work on was changed to: ", num_top_taxa, "\n");
}


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

###############################################################################
# Compute taxonomic correlations and pvalues

uv_fit=list();

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
# Output univariate regression coefficient

uv_coeff_fname=paste(OutputRoot, ".uv_coeff.tsv", sep="");
write.table(t(uv_coeff_mat), file=uv_coeff_fname, quote=F, sep="\t", col.names=NA, row.names=T);

##############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
