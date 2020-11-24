#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
library(car);
options(useFancyQuotes=F);

RM_NA_TRIALS=10000*64;
NUM_TOP_CATEGORIES=30;

params=c(
	"summary_file", "s", 1, "character",
	"factors", "f", 1, "character",
	"num_variables", "p", 2, "numeric",
	"additional_variables_fname", "a", 2, "character",
	"outputroot", "o", 2, "character",

	"reference_levels", "r", 2, "character",
	"model", "m", 2, "character",
	"model_variables_file", "M", 2, "character",
	"contains_remaining", "R", 2, "logical",
	"shorten_category_names", "x", 2, "character",
	"test_run", "T", 2, "logical",
	"rm_na_trials", "N", 2, "numeric",
	"required_var", "q", 2, "character",

	"tag_name", "t", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

script_path=paste(head(strsplit(script_name, "/")[[1]], -1), collapse="/");
source(paste(script_path, "/../../../Metadata/RemoveNAs/Remove_NAs.r", sep=""));

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table>\n",
	"	-f <factors>\n",
	"	[-p <number of variables, in top abundances, default=", NUM_TOP_CATEGORIES, ">]\n",
	"	[-a <additional categories of interest filename>\n",
	"	[-o <output filename root>]\n",
	"\n",
	"	Model building options:\n",
	"	[-r <reference levels file>]\n",
	"	[-m <model formula string>]\n",
	"	[-R (pay attention to 'remaining' category)]\n",
	"	[-M <model variables file>]\n",
	"\n",
	"	NA removal options:\n",
	"	[-N <remove NA trials, trials=", RM_NA_TRIALS, "\n",
	"	[-q <required variables>]\n",
	"\n",
	"	[-x <shorten category names, with separator in double quotes (default=\"\")>]\n",
	"	[-T (test run flag)]\n",
	"\n",
	"	[-t <tag name>]\n",
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
	"\n", sep="");

if(!length(opt$summary_file) || !length(opt$factors)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputRoot=gsub("\\.summary_table\\.xls$", "", opt$summary_file);
	OutputRoot=gsub("\\.summary_table\\.tsv$", "", OutputRoot);
}else{
	OutputRoot=opt$outputroot;
}

if(!length(opt$num_variables)){
	NumVariables=NUM_TOP_CATEGORIES;
}else{
	NumVariables=opt$num_variables;
}

if(!length(opt$additional_variables_fname)){
	AdditionalVariablesFname="";
}else{
	AdditionalVariablesFname=opt$additional_variables_fname;
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

if(length(opt$model_variables_file)){
	ModelVariablesFile=opt$model_variables_file;
}else{
	ModelVariablesFile="";
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

if(length(opt$test_run)){
	TestRun=T;
}else{
	TestRun=F;
}

RequiredFile="";
if(length(opt$required_var)){
	RequiredFile=opt$required_var;
}

Num_Remove_NA_Trials=RM_NA_TRIALS;
if(length(opt$rm_na_trials)){
	Num_Remove_NA_Trials=opt$rm_na_trials;
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
		},
                "append");

}else{
        TagName="";
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
cat("Shorten Category Names: '", ShortenCategoryNames, "'\n", sep="");
cat("\n");

if(ShortenCategoryNames==TRUE){
	cat("Error:  You need to specify a delimitor to split the category names.\n");
	cat("        i.e., this -x option is not a flag, it requires a parameter.\n");
	quit(status=-1);
}

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
	factors=data.frame(read.table(fname,  sep="\t", header=TRUE, row.names=1, check.names=FALSE, comment.char=""));
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

##############################################################################

plot_correl_heatmap=function(mat, title="", noPrintZeros=F, guideLines=F){

	if(is.null(dim(mat))){
		cat(title, " Matrix is NULL.  No heatmap generated.\n");
		return;
	}

	cat("Plotting: ", title, "\n");


	# Generate colors from red to blue
        colors=(rainbow(2^16, start=0, end=0.65));

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

	cscale=min(1, cscale);
	rscale=min(1, rscale);

        max_width=max(nchar(sprintf("%.2f",mat)));
	#cell_cex=sqrt(min(c(cscale, rscale))^2);
	cell_cex=2*(1/max_width)*sqrt(cscale^2 +  rscale^2);

	par(family="Courier");
	par(oma=c(.5, .5, 1.5, .5));

	override_length=10;
	par(mar=c(min(rname_max_len/2, override_length), min(cname_max_len/2, override_length), .5, .5));

	# Remember that rows and columsn are reversed in the image
	image(1:nrow(mat),1:ncol(mat), mat,
                xaxt="n", yaxt="n",
                xlab="", ylab="",
                col=colors
        );

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
	mtext(title, line=0, outer=T, side=3, font=2);

	cat("Done plotting: ", title, "\n");

}

write_top_categorical_effects_by_factor=function(output_fn, coeff_mat, pval_mat, top_n=10){

	cat("Writing top category effects by factor: ", output_fn, "\n");

	top_n=min(top_n, ncol(pval_mat));

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
		if(!is.null(x) && !is.nan(x) && !is.na(x)){
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
		cat("Writing Top ", top_n, "\n");
		ix=1;
		while(ix<=top_n && mat_buf[ix,"ALR"]>0){
			vals=c(paste(ix,"+", sep=""), sort_cat[ix], mat_buf[ix,"ALR"], mat_buf[ix, "p-value"], signif[ix]);
			cat(file=fh, paste(vals, collapse=","), "\n");	
			ix=ix+1;
		}

		# Separator
		cat(file=fh, "...\n");

		# Output Bottom N

		num_cats=nrow(mat_buf);
		cat("Writing Bottom ", top_n, "\n");
		ix=0;
		while((ix<top_n) &&  mat_buf[num_cats-ix,"ALR"]<0 ){
			vals=c(paste(ix+1, "-", sep=""), sort_cat[num_cats-ix], mat_buf[num_cats-ix,"ALR"], 
				mat_buf[num_cats-ix, "p-value"], signif[num_cats-ix]);
			cat(file=fh, paste(vals, collapse=","), "\n");	
			ix=ix+1;
		}

		cat(file=fh, "\n\n");
	}


	close(fh);

}

load_list=function(filename){
        val=scan(filename, what=character(), comment.char="#");
        return(val);
}

##############################################################################

# Open PDF output
pdf(paste(OutputRoot, rnd, ".alr_as_resp.pdf", sep=""), height=11, width=8.5);

# Load summary file table counts 
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

input_info_text=c(
	paste("Summary File: ", SummaryFile, sep=""),
	paste("  Num Samples: ", nrow(counts), sep=""),
	paste("  Num Categories: ", ncol(counts), sep=""),
	"",
	paste("Factor File: ", FactorsFile, sep=""),
	paste("  Num Samples: ", nrow(factors), sep=""),
	paste("  Num Factors: ", ncol(factors), sep=""),
	"",
	paste("Output File Root: ", OutputRoot, sep="") 
);

##############################################################################

# Building Model
if(ModelVariablesFile!=""){
	model_variables_file_list=load_list(ModelVariablesFile);	
	Model=paste(model_variables_file_list, collapse=" + ");
}

if(Model!="All Factors"){
	model_pred_str=Model;
}else{
	model_pred_str=paste(factor_names, collapse=" + ");
}
cat("Model: ", model_pred_str, "\n");


# Remove factors from table that aren't specified in model
model_var_arr=get_var_from_modelstring(model_pred_str);
if(length(model_var_arr)){
	factors=factors[,model_var_arr, drop=F];
	num_factors=ncol(factors);
}

# Load variables to require after NA removal
required_arr=NULL;
if(""!=RequiredFile){
        required_arr=load_list(RequiredFile);
        cat("Required Variables:\n");
        print(required_arr);
        cat("\n");
        missing_var=setdiff(required_arr, factor_names);
        if(length(missing_var)>0){
                cat("Error: Missing required variables from factor file:\n");
                print(missing_var);
        }
}else{
        cat("No Required Variables specified...\n");
}

# Read in additional categories
if(AdditionalVariablesFname!=""){
	cat("Loading Additional ALR Categories...\n");
	additional_categories=load_list(AdditionalVariablesFname);
}else{
	additional_categories=c();
}

plot_text(c(
	"Multivariate ALR Response Regression:",
	"",
	input_info_text,
	"",
	paste("Original Model: ", Model, sep=""),
	"",
	"Required Variables:",
	paste("  File: ", RequiredFile, sep=""),
	paste("  Variable(s): "),
	capture.output(print(required_arr)),
	"",
	"Additional MALR Categories:",
	capture.output(print(additional_categories))
));


##############################################################################
# Remove NAs samples/factors

if(any(is.na(factors))){

	rm_na_res=remove_sample_or_factors_wNA_parallel(factors, required=required_arr, 
		num_trials=Num_Remove_NA_Trials, num_cores=64, outfile=OutputRoot);

	# Update factors, summary table, and model
	factors=rm_na_res$factors;
	num_factors=ncol(factors);
	counts=remove_samples_from_st(rm_na_res, counts);
	model_pred_str=rem_missing_var_from_modelstring(model_pred_str, colnames(factors));
	plot_text(c(
		rm_na_res$summary_text,
		"",
		paste("New Model: ", model_pred_str, sep="")
	));

}

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
##############################################################################


num_factors=ncol(factors);
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

cat("* Reference levels:                                          *\n");
for(f in 1:num_factors){
	cat("*  ", factor_names[f], ":\n*   ", sep="");
	print(levels(factors[,f]));
}
cat("**************************************************************\n");


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
factors=remove_no_variation_factors(factors);
model_pred_str=rem_missing_var_from_modelstring(model_pred_str, colnames(factors)); 
num_factors=ncol(factors);

# Relevel factor levels
if(ReferenceLevelsFile!=""){
        ref_lev_mat=load_reference_levels_file(ReferenceLevelsFile)
        factors=relevel_factors(factors, ref_lev_mat);
}else{
        cat("* No Reference Levels File specified.                        *\n");
}

##############################################################################

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

	normalized_sums=apply(ordered_normalized, 1, sum);
	for(i in 1:num_samples){
		out_mat[i,"Remaining"]=normalized_sums[i]-sum(out_mat[i,]);
	}
	#out_mat[,"Remaining"]=apply(out_mat, 1, function(x){1-sum(x)});

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

plot_overlapping_density=function(mat, title=""){

	cat("Plotting overlapping densities...\n");

	num_cat=ncol(mat);
	range=range(mat);

	cat("Ranges: ", range[1], " - ", range[2], "\n", sep="");

	# Compute the density for each category
	density_list=list();
	max_density=0;
	for(i in 1:num_cat){
		density_list[[i]]=density(mat[,i], n=64);
		max_density=max(max_density, density_list[[i]]$y);
	}
	cat("Max Density: ", max_density, "\n");

	# Open a blank plot
	par(mar=c(5,5,5,1));
	range_span=diff(range);
	plot(0,0, type="n", xlim=c(range[1]-range_span*.3, range[2]+range_span*.3), ylim=c(0, max_density*3), 
		xlab="ALR Value", ylab="Density", main="ALR Density for Extracted Categories (Mode Labelled)");
	title(main=title, outer=F, line=.5, cex.main=.85);
	
	cat_names=colnames(mat);

        colors=rainbow(num_cat, start=0, end=0.65);

	# Plot Densities
	label_pos=numeric(num_cat);
	for(i in 1:num_cat){
		xs=density_list[[i]]$x;
		ys=density_list[[i]]$y;
		max_y=max(ys);
		max_y_ix=max(which(ys==max_y));
		x_at_max_y=xs[max_y_ix];
		label_pos[i]=x_at_max_y;
		points(xs,ys, type="l", col=colors[i], lwd=3);
		points(xs,ys, type="l", col="black", lwd=.5);
		#text(x_at_max_y, max_y, cat_names[i], col=colors[i]);
		points(x_at_max_y, max_y, cex=1, pch=16, col=colors[i]);
	}

	# Tweak label positions so they don't overlap
	sort_ix=order(label_pos);
	label_pos=label_pos[sort_ix]; # Original position
	cat_names=cat_names[sort_ix];
	colors=colors[sort_ix];

	char_size=par()$cxy[1];
	modified=label_pos;	# Tweaked position
	tweaked=T;
	tol=.5;
	while(tweaked){
		tweaked=F;

		max_tweak=max(min(diff(modified)), 0);
		if(max_tweak==0){
			max_tweak=tol/10;
		}
		max_tweak=min(tol/2, max_tweak);

		# Forward adjust
		for(i in 1:(num_cat-1)){
			if(abs(modified[i]-modified[i+1])<tol){
				modified[i+1]=modified[i+1]+max_tweak;	
				tweaked=T;
			}
		}

		# Backward adjust
		for(i in num_cat:2){
			if(abs(modified[i]-modified[i-1])<tol){
				modified[i-1]=modified[i-1]-max_tweak;	
				tweaked=T;
			}
		}

	}	

	# Plot ticks, labels, and connectors
	for(i in 1:num_cat){
		# vertical tick
		points(c(label_pos[i], label_pos[i]), c(max_density*1.05, max_density*1.10), type="l", col=colors[i]);
		# tick to label
		points(c(label_pos[i], modified[i]), c(max_density*1.1, max_density*1.2), type="l", col=colors[i]);
		text(modified[i]-char_size/2, max_density*1.25, cat_names[i], srt=90, pos=4, xpd=T, col=colors[i]);
	}

}


##############################################################################

# Assign 0's to values smaller than smallest abundance across entire dataset
#min_assay=min(normalized[normalized!=0]);
#cat("Lowest non-zero value: ", min_assay, "\n\n", sep="");
#zero_replacment=min_assay/10;
#normalized[normalized==0]=zero_replacment;

##############################################################################

if(num_top_taxa>= num_taxa){
	num_top_taxa = (num_taxa-1);
	cat("Number of taxa to work on was changed to: ", num_top_taxa, "\n");
}


##############################################################################

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

factor_names=colnames(factors);
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

cat("Extracting: ", num_top_taxa, " + 1 (remaining) categories and additional categories.\n", sep="");

# Perform ALR transform
responses=extract_top_categories(normalized, num_top_taxa, additional_cat=additional_categories);
resp_alr_struct=additive_log_rato(responses);
transformed=resp_alr_struct$transformed;

num_cat_to_analyze=ncol(transformed);
sorted_taxa_names=colnames(transformed);
cat("Num ALR Categories to Analyze: ", num_cat_to_analyze, "\n", sep="");

plot_overlapping_density(transformed, title="All");
bottom_half=ceiling(num_cat_to_analyze/2) : num_cat_to_analyze;
top_half=1:floor(num_cat_to_analyze/2);
plot_overlapping_density(transformed[,top_half], title="Top Half by Avg Abundance");
plot_overlapping_density(transformed[,bottom_half], title="Bottom Half by Avg Abundance");

##############################################################################


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

text[1]=paste("Multivariate Regression with ", num_cat_to_analyze, " taxa", sep="");
text[2]=paste("Proportion of top overall mean abundance represented: ", prop_abundance_represented, sep="");
text[3]="";
text=c(text, strsplit(model_string, "(?<=.{80})", perl=T)[[1]]);
text=c(text, "");
text=c(text, manova_txt);
plot_text(text);

###############################################################################
# Compute taxonomic correlations and pvalues

# Compute and Plot Taxonomic correlations
cor_mat=cor(transformed);
print(cor_mat);
plot_correl_heatmap(cor_mat, title="Category Correlations");

# Compute pvalues for correlations, Null Hypothesis is cor=0
cat("Computing P-values for Category Correlations...\n");
pval_matrix=matrix(NA, ncol=num_cat_to_analyze, nrow=num_cat_to_analyze);
colnames(pval_matrix)=colnames(transformed);
rownames(pval_matrix)=colnames(transformed);

pval_vect=numeric(num_cat_to_analyze*(num_cat_to_analyze-1)/2);

num_corr_to_test=0;
for(i in 2:num_cat_to_analyze){
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
fdr_pval_matrix=matrix(NA, ncol=num_cat_to_analyze, nrow=num_cat_to_analyze);
colnames(fdr_pval_matrix)=colnames(transformed);
rownames(fdr_pval_matrix)=colnames(transformed);

num_corr_to_test=0;
for(i in 2:num_cat_to_analyze){
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

uv_pval_mat=matrix(NA, nrow=num_coeff, ncol=num_cat_to_analyze,
		dimnames=list(coeff_names, sorted_taxa_names[1:num_cat_to_analyze]));

uv_coeff_mat=matrix(NA, nrow=num_coeff, ncol=num_cat_to_analyze,
		dimnames=list(coeff_names, sorted_taxa_names[1:num_cat_to_analyze]));


tmp_model_string= paste("transformed[,1] ~", model_pred_str);
test_uv_fit=lm(as.formula(tmp_model_string), data=factors);

anova_factor_names=setdiff(rownames(anova(test_uv_fit)), c("Residuals", "(Intercept)"));
print(anova_factor_names);

uv_anova_pval_mat=matrix(NA, nrow=length(anova_factor_names), ncol=num_cat_to_analyze,
		dimnames=list(anova_factor_names, sorted_taxa_names[1:num_cat_to_analyze]));

uv_model_fit_pval_mat=matrix(NA, ncol=num_cat_to_analyze, nrow=1,
		dimnames=list("p-value", sorted_taxa_names[1:num_cat_to_analyze]));

alr_mean=numeric(num_cat_to_analyze);
alr_stderr=numeric(num_cat_to_analyze);

# Store R^2 for each taxa
rsqrd=numeric(num_cat_to_analyze);
adj_rsqrd=numeric(num_cat_to_analyze);

for(var_ix in 1:num_cat_to_analyze){
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

	alr_mean[var_ix]=mean(ALR_Abundance);
	alr_stderr[var_ix]=sd(ALR_Abundance)/sqrt(length(ALR_Abundance));

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

	# Model p-values
	print(uv_summ$fstatistic);
	uv_model_fit_pval_mat[1, var_ix]=
		1-pf(uv_summ$fstatistic["value"], uv_summ$fstatistic["numdf"], uv_summ$fstatistic["dendf"]);

}

cat("\nUnivariate Regression Coefficients:\n");
print(uv_coeff_mat)
cat("\nUnivariate Regression P-values:\n");
print(uv_pval_mat)
cat("\nUnivariate ANOVA P-values:\n");
print(uv_anova_pval_mat);
cat("\nUnivariate Model P-values:\n");
print(uv_model_fit_pval_mat);

# Plot univariate ANOVA F tests
anova_factors=rownames(uv_anova_pval_mat);
plot_correl_heatmap(uv_anova_pval_mat, title="Univariate F-Tests Pr(>F)", guideLines=T);

# Remove NAs
not_estimable=apply(uv_coeff_mat, 1, function(x){ any(is.na(x))});
coef_names_not_estimable=rownames(uv_coeff_mat)[not_estimable];
uv_coeff_mat=uv_coeff_mat[!not_estimable,, drop=F];
uv_pval_mat=uv_pval_mat[!not_estimable,, drop=F];

# remove intercept
uv_coeff_mat=uv_coeff_mat[-1,, drop=F]; 
uv_pval_mat=uv_pval_mat[-1,, drop=F] 

# Plot pvalues
plot_correl_heatmap(uv_pval_mat, title="Univariate Coefficients Pr(>|t|)");

# Plot log(pvalues, 10)
log_uv_pval_mat=log(uv_pval_mat,10);
plot_correl_heatmap(log_uv_pval_mat, title="Univariate Coefficients Log10[Pr(>|t|)]");


# Plot Heatmap
if(ncol(log_uv_pval_mat)>=2 && nrow(log_uv_pval_mat)>=2 && all(!is.nan(log_uv_pval_mat))){

	# Get current graphic settings so we can restore them.
	par_oma_before=par()$oma;
	par_mar_before=par()$mar;

        cname_max_len=max(nchar(colnames(log_uv_pval_mat)));
        rname_max_len=max(nchar(rownames(log_uv_pval_mat)));

	par(oma=c(1,1,1,1));
	num_hm_colors=20;
	cl_hm_colors=(rainbow(num_hm_colors, start=0, end=0.65));

	override_length=10;
	heatmap(log_uv_pval_mat, col=cl_hm_colors, margins=c(min(rname_max_len/3, override_length), min(cname_max_len/2, override_length)));

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
	cat("No heatmap generated because p-value matrix is not multi-dimensional or all Nan.\n");
}

# Plot log pvalues sorted by most signficiant predictor and taxa
pred_ix=order(apply(log_uv_pval_mat, 1, mean));
taxa_ix=order(apply(log_uv_pval_mat, 2, mean));
plot_correl_heatmap(log_uv_pval_mat[pred_ix, taxa_ix, drop=F], 
	title="Sorted Univariate Coefficients Log10[Pr(>|t|)]", guideLines=T);

# Plot R^2
rsqrd_mat=rbind(rsqrd, adj_rsqrd);
rownames(rsqrd_mat)=c("R^2", "Adjusted R^2");
colnames(rsqrd_mat)=sorted_taxa_names[1:num_cat_to_analyze];
plot_correl_heatmap(rsqrd_mat, title="Univariate R Squared");

# Plot Univariate Model F-stat
plot_correl_heatmap(uv_model_fit_pval_mat, title="Univariate Model Fit F-stat P-values");

# Plot univariate coefficients
not_na_coeff=apply(uv_coeff_mat, 1, function(x){!any(is.na(x))});
plot_correl_heatmap(uv_coeff_mat[not_na_coeff,, drop=F], title="Univariate Coefficients", guideLines=T);

mask_matrix=function(val_mat, mask_mat, mask_thres, mask_val){
        masked_matrix=val_mat;
        masked_matrix[mask_mat>mask_thres]=mask_val;
        return(masked_matrix);
}

# Plot significant coefficients at various pvalue cutoffs
signf_coef=mask_matrix(uv_coeff_mat, uv_pval_mat, .1, 0);
plot_correl_heatmap(signf_coef, "Significant Coefficients at p-value < 0.10", noPrintZeros=T, guideLines=T);
signf_coef=mask_matrix(uv_coeff_mat, uv_pval_mat, .05, 0);
plot_correl_heatmap(signf_coef, "Significant Coefficients at p-value < 0.05", noPrintZeros=T, guideLines=T);
signf_coef=mask_matrix(uv_coeff_mat, uv_pval_mat, .01, 0);
plot_correl_heatmap(signf_coef, "Significant Coefficients at p-value < 0.01", noPrintZeros=T, guideLines=T);


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


# Write coefficient p-values to file
write.table(t(uv_pval_mat), file=paste(OutputRoot, ".alr_as_resp.pvals.tsv", sep=""), 
	sep="\t", quote=F, col.names=NA, row.names=T);

write.table(t(uv_coeff_mat), file=paste(OutputRoot, ".alr_as_resp.coefs.tsv", sep=""),
	sep="\t", quote=F, col.names=NA, row.names=T);

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

setHook("plot.new", function(){mtext(sorted_taxa_names[var_ix], outer=T, line=-.5);}, "prepend");
hooks=getHook("plot.new");

for(var_ix in 1:num_cat_to_analyze){

	if(length(getHook("plot.new"))==0){
		for(hix in 1:length(hooks)){
			setHook("plot.new", hooks[[hix]], "prepend");
		}
	}

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
		"Intercept-only Fit:",
		paste("  ALR Mean: ", round(alr_mean[var_ix], 4)),
		paste("  ALR Standard Error: ", round(alr_stderr[var_ix], 4)),
		"",
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
	h=hist(resp_alr_struct$transformed[,var_ix], breaks=20, plot=F);
	barplot(h$counts, horiz=T, names.arg=signif(h$mids, 3), space=0, las=2, cex.names=.75,
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
# Write MANOVA pvalues to file

if(manova_success){
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
write.table(outmat, file=paste(OutputRoot, ".alr_as_resp.anova.summary.tsv", sep=""), sep="\t", quote=F, col.names=T, row.names=F);

##############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
