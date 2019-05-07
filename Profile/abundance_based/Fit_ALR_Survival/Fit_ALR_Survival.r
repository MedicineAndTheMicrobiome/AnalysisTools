#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
library(survival);

source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');

options(useFancyQuotes=F);

params=c(
	"summary_file", "s", 1, "character",
	"num_top_pred", "v", 2, "numeric",
	"contains_remaining", "R", 2, "logical",
	"shorten_category_names", "x", 2, "character",
	
	"factor_file", "f", 1, "character",
	"model_vars", "m", 2, "character",
	"time_varname", "t", 1, "character",
	"subj_varname", "i", 1, "character",
	"coht_varname", "c", 1, "character",

	"required_varname", "q", 2, "character",
	"reference_levels", "r", 2, "character",
	
	"timeofdeath_file", "d", 1, "character",
	"epoch_file", "e", 2, "character",

	"pvalue_cutoff", "p", 2, "numeric",
	"output_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

NUM_TOP_PRED_CAT=20;
PVAL_CUTOFF=.1;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table>\n",
	"	[-v <number of top predictor (as ALR) categories to include, default=", NUM_TOP_PRED_CAT, ">]\n",
	"	[-R (pay attention to 'remaining' category)]\n",
	"	[-x <shorten category names, with separator in double quotes (default=\"\")>]\n",
	"\n",
	"	-f <factors file, contains covariates and factors>\n",
	"	[-m <list of model variables (file name) to include from the factor file>]\n",
	"	-t <time point variable name (from time 0)>\n",
	"	-i <subject/individual identifier variable name (i.e. grouping of individuals)>\n",
	"	-c <cohort identifier variable name (e.g. control, treatment, disease)>\n",

	"	[-q <required list of variables (file name) to include after NA removal>]\n",
	"	[-r <reference levels (file name) for Y's in factor file>]\n",
	"\n",
	"	-d <time of death/status change file>\n",
	"	[-e <epoch definitions file>\n",
	"\n",
	"	[-p <pvalue cutoff, default=", PVAL_CUTOFF, ">]\n",
	"	-o <output filename root>\n",
	"\n",
	"\n",
	"The format of the 'death/status change file', is:\n",
	"	<subject ID>\\t<time of death>\\n\n",
	"Time of death, should match the time scale used in the -t option.\n",
	"\n",
	"1.) Plot Kaplan-Meier Estimator Curve\n",
	"2.) Calculate Cox Proportional Hazard\n",
	"3.) Pairwise Log-rank test between cohorts (-c option)\n",
	"\n",
	"If the -R flag is set, a 'remaining' category will be be included in the denominator\n",
	"	independent of how large it is.  I.e., do not use it as one of the response variables.\n",
	"\n",
	"The epoch file (-e) should be of the follow format.  Note that all ranges are inclusive.\n",
	"# Begin  End  Name\n",
	"      0    1  Start\n",
	"     1.1   2  Second\n",
	"     ...\n",
	"\n",	
	"\n", sep="");

if(
	!length(opt$summary_file) || 
	!length(opt$factor_file) || 
	!length(opt$time_varname) || 
	!length(opt$subj_varname) || 
	!length(opt$coht_varname) || 
	!length(opt$timeofdeath_file) || 
	!length(opt$output_root)
	){
	cat(usage);
	q(status=-1);
}

# Required
SummaryFile=opt$summary_file;
FactorFile=opt$factor_file;
TimeVarName=opt$time_varname;
SubjVarName=opt$subj_varname;
CohtVarName=opt$coht_varname;
TimeOfDeathFile=opt$timeofdeath_file;
OutputRoot=opt$output_root;

# Optional, i.e. with defaults
NumALRPredictors=NUM_TOP_PRED_CAT;
UseRemaining=F;
ShortenCategoryNames="";
RequiredFile="";
ReferenceLevelsFile="";
ModelVarFile="";
PvalCutoff=PVAL_CUTOFF;
EpochFile="";

if(length(opt$num_top_pred)){
	NumALRPredictors=opt$num_top_pred;
}

if(length(opt$contains_remaining)){
	UseRemaining=T;
}

if(length(opt$shorten_category_names)){
	ShortenCategoryNames=opt$shorten_category_names;
}

if(length(opt$model_vars)){
	ModelVarFile=opt$model_vars;
}

if(length(opt$required_varname)){
	RequiredFile=opt$required_varname;
}

if(length(opt$reference_levels)){
        ReferenceLevelsFile=opt$reference_levels;
}

if(length(opt$pvalue_cutoff)){
	PvalCutoff=opt$pvalue_cutoff;	
}

if(length(opt$epoch_file)){
	EpochFile=opt$epoch_file;
}

###############################################################################

input_param=capture.output({
	cat("\n");
	cat("Summary File: ", SummaryFile, "\n", sep="");
	cat("  Num Top ALR Predictors: ", NumALRPredictors, "\n", sep="");
	cat("  Contains 'Remaining': ", UseRemaining, "\n", sep="");
	cat("  Shorten Categories: ", ShortenCategoryNames, "\n", sep="");
	cat("\n");
	cat("Factor File: ", FactorFile, "\n", sep="");
	cat("  Model Variables File: ", ModelVarFile, "\n", sep="");
	cat("  Time Variable Name: ", TimeVarName, "\n", sep="");
	cat("  Subject Variable Name: ", SubjVarName, "\n", sep="");
	cat("  Cohort Variable Name: ", CohtVarName, "\n", sep="");
	cat("  Required Variables File: ", RequiredFile, "\n", sep="");
	cat("  Reference Levels File: ", ReferenceLevelsFile, "\n", sep="");
	cat("\n");
	cat("Time of Death/Event File: ", TimeOfDeathFile, "\n", sep="");
	cat("Epoch File: ", EpochFile, "\n", sep="");
	cat("P-value Cutoff: ", PvalCutoff, "\n", sep="");
	cat("Output File Root: ", OutputRoot, "\n", sep="");
	cat("\n");
});

cat(paste(input_param, collapse="\n"));


if(ShortenCategoryNames==TRUE){
	cat("Error:  You need to specify a delimitor to split the category names.\n");
	cat("        i.e., this -x option is not a flag, it requires a parameter.\n");
	quit(status=-1);
}

options(width=100);
cat("Text Line Width: ", options()$width, "\n", sep="");

##############################################################################
##############################################################################

load_factors=function(fname){
	factors=data.frame(read.table(fname,  sep="\t", header=TRUE, row.names=1, 
		check.names=FALSE, comment.char=""));
	factor_names=colnames(factors);

	ignore_idx=grep("^IGNORE\\.", factor_names);

	if(length(ignore_idx)!=0){
		return(factors[-ignore_idx]);
	}else{
		return(factors);
	}
}

load_summary_file=function(fname){
	inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, 
		comment.char="", quote="", row.names=1))
	counts_mat=inmat[,2:(ncol(inmat))];

	# Clean category names a little
	cat_names=colnames(counts_mat);
	cat_names=gsub("-", "_", cat_names);
	colnames(counts_mat)=cat_names;
	
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

load_death_times=function(filename){
        inmat=as.matrix(read.table(filename, sep="\t", header=T, check.names=FALSE, comment.char="#", row.names=1))
	dtimes=as.numeric(inmat);
	names(dtimes)=rownames(inmat);
	return(dtimes);
}

load_epochs=function(filename){
        inmat=as.matrix(read.table(filename, sep="\t", header=F, check.names=FALSE, comment.char="#"))
	epochs=list();
	num_rows=nrow(inmat);
	cat("Loading: ", num_rows, " rows from ", filename, "\n");
	for(i in 1:num_rows){
		epochs[[inmat[i,3]]]=as.numeric(inmat[i,c(1,2)]);		
	}
	return(epochs);	
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

        num_row=nrow(mat);
        num_col=ncol(mat);

	if(num_row==1){
		plot_row_dendr=F;
	}
	if(num_col==1){
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
pdf(paste(OutputRoot, ".survival.pdf", sep=""), height=14, width=8.5);

# Load summary file table counts 
cat("\n");
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

num_st_categories=ncol(counts);
num_st_samples=nrow(counts);

cat("Num Categories: ", num_st_categories, "\n");
cat("   Num Samples: ", num_st_samples, "\n");

##############################################################################

# Shorten cateogry names
if(ShortenCategoryNames!=""){
	full_names=colnames(counts);
	splits=strsplit(full_names, ShortenCategoryNames);
	short_names=character();
	for(i in 1:length(full_names)){
		short_names[i]=tail(splits[[i]], 1);
		short_names[i]=gsub("_unclassified$", "_uncl", short_names[i]);
	}
	colnames(counts)=short_names;
	cat("Names have been shortened.\n");
}else{
	cat("Keeping original category names...\n");
}

##############################################################################

# Load factors
cat("Loading Factors...\n");
factors=load_factors(FactorFile);
factor_names=colnames(factors);
num_factors=ncol(factors);
factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);

num_loaded_factors=num_factors;
num_loaded_factor_samp=num_factor_samples;

cat("\n");
cat(num_factors, " Factor(s) Loaded:\n", sep="");
print(factor_names);
cat("\n");

# Load predictors to include in model
model_var_arr=c(TimeVarName, SubjVarName, CohtVarName);
model_pred=c();
if(ModelVarFile!=""){
	model_pred=load_list(ModelVarFile);
	cat("Model Variables:\n");
	print(model_pred);
	cat("\n");
}
model_var_arr=c(model_var_arr, model_pred);

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

plot_text(c(
	"Variables Targeted:",
	"",
	"Predictors:",
	capture.output(print(model_var_arr)),
	"",
	"Required Variables:",
	 capture.output(print(required_arr))
));

# Confirm we can find all the factors
missing_fact=setdiff(model_var_arr, factor_names);
if(length(missing_fact)>0){
	cat("Error: Factors in model, missing Factor File.\n");
	print(missing_fact);
	quit(status=-1);
}else{
	cat("All model variables found in factor file...\n");
}
factors=factors[,model_var_arr, drop=F];
factor_names=colnames(factors);
num_factors=ncol(factors);

# Relevel factor levels
if(ReferenceLevelsFile!=""){
        ref_lev_mat=load_reference_levels_file(ReferenceLevelsFile)
        factors=relevel_factors(kept_factors, ref_lev_mat);
}else{
        cat("No Reference Levels File specified.\n");
}

##############################################################################
# Reconcile factors with samples


reconcile_summarytable_factors=function(st, fc){

	cat("\nReconciling samples between summary table and factor file...\n");
	st_ids=rownames(st);
	fc_ids=rownames(fc);	

	num_st_ids=length(st_ids);
	num_fc_ids=length(fc_ids);

	cat("Num IDs: SummmaryTable: ", num_st_ids, " Factors: ", num_st_ids, "\n");

	shared_ids=intersect(st_ids, fc_ids);

	uniq_to_st=setdiff(st_ids, shared_ids);
	uniq_to_fc=setdiff(fc_ids, shared_ids);

	num_shared=length(shared_ids);
	cat("Shared IDs: ", num_shared, "\n");

	result=list();
	result[["summary_table"]]=st[shared_ids,,drop=F];
	result[["factors"]]=fc[shared_ids,,drop=F];
	result[["unique_to_st"]]=uniq_to_st;
	result[["unique_to_fc"]]=uniq_to_fc;
	result[["shared"]]=shared_ids;
	return(result);
}

##############################################################################
# Remove samples with NAs

recon_res=reconcile_summarytable_factors(counts, factors);
recon_factors=recon_res[["factors"]];

factors_wo_nas_res=remove_sample_or_factors_wNA_parallel(recon_factors, 
	required=required_arr, num_trials=64000, num_cores=64, outfile=paste(OutputRoot, ".noNAs", sep=""));

factors_wo_nas=factors_wo_nas_res$factors;

recon_res=reconcile_summarytable_factors(counts, factors_wo_nas);
counts=recon_res[["summary_table"]];
factors=recon_res[["factors"]];

##############################################################################
# Prepping for ALR calculations

if(NumALRPredictors >= num_st_categories){
	NumALRPredictors= (num_st_categories-1);
	cat("Number of taxa to work on was changed to: ", NumALRPredictors, "\n");
}

##############################################################################

# Normalize
cat("Normalizing counts...\n");
normalized=normalize(counts);

cat("Reordering normalized...\n");
mean_norm=apply(normalized, 2, mean);
ord_ix=order(mean_norm, decreasing=T);
normalized=normalized[,ord_ix, drop=F];
counts=counts[,ord_ix, drop=F];

# Assign 0's to values smaller than smallest abundance across entire dataset
min_assay=min(normalized[normalized!=0]);
cat("Lowest non-zero value: ", min_assay, "\n", sep="");
zero_replacment=min_assay/10;
cat("Substituting 0's with: ", zero_replacment, "\n", sep="");
normalized[normalized==0]=zero_replacment;

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
}else{
	cat("Assuming no categories called 'remainder' or 'remaining'\n");
}

##############################################################################

cat("\n");
cat("Extracting Top categories: ", NumALRPredictors, " from amongst ", ncol(normalized), "\n", sep="");

cat_abundances=extract_top_categories(normalized, NumALRPredictors);
resp_alr_struct=additive_log_rato(cat_abundances);
alr_categories_val=resp_alr_struct$transformed;
alr_cat_names=colnames(alr_categories_val);

plot_text(c(
	paste("ALR Categories (Top ", NumALRPredictors, ")", sep=""),
	capture.output(print(alr_cat_names))
));

##############################################################################

death_times=load_death_times(TimeOfDeathFile);
cat("Death Times:\n");
print(death_times);

##############################################################################

times=factors[,TimeVarName];
subj_ids=factors[,SubjVarName];
coht_ids=factors[,CohtVarName];

uniq_times=sort(unique(times));
uniq_subj_ids=sort(unique(subj_ids));
uniq_coht_ids=sort(unique(coht_ids));

last_measured_time=max(death_times, times, na.rm=T); 

cat("\nTime Points:\n");
print(uniq_times);

cat("\nSubject IDs:\n");
print(uniq_subj_ids);

cat("\nCohort IDs:\n");
print(uniq_coht_ids);

cat("\n");
cat("Last Measured Time: ", last_measured_time, "\n");

###############################################################################

cat("Building data structures for cohort/subject/times...\n");

group_structure=list();

for(sbj_ix in uniq_subj_ids){
	row_ix=(factors[,SubjVarName]==sbj_ix);
	time_ord=order(factors[row_ix, TimeVarName]);
	fact_sub=factors[row_ix,,drop=F];
	ordered=fact_sub[time_ord,,drop=F];
	#print(ordered);

	coht=as.character(unique(ordered[,CohtVarName]));
	if(is.null(group_structure[[as.character(coht)]])){
		group_structure[[coht]]=list();
	}

	group_structure[[coht]][[sbj_ix]]=list();
	group_structure[[coht]][[sbj_ix]][["sample_id"]]=rownames(ordered);
	group_structure[[coht]][[sbj_ix]][["times"]]=ordered[,TimeVarName];
	group_structure[[coht]][[sbj_ix]][["death"]]=death_times[sbj_ix];
}

#print(group_structure);

###############################################################################

plot_alr_over_time=function(ids, times, end_time, alrs, time_range, alr_range, colors, title="", keep=NULL){

	max_m=3;
	min_m=.25;
	mult=max_m-min_m;
	
	par(mar=c(2, 4, .5, 4));

	pad=function(rng, amt){
		pad=abs(diff(rng))*amt;
		return(c(rng[1]-pad, rng[2]+pad));
	}
	
	xrang=pad(time_range,.05);
	yrang=pad(alr_range,.1);

	cat_names=colnames(alrs);

	plot(0,0, type="n", xlim=xrang, ylim=yrang, ylab=title); 
	
	num_alrs=ncol(alrs);
	for(cat_ix in 1:num_alrs){
		if(is.null(keep) || any(cat_names[cat_ix]==keep)){
			points(times, alrs[ids, cat_ix], type="l", 
				col=colors[cat_names[cat_ix]], lwd=max_m-(cat_ix/num_alrs)*mult);
		}
		#points(times, alrs[ids, cat_ix], type="l", col="black", lwd=.05);
	}
	for(cat_ix in 1:num_alrs){
		if(is.null(keep) || any(cat_names[cat_ix]==keep)){
			points(times, alrs[ids, cat_ix], type="p", pch=16, 
				cex=max_m-(cat_ix/num_alrs)*mult, col=colors[cat_ix]);
		}
	}
	
	for(et in end_time){
		if(!is.na(et)){
			abline(v=et, col="grey");
		}
	}

}

plot_death_over_time=function(deaths, time_range, title){


	pad=function(rng, amt){
		pad=abs(diff(rng))*amt;
		return(c(rng[1]-pad, rng[2]+pad));
	}
	xrang=pad(time_range, .05);

	num_subj=length(deaths);
	deaths=deaths[order(deaths)];

	cat("Deaths: \n");
	print(deaths);

	deaths_noNA=deaths[!is.na(deaths)];

	uniq_times=unique(c(0,deaths_noNA));	
	num_utimes=length(uniq_times);

	alive=numeric(num_utimes);
	times=numeric(num_utimes);

	cur_alive=num_subj;
	
	for(i in 1:num_utimes){

		curtime=uniq_times[i];	
		num_died=sum(curtime==deaths_noNA)

		cur_alive=cur_alive-num_died;
		alive[i]=cur_alive;
						
		times[i]=curtime;
	}

	prop_alive=alive/num_subj;

	prop_alive=c(prop_alive, prop_alive[num_utimes], prop_alive[num_utimes]);
	times=c(times, times[num_utimes], time_range[2]);


	par(mar=c(2, 4, .5, 4));
	plot(0,0, type="n", xlim=xrang, ylim=c(-.05, 1.05), ylab=title); 
	axis(4, at=(0:num_subj)/num_subj, labels=0:num_subj, las=2);
	for(i in 1:length(times)){
		points(c(times[i], times[i+1]), c(prop_alive[i], prop_alive[i]), type="l");
	}
	

	
}

plot_category_colors=function(alrs, colors, keep=NULL){
	
	if(!is.null(keep)){
		alrs=alrs[,keep,drop=F];
	}

	cat("Plotting category colors...\n");
	names=colnames(alrs);
	num_cat=ncol(alrs);

	denom_scale=max(num_cat, 10);

	par(mar=c(2, 4, .5, 4));
	plot(0,0, type="n", xlim=c(0,10), ylim=c(0,1), ylab="", xlab="", bty="n", xaxt="n", yaxt="n");
	print(names);
	for(i in 1:num_cat){
		points(0, 1-(i/denom_scale), col=colors[names[i]], pch=16);
		text(0, 1-(i/denom_scale), names[i], pos=4);
	}
	

}

get_colors=function(num_col, alpha=1){
        colors=hsv(seq(0,1,length.out=num_col+1), c(1,.5), c(1,.75,.5), alpha=alpha);
        color_mat_dim=ceiling(sqrt(num_col));
        color_pad=rep("grey", color_mat_dim^2);
        color_pad[1:num_col]=colors[1:num_col];
        color_mat=matrix(color_pad, nrow=color_mat_dim, ncol=color_mat_dim);
        colors=as.vector(t(color_mat));
        colors=colors[colors!="grey"];
}

###############################################################################

alr_range=range(alr_categories_val);

#colors=rainbow(NumALRPredictors, end=0.66);
colors=get_colors(NumALRPredictors);
names(colors)=colnames(alr_categories_val);

plots_per_page=6;

par(oma=c(1,0,4,1));

for(coht_ix in uniq_coht_ids){
	par(mfrow=c(plots_per_page,1));

	cat("Working on Cohort: ", coht_ix, "\n");
	subj_in_coht=names(group_structure[[coht_ix]])
	plot_ix=1;

	for(subj_ix in subj_in_coht){

		cat("  Subject: ", subj_ix, "\n");
		samp_ids=group_structure[[coht_ix]][[subj_ix]][["sample_id"]];
		times=group_structure[[coht_ix]][[subj_ix]][["times"]];
		death=group_structure[[coht_ix]][[subj_ix]][["death"]];
		plot_alr_over_time(samp_ids, times, death, 
			alr_categories_val[samp_ids,,drop=F], 
			c(0, last_measured_time), alr_range, colors,
			title=subj_ix);

		if(plot_ix==1){
			mtext(coht_ix, side=3, line=0, outer=T, font=2, cex=2);
		}

		if(plot_ix==plots_per_page-1 || subj_ix==tail(subj_in_coht,1)){
			plot_category_colors(alr_categories_val, colors);
			plot_ix=1;
		}else{
			plot_ix=plot_ix+1;
		}
	}
}

# Collapse individuals into groups
alr_by_coht=list();
num_alr_cat=ncol(alr_categories_val);

for(coht_ix in uniq_coht_ids){

	cat("Working on Cohort: ", coht_ix, "\n");
	subj_in_coht=names(group_structure[[coht_ix]])
	num_subj_in_coht=length(subj_in_coht);

	alr_by_coht[[coht_ix]]=list();
	alr_by_coht[[coht_ix]][["death"]]=numeric();
	alr_by_coht[[coht_ix]][["alr"]]=list();
	for(time_ix in uniq_times){
		alr_by_coht[[coht_ix]][["alr"]][[as.character(time_ix)]]=numeric();
	}

	# Rearrange samples by time
	for(subj_ix in subj_in_coht){
		times=group_structure[[coht_ix]][[subj_ix]][["times"]];
		samp_ids=group_structure[[coht_ix]][[subj_ix]][["sample_id"]];
		death=group_structure[[coht_ix]][[subj_ix]][["death"]];
		alr_by_coht[[coht_ix]][["death"]]=c(alr_by_coht[[coht_ix]][["death"]], death);
		
		num_times=length(times);

		for(time_ix in 1:num_times){
			cur_time=as.character(times[time_ix]);	
			cur_samp=samp_ids[time_ix];

			old_name=rownames(alr_by_coht[[coht_ix]][["alr"]][[cur_time]]);	

			alr_by_coht[[coht_ix]][["alr"]][[cur_time]]=rbind(
				alr_by_coht[[coht_ix]][["alr"]][[cur_time]],
				alr_categories_val[cur_samp,]
				);

			rownames(alr_by_coht[[coht_ix]][["alr"]][[cur_time]])=c(
				old_name, cur_samp);
			
		}
	}
}

#print(alr_by_coht);

avg_alr_by_coht=list();

num_uniq_time_pts=length(uniq_times);

for(coht_ix in uniq_coht_ids){
	cat("Working on Cohort: ", coht_ix, "\n");
	cur_coht=alr_by_coht[[coht_ix]];

	tmp_avg=matrix(NA, nrow=num_uniq_time_pts, ncol=num_alr_cat);
	rownames(tmp_avg)=uniq_times;
	colnames(tmp_avg)=colnames(alr_categories_val);

	for(time_ix in as.character(uniq_times)){
		tmp_avg[time_ix,]=apply(cur_coht[["alr"]][[time_ix]], 2, mean);
	}

	avg_alr_by_coht[[coht_ix]]=tmp_avg;

}
print(avg_alr_by_coht);


plots_per_age=6;
par(mfrow=c(plots_per_page,1));
plot_ix=0;
for(coht_ix in uniq_coht_ids){

	plot_alr_over_time(
		as.character(uniq_times), 
		uniq_times, 
		alr_by_coht[[coht_ix]][["death"]],
		avg_alr_by_coht[[coht_ix]],
		c(0, last_measured_time), alr_range, colors,
		title=coht_ix);
	
	plot_death_over_time(
		alr_by_coht[[coht_ix]][["death"]], c(0, last_measured_time), coht_ix
	);

	if((plot_ix %% plots_per_page)==0){
		mtext("Combined by Cohort", side=3, line=0, outer=T, font=2, cex=2);
	}

	plot_ix=plot_ix+1;
}

###############################################################################

# Add subject ID and Time to ALR
samp_ids=rownames(factors);
factors_w_alr=cbind(factors, alr_categories_val[samp_ids,]);

# Estimate last study time
last_study_time=max(c(death_times, factors[,TimeVarName]), na.rm=T);

cat("Last Recorded Time of Study: ", last_study_time, "\n");
last_study_time=last_study_time+1;
cat("  Using last study time of: ", last_study_time, "\n");

# Add death 
event_df=as.data.frame(matrix(0, nrow=length(death_times), ncol=3));
colnames(event_df)=c(TimeVarName, SubjVarName, "Status");
rownames(event_df)=names(death_times);

for(sbj in names(death_times)){
	if(is.na(death_times[sbj])){
		dora=0;
		st_time=last_study_time;
	}else{
		dora=1;
		st_time=death_times[sbj];
	}
	event_df[sbj,]=c(st_time, sbj, dora);
}


change_variable_name=function(current, new, df){
	cat("Changing variable name from: ", current, " to: ", new, "\n");
	cur_names=colnames(df);
	ix=which(cur_names==current);
	cur_names[ix]=new;
	colnames(df)=cur_names;
	return(df);
}

times_ge_zero=factors_w_alr[,TimeVarName]>0;
times_eq_zero=factors_w_alr[,TimeVarName]==0;


print(head(event_df));
print(head(factors_w_alr));


merge_events_w_factors=function(fact_alr, events, time_var_name, subj_var_name){
	# Add tstart, tstop, and endpt 
	
	num_rows=nrow(fact_alr);
	
	tstart=rep(0, num_rows);
	tend=rep(0, num_rows);
	endpt=rep(0, num_rows);
	appended_table=cbind(tstart, tend, endpt, fact_alr);

	subjects=unique(events[,subj_var_name]);
	cat("Subjects:\n");
	print(subjects);
	
	out_table=c();

	for(sbj_id in subjects){
		#cat("\n\nWorking on: ", sbj_id, "\n");
		rows_ix=(appended_table[,subj_var_name]==sbj_id);
		sbj_data=appended_table[rows_ix,];
		times=sbj_data[,time_var_name];
		sort_ix=order(times);
		sbj_data=sbj_data[sort_ix,];

		sbj_event=events[sbj_id,];

		num_times=nrow(sbj_data);
		#cat("Num Time Pts: ", num_times, "\n");

		# Set status
		if(sbj_event["Status"]==1){
			sbj_data[num_times, "endpt"]=1;
		}

		# Set tend to the measure time
		for(i in 1:num_times){
			sbj_data[i, "tstart"]=sbj_data[i, time_var_name];
		}

		for(i in 1:(num_times-1)){
			sbj_data[i, "tend"]=sbj_data[i+1, "tstart"];
		}

		sbj_data[num_times, "tend"]=as.numeric(sbj_event[time_var_name]);

		#print(sbj_event);
		#print(sbj_data);

		out_table=rbind(out_table, sbj_data);
	}

	return(out_table);
}

merged_table=merge_events_w_factors(factors_w_alr, event_df, TimeVarName, SubjVarName);

#print(merged_table);

tstarts=merged_table[,"tstart"];
tend=merged_table[,"tend"];
degen_ix=which(tstarts==tend);

if(length(degen_ix)){
	cat("Adjusting degenerate time...\n");

	min_time_step=min(diff(uniq_times));
	time_adj=min_time_step/10;

	cat("Adding: ", time_adj, "\n");

	for(dix in degen_ix){
		merged_table[,"tend"]=merged_table[,"tend"]+time_adj;
	}
}

###############################################################################

cat("Fitting Full Model: Covariates + ALR Categories\n");
full_form_string=paste("Surv(tstart, tend, endpt)~", paste(c(model_pred, CohtVarName, alr_cat_names), collapse="+"));
cat("Cox Proportional Hazards Formula: ", full_form_string, "\n");
full_surv_form=as.formula(full_form_string);
full_coxph_fit=coxph(full_surv_form, data=merged_table);
full_coxph_fit_summ=summary(full_coxph_fit);

cat("Fitting Reduced Model: Covariates (only)\n");
reduced_form_string=paste("Surv(tstart, tend, endpt)~", paste(c(model_pred, CohtVarName), collapse="+"));
cat("Cox Proportional Hazards Formula: ", reduced_form_string, "\n");
reduced_surv_form=as.formula(reduced_form_string);
reduced_coxph_fit=coxph(reduced_surv_form, data=merged_table);
reduced_coxph_fit_summ=summary(reduced_coxph_fit);


full_coxph_res_text=capture.output(summary(full_coxph_fit));
reduced_coxph_res_text=capture.output(summary(reduced_coxph_fit));

###############################################################################

notes=c(
	"Positive 'coef' increase risk, 0 means no change in risk.",
	"exp(coef) are the the hazard ratios. A hazard ratio of 1 means no change.",
	"Likelihood, Wald and Score are global stat signf of model. They should converge for large N."
);

par(mfrow=c(1,1));

par(oma=c(0,0,3,0));
plot_text(c(
	notes,
	"",
	full_coxph_res_text)
);
mtext("Full Model", outer=T, font=2, cex=2);

par(oma=c(0,0,3,0));
plot_text(c(
	notes,
	"",
	reduced_coxph_res_text)
);
mtext("Reduced Model", outer=T, font=2, cex=2);

###############################################################################

reduced_rsq=reduced_coxph_fit_summ$rsq;
full_rsq=full_coxph_fit_summ$rsq;

reduced_concord=reduced_coxph_fit_summ$concordance["C"];
full_concord=full_coxph_fit_summ$concordance["C"];

par(mfrow=c(2,1));
par(oma=c(1,1,3,1));
par(mar=c(4,4,4,1));
barplot(c(full_rsq[1], full_rsq[2], reduced_rsq[1], reduced_rsq[2]), 
	names.arg=c("Full Model", "Full Maximum", "Reduced Model", "Reduced Maximum"),
	ylim=c(0,1),
	main="R^2"
);

barplot(c(full_concord, reduced_concord),
	names.arg=c("Full Model", "Reduced Model"),
	ylim=c(0,1),
	main="Concordance (Ability to predict order of events)"
);

###############################################################################

number_cohts=length(uniq_coht_ids);

par(mfrow=c(number_cohts,2));
par(mar=c(4,4,1,1));

lowess_interval=min(diff(uniq_times))/2;

for(coht_ix in uniq_coht_ids){
	cat("Subsetting for: ", coht_ix, "\n");
	sub_merged=merged_table[merged_table[,CohtVarName]==coht_ix,, drop=F];

	# Survival probability
	full_pred=predict(full_coxph_fit, sub_merged, type="expected");
	reduced_pred=predict(reduced_coxph_fit, sub_merged, type="expected");

	plot(sub_merged[,"tend"], exp(-full_pred), ylim=c(0,1), main="Full", ylab=coht_ix, xlab="");
	points(lowess(sub_merged[,"tend"], exp(-full_pred), delta=lowess_interval), type="l", col="blue");

	plot(sub_merged[,"tend"], exp(-reduced_pred), ylim=c(0,1), main="Reduced", ylab="", xlab="");
	points(lowess(sub_merged[,"tend"], exp(-reduced_pred), delta=lowess_interval), type="l", col="blue");

	mtext("Survival Probability", outer=T, font=2, cex=1.2);
}

for(coht_ix in uniq_coht_ids){

	cat("Subsetting for: ", coht_ix, "\n");
	sub_merged=merged_table[merged_table[,CohtVarName]==coht_ix,, drop=F];

	# Expected number of events
	full_pred=predict(full_coxph_fit, sub_merged, type="expected");
	reduced_pred=predict(reduced_coxph_fit, sub_merged, type="expected");

	plot(sub_merged[,"tend"], full_pred, main="Full", xlab="", ylab=coht_ix, ylim=c(0, 3));
	points(lowess(sub_merged[,"tend"], full_pred, delta=lowess_interval), type="l", col="blue");

	plot(sub_merged[,"tend"], reduced_pred, main="Reduced", xlab="", ylab="", ylim=c(0,3));
	points(lowess(sub_merged[,"tend"], reduced_pred, delta=lowess_interval), type="l", col="blue");

	mtext("Expected Number of Events", outer=T, font=2, cex=1.2);

}


###############################################################################

full_coef_table=full_coxph_fit_summ$coefficients;

signif_coefficients=(full_coef_table[,"Pr(>|z|)"]<=PvalCutoff)
signif_categories=full_coef_table[signif_coefficients,,drop=F];
signif_cat_coef_table=signif_categories[intersect(rownames(signif_categories), alr_cat_names),,drop=F];
print(signif_cat_coef_table);

signif_cat=rownames(signif_cat_coef_table);

par(mfrow=c(1,1));
plot_text(c(
	paste("Significant ALR Categories (alpha = ", PvalCutoff, "):", sep=""),
	capture.output(print(signif_cat_coef_table)),
	"",
	"Subsequent time series plots only illustrate significant categories."
));

par(mfrow=c(plots_per_page,1));

for(coht_ix in uniq_coht_ids){
	par(mfrow=c(plots_per_page,1));

	cat("Working on Cohort: ", coht_ix, "\n");
	subj_in_coht=names(group_structure[[coht_ix]])
	plot_ix=1;

	for(subj_ix in subj_in_coht){

		cat("  Subject: ", subj_ix, "\n");
		samp_ids=group_structure[[coht_ix]][[subj_ix]][["sample_id"]];
		times=group_structure[[coht_ix]][[subj_ix]][["times"]];
		death=group_structure[[coht_ix]][[subj_ix]][["death"]];
		plot_alr_over_time(samp_ids, times, death, 
			alr_categories_val[samp_ids, ,drop=F], 
			c(0, last_measured_time), alr_range, colors,
			title=subj_ix, keep=signif_cat);

		if(plot_ix==1){
			mtext(coht_ix, side=3, line=0, outer=T, font=2, cex=2);
		}

		if(plot_ix==plots_per_page-1 || subj_ix==tail(subj_in_coht,1)){
			plot_category_colors(alr_categories_val, colors, keep=signif_cat);
			plot_ix=1;
		}else{
			plot_ix=plot_ix+1;
		}
	}
}

###############################################################################

print(last_measured_time);
print(death_times);
death_time_nona=death_times[!is.na(death_times)];
print(death_time_nona);
num_death_time_nona=length(death_time_nona);
death_time_95pi=quantile(death_time_nona, c(.025, .5, .975));

jitter=rnorm(num_death_time_nona, 0, .1);

par(mar=c(4,4,4,4));
par(mfrow=c(2,2));
plot(0,0, type="n", 
	xlim=c(-.5,.5), ylim=c(0, last_measured_time), 
	ylab="Event Times", xlab="", xaxt="n", main="Event Prediction Intervals");
abline(h=death_time_95pi[c(1,3)], col="blue", lty=2);
abline(h=death_time_95pi[2], col="blue", lwd=2);
axis(4, death_time_95pi, labels=c("95% LB", "Median", "95% UB"));
points(jitter, death_time_nona);

cat("Prediction Intervals around Death Times:\n");
print(death_time_95pi);

if(EpochFile!=""){
	epochs=load_epochs(EpochFile);
}else{
	epochs=list();
	epochs[["Start"]]=c(0,0);
	epochs[["Before"]]=c(0, death_time_95pi[1]);
	epochs[["During"]]=c(death_time_95pi[1], death_time_95pi[3]);
	epochs[["After"]]=c(death_time_95pi[3], last_measured_time);
	epochs[["End"]]=c(last_measured_time, last_measured_time);
}

print(epochs);
epoch_names=names(epochs);
num_epochs=length(epoch_names);

###############################################################################

plot(0,0, type="n", 
	xlim=c(-.5,.5+num_epochs+2), ylim=c(0, last_measured_time), 
	ylab="Event Times", xlab="", xaxt="n", main="Defined Epochs");

for(i in 1:num_epochs){
	ep_nm=epoch_names[i];
	points(
		c(i,i),
		c(epochs[[ep_nm]][1], epochs[[ep_nm]][2]), 
		lwd=5, type="l",
		col="darkgreen");
}

for(i in 1:num_epochs){
	ep_nm=epoch_names[i];
	text(i, mean(epochs[[ep_nm]]), labels=ep_nm, font=2, pos=4 );
}

points(jitter, death_time_nona);

###############################################################################

hist(death_time_nona, xlim=c(0, last_measured_time), xlab="Event Times",
	main="Event Times");

###############################################################################

# Plot number of samples over time

hist_res=hist(factors[,TimeVarName], plot=F, breaks=seq(-.5, last_measured_time+.5, 1));
print(hist_res);

max_count=max(hist_res$counts);
plot(0,0, type="n", xlim=c(0, max_count*2), ylim=c(0, last_measured_time),
	yaxt="n", xlab="Counts", ylab="Time", main="Number of Samples Over Time");

num_bins=length(hist_res$counts);
for(i in 1:num_bins){
	rect(0, hist_res$breaks[i], hist_res$counts[i], hist_res$breaks[i+1], col="grey");
}

num_labels=10;
axis_ticks=round(seq(0, last_measured_time, length.out=num_labels), 0);
axis(side=2, 
	at=axis_ticks, 
	labels=axis_ticks,
	las=2);

for(i in 1:num_epochs){
	ep_nm=epoch_names[i];
	points(
		c(i*2+max_count+3,i*2+max_count+3),
		c(epochs[[ep_nm]][1], epochs[[ep_nm]][2]), 
		lwd=5, type="l", 
		col="darkgreen");
}

for(i in 1:num_epochs){
	ep_nm=epoch_names[i];
	text(i*2+max_count+3, mean(epochs[[ep_nm]]), labels=ep_nm, 
		font=2, pos=4, cex=.7 );
}

###############################################################################

avg_by_subj=function(comb_fact, subj_id_colname, coht_id_colname, alr_id_colname){

	uniq_subj_ids=unique(comb_fact[,subj_id_colname]);
	out_mat=as.data.frame(matrix(NA, nrow=length(uniq_subj_ids), ncol=2));
	colnames(out_mat)=c(coht_id_colname, alr_id_colname);
	rownames(out_mat)=uniq_subj_ids;

	for(sbj in uniq_subj_ids){
		subj_rec_ix=comb_fact[,subj_id_colname]==sbj;
		subj_rec=comb_fact[subj_rec_ix,];
		mean=mean(subj_rec[,alr_id_colname]);
		out_mat[sbj, coht_id_colname]=as.character(subj_rec[1,coht_id_colname]);
		out_mat[sbj, alr_id_colname]=mean;
	}

	return(out_mat);
}

###############################################################################

plot_alr_diff=function(alrA, alrB, nameA, nameB, title, y_range){

	mean_alra=mean(alrA[,1]);
	mean_alrb=mean(alrB[,1]);

	tt_res=t.test(alrA[,1], alrB[,1]);

	plot(0,0, type="n", xlim=c(0,1), ylim=y_range, main="",
		xlab="", ylab="ALR", xaxt="n");

	if(tt_res$p.value<.001){
		sig_char="***";
	}else if(tt_res$p.value<.01){
		sig_char="**";
	}else if(tt_res$p.value<.05){
		sig_char="*";
	}else{
		sig_char="";
	}

	if(nchar(title)>18){
		title_size=1.2*18/nchar(title);
	}else{
		title_size=1.2;
	}
	mtext(text=paste(title, sig_char, sep=""), line=1.8, cex=title_size, font=2);

	mtext(text=paste(nameA, ": u=", round(mean_alra,4), sep=""), line=1.2, cex=.6);
	mtext(text=paste(nameB, ": u=", round(mean_alrb,4), sep=""), line=.6, cex=.6);
	mtext(text=paste("T-Test p-value: ", round(tt_res$p.value, 4), sep=""), line=0, cex=.6, font=3);

	axis(1, at=c(.25, .75), labels=c(nameA, nameB));

	alrA_len=nrow(alrA);
	alrB_len=nrow(alrB);

	points(c(.15,.35), c(mean_alra, mean_alra), type="l", col="blue");
	points(c(.65,.85), c(mean_alrb, mean_alrb), type="l", col="blue");

	points(rep(.25, alrA_len), alrA[,1]);
	points(rep(.75, alrB_len), alrB[,1]);

}

num_uniq_cohts=length(uniq_coht_ids);
num_epochs=length(epoch_names);

samp_tpt=factors[,TimeVarName];
samp_cht=factors[,CohtVarName];
samp_names=rownames(factors);

par(mfrow=c(1,1));
par(mar=c(0,0,0,0));
plot(0,0, type="n", xlab="", ylab="", xaxt="n", yaxt="n", main="", bty="n");
text(0,0, paste("Comparisons of ", CohtVarName,"\nby Epoch",sep=""), font=2, cex=2);

# Compare the cohorts by epochs
par(mar=c(1,3,5,1));
par(oma=c(1,1,3,1));

for(cat_ix in signif_cat){

	alr_range=range(alr_categories_val[, cat_ix]);

	for(ep_nm in epoch_names){

		par(mfrow=c(3,3));
		cat("Look at epoch: ", ep_nm, "\n");

		ep_range=epochs[[ep_nm]];
		ep_ix=(samp_tpt>=ep_range[1] & samp_tpt<=ep_range[2]);
		
		for(chtA_ix in 1:num_uniq_cohts){
			chtA=uniq_coht_ids[chtA_ix];
			a_pts_ix=(samp_cht==chtA) & ep_ix;
			fact_A=factors[a_pts_ix,];
			cur_alr_cat=alr_categories_val[rownames(fact_A), cat_ix];
			combined_A=cbind(fact_A, cur_alr_cat);
			collapsed_A=avg_by_subj(combined_A, SubjVarName, CohtVarName, "cur_alr_cat");

			for(chtB_ix in 1:num_uniq_cohts){

				if(chtA_ix<chtB_ix){

					chtB=uniq_coht_ids[chtB_ix];
					b_pts_ix=(samp_cht==chtB) & ep_ix;
					fact_B=factors[b_pts_ix,];
					cur_alr_cat=alr_categories_val[rownames(fact_B), cat_ix];
					combined_B=cbind(fact_B, cur_alr_cat);
					collapsed_B=avg_by_subj(combined_B, SubjVarName, CohtVarName, "cur_alr_cat");

					cat("Comparing: ", chtA, " vs ", 
						chtB, "\n");

					plot_alr_diff(
						#alr_categories_val[a_pts_ix, cat_ix, drop=F], 
						#alr_categories_val[b_pts_ix, cat_ix, drop=F],
						collapsed_A[,"cur_alr_cat",drop=F],
						collapsed_B[,"cur_alr_cat",drop=F],
						nameA=as.character(chtA),
						nameB=as.character(chtB),
						title=cat_ix,
						y_range=alr_range
					);
				}
			}
		}

		mtext(paste(ep_nm, ": [", epochs[[ep_nm]][1], ", ", 
			epochs[[ep_nm]][2], "]", sep="") , side=3, outer=T, font=2);

	}
}


###############################################################################

# Compare epochs by cohorts
#print(factors)

plot_epoch_comp=function(alr_a_table, alr_b_table, nameA, nameB, 
	alr_colname, cht_colname, cht_colors, alr_ranges, mtitle){

	laymat=matrix(c(
		1,1,1,1,
		2,2,3,3,
		4,4,5,5), ncol=4, byrow=T);
	layout(laymat);

	orig_par=par(no.readonly=T);

	subjects=sort(unique(rownames(alr_a_table), rownames(alr_b_table)));
	num_subjects=length(subjects);

	cat("Plotting: ", nameA, " vs. ", nameB, "\n", sep="");

	shrd_sbj=intersect(rownames(alr_a_table), rownames(alr_b_table));
	num_shrd_sbj=length(shrd_sbj);
	shrd_a=alr_a_table[shrd_sbj,,drop=F];
	shrd_b=alr_b_table[shrd_sbj,,drop=F];
	num_censored=num_subjects-num_shrd_sbj;
	alr_range=range(alr_a_table[,alr_colname], alr_b_table[,alr_colname]);

	excl_a_sbj=setdiff(rownames(alr_a_table),shrd_sbj);

	# Compute p-values within each cohort
	chts=sort(unique(alr_a_table[,cht_colname]));
	num_chts=length(chts);
	legend_info=character(num_chts);
	names(legend_info)=chts;
	for(i in 1:num_chts){
		cc=chts[i];
		c_ix=(shrd_a[,cht_colname]==cc);
		tres=t.test(shrd_a[c_ix,alr_colname], shrd_b[c_ix,alr_colname], paired=T);
		pval=tres$p.value;

		sigchar="";
		if(!is.na(pval)){
			if(pval<.001){
				sigchar="***";
			}else if(pval<.01){
				sigchar="**";
			}else if(pval<.05){
				sigchar="*";
			}
		}

		legend_info[cc]=paste(cc, sigchar, ", p=", round(pval,3), sep="");
	}

	# Setup plot area
	alr_diff=diff(alr_ranges);
	pad=alr_diff*.05;
	x_plot_rang=c(alr_ranges[1]-pad, alr_ranges[2]+pad);
	y_plot_rang=c(alr_ranges[1]-pad, alr_ranges[2]+4*pad);


	# Plot each subject
	a_only=c();
	b_only=c();
	a_compl=c();
	b_compl=c();
	
	# lines through centroids
	for(i in 1:num_subjects){
		subj_id=subjects[i];

		a_alr=alr_a_table[subj_id, alr_colname];
		b_alr=alr_b_table[subj_id, alr_colname];

		pt_col=cht_colors[alr_a_table[subj_id, cht_colname]];

		if(!is.na(a_alr) && is.na(b_alr)){
			#abline(v=a_alr, col=pt_col, lwd=.5);
			a_only=c(a_only, a_alr);
		}else if(is.na(a_alr) && !is.na(b_alr)){
			#abline(h=b_alr, col=pt_col, lwd=.5);
			b_only=c(b_only, b_alr);
		}else{
			a_compl=c(a_compl, a_alr);
			b_compl=c(b_compl, b_alr);
		}

	}

	##############################################################################
	# Generate pairwise scatter plots
	plot(0,0, type="n", xlim=x_plot_rang, ylim=y_plot_rang, xlab="", ylab="");

	title(xlab=nameA, ylab=nameB, line=2, font.lab=2);
	abline(a=0, b=1, col="grey", lty=1, lwd=2);

	legend(x_plot_rang[1], y_plot_rang[2], 
		legend=legend_info[chts],
		fill=cht_colors[chts], bty="n");

	abline(v=mean(a_only), lty=3, col="gray50");
	abline(h=mean(b_only), lty=3, col="gray50");
	abline(v=mean(a_compl), lty=2, lwd=1,col="grey40");
	abline(h=mean(b_compl), lty=2, lwd=1, col="grey40");

	
	# Draw Points
	for(i in 1:num_subjects){
		subj_id=subjects[i];

		a_alr=alr_a_table[subj_id, alr_colname];
		b_alr=alr_b_table[subj_id, alr_colname];

		pt_col=cht_colors[alr_a_table[subj_id, cht_colname]];

		if(!is.na(a_alr) && is.na(b_alr)){
			#abline(v=a_alr, col=pt_col, lwd=.5);
			points(a_alr, alr_ranges[1]-pad, pch=4, col=pt_col);
			a_only=c(a_only, a_alr);
		}else if(is.na(a_alr) && !is.na(b_alr)){
			#abline(h=b_alr, col=pt_col, lwd=.5);
			points(alr_ranges[1]-pad, b_alr, pch=4, col=pt_col);
			b_only=c(b_only, b_alr);
		}else{
			points(a_alr, b_alr, pch=19, col=pt_col);
			points(a_alr, b_alr, cex=.2, col="black");
			a_compl=c(a_compl, a_alr);
			b_compl=c(b_compl, b_alr);
		}

	}

	##############################################################################
	# Generate difference bar plots

	ci95=function(values){
		stdev=sd(values);
		mean=mean(values);

		num_smp=length(values);
		intvl=stdev/sqrt(num_smp)*1.96;
		ub=mean+intvl;
		lb=mean-intvl;
		return(c(lb,ub));
	}

	shrd_diff=as.data.frame(matrix(NA, nrow=num_shrd_sbj, ncol=2));
	alr_diff_colname=paste(alr_colname, "_diff", sep="");
	colnames(shrd_diff)=c(cht_colname, alr_diff_colname);
	rownames(shrd_diff)=rownames(shrd_a);
	shrd_diff[,cht_colname]=shrd_a[,cht_colname];
	shrd_diff[,alr_diff_colname]=shrd_b[,alr_colname]-shrd_a[,alr_colname];
	num_shrd_diff=length(shrd_diff[,alr_diff_colname]);
	cat("Number of Shared Differences: ", num_shrd_diff, "\n");
	print(shrd_diff[,alr_diff_colname]);

	cht_shrd=shrd_a[,cht_colname];

	mean_comb_diff=mean(shrd_diff[,alr_diff_colname]);
	ci95_comb_diff=ci95(shrd_diff[,alr_diff_colname]);
	max_mag_diff=max(abs(shrd_diff[,alr_diff_colname]));
	cat("Maximum difference: ", max_mag_diff, "\n");
	mean_comb_barwidth=2;
	spacer_width=1;

	orig_mar=orig_par$mar;
	cur_mar=orig_mar;
	cur_mar[1]=orig_mar[1]+5;
	par(mar=cur_mar);
	plot(0,0, type="n", 
		xlim=c(0, num_chts+spacer_width+mean_comb_barwidth+1),
		ylim=c(-max_mag_diff, max_mag_diff),
		xaxt="n",
		xlab="", ylab="", main="");
	title(ylab=paste("Difference: ", nameB, " - ", nameA), line=2, font.lab=2);
	abline(h=0, col="grey", lwd=2);

	bar_pos=1;
	mean_line_width=.75;
	for(ch_id in chts){
		ch_ix=(cht_shrd==ch_id);
		ch_diff=shrd_diff[ch_ix,alr_diff_colname,drop=F];
		num_ch_smp=nrow(ch_diff);
		mean_diff=mean(ch_diff[,alr_diff_colname]);

		pt_col=cht_colors[ch_id];

		# mean line
		points(c(bar_pos-mean_line_width/2, bar_pos+mean_line_width/2), rep(mean_diff, 2),
			type="l", lwd=1.5, col=pt_col, lend="square");
		points(c(bar_pos-mean_line_width/2, bar_pos-mean_line_width/2), c(mean_diff, 0),
			type="l", lwd=1.2, col=pt_col, lend="square");
		points(c(bar_pos+mean_line_width/2, bar_pos+mean_line_width/2), c(mean_diff, 0),
			type="l", lwd=1.2, col=pt_col, lend="square");

		points(rep(bar_pos, num_ch_smp), ch_diff[,alr_diff_colname], pch=19, col=pt_col);
		points(rep(bar_pos, num_ch_smp), ch_diff[,alr_diff_colname], cex=.2, col="black");

		axis(side=1, at=bar_pos, labels=ch_id, las=2, font.axis=2);
		bar_pos=bar_pos+1;
	}

	# overall mean line
	bar_pos=num_chts+spacer_width+1;
	jitter=rnorm(num_shrd_diff, 0, .1*mean_comb_barwidth);
	points(c(bar_pos-mean_comb_barwidth/2, bar_pos+mean_comb_barwidth/2), rep(mean_comb_diff, 2),
		type="l", lwd=1.5, col="black", lty=2, lend="square");
	points(c(bar_pos-mean_comb_barwidth/2, bar_pos-mean_comb_barwidth/2), c(mean_comb_diff, 0),
		type="l", lwd=1.2, col="black", lty=2, lend="square");
	points(c(bar_pos+mean_comb_barwidth/2, bar_pos+mean_comb_barwidth/2), c(mean_comb_diff, 0),
		type="l", lwd=1.2, col="black", lty=2, lend="square");
	points(bar_pos+jitter, 
		shrd_diff[,alr_diff_colname], pch=1, col="grey40");
	axis(side=1, at=bar_pos, labels="Combined", las=2, font.axis=2);
	mtext(text=paste("Num Censored Meas.: ", num_censored, sep=""), side=3, line=0, cex=.5, font=3);

	##############################################################################
	# Generate bar plot without CI bars only

	orig_mar=orig_par$mar;
	cur_mar=orig_mar;
	cur_mar[1]=orig_mar[1]+5;
	par(mar=cur_mar);
	plot(0,0, type="n", 
		xlim=c(0, num_chts+spacer_width+mean_comb_barwidth+1),
		ylim=c(-max_mag_diff, max_mag_diff),
		xaxt="n",
		xlab="", ylab="", main="");
	title(ylab=paste("Difference: ", nameB, " - ", nameA), line=2, font.lab=2);
	abline(h=0, col="grey", lwd=2);

	bar_pos=1;
	mean_line_width=.75;
	for(ch_id in chts){
		ch_ix=(cht_shrd==ch_id);
		ch_diff=shrd_diff[ch_ix,alr_diff_colname,drop=F];
		num_ch_smp=nrow(ch_diff);

		mean_diff=mean(ch_diff[,alr_diff_colname]);
		ch_diff_ci95=ci95(ch_diff[,alr_diff_colname]);

		pt_col=cht_colors[ch_id];

		# mean line
		points(c(bar_pos-mean_line_width/2, bar_pos+mean_line_width/2), rep(mean_diff, 2),
			type="l", lwd=1.5, col=pt_col, lend="square");
		points(c(bar_pos-mean_line_width/2, bar_pos-mean_line_width/2), c(mean_diff, 0),
			type="l", lwd=1.2, col=pt_col, lend="square");
		points(c(bar_pos+mean_line_width/2, bar_pos+mean_line_width/2), c(mean_diff, 0),
			type="l", lwd=1.2, col=pt_col, lend="square");

		# 95% CI bar
		points(c(bar_pos-mean_line_width/4, bar_pos+mean_line_width/4), rep(ch_diff_ci95[1], 2),
			type="l", lwd=1, col="grey25", lend="square");
		points(c(bar_pos-mean_line_width/4, bar_pos+mean_line_width/4), rep(ch_diff_ci95[2], 2),
			type="l", lwd=1, col="grey25", lend="square");
		points(c(bar_pos, bar_pos), ch_diff_ci95,
			type="l", lwd=1, col="grey25", lend="square");
		

		axis(side=1, at=bar_pos, labels=ch_id, las=2, font.axis=2);
		bar_pos=bar_pos+1;
	}

	# overall mean line
	bar_pos=num_chts+spacer_width+1;

	# Plot bars
	points(c(bar_pos-mean_comb_barwidth/2, bar_pos+mean_comb_barwidth/2), rep(mean_comb_diff, 2),
		type="l", lwd=1.5, col="black", lty=2, lend="square");
	points(c(bar_pos-mean_comb_barwidth/2, bar_pos-mean_comb_barwidth/2), c(mean_comb_diff, 0),
		type="l", lwd=1.2, col="black", lty=2, lend="square");
	points(c(bar_pos+mean_comb_barwidth/2, bar_pos+mean_comb_barwidth/2), c(mean_comb_diff, 0),
		type="l", lwd=1.2, col="black", lty=2, lend="square");

	# Annotate 95% CI
	points(c(bar_pos-mean_line_width/4, bar_pos+mean_line_width/4), rep(ci95_comb_diff[1], 2),
		type="l", lwd=1, col="grey25", lend="square");
	points(c(bar_pos-mean_line_width/4, bar_pos+mean_line_width/4), rep(ci95_comb_diff[2], 2),
		type="l", lwd=1, col="grey25", lend="square");
	points(c(bar_pos, bar_pos), ci95_comb_diff,
		type="l", lwd=1, col="grey25", lend="square");

	axis(side=1, at=bar_pos, labels="Combined", las=2, font.axis=2);
	mtext(text=paste("Num Censored Meas.: ", num_censored, sep=""), side=3, line=0, cex=.5, font=3);
	mtext(text=paste("Mean Difference\nw/ 95% CI around Means"), side=3, line=-2.5, cex=.75, font=2);

	# Reset original margins
	par(mar=orig_mar);


	##############################################################################
	# Generate connected lines plot
	
	# Plot positions and widths
	censor_colwid=1;
	a_colwid=1;
	b_colwid=1;
	xpos=.75;
	apos=2;
	bpos=4;

	# Plot specific margins
	cur_mar=orig_mar;
	cur_mar[2]=2;
	cur_mar[1]=2.5;
	par(mar=cur_mar);

	# Extract out samples in A, censored in B
	excl_a_sbj=setdiff(rownames(alr_a_table),shrd_sbj);
	plot_excl_a=F;
	if(length(excl_a_sbj)){
		excl_a=alr_a_table[excl_a_sbj,];
		mean_excl_a=mean(excl_a[,alr_colname]);
		plot_excl_a=T;
	}

	# Setup plot
	plot_wid=censor_colwid+a_colwid+b_colwid+.5;
	plot(0,0, type="n",
		xlim=c(0, plot_wid+1),
		ylim=alr_ranges,
		xaxt="n",
		xlab="", ylab=""
		);

	# Label B column
	axis(side=1, at=c(bpos), labels=c(nameB), 
		las=1, font.axis=2, tick=F, line=0);
	# Highlight A and B column
	abline(v=c(apos, bpos), lwd=30, col="grey90", lend=2);

	a_shrd_mean=mean(shrd_a[,alr_colname]);
	b_shrd_mean=mean(shrd_b[,alr_colname]);
	abline(h=c(a_shrd_mean, b_shrd_mean), lwd=.5, col="grey90");

	# Label A column depending on whether there any censored samples
	if(plot_excl_a){
		axis(side=1, at=xpos, labels=paste("died before\n",nameB,sep=""), 
			las=1, font.axis=3, cex.axis=.60, line=0);
		axis(side=1, at=apos, labels=paste("survived to\n",nameB,sep=""), 
			las=1, font.axis=3, cex.axis=.60, line=0);
		abline(v=c(xpos), lwd=20, col="grey90", lend=1);
		axis(side=1, at=c((xpos+apos)/2), labels=c(nameA), 
			las=1, font.axis=2, tick=F, line=1.1);

		abline(h=c(mean_excl_a), lwd=.5, col="grey90");

	}else{
		axis(side=1, at=c(apos), labels=c(nameA), 
			las=1, font.axis=2, tick=F, line=0);
	}

	# Draw line through means
	points(c(apos-a_colwid/2, apos+a_colwid/2), c(a_shrd_mean, a_shrd_mean), 
		lty=1, lwd=2, type="l");
	points(c(bpos-b_colwid/2, bpos+b_colwid/2), c(b_shrd_mean, b_shrd_mean), 
		lty=1, lwd=2, type="l");

	# Label mean of censored
	if(plot_excl_a){
		points(c(xpos-censor_colwid/2, xpos+censor_colwid/2), c(mean_excl_a, mean_excl_a), 
			lty=1, lwd=2, type="l");
	}
	
	# Draw points by color	
	bar_pos=1;
	for(ch_id in chts){

		ch_ix=(cht_shrd==ch_id);
		a_pts=shrd_a[ch_ix,alr_colname];
		b_pts=shrd_b[ch_ix,alr_colname];
		pt_col=cht_colors[ch_id];
		num_ch_smp=length(a_pts);

		# Draw lines between A and B points 
		jitter=rnorm(num_ch_smp, 0, .10);
		for(i in 1:num_ch_smp){
			points(
				c(apos+jitter[i], bpos+jitter[i]), 
				c(a_pts[i], b_pts[i]),
				type="l", lwd=.2, col=pt_col
			); 
		}

		# Draw points for A and B
		points(apos+jitter, a_pts, pch=19, col=pt_col);
		points(apos+jitter, a_pts, cex=.1, col="black");

		points(bpos+jitter, b_pts, pch=19, col=pt_col);
		points(bpos+jitter, b_pts, cex=.1, col="black");

		# Draw censored points if there are any
		if(plot_excl_a){
			ch_ix=(ch_id==excl_a[,cht_colname]);
			excl=excl_a[ch_ix, alr_colname];
			points(rep(xpos, length(excl)), excl, col=pt_col, pch=4);
		}


	}

	##############################################################################
	# Plot means and CI only
	
	# Setup plot
	plot_wid=censor_colwid+a_colwid+b_colwid+.5;
	plot(0,0, type="n",
		xlim=c(0, plot_wid+1),
		ylim=alr_ranges,
		xaxt="n",
		xlab="", ylab=""
		);

	# Label B column
	axis(side=1, at=c(bpos), labels=c(nameB), 
		las=1, font.axis=2, tick=F, line=0);
	# Highlight A and B column
	abline(v=(apos+bpos)/2, lwd=2, col="grey50", lend=2);

	a_shrd_mean=mean(shrd_a[,alr_colname]);
	b_shrd_mean=mean(shrd_b[,alr_colname]);
	#abline(h=c(a_shrd_mean, b_shrd_mean), lwd=.5, col="grey90");

	# Label A column depending on whether there any censored samples
	if(plot_excl_a){
		axis(side=1, at=xpos, labels=paste("died before\n",nameB,sep=""), 
			las=1, font.axis=3, cex.axis=.60, line=0);
		axis(side=1, at=apos, labels=paste("survived to\n",nameB,sep=""),
			las=1, font.axis=3, cex.axis=.60, line=0);
		abline(v=(xpos+apos)/2, lwd=1, col="grey60", lend=1);
		axis(side=1, at=c((xpos+apos)/2), labels=c(nameA), 
			las=1, font.axis=2, tick=F, line=1.1);

		#abline(h=c(mean_excl_a), lwd=.5, col="grey90");

	}else{
		axis(side=1, at=c(apos), labels=c(nameA), 
			las=1, font.axis=2, tick=F, line=0);
	}

	# Draw line through means
	#points(c(apos-a_colwid/2, apos+a_colwid/2), c(a_shrd_mean, a_shrd_mean), 
	#	lty=1, lwd=2, type="l");
	#points(c(bpos-b_colwid/2, bpos+b_colwid/2), c(b_shrd_mean, b_shrd_mean), 
	#	lty=1, lwd=2, type="l");

	# Label mean of censored
	#if(plot_excl_a){
	#	points(c(xpos-censor_colwid/2, xpos+censor_colwid/2), c(mean_excl_a, mean_excl_a), 
	#		lty=1, lwd=2, type="l");
	#}
	
	colwid=a_colwid;
	cht_pos=seq(0, colwid, length.out=num_chts+4)[2:(num_chts+3)];
	

	# Draw points by color	
	cht_ix=1;
	for(ch_id in chts){

		ch_ix=(cht_shrd==ch_id);
		a_pts=shrd_a[ch_ix,alr_colname];
		b_pts=shrd_b[ch_ix,alr_colname];
		pt_col=cht_colors[ch_id];
		num_ch_smp=length(a_pts);
	
		ci95a_pts=ci95(a_pts);
		ci95b_pts=ci95(b_pts);

		#points(rep(apos-colwid/2+cht_pos[cht_ix],2), ci95a_pts, col=pt_col, type="l", lwd=2);
		#points(rep(bpos-colwid/2+cht_pos[cht_ix],2), ci95b_pts, col=pt_col, type="l", lwd=2);

		if(length(a_pts)>=2){
			ci95a_pts=ci95(a_pts);
			points(rep(apos-colwid/2+cht_pos[cht_ix],2), ci95a_pts, col=pt_col, type="l", lwd=2);
			points(apos-colwid/2+cht_pos[cht_ix], mean(ci95a_pts), col=pt_col, type="p");
			
		}else if(length(a_pts)==1){
			points(apos-colwid/2+cht_pos[cht_ix], a_pts, col=pt_col, type="p");
		}

		if(length(b_pts)>=2){
			ci95b_pts=ci95(b_pts);
			points(rep(bpos-colwid/2+cht_pos[cht_ix],2), ci95b_pts, col=pt_col, type="l", lwd=2);
			points(bpos-colwid/2+cht_pos[cht_ix], mean(ci95b_pts), col=pt_col, type="p");
		}else if(length(b_pts)==1){
			points(bpos-colwid/2+cht_pos[cht_ix], b_pts, col=pt_col, type="p");
		}

		# Draw censored points if there are any
		if(plot_excl_a){
			ch_ix=(ch_id==excl_a[,cht_colname]);
			excl=excl_a[ch_ix, alr_colname];

			if(length(excl)>=2){
				ci95x_pts=ci95(excl);
				points(rep(xpos-colwid/2+cht_pos[cht_ix],2), ci95x_pts, col=pt_col, type="l", lwd=2);
				points(xpos-colwid/2+cht_pos[cht_ix], mean(ci95x_pts), col=pt_col, type="p");
			}else if(length(excl)==1){
				points(xpos-colwid/2+cht_pos[cht_ix], excl, col=pt_col, type="p");
			}
		}

		cht_ix=cht_ix+1;

	}

	cht_ix=num_chts+2;
	a_pts=shrd_a[,alr_colname];

	pt_col="black";
	if(length(a_pts)>=2){
		ci95a_pts=ci95(a_pts);
		points(rep(apos-colwid/2+cht_pos[cht_ix],2), ci95a_pts, col=pt_col, type="l", lwd=3);
		points(apos-colwid/2+cht_pos[cht_ix], mean(ci95a_pts), col=pt_col, type="p");
		
	}else if(length(a_pts)==1){
		points(apos-colwid/2+cht_pos[cht_ix], a_pts, col=pt_col, type="p");
	}

	b_pts=shrd_b[,alr_colname];
	if(length(b_pts)>=2){
		ci95b_pts=ci95(b_pts);
		points(rep(bpos-colwid/2+cht_pos[cht_ix],2), ci95b_pts, col=pt_col, type="l", lwd=3);
		points(bpos-colwid/2+cht_pos[cht_ix], mean(ci95b_pts), col=pt_col, type="p");
	}else if(length(b_pts)==1){
		points(bpos-colwid/2+cht_pos[cht_ix], b_pts, col=pt_col, type="p");
	}

	# Draw censored points if there are any
	if(plot_excl_a){
		excl=excl_a[, alr_colname];

		if(length(excl)>=2){
			ci95x_pts=ci95(excl);
			points(rep(xpos-colwid/2+cht_pos[cht_ix],2), ci95x_pts, col=pt_col, type="l", lwd=3);
			points(xpos-colwid/2+cht_pos[cht_ix], mean(ci95x_pts), col=pt_col, type="p");
		}else if(length(excl)==1){
			points(xpos-colwid/2+cht_pos[cht_ix], excl, col=pt_col, type="p");
		}
	}

	mtext(text="Mean ALR Abundance\nw/ 95% CI's", side=3, line=0, cex=.9, font=2);

	##############################################################################

	#plot(0,0, type="n", xlab="", ylab="", main="", xaxt="n", yaxt="n", bty="n");
	
	par(mar=orig_mar);
	mtext(mtitle, side=3, outer=T, font=2);

}

###############################################################################

par(mfrow=c(1,1));
par(mar=c(0,0,0,0));
plot(0,0, type="n", xlab="", ylab="", xaxt="n", yaxt="n", main="", bty="n");
text(0,0, paste("Comparisons of Epochs\nby ", CohtVarName, "\nMatched Paired by ", SubjVarName, sep=""), 
	font=2, cex=2);

cht_colors=rainbow(num_uniq_cohts, start=0, end=2/3);
names(cht_colors)=uniq_coht_ids;

for(cat_ix in signif_cat){
	cat("\nWorking on: ", cat_ix, "\n");	

	par(mfrow=c(3,3));
	par(mar=c(4,4,3,1));

	alr_ranges=range(alr_categories_val[,cat_ix]);

	for(ep_ix_A in 1:(num_epochs-1)){


		ep_nm_A=epoch_names[ep_ix_A];
		ep_range=epochs[[ep_nm_A]];
                ep_ix=(samp_tpt>=ep_range[1] & samp_tpt<=ep_range[2]);
		fact_A=factors[ep_ix,];

		cur_alr_cat=alr_categories_val[rownames(fact_A), cat_ix];
		combined_A=cbind(fact_A, cur_alr_cat);
		#print(combined_A);

		collapsed_A=avg_by_subj(combined_A, SubjVarName, CohtVarName, "cur_alr_cat");
	
		for(ep_ix_B in (ep_ix_A+1):num_epochs){

			ep_nm_B=epoch_names[ep_ix_B];

			cat("Comparing: ", ep_nm_A, " and ", ep_nm_B, "\n", sep="");

			ep_range=epochs[[ep_nm_B]];
			ep_ix=(samp_tpt>=ep_range[1] & samp_tpt<=ep_range[2]);
			fact_B=factors[ep_ix,];

			cur_alr_cat=alr_categories_val[rownames(fact_B), cat_ix];
			combined_B=cbind(fact_B, cur_alr_cat);
			#print(combined_B);

			collapsed_B=avg_by_subj(combined_B, SubjVarName, CohtVarName, "cur_alr_cat");

			plot_epoch_comp(collapsed_A, collapsed_B, ep_nm_A, ep_nm_B, 
				"cur_alr_cat", CohtVarName, cht_colors, alr_ranges,
				mtitle=cat_ix);
		}
	}



}



if(0){

#pred_res=predict(coxph_fit, type="response", se.fit=T);
#print(pred_res);
lp_red=predict(coxph_fit,type="lp")
exp_pred=predict(coxph_fit,type="expected")
risk_pred=predict(coxph_fit,type="risk",se.fit=TRUE)
terms_pred=predict(coxph_fit,type="terms",se.fit=TRUE)


print(coxph_fit$y);
print(coxph_fit$linear.predictors);

print(summary(lp_red));

print(summary(exp_pred));

print(summary(risk_pred));

print(summary(terms_pred));
}


###############################################################################


cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
