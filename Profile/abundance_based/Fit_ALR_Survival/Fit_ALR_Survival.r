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

	"output_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

NUM_TOP_PRED_CAT=20;

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
	"\n",
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
pdf(paste(OutputRoot, ".surival.pdf", sep=""), height=14, width=8.5);

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
if(ModelVarFile!=""){
	model_var_arr=c(model_var_arr, load_list(ModelVarFile));
	cat("Model Variables:\n");
	print(model_var_arr);
	cat("\n");
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

plot_alr_over_time=function(ids, times, end_time, alrs, time_range, alr_range, colors, title=""){

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

	plot(0,0, type="n", xlim=xrang, ylim=yrang, ylab=title); 
	
	num_alrs=ncol(alrs);
	for(cat_ix in 1:num_alrs){
		points(times, alrs[ids, cat_ix], type="l", col=colors[cat_ix], lwd=max_m-(cat_ix/num_alrs)*mult);
		#points(times, alrs[ids, cat_ix], type="l", col="black", lwd=.05);
	}
	for(cat_ix in 1:num_alrs){
		points(times, alrs[ids, cat_ix], type="p", pch=16, cex=max_m-(cat_ix/num_alrs)*mult, col=colors[cat_ix]);
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

plot_category_colors=function(alrs, colors){
	
	cat("Plotting category colors...\n");
	names=colnames(alrs);
	num_cat=ncol(alrs);

	par(mar=c(2, 4, .5, 4));
	plot(0,0, type="n", xlim=c(0,10), ylim=c(0,1), ylab="", xlab="", bty="n", xaxt="n", yaxt="n");
	print(names);
	for(i in 1:num_cat){
		points(0, 1-(i/num_cat), col=colors[i], pch=16);
		text(0, 1-(i/num_cat), names[i], pos=4);
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

#print(rownames(alr_categories_val));
#print(alr_categories_val);


cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
