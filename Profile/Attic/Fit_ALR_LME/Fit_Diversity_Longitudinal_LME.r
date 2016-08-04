#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
library(nlme);
options(useFancyQuotes=F);

params=c(
	"summary_file", "s", 1, "character",
	"factors", "f", 1, "character",
	"groups", "r", 1, "character",
	"series", "e", 1, "character",
	"control", "c", 1, "character",
	"linear_model", "m", 1, "character",
	"trt_factor", "t", 1, "character",

	"outputroot", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

# Source the library
ALR_LIB_NAME="ALR_Library.r";

script_name_parts=strsplit(script_name, "/")[[1]];
if(length(script_name_parts)==1){
	script_name_parts=c(".",script_name_parts);
}
script_name_parts=script_name_parts[-length(script_name_parts)];
script_path=paste(script_name_parts, collapse="/");

lib_path=paste(script_path, "/", ALR_LIB_NAME, sep="");
source(lib_path);

# Describe usage
usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <Summary file table>\n",
	"	-f <Factors/metadata file>\n",
	"\n",
	"	-r <factor name for gRouping, e.g. patient ID>\n",
	"	-e <factor name for sEries, e.g. time>\n",
	"	-t <Treatment factor name>\n",
	"	-c <Control factor level, e.g. healthy>\n",
	"	-m \"<formula for linear Model, e.g. A + B + D*E>\"\n",
	"\n",
	"	[-o <Output filename root>]\n",
	"\n",
	"This script will fit a linear mixed effect model to the various\n",
	"diversity measures.  \n",
	"\n",
	"The random effect will be specified by the Grouping (-r) and\n",
	"the correlation structure will be based on AR1 using the information\n",
	"in the Series (-e).  The reference treatment level needs to be specified\n",
	"by indicating the name of the treatment (-t) and the reference or control level\n",
	"(-c).\n",
	"\n",
	"Specify the RHS (right hand side) of the model (-m) as a formula.  The response (LHS)\n",
	"will be dynamically specified by cycling through each of the diversity indices.\n",
	"\n",
	"\n");

# Prepare parameters based on options
if(!length(opt$summary_file) || !length(opt$factors) || 
   !length(opt$groups) || !length(opt$series) ||
   !length(opt$trt_factor) ||
   !length(opt$control) || !length(opt$linear_model)
){
	cat(usage);
	q(status=-1);
}else{
	SummaryFile=opt$summary_file;
	FactorsFile=opt$factors;
	GroupFactorName=opt$groups;
	SeriesFactorName=opt$series;
	LinearModel=opt$linear_model;
	TrtFactorName=opt$trt_factor;
	ControlFactorLevel=opt$control;
}

if(!length(opt$outputroot)){
	OutputRoot=gsub(".summary_table.xls$", "", opt$summary_file);
	OutputRoot=gsub(".summary_table.tsv$", "", OutputRoot);
}else{
	OutputRoot=opt$outputroot;
}

##############################################################################

cat("\n");
cat("Summary File: ", SummaryFile, "\n", sep="");
cat("Factors File: ", FactorsFile, "\n", sep="");
cat(" Output File: ", OutputRoot, "\n", sep="");
cat("\n");
cat(" Grouping Factor Name: ", GroupFactorName, "\n");
cat("   Series Factor Name: ", SeriesFactorName, "\n");
cat("Treatment Factor Name: ", TrtFactorName, "\n");
cat("\n");
cat("Linear Model: ", LinearModel, "\n");
cat("\n");

options(width=80);
cat("Text Line Width: ", options()$width, "\n\n", sep="");

##############################################################################

tail_statistic=function(x){
        sorted=sort(x, decreasing=TRUE);
        norm=sorted/sum(x);
        n=length(norm);
        tail=0;
        for(i in 1:n){
                tail=tail + norm[i]*((i-1)^2);
        }
        return(sqrt(tail));
}

compute_diversity=function(norm_cat_prof){
        zf_norm=norm_cat_prof[norm_cat_prof>0];

	div=numeric(5);
	names(div)=c("Tail", "Shannon", "Simpson", "Evenness", "SimpsonsRecip");

	div["Tail"]=tail_statistic(zf_norm);
        div["Shannon"]=-sum(zf_norm*log(zf_norm));
        div["Simpson"]=1-sum(zf_norm^2);
        div["Evenness"]=div["Shannon"]/log(length(zf_norm));
        div["SimpsonsRecip"]=1/sum(zf_norm^2);

	return(div);
}

##############################################################################
# Load summary table and metadata

counts_mat=load_summary_file(SummaryFile);
factors_mat=load_factors(FactorsFile);

factors_mat[,TrtFactorName]=relevel(factors_mat[,TrtFactorName], ControlFactorLevel);

reconciled=reconcile_samples(factors_mat, counts_mat);

trt_counts=table(factors_mat[, TrtFactorName]);	
cat("Samples by Trt Level:\n");
print(trt_counts);
cat("\n");

##############################################################################

model_string=LinearModel;
model_formula=as.formula(paste("response_variable ~ ", model_string));

# Identify the factors in the linear model
model_factors=rownames(attributes(terms(model_formula))$factors);
model_factors=setdiff(model_factors, "response_variable");
num_model_factors=length(model_factors);

# Identify the coefficients that will be estimated in the linear model
model_coefficients=colnames(attributes(terms(model_formula))$factors);
model_coefficients=setdiff(model_coefficients, "(Intercept)");
num_model_coefficients=length(model_coefficients);

##############################################################################

# Normalize the counts
norm_mat=normalize(reconciled$counts);
num_samples=nrow(norm_mat);

# Prepare space for diversity calculation results
diversity_index_names=c("Tail", "Shannon", "Simpson", "Evenness", "SimpsonsRecip");
num_indices=length(diversity_index_names);

div_mat=matrix(NA, nrow=num_samples, ncol=num_indices);
rownames(div_mat)=rownames(norm_mat);
colnames(div_mat)=diversity_index_names;

# Compute diversity for each sample and store in matrix
for(samp_ix in 1:num_samples){
	div_res=compute_diversity(norm_mat[samp_ix,]);
	div_mat[samp_ix, diversity_index_names]=div_res[diversity_index_names];
}

##############################################################################
# Compute diversity means across controls

# Extract columns we need
grouping_col=reconciled$factors[,GroupFactorName, drop=F];
series_col=reconciled$factors[,SeriesFactorName, drop=F];
treatment_col=reconciled$factors[,TrtFactorName, drop=F];

# Get number of groups and treatments
uniq_grps=sort(unique(grouping_col[,1]));
num_grps=length(uniq_grps);
uniq_trt=unique(treatment_col[,1]);

# Allocate space for per group means
grp_mean=matrix(NA, nrow=num_grps, ncol=num_indices);
rownames(grp_mean)=uniq_grps;
colnames(grp_mean)=diversity_index_names;

# Allocate space to store per group treatment level (i.e. is Patient control? H1N1? H3N2?)
trt_by_group=factor(rep(0, num_grps), levels=uniq_trt);
names(trt_by_group)=uniq_grps;

# For each group (samples from the same individual), pull sample diversities and compute mean
for(grp_id in uniq_grps){

	# Yank the samples from each group
	samp_ix=which(grouping_col==grp_id);	

	# Compute mean for each diversity index
	for(div_ix in diversity_index_names){
		grp_mean[grp_id, div_ix]=mean(div_mat[samp_ix, div_ix]);
	}

	# Keep track of the treatment level of each group
	trt_by_group[grp_id]=unique(treatment_col[samp_ix,1]);
}

# Compute mean across all the controls
control_mean=apply(grp_mean[trt_by_group==ControlFactorLevel,], 2, mean);
cat("Control means:\n");
print(control_mean);
cat("\n\n");

##############################################################################

plot_lowess=function(data, groups, series, treatments){
	#print(data);
	#print(groups);
	#print(series);
	#print(treatments);

	orig_par=par(c("mfrow"));

	# Story the x by y points by treatment
	xy_by_trt=list();

	# Identify the number of treatments
	uniq_trt=sort(unique(treatments[,1]));
	num_trt=length(uniq_trt);
	par(mfrow=c(num_trt, 1));

	# Get data name for the Y axis label
	data_name=colnames(data);
	#cat("Data name: ", data_name, "\n");

	for(cur_trt in uniq_trt){

		cat("Current Trt: ", cur_trt, "\n");

		# Get all samples for each treatment
		trt_ix=(cur_trt==treatments);	
		y=data[trt_ix];
		x=series[trt_ix];
		
		# Build 2 columns matrix for x, y
		mat=cbind(x,y);
		colnames(mat)=c("x","y");

		# Fit a lowess curve to the data
		lowess_fit=lowess(x,y);

		# Plot the lowess curve
		plot(lowess_fit, xlim=range(series), ylim=range(data), 
			type="l", col="blue", main=cur_trt,
			xlab=colnames(series), ylab=data_name
		);
		
		# Plot the observed points
		points(mat);

		# Save results
		xy_by_trt[[cur_trt]]=list(
			observed=list(x=mat[,"x"], y=mat[,"y"]), 
			fitted=list(x=lowess_fit$x, y=lowess_fit$y)
		);
	}

	par(orig_par);
	return(xy_by_trt);
}

##############################################################################


fit_lme=function(data, groups, series, treatments, formula, factors, model){

	#print(dim(data));
	#print(dim(groups));
	#print(dim(series));
	#print(dim(treatments));
	#print(dim(factors));

	data_name=colnames(data);
	groups_name=colnames(groups);
	series_name=colnames(series);

	cat("\n");
	cat("Data name: ", data_name, "\n", sep="");
	cat("Groups name: ", groups_name, "\n", sep="");
	cat("Series name: ", series_name, "\n", sep="");
	cat("\n");

	# Place columns into dataframe
	dtfr=as.data.frame(cbind(data, groups, series, treatments, factors));

	# Fixed effects formula construction
	fixed_forml=as.formula(paste(data_name, " ~ ", model, sep=""));

	# Random effects formula construction
	rand_model=paste("~1 |", groups_name);
	rand_forml=as.formula(rand_model);

	# Fit the model with lme 
	lme_fit=lme(fixed=fixed_forml, random=rand_forml, correlation=corAR1(), data=dtfr);
	print(ACF(lme_fit));

	return(lme_fit);
}


##############################################################################

pdf(paste(OutputRoot, ".div_lme.pdf", sep=""), height=11, width=8.5);

# Matrices to store per taxa results
fixed_mat=numeric();
fixed_pval_mat=numeric();

for(div_idx in diversity_index_names){

	cat("Working on Diversity Index: ", div_idx, "\n");

	#grouping_col=reconciled$factors[,GroupFactorName, drop=F];
	#series_col=reconciled$factors[,SeriesFactorName, drop=F];

	series_col=series_col+1;
	series_col[series_col==0]=1;
	
#	means_series_list=c();
	#if(TrtFactorName!=""){
	#	treatment_col=reconciled$factors[,TrtFactorName, drop=F];
	#	means_series_list=compute_means_over_series_by_treatment(
	#			div_mat[,div_idx, drop=F],	
	#			series_col,
	#			treatment_col
	#		);	
	#}

	#fits=fit_nlm(alr_data, grouping_col, series_col, treatment_col, control=ControlFactorLevel);

	#-----------------------------------------------------------------------------

	#nlm_fits=fits[["predictions"]];

	#print(div_mat);
	#print(grouping_col);
	#print(series_col);
	#print(treatment_col);

	plot_trellis(
		data=div_mat[, div_idx, drop=F],
		group_data=grouping_col,
		series_data=series_col,
		trt_levels=treatment_col,
		fits=NULL,
		means=control_mean[div_idx]
	);

	lowess_fits=plot_lowess(data=div_mat[,div_idx, drop=F], groups=grouping_col, 
		series=series_col, treatments=treatment_col);

	#print(lowess_fits);


	lme_fit=fit_lme(data=div_mat[,div_idx, drop=F], 
		groups=grouping_col,
                series=series_col, 
		treatments=treatment_col,
		factors=reconciled$factors[,model_factors],
		model=LinearModel
	);
	sum_fit=summary(lme_fit);

	plot_text(c(
		capture.output(print(lme_fit)),
		"\n",
		"\n",
		capture.output(print(sum_fit))
	), "P", title=div_idx);

	fix_eff_coeff=sum_fit$tTable;
	fixed_mat=cbind(fixed_mat, fix_eff_coeff[,"Value"]);
	fixed_pval_mat=cbind(fixed_pval_mat, fix_eff_coeff[,"p-value"]);
	
}


colnames(fixed_mat)=diversity_index_names;
colnames(fixed_pval_mat)=diversity_index_names;

print(fixed_mat);
print(fixed_pval_mat);

# Remove intercept row from table
fixed_eff_names=rownames(fixed_mat);
intercept_ix=which(fixed_eff_names=="(Intercept)");
fixed_mat=fixed_mat[-intercept_ix,];
fixed_pval_mat=fixed_pval_mat[-intercept_ix,];

## Plot heatmap
par(mfrow=c(1,1));
plot_heatmap(fixed_mat, pval_mat=fixed_pval_mat, "Diversity");
plot_heatmap(fixed_pval_mat, title="Diversity p-values");

##############################################################################

cat("\n\nDone.\n");
dev.off();
print(warnings());
q(status=0);
