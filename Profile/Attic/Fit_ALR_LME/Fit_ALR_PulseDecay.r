#!/usr/bin/env Rscript

###############################################################################

NUM_VAR=25;
SEARCH_GRANULARITY=5;

library(MASS);
library('getopt');
options(useFancyQuotes=F);

params=c(
	"summary_file", "s", 1, "character",
	"factors", "f", 1, "character",
	"groups", "r", 1, "character",
	"series", "e", 1, "character",
	"control", "c", 1, "character",
	"linear_model", "m", 1, "character",
	"trt_factor", "t", 1, "character",

	"num_variables", "p", 2, "numeric",
	"outputroot", "o", 2, "character",
	"search_granularity", "n", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

# Source the library
ALR_LIB_NAME="ALR_Library.r";
#SS_PULSE_DECAY="SSpulseDecay.r";

script_name_parts=strsplit(script_name, "/")[[1]];
if(length(script_name_parts)==1){
	script_name_parts=c(".",script_name_parts);
}
script_name_parts=script_name_parts[-length(script_name_parts)];
script_path=paste(script_name_parts, collapse="/");

lib_path=paste(script_path, "/", ALR_LIB_NAME, sep="");
source(lib_path);
#lib_path=paste(script_path, "/", SS_PULSE_DECAY, sep="");
#source(lib_path);

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
	"	[-p <number of variables, default = ", NUM_VAR, ">]\n",
	"	[-o <Output filename root>]\n",
	"	[-n <search 1/graNularity, default = ", SEARCH_GRANULARITY, ">]\n",
	"\n",
	"This script fits a 3 parameter pulse decay model to each of the individuals\n",
	"and then fits a linear model to test for whether the factors had an effect\n",
	"on the pulse decay parameters.\n",
	"\n",
	"The algorithm goes:\n",
	"	1.) Convert the summary file table (-s) into ALR values per taxa\n",
	"	2.) For the control factor level (-c and -t) compute the average\n",
	"		ALR and use that as the base line.\n",
	"	3.) Generate a trellis plot across all individual over time.\n",
	"	4.) For each group of samples belonging to an individual over time\n",
	"		fit the pulse decay model by fitting the pulse decay\n",
	"		curve over the ALR values for each grouping (i.e. individual).\n",
	"	5.) Given the 3 parameters per individual, fit a linear regression\n",
	"		model using the factors specified in the linear model (-m)\n",
	"	6.) Perform an MANOVA across the 3 parameters.\n",
	"\n",
	"Search granularity specifies how many increments to use as starting positions\n",
	"for search space in each of the 3 pulse decay parameters.  Smaller numbers\n",
	"will be faster, but the solution could be a local minina.\n",
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

if(!length(opt$num_variables)){
	NumVariables=NUM_VAR;
}else{
	NumVariables=opt$num_variables;
}

if(length(opt$contains_remaining)){
	UseRemaining=T;
}else{
	UseRemaining=F;
}

if(length(opt$search_granularity)){
	SearchGranularity=opt$search_granularity;
}else{
	SearchGranularity=SEARCH_GRANULARITY;
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
cat("Search Granularity: ", SearchGranularity, "\n");
cat("\n");

options(width=120);
cat("Text Line Width: ", options()$width, "\n\n", sep="");

##############################################################################

compute_means_over_series_by_treatment=function(data, series, treatments){

	data=data[,1];
	series_name=colnames(series);
	series=series[,1];
	treatments=treatments[,1];

	uniq_trts=sort(unique(treatments));
	num_trts=length(uniq_trts);

	uniq_series_pts=sort(unique(series));
	num_series_pts=length(uniq_series_pts);

	results=list();
	for(trt in uniq_trts){

		trt_ix=(trt==treatments);
		results[[trt]]=numeric();

		for(ser in uniq_series_pts){
		
			ser_ix=(ser==series);
			tgts_ix=(ser_ix & trt_ix);

			mean_val=mean(data[tgts_ix]);

			if(!is.nan(mean_val)){
				results[[trt]]=rbind(results[[trt]], c(ser, mean_val));
			}
		}
	}

	# Returns an associative matrix of (series, mean)
	return(results);

}

#-----------------------------------------------------------------------------

pulseDecay_nobase=function(t, magnitude, peak_time, peak_width){
	val= magnitude*exp(-(log(peak_time)-log(t))^2/log(peak_width)^2);

	#cat("\n");
	#cat("    Magnitude: ", magnitude, " ");
	#cat("    Peak Time: ", peak_time, " ");
	#cat("   Peak Width: ", peak_width, "\n");
	#cat("            t: ", t, "\n");
	#cat("            Y: ", val, "\n");
	return(val);
}

MAG_MARGIN=0.25;

MIN_TIME=1/3600;
MAX_TIME=15;

MIN_PKW=1/5;
MAX_PKW=15;
#MAX_PKW=5;

fit_pulseDecay_nobase=function(obs_y, obs_time){
	
	ssd=function(param){
		magnitude=param[1];
		peak_time=param[2];
		peak_width=param[3];
		pred_y=pulseDecay_nobase(obs_time, magnitude, peak_time, peak_width);
		#ssd=(sum((pred_y-obs_y)^4));
		ssd=(sum((pred_y-obs_y)^2));
		return(ssd);
	}
	
	range_obs_y=range(obs_y);
	range_obs_time=range(obs_time);

	# Estimate magnitude and peak position
	abs_obs_y=abs(obs_y);
	max_abs=max(abs_obs_y);
	max_abs_ix=min(which(max_abs==abs_obs_y));

	mag_est=abs_obs_y[max_abs_ix];
	pkt_est=obs_time[max_abs_ix];
	pkw_est=2;

	par_init=c(mag_est, pkt_est, pkw_est);
	names(par_init)=c("magnitude", "peak_time", "peak_width");

	# Limit magnitude search range by constant margin
	spread=range_obs_y[2]-range_obs_y[1];
	mag_search_range=c(
		range_obs_y[1]-spread*MAG_MARGIN,
		range_obs_y[2]+spread*MAG_MARGIN
	);

	# Search through 3 parameters at various start positions
	best_optim_val=Inf;
	for(mag_est in seq(range_obs_y[1]-spread*.1, range_obs_y[2]+spread*.1, length.out=SearchGranularity)){
	for(pkt_est in seq(MIN_TIME, MAX_TIME, length.out=SearchGranularity)){
	for(pkw_est in seq(MIN_PKW, MAX_PKW, length.out=SearchGranularity*2)){

		par_init=c(mag_est, pkt_est, pkw_est);

		optim_res=optim(par=par_init, fn=ssd, method="L-BFGS-B", 
			lower=c(mag_search_range[1], MIN_TIME, MIN_PKW),
			upper=c(mag_search_range[2], MAX_TIME, MAX_PKW),
			control=list(trace=0)
		);

		if(best_optim_val>optim_res$value){
			best_optim_val=optim_res$value;
			best_optim_res=optim_res;
		}

	}
	}
	}

	return(best_optim_res);
}

#-----------------------------------------------------------------------------

fit_nlm=function(y_val, grouping, series, treatment, control){

	cat("Working on fitting non-linear model...\n");

	y_val=y_val[,1];
	grouping=grouping[,1];
	series=series[,1];
	treatment=treatment[,1];

	series_range=range(series);
	pred_series=seq(series_range[1], series_range[2], length.out=30);

	unique_groupings=sort(unique(grouping));
	unique_series=sort(unique(as.numeric(series)));
	unique_treatments=sort(unique(treatment));

	min_time_point=min(unique_series);
	max_time_point=max(unique_series);

	control_ix=(treatment==control);
	baseline=mean(y_val[control_ix]);
	cat("Control baseline: ", baseline, "\n");
	
	baselined_yval = y_val-baseline;

	param=numeric();
	predictions=list();

	for(cur_grp in unique_groupings){

		#cat("\n\n##############################################\n");
		#cat("Fitting ", cur_grp, " :\n", sep="");
		
		grp_ix=(cur_grp==grouping);
		cur_series=series[grp_ix];
		cur_yval=baselined_yval[grp_ix];

		#print(cur_yval);
		#print(cur_series);
		cur_param=fit_pulseDecay_nobase(cur_yval, cur_series);	
		#cat("Optimized parameter estimates:\n");
		#print(cur_param$par);
		param=rbind(param,cur_param$par);

		predictions[[cur_grp]]=cbind(
			pred_series,
			pulseDecay_nobase(pred_series, cur_param$par[1], cur_param$par[2], cur_param$par[3])+baseline
		);
		cat(".");
	}	
	cat("\n");

	rownames(param)=unique_groupings;
	colnames(param)=c("magnitude", "peak_time", "peak_width");

	results=list();
	results[["parameters"]]=param;
	results[["predictions"]]=predictions;

	return(results);
}

#-----------------------------------------------------------------------------

plot_parameter_histograms=function(nlm_parameters, grouping, treatment, title=""){

	orig_par=par(c("oma","mar","mfrow"));
	
	grouping=grouping[,1];
	treatment=treatment[,1];

	trt_names=names(treatment);
	param_names=colnames(nlm_parameters);

	#print(nlm_parameters);
	#print(grouping);
	#print(treatment);

	# Number of parameters estimated
	num_param=ncol(nlm_parameters);

	unique_treatments=sort(unique(treatment));
	num_treatments=length(unique_treatments);


	# Param Ranges:
	param_range=matrix(0, nrow=num_param, ncol=2);
	nclass=numeric();
	for(param_ix in 1:num_param){
		param_val=nlm_parameters[,param_ix];
		param_range[param_ix,]=range(param_val);
		nclass[param_ix]=nclass.Sturges(param_val)*2;
	}
	

	par(oma=c(1,1,1,1));
	par(mar=c(3,4,3,3));
	par(mfrow=c(num_treatments, num_param));

	for(cur_trt in unique_treatments){

		trt_ix=(treatment==cur_trt);
		grp_id=grouping[trt_ix];

		for(param_ix in 1:num_param){
			
			hist(nlm_parameters[grp_id, param_ix],  
				breaks=seq(param_range[param_ix, 1], param_range[param_ix, 2], length.out=nclass[param_ix]),
				xlim=param_range[param_ix,],
				ylab=cur_trt, xlab=param_names[param_ix], main=param_names[param_ix] );
	
		}
	}

	mtext(title, line=-.5, side=3, font=2, outer=T, cex=.75);

	par(orig_par);

}

##############################################################################

counts_mat=load_summary_file(SummaryFile);
factors_mat=load_factors(FactorsFile);
reconciled=reconcile_samples(factors_mat, counts_mat);

if(TrtFactorName!=""){
	trt_counts=table(factors_mat[, TrtFactorName]);	

	cat("Samples by Trt Level:\n");
	print(trt_counts);
}

##############################################################################

collapsed_res=collapse_to_top_matrix(reconciled$counts, NumVariables, remove_remaining_cat_from_top=T);
norm_mat=normalize(collapsed_res$matrix);
nz_norm_mat=reset_zeros(norm_mat);
alr_results=additive_log_ratio(nz_norm_mat, collapsed_res$remaining);

num_samples_clps_mat=nrow(alr_results$transformed);
num_categories_clps_mat=ncol(alr_results$transformed);
num_alr_categories=num_categories_clps_mat;

cat("ALR Matrix: \n");
cat("         Num Samples: ", num_samples_clps_mat, "\n");
cat("      Num Categories: ", num_categories_clps_mat, "\n");
cat("  Num ALR Categories: ", num_alr_categories, "\n");
cat("\n");

##############################################################################
# Collapse the factors by grouping, since after fitting the pulse decay curve, all the
# samples from a group (individual) will be represented by the estimated parameters

# Specify linear model where the predictors are in the meta file and response are the 
# estimated parameters.
#model_string="Host_Age + Host_Gender_isF + Sudden_onset + Symptom_cough + Flu_Type";
model_string=LinearModel;
model_formula=as.formula(paste("nlm_parameters ~ ", model_string));

# Identify the factors in the linear model
model_factors=rownames(attributes(terms(model_formula))$factors);
model_factors=setdiff(model_factors, "nlm_parameters");
num_model_factors=length(model_factors);

# Identify the coefficients that will be estimated in the linear model
model_coefficients=colnames(attributes(terms(model_formula))$factors);
model_coefficients=setdiff(model_coefficients, "(Intercept)");
num_model_coefficients=length(model_coefficients);

# Extract out the factors that are specified in the model
factors_by_grouping=reconciled$factors[, GroupFactorName];
unique_grouping=sort(unique(factors_by_grouping));
selected_factors=reconciled$factors[,c(GroupFactorName, model_factors)];
factors_by_grouping=numeric();
for(grp in unique_grouping){
	row_ix=min(which(grp==selected_factors[,GroupFactorName]));
	factors_by_grouping=rbind(factors_by_grouping, selected_factors[row_ix,]);
}
rownames(factors_by_grouping)=factors_by_grouping[,GroupFactorName];
factors_by_grouping=factors_by_grouping[,-1];

# Relevel to make control the reference
factors_by_grouping[,TrtFactorName]=relevel(factors_by_grouping[,TrtFactorName], ControlFactorLevel);

# factors_by_grouping contains all the metadata by group id, instead of sample id

##############################################################################

pdf(paste(OutputRoot, ".alr_pdm.pdf", sep=""), height=11, width=8.5);
alr_cat_names=colnames(alr_results$transformed);

# Matrices to store per taxa results
magnitude_mat=numeric();
peak_time_mat=numeric();
peak_width_mat=numeric();
magnitude_pval_mat=numeric();
peak_time_pval_mat=numeric();
peak_width_pval_mat=numeric();
manova_pval_mat=numeric();

#NumVariables=4;
for(i in 1:NumVariables){

	cur_cat_name=alr_cat_names[i];
	cat(i, ": Working on: ", cur_cat_name, "\n");

	grouping_col=reconciled$factors[,GroupFactorName, drop=F];
	series_col=reconciled$factors[,SeriesFactorName, drop=F];
	alr_data=alr_results$transformed[,i, drop=F];

	series_col=series_col+1;
	series_col[series_col==0]=1;
	
	means_series_list=c();
	treatment_col=NULL;
	if(TrtFactorName!=""){
		treatment_col=reconciled$factors[,TrtFactorName, drop=F];
		means_series_list=compute_means_over_series_by_treatment(
				alr_data,	
				series_col,
				treatment_col
			);	
	}

	fits=fit_nlm(alr_data, grouping_col, series_col, treatment_col, control=ControlFactorLevel);

	#-----------------------------------------------------------------------------

	nlm_fits=fits[["predictions"]];

	plot_trellis(
		data=alr_data,
		group_data=grouping_col,
		series_data=series_col,
		trt_levels=treatment_col,
		fits=nlm_fits
	);

	#-----------------------------------------------------------------------------

	nlm_parameters=fits[["parameters"]];

	plot_parameter_histograms(
		nlm_parameters,
		grouping_col,
		treatment_col,
		title=cur_cat_name
	);

	#-----------------------------------------------------------------------------

	# For the grouping_col ID's, extract out metadata from factor file
	
	# Perform regression and ANOVA
	fit=lm(model_formula, data=factors_by_grouping);
	sum_fit=summary(fit);
	manova_res=anova(fit);
	
	# Output Univariate regression
	plot_text(c(
		capture.output(print(fit)),
		"\n",
		"\n",
		capture.output(print(sum_fit))
	), "P", title=alr_cat_names[i]);

	# Output MANOVA
	plot_text(c(
		"Multivariate",
		capture.output(print(manova_res))
	), "P", title=alr_cat_names[i]);
	
	#-----------------------------------------------------------------------------
	
	mag_coef=sum_fit[[1]]$coefficients;
	pkt_coef=sum_fit[[2]]$coefficients;
	pkw_coef=sum_fit[[3]]$coefficients;

	mag_coef=mag_coef[-1,];
	pkt_coef=pkt_coef[-1,];
	pkw_coef=pkw_coef[-1,];

	magnitude_mat=rbind(magnitude_mat, mag_coef[,"Estimate"]);
	peak_time_mat=rbind(peak_time_mat, pkt_coef[,"Estimate"]);
	peak_width_mat=rbind(peak_width_mat, pkw_coef[,"Estimate"]);

	magnitude_pval_mat=rbind(magnitude_pval_mat, mag_coef[,"Pr(>|t|)"]);
	peak_time_pval_mat=rbind(peak_time_pval_mat, pkt_coef[,"Pr(>|t|)"]);
	peak_width_pval_mat=rbind(peak_width_pval_mat, pkw_coef[,"Pr(>|t|)"]);

	# Save the manova p-values
	rn=rownames(manova_res);
	intercept_residuals_ix=which((rn=="Residuals") | (rn=="(Intercept)"));
	manova_res=manova_res[-intercept_residuals_ix,];
	manova_pval_mat=rbind(manova_pval_mat, manova_res[, "Pr(>F)"]);
	colnames(manova_pval_mat)=rownames(manova_res);
	
}

if(NumVariables != length(alr_cat_names)){
	cat("Warning! NumVariables != length(alr_cat_names)\n");
	alr_cat_names=alr_cat_names[1:NumVariables];
}

# Save category names to matrix rows
rownames(magnitude_mat)=alr_cat_names;
rownames(peak_time_mat)=alr_cat_names;
rownames(peak_width_mat)=alr_cat_names;
rownames(magnitude_pval_mat)=alr_cat_names;
rownames(peak_time_pval_mat)=alr_cat_names;
rownames(peak_width_pval_mat)=alr_cat_names;
rownames(manova_pval_mat)=alr_cat_names;

# Plot heatmap
plot_heatmap(magnitude_mat, pval_mat=magnitude_pval_mat, "Amplitude");
plot_heatmap(magnitude_pval_mat, title="Amplitude p-values");

plot_heatmap(peak_time_mat, pval_mat=peak_time_pval_mat, "Peak Time");
plot_heatmap(peak_time_pval_mat, title="Peak Time p-values");

plot_heatmap(peak_width_mat, pval_mat=peak_width_pval_mat, "Peak Decay");
plot_heatmap(peak_width_pval_mat, title="Peak Decay p-values");

plot_heatmap(manova_pval_mat, title="MANOVA p-values");



##############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
