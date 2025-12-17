#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
library('glmnet');
library('doMC');

source('~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r');
source('~/git/AnalysisTools/Metadata/OutputFileLibrary/OutputFileLibrary.r');

options(useFancyQuotes=F, width=80);

CVSEED=1;

params=c(
	"factor_fn", "f", 1, "character",
	"subjectid_cn", "s", 2, "character",
	"responses_fn", "r", 1, "character",
	"covariates_fn", "c" , 1, "character",
	"target_var_fn", "t", 1, "character",
	"outputroot", "o", 1, "character",
	"cvseed", "S", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"\n",
	"	-f <Factor Filename>\n",
	"	[-s <Subject ID Column Name>]\n",
	"	-r <Response List>\n",
	"	-c <Covariates List>\n",
	"	-t <Target Predictors>\n",
	"	-o <Output Filename Root>\n",
	"	-S <CrossValidation Seed, default=", CVSEED, "\n",
	"\n",
	"This script will LASSO the following logistic regression model:\n",
	"\n",
	"	[response] = [covariates] + var1 + var2 + ... + varn\n",
	"\n",
	"\n",
	"The Factor file should contain all the variables for all the subjects.\n",
	"The Response Column Name should be the name of the response you are looking\n",
	"  for targets predictors to predict.\n",
	"\n",
	"The Covariates variables are always included in the model as predictors.\n",
	"The Target Predictors should contain the list of predictors LASSO can manipulate.\n",
	"\n",
	"\n", sep="");

if(!(length(opt$factor_fn) && 
	length(opt$responses_fn) &&
	length(opt$covariates_fn) && 
	length(opt$target_var_fn) && 
	length(opt$outputroot)
)){
	cat(usage);
	q(status=-1);
}

FactorsFile=opt$factor_fn;
ResponsesFile=opt$responses_fn;
CovariatesFile=opt$covariates_fn;
TargetVarFile=opt$target_var_fn;
OutputRoot=opt$outputroot;
SubjectIDColname=opt$subjectid_cn;

if(length(opt$cvseed)){
	CVSeed=opt$cvseed;
}else{
	CVSeed=CVSEED;
}

params=capture.output({
cat("\n");
cat(script_name, "\n", sep="");
cat("\n");
cat("      Factors File: ", FactorsFile, "\n", sep="");
cat("Subject ID Colname: ", SubjectIDColname, "\n", sep="");
cat("    Responses File: ", ResponsesFile, "\n", sep="");
cat("   Covariates File: ", CovariatesFile, "\n", sep="");
cat("  Targ. Var.  File: ", TargetVarFile, "\n", sep="");
cat("       Output Root: ", OutputRoot, "\n", sep="");
cat("  Cross Valid Seed: ", CVSeed, "\n", sep="");
cat("\n");
});

print(params);
options(width=120);

##############################################################################

load_list=function(filename){
	cat("Loading List: ", filename, "\n", sep="");
	val=scan(filename, what=character(), comment.char="#");
	return(val);
}

##############################################################################
##############################################################################

factors_loaded=load_factors_file(FactorsFile, SubjectIDColname);
available_variables=colnames(factors_loaded);

targets_arr=load_list(TargetVarFile);
covariates_arr=load_list(CovariatesFile);
responses_arr=load_list(ResponsesFile);

num_targets=length(targets_arr);
num_covariates=length(covariates_arr);
num_responses=length(responses_arr);

pdf(paste(OutputRoot, ".uni_lasso.pdf", sep=""), height=11, width=9.5);

plot_text(params);

plot_text(c(
	paste("Responses [y] (", num_responses, "): ", sep=""),
	capture.output(print(responses_arr)),
	"",
	paste("Covariates [Required x] (", num_covariates, "):", sep=""),
	capture.output(print(covariates_arr)),
	"",
	paste("Target Vars [Selectable x] (", num_targets, "):", sep=""),
	capture.output(print(targets_arr))
));

all_used_variables=c(targets_arr, responses_arr, covariates_arr);

missing_variables=setdiff(all_used_variables, available_variables);
if(length(missing_variables)>0){
	cat("Error: Some variables missing from factor file:\n");
	print(missing_variables);
	quit(status=-1);
}else{
	cat("All variables found.\n");
}

factors_used=factors_loaded[,all_used_variables];
nonas=apply(factors_used, 1, function(x){ all(!is.na(x));});
factors_used=factors_used[nonas,,drop=F];

###############################################################################

cat("\nResponses:\n");
print(responses_arr);

cat("\nCovariates:\n");
print(covariates_arr);

cat("\nTargets:\n");
print(targets_arr);
cat("\n");

###############################################################################
# Generate correlation heat map 

response_matrix=factors_used[,responses_arr];
resp_corr=cor(response_matrix);
paint_matrix(resp_corr, "Response Correlation Matrix", plot_min=-1, plot_max=1,
	deci_pts=2, label_zeros=F, show_leading_zero=F);

###############################################################################

num_targets=length(targets_arr);

par(mfrow=c(3,1));
par(mar=c(5,5,7,1));

num_folds=10;
registerDoMC(num_folds);

#-----------------------------------------------------------------------------

num_predictors=num_covariates+num_targets;
cat("Num Predictors (Cov + Targets): ", num_predictors, "\n");

# Set up penalty factor for covariates, so no shrinkage is possible.
penalty_fact_arr=rep(1, num_predictors);
if(num_covariates>0){
	penalty_fact_arr[1:num_covariates]=0.0;
}

cat("Penalty Factors: 0's for covariates (i.e. no shrinkage allowed):\n");
print(penalty_fact_arr);

#-----------------------------------------------------------------------------

reorder_coefficients=function(coef_tab, cov_arr){

	coef_names=rownames(coef_tab);
	other_preds=setdiff(coef_names, c("(Intercept)", cov_arr));
	num_other_preds=length(other_preds);
	num_cov=length(cov_arr);
	
	out_tab=matrix(NA, nrow=0, ncol=ncol(coef_tab));
	colnames(out_tab)=colnames(coef_tab);

	out_tab=rbind(out_tab, coef_tab["(Intercept)",,drop=F]);

	if(num_cov>0){
		out_tab=rbind(out_tab, coef_tab[cov_arr,,drop=F]);
	}

	if(num_other_preds>0){
		other_preds_tab=coef_tab[other_preds,,drop=F];
		other_preds_tab_sorted=other_preds_tab[order(other_preds_tab[,"Pr(>|t|)"]),,drop=F];
		out_tab=rbind(out_tab, other_preds_tab_sorted);
	}

	return(out_tab);

}

#-----------------------------------------------------------------------------
# Run cv.glmnet

lasso_family="gaussian";

lasso_stats_hdr=c("NumSgnfCov", "NumSgnfSelected");
summary_matrix=matrix(0, ncol=length(lasso_stats_hdr), nrow=num_responses);

cumulative_AIC=list();
cumulative_num_pred_sel=list();

for(response_name in responses_arr){

	par(oma=c(0,0,0,0));
	par(mar=c(10,0,10,0));
	plot_title_page(response_name);

	y=factors_used[,response_name];

	#-----------------------------------------------------------------------------

	yexp=exp(y);
	ylog=log(y);
	ysqrt=sqrt(y);

	y_normp=shapiro.test(y)$p.val;
	yexp_normp=shapiro.test(yexp)$p.val;
	ylog_normp=shapiro.test(ylog)$p.val;
	ysqrt_normp=shapiro.test(ysqrt)$p.val;

	layout_mat=matrix(c(1,1,2,3,4), ncol=1);
	layout(layout_mat);

	par(mar=c(5,5,5,2));
	hist(y, main=response_name, breaks=20, cex.main=2, 
		col=ifelse(y_normp<.1, "red", "green"));
	title(main=sprintf("Shapiro-Wilks P-value: %3.4f", y_normp), line=0, cex.main=.95, font.main=3);

	par(mar=c(6,8,6,8));
	hist(yexp, main=paste("exp(", response_name, ")", sep=""), breaks=20, cex.main=1.2,
		col=ifelse(yexp_normp<.1, "pink", "lightgreen"));
	title(main=sprintf("Shapiro-Wilks P-value: %3.4f", yexp_normp), line=1.2, cex.main=.95, font.main=3);

	hist(ylog, main=paste("log(", response_name, ")", sep=""), breaks=20, cex.main=1.2,
		col=ifelse(ylog_normp<.1, "pink", "lightgreen"));
	title(main=sprintf("Shapiro-Wilks P-value: %3.4f", ylog_normp), line=1.2, cex.main=.95, font.main=3);

	hist(ysqrt, main=paste("sqrt(", response_name, ")", sep=""), breaks=20, cex.main=1.2,
		col=ifelse(ysqrt_normp<.1, "pink", "lightgreen"));
	title(main=sprintf("Shapiro-Wilks P-value: %3.4f", ysqrt_normp), line=1.2, cex.main=.95, font.main=3);

	#-----------------------------------------------------------------------------

	cat("Set Cross Validation Seed to: ", CVSeed, "\n");
	set.seed(CVSeed);

	cat("Run CV GLMNet:\n");
	y=factors_used[,response_name];
	x=as.matrix(factors_used[,c(covariates_arr, targets_arr)]);

	cvfit=cv.glmnet(x, y, family=lasso_family, 
		parallel=TRUE,
		nfolds=num_folds,
		keep=TRUE, # Keep fold information
		penalty.factor=penalty_fact_arr);
	cat("CV GLMNet Done.\n");
	
	#----------------------------------------------------------------------

	# Gives minimum mean CV error
	mindev_lambda=cvfit$lambda.min;
	log_mindev_lambda=log(mindev_lambda);

	# Gives less regularized/penalized (more selected variables) within 1se of lowest 
	stricter_lambda=cvfit$lambda.1se;
	log_stricter_lambda=log(stricter_lambda);

	# Looser lambda
	log_looser_lambda=log_mindev_lambda - (log_stricter_lambda-log_mindev_lambda);
	looser_lambda=exp(log_looser_lambda);

	msg_lambda=capture.output({
		cat("Cross Validation Randomization Seed: ", CVSeed, "\n");
		cat("Min Deviation Lambda: ", mindev_lambda, 
			"  Log():", log_mindev_lambda, "\n", sep="");
		cat("Stricter Deviation Lambda: ", stricter_lambda, 
			"  Log():", log_stricter_lambda, "\n", sep="");
		cat("Looser Deviation Lambda: ", looser_lambda, 
			"  Log():", log_looser_lambda, "\n", sep="");
	});
	print(msg_lambda);
	plot_text(msg_lambda);

	#-----------------------------------------------------------------------------
	# Get the coefficients (selected predictors will have non-zero coefficients)

	coef_tabs=list();
	coef_tabs[["Loose"]]=coef(cvfit, looser_lambda);
	coef_tabs[["Mindev"]]=coef(cvfit, mindev_lambda);
	coef_tabs[["Strict"]]=coef(cvfit, stricter_lambda);

	cutoffs=c("Loose", "Mindev", "Strict");

	coef_mat=matrix(numeric(), nrow=length(targets_arr), ncol=length(cutoffs));
	colnames(coef_mat)=cutoffs;
	rownames(coef_mat)=targets_arr;

	for(cutoff in cutoffs){
		coefs=coef_tabs[[cutoff]];
		#cat("Cutoff: ", cutoff, "\n");
		#print(coefs);
		coef_mat[targets_arr,cutoff]=coefs[targets_arr,1];

	}

	cat("\n");
	#cat("Coefficients Matrix:\n");
	#print(coef_mat);

	selected_var=apply(coef_mat, 2, function(x){sum(x!=0);});
	cat("\n");
	cat("Num of Variables Selected:\n");
	print(selected_var);


	par(mfrow=c(2,1));
	#-----------------------------------------------------------------------------
	# Plot and annotate cross validation deviances
	cat("Plotting CV: Lambda vs. Deviance...\n");
	plot(cvfit);
	title(main=paste("Response: ", response_name, sep=""), line=4);
	title(main=paste("Num Variables Selected (including ", num_covariates, " required covariates)", 
		sep=""), line=2.1, cex.main=1, font.main=1);

	abline(v=log_looser_lambda, lty="dotted");

	axis(side=3, at=log_stricter_lambda, labels="Stricter", 
		cex=.3, line=-1, col.axis="blue", tick=F);
	axis(side=3, at=log_mindev_lambda, labels="MinDev", 
		cex=.3, line=-1, col.axis="black", tick=F);
	axis(side=3, at=log_looser_lambda, labels="Looser", 
		cex=.3, line=-1, col.axis="green", tick=F);

	#-----------------------------------------------------------------------------
	# Variable counts per response
	cat("Plotting Selected Variables Per Level/Categories...\n");
	NumStrict=selected_var["Strict"];
	NumMedium=selected_var["Mindev"];
	NumLoose=selected_var["Loose"];

	cumulative_num_pred_sel[[response_name]]=selected_var;
	

	mids=barplot(c(NumLoose, NumMedium, NumStrict),
		ylim=c(0, num_targets+1),
		main=paste("Response: ", response_name, sep=""),
		las=2,
		col=c("green", "black", "blue"));
	text(mids, 
		c(NumLoose, NumMedium, NumStrict), 
		labels=c(
			sprintf("%i (%3.1f %%)", NumLoose, NumLoose/num_targets*100),
			sprintf("%i (%3.1f %%)", NumMedium, NumMedium/num_targets*100), 
			sprintf("%i (%3.1f %%)", NumStrict, NumStrict/num_targets*100)), 
		pos=3);
	abline(h=num_targets, col="black", lwd=1);
	axis(side=2, at=num_targets, labels=paste("Targets: ", num_targets, sep=""),
		las=2, cex.axis=.7);


	#-----------------------------------------------------------------------------
	# Accumulate across all responses

	stict_var_arr=names(coef_mat[coef_mat[,"Strict"]!=0,"Strict"]);
	medium_var_arr=names(coef_mat[coef_mat[,"Mindev"]!=0,"Mindev"]);
	loose_var_arr=names(coef_mat[coef_mat[,"Loose"]!=0,"Loose"]);

	varlist=list("Strict"=stict_var_arr, "Mindev"=medium_var_arr, "Loose"=loose_var_arr);

	par(mfrow=c(1,1));
	summary_lines_list=list();
	summary_lines_array=c();
	num_summary_lines=0;

	cumulative_AIC[[response_name]]=list();
	for(c_ix in cutoffs){

		# Get variable list by cutoff
		var_arr=varlist[[c_ix]];
		cumulative_AIC[[response_name]][[c_ix]]=list();

		# Fit covariates only model
		reduced_model_str=paste("y ~ ", paste(c(covariates_arr), collapse="+"), sep="");
		print(reduced_model_str);
		reduced_fit=glm(as.formula(reduced_model_str), data=as.data.frame(cbind(y,x)));
		reduced_sumfit=summary(reduced_fit);
		cumulative_AIC[[response_name]][[c_ix]][["reduced"]]=reduced_sumfit$aic;
		
		# Fit model to confirm selected variables
		model_str=paste("y ~ ", paste(c(covariates_arr, var_arr), collapse="+"), sep="");
		print(model_str);
		fit=glm(as.formula(model_str), data=as.data.frame(cbind(y,x)));
		sumfit=summary(fit);
		cumulative_AIC[[response_name]][[c_ix]][["full"]]=sumfit$aic;

		#ordered_sumfit_table=sumfit$coefficients[order(sumfit$coefficients[,"Pr(>|t|)"]),];
		ordered_sumfit_table=reorder_coefficients(sumfit$coefficients, covariates_arr);
		printCoefmat(ordered_sumfit_table);
		sumfit_txt=capture.output({printCoefmat(ordered_sumfit_table, quotes=F)});	

		# Output Coefficients/Pvalues Table
		summary_lines_list[[c_ix]]=c(
			paste("Response: ", response_name, sep=""),
			paste("Cutoff: ", c_ix, sep=""),
			"",
			sumfit_txt
			);

		summary_lines_array=c(summary_lines_array, rep(" ", 5), summary_lines_list[[c_ix]]);
	}

	if(length(summary_lines_array)<70){
		plot_text(summary_lines_array);	
	}else{
		for(c_ix in cutoffs){
			plot_text(summary_lines_list[[c_ix]]);
		}	
	}

	#-----------------------------------------------------------------------------

	par(oma=c(2,2,4,2));
	par(mfrow=c(1,3));
	par(mar=c(24,3,24,1));
	for(c_ix in cutoffs){

		# Get variable list by cutoff
		var_arr=varlist[[c_ix]];

		# Fit model to confirm selected variables
		model_str=paste("y ~ ", paste(c(covariates_arr, var_arr), collapse="+"), sep="");
		print(model_str);
		fit=lm(as.formula(model_str), data=as.data.frame(cbind(y,x)));
		sumfit=summary(fit);

		# Plot obs/pred points
		obspred_rng=range(c(y, fit$fitted.values));
		plot(y, fit$fitted.values, xlab="Observed", ylab="Predicted",
			xlim=obspred_rng, ylim=obspred_rng, main=""
			);
		title(main=c_ix, line=3, cex.main=3);

		# Draw reference and obs/pred lines
		model_fit=lm(fit$fitted.values~y);
		abline(a=0, b=1, col="blue", lty="dashed");
		abline(model_fit, col="black");

		# Output selected variable list
		if(length(var_arr)){

			outvarlist_fn=paste(OutputRoot, ".", response_name, ".", c_ix, ".lst", sep="");
			writeLines(
				var_arr, 
				outvarlist_fn
				);

			outvarfactors_fn=paste(OutputRoot, ".", response_name, ".",  c_ix, ".tsv", sep="");

			
			out_table=cbind(rownames(factors_loaded), factors_loaded[,var_arr,drop=F]);
			colnames(out_table)=c(SubjectIDColname, var_arr);
			
			# Output factors with selected variables only
			write.table(
				out_table,
				outvarfactors_fn,
				sep="\t", quote=F,
				row.names=F, col.names=T
				);
		}
	}
	mtext(response_name, side=3, line=0, outer=T, cex=2, font=2);
	
}

#------------------------------------------------------------------------------
# Output Summary of Responses and Fits 

par(oma=c(0,0,0,0));
par(mar=c(10,0,10,0));
plot_title_page("Summary across Responses");

# Convert num_pred_sel to matrix
#print(cumulative_num_pred_sel);
num_pred_sel_mat=matrix(unlist(cumulative_num_pred_sel), ncol=3, byrow=T);
rownames(num_pred_sel_mat)=names(cumulative_num_pred_sel);
colnames(num_pred_sel_mat)=names(cumulative_num_pred_sel[[1]]);
print(num_pred_sel_mat);

nonzero=apply(num_pred_sel_mat, 1, function(x){sum(x)>0;});
nz_num_pred_sel_mat=num_pred_sel_mat[nonzero,,drop=F];
print(nz_num_pred_sel_mat);

plot_text(c(
	"Num Predictors Selected:",
	capture.output({print(num_pred_sel_mat, quotes=F)}),
	"",
	"",
	"Nonzero Filtered:",
	capture.output({print(nz_num_pred_sel_mat, quotes=F)})
));

cat("\n\n");

#------------------------------------------------------------------------------

# Convert AIC to matrix
#print(cumulative_AIC);
aic_mat=matrix(unlist(cumulative_AIC), ncol=6, byrow=T);
rownames(aic_mat)=names(cumulative_AIC);
aic_mat=aic_mat[,c(1,2,4,6)];
colnames(aic_mat)=c("Reduced", "Loose", "MinDev", "Strict");
print(aic_mat);

wdiff=apply(aic_mat, 1, function(x){ any(x!=x[1])});
aic_mat_wdiff=aic_mat[wdiff,,drop=F];
print(aic_mat_wdiff);

plot_text(c(
	"AIC Summary:",
	capture.output({print(aic_mat, quotes=F)}),
	"",
	"",
	"AIC w/ Diff From Reduced",
	capture.output({print(aic_mat_wdiff, quotes=F)})
));

###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
