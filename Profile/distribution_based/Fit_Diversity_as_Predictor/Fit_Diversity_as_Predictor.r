#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');

source('~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r');
source('~/git/AnalysisTools/Metadata/OutputFileLibrary/OutputFileLibrary.r');
source('~/git/AnalysisTools/Profile/distribution_based/DiversityLibrary/DiversityLibrary.r');
source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');

options(useFancyQuotes=F);


params=c(
	"summary_file", "s", 1, "character",
	"factors", "f", 1, "character",
	"covariates", "c", 1, "character",
	"responses", "y", 1, "character",
	"required", "q", 2, "character",

	"reference_levels", "r", 2, "character",
	"outputroot", "o", 2, "character",

	"tag_name", "t", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table for taxa/function (used as Predictor/X's)>\n",
	"	-f <factors file, contains covariates and multivariate Y>\n",
	"	-c <list of covariate X's names to select from factor file (filename)>\n",
	"	-y <list of response Y's names to select from factor file (filename)>\n",
	"	[-q <required list of variables to include after NA removal>]\n",
	"\n",
	"	[-r <reference levels file for variables in factor file>]\n",
	"	[-o <output filename root>]\n",
	"\n",
	"	[-t <tag name>]\n",
	"\n",
	"This script will fit the following model:\n",
	"\n",
	" Multivariate Response = covariates + <diversity index>\n",
	"\n", sep="");

if(!length(opt$summary_file) || !length(opt$factors) || !length(opt$covariates) || 
	!length(opt$responses)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputRoot=gsub(".summary_table.xls", "", opt$summary_file);
	OutputRoot=gsub(".summary_table.tsv", "", OutputRoot);
}else{
	OutputRoot=opt$outputroot;
}

if(!length(opt$reference_levels)){
        ReferenceLevelsFile="";
}else{
        ReferenceLevelsFile=opt$reference_levels;
}

if(length(opt$required)){
	RequiredFile=opt$required;
}else{
	RequiredFile="";
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
CovariatesFile=opt$covariates;
ResponseFile=opt$responses;

cat("\n");
cat("   Summary File: ", SummaryFile, "\n", sep="");
cat("   Factors File: ", FactorsFile, "\n", sep="");
cat("Covariates File: ", CovariatesFile, "\n", sep="");
cat("  Response File: ", ResponseFile, "\n", sep="");
cat("  Required File: ", RequiredFile, "\n", sep="");
cat("\n");
cat("Output File: ", OutputRoot, "\n", sep="");
cat("\n");
cat("Reference Levels File: ", ReferenceLevelsFile, "\n", sep="");
cat("\n");

options(width=100);
cat("Text Line Width: ", options()$width, "\n", sep="");

##############################################################################

plot_fit_scatter=function(div_name, resp_name, fit, full_mcf, anpv, null_aic, full_aic){

	#print(names(fit));

	par.orig=par(no.readonly=T);

	observed=fit[["model"]][,resp_name];
	predicted=fit[["fitted.values"]];
	model_fam=fit[["family"]][[1]];

	par(mar=c(15, 15, 15, 15));

	obs_pred_range=range(c(observed, predicted));
	span=diff(obs_pred_range);
	pad=span*.07;

	plot(observed, predicted, 
		xlim=c(obs_pred_range[1]-pad, obs_pred_range[2]+pad),
		ylim=c(obs_pred_range[1]-pad, obs_pred_range[2]+2*pad),
		main="", xlab="Observed", ylab="Predicted");

	abline(a=0, b=1, col="blue");
	points(lowess(observed, predicted), type="l", col="red");
	legend(obs_pred_range[1]+pad, obs_pred_range[2]+pad, legend=c("Model", "Ideal"), 
		fill=c("red", "blue"), bty="n");

	mtext(resp_name, line=8.5, font=2, cex=3);
	mtext("Predicted vs. Observed", line=7, font=2, cex=1.1);

	# Annotate with some stats
	mtext(sprintf("Model Family: %s", model_fam), line=6, cex=.8);
	mtext(sprintf("McFadden's R^2 = %4.3f", full_mcf), line=5, cex=.9);
	mtext(sprintf("Full Model (vs. Null Model) P-val: %6.3g",
		anpv[["null_full"]]), line=3, cex=.9);
	mtext(sprintf("Diversity's Contrib to Model Improvement (Full vs. Reduced Model): P-val = %6.3g", 
		anpv[["reduced_full"]]), line=2, cex=.9);
	mtext(sprintf("Covariates-Only Model (Reduced vs. Null Model): P-val = %6.3g", 
		anpv[["null_reduced"]]), line=1, cex=.9);

	par(par.orig);
}

plot_predictor_barplot=function(full_coef_tab, reduced_coef_tab, resp_name){
	
	par.orig=par(no.readonly=T);

	# Remove intercept from tables
	predictors=setdiff(rownames(full_coef_tab), "(Intercept)");
	full_coef_tab=full_coef_tab[predictors,,drop=F];

	predictors=setdiff(rownames(reduced_coef_tab), "(Intercept)");
	reduced_coef_tab=reduced_coef_tab[predictors,,drop=F];

	# diversity name
	div_name=setdiff(rownames(full_coef_tab), rownames(reduced_coef_tab));
	cat("Diversity Name: ", div_name, "\n");
	orig_rownames=rownames(reduced_coef_tab);
	reduced_coef_tab=rbind(c(0, 0, 0, 1), reduced_coef_tab);
	rownames(reduced_coef_tab)=c(div_name, orig_rownames);

	#cat("-----------------------------------------------\n");
	#cat("Full:\n");
	#print(full_coef_tab);
	#cat("\nReduced:\n");
	#print(reduced_coef_tab);
	
	# Get signf preds
	full_signf_preds=rownames(full_coef_tab[full_coef_tab[,4]<.1,,drop=F]);
	reduced_signf_preds=rownames(reduced_coef_tab[reduced_coef_tab[,4]<.1,,drop=F]);
	combined_signf_preds=unique(c(reduced_signf_preds, full_signf_preds));
	
	# Subset only signf from both tabs	
	full_coef_tab_signf=full_coef_tab[combined_signf_preds,,drop=F];
	reduced_coef_tab_signf=reduced_coef_tab[combined_signf_preds,,drop=F];
	
	#cat("-----------------------------------------------\n");
	#cat("Only Significant:\n\n");
	#cat("Full:\n");
	#print(full_coef_tab_signf);
	#cat("\nReduced:\n");
	#print(reduced_coef_tab_signf);


	# Sort tables by pval of full model
	pval=full_coef_tab_signf[,4];
	order_inc_pval_ix=order(pval, decreasing=F);
	full_coef_tab_signf=full_coef_tab_signf[order_inc_pval_ix,,drop=F];
	reduced_coef_tab_signf=reduced_coef_tab_signf[order_inc_pval_ix,,drop=F];
	sorted_pred_names=rownames(full_coef_tab_signf);
	num_signf_preds=length(sorted_pred_names);

	if(num_signf_preds==0){
		plot_text(
			paste("No Significant Predictors for ", resp_name, sep="")
		);
		return();
	}

	#cat("-----------------------------------------------\n");
	#cat("Only Significant, sorted:\n\n");
	#cat("Full:\n");
	#print(full_coef_tab_signf);
	#cat("\nReduced:\n");
	#print(reduced_coef_tab_signf);

	# Add nlp and colors
	orig_colnames=colnames(full_coef_tab_signf);

	full_coef_tab_signf=cbind(full_coef_tab_signf, -log10(full_coef_tab_signf[,4]));
	reduced_coef_tab_signf=cbind(reduced_coef_tab_signf, -log10(reduced_coef_tab_signf[,4]));

	full_coef_tab_signf=cbind(full_coef_tab_signf, 
		sapply(full_coef_tab_signf[,"Estimate"], function(x){ifelse(x>=0, 1, -1)}));
	reduced_coef_tab_signf=cbind(reduced_coef_tab_signf, 
		sapply(reduced_coef_tab_signf[,"Estimate"], function(x){ifelse(x>=0, 1, -1)}));

	append_colnames=c(orig_colnames, "NegLogPval", "PosNeg");
	colnames(full_coef_tab_signf)=append_colnames;
	colnames(reduced_coef_tab_signf)=append_colnames;

	#cat("-----------------------------------------------\n");
	#cat("Only Significant, sorted, extra columns:\n\n");
	#cat("Full:\n");
	#print(full_coef_tab_signf);
	#cat("\nReduced:\n");
	#print(reduced_coef_tab_signf);
	
	# Generate matrix for "beside" barplots
	stacked=matrix(0, ncol=length(combined_signf_preds), nrow=2);
	colnames(stacked)=sorted_pred_names;
	rownames(stacked)=c("Full", "Reduced");
	stacked["Full", combined_signf_preds]=
		full_coef_tab_signf[combined_signf_preds, "NegLogPval"];
	stacked["Reduced", combined_signf_preds]=
		reduced_coef_tab_signf[combined_signf_preds, "NegLogPval"];

	min_pval=min(full_coef_tab_signf[,4], reduced_coef_tab_signf[,4]);

	#cat("-----------------------------------------------\n");
	#cat("Stacked Matrix for Barplot:\n");
	#print(stacked);
	#cat("-----------------------------------------------\n");

	# Generate barplot
	ref_cutoffs=c(1, 0.1, 0.05, 0.01, 0.001);
	nlp_ref_cutoffs=-log10(ref_cutoffs);
	nlp_min_pval=-log10(min(min_pval, 0.001));

	barcols=as.vector(matrix(c(
		ifelse(full_coef_tab_signf[,"PosNeg"]==1, "darkgreen", "darkred"),  
		ifelse(reduced_coef_tab_signf[,"PosNeg"]==1, "green", "red")
		), byrow=T, nrow=2));

	cat("barcols:\n");
	print(barcols);
	

	par(mfrow=c(1,1));
	par(mar=c(15,5,10,5));
	mids=barplot(stacked, names.arg=rep(c("Full", "Reduced"), num_signf_preds), 
		main="", ylab="-Log10(p-value)",
		col=barcols, cex.names=.9,
		ylim=c(0, nlp_min_pval), beside=T);
	#print(mids);

	mids_of_mids=apply(mids, 2, mean);

	abline(h=nlp_ref_cutoffs, lty="dashed", col="blue");
	axis(1, at=mids_of_mids, labels=sorted_pred_names, las=2, line=4);
	axis(4, at=nlp_ref_cutoffs, labels=ref_cutoffs, las=2);
	mtext("P-values", side=4, line=3.5, cex=.8);	

	title(main="most significant predictors", line=1.5, cex.main=1.5, font.main=1);
	title(main=resp_name, line=3, cex.main=3, font.main=2);

	par(par.orig);

}

##############################################################################
##############################################################################

pdf(paste(OutputRoot, ".div_as_pred.pdf", sep=""), height=11, width=9.5);

##############################################################################

loaded_files=load_and_reconcile_files(
	sumtab=list(fn=SummaryFile),
	factors=list(fn=FactorsFile),
	covariates=list(fn=CovariatesFile),
	grpvar=list(fn=ResponseFile),
	reqvar=list(fn=RequiredFile));

cat("Loaded Files:\n");
print(names(loaded_files));

write_file_report(loaded_files[["Report"]]);

counts=loaded_files[["SummaryTable_counts"]];
num_taxa=ncol(counts);
num_samples=nrow(counts);

diversity_mat=sample_counts_to_diversity_matrix(counts);
diversity_names=colnames(diversity_mat);

plot_histograms(diversity_mat, "Diversity");

##############################################################################

factors=loaded_files[["Factors"]];

factor_names=colnames(factors);
num_factors=ncol(factors);
factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);

cat("\n");
cat("Num Factors/Columns Loaded: ", num_factors, "\n", sep="");
cat("Num Samples for Factors Loaded: ", num_factor_samples, "\n", sep="");
print(factor_names);
cat("\n");

cat("Covariates Variables:\n");
covariates_arr=loaded_files[["Covariates"]];
print(covariates_arr);

cat("\n");

cat("Response Variables:\n");
responses_arr=loaded_files[["GroupVariables"]];
print(responses_arr);
if(length(responses_arr)==0){
	msg=c(
		"Response variables not specified.",
		"Could not run diversity as a predictor analysis."
	);
	plot_text(msg, echo=T);
	quit(status=0);
}

cat("\n");

cat("Required Variables:\n");
required_arr=loaded_files[["RequiredVariables"]];
print(required_arr);

cat("\n");

#------------------------------------------------------------------------------

responses_factors_mat=factors[,responses_arr,drop=F];
covariates_factors_mat=factors[,covariates_arr,drop=F];

###############################################################################

plot_histograms(responses_factors_mat, "Response Variables");

plot_histograms(covariates_factors_mat, "Covariates Variables");

###############################################################################

# Determine whether response is normal or binomial, or other?
cat("Determining whether response variable is continuous or boolean...\n");

assign_model_families=function(resp_var_mat){

	num_resp=ncol(resp_var_mat);

	mod_fam=character(num_resp);
	names(mod_fam)=colnames(resp_var_mat);

	for(i in 1:num_resp){
		cur_resp=resp_var_mat[,i];

		uniq_resp=unique(cur_resp);
		num_uniq_resp=length(uniq_resp);
		min_resp=min(uniq_resp);
		max_resp=max(uniq_resp);

		if(num_uniq_resp==2 && min_resp==0 && max_resp==1){
			mod_fam[i]="binomial";
		}else{
			mod_fam[i]="gaussian";
		}
	}

	return(mod_fam);
}

model_families=assign_model_families(responses_factors_mat);
num_bool_resp=sum(model_families=="binomial");
all_gaussian_resp=all(model_families=="gaussian");

cat("Model Families:\n");
print(model_families);

if(!all_gaussian_resp){
	cat("Gaussian family responses identified.\n");
}

cat("\n");

###############################################################################
# Set up and run univariate and manova


build_model_formulas=function(div_nm, cov_arr, resp_arr){

	cov_str=paste(cov_arr, collapse=" + ");

	div_mod_list=list();
	for(div in div_nm){
		resp_list=list();
		for(resp in resp_arr){

			null_model_str=paste(resp, " ~ ", 1, sep="");
			if(length(cov_arr)){
				full_model_str=paste(resp, " ~ ", div, " + ", cov_str, sep="");
				reduced_model_str=paste(resp, " ~ ", cov_str, sep="");
			}else{
				full_model_str=paste(resp, " ~ ", div, sep="");
				reduced_model_str=null_model_str;
			}

			mod_typ=list();
			mod_typ[["Null"]]=null_model_str;
			mod_typ[["Full"]]=full_model_str;
			mod_typ[["Reduced"]]=reduced_model_str;
	
			resp_list[[resp]]=mod_typ;
		}
		div_mod_list[[div]]=resp_list;
	}
	return(div_mod_list);
}

diversity_model_list=build_model_formulas(diversity_names, covariates_arr, responses_arr);
#print(diversity_model_list);

###############################################################################

compute_mod_stats=function(fit, sumfit){

	ms=list();

	ms[["mcfadden"]]=1-sumfit$deviance/sumfit$null.deviance;
	ms[["aic"]]=fit$aic;

	model_family=fit$family$family;
	if(model_family=="gaussian"){

		ss_res=sum(fit$residuals^2);
		obs_resp=fit$y;
		ss_tot=sum((obs_resp - mean(obs_resp))^2);
		unadj_r2=1-(ss_res/ss_tot);

		ms[["unadjusted"]]=unadj_r2;

		num_obs=length(obs_resp);
		num_pred_var=length(fit$coefficients)-1; # exclude intercept
		adj_r2=1-((1-unadj_r2)*(num_obs-1)/(num_obs-num_pred_var-1));

		ms[["adjusted"]]=adj_r2;

	}else{
		
		ms[["unadjusted"]]=NA;
		ms[["adjusted"]]=NA;
	}

	return(ms);

}

#------------------------------------------------------------------------------

compute_glm_manovas=function(cov, div, all_data, resp_fam){

	# cov: list of covariates
	# div: current diversity index name
	# all_data: all data cov, grp, and diversity indices
	# resp_fam: binomial, gaussian, etc.

	# See if we have enough gaussian response variables to run MANOVA
	gausfam_ix=(resp_fam=="gaussian");
	num_gaussian_resp=sum(gausfam_ix);

	cat("\n");
	if(num_gaussian_resp>2){

		gaus_resp_names=names(resp_fam[gausfam_ix]);

		cat("Running MANOVA on Response Variables: \n");
		print(gaus_resp_names);

		# Copy all the residuals from the glm fit to a single matrix
		cat("Copying resp values to single matrix.\n");
		num_samples=nrow(all_data);
		gaus_resp_mat=matrix(NA, nrow=num_samples, ncol=num_gaussian_resp);
		colnames(gaus_resp_mat)=gaus_resp_names;
		for(grn in gaus_resp_names){
			gaus_resp_mat[,grn]=all_data[,grn];
		}
		#print(gaus_resp_mat);

		full_pred_str=paste(c(div, cov), collapse=" + ");
		cat("Full Predictors String: \n");
		print(full_pred_str);

		tryCatch({

			manova_res=manova(
				as.formula(paste("gaus_resp_mat ~ ", full_pred_str, sep="")), 
				data=all_data);
			manova_sumres=summary(manova_res);

		}, error=function(e){
			cat("Error: MANOVA was attempted but not successful.\n");
			print(e);
		});

		#print(manova_res);
		print(summary(manova_res));

	}else{
		cat("Insufficient Gaussian family responses for MANOVA.\n");
		manova_sumres=NA;
	}

	return(manova_sumres);
}

#------------------------------------------------------------------------------

fit_glms=function(div_mod_lst, div_mat, fact_mat, rsp_mod_fam){

	diversity_names=names(div_mod_lst);

	div_fits_lst=list();

	samp_id=rownames(fact_mat);
	resp_pred_div_mat=cbind(div_mat[samp_id,,drop=F], fact_mat[samp_id,,drop=F]); 
	#print(resp_pred_div_mat);

	for(div_nm in diversity_names){
		

		cat("Fitting ", div_nm, " models...\n", sep="");

		rsp_names=names(div_mod_lst[[div_nm]]);
		
		glm_fits_by_rsp=list();
		for(rsp_nm in rsp_names){
			rsp_fam=rsp_mod_fam[rsp_nm];
			cat("  Fitting Resp: ", rsp_nm, " (family=", rsp_fam, ")\n", sep="");
			
			# Compute GLMs
			glm_fits=list();
			glm_sumfits=list();
			glm_r2s=list();
			for(models in c("Null", "Full", "Reduced")){

				cat("    Fitting Model: ", models, "\n", sep="");
				glm_fits[[models]]=glm(
					as.formula(div_mod_lst[[div_nm]][[rsp_nm]][[models]]),
					data=resp_pred_div_mat,
					family=rsp_fam
					);
				glm_sumfits[[models]]=summary(glm_fits[[models]]);
				#print(glm_sumfits[[models]]);

				glm_r2s[[models]]=compute_mod_stats(
					glm_fits[[models]], glm_sumfits[[models]]);
			}
			
			glm_fits_by_rsp[[rsp_nm]][["fits"]]=glm_fits;	
			glm_fits_by_rsp[[rsp_nm]][["sumfits"]]=glm_sumfits;	
			glm_fits_by_rsp[[rsp_nm]][["r2"]]=glm_r2s;	

			# Compute ANOVAs
			anova_null_reduced_res=anova(
				glm_fits[["Null"]], 
				glm_fits[["Reduced"]], 
				test="Chi");
			anova_null_reduced_pval=anova_null_reduced_res[["Pr(>Chi)"]][2];

			anova_null_full_res=anova(
				glm_fits[["Null"]], 
				glm_fits[["Full"]], 
				test="Chi");
			anova_null_full_pval=anova_null_full_res[["Pr(>Chi)"]][2];

			anova_reduced_full_res=anova(
				glm_fits[["Reduced"]], 
				glm_fits[["Full"]], 
				test="Chi");
			anova_reduced_full_pval=anova_reduced_full_res[["Pr(>Chi)"]][2];

			anovas=list();
			anovas[["null_reduced"]]=anova_null_reduced_pval;
			anovas[["null_full"]]=anova_null_full_pval;
			anovas[["reduced_full"]]=anova_reduced_full_pval;
			glm_fits_by_rsp[[rsp_nm]][["anova_pvals"]]=anovas;

		}

		div_fits_lst[[div_nm]]=list();
		div_fits_lst[[div_nm]][["glms"]]=glm_fits_by_rsp;
		div_fits_lst[[div_nm]][["manova"]]=
			compute_glm_manovas(
				covariates_arr,
				div_nm,
				resp_pred_div_mat,
				rsp_mod_fam);
	}

	return(div_fits_lst);

}

fit_list=fit_glms(diversity_model_list, diversity_mat, factors, model_families);


#print(fit_list);

###############################################################################

print_list_as_tree=function(lst, lev=0, max_lev=Inf){
	
	if(!is.list(lst)){
		return();
	}else{
		if(lev>max_lev){
			return();
		}

		for(n in names(lst)){
			cat(paste(rep("   ", lev), collapse=""), n, "\n", sep="");
			print_list_as_tree(lst[[n]], lev=lev+1, max_lev=max_lev);
		}
	}

}

###############################################################################

#print_list_as_tree(fit_list, max_lev=4);


accumulate_coef_pval_matrices=function(fit_lst){

	cat("\n");
	cat("Accumulating matrices from across list.\n");

	div_names=names(fit_lst);
	example_glm=fit_lst[[div_names[1]]][["glms"]]
	resp_names=names(example_glm);

	cat("Diversity Indices:\n");
	print(div_names);
	cat("\n");
	cat("Response Names:\n");
	print(resp_names);
	cat("\n");

	# Allocate matrices
	num_div_ix=length(div_names);
	num_resp=length(resp_names);
	pval_mat=matrix(NA, nrow=num_div_ix, ncol=num_resp, 
			dimnames=list(div_names, resp_names));
	coef_mat=pval_mat;

	for(div_ix in div_names){
		#cat("Extracting: ", div_ix, "\n", sep="");

		glms=fit_list[[div_ix]][["glms"]];

		for(rsp_ix in resp_names){

			cat("  Response: ", rsp_ix, "\n");

			full_sumfits=glms[[rsp_ix]][["sumfits"]][["Full"]];
			print(full_sumfits);

			regress_tab=full_sumfits$coefficients;
			print(regress_tab);

			coef_mat[div_ix, rsp_ix]=regress_tab[div_ix, "Estimate"];
			pval_mat[div_ix, rsp_ix]=regress_tab[div_ix, 4];
		}
	}	
	
	matrices=list();
	matrices[["pval"]]=pval_mat;
	matrices[["coef"]]=coef_mat;

	return(matrices);
}

pv_cf_mat_full_list=accumulate_coef_pval_matrices(fit_list);


###############################################################################

accumulate_coef_pval_predictors_matrices_by_diversity=function(fit_lst){

	div_names=names(fit_lst);
	example_glm=fit_lst[[div_names[1]]][["glms"]]
	resp_names=names(example_glm);
	pred_names=setdiff(
		rownames(example_glm[[1]][["sumfits"]][["Full"]][["coefficients"]]),
		c("(Intercept)", div_names[1])
		);

	cat("Diversity Indices:\n");
	print(div_names);
	cat("\n");
	cat("Response Names:\n");
	print(resp_names);
	cat("\n");
	cat("Predictor Names:\n");
	print(pred_names);
	cat("\n");
	

	# Allocate matrices
	num_div_ix=length(div_names);
	num_resp=length(resp_names);
	num_pred=length(pred_names);
	mat_tmp=matrix(NA, nrow=num_pred+1, ncol=num_resp, 
			dimnames=list(c("",pred_names), resp_names));

	cfpv_list=list();
	for(div_ix in div_names){
		cat("Extracting: ", div_ix, "\n", sep="");

		glms=fit_list[[div_ix]][["glms"]];
		coef_mat=mat_tmp;
		pval_mat=mat_tmp;

		rownames(coef_mat)=c(div_ix, pred_names);
		rownames(pval_mat)=c(div_ix, pred_names);

		for(rsp_ix in resp_names){

			cat("  Response: ", rsp_ix, "\n");

			full_sumfits=glms[[rsp_ix]][["sumfits"]][["Full"]];
			regress_tab=full_sumfits$coefficients;

			for(prd_ix in c(div_ix, pred_names)){
				coef_mat[prd_ix, rsp_ix]=regress_tab[prd_ix, "Estimate"];
				pval_mat[prd_ix, rsp_ix]=regress_tab[prd_ix, 4];
			}
		}

		mat_lst=list();
		mat_lst[["coef"]]=coef_mat;
		mat_lst[["pval"]]=pval_mat;
		cfpv_list[[div_ix]]=mat_lst;
	}	
	
	return(cfpv_list);

}

cf_pv_mat_list=accumulate_coef_pval_predictors_matrices_by_diversity(fit_list);


###############################################################################

accumulate_model_improv_stats=function(fit_lst){

	cat("\n");
	cat("Accumulating Model Improvement (R^2 and AIC) Stats across list.\n");

	div_names=names(fit_lst);
	example_glm=fit_lst[[div_names[1]]][["glms"]]
	resp_names=names(example_glm);

	cat("Diversity Indices:\n");
	print(div_names);
	cat("\n");
	cat("Response Names:\n");
	print(resp_names);
	cat("\n");

	# Allocate matrices
	num_div_idx=length(div_names);
	num_resp=length(resp_names);

	# AIC Matrix Template
	aic_coln=c("Null", "Reduced", "Full");
	aic_mat=matrix(NA, nrow=num_resp, ncol=length(aic_coln),
			dimnames=list(resp_names, aic_coln));

	# R2 Matrix Template
	r2_coln=c("McF Reduced", "McF Full", "Unadj Reduced", "Unadj Full", "Adj Reduced", "Adj Full");
	r2_mat=matrix(NA, nrow=num_resp, ncol=length(r2_coln),
			dimnames=list(resp_names, r2_coln));

	# Anova Matrix Template
	anova_coln=c("Null Full", "Null Reduced", "Reduced Full");
	anova_mat=matrix(NA, nrow=num_resp, ncol=length(anova_coln),
			dimnames=list(resp_names, anova_coln));

	stats_by_div=list();

	for(div_ix in div_names){
		cat("Extracting: ", div_ix, "\n", sep="");

		glms=fit_list[[div_ix]][["glms"]];

		for(rsp_ix in resp_names){

			cat("  Response: ", rsp_ix, "\n");


			aic_mat[rsp_ix, "Null"]=
				glms[[rsp_ix]][["fits"]][["Null"]][["aic"]];
			aic_mat[rsp_ix, "Reduced"]=
				glms[[rsp_ix]][["fits"]][["Reduced"]][["aic"]];
			aic_mat[rsp_ix, "Full"]=
				glms[[rsp_ix]][["fits"]][["Full"]][["aic"]];


			r2_mat[rsp_ix, "McF Reduced"]=
				glms[[rsp_ix]][["r2"]][["Reduced"]][["mcfadden"]]; 
			r2_mat[rsp_ix, "McF Full"]=
				glms[[rsp_ix]][["r2"]][["Full"]][["mcfadden"]]; 
			r2_mat[rsp_ix, "Unadj Reduced"]=
				glms[[rsp_ix]][["r2"]][["Reduced"]][["unadjusted"]]; 
			r2_mat[rsp_ix, "Unadj Full"]=
				glms[[rsp_ix]][["r2"]][["Full"]][["unadjusted"]]; 
			r2_mat[rsp_ix, "Adj Reduced"]= 
				glms[[rsp_ix]][["r2"]][["Reduced"]][["adjusted"]]; 
			r2_mat[rsp_ix, "Adj Full"]=
				glms[[rsp_ix]][["r2"]][["Full"]][["adjusted"]]; 


			anova_mat[rsp_ix, "Null Full"]=
				nf_anova_pv=glms[[rsp_ix]][["anova_pvals"]][["null_full"]];
			anova_mat[rsp_ix, "Null Reduced"]=
				nr_anova_pv=glms[[rsp_ix]][["anova_pvals"]][["null_reduced"]];
			anova_mat[rsp_ix, "Reduced Full"]=
				rf_anova_pv=glms[[rsp_ix]][["anova_pvals"]][["reduced_full"]];


		}

		mod_stats=list();
		mod_stats[["aic"]]=aic_mat;
		mod_stats[["r2"]]=r2_mat;
		mod_stats[["anova"]]=anova_mat;
		stats_by_div[[div_ix]]=mod_stats;
	}	
	
	
	return(stats_by_div);

}

stats_by_div=accumulate_model_improv_stats(fit_list);
print(stats_by_div);

###############################################################################

format_pval=function(x){
	res=ifelse(x>=.001,
		sprintf("%11.3f", x),
		sprintf("%11.03e", x));
	return(res);
}

###############################################################################
#
# Report by diversity index:

#------------------------------------------------------------------------------

reformat_coef_tab=function(coef_tab){

	num_col=ncol(coef_tab);
	num_coef=nrow(coef_tab);

	out_tab=matrix(character(), nrow=num_coef, ncol=num_col+1);
	rownames(out_tab)=rownames(coef_tab);
	colnames(out_tab)=c(colnames(coef_tab), "Signf");

	for(i in 1:3){
		out_tab[,i]=sprintf("%8.04f", coef_tab[,i]);
	}

	out_tab[,4]=ifelse(coef_tab[,4]>=0.001,
		sprintf("%11.3f", coef_tab[,4]),
		sprintf("%11.03e", coef_tab[,4]));

	out_tab[,"Signf"]=sapply(coef_tab[,4], signf_char);

	return(out_tab);
}

report_fit=function(fit_rec, div_name, rsp_name){

	cat("Reporting GLMs: ", div_name, " : ", rsp_name, "\n");
	#print(fit_rec);

	# Set up coefficients / p-val table
	full_coef_tab=fit_rec[["sumfits"]][["Full"]][["coefficients"]];
	reduced_coef_tab=fit_rec[["sumfits"]][["Reduced"]][["coefficients"]];

	cat("Full Model Coefficients:\n");
	print(full_coef_tab);

	cat("Reduced Model Coefficients:\n");
	print(reduced_coef_tab);

	reduced_predictors=rownames(reduced_coef_tab);
	cat("Reduced Predictors:\n");
	print(reduced_predictors);

	full_coef_tab_char=reformat_coef_tab(full_coef_tab);
	reduced_coef_tab_char=reformat_coef_tab(reduced_coef_tab);

	print(full_coef_tab_char, quote=F);
	print(reduced_coef_tab_char, quote=F);

	# Pull together AIC and R^2's
	r2=fit_rec[["r2"]];

	full_mcfadden=r2[["Full"]][["mcfadden"]];
	full_adj=r2[["Full"]][["adjusted"]];
	full_unadj=r2[["Full"]][["unadjusted"]];

	reduced_mcfadden=r2[["Reduced"]][["mcfadden"]];
	reduced_adj=r2[["Reduced"]][["adjusted"]];
	reduced_unadj=r2[["Reduced"]][["unadjusted"]];

	null_aic=r2[["Null"]][["aic"]];
	full_aic=r2[["Full"]][["aic"]];
	reduced_aic=r2[["Reduced"]][["aic"]];
	nf_aic_diff=null_aic-full_aic;
	rf_aic_diff=reduced_aic-full_aic;

	# Pull together ANOVA p-values
	anpv=fit_rec[["anova_pvals"]];
	anpv_null_reduced=anpv[["null_reduced"]];
	anpv_null_full=anpv[["null_full"]];
	anpv_reduced_full=anpv[["reduced_full"]];

	#----------------------------------------------------------------------

	full_reduced_coef_comp=c(
		paste("Univariate Response: ", rsp_name, sep=""),
		"",
		paste("Comparison of the Covariates portion of w/ and w/o ", div_name, ":", sep=""),
		"",
		"---------------------------------------------------------------------------------------",
		"",
		"Predictors of Reduced Model (Covariates-only Model):",
		sprintf("McFadden's R^2: %4.3f, AIC: %.1f", reduced_mcfadden, reduced_aic),
		"",
		capture.output(print(reduced_coef_tab_char, quote=F)),
		"",
		"---------------------------------------------------------------------------------------",
		"",

		"Predictors Portion of Full Model (Covariates + Diversity Model):",
		paste("  (Contribution of Covariates when controlling for ", div_name, ".)", sep=""),
		sprintf("McFadden's R^2: %4.3f, AIC: %.1f", full_mcfadden, full_aic),
		"",
		capture.output(print(full_coef_tab_char, quote=F)),
		"",
		"---------------------------------------------------------------------------------------",
		"",
		paste("Model Improvement ANOVA (Full vs. Reduced): ", 
			sprintf("%6.02g", anpv_reduced_full), 
			" ", signf_char(anpv_reduced_full), sep=""),
		"",
		paste("Diff. in AIC: ", round(reduced_aic,2), " [reduced] - ", round(full_aic,2), 
			" [full] = ", round(rf_aic_diff,2), " [dAIC]", sep=""),
		"If dAIC is > 2, then full model is substantially better than reduced model.",
		paste("Full model is ", 
			ifelse(rf_aic_diff>2, "", "NOT "), 
			"substantially better then reduced model.", sep=""),
		"",
		sprintf("McFadden's R^2: %4.3f [full] - %4.3f [reduced] = %4.3f [delta]",
			full_mcfadden, reduced_mcfadden, (full_mcfadden-reduced_mcfadden))	
		
	);

	plot_text(full_reduced_coef_comp, echo=T);
	mtext(div_name, side=1, line=0, outer=T, font=2, col="darkred");

	#----------------------------------------------------------------------

	plot_predictor_barplot(full_coef_tab, reduced_coef_tab, rsp_name);
	mtext(div_name, side=1, line=0, outer=T, font=2, col="darkred");

	plot_fit_scatter(
		div_name, rsp_name,
		fit_rec[["fits"]][["Full"]], 
		full_mcfadden, anpv, null_aic, full_aic);
	mtext(div_name, side=1, line=0, outer=T, font=2, col="darkred");

}

###############################################################################
# Manova output

format_manova_tab=function(mt){

	format_pval=function(x){
		res=ifelse(x>=.001,
			sprintf("%11.3f", x),
			sprintf("%11.03e", x));
		return(res);
	}

	manova_tab_char=matrix("", ncol=ncol(mt), nrow=nrow(mt));
	colnames(manova_tab_char)=colnames(mt);
	rownames(manova_tab_char)=rownames(mt);

	manova_tab_char[,"Df"]=mt[,"Df"];
	manova_tab_char[,"Pillai"]=sprintf("%7.4f", mt[,"Pillai"]);
	manova_tab_char[,"approx F"]=sprintf("%8.4f", mt[,"approx F"]);
	manova_tab_char[,"num Df"]=mt[,"num Df"];
	manova_tab_char[,"den Df"]=mt[,"den Df"];
	manova_tab_char[,"Pr(>F)"]=format_pval(mt[,"Pr(>F)"]);

	manova_tab_char=cbind(manova_tab_char, sapply(mt[,"Pr(>F)"], signf_char));
	colnames(manova_tab_char)=c(colnames(mt), "signf");

	return(manova_tab_char);

}

plot_manova_barplot=function(manova_tab){

	ref_cutoffs=c(1, 0.1, 0.05, 0.01, 0.001);

	pval=manova_tab[,"Pr(>F)"];
	pred=rownames(manova_tab);
	min_pval=min(pval);

	nlp_ref_cutoffs=-log10(ref_cutoffs);
	nlp_min_pval=-log10(min(min_pval, 0.001));
	nlp_pvals=-log10(pval);

	signf_col=sapply(pval, function(x){ ifelse(x<0.1, "blue", "grey")});

	par(mfrow=c(1,1));
	par(mar=c(15,5,10,5));
	mids=barplot(nlp_pvals, names.arg="", main="", ylab="-Log10(p-value)",
		col=signf_col,
		ylim=c(0, nlp_min_pval));
	abline(h=nlp_ref_cutoffs, lty="dashed", col="orange");
	axis(1, at=mids, labels=pred, las=2);
	axis(4, at=nlp_ref_cutoffs, labels=ref_cutoffs, las=2);
	mtext("P-values", side=4, line=3.5, cex=.8);	

	title(main="MANOVA: Most Significant Predictors", line=1.5, cex.main=1.5, font.main=1);
}

report_manova=function(manova_res_sum, title){

	print(names(manova_res_sum));
	resp_var=colnames(manova_res_sum$SS$Residuals);

	if(length(manova_res_sum)>1){

		manova_tab=manova_res_sum$stats;
		pred_names=rownames(manova_tab);
		pred_names=setdiff(pred_names, "Residuals");
		num_pred_var=length(pred_names);
		manova_tab=manova_tab[pred_names,,drop=F];

		ord_pval_ix=order(manova_tab[,"Pr(>F)"], decreasing=F);
		manova_tab=manova_tab[ord_pval_ix,,drop=F];

		manova_tab_char=format_manova_tab(manova_tab);
		print(manova_tab_char, quote=F);

		manout=c(
			title,
			"",
			"MANOVA on Gaussian Family Variables",
			capture.output(print(resp_var, quote=F)),
			"",
			paste("Num Samples: ", num_samples, sep=""),
			paste("Num Covariates+",title, " Predictors: ", num_pred_var, sep=""),
			"",
			capture.output(print(manova_tab_char, quote=F))
		);

		plot_text(manout, echo=T);
		mtext(title, side=1, line=0, outer=T, font=2, col="darkred");

		#--------------------------------------------------------------
		
		plot_manova_barplot(manova_tab);
		mtext(title, side=1, line=0, outer=T, font=2, col="darkred");


	}else{
		plot_text(c(
			"MANOVA could not be performed.",
			"",
			"If <2 response variables were assumed to be Gaussian Family",
			"or there were insufficient samples for the number of predictors",
			"and response variables, then the MANOVA could not be performed.",
			"",
			paste("Num Samples: ", num_samples, sep=""),
			paste("Num Covariates+",title, " Predictors: ", num_pred_var, sep="")
		));
		mtext(title, side=1, line=0, outer=T, font=2, col="darkred");
	}

}


#------------------------------------------------------------------------------

plot_coef_pval_heatmaps=function(cfpv_mat, divn){

	coef_mat=cfpv_mat[["coef"]];
	pval_mat=cfpv_mat[["pval"]];

	par.orig=par(no.readonly=T);
	par(oma=c(1,0,0,0));

	paint_matrix(coef_mat,
		title="Full Models: Predictor Coefficients: All",value.cex=2, deci_pts=2);
	mtext(divn, side=1, line=0, outer=T, font=2, col="darkred");

	paint_matrix(pval_mat,
		title="Full Models: Predictor P-values: All",value.cex=2, deci_pts=2,
		plot_min=0, plot_max=1, high_is_hot=F);
	mtext(divn, side=1, line=0, outer=T, font=2, col="darkred");

	for(cutoff in c(0.10, 0.05, 0.01, 0.001)){
		masked_coef_mat=mask_matrix(coef_mat, pval_mat, cutoff, 0);
		paint_matrix(masked_coef_mat, 
			title=paste("Full Models: Coefficients: p-val < ", cutoff, sep=""),
			value.cex=2, deci_pts=2, label_zeros=F);
		mtext(divn, side=1, line=0, outer=T, font=2, col="darkred");
	}

	par(par.orig);
		
}

plot_model_fitness_heatmaps=function(mod_stats, divn){

	aic_mat=mod_stats[["aic"]];
	r2_mat=mod_stats[["r2"]][,c("McF Reduced","McF Full")];
	anova_mat=mod_stats[["anova"]];

	
	par.orig=par(no.readonly=T);
	par(oma=c(1,0,0,0));

	# AIC
	orig_cnames=colnames(aic_mat);
	aic_mat=cbind(aic_mat, aic_mat[,"Full"]-aic_mat[,"Reduced"]);
	colnames(aic_mat)=c(orig_cnames, "\"Difference\"");
	sort_best_full_ix=order(aic_mat[,"\"Difference\""]);
	paint_matrix(aic_mat[sort_best_full_ix,,drop=F],
		title="Full Model Quality: AIC",
		subtitle="Sorted by Full-Reduced Difference (>2 is Signf, Lower is Better)",
		 value.cex=2, deci_pts=1);
	mtext(divn, side=1, line=0, outer=T, font=2, col="darkred");

	# R2
	orig_cnames=colnames(r2_mat);
	r2_mat=cbind(r2_mat, r2_mat[,"McF Full"]-r2_mat[,"McF Reduced"]);
	colnames(r2_mat)=c(orig_cnames, "\"Difference\"");
	sort_best_full_ix=order(r2_mat[, "\"Difference\""], decreasing=T);
	paint_matrix(r2_mat[sort_best_full_ix,,drop=F],
		plot_min=0, plot_max=1,
		title="Full Model Quality: McFadden's R^2", 
		subtitle="Sorted by Full-Reduced McF Difference (Large Difference is Better)",
		value.cex=2, deci_pts=3);
	mtext(divn, side=1, line=0, outer=T, font=2, col="darkred");

	# ANOVA
	sort_best_full_ix=order(anova_mat[,"Reduced Full"]);
	paint_matrix(anova_mat[sort_best_full_ix,,drop=F],
		log_col=T, high_is_hot=F,
		title="Model Improvement: ANOVA P-values", 
		subtitle="Sorted by Reduced-Full P-Values (Smaller P-val is Better)",
		value.cex=2, deci_pts=3);
	mtext(divn, side=1, line=0, outer=T, font=2, col="darkred");

	par(par.orig);

}


#------------------------------------------------------------------------------

generate_fit_details=function(fit_lst){

	div_names=names(fit_lst);

	for(div in div_names){

		cur_glm=fit_lst[[div]][["glms"]];		

		for(rsp in names(cur_glm)){
			report_fit(cur_glm[[rsp]], div, rsp);

		}

		report_manova(manova_res_sum=fit_lst[[div]][["manova"]], div);

		plot_coef_pval_heatmaps(cf_pv_mat_list[[div]], div);
		plot_model_fitness_heatmaps(stats_by_div[[div]], div);

	}

}


generate_fit_details(fit_list);


###############################################################################
###############################################################################
# Combine diversity indices together for each response variable

combine_fits_across_div_into_matrices=function(fit_lst){

	#print(fit_lst);

	cat("Combining Fits Across Diversities:\n");
	div_names=names(fit_lst);
	
	resp_names=names(fit_lst[[1]][["glms"]]);
	
	cat("Diversity Names:\n");
	print(div_names);
	num_div=length(div_names);

	cat("Response Names:\n");
	print(resp_names);
	num_rsp=length(resp_names);

	fit_mat=matrix(NA, nrow=num_rsp, ncol=num_div, 
		dimnames=list(resp_names, div_names));

	aic_mat=fit_mat;
	mcf_mat=fit_mat;
	pvl_mat=fit_mat;

	for(div in div_names){
		for(rsp in resp_names){
			full_r2_rec=fit_lst[[div]][["glms"]][[rsp]][["r2"]][["Full"]];
			full_aic=full_r2_rec[["aic"]];
			full_mcf=full_r2_rec[["mcfadden"]];

			aic_mat[rsp, div]=full_aic;
			mcf_mat[rsp, div]=full_mcf;

			red_ful_pv=fit_lst[[div]][["glms"]][[rsp]][["anova_pvals"]][["reduced_full"]];
			pvl_mat[rsp, div]=red_ful_pv;
			
		}
	}

	#print(aic_mat);
	#print(mcf_mat);

	rec=list();
	rec[["aic"]]=aic_mat;
	rec[["mcf"]]=mcf_mat;
	rec[["red_ful_pval"]]=pvl_mat;
	
	return(rec);

}

fit_mat_list=combine_fits_across_div_into_matrices(fit_list);

plot_div_fit_heatmaps=function(fit_mat_lst){

	add_margins=function(mat, rowm, colm){
		
		all_mean=mean(mat);	
		orig_rn=rownames(mat);
		orig_cn=colnames(mat);

		mat_wmar=cbind(mat, rowm);
		mat_wmar=rbind(mat_wmar, c(colm, all_mean));

		colnames(mat_wmar)=c(orig_cn, "\"Mean\"");
		rownames(mat_wmar)=c(orig_rn, "\"Mean\"");

		return(mat_wmar);
	
	}

	#----------------------------------------------------------------------
	
	mat=fit_mat_lst[["aic"]];

	resp_margins=apply(mat, 1, mean);
	div_margins=apply(mat, 2, mean);

	resp_margins_sort_ix=order(resp_margins, decreasing=F);
	div_margins_sort_ix=order(div_margins, decreasing=F);

	sorted_mat=mat[resp_margins_sort_ix, div_margins_sort_ix, drop=F];
	resp_margins=resp_margins[resp_margins_sort_ix];
	div_margins=div_margins[div_margins_sort_ix];

	sorted_wmargins=add_margins(mat=sorted_mat, rowm=resp_margins, colm=div_margins);

	paint_matrix(sorted_wmargins, 
		title="Comparison of Full AIC Across Diversity Ind.",
		subtitle="(Lower AIC is a better model)",
		deci_pts=1,
	);

	#----------------------------------------------------------------------

	mat=fit_mat_lst[["mcf"]];

	resp_margins=apply(mat, 1, mean);
	div_margins=apply(mat, 2, mean);

	resp_margins_sort_ix=order(resp_margins, decreasing=T);
	div_margins_sort_ix=order(div_margins, decreasing=T);

	sorted_mat=mat[resp_margins_sort_ix, div_margins_sort_ix,drop=F];
	resp_margins=resp_margins[resp_margins_sort_ix];
	div_margins=div_margins[div_margins_sort_ix];

	sorted_wmargins=add_margins(mat=sorted_mat, rowm=resp_margins, colm=div_margins);

	paint_matrix(sorted_wmargins, 
		title="Comparison of Full McFadden's R^2 Across Diversity Ind.",
		subtitle="(Greater R^2 is a better model)",
		deci_pts=3,
		plot_min=0, plot_max=1
	);

	#----------------------------------------------------------------------

	mat=fit_mat_lst[["red_ful_pval"]];

	comb_pv=function(x){
		nlog=-log(x);
		mean_nlog=mean(nlog);
		cpv=exp(-mean_nlog);
		return(cpv);
	}

	resp_margins=apply(mat, 1, comb_pv);
	div_margins=apply(mat, 2, comb_pv);

	resp_margins_sort_ix=order(resp_margins, decreasing=F);
	div_margins_sort_ix=order(div_margins, decreasing=F);

	sorted_mat=mat[resp_margins_sort_ix, div_margins_sort_ix,drop=F];
	resp_margins=resp_margins[resp_margins_sort_ix];
	div_margins=div_margins[div_margins_sort_ix];

	sorted_wmargins=add_margins(mat=sorted_mat, rowm=resp_margins, colm=div_margins);

	paint_matrix(sorted_wmargins, 
		title="Reduced-Full Improv. ANOVA P-Vals Across Div. Ind.",
		subtitle="(More significant p-values, great contribution of Div. Ind.)",
		deci_pts=3,
		high_is_hot=F,
		plot_min=0, plot_max=1
	);
		
}

plot_title_page("Comparison of Models\nAcross Diversity Indices");
plot_div_fit_heatmaps(fit_mat_list);

###############################################################################




###############################################################################

pull_div_predictors_tab=function(fit_lst){

	# Pull predictors coef and pval into table
	#print(fit_lst);

	div_names=names(fit_lst);
	resp_names=names(fit_lst[[1]][["glms"]]);
	
	cat("Diversity Names:\n");
	print(div_names);
	num_div=length(div_names);

	cat("Response Names:\n");
	print(resp_names);
	num_rsp=length(resp_names);

	fit_mat=matrix(NA, nrow=num_rsp, ncol=num_div, 
		dimnames=list(resp_names, div_names));

	coef_mat=fit_mat;
	pval_mat=fit_mat;

	for(div in div_names){
		for(rsp in resp_names){
			full_coef_tab=
				fit_lst[[div]][["glms"]][[rsp]][["sumfits"]][["Full"]][["coefficients"]];

			coef_mat[rsp, div]=full_coef_tab[div, "Estimate"];
			pval_mat[rsp, div]=full_coef_tab[div, 4];
			
		}
	}

	rec=list();
	rec[["pval"]]=pval_mat;
	rec[["coef"]]=coef_mat;
	
	return(rec);

}

div_coef_pval_list=pull_div_predictors_tab(fit_list);	
print(div_coef_pval_list);	

#------------------------------------------------------------------------------

generate_div_resp_heatmaps=function(div_cf_pv_lst){

	pval_mat=div_cf_pv_lst[["pval"]];
	coef_mat=div_cf_pv_lst[["coef"]];

	coef_mag=max(abs(coef_mat));

	for(cutoff in c(1, .1, .05, .01, .001)){

		masked_coef_mat=mask_matrix(coef_mat, pval_mat, cutoff, 0);
		paint_matrix(masked_coef_mat, 
			title=paste("[Full Models] Diversity Coef w/ p-val < ", cutoff, sep=""),
			plot_max=coef_mag, plot_min=-coef_mag,
			value.cex=2, deci_pts=2, label_zeros=F);

	}

}

plot_title_page("Comparisons of\nDiversity-to-Response\nAssociations");
generate_div_resp_heatmaps(div_coef_pval_list);

#------------------------------------------------------------------------------

draw_list=function(abbr_tab, coloffset){

	# the abbr_tab has is a matrix with: col1=estimate, col2=pval
	# rownames are the predictors

	chw=par()$cxy[1]/1.1; # width is over estimated
	chh=par()$cxy[2];	
	
	ch_rows=1/chh;
	ch_cols=1/chw;

	#cat("Rows: ", ch_rows, " / Cols: ", ch_cols, "\n");

	# Include cutoffs separators into table
	cutoffs=c(-1, 0.001, 0.01, 0.05, 0.10);
	cutoff_entries=cbind(0, cutoffs);
	rownames(cutoff_entries)=c(" < 0.001 :", " < 0.01 :", " < 0.05 :", " < 0.1 :", " n.s. :");
	abbr_tab_wcutoffs=rbind(abbr_tab, cutoff_entries);
	order_ix=order(abbr_tab_wcutoffs[,2]);
	abbr_tab_wcutoffs=abbr_tab_wcutoffs[order_ix,,drop=F];

	#print(abbr_tab_wcutoffs);

	# Calculate font size adjustments based on number of items in list
	num_items=nrow(abbr_tab_wcutoffs);
	vadj=1;
	# If too many items, reduce size
	if(num_items>ch_rows){
		ch_rows=num_items;
		vadj=num_items/ch_rows;
	}

	item_names=rownames(abbr_tab_wcutoffs);

	row_id=0;
	for(i in 1:num_items){

		item_len=nchar(item_names[i]);
		hadj=1;
		# If size was adjust and the item is still to long, apply another adjustment
		if(vadj*item_len*chw > 1){
			hadj=1/(vadj*item_len*chw);
		}

		estimate=abbr_tab_wcutoffs[i,"Estimate"];
		signf=abbr_tab_wcutoffs[i,2]<=0.1;
		if(estimate>0){
			draw_bullets=T;
			est_bgcol="green";
			est_sycol="black";
			est_ch="+";
			est_labfont=1;
			est_labcol=ifelse(signf, "black", "grey35");
			co_adj=1;
		}else if(estimate<0){
			draw_bullets=T;
			est_bgcol="red";
			est_sycol="white";
			est_ch="-";
			est_labfont=1;
			est_labcol=ifelse(signf, "black", "grey35");
			co_adj=1;
		}else{
			# If cutoff label
			draw_bullets=F;
			est_ch="";
			est_labfont=3;
			est_labcol="blue";
			co_adj=.7;
		}


		# Mark +/- bullet glyphs
		if(draw_bullets){
			points(coloffset+chw*.6, 1-chh*row_id, pch=21, cex=1.3*vadj, 
				col="black", bg=est_bgcol);
			points(coloffset+chw*.6, 1-chh*row_id, pch=est_ch, cex=1*vadj, 
				col=est_sycol);
		}

		# Draw response name
		text(x=coloffset+chw*.6, y=1-chh*row_id, item_names[i], cex=hadj*vadj*co_adj, 
			pos=4, font=est_labfont, col=est_labcol);

		row_id=row_id+1;
	}
	
}

#------------------------------------------------------------------------------

plot_div_summary_list=function(div_coef_pval_lst){

	par.orig=par(no.readonly=T);

	#----------------------------------------------------------------------

	par(mar=c(.5, .5, 3, .5));
	par(mfrow=c(1,1));

	item_cex=1.1;
	
	pval_mat=div_coef_pval_lst[["pval"]];
	coef_mat=div_coef_pval_lst[["coef"]];

	div_names=colnames(pval_mat);
	num_divs=length(div_names);

	plot(0, type="n", xlim=c(0,num_divs), ylim=c(0,1), bty="n",
		xaxt="n", yaxt="n", xlab="", ylab="");

	title(ylab="Responses", font.lab=3, line=-1);
	abline(v=seq(0,num_divs));
	

	div_ix=0;
	for(div_nm in div_names){

		cat("Drawing: ", div_nm, "\n", sep="");
		
		mtext("predictor", side=3, line=1.1, at=div_ix+.5, cex=.6, font=3);
		mtext(div_nm, side=3, line=.1, at=div_ix+.5, cex=item_cex, font=2);

		abtab=cbind(coef_mat[,div_nm], pval_mat[,div_nm]);
		rownames(abtab)=rownames(coef_mat);
		colnames(abtab)=c("Estimate", "Pr(>|.|)");

		print(abtab);

		draw_list(abtab, div_ix);
		
		div_ix=div_ix+1;
	}

	par(par.orig);

}

plot_div_summary_list(div_coef_pval_list);

###############################################################################
# Exports to text files
###############################################################################

###############################################################################
# Required for predictor/response downstream analysis

write_coef_pval_matrices=function(coef_mat, pval_mat, pred_names, out_rootfn){
	# Format: Factors as columns, alr (predictor) as rows

	cat("Writing coef and pvals matrices to: ", out_rootfn, " (root).\n", sep="");

	coef_mat=coef_mat[pred_names,,drop=F];
	pval_mat=pval_mat[pred_names,,drop=F];

	write.table(coef_mat, file=paste(out_rootfn, ".alr_as_pred.coefs.tsv", sep=""), 
		sep="\t", quote=F, col.names=NA, row.names=T);

	write.table(pval_mat, file=paste(out_rootfn, ".alr_as_pred.pvals.tsv", sep=""), 
		sep="\t", quote=F, col.names=NA, row.names=T);
}

write_coef_pval_matrices(t(div_coef_pval_list[["coef"]]), t(div_coef_pval_list[["pval"]]), 
	diversity_names, OutputRoot);

#------------------------------------------------------------------------------

write_alr_as_pred_summary=function(fn, div_cfpv_lst){

	fmt_cf=function(x){
		fmt=sapply(x, function(x){sprintf("%7.4f", x)});
		return(fmt)
	}

	fmt_pv=function(x){
		fmt=sapply(x, function(x){
			if(x>=0.001){
				sprintf("%9.3f", x);
			}else{
				sprintf("%9.3e", x);
			}});
		return(fmt);
	}

	# Just write out shannon and tail
	cat("Writing summary to: ", fn, "\n");

	coef_mat=div_cfpv_lst[["coef"]];
	pval_mat=div_cfpv_lst[["pval"]];

	shan_tab=cbind(coef_mat[,"Shannon"], pval_mat[,"Shannon"]);
	tail_tab=cbind(coef_mat[,"Tail"], pval_mat[,"Tail"]);

	# Order by decreasing significance
	shan_ord_ix=order(shan_tab[,2]);
	tail_ord_ix=order(tail_tab[,2]);

	shan_tab=shan_tab[shan_ord_ix,,drop=F];
	tail_tab=tail_tab[tail_ord_ix,,drop=F];

	# Format the numbers
	shan_txt_tab=cbind(fmt_cf(shan_tab[,1]), fmt_pv(shan_tab[,2]));
	tail_txt_tab=cbind(fmt_cf(tail_tab[,1]), fmt_pv(tail_tab[,2]));

	shan_txt_tab=cbind(rownames(shan_txt_tab), shan_txt_tab);
	tail_txt_tab=cbind(rownames(tail_txt_tab), tail_txt_tab);

	cn=c("Coefficients", "P-values");

	comb=rbind(
		c("Shannon", cn),
		shan_txt_tab, 
		c("","",""), 
		c("Tail", cn),
		tail_txt_tab);

	#print(comb);

	write.table(comb, fn, sep="\t", quote=F, col.names=F, row.names=F);
}

write_alr_as_pred_summary(paste(OutputRoot, ".alr_as_pred.summary.tsv", sep=""), div_coef_pval_list);

###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
