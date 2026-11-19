#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
library(tidymodels);
library(correctR);
#library(glmnet);
#library(logistf);
#library(hdi);

source("~/git/AnalysisTools/Metadata/OutputFileLibrary/OutputFileLibrary.r");

options(useFancyQuotes=F);
options(width=200);

params=c(
	"factor_fn", "f", 1, "character",
	"responses_fn", "r", 1, "character",
	"covariates_fn", "c" , 1, "character",
	"test_groups_fn", "t", 1, "character",
	"outputroot", "o", 1, "character",
	"skip_full_less_target", "R", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"\n",
	"	-f <Factor Filename>\n",
	"	-r <Response Variables List>\n",
	"	-c <Covariates List>\n",
	"	-t <Test Predictor Groups Map>\n",
	"	-o <Output Filename Root>\n",
	"\n",
	"\n", sep="");

if(!(length(opt$factor_fn) && 
	length(opt$responses_fn) &&
	length(opt$covariates_fn) && 
	length(opt$test_groups_fn) && 
	length(opt$output)
)){
	cat(usage);
	q(status=-1);
}

FactorsFile=opt$factor_fn;
ResponsesFile=opt$responses_fn;
CovariatesFile=opt$covariates_fn;
TestGroupsFile=opt$test_groups_fn;
OutputRoot=opt$outputroot;

SkipFullLessTarget=F;
if(length(opt$skip_full_less_target)){
	SkipFullLessTarget=T;
}

params=capture.output({
cat("\n");
cat("         Factors File: ", FactorsFile, "\n", sep="");
cat("      Responses File : ", ResponsesFile, "\n", sep="");
cat("      Covariates File: ", CovariatesFile, "\n", sep="");
cat("     Test Groups File: ", TestGroupsFile, "\n", sep="");
cat("          Output Root: ", OutputRoot, "\n", sep="");
cat("Skip Full Less Target: ", SkipFullLessTarget, "\n", sep="");
cat("\n");
});

print(params);
options(width=80);

##############################################################################

load_factors=function(fname){
	cat("Loading Factors/Metadata: ", fname, "\n", sep="");
	factors=data.frame(read.table(fname,  sep="\t", header=TRUE, row.names=NULL, check.names=FALSE));
	return(factors);
}

load_map=function(fname){
	# Variable name / group
	cat("Loading Map: ", fname, "\n", sep="");
	map=read.delim(fname, header=F, sep="\t", comment.char="#", stringsAsFactors=F);

	if(ncol(map)<2){
		cat("Error:  Map file needs two columns.\n");
		quit(status=-1);
	}

	map_list=list();
	grps=sort(unique(map[,2]));
	for(g in grps){
		map_list[[g]]=map[map[,2]==g, 1];
	}

	return(map_list);
}

load_list=function(filename){
	cat("Loading List: ", filename, "\n", sep="");
	val=scan(filename, what=character(), comment.char="#");
	return(val);
}

plot_title_page=function(title, subtitle=""){

        orig.par=par(no.readonly=T);
        par(family="serif");
        par(mfrow=c(1,1));

        plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n",  xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                );

        # Title
        title_cex=3;
        title_line=1;
        text(0.5, title_line, title, cex=title_cex, font=2, adj=c(.5,1));

        # Subtitle
        num_subt_lines=length(subtitle);
        cxy=par()$cxy;
        for(i in 1:num_subt_lines){
                text(.5, title_line -title_cex*cxy[2] -i*cxy[2], subtitle[i], adj=.5);
        }

        par(orig.par);
}


plot_page_separator=function(title, subtitle="", notes=""){

	par(mfrow=c(1,1));
	plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n",  xaxt="n", yaxt="n",
		xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
		);

	text(.5, .8, title, adj=c(.5, -1), cex=4, font=2); 
	text(.5, .5, subtitle, adj=c(.5, 0), cex=2, font=1); 
	text(.5, .4, notes, adj=c(.5, 1), cex=1, font=3); 
}

calc_guide_lines=function(num_cells, min_cuts=5, max_cuts=10){

	if(num_cells<max_cuts){return(c())}

	cuts=min_cuts:max_cuts;
	cut_remainders=num_cells %% cuts;
	names(cut_remainders)=cuts;
	sorted_remainders=sort(cut_remainders, method="shell");
	recommended_cut=as.numeric(names(sorted_remainders[1]));
	guides=seq(0, num_cells, recommended_cut);
	num_guides=length(guides);
	return(guides[2:(num_guides-1)]);
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

load_test_groups=function(tg_fname, include_null=F){

	test_group_map=load_map(tg_fname);

	if(include_null){
		test_group_map[["covariates_only"]]=character();
	}

	test_group_names=names(test_group_map);

	test_group_variables=c();
	for(tgrp in test_group_names){ 
		test_group_variables=c(test_group_variables, test_group_map[[tgrp]]);
	}
	num_test_groups=length(test_group_names);

	res=list();
	res[["map"]]=test_group_map;
	res[["names"]]=test_group_names;
	res[["num_grps"]]=num_test_groups;
	res[["all_var"]]=test_group_variables;

	return(res);
}

##############################################################################

cat("----------------------------------------\n");
cat("  Loading Files...\n");
cat("----------------------------------------\n");

factors_loaded=load_factors(FactorsFile);
available_variables=colnames(factors_loaded);
responses_arr=load_list(ResponsesFile);
covariates_arr=load_list(CovariatesFile);
test_group_rec=load_test_groups(TestGroupsFile, include_null=T);

cat("----------------------------------------\n");
cat("  Files Loaded...\n");
cat("----------------------------------------\n");

pdf(paste(OutputRoot, ".compare_models.pdf", sep=""), height=11, width=16);

plot_text(params);

plot_text(c(
	"Test Groups:",
	capture.output(print(test_group_rec[["map"]])),
	"",
	"Covariates:",
	capture.output(print(covariates_arr)),
	"",
	"Responses:",
	capture.output(print(responses_arr))
));


#------------------------------------------------------------------------------
# Extract variables we need, make sure they are in the factor file

all_used_variables=c(test_group_rec[["all_var"]], responses_arr, covariates_arr);

missing_variables=setdiff(all_used_variables, available_variables);
if(length(missing_variables)>0){
	cat("Error: Some variables missing from factor file:\n");
	print(missing_variables);
	quit(status=-1);
}else{
	cat("All variables found.\n");
}

factors_used=factors_loaded[,all_used_variables];

#------------------------------------------------------------------------------
# Changing the responses

change_responses_to_factors=function(factors_dataframe, rsp_arr){
	for(rsp_ix in rsp_arr){
		factors_dataframe[,rsp_ix]=
			factor(factors_dataframe[,rsp_ix], 
				levels=c(1,0),
				labels=c("1","0"));
	}
	return(factors_dataframe);
}

factors_used=change_responses_to_factors(factors_used, responses_arr);



cat("--------------------------------------------------------------------------\n");
cat("Split into Training and Test Set.\n");
#ModelingSetProp=.8;
ModelingSetProp=.99;

cat("initial_split Parameters:\n");
cat("\tprop = proportion retrained for modeling = ", ModelingSetProp, "\n");
data_split=initial_split(factors_used, prop=ModelingSetProp);

train_data=training(data_split);
test_data=testing(data_split);

num_train_samples=nrow(train_data);
num_test_samples=nrow(test_data);

cat("Number of Training Samples: ", num_train_samples, " (For CV)\n");
cat("Number of Testing Samples: ", num_test_samples, " (Held Out)\n");


cat("--------------------------------------------------------------------------\n");
cat("Split Training Set into multiple folds/partitions for Cross Validation.\n");
CVSeed=1234;
cat("Seed: ", CVSeed, "\n");
set.seed(CVSeed);

NumPartitions=10; # Number of folds
NumPartitionRepeats=5; # Number of times to repeat partitioning

#NumPartitions=3; # Number of folds
#NumPartitionRepeats=2; # Number of times to repeat partitioning

cat("vfold_cv Parameters:\n");
cat("\tv = num partitions = ", NumPartitions, "\n");
cat("\trepeats = num times to repeat partitioning = ", NumPartitionRepeats, "\n");
folds=vfold_cv(train_data, v=NumPartitions, repeats=NumPartitionRepeats);
print(folds);

cat("--------------------------------------------------------------------------\n");
cat("Building model for each response.\n");

build_workflows=function(rsp_arr, tst_gp_rec, data_mat){

	workflow_byRsp_list=list();

	# Right now, only looking at boolean responses
	lr_spec=logistic_reg() |>
		set_engine("glm") |>
		set_mode("classification");

	# For each response, build a workflow for each model
	for(rsp_ix in rsp_arr){

		cat("Building workflow for: ", rsp_ix, "\n");

		workflow_list=list();

		for(tst_gp in tst_gp_rec[["names"]]){
			cat("\tBuilding model for: ", tst_gp, "\n");

			wcov=c(covariates_arr, tst_gp_rec[["map"]][[tst_gp]]);

			rcp=recipe(
				reformulate(wcov, response=rsp_ix),
				data=factors_used) |>
				step_zv(all_predictors());

			wkfl=workflow() |>
				add_recipe(rcp) |>
				add_model(lr_spec);

			workflow_list[[tst_gp]]=wkfl

		}

		workflow_byRsp_list[[rsp_ix]]=workflow_list;
	}

	return(workflow_byRsp_list);
}


wkflw_list=build_workflows(responses_arr, test_group_rec, factors_used);

#print(wkflw_list);

cat("--------------------------------------------------------------------------\n");
cat("Fit models using the same folds\n");

fit_workflows=function(wkflw_lst, flds){

	ctrl_params=control_resamples(save_pred=T);

	ft_lst=list();
	for(rsp_ix in names(wkflw_lst)){
	
		cat("Response: ", rsp_ix, "\n");	
		gp_lst=list();
		for(grp_ix in names(wkflw_lst[[rsp_ix]])){

			cat("\tWorking on Fitting: ", grp_ix, "\n");	

			fit=fit_resamples(
				wkflw_lst[[rsp_ix]][[grp_ix]],
				resamples=flds,
				metrics=metric_set(mn_log_loss, roc_auc),
				control=ctrl_params);

			gp_lst[[grp_ix]]=fit;
		}		

		ft_lst[[rsp_ix]]=gp_lst;
	}
	
	return(ft_lst);
}

fit_list=fit_workflows(wkflw_list, folds);
print(fit_list);

cat("--------------------------------------------------------------------------\n");
cat("Summarize CV performance\n");

summarize_cv=function(ft_lst){

	cv_tab_list=list();
	for(rsp_ix in names(ft_lst)){
		cat("Response: ", rsp_ix, "\n");	

		fits=ft_lst[[rsp_ix]];

		cv_summary=fits |>
			map(collect_metrics) |>
			bind_rows(.id="model");

		cv_tab=cv_summary |>
			select(model, .metric, mean, std_err, n);

		cv_tab_list[[rsp_ix]]=cv_tab;
	}

	return(cv_tab_list);
}

cv_sumtabs_list=summarize_cv(fit_list);
print(cv_sumtabs_list);

cat("--------------------------------------------------------------------------\n");
cat("Calculate fold-level log loss\n");

parse_num=function(x){
	as.integer(gsub("[^0-9]", "", x));
}

calc_fold_logloss=function(ft_lst){

	flloss_list=list();
	for(rsp_ix in names(ft_lst)){

		cat("Response: ", rsp_ix, "\n");	
		fits=ft_lst[[rsp_ix]];

		gp_list=list();
		for(grp_ix in names(ft_lst[[rsp_ix]])){

			fit=ft_lst[[rsp_ix]][[grp_ix]];
			gp_list[[grp_ix]]=
				collect_predictions(fit) |>
				mutate(
					p_true=if(rsp_ix=="1") .pred_1 else .pred_0,
					log_loss=-log(pmax(pmin(p_true, 1-1e-15),1e-15))) |>
				group_by(id, id2) |>
				summarise(
					log_loss=mean(log_loss),
					.groups="drop") |>
				transmute(
					model=grp_ix,
					log_loss,
					k=parse_num(id2), # fold number
					r=parse_num(id))  # repeat number


		}

		flloss_list[[rsp_ix]]=gp_list;
	}

	return(flloss_list);
}

cat("Log Loss = -log(Pr(observed outcome)), so if the outcome is supposd to be 1,\n");
cat(" and the model predicts it should be 1 with a probility of .4, the log loss is -log(.4)\n");
cat(" The lower the log loss, the better the score.\n");
fold_ll_list=calc_fold_logloss(fit_list);

print(fold_ll_list);

cat("--------------------------------------------------------------------------\n");
cat("Test full vs covariates only\n");


compare_against_covariatesonly=function(fll_lst, train_n, v, r){
# For each response, compare the test group against each covariate_only

	ttest_tab_list=list();

	for(rsp_ix in names(fll_lst)){

		cat("Response: ", rsp_ix, "\n");

		grps=names(fll_lst[[rsp_ix]]);
		cat("Test Groups:\n");
		print(grps);
		comp_grp=setdiff(grps, "covariates_only");

		ttest_tab=data.frame();

		for(cmpgrp_ix in comp_grp){
			paired=inner_join(
				fll_lst[[rsp_ix]][["covariates_only"]] |>
					select(k, r, ll_reduced=log_loss),
				fll_lst[[rsp_ix]][[cmpgrp_ix]] |>
					select(k, r, ll_full=log_loss),
				by=c("k", "r"));

			#print(paired);

			delta=paired$ll_reduced - paired$ll_full;
			mean_del=mean(delta);
	
			#print(delta);

			test_data=bind_rows(
				fll_lst[[rsp_ix]][["covariates_only"]] |>
					transmute(
						model="reduced",
						values=-log_loss,
						k,r),
				fll_lst[[rsp_ix]][[cmpgrp_ix]] |>
					transmute(
						model="full",
						values=-log_loss,
						k,r)
				);

			ttest_res=correctR::repkfold_ttest(
				data=as.data.frame(test_data),
				n1=floor(train_n * (v-1)/v),
				n2=ceiling(train_n / v),
				k=v, r=r, tailed="two");
			# list of "statistic" and "p.value"

			ttest_tab=rbind(ttest_tab, 
				c(rsp_ix, cmpgrp_ix, mean_del, unlist(ttest_res), sig_char(ttest_res[["p.value"]])));

		}
		
		colnames(ttest_tab)=c("Response", "Group", "mean delta", "t-statistic", "p.value", "signf");
		ttest_tab_list[[rsp_ix]]=ttest_tab;
	
	}

	return(ttest_tab_list);

}

ttest_list=compare_against_covariatesonly(fold_ll_list, num_train_samples, v=NumPartitions, r=NumPartitionRepeats);

print(ttest_list);

# Rbind all the matrices in the list together

ttest_matrix=do.call(rbind, ttest_list);

print(ttest_matrix);

ttest_matrix[,"mean delta"]=round(as.numeric(ttest_matrix[,"mean delta"]),5);
ttest_matrix[,"t-statistic"]=round(as.numeric(ttest_matrix[,"t-statistic"]),3);
ttest_matrix[,"p.value"]=signif(as.numeric(ttest_matrix[,"p.value"]),5);

cat("\n\n");
cat("By Response and Group:\n");
print(ttest_matrix);

cat("\n\n");
cat("By Significance:\n");
by_signf=order(ttest_matrix[,"p.value"], decreasing=F);

ttest_matrix_by_signf=ttest_matrix[by_signf,];
print(ttest_matrix_by_signf);

###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
