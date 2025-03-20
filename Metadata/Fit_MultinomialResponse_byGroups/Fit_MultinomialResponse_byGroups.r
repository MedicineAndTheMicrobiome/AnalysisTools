#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

source("~/git/AnalysisTools/Metadata/OutputFileLibrary/OutputFileLibrary.r");

options(useFancyQuotes=F);

params=c(
	"factor_fn", "f", 1, "character",
	"responses_fn", "r", 1, "character",
	"covariates_fn", "c" , 2, "character",
	"test_groups_fn", "t", 2, "character",
	"outputroot", "o", 1, "character"
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
	"This script will fit the following multinomial logistic regression model:\n",
	"\n",
	"	[response] = [covariates] + [test group1] + [test group2] + [test groupn]\n",
	"\n",
	"\n",
	"The Factor file should contain all the variables for all the subjects.\n",
	"The Response variables is a list of all the variables (columns) that should be\n",
	"   independently analyzed.\n",
	"The Covariates variables are always included in the model as predictors.\n",
	"The Test Predictor Groups Map should contain two colums\n",
	"   <variable name>\\t<group name>\n",
	"\n",
	"A full model will be fit with all the variables, then a reduced model excluding the\n",
	"   the test group will then be fit.  Since the response is multinomial, a model will\n",
	"   be fit for each category.\n",
	"\n",
	"This will be done for each of the multinomial response variables.\n",
	"A separate pdf file result will be generated for each multinomial response.\n",
	"\n",
	"This script was designed to be run on the results of a clustering algorithm, where each\n",
	"   cut is a separate multinomial response.\n",
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

params=capture.output({
cat("\n");
cat("    Factors File: ", FactorsFile, "\n", sep="");
cat(" Responses File : ", ResponsesFile, "\n", sep="");
cat(" Covariates File: ", CovariatesFile, "\n", sep="");
cat("Test Groups File: ", TestGroupsFile, "\n", sep="");
cat("     Output Root: ", OutputRoot, "\n", sep="");
cat("\n");
});

print(params);
options(width=80);

##############################################################################

load_factors=function(fname){
	cat("Loading Factors/Metadata: ", fname, "\n", sep="");
	factors=data.frame(read.table(fname,  sep="\t", header=TRUE, row.names=1, check.names=FALSE));
	return(factors);
}

load_map=function(fname){
	# Variable name / group
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

plot_text=function(strings, max_lines_pp=Inf, oma_tag=""){

        orig.par=par(no.readonly=T);

        par(mfrow=c(1,1));
        par(family="Courier");

	if(oma_tag==""){
		par(oma=rep(.5,4));
	}else{
		par(oma=c(1, .5, .5, .5));
	}

        par(mar=rep(0,4));

        num_lines=length(strings);
        num_pages=max(1, ceiling(num_lines/max_lines_pp));

        cat("Num Pages for ", num_lines, " lines: ", num_pages, "\n", sep="");

        lines_pp=min(num_lines, max_lines_pp);
        for(p in 1:num_pages){

                top=max(as.integer(lines_pp), 52);
	
                plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
                        xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                        );

		if(oma_tag!=""){
			mtext(paste("[", oma_tag, "]", sep=""), side=1, col="grey25");
		}

                text_size=max(.01, min(.7, .7 - .003*(lines_pp-52)));
                #print(text_size);

                start=(p-1)*lines_pp+1;
                end=start+lines_pp-1;
                end=min(end, num_lines);
                line=1;
                for(i in start:end){
                        #cat(strings[i], "\n", sep="");
                        strings[i]=gsub("\t", "", strings[i]);
                        text(0, top-line, strings[i], pos=4, cex=text_size);
                        line=line+1;
                }

        }

        par(orig.par);
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
##############################################################################


factors_loaded=load_factors(FactorsFile);
available_variables=colnames(factors_loaded);

test_group_map=load_map(TestGroupsFile);
test_group_names=names(test_group_map);
test_group_variables=c();
for(tgrp in test_group_names){ 
	test_group_variables=c(test_group_variables, test_group_map[[tgrp]]);
}
num_test_groups=length(test_group_names);

responses_arr=load_list(ResponsesFile);
covariates_arr=load_list(CovariatesFile);


pdf(paste(OutputRoot, ".multn_resp.pdf", sep=""), height=11, width=12);

plot_text(params);

plot_text(c(
	"Test Groups:",
	capture.output(print(test_group_map)),
	"",
	"Covariates:",
	capture.output(print(covariates_arr)),
	"",
	"Responses:",
	capture.output(print(responses_arr))
));



all_used_variables=c(test_group_variables, responses_arr, covariates_arr);

missing_variables=setdiff(all_used_variables, available_variables);
if(length(missing_variables)>0){
	cat("Error: Some variables missing from factor file:\n");
	print(missing_variables);
	quit(status=-1);
}else{
	cat("All variables found.\n");
}

factors_used=factors_loaded[,all_used_variables];

###############################################################################

covariates_formula_string=paste(covariates_arr, collapse=" + ");

test_group_strings_list=character(num_test_groups);
names(test_group_strings_list)=test_group_names;
for(tg in test_group_names){
	test_group_strings_list[tg]=paste(test_group_map[[tg]], collapse=" + ");
}

cat("--------------------------------------------------------------------------\n");

cat("Covariates Components:\n");
print(covariates_formula_string);

cat("Test Group Components:\n");
print(test_group_strings_list);

cat("--------------------------------------------------------------------------\n");
cat("--------------------------------------------------------------------------\n");

cat("Models:\n\n");

full_model=paste(covariates_formula_string,
	paste(test_group_strings_list, collapse=" + "), 
	sep=" + ");

covariates_model=paste(covariates_formula_string);

targeted_reduced_model_list=character(num_test_groups);
names(targeted_reduced_model_list)=test_group_names;
for(tg  in test_group_names){
	reduced_list=setdiff(test_group_names, tg);
	targeted_reduced_model_list[tg]=paste(covariates_formula_string,
		paste(test_group_strings_list[reduced_list], collapse=" + "), 
		sep=" + ");
}

targets_cov_model_list=paste(covariates_formula_string, test_group_strings_list, sep=" + ");
names(targets_cov_model_list)=test_group_names;


cat("Full Model:\n");
cat(full_model, "\n");
cat("\n\n");
cat("Targeted Reduced Model:\n");
cat(paste(targeted_reduced_model_list, collapse="\n\n"));
cat("\n\n");
cat("Covariates Model:\n");
cat(covariates_model, "\n");
cat("\n\n");
cat("Targets Model:\n");
cat(paste(targets_cov_model_list, collapse="\n\n"));
cat("\n\n");

cat("--------------------------------------------------------------------------\n");
cat("--------------------------------------------------------------------------\n");
cat("\n\n");

mm=model.matrix(as.formula(paste(responses_arr[1], " ~ ", full_model)), factors_loaded);
expected_variables=setdiff(colnames(mm), "(Intercept)");
num_expected_variables=length(expected_variables);
cat("Expected Variables:\n");
print(expected_variables);
cat("\n");

###############################################################################

plot_title_page("Multinomial Response Analysis", c(
	"In this analysis, each of the columns specified as a targeted response column",
	"will be modeled as a response against the specified covariates and test predictor",
	"groups.",
	"",
	"In all models, the covariates will be included, but specific combinations of",
	"the test predictor groups will be included or excluded, in order to identify whether the",
	"group as a whole contributes to the prediction of the response.",
	"",
	"The following models will be computed:",
	"Full: Covariates + All Test Predictor Groups",
	"Covariates: Covariates Only",
	"Covariates+[Group]: Covariates and a single test predictor group",
	"Full-[Group]: Full model except a single test predictor group",
	"",
	"(If only 2 test groups are supplied, then the Full-[Group] model is not fit since it is",
	"redundant to the Covariates+[Group] models.)",
	"",
	"For each of the specified response columns, a battery of analyses will be generated.",
	"Procedurally, each of the categories in the multivariate response will be fit independently",
	"(in category vs remaining) with a generalized linear model (GLM) using the binomial",
	"family as the error distribution and link function." 
));

plot_title_page("Associations From Model Fitting", c(
	"The following plots will be generated for each targeted response column:",
	"",
	"1.) Full Model Heatmap (All Associations):",
	"The coefficients for each of the predictors (rows) and for the categories (columns)",
	"are represented in the heat map.  Cells are colored red and blue for positive and",
	"negative associations, respectively.",
	"",
	"2.) Full Model Significant Associations:",
	"Similar to the previous figure, but only the associations with p-values below the",
	"specified cutoff are displayed.",
	"",
	"3.) Formatted Text List of Significant Assocations by Response Category",
	"",
	"4.) Predictor Associations by Response Category",
	"In these barplots, the number of significant associations between each response category",
	"and the color coded predictor groups are displyed.  Response categories are sorted with the",
	"greatest number of significant predictor associations first.  The dashed blue line",
	"represents the maximum number of predictors that can be associated with each category.",
	"",
	"5.) Response Categories by Predictor", 
	"In these barplots, the number of response categories each of the predictors is",
	"significantly associated with is plotted.  The predictor variables are sorted in",
	"decreasing order by the number of response categories each is associated with.",
	"The dashed blue line represents the maximum number of response categories each predictor",
	"can be associated with."
));	

plot_title_page("Model Fitness Evaluations", c(
	"After the figures for the model fitting, an analysis of model fit using the AIC",
	"statistic is performed.  The greater the AIC, the better the model.  If two models are",
	"greater than 2 AIC values apart, then the model with the greater AIC can be regarded as",
	"signficantly better than the lower one.",
	"",
	"1.) Model Fit AIC Values:",
	"For all the models fit and each of the response categories, the AIC values are",
	"reported in a table.",
	"",
	"2.) Relative AIC: All Models",
	"For each model and response category, the models are plotted by their AIC values.",
	"so their relative fitness within and between the response categories can be visualized.",
	"",
	"3.) Top AICs: Only Models within 2 from Top Model",
	"Only the models that are within 2 from the best fitting model for each category",
	"are plotted, in order to reduce some clutter.",
	"",
	"4.) Top Models Table",
	"In this table, the top models for each response category are colored blue.",
	"If there is a single best model, then the color will be dark blue.  If there are",
	"multiple best models then all competing best models will share a lighter blue.",
	"The model with the numerically greatest AIC value will be demarked with an",
	"orange asterisk."
));

###############################################################################

accumulate_sumfit_to_matrix=function(sumfit, mat, cname){

	coef_tab=sumfit[["coefficients"]]	
	varnames=setdiff(rownames(coef_tab), "(Intercept)"); 
	
	mat[["coef"]][varnames, cname]=coef_tab[varnames, "Estimate"];
	mat[["pval"]][varnames, cname]=coef_tab[varnames, "Pr(>|z|)"];
	mat[["aic"]][cname]=sumfit[["aic"]];

	return(mat);
}

responses_results=list();

for(resp_ix in responses_arr){

	cat("Working on Multinomial Response: ", resp_ix, "\n");

	resp_categories_values=as.character(factors_loaded[,resp_ix]);
	resp_categories=sort(unique(resp_categories_values));
	cat("Response Categories: \n");
	print(resp_categories);
	num_resp_categories=length(resp_categories);
	cat("\n");

	# pattern the matrix depending on number of response categories
	empty_matrix=matrix(NA, nrow=num_expected_variables, ncol=num_resp_categories);
	colnames(empty_matrix)=resp_categories;
	rownames(empty_matrix)=expected_variables;

	empty_arr=numeric(num_resp_categories);
	names(empty_arr)=resp_categories;

	empty_pval_coef=list();
	empty_pval_coef[["coef"]]=empty_matrix;
	empty_pval_coef[["pval"]]=empty_matrix;
	empty_pval_coef[["aic"]]=empty_arr;;
	

	responses_results[[resp_ix]]=list();

	responses_results[[resp_ix]][["full"]]=empty_pval_coef;
	responses_results[[resp_ix]][["covariates"]]=empty_pval_coef;
	responses_results[[resp_ix]][["target_reduced"]]=list()
	responses_results[[resp_ix]][["target_covar"]]=list();

	for(tg in test_group_names){
		if(num_test_groups>2){
			responses_results[[resp_ix]][["target_reduced"]][[tg]]=empty_pval_coef;
		}
		responses_results[[resp_ix]][["target_covar"]][[tg]]=empty_pval_coef;
	}

	#print(responses_results);

	for(category_ix in resp_categories){

		cat("Fitting Model: ", category_ix, " vs. Other.\n");
		cur_responses=factors_used[,resp_ix];
		responses=(category_ix==cur_responses);

		tmp_factors=cbind(responses, factors_used);

		# Test Full (All Targets + Covariates)
		cat("\tFitting Full Model...\n");
		full_formula=as.formula(paste("responses ~ ", full_model));
		full_fit=glm(full_formula, data=as.data.frame(tmp_factors), family="binomial");
		full_sum_fit=summary(full_fit);
		#print(full_sum_fit);

		responses_results[[resp_ix]][["full"]]=
			accumulate_sumfit_to_matrix(full_sum_fit, 
				responses_results[[resp_ix]][["full"]], category_ix);


		# Test Covariates only
		cat("\tFitting Covariates-only Model...\n");
		cov_formula=as.formula(paste("responses ~ ", covariates_model));
		cov_fit=glm(cov_formula, data=as.data.frame(tmp_factors), family="binomial");
		cov_sum_fit=summary(cov_fit);
		#print(cov_sum_fit);

		responses_results[[resp_ix]][["covariates"]]=
			accumulate_sumfit_to_matrix(cov_sum_fit, 
				responses_results[[resp_ix]][["covariates"]], category_ix);
	
		for(tg in test_group_names){
	
			cat("\tFitting Targeted Models (", tg, "):\n");

			if(num_test_groups>2){

				# Test Targeted Reduced (All Targets except Target, + Covariates)
				cat("\t\tFitting Targeted Reduced Model...\n");
				tar_red_formula=as.formula(paste("responses ~ ", targeted_reduced_model_list[tg]));
				tar_red_fit=glm(tar_red_formula, data=as.data.frame(tmp_factors), family="binomial");
				tar_red_sum_fit=summary(tar_red_fit);
				#print(tar_red_sum_fit);

				responses_results[[resp_ix]][["target_reduced"]][[tg]]=
					accumulate_sumfit_to_matrix(tar_red_sum_fit, 
						responses_results[[resp_ix]][["target_reduced"]][[tg]], category_ix);
			}


			# Test Target (Target + Covariates)
			cat("\t\tFitting Targeted+Covariates Model...\n");
			tar_formula=as.formula(paste("responses ~ ", targets_cov_model_list[tg]));
			tar_fit=glm(tar_formula, data=as.data.frame(tmp_factors), family="binomial");
			tar_sum_fit=summary(tar_fit);
			#print(tar_sum_fit);

			responses_results[[resp_ix]][["target_covar"]][[tg]]=
				accumulate_sumfit_to_matrix(tar_sum_fit, 
					responses_results[[resp_ix]][["target_covar"]][[tg]], category_ix);


		}
	}

	cat("\n");
}

###############################################################################

mask_by_cutoff=function(coef_mat, pval_mat, pval_cutoff=1, mask_val=0){
        masked_mat=coef_mat;
        masked_mat[pval_mat>pval_cutoff]=mask_val;
	return(masked_mat);
}

plot_signif_summary=function(coef_mat, pval_mat){

	print(coef_mat);

	dir_mat=apply(coef_mat, c(1,2), function(x){ 
		ifelse(x==0, 0, ifelse(x>0, 1, -1));}
		);

	signf_mat=apply(pval_mat, c(1,2), function(x){
		if(x<.001){
			s=4;
		}else if(x<0.01){
			s=3;
		}else if(x<0.05){
			s=2;
		}else if(x<0.10){
			s=1;
		}else{
			s=0;
		}
	});

	comp_mat=dir_mat * signf_mat;	
		
	paint_matrix(comp_mat, title="Significant Associations", 
		deci_pts=2, label_zeros=F, value.cex=1);
	
}


get_aic=function(fit_results){
	
	# Pull AIC values from results
	model_aics=list();
	model_aics[["full"]]=fit_results[["full"]]$aic;
	model_aics[["covar_only"]]=fit_results[["covariates"]]$aic;

	tr_names=names(fit_results[["target_reduced"]]);
	for(tr in tr_names){
		trname=paste("full-", tr, sep="");
		model_aics[[trname]]=fit_results[["target_reduced"]][[tr]]$aic;
	}

	tr_names=names(fit_results[["target_covar"]]);
	for(tr in tr_names){
		trname=paste("covar+", tr, sep="");
		model_aics[[trname]]=fit_results[["target_covar"]][[tr]]$aic;
	}

	# Combine into matrix
	results=numeric();
	mnames=names(model_aics);
	for(m in mnames){
		results=rbind(results, model_aics[[m]]);
	}
	rownames(results)=mnames;
	colnames(results)=names(model_aics[["full"]]);

	return(results);

}

#------------------------------------------------------------------------------

adj_spacing=function(aicm){
	param=par()
	ch=param$cxy[2];	
	min_pos=param$usr[3]+.5*ch;
	max_pos=param$usr[4]-.5*ch;
	num_labels=nrow(aicm)*ncol(aicm);

	aic_arr=as.numeric(aicm);
	names(aic_arr)=1:num_labels;
	aic_arr=sort(aic_arr);
	num_aics=length(aic_arr);

	ch=min(ch, (max_pos-min_pos)/(num_aics+2));
	adj=.5*ch;
	cat("Targeted Spacings: ", ch, " (Ideal: ", param$cxy[2], ")\n");
	cat("Allowed Adjustment Range: ", min_pos, " - ", max_pos, "\n");

	adj_made=T;
	max_iter=1000000;

	tweak_iter=0
	iter=0;
	
	while(adj_made && iter<max_iter){
		adj_made=F;

		for(i in 1:(num_aics-1)){
			if((aic_arr[i+1]-aic_arr[i])<ch){
				aic_arr[i+1]=min(aic_arr[i+1]+adj, max_pos);
				adj_made=T;
			}
		}

		for(i in (num_aics:2)){
			if((aic_arr[i]-aic_arr[i-1])<ch){
				aic_arr[i-1]=max(aic_arr[i-1]-adj, min_pos);
				adj_made=T;
			}
		}

		if(tweak_iter>2*num_labels){
			ch=ch*.95;	
			adj=adj*.95;
			tweak_iter=0;
			cat("Reducing spacing requirements: ", ch, "\n");
		}

		iter=iter+1;
		tweak_iter=tweak_iter+1;
	}

	if(iter==max_iter){
		cat("WARNING: Hit adjustment iteration limit.\n");
	}

	reorder=order(as.numeric(names(aic_arr)));
	aic_arr=aic_arr[reorder];

	na_padded=rep(NA, num_labels);
	names(na_padded)=1:num_labels;
	na_padded[as.numeric(names(aic_arr))]=aic_arr;	

	out_aicm=matrix(na_padded, nrow=nrow(aicm));
	colnames(out_aicm)=colnames(aicm);
	rownames(out_aicm)=rownames(aicm);

	return(out_aicm);

}

#------------------------------------------------------------------------------

plot_top_aics_table=function(aic_matrix, aic_thres=2, title=""){
	# Colored Cells
	# The "best" AIC is the lowest.

	cat("Generating table for top AICs Models...\n");

	num_col=ncol(aic_matrix);
	num_models=nrow(aic_matrix);
	model_names=rownames(aic_matrix);
	category_names=colnames(aic_matrix);

	# Remove models not within aic_thres of top
	min_aics=apply(aic_matrix, 2, min);
	for(i in 1:num_col){
		val=aic_matrix[,i];
		val[(val-min_aics[i])>aic_thres]=NA;
		aic_matrix[,i]=val;
	}


	par(mar=c(1,26,14,5));
	par(oma=c(0,0,2,0));

	plot(0,0, type="n", xlim=c(0, num_col), ylim=c(0, num_models),
		ylab="", xlab="Response Categories", bty="n",
		xaxt="n", yaxt="n"
		);

	text((1:num_col)-1, num_models+.05, pos=4, srt=45, xpd=T, cex=.8, category_names);

	text(0, (1:num_models)-.5, pos=2, srt=0, xpd=T, model_names);

	num_best=apply(aic_matrix, 2, function(x){ sum(!is.na(x))});

	for(cix in 1:num_col){
		for(mix in 1:num_models){

			best=F;

			if(!is.na(aic_matrix[mix, cix])){

				fill=ifelse(num_best[cix]==1, "blue", "cornflowerblue");

				otherval=aic_matrix[,cix];
				otherval=otherval[!is.na(otherval)];
				best=T;
				if(any(aic_matrix[mix, cix]>otherval)){
					best=F;
				}

			}else{
				fill="white";
			}

			# Fill in cell
			rect(xleft=cix-1, ybottom=mix-1, xright=cix, ytop=mix, col=fill);

			# Mark as best
			if(best){
				points(cix-.5, mix-.5, type="p", col="darkorange3", 
					font=2, cex=3, pch="*")
			}

		}	
	}

	mtext(title, side=3, outer=T, font=2, cex=1.5);
	mtext("(Lower AIC is Better Model)", side=3, outer=T, line=-1, font=1, cex=1);

}

#------------------------------------------------------------------------------

plot_top_aics=function(aic_matrix, aic_thres=2, title=""){

	cat("Plotting top AICs Models...\n");
	#print(aic_matrix);

	num_col=ncol(aic_matrix);
	num_models=nrow(aic_matrix);
	model_names=rownames(aic_matrix);
	category_names=colnames(aic_matrix);

	# Remove models not within aic_thres of top
	min_aics=apply(aic_matrix, 2, min);
	for(i in 1:num_col){
		val=aic_matrix[,i];
		val[(val-min_aics[i])>aic_thres]=NA;
		aic_matrix[,i]=val;
	}

	aic_range=range(as.numeric(aic_matrix), na.rm=T);
	min_aic=aic_range[1];
	max_aic=aic_range[2];
	cat("Input Range: ", aic_range, "\n");

	aic_range=max_aic-min_aic;

	par(mar=c(1,4,14,14));
	par(oma=c(0,0,2,0));

	plot(0,0, type="n", xlim=c(0, num_col), ylim=c(min_aic, max_aic),
		ylab="AIC", xlab="Response Categories", bty="n",
		xaxt="n"
		);
	plot_top=par()$usr[4];
	plot_bottom=par()$usr[3];
	text((1:num_col)-1, plot_top, pos=4, srt=45, xpd=T, font=2, cex=.8, category_names);
	abline(v=0:num_col, col="grey");

	aic_matrix_spc_adj=adj_spacing(aic_matrix);

	chw=par()$cxy[1];

	for(cix in 1:num_col){
		for(mix in 1:num_models){
			points(c(cix-.95, cix-1.05), rep(aic_matrix[mix, cix], 2), 
				type="l", col="blue");
			points(c(cix-1+.05, cix-1+.1), 
				c(aic_matrix[mix, cix], aic_matrix_spc_adj[mix, cix]), 
				type="l", col="blue");

			text(cix-1+.1-chw*.70, aic_matrix_spc_adj[mix, cix], model_names[mix], 
				pos=4, xpd=T);

		}	
	}

	mtext(title, side=3, outer=T, font=2, cex=1.5);
	mtext("(Lower AIC is Better Model)", side=3, outer=T, line=-1, font=1, cex=1);

}

#------------------------------------------------------------------------------

list_signf_pred_by_category=function(coef_mat, pval_mat, pval_cutoff=.1, oma_tag=""){

	category_names=colnames(coef_mat);
	num_categories=ncol(coef_mat);
	predictor_names=rownames(coef_mat);

	orig_width=options()$width;

	out_text=c();
	options(width=1000);
	for(cat in category_names){

		# Extract signf rows for this category
		signif=pval_mat[,cat]<=pval_cutoff;
		signf_coef=coef_mat[signif, cat];
		signf_pval=pval_mat[signif, cat];
		pred_names=predictor_names[signif];

		names(signf_coef)=pred_names;
		names(signf_pval)=pred_names;

		# Sort by pvalue
		sort_pval_order=order(signf_pval, decreasing=F);
		cat_pval_sort=signf_pval[sort_pval_order];
		cat_coef_sort=signf_coef[sort_pval_order];
		pred_names=pred_names[sort_pval_order];

		num_signf=length(signf_coef);
	
		pos_assoc=character();
		neg_assoc=character();
		
		if(length(pred_names)){
			for(pred in pred_names){
				formatted=paste(
					pred, 
					" (coef=", sprintf("%3.4f", cat_coef_sort[pred]), ", ",
					"p-val=", sprintf("%3.4f", cat_pval_sort[pred]), ")", 
					sep="");

				if(!is.na(cat_coef_sort[pred])){
					if(cat_coef_sort[pred]>0){
						pos_assoc=c(pos_assoc, formatted);
					}else{
						neg_assoc=c(neg_assoc, formatted);
					}
				}
			}
		}	

		num_pos=length(pos_assoc);
		num_neg=length(neg_assoc);

		max_len=max(num_pos, num_neg);

		pos_neg_mat=matrix("", nrow=max_len, ncol=2);
		colnames(pos_neg_mat)=c("[Negative (-)]", "[Positive (+)]");

		if(num_neg>0){
			pos_neg_mat[1:num_neg,1]=neg_assoc;
		}
		if(num_pos>0){
			pos_neg_mat[1:num_pos,2]=pos_assoc;
		}

		out_text=c(out_text,
			paste("Response Category: ", cat, sep=""),
			capture.output(print(pos_neg_mat, quote=F, width=1000)),
			""
			);

	}	

	#print(out_text);
	plot_text(c(
		paste("Number of Categories: ", num_categories),
		"",
		out_text,
		"----------------------------------------------End----------------------------------------------"
	), max_lines_pp=50, oma_tag=oma_tag);

	options(width=orig_width);

}

list_signf_category_by_pred=function(coef_mat, pval_mat, pval_cutoff=.1, oma_tag=""){

	category_names=colnames(coef_mat);
	num_categories=ncol(coef_mat);
	num_predictors=nrow(coef_mat);
	predictor_names=rownames(coef_mat);

	orig_width=options()$width;

	out_text=c();
	options(width=1000);

	for(pred in predictor_names){

		# Extract signf rows for this category
		signif=pval_mat[pred,]<=pval_cutoff;
		signf_coef=coef_mat[pred, signif];
		signf_pval=pval_mat[pred, signif];
		cat_names=category_names[signif];

		names(signf_coef)=cat_names;
		names(signf_pval)=cat_names;

		# Sort by pvalue
		sort_pval_order=order(signf_pval, decreasing=F);
		pred_pval_sort=signf_pval[sort_pval_order];
		pred_coef_sort=signf_coef[sort_pval_order];
		cat_names=cat_names[sort_pval_order];


		num_signf=length(signf_coef);
	
		pos_assoc=character();
		neg_assoc=character();

		# Split the pos and neg associations, so we can make a two column
		# report
		if(length(cat_names)){		
			for(cat in cat_names){
				formatted=paste(
					cat, 
					" (coef=", sprintf("%3.4f", pred_coef_sort[cat]), ", ",
					"p-val=", sprintf("%3.4f", pred_pval_sort[cat]), ")", 
					sep="");

				if(!is.na(pred_coef_sort[cat])){
					if(pred_coef_sort[cat]>0){
						pos_assoc=c(pos_assoc, formatted);
					}else{
						neg_assoc=c(neg_assoc, formatted);
					}
				}
			}	
		}

		num_pos=length(pos_assoc);
		num_neg=length(neg_assoc);

		max_len=max(num_pos, num_neg);

		pos_neg_mat=matrix("", nrow=max_len, ncol=2);
		colnames(pos_neg_mat)=c("[Negative (-)]", "[Positive (+)]");

		if(num_neg>0){
			pos_neg_mat[1:num_neg,1]=neg_assoc;
		}
		if(num_pos>0){
			pos_neg_mat[1:num_pos,2]=pos_assoc;
		}

		out_text=c(out_text,
			paste("Predictor: ", pred, sep=""),
			capture.output(print(pos_neg_mat, quote=F, width=1000)),
			""
			);

	}	

	#print(out_text);
	plot_text(c(
		paste("Number of Predictors: ", num_predictors),
		"",
		out_text,
		"----------------------------------------------End----------------------------------------------"
	), max_lines_pp=50, oma_tag=oma_tag);

	options(width=orig_width);

}

#------------------------------------------------------------------------------

plot_contributors=function(pval_mat, pval_cutoff, covariates, group_map, resp_name){

	group_map[["covariates"]]=covariates;

	num_predictors=nrow(pval_mat);
	num_categories=ncol(pval_mat);

	num_groups=length(group_map);
	grp_names=names(group_map);

	num_resp_cat=ncol(pval_mat);
	resp_names=colnames(pval_mat);	

	#----------------------------------------------------------------------
	# Get group colors
	grp_col=rainbow(num_groups, end=2/3);
	names(grp_col)=grp_names;

	# Assign group colors to underlying variable names
	pred_cols=character();
	for(g_ix in grp_names){
		grp_var=rep(grp_col[[g_ix]], length(group_map[[g_ix]]));
		names(grp_var)=group_map[[g_ix]];
		pred_cols=c(pred_cols, grp_var);
	}
	print(pred_cols);

	#----------------------------------------------------------------------
	#----------------------------------------------------------------------
	# Split predictors into their groups

	contrib_matrix=matrix(NA, nrow=num_groups, ncol=num_resp_cat);
	rownames(contrib_matrix)=grp_names;
	colnames(contrib_matrix)=resp_names;

	for(g_ix in grp_names){
		gr_var_list=group_map[[g_ix]];
		grp_pval_mat=pval_mat[gr_var_list,,drop=F];
		num_predict_contrib=apply(grp_pval_mat, 2, function(x){
			x=x[!is.na(x)];
			sum(x<=pval_cutoff)});
		contrib_matrix[g_ix,]=num_predict_contrib;
	}

	# Sort by total contributions
	all_contrib_arr=apply(contrib_matrix, 2, sum);

	sort_ix=order(all_contrib_arr, decreasing=T);
	contrib_matrix=contrib_matrix[,sort_ix, drop=F];
	all_contrib_arr=all_contrib_arr[sort_ix];

	cat("\nContributions by Group:\n");
	print(contrib_matrix);
	cat("\n\nTotal Contributions across Groups:\n");
	print(all_contrib_arr);
	cat("\n");

	#----------------------------------------------------------------------
	# Generate Predictors by Category plot

	par(mfrow=c(2,1));
	par(mar=c(15, 4, 2, 10));
	par(oma=c(1,0,0,0));

	plot(0,0, xlim=c(0, num_resp_cat+1), ylim=c(0, num_predictors*1.1),
		main=paste("Num. of Assoc. Predictors by Resp. Category: p-val < ", pval_cutoff, sep=""),
		bty="n",
		xaxt="n", xlab="", ylab="Num Predictors", type="n");

	cumsums=apply(contrib_matrix, 2, cumsum);
	print(cumsums);
	for(g_ix in rev(grp_names)){
		for(cat_ix in 1:num_categories){
			rect(
				xleft=cat_ix-1+.05, 
				xright=cat_ix-.05, 
				ybottom=0, 
				ytop=cumsums[g_ix, cat_ix], 
			col=grp_col[g_ix]);
		}
	}

	bmids=(1:num_resp_cat)-.5;

	abline(h=num_predictors, lty="dashed", col="blue");
	chx=par()$cxy[1];
	chy=par()$cxy[2]
	text(bmids-chx, -.75*chy, names(all_contrib_arr), pos=4, srt=-45, cex=.7, xpd=T);

	legend(num_resp_cat*3/4, num_predictors, legend=grp_names, fill=grp_col, bty="n");

	#----------------------------------------------------------------------
	# Generate Categories by Predictor plot

	pred_contrib=apply(pval_mat, 1, function(x){
		x=x[!is.na(x)];
		sum(x<=pval_cutoff)});

	pred_contrib=sort(pred_contrib, decreasing=T);

	bmids=barplot(pred_contrib, 
		main=paste("Num. of Assoc. Categories by Predictor: p-val < ", pval_cutoff, sep=""),
		xaxt="n", ylim=c(0, num_categories*1.1), ylab="Num Categories",
		col=pred_cols[names(pred_contrib)]
		);
	abline(h=num_categories, lty="dashed", col="blue");
	chx=par()$cxy[1];
	chy=par()$cxy[2]

	labelcex=min(1, 45/num_predictors);
	text(bmids-chx, -.75*chy, names(pred_contrib), pos=4, srt=-45, cex=labelcex, xpd=T);

	legend(num_predictors*3/4, num_resp_cat, legend=grp_names, fill=grp_col, bty="n");

	mtext(text=paste("[", resp_name, "]", sep=""), side=1, family="Courier", col="grey25", outer=T);

}

#------------------------------------------------------------------------------

plot_results=function(resp_rec){
	
	cat("Plotting Results:\n");
	resp_var_names=names(resp_rec);
	num_responses=length(resp_var_names);

	cat("Number of Response Variables:", num_responses, "\n");

	plot_text(c(
		paste("Number of Respose Variables:", num_responses),
		"",
		print(resp_var_names)
	));

	for(resp_var_ix in resp_var_names){

		cat("Working on Response Variable: ", resp_var_ix, "\n");
		plot_page_separator(resp_var_ix);

		resp_var_res=resp_rec[[resp_var_ix]];
		models=names(resp_var_res);
		#print(models);

		#print(resp_var_res[["full"]]);
		full_coef_mat=resp_var_res[["full"]][["coef"]];
		full_pval_mat=resp_var_res[["full"]][["pval"]];

		resp_guide_lines=calc_guide_lines(ncol(full_coef_mat), min_cuts=5, max_cuts=10);
		pred_guide_lines=calc_guide_lines(nrow(full_coef_mat), min_cuts=5, max_cuts=10);

		masked_coef=mask_by_cutoff(full_coef_mat, full_pval_mat, pval_cutoff=1);
		paint_matrix(masked_coef, 
			title=paste(resp_var_ix, ":  Full Model, All Assoc", sep=""), 
			deci_pts=2, label_zeros=F, value.cex=1);

		masked_coef=mask_by_cutoff(full_coef_mat, full_pval_mat, pval_cutoff=.1);
		paint_matrix(masked_coef, 
			title=paste(resp_var_ix, ":  Full Model, Signf Assoc (p-val<0.1)", sep=""), 
			deci_pts=2, label_zeros=F, value.cex=1);

		list_signf_pred_by_category(full_coef_mat, full_pval_mat, pval_cutoff=.1, 
			oma_tag=resp_var_ix);
		list_signf_category_by_pred(full_coef_mat, full_pval_mat, pval_cutoff=.1, 
			oma_tag=resp_var_ix);

		plot_contributors(full_pval_mat, pval_cutoff=.1, covariates_arr, 
			test_group_map, resp_var_ix);
	
		aic_values=get_aic(resp_var_res);	

		par(mfrow=c(1,1));		

		options(width=120);
		plot_text(c(
			paste("Response Name: ", resp_var_ix),
			"Model Fit AIC Values:",
			"(Remember: If two models are within 2 AIC units of each other,",
			"  then they can not be considered significantly better/worse than",
			"  each other.)", 
			"",
			capture.output(print(aic_values))
		), oma_tag=resp_var_ix);

		plot_top_aics(aic_values, Inf, 
			paste(resp_var_ix, ": Relative AIC, All Models", sep=""));
		plot_top_aics(aic_values, 2, 
			paste(resp_var_ix, ": Models with AIC within 2 of Top Model", sep=""));

		plot_top_aics_table(aic_values, 2, paste(resp_var_ix, ": Top Models Table", sep=""));

	}	
}

#------------------------------------------------------------------------------

plot_results(responses_results);

#------------------------------------------------------------------------------

par(mfrow=c(1,1));
par(mar=c(1,1,1,1));

plot_title_page("Model Improvement through Group Inclusion Measured by Num Factor Levels", c(
	"The following matrix:",
	"    Rows: Group Variable Names:",
	"    Cols: Response Variable Names",
	"",
	"Each cell represents the number of models where including the",
	"Group's variables, improved the model for the factor level:",
	"",
	"    AIC(Group Variable + Covariates) < AIC(Covariates)",
	"           (and the difference is greater than 2)",
	"",
	"For example if a response variable has 6 levels (e.g. K6), and 5 of the levels",
	"have models with an improved AIC when a group's variables is included",
	"then 5 will be recorded in the cell.",
	"",
	"A response level is considered improved if the 'Group w/ Covariates' model",
	"is better than the 'Covariates-only' model"
));

summarize_groups_by_response=function(resp_res, aic_diff_cutoff=2){

	# Generate a num_groups by num_responses matrix
	# For each group, count up the the number of factor levels that have improved.

	response_names=names(resp_res);
	print(response_names);

	impr_mat_as_list=list();

	for(resp_ix in response_names){
		cat("Extracting: ", resp_ix, "\n");
		covariates_res=resp_res[[resp_ix]][["covariates"]];
		target_covar_res=resp_res[[resp_ix]][["target_covar"]];

		covariates_only_aic=covariates_res$aic;
		cat("Covariates Only\n");	
		print(covariates_only_aic);

		target_names=names(target_covar_res);

		impr_list_by_resp=list();
		for(targ_ix in target_names){
			cat("Target Group: ", targ_ix, "\n");
			target_covar_aic=target_covar_res[[targ_ix]][["aic"]];
			cat("Target w/ Covariates:\n");
			print(target_covar_aic);

			signf_impr=((covariates_only_aic-target_covar_aic)>aic_diff_cutoff);
			num_sigimp=sum(signf_impr);
			impr_list_by_resp[[targ_ix]]=num_sigimp;
		}

		impr_mat_as_list[[resp_ix]]=impr_list_by_resp;
	}

	print(impr_mat_as_list);

	grp_names=names(impr_mat_as_list[[1]]);
	rsp_names=names(impr_mat_as_list);

	num_rsp=length(rsp_names);
	num_grps=length(grp_names);

	grp_rsp_matrix=matrix(0, nrow=num_grps, ncol=num_rsp);
	colnames(grp_rsp_matrix)=rsp_names;
	rownames(grp_rsp_matrix)=grp_names;

	for(g_ix in grp_names){
		for(r_ix in rsp_names){
			grp_rsp_matrix[g_ix, r_ix]=impr_mat_as_list[[r_ix]][[g_ix]];
		}
	}
	print(grp_rsp_matrix);
	return(grp_rsp_matrix);

}

summary_matrix=summarize_groups_by_response(responses_results);

par(mfrow=c(1,1));
par(mar=c(1,1,1,1));
paint_matrix(summary_matrix, title="Num Factor Levels with Model Improvement when Group Incl", 
		deci_pts=0, label_zeros=F, value.cex=1);

out_sum_mat=cbind(rownames(summary_matrix), summary_matrix);
colnames(out_sum_mat)=c("Group", colnames(summary_matrix));
outfn=paste(OutputRoot, ".multn_resp.grp_sum.tsv", sep="");
write.table(out_sum_mat, outfn, quote=F, sep="\t", row.names=F, col.names=T);

###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
