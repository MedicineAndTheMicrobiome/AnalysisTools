#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

options(useFancyQuotes=F);
options(width=200);

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
options(width=200);

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

title_page=function(title, subtitle="", notes=""){

	par(mfrow=c(1,1));
	plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n",  xaxt="n", yaxt="n",
		xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
		);

	text(.5, .8, title, adj=c(.5, -1), cex=4, font=2); 
	text(.5, .5, subtitle, adj=c(.5, 0), cex=2, font=1); 
	text(.5, .4, notes, adj=c(.5, 1), cex=1, font=3); 
}

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, deci_pts=4, 
	label_zeros=T, counts=F, value.cex=2, 
	plot_col_dendr=F,
	plot_row_dendr=F
){

        num_row=nrow(mat);
        num_col=ncol(mat);

	row_names=rownames(mat);
	col_names=colnames(mat);

	orig.par=par(no.readonly=T);

	cat("Painting Matrix: ", title, "\n");
        cat("Num Rows: ", num_row, "\n");
        cat("Num Cols: ", num_col, "\n");


	if(num_row==0 || num_col==0){
		plot(0, type="n", xlim=c(-1,1), ylim=c(-1,1), xaxt="n", yaxt="n", bty="n", xlab="", ylab="",
			main=title);
		text(0,0, "No data to plot...");
		return();
	}

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
                plot_min=min(c(mat,Inf), na.rm=T);
        }
        if(is.na(plot_max)){
                plot_max=max(c(mat,-Inf), na.rm=T);
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

	if(num_row==1){
		plot_row_dendr=F;
	}
	if(num_col==1){
		plot_col_dendr=F;
	}

	# Don't plot dendrogram if there are any NAs in the matrix
	#if(any(is.na(mat))){
	#	plot_col_dendr=F;
	#	plot_row_dendr=F;
	#}

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
		mat=mat[row_dendr[["names"]],,drop=F];
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

                        if(mat[y,x]!=0 || label_zeros){
                                if(counts){
                                        text_lab=sprintf("%i", mat[y,x]);
                                }else{
                                        text_lab=sprintf(paste("%0.", deci_pts, "f", sep=""), mat[y,x]);
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


pdf(paste(OutputRoot, ".multn_resp.pdf", sep=""), height=11, width=9.5);

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

accumulate_sumfit_to_matrix=function(sumfit, mat, cname){

	coef_tab=sumfit[["coefficients"]]	
	varnames=setdiff(rownames(coef_tab), "(Intercept)"); 
	
	mat[["coef"]][varnames, cname]=coef_tab[varnames, "Estimate"];
	mat[["pval"]][varnames, cname]=coef_tab[varnames, "Pr(>|z|)"];

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

	empty_pval_coef=list();
	empty_pval_coef[["coef"]]=empty_matrix;
	empty_pval_coef[["pval"]]=empty_matrix;

	responses_results[[resp_ix]]=list();

	responses_results[[resp_ix]][["full"]]=empty_pval_coef;
	responses_results[[resp_ix]][["covariates"]]=empty_pval_coef;
	responses_results[[resp_ix]][["target_reduced"]]=list()
	responses_results[[resp_ix]][["target_covar"]]=list();

	for(tg in test_group_names){
		responses_results[[resp_ix]][["target_reduced"]][[tg]]=empty_pval_coef;
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
			accumulate_sumfit_to_matrix(full_sum_fit, 
				responses_results[[resp_ix]][["covariates"]], category_ix);
	
		for(tg in test_group_names){
	
			cat("\tFitting Targeted Models (", tg, "):\n");

			# Test Targeted Reduced (All Targets except Target, + Covariates)
			cat("\t\tFitting Targeted Reduced Model...\n");
			tar_red_formula=as.formula(paste("responses ~ ", targeted_reduced_model_list[tg]));
			tar_red_fit=glm(tar_red_formula, data=as.data.frame(tmp_factors), family="binomial");
			tar_red_sum_fit=summary(tar_red_fit);
			#print(tar_red_sum_fit);

			responses_results[[resp_ix]][["target_reduced"]][[tg]]=
				accumulate_sumfit_to_matrix(full_sum_fit, 
					responses_results[[resp_ix]][["target_reduced"]][[tg]], category_ix);



			# Test Target (Target + Covariates)
			cat("\t\tFitting Targeted+Covariates Model...\n");
			tar_formula=as.formula(paste("responses ~ ", targets_cov_model_list[tg]));
			tar_fit=glm(tar_formula, data=as.data.frame(tmp_factors), family="binomial");
			tar_sum_fit=summary(tar_fit);
			#print(tar_sum_fit);

			responses_results[[resp_ix]][["target_covar"]][[tg]]=
				accumulate_sumfit_to_matrix(full_sum_fit, 
					responses_results[[resp_ix]][["target_covar"]][[tg]], category_ix);


		}
	}

	cat("\n");
}

print(responses_results);


###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
