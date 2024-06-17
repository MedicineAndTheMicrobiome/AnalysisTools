#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
library('glmnet');

options(useFancyQuotes=F);

params=c(
	"factor_fn", "f", 1, "character",
	"responses_fn", "r", 1, "character",
	"covariates_fn", "c" , 2, "character",
	"target_var_fn", "t", 2, "character",
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
	"	-t <Target Predictors>\n",
	"	-o <Output Filename Root>\n",
	"\n",
	"This script will LASSO the following multinomial logistic regression model:\n",
	"\n",
	"	[multi-nomial response] = [covariates] + var1 + var2 + ... + varn\n",
	"\n",
	"\n",
	"The Factor file should contain all the variables for all the subjects.\n",
	"The Response variables is a list of all the variables (columns) that should be\n",
	"   independently analyzed.\n",
	"\n",
	"The Covariates variables are always included in the model as predictors.\n",
	"\n",
	"\n", sep="");

if(!(length(opt$factor_fn) && 
	length(opt$responses_fn) &&
	length(opt$covariates_fn) && 
	length(opt$target_var_fn) && 
	length(opt$output)
)){
	cat(usage);
	q(status=-1);
}

FactorsFile=opt$factor_fn;
ResponsesFile=opt$responses_fn;
CovariatesFile=opt$covariates_fn;
TargetVarFile=opt$target_var_fn;
OutputRoot=opt$outputroot;

params=capture.output({
cat("\n");
cat(script_name, "\n", sep="");
cat("\n");
cat("    Factors File: ", FactorsFile, "\n", sep="");
cat(" Responses File : ", ResponsesFile, "\n", sep="");
cat(" Covariates File: ", CovariatesFile, "\n", sep="");
cat("Targ. Var.  File: ", TargetVarFile, "\n", sep="");
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

##############################################################################
##############################################################################

factors_loaded=load_factors(FactorsFile);
available_variables=colnames(factors_loaded);

targets_arr=load_list(TargetVarFile);
responses_arr=load_list(ResponsesFile);
covariates_arr=load_list(CovariatesFile);


pdf(paste(OutputRoot, ".lasso.pdf", sep=""), height=11, width=9.5);

plot_text(params);

plot_text(c(
	"Target Vars:",
	capture.output(print(targets_arr)),
	"",
	"Covariates:",
	capture.output(print(covariates_arr)),
	"",
	"Responses:",
	capture.output(print(responses_arr))
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

###############################################################################

cat("\nTargets:\n");
print(targets_arr);
cat("\nCovariates:\n");
print(covariates_arr);
cat("\nResponses:\n");
print(responses_arr);
cat("\n");

#data(MultinomialExample);
#x=MultinomialExample$x;
#y=MultinomialExample$y;
#print(x);
#print(y);
#fit=glmnet(x,y);

#quit();

###############################################################################

num_responses=length(responses_arr);
num_targets=length(targets_arr);
all_lasso_matrix=matrix(FALSE, ncol=num_responses, nrow=num_targets);
colnames(all_lasso_matrix)=responses_arr;
rownames(all_lasso_matrix)=targets_arr;

coef_list=list();

par(mfrow=c(3,1));
par(mar=c(5,5,7,1));
for(resp in responses_arr){

	cat("Working on Response: ", resp, "\n");
	x=factors_loaded[,c(covariates_arr, targets_arr),drop=F];
	y=factors_loaded[,resp, drop=F];

	nonas=apply(cbind(x,y), 1, function(x){ all(!is.na(x));});
	x=as.matrix(x[nonas,,drop=F]);
	y=as.factor(y[nonas,]);

	incl_cov=intersect(covariates_arr, colnames(x));
	num_incl_cov=length(incl_cov);
	
	#cat("X: Predictors:\n");
	#print(x);
	#cat("\n");
	#cat("Y: Response: \n");
	#print(y);

	num_predictors=ncol(x);

	cat("Num Predictors: ", num_predictors, "\n");

	# Set up penalty factor for covariates, so no shrinkage is possible.
	penalty_fact_arr=rep(1, num_predictors);
	if(num_incl_cov>0){
		penalty_fact_arr[1:num_incl_cov]=0.0;
	}

	cat("Penalty Factors: 0's for covariates (i.e. no shrinkage allowed):\n");
	print(penalty_fact_arr);

	# Run GLMNET
	#fit=glmnet(x,y,family="multinomial", penalty.factor=penalty_fact_arr);
	#print(fit);
	#plot(fit, label=T);
	#print(summary(fit));
	#print(coef(fit));

	cat("CrossValidation GLMNet:\n");
	cvfit=cv.glmnet(x,y, family="multinomial", penalty.factor=penalty_fact_arr);
	min_lambda=cvfit$lambda.min;
	max_lambda=cvfit$lambda.1se;
	plot(cvfit);
	title(main=paste("Response: ", resp, sep=""), line=4);
	title(main=paste("Num Variables Selected (including ", num_incl_cov, " required covariates)", 
		sep=""), line=2.1, cex.main=1, font.main=1);
	axis(side=3, at=log(min_lambda), labels="More", cex=.3, line=-1, col.axis="blue", tick=F);
	axis(side=3, at=log(max_lambda), labels="Less", cex=.3, line=-1, col.axis="red", tick=F);
	

	# Gives minimum mean CV error
	#cvfit$lambda.min;
	# Gives most regularized (fewest selected variables) within 1se of 
	#cvfit$lambda.1se;

	coef_tabs=coef(cvfit, cvfit$lambda.min);

	print(coef_tabs[[1]]);
	
	first_coef_tab=coef_tabs[[1]];
	nonzero_ix=first_coef_tab[,1]!=0;

	sel_pred_names=names(first_coef_tab[nonzero_ix,1]);
	
	cat("Selected Predictors:\n");
	print(sel_pred_names);

	selected_targets=intersect(sel_pred_names, targets_arr);
	all_lasso_matrix[selected_targets,resp]=TRUE;

	cat("\n\n");
}

###############################################################################

# all_lasso_matrix is T/F boolean matrix
print(all_lasso_matrix);

sel_var_across_all_resp=apply(all_lasso_matrix, 1, function(x){ any(x);});

###############################################################################

# Generate output summary matrix

output_variable_selection_table=function(alm, title=""){
	out_matrix=apply(alm, 1:2, function(x){ ifelse(x, "X", ".")});
	out_matrix=cbind(out_matrix, " ");

	NumResp=apply(alm, 1, sum);
	NumUniqResp=apply(all_lasso_matrix, 1, function(x){ any(x)});
	out_matrix=cbind(out_matrix, NumResp);

	out_matrix=rbind(out_matrix, "");
	NumSelected=c(apply(alm, 2, sum), " ", "");
	out_matrix=rbind(out_matrix, NumSelected);

	out_mat_txt=capture.output(print(out_matrix, quote=F));

	num_unique_responses=sum(NumUniqResp);

	print(out_mat_txt);
	cat("\n");
	cat("Total Unique Selected Variables: ", num_unique_responses, "\n");
	plot_text(c(
		title,
		"Selected Variables (Min Lambda / Max Variables):",
		"",
		out_mat_txt,
		"",
		paste("Total (Unique) Selected Variables: ", num_unique_responses, sep="")
	));
}

output_variable_selection_table(all_lasso_matrix, title="All Variables");

selected_lasso_matrix=all_lasso_matrix[sel_var_across_all_resp,,drop=F];

print(selected_lasso_matrix);

output_variable_selection_table(selected_lasso_matrix, title="Selected Variables");

###############################################################################

# Export selected

select_var_arr=rownames(selected_lasso_matrix);
write.table(
	x=select_var_arr, 
	file=paste(OutputRoot, ".lasso.selected.lst", sep=""),
	quote=F, row.names=F, col.names=F);

###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
