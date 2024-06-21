#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
library('glmnet');
library('doMC');

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

output_variable_selection_table=function(alm, title="", abridge=F, lookup_tab=F){

	if(abridge){
		hits=apply(alm, 1, any);
		alm=alm[hits,,drop=F];
	}

	categories=colnames(alm);
	colnames(alm)=1:ncol(alm);

	out_matrix=apply(alm, 1:2, function(x){ ifelse(x, "X", ".")});

	out_matrix=cbind(out_matrix, rep(" ", nrow(out_matrix)));

	NumResp=apply(alm, 1, sum);
	NumUniqResp=apply(alm, 1, function(x){ any(x)});
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
		ifelse(abridge, "[Abridged/Only Selected]", "[All Predictors]"),
		"",
		out_mat_txt,
		"",
		paste("Total (Unique) Selected Variables: ", num_unique_responses, sep="")
	));

	if(lookup_tab){
		plot_text(c(
			title,
			"Category/Level Names",
			paste(1:ncol(alm), ". ", categories, sep="")
		));
	}
}

###############################################################################

num_responses=length(responses_arr);
num_targets=length(targets_arr);
all_lasso_matrix=matrix(FALSE, ncol=num_responses, nrow=num_targets);
colnames(all_lasso_matrix)=responses_arr;
rownames(all_lasso_matrix)=targets_arr;

coef_list=list();

par(mfrow=c(3,1));
par(mar=c(5,5,7,1));

num_folds=20;
registerDoMC(num_folds);

for(resp in responses_arr){
#for(resp in responses_arr[1]){
	plot_title_page(resp, "response");

	#-----------------------------------------------------------------------------
	# Remove NAs
	cat("Working on Response: ", resp, "\n");
	x=factors_loaded[,c(covariates_arr, targets_arr),drop=F];
	y=factors_loaded[,resp, drop=F];

	nonas=apply(cbind(x,y), 1, function(x){ all(!is.na(x));});
	x=as.matrix(x[nonas,,drop=F]);
	y=as.factor(y[nonas,]);

	#-----------------------------------------------------------------------------
	# Set up covariates as penalty free

	incl_cov=intersect(covariates_arr, colnames(x));
	num_incl_cov=length(incl_cov);

	num_predictors=ncol(x);

	cat("Num Predictors: ", num_predictors, "\n");

	# Set up penalty factor for covariates, so no shrinkage is possible.
	penalty_fact_arr=rep(1, num_predictors);
	if(num_incl_cov>0){
		penalty_fact_arr[1:num_incl_cov]=0.0;
	}

	cat("Penalty Factors: 0's for covariates (i.e. no shrinkage allowed):\n");
	print(penalty_fact_arr);

	#-----------------------------------------------------------------------------
	# Run cv.glmnet

	cat("CrossValidation GLMNet:\n");
	cvfit=cv.glmnet(x,y, family="multinomial", 
		parallel=TRUE,
		nfolds=num_folds,
		keep=TRUE, # Keep fold information
		penalty.factor=penalty_fact_arr);

	# Gives minimum mean CV error
	mindev_lambda=cvfit$lambda.min;
	log_mindev_lambda=log(mindev_lambda);

	# Gives less regularized/penalized (more selected variables) within 1se of lowest 
	stricter_lambda=cvfit$lambda.1se;
	log_stricter_lambda=log(stricter_lambda);

	# Looser lambda
	log_looser_lambda=log_mindev_lambda - (log_stricter_lambda-log_mindev_lambda);
	looser_lambda=exp(log_looser_lambda);

	cat("Min Deviation Lambda: ", mindev_lambda, "Log():", log_mindev_lambda, "\n", sep="");
	cat("Stricter Deviation Lambda: ", stricter_lambda, "Log():", log_stricter_lambda, "\n", sep="");
	cat("Looser Deviation Lambda: ", looser_lambda, "Log():", log_looser_lambda, "\n", sep="");

	#-----------------------------------------------------------------------------
	# Get the coefficients (selected predictors will have non-zero coefficients)

	coef_tabs=coef(cvfit, cvfit$lambda.min);

	num_levels=length(coef_tabs);

	coef_mat=matrix(NA, nrow=length(targets_arr), ncol=num_levels);
	colnames(coef_mat)=names(coef_tabs);
	rownames(coef_mat)=targets_arr;
	for(i in 1:num_levels){
		# Skip over covariates AND intercept
		coef_mat[targets_arr,i]=(0!=coef_tabs[[i]][targets_arr,1]);
	}

	print(coef_mat);

	output_variable_selection_table(coef_mat, resp, abridge=F, lookup_tab=F);
	output_variable_selection_table(coef_mat, resp, abridge=T, lookup_tab=T);

	#-----------------------------------------------------------------------------
	# Plot and annotate cross validation deviances
	plot(cvfit);
	title(main=paste("Response: ", resp, sep=""), line=4);
	title(main=paste("Num Variables Selected (including ", num_incl_cov, " required covariates)", 
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
	var_counts=apply(coef_mat, 2, sum);
	barplot(var_counts, 
		main=paste(resp, ": Factor Level/Category Specific Variables", sep=""),
		las=2);

	#-----------------------------------------------------------------------------
	# Accumulate across all responses
	nonzero_targets=apply(coef_mat, 1, any);
	targ_names=rownames(coef_mat[nonzero_targets,,drop=F]);
	all_lasso_matrix[targ_names, resp]=TRUE;
	#output_variable_selection_table(all_lasso_matrix, paste("ALL"), abridge=F, lookup_tab=F);

	cat("\n\n");
}


plot_cdf_of_sel_var=function(sel_mat){

	num_resp=ncol(sel_mat);
	num_uniq_sel_var=numeric(num_resp);

	for(i in 1:num_resp){
		num_uniq_sel_var[i]=sum(apply(sel_mat[,1:i,drop=F], 1, any));
	}

	names(num_uniq_sel_var)=colnames(sel_mat);

	max_sel=max(num_uniq_sel_var);
	mids=barplot(num_uniq_sel_var, las=2, main="Cumulative Uniquely Selected Variables",
		ylab="Num Uniq Sel Var", ylim=c(0, max_sel*1.1)
		);
	text(mids, num_uniq_sel_var, num_uniq_sel_var, font=3, pos=3);
	abline(h=max_sel, col="blue", lty="dashed");

}



par(mfrow=c(3,1));

plot_title_page("Selected Variables Across\nALL\nResponses");

var_per_resp=apply(all_lasso_matrix, 2, sum);
barplot(var_per_resp, main="Num Selected Var Per Response");

plot_cdf_of_sel_var(all_lasso_matrix);

output_variable_selection_table(all_lasso_matrix, paste("ALL"), abridge=F, lookup_tab=T);
output_variable_selection_table(all_lasso_matrix, paste("Selected"), abridge=T, lookup_tab=F);

###############################################################################

# Export selected
selected_across_all=apply(all_lasso_matrix, 1, any);
select_var_arr=rownames(selected_across_all);
write.table(
	x=select_var_arr, 
	file=paste(OutputRoot, ".lasso.selected.lst", sep=""),
	quote=F, row.names=F, col.names=F);

###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
