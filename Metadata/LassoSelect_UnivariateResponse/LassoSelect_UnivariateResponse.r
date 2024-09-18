#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
library('glmnet');
library('doMC');

source('~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r');

options(useFancyQuotes=F);

params=c(
	"factor_fn", "f", 1, "character",
	"subjectid_cn", "s", 1, "character",
	"responses_cn", "r", 1, "character",
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
	"	-s <Subject ID Column Name>\n",
	"	-r <Response Column Name>\n",
	"	-c <Covariates List>\n",
	"	-t <Target Predictors>\n",
	"	-o <Output Filename Root>\n",
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
	length(opt$responses_cn) &&
	length(opt$covariates_fn) && 
	length(opt$target_var_fn) && 
	length(opt$outputroot)
)){
	cat(usage);
	q(status=-1);
}

FactorsFile=opt$factor_fn;
ResponseColname=opt$responses_cn;
CovariatesFile=opt$covariates_fn;
TargetVarFile=opt$target_var_fn;
OutputRoot=opt$outputroot;
SubjectIDColname=opt$subjectid_cn;

params=capture.output({
cat("\n");
cat(script_name, "\n", sep="");
cat("\n");
cat("      Factors File: ", FactorsFile, "\n", sep="");
cat("Subject ID Colname: ", SubjectIDColname, "\n", sep="");
cat(" Responses Colname: ", ResponseColname, "\n", sep="");
cat("   Covariates File: ", CovariatesFile, "\n", sep="");
cat("  Targ. Var.  File: ", TargetVarFile, "\n", sep="");
cat("       Output Root: ", OutputRoot, "\n", sep="");
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

factors_loaded=load_factors_file(FactorsFile, SubjectIDColname);
available_variables=colnames(factors_loaded);

targets_arr=load_list(TargetVarFile);
covariates_arr=load_list(CovariatesFile);

num_targets=length(targets_arr);
num_covariates=length(covariates_arr);

pdf(paste(OutputRoot, ".uni_lasso.pdf", sep=""), height=11, width=9.5);

plot_text(params);

plot_text(c(
	paste("Target Vars (", num_targets, "):", sep=""),
	capture.output(print(targets_arr)),
	"",
	paste("Covariates (", num_covariates, "):", sep=""),
	capture.output(print(covariates_arr)),
	"",
	paste("Responses: ", ResponseColname, sep="")
));

all_used_variables=c(targets_arr, ResponseColname, covariates_arr);

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
cat("\nResponse:\n");
print(ResponseColname);
cat("\n");

###############################################################################

num_targets=length(targets_arr);

coef_list=list();

par(mfrow=c(3,1));
par(mar=c(5,5,7,1));

num_folds=10;
registerDoMC(num_folds);

#-----------------------------------------------------------------------------
# Remove NAs
x=factors_loaded[,c(covariates_arr, targets_arr),drop=F];
y=factors_loaded[,ResponseColname, drop=F];

#cat("X:\n");
#print(x);

#cat("Y:\n");
#print(y);

nonas=apply(cbind(x,y), 1, function(x){ all(!is.na(x));});
x=x[nonas,,drop=F];

y=unlist(y[nonas, ResponseColname, drop=F]);


#-----------------------------------------------------------------------------

layout_mat=matrix(c(1,1,2,3,4), ncol=1);
layout(layout_mat);

yexp=exp(y);
ylog=log(y);
ysqrt=sqrt(y);

y_normp=shapiro.test(y)$p.val;
yexp_normp=shapiro.test(yexp)$p.val;
ylog_normp=shapiro.test(ylog)$p.val;
ysqrt_normp=shapiro.test(ysqrt)$p.val;

par(mar=c(5,5,5,2));
hist(y, main=ResponseColname, breaks=20, cex.main=2, 
	col=ifelse(y_normp<.1, "red", "green"));
title(main=sprintf("Shapiro-Wilks P-value: %3.4f", y_normp), line=0, cex.main=.95, font.main=3);

par(mar=c(6,8,6,8));
hist(yexp, main=paste("exp(", ResponseColname, ")", sep=""), breaks=20, cex.main=1.2,
	col=ifelse(yexp_normp<.1, "pink", "lightgreen"));
title(main=sprintf("Shapiro-Wilks P-value: %3.4f", yexp_normp), line=1.2, cex.main=.95, font.main=3);

hist(ylog, main=paste("log(", ResponseColname, ")", sep=""), breaks=20, cex.main=1.2,
	col=ifelse(ylog_normp<.1, "pink", "lightgreen"));
title(main=sprintf("Shapiro-Wilks P-value: %3.4f", ylog_normp), line=1.2, cex.main=.95, font.main=3);

hist(ysqrt, main=paste("sqrt(", ResponseColname, ")", sep=""), breaks=20, cex.main=1.2,
	col=ifelse(ysqrt_normp<.1, "pink", "lightgreen"));
title(main=sprintf("Shapiro-Wilks P-value: %3.4f", ysqrt_normp), line=1.2, cex.main=.95, font.main=3);

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

lasso_family="gaussian";

#set.seed(100);

cat("CrossValidation GLMNet:\n");
x=as.matrix(x);
y=as.matrix(y);
cvfit=cv.glmnet(x, y, family=lasso_family, 
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

msg_lambda=capture.output({
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
cat("Coefficients Matrix:\n");
print(coef_mat);

selected_var=apply(coef_mat, 2, function(x){sum(x!=0);});
cat("\n");
cat("Num of Variables Selected:\n");
print(selected_var);


par(mfrow=c(2,1));
#-----------------------------------------------------------------------------
# Plot and annotate cross validation deviances
cat("Plotting CV: Lambda vs. Deviance...\n");
plot(cvfit);
title(main=paste("Response: ", ResponseColname, sep=""), line=4);
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
cat("Plotting Selected Variables Per Level/Categories...\n");
NumStrict=selected_var["Strict"];
NumMedium=selected_var["Mindev"];
NumLoose=selected_var["Loose"];

mids=barplot(c(NumLoose, NumMedium, NumStrict),
	ylim=c(0, num_targets+1),
	main=paste("Response: ", ResponseColname, sep=""),
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

for(c_ix in cutoffs){

	# Get variable list by cutoff
	var_arr=varlist[[c_ix]];
	
	# Fit model to confirm selected variables
	model_str=paste("y ~ ", paste(c(covariates_arr, var_arr), collapse="+"), sep="");
	print(model_str);
	fit=lm(as.formula(model_str), data=as.data.frame(cbind(y,x)));
	sumfit=summary(fit);

	#cat("[FIT]------------------------------------------------------>>\n");
	#print(fit);
	#cat("[FIT]------------------------------------------------------<<\n");
	#print(names(fit));

	#cat("[SUM_FIT]------------------------------------------------------>>\n");
	#print(sumfit);
	#cat("[SUM_FIT]------------------------------------------------------<<\n");
	#print(names(sumfit));

	#print(sumfit$coefficients);
	sumfit$coefficients=sumfit$coefficients[order(sumfit$coefficients[,"Pr(>|t|)"]),];
	#print(sumfit);
	print(sumfit);

	sumfit_txt=capture.output({print(sumfit)});
	# Output Coefficients/Pvalues Table
	plot_text(c(
		paste("Cutoff: ", c_ix, sep=""),
		sumfit_txt
		));

	# Plot obs/pred points
	obspred_rng=range(c(y, fit$fitted.values));
	plot(y, fit$fitted.values, xlab="Observed", ylab="Predicted",
		xlim=obspred_rng, ylim=obspred_rng, main=c_ix
		);

	# Draw reference and obs/pred lines
	model_fit=lm(fit$fitted.values~y);
	abline(a=0, b=1, col="blue", lty="dashed");
	abline(model_fit, col="black");

	# Output selected variable list
	if(length(var_arr)){

		outvarlist_fn=paste(OutputRoot, ".", c_ix, ".lst", sep="");
		writeLines(
			var_arr, 
			outvarlist_fn
			);

		outvarfactors_fn=paste(OutputRoot, ".", c_ix, ".tsv", sep="");

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


###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
