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
options(width=120);

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

num_targets=length(targets_arr);
num_responses=length(responses_arr);
num_covariates=length(covariates_arr);


pdf(paste(OutputRoot, ".lasso.pdf", sep=""), height=11, width=9.5);

plot_text(params);

plot_text(c(
	paste("Target Vars (", num_targets, "):", sep=""),
	capture.output(print(targets_arr)),
	"",
	paste("Covariates (", num_covariates, "):", sep=""),
	capture.output(print(covariates_arr)),
	"",
	paste("Responses (", num_responses, "):", sep=""),
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
#----------------------------------------------------------------------	 

nstrict=function(x){
	sum(x=="S");
}

nmed=function(x){
	sum(x=="S" | x=="M");
}

nloose=function(x){
	sum(x!=".");
}

###############################################################################

output_variable_selection_table=function(alm, title="", abridge=F, lookup_tab=F){

	# The alm table contains ., L, M, and S.	
	# print(alm);

	#----------------------------------------------------------------------	 

	if(abridge){
		hits=apply(alm, 1, function(x){ any(x!=".")});
		alm=alm[hits,,drop=F];
	}

	categories=colnames(alm);
	colnames(alm)=1:ncol(alm);

	out_matrix=alm;
	
	#----------------------------------------------------------------------
	# Calculate row margin values
	NumStrict=apply(alm, 1, nstrict);
	NumMedium=apply(alm, 1, nmed);
	NumLoose=apply(alm, 1, nloose);

	out_matrix=cbind(out_matrix, 
		rep(" ", nrow(out_matrix)), 
		NumStrict, NumMedium, NumLoose);

	num_uniq_strict=sum(NumStrict>0);
	num_uniq_medium=sum(NumMedium>0);
	num_uniq_loose=sum(NumLoose>0);

	#----------------------------------------------------------------------
	# Calculate column margin values
	NumStrict=apply(alm, 2, nstrict);
	NumMedium=apply(alm, 2, nmed);
	NumLoose=apply(alm, 2, nloose);

	out_ncol=ncol(out_matrix);
	pad=rep("  ", out_ncol-ncol(alm));
	out_matrix=rbind(out_matrix, rep("  ", out_ncol));
	out_matrix=rbind(out_matrix, rep("  ", out_ncol));

	margin=rbind(
		c(NumStrict, pad),	
		c(NumMedium, pad),	
		c(NumLoose, pad));

	rownames(margin)=c("NumStrict", "NumMedium", "NumLoose");

	out_matrix=rbind(out_matrix, margin);

	#----------------------------------------------------------------------
	# Output 
	out_mat_txt=capture.output(print(out_matrix, quote=F));
	print(out_mat_txt);

	cat("\n");
	plot_text(c(
		title,
		"Selected Variables:",
		ifelse(abridge, "[Abridged/Only Selected]", "[All Predictors]"),
		"",
		out_mat_txt,
		"",
		"Total (Unique) Selected Variables: ",
		paste("	Strict: ", num_uniq_strict, sep=""),
		paste("	Medium: ", num_uniq_medium, sep=""),
		paste("	Loose:  ", num_uniq_loose, sep="")
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
all_responses_matrix=matrix(FALSE, ncol=num_responses, nrow=num_targets);
colnames(all_responses_matrix)=responses_arr;
rownames(all_responses_matrix)=targets_arr;

coef_list=list();

par(mfrow=c(3,1));
par(mar=c(5,5,7,1));

num_folds=10;
registerDoMC(num_folds);

for(resp in responses_arr){
#for(resp in responses_arr[1]){

	#-----------------------------------------------------------------------------
	# Remove NAs
	cat("Working on Response: ", resp, "\n");
	x=factors_loaded[,c(covariates_arr, targets_arr),drop=F];
	y=factors_loaded[,resp, drop=F];

	nonas=apply(cbind(x,y), 1, function(x){ all(!is.na(x));});
	x=as.matrix(x[nonas,,drop=F]);
	y=as.factor(y[nonas,]);

	resp_levels=levels(y);
	num_resp_levels=length(resp_levels);

	#-----------------------------------------------------------------------------

	plot_title_page(
		paste("Response: ", resp, sep=""),
		c(
			"",
			"",
			"Levels:",
			"",	
			paste("[", 1:num_resp_levels, "] ", resp_levels, "\n", sep="")
		)
	);

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

	cat("Min Deviation Lambda: ", mindev_lambda, 
		"  Log():", log_mindev_lambda, "\n", sep="");
	cat("Stricter Deviation Lambda: ", stricter_lambda, 
		"  Log():", log_stricter_lambda, "\n", sep="");
	cat("Looser Deviation Lambda: ", looser_lambda, 
		"  Log():", log_looser_lambda, "\n", sep="");

	#-----------------------------------------------------------------------------
	# Get the coefficients (selected predictors will have non-zero coefficients)

	coef_tabs=list();
	coef_tabs[["Loose"]]=coef(cvfit, looser_lambda);
	coef_tabs[["Mindev"]]=coef(cvfit, mindev_lambda);
	coef_tabs[["Strict"]]=coef(cvfit, stricter_lambda);

	cutoffs=c("Loose", "Mindev", "Strict");

	coef_mat=matrix(".", nrow=length(targets_arr), ncol=num_resp_levels);
	colnames(coef_mat)=resp_levels;
	rownames(coef_mat)=targets_arr;
	for(cutoff in cutoffs){
		for(resp_ix in 1:num_resp_levels){

			coefs_list=coef_tabs[[cutoff]];

			#cat("Cutoff: ", cutoff, "\n");
			#print(coefs_list[[resp_ix]]);

			# Skip over covariates AND intercept
			selected_ix=(0!=coefs_list[[resp_ix]][targets_arr,1]);

			coef_mat[selected_ix[targets_arr],resp_ix]=substr(cutoff,1,1);

		}
	}

	output_variable_selection_table(coef_mat, resp, abridge=F, lookup_tab=F);
	output_variable_selection_table(coef_mat, resp, abridge=T, lookup_tab=T);

	#-----------------------------------------------------------------------------
	# Plot and annotate cross validation deviances
	cat("Plotting CV: Lambda vs. Deviance...\n");
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
	cat("Plotting Selected Variables Per Level/Categories...\n");
	NumStrict=apply(coef_mat, 2, nstrict);
	NumMedium=apply(coef_mat, 2, nmed);
	NumLoose=apply(coef_mat, 2, nloose);
	barplot(NumLoose, 
		ylim=c(0, num_targets+1),
		main=paste(resp, ": Factor Level/Category Specific Variables", sep=""),
		las=2,
		col="green");
	barplot(NumMedium, col="black", xaxt="n", yaxt="n", add=T);
	barplot(NumStrict, col="blue", xaxt="n", yaxt="n", add=T);
	abline(h=num_targets, col="black", lwd=1);
	axis(side=2, at=num_targets, labels=paste("Targets: ", num_targets, sep=""),
		las=2, cex.axis=.7);

	#-----------------------------------------------------------------------------
	# Accumulate across all responses

	all_responses_matrix[targets_arr, resp]=
		apply(coef_mat, 1, function(x){
			if(any(x=="S")){
				return("S");
			}else if(any(x=="M")){
				return("M");
			}else if(any(x=="L")){
				return("L");
			}else{
				return(".");
			}
		});

	#output_variable_selection_table(all_responses_matrix, paste("ALL"), abridge=F, lookup_tab=F);

	cat("\n\n");
}

#------------------------------------------------------------------------------

cat("----------------------------------------------------------------------------------\n");
cat("Across ALL Responses:\n");
cat("----------------------------------------------------------------------------------\n");

plot_title_page("Selected Variables Across\nALL\nResponses");

output_variable_selection_table(all_responses_matrix, paste("ALL"), abridge=F, lookup_tab=F);
output_variable_selection_table(all_responses_matrix, paste("Selected"), abridge=T, lookup_tab=F);

par(mfrow=c(3,1));

plot_stacked_all_responses=function(arm){

	num_targets=nrow(arm);
	NumStrict=apply(arm, 2, nstrict);
	NumMedium=apply(arm, 2, nmed);
	NumLoose=apply(arm, 2, nloose);

	barplot(NumLoose,
		main="Selected Variables: All Responses",
		ylab="Num Uniq Sel Var",
		las=2,
		ylim=c(0, num_targets+1),
		col="green");
	barplot(NumMedium, col="black", xaxt="n", yaxt="n", add=T);
	barplot(NumStrict, col="blue", xaxt="n", yaxt="n", add=T);

	abline(h=num_targets, col="black", lwd=1);
	axis(side=2, at=num_targets, labels=paste("Targets: ", num_targets, sep=""),
		las=2, cex.axis=.7);

}

cat("Plotting Stacked Barplots Across All Responses...\n");
plot_stacked_all_responses(all_responses_matrix);

#------------------------------------------------------------------------------

plot_cdf_of_sel_var=function(arm){

	num_resp=ncol(arm);
	num_targets=nrow(arm);
	cutoffs=c("S", "M", "L");
	funcs=list("S"=nstrict, "M"=nmed, "L"=nloose);
	num_uniq_sel_var=list();
	max_sel=rep(0,3);
	names(max_sel)=cutoffs;

	for(co in cutoffs){
	
		num_uniq_sel_var[[co]]=numeric(num_resp);
		names(num_uniq_sel_var[[co]])=colnames(arm);

		for(i in 1:num_resp){
			counts=apply(arm[,1:i,drop=F], 1, funcs[[co]]);
			selected=counts>0
			num_selected=sum(selected);
		
			#cat("num resp: ", i, " / cutoff: ", co, "\n");
			#print(selected);

			num_uniq_sel_var[[co]][i]=num_selected;
			max_sel[co]=max(c(max_sel[co], num_selected));
		}

	}

	#print(num_uniq_sel_var);

	mids=barplot(num_uniq_sel_var[["L"]], col="green",
		ylim=c(0, num_targets+1),
		las=2,
		main="Cumulative Uniquely Selected Variables Across All Responses",
		ylab="Num Uniq Sel Var");

	title(main="(Assuming the response variables have a meaningful order.)",
		cex.main=.7, font.main=3, line=2);

	mids=barplot(num_uniq_sel_var[["M"]], col="black", xaxt="n", yaxt="n", add=T);
	mids=barplot(num_uniq_sel_var[["S"]], col="blue", xaxt="n", yaxt="n", add=T);

	abline(h=num_targets, col="black", lwd=1);
	abline(h=max_sel, col="purple", lty="dashed");

	axis(side=2, at=num_targets, labels=paste("Targets: ", num_targets, sep=""),
		las=2, cex.axis=.7);

	#----------------------------------------------------------------------

	# Plot legend
	plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", bty="n",
		xaxt="n", yaxt="n", main="");
	legend(0,1, legend=c(
		sprintf("Strict: %i (%2.1f %%)", max_sel["S"], 100*max_sel["S"]/num_targets),
		sprintf("MinDev: %i (%2.1f %%)", max_sel["M"], 100*max_sel["M"]/num_targets),
		sprintf("Loose : %i (%2.1f %%)", max_sel["L"],  100*max_sel["L"]/num_targets)
		), fill=c("blue", "black", "green"),
		cex=2
		);

}

cat("Plotting Cumulative Stacked Barplots Across All Responses...\n");
plot_cdf_of_sel_var(all_responses_matrix);

###############################################################################

cutoffs=c("S", "M", "L");
funcs=list("S"=nstrict, "M"=nmed, "L"=nloose);
num_selected=numeric(3);
names(num_selected)=cutoffs;

# Export selected
for(co in cutoffs){
	cutoff_counts=apply(all_responses_matrix, 1, funcs[[co]] );
	selected=cutoff_counts>0;
	select_var_arr=names(cutoff_counts[selected]);
	num_selected[co]=length(select_var_arr);
	write.table(
		x=select_var_arr, 
		file=paste(OutputRoot, ".lasso.sel.", co, ".lst", sep=""),
		quote=F, row.names=F, col.names=F);
}

###############################################################################

outroot=tail(strsplit(OutputRoot, "/")[[1]], 1);

sumfile=paste(OutputRoot, ".lasso.summary.tsv", sep="");
fh=file(sumfile, "w");
cat(file=fh, 
	c("Name", "NumTargets", "nStrict", "pStrict", "nMinDev", "pMinDev", "nLoose", "pLooose"),
	sep="\t");
cat(file=fh, "\n");

print(num_selected);
perc_sel=paste(round(num_selected/num_targets*100, 1), "%");
names(perc_sel)=cutoffs;

cat(file=fh,
	c(outroot,
	num_targets,
		num_selected["S"], perc_sel["S"],
		num_selected["M"], perc_sel["M"],
		num_selected["L"], perc_sel["L"]
	),
	sep="\t");	
		
cat(file=fh, "\n");

close(fh);

###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
