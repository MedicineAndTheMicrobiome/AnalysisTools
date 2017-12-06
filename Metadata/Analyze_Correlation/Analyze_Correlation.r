#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

options(useFancyQuotes=F);
options(digits=5)

params=c(
	"factors", "f", 1, "character",
	"outputroot", "o", 1, "character",
	"response", "r", 1, "character",
	"predictor", "p", 2, "character",
	"pc_coverage", "c", 2, "numeric"
);

PCA_COVERAGE=.95;

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <factors/metadata file name>\n",
	"	-o <output filename root.\n",
	"\n",
	"	-r <response variable list>\n",
	"	[-p <predictor variable list>]\n",
	"\n",
	"	[-c PC Coverage, default=", PCA_COVERAGE, "\n",
	"\n",
	"This script will perform a PCA analysis\n",
	"across the response variables and\n",
	"recommend a set of principal coordinates.\n",
	"\n",
	"If the predictor variable list is specified\n",
	"then plots between each predictor (x) and\n",
	"all the responses (y) will be made.\n",
	"\n",
	"The predictors will be clustered based on\n",
	"their correlation magnitude.\n",
	"\n");

if(
	!length(opt$factors) || 
	!length(opt$outputroot) || 
	!length(opt$response)
){
	cat(usage);
	q(status=-1);
}

FactorsFname=opt$factors;
OutputFnameRoot=opt$outputroot;
ResponseListName=opt$response;

PredictorListName=opt$predictor;
PCCoverage=PCA_COVERAGE;

cat("\n");
cat("Factor File Name: ", FactorsFname, "\n");
cat("Output File Name Root: ", OutputFnameRoot, "\n");
cat("Response List Name: ", ResponseListName, "\n");
cat("Predictor List Name: ", PredictorListName, "\n");
cat("PC Min Coverage: ", PCCoverage, "\n");
cat("\n");

###############################################################################

load_factors=function(fname){
	factors=as.data.frame(read.delim(fname,  header=TRUE, row.names=1, check.names=FALSE, sep="\t"));
	return(factors);
}

load_list=function(fname){
	lst=scan(fname, what=character());
	return(lst);	
}

test_log=function(mat_val){
	nrows=nrow(mat_val);
	ncols=ncol(mat_val);

	trans_mat=mat_val;
	orig_colnames=colnames(mat_val);
	new_colnames=character();

	for(var in orig_colnames){
		values=mat_val[,var];
		test_res=shapiro.test(values);

		if(test_res$p.value<.80){
			log_values=log(values+1);
			test_log_res=shapiro.test(log_values);

			if(test_log_res$p.value < test_res$p.value){
				# Keep original
				new_colnames=c(new_colnames, var);
			}else{
				# Keep log transformed
				new_colnames=c(new_colnames, paste("log(", var, ")", sep=""));
				trans_mat[, var]=log_values;
			}
		}else{
			new_colnames=c(new_colnames, var);
		}
	}
	colnames(trans_mat)=new_colnames;
	return(trans_mat);
}

##############################################################################
# Main Program Starts Here!
##############################################################################

pdf(paste(OutputFnameRoot, ".pdf", sep=""), height=11, width=8.5);

# Load factors
cat("Loading Factors...\n");
loaded_factors=load_factors(FactorsFname);
loaded_factor_names=colnames(loaded_factors);
loaded_sample_names=rownames(loaded_factors);

cat("Loaded factors:\n");
print(loaded_factor_names);
cat("\n");
cat("Loaded sample ids:\n");
print(loaded_sample_names);
cat("\n");

# Subset factors
responses_arr=load_list(ResponseListName);
predictors_arr=load_list(PredictorListName);

cat("\n");
cat("Targeted Responses:\n");
print(responses_arr);
missing=setdiff(responses_arr, loaded_factor_names);
if(length(missing)>0){
	cat("Missing:\n");
	print(missing);
	quit(status=-1);
}
cat("\n");

cat("Targeted Predictors:\n");
print(predictors_arr);
missing=setdiff(predictors_arr, loaded_factor_names);
if(length(missing)>0){
	cat("Missing:\n");
	print(missing);
	quit(status=-1);
}
cat("\n");

num_pred=length(predictors_arr);
num_resp=length(responses_arr);

pred_mat=loaded_factors[, predictors_arr, drop=F];
resp_mat=loaded_factors[, responses_arr, drop=F];

pred_mat=test_log(pred_mat);
resp_mat=test_log(resp_mat);

predictors_arr=colnames(pred_mat);
responses_arr=colnames(resp_mat);

##############################################################################

par(oma=c(1,1,4,1));

for(pred in predictors_arr){
	par(mfrow=c(4,3));
	for(resp in responses_arr){
		cat("Working on: ", resp, " vs ", pred, "\n");
		pred_val=pred_mat[,pred];
		resp_val=resp_mat[,resp];

		plot(pred_val, resp_val, xlab=pred, ylab=resp, main=resp);
		fit=lm(resp_val~pred_val);
		sumfit=summary(fit);

		# Draw regression line
		abline(fit, col="blue");

		# Pull regression stats
		coeff_pval=sumfit$coefficients["pred_val", "Pr(>|t|)"];
		coeff_estm=sumfit$coefficients["pred_val", "Estimate"];
		rsquared=sumfit$r.squared;

		stat_info=paste(
			"coeff=", signif(coeff_estm,2), 
			"  p-val=", sprintf("%3.3f", coeff_pval),
			"  R^2=", signif(rsquared,2), 
			sep="");

		mtext(stat_info, side=3, outer=F, cex=.5, col="blue");
		mtext(pred, side=3, outer=T, cex=2);
	}
}

##############################################################################

compute_correlations=function(mat){
	num_col=ncol(mat);
	cor_mat=matrix(0, nrow=num_col, ncol=num_col);
	pval_mat=matrix(0, nrow=num_col, ncol=num_col);
	rownames(cor_mat)=colnames(mat);
	colnames(cor_mat)=colnames(mat);
	rownames(pval_mat)=colnames(mat);
	colnames(pval_mat)=colnames(mat);
	for(i in 1:num_col){
		for(j in 1:i){
			v1=mat[,i];
			v2=mat[,j];
			notna=!(is.na(v1) | is.na(v2));
			#cor_mat[i,j]=cor(v1[notna], v2[notna]);
			test=cor.test(v1[notna], v2[notna]);
			pval_mat[i,j]=test$p.value;
			pval_mat[j,i]=test$p.value;
			cor_mat[i,j]=test$estimate;
			cor_mat[j,i]=test$estimate;
		}
	}
	res=list();
	res[["val"]]=cor_mat;
	res[["pval"]]=pval_mat;;
	res[["dist"]]=as.dist(1-abs(cor_mat));
	return(res);
}

correl=compute_correlations(cbind(pred_mat, resp_mat));
print(correl$val);
print(correl$pval);
print(correl$dist);

par(mfrow=c(1,1));
par(mar=c(15,2,1,2));
hcl=hclust(correl$dist, method="ward.D2");
dend=as.dendrogram(hcl);

highlight_predictors=function(x){
	if(is.leaf(x)){
		leaf_attr=attributes(x);
		label=leaf_attr$label;
		print(label);
		if(any(label==predictors_arr)){
			color="red";
			font=2;
		}else{
			color="black";
			font=1;
		}
		attr(x, "nodePar")=c(leaf_attr$nodePar, list(lab.font=font, lab.col=color, cex=0));
	}
	return(x);
}
dend=dendrapply(dend, highlight_predictors);

plot(dend, main="Ward's Minimum Variance: dist(1-abs(cor))");

pca=princomp(correl$val, cor=T);
print(pca);

##############################################################################


cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
