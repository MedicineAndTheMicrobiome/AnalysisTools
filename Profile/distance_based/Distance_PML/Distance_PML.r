#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(glmnet);
library('getopt');
library(plotrix);
library(doMC);

registerDoMC(cores=8);

options(useFancyQuotes=F);
options(digits=5)

DEF_COR_EFF_CUTOFF=0.5;
DEF_COR_PVL_CUTOFF=0.05;
DEF_MAX_INTERACT=10;
DEF_ADD_TRANS=T;

params=c(
	"distmat", "d", 1, "character",
	"factors", "f", 1, "character",
	"replaceNAs", "n", 2, "logical",
	"model_formula", "m", 2, "character",
	"outputroot", "o", 2, "character",
	"cor_pval_cutoff", "c", 2, "numeric",
	"cor_eff_cutoff", "e", 2, "numeric",
	"num_max_interact", "i", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-d <distance matrix>\n",
	"	-f <factors>\n",
	"	[-n (replace NAs with something)]\n",
	"	[-m \"model formula string\"]\n",
	"	[-o <output filename root>]\n",
	"\n",
	"	[-c <correlation magnitude effect cutoff, default=", DEF_COR_EFF_CUTOFF, ">]\n",
	"	[-e <correlation pvalue cutoff, default=", DEF_COR_PVL_CUTOFF, ">]\n",
	"	[-i <maximum interactions terms to add, default=", DEF_MAX_INTERACT,">]\n",
	"	[-t <don't add recommended transforms, default=", DEF_ADD_TRANS,">]\n",
	"\n",
	"This script will utilize Penalized Maximal Likelihood regresssion\n",
	"in order to perform variable selection on a set of factors\n",
	"with the response being a distance matrix.\n",
	"\n",
	"The interactions terms to include will be based on selecting factors in the order\n",
	"of the pairs that are the least correlated to each other.\n",
	"\n",
	"Essentially, it is using LASSO to select variables PERMANOVA.\n",
	"\n");

if(!length(opt$distmat) || !length(opt$factors)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputFnameRoot=gsub(".distmat", "", opt$distmat);
}else{
	OutputFnameRoot=opt$outputroot;
}
OutputFnameRoot=paste(OutputFnameRoot, ".dist_pml", sep="");

if(!length(opt$model_formula)){
	ModelFormula="";
}else{
	ModelFormula=opt$model_formula;
}

if(!length(opt$replaceNAs)){
	RemoveSamples_wNA_Factors=F;
}else{
	RemoveSamples_wNA_Factors=T
}

CorPvalCutoff=DEF_COR_PVL_CUTOFF;
if(length(opt$cor_pval_cutoff)){
	CorPvalCutoff=opt$cor_pval_cutoff;
}

CorEffCutoff=DEF_COR_EFF_CUTOFF;
if(length(opt$cor_eff_cutoff)){
	CorEffCutoff=opt$cor_eff_cutoff;
}

NumMaxInteractions=DEF_MAX_INTERACT;
if(length(opt$num_max_interact)){
	NumMaxInteractions=opt$num_max_interact;
}

DistmatFname=opt$distmat;
FactorsFname=opt$factors;

cat("Distance Matrix Filename: ", DistmatFname, "\n", sep="");
cat("Factors Filename: ", FactorsFname, "\n", sep="");
cat("Output Filename Root: ", OutputFnameRoot, "\n", sep="");
cat("Remove Samples with NAs in Factors: ", RemoveSamples_wNA_Factors, "\n", sep=""); 
cat("\n");

if(ModelFormula!=""){
	cat("Model Formula specified: ", ModelFormula, "\n\n");
}

###############################################################################

load_distance_matrix=function(fname){
	distmat=as.matrix(read.delim(fname, sep=" ",  header=TRUE, row.names=1, check.names=FALSE, comment.char="", quote=""));
	#print(distmat);
	mat_dim=dim(distmat);
	cat("Read in distance matrix: \n");
	cat("  Rows: ", mat_dim[1], "\n");
	cat("  Cols: ", mat_dim[2], "\n");
	if(mat_dim[1]!=mat_dim[2]){
		cat("Error: Distance Matrix is not squared.\n");
		print(colnames(distmat));
		print(rownames(distmat));
		q(status=-1);
	}
	return(distmat);
}

##############################################################################

load_factors=function(fname){
	factors=as.data.frame(read.delim(fname,  header=TRUE, row.names=1, check.names=FALSE, sep="\t"));
	return(factors);
}

##############################################################################

plot_text=function(strings){
	orig_par=par(no.readonly=T);

        par(mfrow=c(1,1));
        par(family="Courier");
        par(oma=rep(.5,4));
        par(mar=rep(0,4));

        num_lines=length(strings);

        top=max(as.integer(num_lines), 40);

        plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                );
        for(i in 1:num_lines){
                #cat(strings[i], "\n", sep="");
                text(0, top-i, strings[i], pos=4, cex=.8);
        }

	par(orig_par);
}

##############################################################################

remove_samples_wNA=function(factors){
	num_factors=ncol(factors);
	num_samples=nrow(factors);
	
	keepers=numeric();
	for(i in 1:num_samples){
		if(!any(is.na(factors[i,]))){
			keepers=c(keepers, i);
		}
	}
	num_keepers=length(keepers);
	cat("Kept ", num_keepers, " of ", num_samples, " samples.\n");
	return(factors[keepers,, drop=F]);
}

check_factors=function(factor){
	# This function will try to transform the data to be more normally distributed	

	num_factors=ncol(factor);
	num_samples=nrow(factor);
	factor_names=colnames(factor);

	trans_types=c("Avg", "Stdev", "Med", "Min", "Max", "LB95", "UB95", "PropNonNAs", "NormDist", "Proportion", "Lognormal", "Count", "PropUnique", "RecTrans");
	num_types=length(trans_types);

	factor_info=matrix(0, nrow=num_factors, ncol=num_types, dimnames=list(factor_names, trans_types));
	
	shrink=function(x){
		shrink_factor=min(0.0001, min(diff(sort(x))));
		shrunk=x*(1-shrink_factor)+shrink_factor/2;
		return(shrunk);
	}

	transformed_matrix=matrix(NA, nrow=num_samples, ncol=num_factors);
	transform_name=character(num_factors);

	for(i in 1:num_factors){
		cur_fact=factor[,i, drop=F];

		# Find NAs
		blank_ix=(cur_fact=="");
		cur_fact[blank_ix]=NA;
		NA_ix=is.na(cur_fact);
		num_NAs=sum(NA_ix)
		cur_fact=cur_fact[!NA_ix];
		num_non_na_samples=length(cur_fact);
		factor_info[i, "PropNonNAs"]=1-(num_NAs/num_samples);

		# Get descriptive statistics
		factor_info[i, "Avg"]=mean(cur_fact);
		factor_info[i, "Stdev"]=sd(cur_fact);

		qs=quantile(cur_fact, c(0, .025, .5, .975, 1));
		factor_info[i, "Min"]=qs[1];
		factor_info[i, "LB95"]=qs[2];
		factor_info[i, "Med"]=qs[3];
		factor_info[i, "UB95"]=qs[4];
		factor_info[i, "Max"]=qs[5];

		# Test for normality
		st=shapiro.test(cur_fact);
		if(st$p.value>.05){
			factor_info[i, "NormDist"]=TRUE;
		}

		# If factors are not normally distributed, try log transforming it
		if(!factor_info[i, "NormDist"]){
			if(factor_info[i, "Min"]>0){
				transf=log(cur_fact+1);
				st=shapiro.test(transf);
				if(st$p.value>.05){
					factor_info[i, "Lognormal"]=1;
					factor_info[i, "RecTrans"]=TRUE;
					transform_name[i]=paste("log_", factor_names[i], sep="");
					transformed_matrix[!NA_ix, i]=transf;
				}
			}
		}

		# Are variables discontinuous?
		unique_val=unique(cur_fact);
		num_unique=length(unique_val);
		factor_info[i, "PropUnique"]=num_unique/num_samples;

		# If factors are proportions, try logit transform
		if(factor_info[i, "Min"]>=0 && factor_info[i, "Max"]<=1 && num_unique>=3){
			factor_info[i, "Proportion"]=TRUE;
			
			shrunk=shrink(cur_fact);
			logit=log(shrunk/(1-shrunk));
			st=shapiro.test(logit);
			if(st$p.value>.05){
				factor_info[i, "RecTrans"]=TRUE;
				transform_name[i]=paste("logit_", factor_names[i], sep="");
				transformed_matrix[!NA_ix, i]=logit;
			}
		}

		# If factors are counts, try sqrt 
		if(as.integer(cur_fact)==cur_fact && any(cur_fact>1)){
			factor_info[i, "Count"]=TRUE;
			transf=sqrt(cur_fact);
			factor_info[i, "RecTrans"]=TRUE;
			transform_name[i]=paste("sqrt_", factor_names[i], sep="");
			transformed_matrix[!NA_ix, i]=transf;
		}

	}
	
	rownames(transformed_matrix)=rownames(factor);
	colnames(transformed_matrix)=transform_name;
	transformed_factors=which(factor_info[,"RecTrans"]==1);

	results=list();
	results[["info"]]=factor_info;
	results[["transf"]]=transformed_matrix[,transformed_factors];
	return(results);	

}

compute_correl=function(factors, pval_cutoff=0.05, abs_correl_cutoff=0.5){
# This function will calculate the correlation and pvalue between all
# factors and then recommend interaction terms 

	num_factors=ncol(factors);
	factor_names=colnames(factors);
	pvalue_mat=matrix(1, nrow=num_factors, ncol=num_factors, dimnames=list(factor_names,factor_names));
	correl_mat=matrix(0, nrow=num_factors, ncol=num_factors, dimnames=list(factor_names,factor_names));
	corrltd_mat=matrix("", nrow=num_factors, ncol=num_factors, dimnames=list(factor_names,factor_names));

	for(i in 1:num_factors){
		for(j in 1:num_factors){
			if(i<j){
				# Remove NAs and then compute correl 
				not_na_ix=!(is.na(factors[,i] | is.na(factors[,j])));
				test_res=cor.test(factors[not_na_ix,i], factors[not_na_ix,j]);
				correl_mat[i,j]=test_res$estimate;
				pvalue_mat[i,j]=test_res$p.value;
				is_cor=abs(test_res$estimate)>abs_correl_cutoff && test_res$p.value<=pval_cutoff;
				corrltd_mat[i,j]=ifelse(is_cor, "X", ".");

				# Copy over symmetric values
				correl_mat[j,i]=correl_mat[i,j];
				pvalue_mat[j,i]=pvalue_mat[i,j];
				corrltd_mat[j,i]=corrltd_mat[i,j];
			}
		}
	}	

	results=list();
	results$pval_cutoff=pval_cutoff;
	results$correl_cutoff=abs_correl_cutoff;
	results$corrltd_mat=corrltd_mat;
	results$correl_mat=correl_mat;
	results$pvalue_mat=pvalue_mat;

	return(results);

}

add_interactions=function(factors, correl_results, max_interactions=10){
	
	# Convert correlations to magnitudes
	asymmetric=abs(correl_results$correl_mat);

	# Set the other side of the correlation to 1, so we don't count it
	num_factors=ncol(asymmetric);
	for(i in 1:num_factors){
		for(j in 1:num_factors){
			if(i<=j){
				asymmetric[i,j]=1;
			}
		}
	}
	
	print(asymmetric);

	# Order the correlations from smallest to 1
	correl_ord=order(as.vector(asymmetric));

	# Convert 1D position to 2D correlation matrix position
	to2D=function(x, side_len){
		x=x-1;
		return(c(x %% side_len, x %/% side_len)+1);
	}	

	# Compute max pos correl values
	max_correl_val=num_factors*(num_factors-1)/2;
	max_interactions=min(max_interactions, max_correl_val);

	# Store new interaction values:
	num_samples=nrow(factors);
	interact_mat=matrix(0, ncol=max_interactions, nrow=num_samples);
	rownames(interact_mat)=rownames(factors);

	# Compute interactions
	interact_names=character(max_interactions);
	factor_names=colnames(factors);
	least_correl_mat=correl_results$corrltd_mat;
	for(i in 1:max_interactions){
		#print(to2D(correl_ord[i], num_factors));
		fact_ix=to2D(correl_ord[i], num_factors);
		
		a=fact_ix[1];
		b=fact_ix[2];
		interact_mat[,i]=factors[,a]*factors[,b];
		interact_names[i]=paste(factor_names[a], "_x_", factor_names[b], sep="");
		least_correl_mat[a,b]=i;
		least_correl_mat[b,a]=i;

	}
	colnames(interact_mat)=interact_names;

	results=list();
	results$least_corrltd_mat=least_correl_mat;
	results$interaction_mat=interact_mat;
	
	return(results);

}

##############################################################################

pdf(paste(OutputFnameRoot, ".pdf", sep=""), height=8.5, width=11);

# Load distance matrix
distmat=load_distance_matrix(DistmatFname);
num_distmat_samples=ncol(distmat);
distmat_sample_names=colnames(distmat);
cat("\n");
#print(distmat):

# Load factors
factors=load_factors(FactorsFname);
factor_names=colnames(factors);
num_factors=ncol(factors);
num_factor_orig_samples=nrow(factors);
cat(num_factors, " Factor(s) Loaded:\n", sep="");
print(factor_names);
cat("\n");

if(ModelFormula!=""){
	# Based on factors in model string, identity which factors are used
	model_var=gsub("\\+", " ", ModelFormula);
	model_var=gsub("\\:", " ", model_var);
	model_var=gsub("\\*", " ", model_var);
	model_var=unique(strsplit(model_var, " ")[[1]]);
	model_var=intersect(model_var, factor_names);
	factors=factors[,model_var, drop=F];
	num_factors=ncol(factors);
	
}else{
	model_var=factor_names;
}

if(RemoveSamples_wNA_Factors){
	factors=remove_samples_wNA(factors);
}

factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);

# Confirm/Reconcile that the samples in the matrix and factors file match
common_sample_names=intersect(distmat_sample_names, factor_sample_names);
num_common_samples=length(common_sample_names);
if(num_common_samples < num_distmat_samples || num_common_samples < num_factor_samples){
	cat("\n");
	cat("*** Warning: The number of samples in factors file does not match those in your distance matrix. ***\n");
	cat("Taking intersection (common) sample IDs between both.\n");
	cat("Please confirm this is what you want.\n");
	cat("\tNum Distmat Samples: ", num_distmat_samples, "\n");
	cat("\tNum Factor  Samples: ", num_factor_samples, "\n");
	cat("\tNum Common  Samples: ", num_common_samples, "\n");	
	cat("\n");
}

# Set the working distance matrix to the same order
distmat=distmat[common_sample_names, common_sample_names];
num_samples=ncol(distmat);
sample_names=colnames(distmat);
cat("Num Samples used: ", num_samples, "\n\n");

factors=factors[common_sample_names, , drop=F];
for(i in 1:num_factors){
	categories=unique(unique(factors[,i]));
	cat("'", factor_names[i], "' has ", length(categories), " categories.\n", sep="");
	cat("\t", paste(categories, collapse=", "), sep="");
	cat("\n");
}

##############################################################################

intro_text=c(
	"Distance Matrix Penalized Maximum Likelihood",
	"",
	paste("Distance Matrix Filename: ", DistmatFname, sep=""),
	paste("Factors Filename: ", FactorsFname, sep=""),
	paste("Output Filename Root: ", OutputFnameRoot, sep=""),
	"",
	paste("Num Samples in Dist Mat: ", num_distmat_samples, sep=""),
	paste("Num Samples in Factor File: ", num_factor_orig_samples, sep=""),
	paste("Num Shared/Common Samples: ", num_common_samples, sep=""),
	"",
	paste("Num Factors in Factor File: ", num_factors, sep="")

);

#plot_text(intro_text);

##############################################################################
# Describe factors and recommend transfomrations

check_res=check_factors(factors);

#print(check_res);

factor_file_desc_text=c(
	"Normal Summary:",
	capture.output(print(check_res$info[,c("Avg", "Stdev", "Min", "Max", "NormDist")])),
	"",
	"Non-Parametric Summary:",
	capture.output(print(check_res$info[,c("Med", "LB95", "UB95", "PropNonNAs")])),
	"",
	"Inferred Data Types and Transformation Recommendation:",
	capture.output(print(check_res$info[,c("Proportion", "Lognormal", "PropUnique", "RecTrans")]))
);

plot_text(factor_file_desc_text);

##############################################################################
# Add recommended transformations to factors

factors=cbind(factors, check_res$transf);
#print(factors);

##############################################################################
# Compute correlation matrices even with NAs

correl_res=compute_correl(factors, CorPvalCutoff, CorEffCutoff);
plot_text(c(
	"Significantly Correlated Factors:",
	paste("P-value Cutoff: ", CorPvalCutoff, "  Effect Size Cutoff: ", CorEffCutoff, sep=""),
	"",
	capture.output(print(correl_res$corrltd_mat, quote=F))
));

##############################################################################
# Insert interaction terms for those factors least correlated

interactions_res=add_interactions(factors, correl_res, NumMaxInteractions);
num_interaction_terms_added=ncol(interactions_res$interaction_mat);

plot_text(c(
	"Automatically Generate Interations:",
	"(Based on pairs of factors with the least correlation.)",
	"",
	paste("Num of Terms Added: ", num_interaction_terms_added, sep=""),
	"",
	paste("    ", 1:num_interaction_terms_added, ".) ", colnames(interactions_res$interaction_mat), sep="")
));

plot_text(c(
	"Least Correlated Factors Matrix:",
	"(1 least correlated ... N most correlated)",
	"",
	capture.output(print(interactions_res$least_corrltd_mat, quote=F))
));

factors=cbind(factors, interactions_res$interaction_mat);

##############################################################################

# alpha=1 is lasso (i.e. solve for minimizing l1-norm) "sum of abs values"
# alpha=0 is ridge (i.e. solve for minimizing l2-norm) "sum of squares"
# alpha is the weighting between ridge and lasso behavior
# lambda is the weighting between ML and l-norm penalty
# ||b||p is the l-normal penalty, it is a function of the coefficients

factors_mod_matrix=as.matrix(factors);
factor_names=colnames(factors_mod_matrix);
print(factors_mod_matrix);
num_xs=ncol(factors_mod_matrix);

plot_text(c(
	paste(num_xs, " predictors/factors/variables entering into LASSO:", sep=""),
	paste("   ", 1:num_xs, ".) ", factor_names, sep="")
));

cat("Computing GLMNET Fit:\n");
fit=glmnet(x=factors_mod_matrix, y=distmat, family="mgaussian", standardize=T, alpha=1);
print(fit);

# Plot of Coefficients vs. L1Norm(Coefficients)
par(mfrow=c(5,3));
par(mar=c(0,1,0,1));
plot(fit, xvar="norm", label=TRUE, ylim=c(-10,30))

# Plot of Coefficients vs. penalty strength>
#plot(fit, xvar="lambda", label=TRUE, ylim=c(-.3,.3))

# Plot of Coefficients vs R^2?
#plot(fit, xvar="dev", label=TRUE, ylim=c(-2,2))

cat("Running Cross Validation...\n");
cvfit=cv.glmnet(x=factors_mod_matrix, y=distmat, family="mgaussian", standardize=T, alpha=1, parallel=T);

cv_num_var=cvfit$nzero;
cv_mean_cv_err=cvfit$cvm;

coefficients=fit$beta;
df=fit$df;
lambdas=fit$lambda;

y_idx_names=names(coefficients);
cat("Response Names: ");
print(y_idx_names);
num_lambdas=fit$dim[2];

# Compute median coefficient across samples at each Lambda for each x
median_coeff=matrix(0, nrow=num_lambdas, ncol=num_xs);
colnames(median_coeff)=factor_names;
# rows =x's
# col  =lambda's

#print(coefficients);
# The number of y's is the number of samples
for(lamb_ix in 1:num_lambdas){
	for(x_ix in 1:num_xs){
	
		# Accumulate coefficients across samples before calculating median
		across_samp=numeric(num_samples);
		for(y_ix in 1:num_samples){
			across_samp[y_ix]=coefficients[[y_idx_names[y_ix]]][x_ix, lamb_ix];
		}
		median_coeff[lamb_ix, x_ix]=median(abs(across_samp));
	}
}

plot_coefficients=function(coeff_mat, lambdas, mark_lambda_ix=NA, title=""){

	num_lambdas=nrow(coeff_mat);
	num_xs=ncol(coeff_mat);

	coef_range=range(coeff_mat);
	coef_span=diff(coef_range);
	extra_buf=coef_span/10;

	# Set up plot
	plot(0, xlim=c(0, num_lambdas), ylim=c(coef_range[1]-extra_buf, coef_range[2]+extra_buf), type="n",
		xaxt="n",
		ylab="Coefficients of Standardized Predictors",
		xlab="ML Penalty: Log10(Lambda)",
		bty="c"
	);

	# Mark the best lambda value
	if(!is.na(mark_lambda_ix)){
		abline(v=mark_lambda_ix, lty=2, col="black");
	}

	# Plot curves
	for(x_ix in 1:num_xs){
		points(coeff_mat[,x_ix], col=x_ix, type="l");
	} 

	# Label lambda/DFs positions
	x_axis_pos=floor(seq(1, num_lambdas, length.out=20));

	axis(side=3, at=x_axis_pos, labels=df[x_axis_pos], cex.axis=.5);
	axis(side=1, at=x_axis_pos, labels=round(log10(lambdas[x_axis_pos]),2), las=2, cex.axis=.5);
	axis(side=4, at=coeff_mat[num_lambdas, ], labels=factor_names, las=2, cex.axis=.5, lwd=0, lwd.tick=1, line=-1);
	title(main="Number of Variables", cex.main=1, font.main=1, line=2)
	title(main=title, cex.main=2, line=4)
}

par(mfrow=c(1,1));
par(mar=c(5, 5, 7, 8));

# Plot median coefficients across samples (y's) across all variables
plot_coefficients(median_coeff, lambdas, title="Median Magnitude of Coefficients Across All Samples");


# Plot cross validation error vs num variables
cv_num_lambdas=length(cvfit$lambda);
cv_num_var=cvfit$nzero;
cv_mean_err=cvfit$cvm;
cv_min_err=min(cv_mean_err);
cv_min_err_ix=which(cv_min_err==cv_mean_err);
cv_min_err_num_var=cvfit$nzero[cv_min_err_ix];
cv_min_err_lambda=cvfit$lambda[cv_min_err_ix];

plot(cv_num_var, cv_mean_cv_err, main="Influence of Variable Inclusion on Prediction Error", xlab="Number of Variables Included", ylab="Mean CV Error");

plotCI(log10(cvfit$lambda), cvfit$cvm, ui=cvfit$cvup, li=cvfit$cvlo, col="red", scol="grey",
	pch=16, 
	xlab="Log10(Lambda)",
	ylab="Mean Cross-Validated Error"	
);

abline(h=cv_min_err, col="blue", lty=2);
abline(v=log10(cv_min_err_lambda), col="blue", lty=2);
x_axis_pos=floor(seq(1, cv_num_lambdas, length.out=20));
axis(side=3, at=log10(cvfit$lambda[x_axis_pos]), labels=cvfit$nzero[x_axis_pos], cex.axis=.5);
title(main="Number of Variables", cex.main=1, font.main=1, line=2)
title(main="Influence of ML Penalty on Prediction Error ", cex.main=2, line=4)

# Get Variables at min error
all_min_error_coeff=numeric();
for(samp_ix in 1:num_samples){
	cv_min_err_coeff=cvfit$glmnet.fit$beta[[y_idx_names[samp_ix]]][,cv_min_err_ix];
	all_min_error_coeff=rbind(all_min_error_coeff, cv_min_err_coeff);
}
rownames(all_min_error_coeff)=sample_names;
print(all_min_error_coeff);

# Extract the predictor names based on the coefficients that were non zero
non_zero_xs=apply(all_min_error_coeff, 2, function(x){ return(!all(x==0))});
non_zero_coeff=all_min_error_coeff[,non_zero_xs, drop=F];
num_nonzero_coeff=ncol(non_zero_coeff);
print(num_nonzero_coeff);
non_zero_x_names=colnames(non_zero_coeff);
print(non_zero_x_names);

# Output which coefficients were selected:
plot_text(c(
	"Top Selected Predictors Based on Mimimum Mean CV Error:",
	paste("   ", 1:num_nonzero_coeff, ".)", non_zero_x_names, sep="")
));


# Plot median coefficients across samples (y's) across all variables zoomed
zoom_ix=df <= (num_nonzero_coeff+1);
plot_coefficients(
	median_coeff[zoom_ix,], 
	lambdas,
	mark_lambda_ix=cv_min_err_ix,
	title=paste("Median Magnitude of Coefficients Across All Samples (DF < ", num_nonzero_coeff, "+1 )", sep="")
);







par(mfrow=c(3,3));
for(x_ix in 1:num_nonzero_coeff){
	hist(non_zero_coeff[,x_ix], xlim=c(-2, 2), main=non_zero_x_names[x_ix]);
}



# For each sample, plot coefficients

# For each x, plot coefficients across all samples

#

##############################################################################

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
