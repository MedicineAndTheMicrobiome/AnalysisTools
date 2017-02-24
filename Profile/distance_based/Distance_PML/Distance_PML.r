#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(glmnet);
library('getopt');
options(useFancyQuotes=F);
options(digits=5)

params=c(
	"distmat", "d", 1, "character",
	"factors", "f", 1, "character",
	"replaceNAs", "n", 2, "logical",
	"model_formula", "m", 2, "character",
	"outputroot", "o", 2, "character"
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
	"This script will utilize Penalized Maximal Likelihood regresssion\n",
	"in order to perform variable selection on a set of factors\n",
	"with the response being a distance matrix.\n",
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

check_res=check_factors(factors);

#print(check_res);

factor_file_desc_text=c(
	capture.output(print(check_res$info[,c("Avg", "Stdev", "Min", "Max", "NormDist")])),
	"",
	capture.output(print(check_res$info[,c("Med", "LB95", "UB95", "PropNonNAs")])),
	"",
	capture.output(print(check_res$info[,c("Proportion", "Lognormal", "PropUnique", "RecTrans")]))
);

plot_text(factor_file_desc_text);

##############################################################################
# Add recommended transformations to factors

factors=cbind(factors, check_res$transf);
#print(factors);

##############################################################################

# alpha=1 is lasso (i.e. solve for minimizing l1-norm) "sum of abs values"
# alpha=0 is ridge (i.e. solve for minimizing l2-norm) "sum of squares"
# alpha is the weighting between ridge and lasso behavior
# lambda is the weighting between ML and l-norm penalty
# ||b||p is the l-normal penalty, it is a function of the coefficients

factors_mod_matrix=as.matrix(factors);
print(factors_mod_matrix);

fit=glmnet(x=factors_mod_matrix, y=distmat, family="mgaussian", standardize=T, alpha=1);
print(fit);

# Plot of Coefficients vs. L1Norm(Coefficients)
plot(fit, xvar="norm", label=TRUE, ylim=c(-.3,.3))

# Plot of Coefficients vs. penalty strength>
#plot(fit, xvar="lambda", label=TRUE, ylim=c(-.3,.3))

# Plot of Coefficients vs R^2?
plot(fit, xvar="dev", label=TRUE, ylim=c(-.3,.3))

cvfit=cv.glmnet(x=factors_mod_matrix, y=distmat, family="mgaussian", standardize=T, alpha=1);
print(cvfit)
plot(cvfit);


##############################################################################

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
