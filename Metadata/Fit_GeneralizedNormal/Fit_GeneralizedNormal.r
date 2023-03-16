#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
library('gnorm'); # for gnrom functions
library('stats4'); # For mle

options(useFancyQuotes=F);
options(digits=5)

params=c(
	"factors", "f", 1, "character",
	"targets", "t", 1, "character",
	"outputroot", "o", 1, "character",
	"debug_targets", "d", 2, "character",
	"fast_fit", "F", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

NORM_PVAL_CUTOFF=.2;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <factors/metadata file name>\n",
	"	-t <list of variable targets>\n",
	"	-o <output filename root>\n",
	"\n",
	"	[-d <debug list of variable targets>]\n",
	"	[-F (fast fit)>]\n",
	"\n",
	"This script will try to fit the variables based on\n",
	"how much 'information' the variables have.  Since\n",
	"one can't strictly compare entropy between continuous\n",
	"distributions, and the variance of a variable is\n",
	"dependent on scale, we fit the generalized normal distribution\n",
	"and then later preferentially analyze variables with a larger\n",
	"shape parameters (flatter).\n",
	"\n",
	"This script only performs the fit and exports the beta for\n",
	"each variable.\n",
	"\n",
	"The -d debug list lets you specify the names of variables\n",
	"for subsetting.  In case some variables keep failing the fit\n",
	"optimization.\n",
	"\n",
	"The -F performs a less thorough, but not bad fit, by restricting\n",
	"the search space.  Use this for testing.\n", 
	"\n",
	"The fit_gnd.raw.pdf contains the log/sqrt transformations to normal if necessary,\n",
	"as well as the fits based on the order in the input factors/metadata file.\n",
	"\n",
	"The nLL_sorted and beta_sorted pdf files are sorted based on quality of fit and\n",
	"beta values.\n",
	"\n");

if(
	!length(opt$factors) || 
	!length(opt$targets) || 
	!length(opt$outputroot)
){
	cat(usage);
	q(status=-1);
}

FactorsFname=opt$factors;
OutputFnameRoot=opt$outputroot;
TargetsFname=opt$targets;

DebugVarList="";
if(length(opt$debug_targets)){
	DebugVarList=opt$debug_targets;	
}

FastFit=F;
if(length(opt$fast_fit)){
	FastFit=T;
}

if(FastFit){
	OutputFnameRoot=paste(OutputFnameRoot, ".fast", sep="");
}

param_text=capture.output({
	cat("\n");
	cat("Factor File Name: ", FactorsFname, "\n");
	cat("Output File Name Root: ", OutputFnameRoot, "\n");
	cat("Targets: ", TargetsFname, "\n");
	cat("Fast Fit? ", FastFit, "\n");
	cat("\n");
});

print(param_text, quote=F);

###############################################################################

load_factors=function(fname){
	factors=as.data.frame(read.delim(fname,  header=TRUE, row.names=1, check.names=FALSE, sep="\t"));
	return(factors);
}

load_list=function(fname){
	cat("Loading: ", fname, "\n");
	lst=read.delim(fname, header=F, check.names=F, comment.char="#", as.is=T);
	return(lst[,1]);	
}

test_and_apply_log_transform=function(mat_val, pval_cutoff=.2, plot_before_after=T){
	nrows=nrow(mat_val);
	ncols=ncol(mat_val);

	trans_mat=mat_val;
	orig_colnames=colnames(mat_val);
	new_colnames=character();

	if(plot_before_after){
		orig_par=par(no.readonly=T);
		par(mfrow=c(5,2));
	}

	delete_list=c();
	for(var in orig_colnames){
		values=mat_val[,var];


		log_transformed=F;
		sqrt_transformed=F;

		num_unique_val=length(setdiff(unique(values), NA));
		values_nona=values[!is.na(values)];
		num_nona=length(values_nona);

		if(!is.numeric(values_nona)){
			cat("Error: Values not numeric for: ", var, "\n", sep="");
			print(values_nona);
		}

		if(num_nona<=3){
			cat("Not enough non NA values to measure normality.\n");
			new_colnames=c(new_colnames, var);
			next;
		}

		if(any(values_nona<0)){
			cat("Negative values.  Skipping transformation.\n");
			new_colnames=c(new_colnames, var);
			next;
		}

		cat("\n", var, ": Num Unique Values: ", num_unique_val, "\n");

		if(num_unique_val>1){

			test_res=shapiro.test(values);
			test_log_res=NULL;

			if(test_res$p.value<=pval_cutoff && num_unique_val>2){
				cat(" Not normal: ", test_res$p.value, "\n");

				if(any(values_nona==0)){
					log_values=log(values+1);
				}else{
					log_values=log(values);
				}
				sqrt_values=sqrt(values);

				test_log_res=shapiro.test(log_values);
				test_sqrt_res=shapiro.test(sqrt_values);

				if(test_log_res$p.value < test_res$p.value && 
					test_sqrt_res$p.value < test_res$p.value){
					# Keep original
					cat("  No Improvement: ", test_log_res$p.value, "\n");
					new_colnames=c(new_colnames, paste("orig_", var, sep=""));
				}else{

					cat("     Log p-val : ", test_log_res$p.value, "\n");
					cat("    Sqrt p-val : ", test_sqrt_res$p.value, "\n");
					
					if(test_log_res$p.value > test_sqrt_res$p.value){
						# Keep log transformed
						cat("  Log Transformation Effective: ", test_log_res$p.value, "\n");
						new_colnames=c(new_colnames, paste("log_", var, sep=""));
						trans_mat[, var]=log_values;
						log_transformed=T;		
					}else{
						# Keep sqrt transformed
						cat("  Sqrt Transformation Effective: ", test_sqrt_res$p.value, "\n");
						new_colnames=c(new_colnames, paste("sqrt_", var, sep=""));
						trans_mat[, var]=sqrt_values;
						sqrt_transformed=T;		
					}
				}
			}else{
				cat(" Normal enough. ", test_res$p.value, "\n");
				new_colnames=c(new_colnames, var);
			}

		}else{
			cat("  All values identical, removing...\n");
			new_name=paste("all_ident_", var, sep="");
			new_colnames=c(new_colnames, new_name);
			delete_list=c(delete_list, new_name);
		}

		if(plot_before_after){
			nclass=nclass.Sturges(values)*4;

			hist(values, main=var, breaks=nclass);
			title(main=sprintf("p-value: %4.4f", test_res$p.value), cex.main=.8, line=.5);
			
			if(log_transformed){
				hist(log_values, breaks=nclass, main=paste("log(", var,")", sep=""));
				title(main=sprintf("p-value: %4.4f", test_log_res$p.value), cex.main=.8, line=.5);
			}else if(sqrt_transformed){
				hist(sqrt_values, breaks=nclass, main=paste("sqrt(", var,")", sep=""));
				title(main=sprintf("p-value: %4.4f", test_sqrt_res$p.value), cex.main=.8, line=.5);
			}else{
				plot(0,0, xlab="", ylab="", main="", xaxt="n", yaxt="n", bty="n", type="n");

				if(test_res$p.value>pval_cutoff){
					text(0,0, "Transform not necessary");
				}else{
					text(0,0, "Transform not beneficial");
				}
			}
		}

	}

	colnames(trans_mat)=new_colnames;

	trans_mat=trans_mat[,setdiff(new_colnames, delete_list),drop=F];

	if(plot_before_after){
		par(orig_par);
	}

	return(trans_mat);
}

plot_text=function(strings){

	orig.par=par(no.readonly=T);

	par(mfrow=c(1,1));
	par(family="Courier");
	par(oma=rep(.5,4));
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

	par(orig.par);
}


##############################################################################
# Main Program Starts Here!
##############################################################################

pdf(paste(OutputFnameRoot, ".fit_gnd.raw.pdf", sep=""), height=11, width=8.5);

plot_text(param_text);

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

target_var_arr=load_list(TargetsFname);
cat("Num of Targets: ", length(target_var_arr));
target_variables=intersect(loaded_factor_names, load_list(TargetsFname));
targeted_mat=loaded_factors[,target_variables, drop=F];

cat("Testing Predictor Variables for normality.\n");
curated_targets_mat=test_and_apply_log_transform(targeted_mat, NORM_PVAL_CUTOFF);

if(DebugVarList!=""){
	debug_tar_list=scan(DebugVarList, what=character(), comment.char="#");
	cat("Debug target variable list:\n");
	print(debug_tar_list);
	curated_targets_mat=curated_targets_mat[,debug_tar_list,drop=F];
}


curated_target_names=colnames(curated_targets_mat);


##############################################################################

num_variables=length(curated_target_names);

params_matrix=matrix(NA, nrow=num_variables, ncol=6);
colnames(params_matrix)=c("nLL", "pert_ix", "method", "mu", "alpha", "beta");
rownames(params_matrix)=curated_target_names;

gn_density=list();

normal_estim_matrix=matrix(NA, nrow=num_variables, ncol=2);
colnames(normal_estim_matrix)=c("mu", "sd");
rownames(normal_estim_matrix)=curated_target_names;

norm_density=list();


failed_variables=c();
par(mfrow=c(2,2));

# If using fast fitting, use smaller search space and fewer repeats
if(!FastFit){
	# Complete
	mu_pert=3;
	alpha_pert=5;
	beta_pert=7;

	pert_range=4;
	max_successes=20;
}else{
	# Fast
	mu_pert=3;
	alpha_pert=3;
	beta_pert=3;

	pert_range=3;
	max_successes=3
}

# Compute range and steps for perturbation
mu_adj=2^(seq(-pert_range, pert_range, length.out=mu_pert));
alpha_adj=2^(seq(-pert_range, pert_range, length.out=alpha_pert));
beta_adj=2^(seq(-pert_range, pert_range, length.out=beta_pert));

total_perts=mu_pert*alpha_pert*beta_pert;

# Compute how far each perturbation is from 'no perturbation', and order.
pert_matrix=matrix(0, nrow=total_perts, ncol=3);
colnames(pert_matrix)=c("mu", "alpha", "beta");
pert_distance=numeric(total_perts);
pert_ix=1;
for(i in 1:mu_pert){
	for(j in 1:alpha_pert){
		for(k in 1:beta_pert){

			pert_matrix[pert_ix,]=c(mu_adj[i], alpha_adj[j], beta_adj[k]);

			pert_distance[pert_ix]=
				sqrt(log2(mu_adj[i])^2+log2(alpha_adj[j])^2+log2(beta_adj[k])^2);

			pert_ix=pert_ix+1;
		}
	}
}

pert_dist_sort_ix=order(pert_distance);
pert_matrix=pert_matrix[pert_dist_sort_ix,];

cat("Pertubration Matrix (Sorted):\n");
print(pert_matrix);

##############################################################################

MIN_REAL=1e-323;

num_curated_targets=length(curated_target_names);
counter=1;
for(ctn in curated_target_names){
	data=curated_targets_mat[,ctn];

	cat("********************************************************************\n");
	cat("Working on: ", ctn, " (", counter, "/", num_curated_targets,")\n");

	mu0=mean(data);
	variance=var(data);
	alpha0=sqrt(variance*2);
	beta0=2;

	normal_estim_matrix[ctn, "mu"]=mu0;
	normal_estim_matrix[ctn, "sd"]=sqrt(variance);

	nLL=function(mu, alpha, beta){
		if(alpha<MIN_REAL){alpha=MIN_REAL;};
		if(beta<MIN_REAL){beta=MIN_REAL};
		neg_loglik=-sum(dgnorm(data, mu, alpha, beta, log=T));
		#neg_loglik=ifelse(neg_loglik==Inf, 1e200, neg_loglik);
		neg_loglik=ifelse(neg_loglik==Inf, 1e308, neg_loglik);
		if(!is.finite(neg_loglik)){
			cat(mu, ", ", alpha, ", ", beta, ": ", neg_loglik, "\n", sep="");
		}
		return(neg_loglik);
	};

	cat("Initial Values: mean=",mu0, " alpha=", alpha0, 
		" (var=", variance, ") beta=", beta0, "\n", sep="");
	muf=NA; alphaf=NA; betaf=NA;

	pert_ix=1;

	pert_res_matrix=matrix(NA, nrow=total_perts, ncol=6);
	colnames(pert_res_matrix)=c("nLL", "pert_ix", "method", "mu", "alpha", "beta");
	
	num_success=0;

	while(num_success<max_successes && pert_ix<=total_perts){
		success=F;

		#cat("Perturbations: ", pert_ix, "  Trials: ", tries, "\n");

		muTry=mu0*pert_matrix[pert_ix, "mu"];
		alphaTry=alpha0*pert_matrix[pert_ix, "alpha"];
		betaTry=beta0*pert_matrix[pert_ix, "beta"];

		method_ix=1
		for(method in c("Nelder-Mead", "SANN")){
			method_success=F;
			tryCatch(
				{

					cat(method, ": Adj Init: (",
						muTry, ", ",
						alphaTry, ", ",
						betaTry,"): \n",sep="");

					names(muTry)=NULL;
					names(alphaTry)=NULL;
					names(betaTry)=NULL;

					start_param=list(mu=muTry, alpha=alphaTry, beta=betaTry);

					#parscale=c(1,1,1);
					parscale=c(abs(mu0)/5, alpha0/5, beta0/10);

					if(method=="SANN"){
						control_var=list(parscale=parscale);
					}else{
						control_var=list(parscale=parscale);
					}

					#fit=mle(nLL, start=start_param, method);
					fit=mle(nLL, start=start_param, method, control=control_var);

					#print(fit);
					fit_res=attributes(fit);
					#print(fit_res);

					muf=fit_res$coef[1];	
					alphaf=max(fit_res$coef[2], MIN_REAL);
					betaf=max(fit_res$coef[3], MIN_REAL);
					num_success=num_success+1;
					method_success=T;
				}, 
				error=function(e){
					cat("Error:\n");
					print(e);
				}
			);
			if(method_success){
				cat(method, ": Success!\n");
				break;
			}else{
				cat(method, ": Failure.\n");
			}

			method_ix=method_ix+1;
		}

		if(method_success){
			nll_pert=nLL(mu=muf, alpha=alphaf, beta=betaf);

			pert_res_matrix[pert_ix,]=c(nll_pert, pert_ix, method_ix, muf, alphaf, betaf);
		}

		pert_ix=pert_ix+1;		

		cat("\n");
	}

	# Find perturbation with min nll_pert
	nona_ix=!is.na(pert_res_matrix[,"nLL"]);
	nona_pert_res_matrix=pert_res_matrix[nona_ix,,drop=F];

	
	print(nona_pert_res_matrix);

	min_nLL=min(nona_pert_res_matrix[,"nLL"]);
	min_ix=min(which(nona_pert_res_matrix[,"nLL"]==min_nLL));

	nll_best=nona_pert_res_matrix[min_ix,"nLL"];
	pert_ix_best=nona_pert_res_matrix[min_ix, "pert_ix"];
	method_best=nona_pert_res_matrix[min_ix,"method"];
	mu_best=nona_pert_res_matrix[min_ix,"mu"];
	alpha_best=nona_pert_res_matrix[min_ix,"alpha"];
	beta_best=nona_pert_res_matrix[min_ix,"beta"];

	cat("Best: (", mu_best, ", ", alpha_best, ", ", beta_best, ")\n", sep="");

	hist(data, breaks=30, main=paste(counter, ".) ", ctn, sep=""), freq=F, xlab="");

	if(num_success>0){
		cat("Overall Success.\n")
		title(main=paste("Mean =", round(mu_best,4)), line=-1);
		title(main=paste("Alpha =", round(alpha_best,4)), line=-2);
		title(main=paste("Beta =", round(beta_best,4)), line=-3);

		range=range(data);
		span=diff(range);
		margin=span*.15;
		x=seq(range[1]-margin, range[2]+margin, length.out=100);
		#print(c(muf, alphaf, betaf));
		y=dgnorm(x, mu_best, alpha_best, beta_best);
		#print(x);
		#print(y);
		points(x,y, col="blue");
	
		params_matrix[ctn,]=c(nll_best, pert_ix_best, method_best, mu_best, alpha_best, beta_best);
		gn_density[[ctn]]=list();
		gn_density[[ctn]][["x"]]=x;
		gn_density[[ctn]][["y"]]=y;

		y=dnorm(x, normal_estim_matrix[ctn, "mu"], normal_estim_matrix[ctn, "sd"]);
		points(x,y, col="red", cex=.5);
		norm_density[[ctn]]=list();
		norm_density[[ctn]][["x"]]=x;
		norm_density[[ctn]][["y"]]=y;
		norm_density[[ctn]][["mu"]]=mu0;
		norm_density[[ctn]][["sd"]]=sqrt(variance);
		norm_density[[ctn]][["nLL"]]=-sum(dnorm(data, mu0, sqrt(variance), log=T));

	}else{
		cat("Overall Failure.\n");
		title(main="Unsuccessful Fit.", line=-1);
		failed_variables=c(failed_variables, ctn);	
	}


	cat("\n\n");
	counter=counter+1;
}

##############################################################################

dev.off();

##############################################################################
# Report Failed Variables

if(length(failed_variables)){

	failed_fn=paste(OutputFnameRoot, ".failed.meta.tsv", sep="");
	fh=file(failed_fn, "w");
	cat(file=fh, "SampleID\t");
	close(fh);

	write.table(
		curated_targets_mat[,failed_variables,drop=F], 
		file=failed_fn, append=T,
		quote=F, sep="\t");

	#----------------------------------------------------------------------

	failed_fn=paste(OutputFnameRoot, ".failed.lst", sep="");
	fh=file(failed_fn, "w");
	cat(file=fh, paste(failed_variables, collapse="\n"), "\n", sep="");
	close(fh);
}

##############################################################################
##############################################################################
# Write out the variable histograms sorted by decreasing beta

plot_hist=function(d, var_name, params, var_ix, gden, nden){

	nLL=params[var_name, "nLL"];
	pert_ix=params[var_name, "pert_ix"];
	method=params[var_name, "method"];
	mu=params[var_name,"mu"];
	alpha=params[var_name,"alpha"];
	beta=params[var_name,"beta"];

	hist(data, main=paste(var_ix, ".) ", var_name, sep=""), xlab="", freq=F, breaks=30,
		border="grey50", col="grey75" );

	xy_bounds=par()$usr;
	xspan=xy_bounds[2]-xy_bounds[1];
	left=xy_bounds[1]+xspan*3/8;
	right=xy_bounds[1]+xspan*5/8;
	top=xy_bounds[4];

	text(left, top,
		paste(
			"[Generalized Norm Dist]\n",
			"nLL = ", sprintf("%g", nLL), "\n",
			"Mean = ", sprintf("%g", mu), "\n",
			"(StDev = ", sprintf("%g", alpha/sqrt(2)), ")\n",
			"Alpha = ", sprintf("%g", alpha), "\n",
			"Beta = ", sprintf("%g", beta), "\n",
			"Method = ", ifelse(method==1, "Nelder-Mead", "SANN"), "\n", 
			"pert iter = ", pert_ix, "\n",
		sep=""),
		cex=.7, col="blue", adj=c(1,1));

	text(right, top,
		paste(
			"[Standard Norm Dist]\n",
			"nLL = ", sprintf("%g", nden[[var_name]][["nLL"]]), "\n",
			"Mean = ", sprintf("%g", nden[[var_name]][["mu"]]), "\n",
			"StDev = ", sprintf("%g", nden[[var_name]][["sd"]]), "\n",
			"Method = MoM\n", sep=""),
		cex=.7, col="red", adj=c(0,1));

	points(gden[[var_name]][["x"]], gden[[var_name]][["y"]], col="blue");
	points(nden[[var_name]][["x"]], nden[[var_name]][["y"]], col="red", cex=.5);
}

#------------------------------------------------------------------------------

pdf(paste(OutputFnameRoot, ".beta_sorted.pdf", sep=""), height=11, width=8.5);
par(mfrow=c(3,2));

beta_arr=params_matrix[,"beta"];
beta_arr_sort_ix=order(beta_arr, decreasing=F);

params_matrix_beta_sort=params_matrix[beta_arr_sort_ix,, drop=F];
beta_arr_sorted=params_matrix_beta_sort[, "beta"];
mean_arr_sorted=params_matrix_beta_sort[, "mu"];
var_sorted=rownames(params_matrix_beta_sort);

counter=1;
for(varname in var_sorted){
	data=curated_targets_mat[,varname];
	plot_hist(data, varname, params_matrix, counter, gn_density, norm_density);	
	counter=counter+1;
}

dev.off();

#------------------------------------------------------------------------------

pdf(paste(OutputFnameRoot, ".nLL_sorted.pdf", sep=""), height=11, width=8.5);
par(mfrow=c(3,2));

nLL_arr=params_matrix[,"nLL"];
nLL_arr_sort_ix=order(nLL_arr, decreasing=F);

params_matrix_nLL_sort=params_matrix[nLL_arr_sort_ix,, drop=F];
nLL_var_sorted=rownames(params_matrix_nLL_sort);

counter=1;
for(varname in nLL_var_sorted){
	data=curated_targets_mat[,varname];
	plot_hist(data, varname, params_matrix, counter, gn_density, norm_density);	
	counter=counter+1;
}

dev.off();

##############################################################################
# Output transformed variables and beta values

transformed_factors_fn=paste(OutputFnameRoot, ".transf_var.tsv", sep="");

fh=file(transformed_factors_fn, "w");
cat(file=fh, "SampleID\t");
close(fh);
write.table(
	curated_targets_mat, 
	file=transformed_factors_fn, append=T,
	quote=F, sep="\t");

#------------------------------------------------------------------------------

estimated_betas_fn=paste(OutputFnameRoot, ".betas.tsv", sep="");
fh=file(estimated_betas_fn, "w");
cat(file=fh, "Variable\t");
close(fh);

write.table(
	params_matrix[,"beta", drop=F],
	file=estimated_betas_fn, append=T,
	quote=F, sep="\t");


##############################################################################

cat("Done.\n");

print(warnings());
q(status=0);
