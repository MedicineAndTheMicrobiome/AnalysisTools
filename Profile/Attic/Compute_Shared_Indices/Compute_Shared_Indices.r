#!/usr/local/bin/Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	
	"bootstrap_iterations", "b", 2, "integer",
	"bootstrap_resample_size", "n", 2, "integer",

	"job_index", "j", 2, "integer",
	"total_jobs", "J", 2, "integer"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];
base_name=dirname(script_name);
source(paste(base_name, '/WeightedRankDifference.r', sep=""));

if(length(opt$input_file)==0){
	usage <- paste (
		"\nUsage:\n\t", script_name, "\n",
	
		"		-i <input summary_table.xls file>\n",
		"\n",
		"	Bootstrap related:\n",
		"		[-b <bootstrap iterations, 40 for 95% CI>]\n",
		"		[-n <bootstrap resample size, 2000>]\n",
		"\n",
		"	Grid related:\n",
		"		[-j <task id, usually $SGE_TASK_ID>]\n",
		"		[-J <total tasks, usually $SGE_TASK_LAST>]\n",
		"\n",

		"Computes a variety of beta-diversity indices pairwise.\n",
		"Note for every sample in your input, an output file will be generated.\n",
		"\n",
		"\n")

	writeLines(usage)
	writeLines("Input FileName not defined.\n")
	quit(status=0)
}

###############################################################################
# Bootstrapping parameters

do_bootstrap=FALSE;
RESAMPLE_SIZE=2000;
NUM_BOOTSTRAPS=40;

if(!is.null(opt$bootstrap_iterations) & !is.null(opt$bootstrap_resample_size)){
	do_bootstrap=TRUE;
	RESAMPLE_SIZE=opt$bootstrap_resample_size;
	NUM_BOOTSTRAPS=opt$bootstrap_iterations;
}

if(sum(c(is.null(opt$bootstrap_iterations),is.null(opt$bootstrap_resample_size)))==1){
	cat("Error, you need to specify both -b and -n to perform bootstrapping.\n");	
	q(status=-1)
}

if(do_bootstrap){
	cat("Bootstrapping is ON.\n");
	cat("  BS Trials: ", NUM_BOOTSTRAPS, "\n");
	cat("  BS Resample Size: ", RESAMPLE_SIZE, "\n");
}else{
	cat("Bootstrapping is OFF.\n");
}

###############################################################################
# Grid parameters

on_grid=FALSE;

if(!is.null(opt$job_index) & !is.null(opt$total_jobs)){
	on_grid=TRUE;	
}
if(sum(c(is.null(opt$job_index),is.null(opt$total_jobs)))==1){
	cat("Error, you need to specify both -t and -T to for grid mode.\n");	
	q(status=-1)
}

if(on_grid){
	cat("Mode: On grid parameters ON.\n");
	cat("  job index:  ", opt$job_index, "\n");
	cat("  total jobs: ", opt$total_jobs, "\n");
	opt$job_index=opt$job_index %% opt$total_jobs;
}else{
	cat("Mode: On grid parameters OFF.\n");
}

###############################################################################

compute_chaobased = function (X, Y){

	#print(X);
	#print(Y);

	n=sum(X);
	m=sum(Y);
	#cat("X count: ", n, "\n", sep="");
	#cat("Y count: ", m, "\n", sep="");
	
	# Shared OTUs
	D_12= X & Y;
	# Index of shared OTUs
	which_shared=which(D_12);

	# f_1+ = obs shared singletons, in X
	f_1plus=sum(X==1 & Y>=1);
	# f_2+ = obs shared doubletons, in X
	f_2plus=sum(X==2 & Y>=1);

	# f_+1 = obs shared singletons, in Y
	f_plus1=sum(Y==1 & X>=1)
	# f_+2 = obs shared doubletons, in Y
	f_plus2=sum(Y==2 & X>=1)
	
	#cat("X shared singletons: ", f_1plus, "\n", sep="");
	#cat("X shared doubletons: ", f_2plus, "\n", sep="");
	#cat("Y shared singletons: ", f_plus1, "\n", sep="");
	#cat("Y shared doubletons: ", f_plus2, "\n", sep="");

	# Compute observed portion
	U_observed=sum(X[which_shared])/n;
	V_observed=sum(Y[which_shared])/m;

	#cat("X observed shared (proportion): ", U_observed, "\n", sep="");
	#cat("Y observed shared (proportion): ", V_observed, "\n", sep="");

	# Compute unseen based on chao estimator
	U_unseen=	((m-1)/m) *
			 (f_plus1/(2*f_plus2))*
			sum(X[Y==1])/n;
	if(is.infinite(U_unseen)){
		U_unseen=0;
	}
	U=U_observed+U_unseen;

	V_unseen=	((n-1)/n) *
			 (f_1plus/(2*f_2plus))*
			sum(Y[X==1])/m;
	if(is.infinite(V_unseen)){
		V_unseen=0;
	}
	V=V_observed+V_unseen;
	
	#cat("X unseen (proportion): ", U_unseen, "\n", sep="");
	#cat("Y unseen (proportion): ", V_unseen, "\n", sep="");

	#cat("X combined (proportion): ", U, "\n", sep="");
	#cat("Y combined (proportion): ", V, "\n", sep="");

	jacc=(U*V)/(U+V-(U*V));
	sore=(2*U*V)/(U+V);
	
	out=c(jacc, sore);
	names(out)=c("chao jaccard", "chao sorensen");

	out;

}

###############################################################################

compute_setbased = function(X, Y){
		
	# total
	x_pres=(X>0);
	y_pres=(Y>0);

	# computed sharedness
	shared=sum(x_pres & y_pres);
	x=sum(x_pres);
	y=sum(y_pres);

	# compute indices
	jacc=shared/(x+y-shared);
	sore=(2*shared)/(x+y);

	# return indices
	out=c(jacc, sore);
	names(out)=c("set jaccard", "set sorensen");
	out;

}

###############################################################################

compute_probbased = function(X_prob, Y_prob){
	
	# compute shared
	shared=numeric(0);
	for(i in 1:length(X_prob)){
		shared[i]=min(X_prob[i], Y_prob[i]);
	}

	# compute sharedness, x_only and y_only are the same because we normalized.
	shared = sum(shared);
	x_only = 1-shared;
	y_only = 1-shared;

	# compute indices
	jacc=shared/(shared + x_only + y_only);
	sore=(2*shared)/(2*shared + x_only + y_only);

	# return indices
	out=c(jacc, sore);
	names(out)=c("prob jaccard", "prob sorensen");
	out;

}

###############################################################################

compute_eucliddist = function (X_prob, Y_prob){
	return(sqrt(sum((X_prob-Y_prob)^2)));
}

###############################################################################

compute_thetayc = function(X_prob, Y_prob){
	sum_xy=sum(X_prob*Y_prob);
	sum_sqrdiff=sum((X_prob-Y_prob)^2);
	return(sum_xy/(sum_sqrdiff + sum_xy));
}

###############################################################################

compute_morhor = function(X_prob, Y_prob){
	sum_xy=sum(X_prob*Y_prob);
	sum_sqrx=sum(X_prob^2);
	sum_sqry=sum(Y_prob^2);
	return(2*(sum_xy/(sum_sqrx+sum_sqry)));
}

###############################################################################

samples_to_summary_vector = function (samples, NumCategories){

	# Converts a list of samples, into category counts
	vector=rep(0, NumCategories);
	num_samples=length(samples);
	for(i in 1:num_samples){
		vector[samples[i]]=vector[samples[i]]+1;
	}
	vector;
}

###############################################################################

get_CI = function(x){

	# Computes the maximum confidence interval given then number of samples
	len=length(x);
	xsorted=sort(x);
	mean=mean(x)

	ci_range=(len-2)*100/len;
	# Lower bound
	lb=xsorted[2];
	# Upper bound
	up=xsorted[len-1];
	
	c(mean, lb, up);
}

###############################################################################

# Get filename of summary table
InputFileName=opt$input_file;
cat("Working on ", InputFileName, "\n", sep="");

# Load summary_table.xls
mat<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, row.names=1))

# Exclude total counts column
count_mat=mat[,2:ncol(mat)]

categories=as.vector(colnames(count_mat));
sample_names=as.vector(rownames(count_mat));
num_categories=length(categories);

#print("Original category names:\n");
#print(categories);

# Compute shorted names
short_names=character(num_categories);
for(i in 1:num_categories){
	taxonomy=unlist(strsplit(categories[i], " "));
	short_names[i]=taxonomy[length(taxonomy)];
}

#print("Shortened category names:\n");
#print(short_names);

NumSamples=nrow(count_mat);
NumCategories=ncol(count_mat);

cat("\n");
cat("Num Samples: ", NumSamples, "\n");
cat("Num Categories: ", NumCategories, "\n");
cat("\n");

indexnames=c(
	"jaccard-set", "sorensen-set", 
	"jaccard-chao", "sorensen-chao",
	"jaccard-prob", "sorensen-prob", "euclid",
	"theta-YC", 
	"morisita-horn",
	"weighted_rank_diff-p4"
);
if(!do_bootstrap){
	header_names=indexnames;
}else{
	header_names=paste(indexnames, "lb", "ub", sep=",");
}


for(i in 1:NumSamples){

	if(on_grid){
		sample_owner=i %% opt$total_jobs; # The % returns 0 to total_jobs-1, and the job_index goes from 1 to total_jobs
		if(opt$job_index != sample_owner){
			next;
		}
	}
	
	fc=file(paste(sample_names[i], ".beta_diversity",sep=""), "w");

	# Write out header lines
	outline=paste(header_names, collapse=",");
	write(paste("# SampleA,SampleB,", outline, sep=""),file=fc);

	for(j in 1:NumSamples){
		cat("Working on: [", i, "] ", sample_names[i], " vs. [", j, "] ", sample_names[j], ".\n", sep="");

		X=as.vector(count_mat[i,]);
		Y=as.vector(count_mat[j,]);

		#print(X);
		#print(Y);

		# Normalize once, instead of inside of function
		X_prob=X/sum(X);
		Y_prob=Y/sum(Y);

		if(do_bootstrap){

			idx=1:NumCategories;


			all_indices=matrix(0,nrow=NUM_BOOTSTRAPS, ncol=length(indexnames));

			for(b in 1:NUM_BOOTSTRAPS){
				Xsamples=sample(idx, RESAMPLE_SIZE, replace=TRUE, prob=X_prob);
				Ysamples=sample(idx, RESAMPLE_SIZE, replace=TRUE, prob=Y_prob);

				Xvector=samples_to_summary_vector(Xsamples, NumCategories);
				Yvector=samples_to_summary_vector(Ysamples, NumCategories);

				#print(Xvector);
				#print(Yvector);
				X_bs_prob=Xvector/sum(Xvector);
				Y_bs_prob=Yvector/sum(Yvector);

				indicesSet=compute_setbased(Xvector, Yvector);
				indicesChao=compute_chaobased(Xvector, Yvector);
				indicesProb=compute_probbased(X_bs_prob, Y_bs_prob);
				indicesEuclid=compute_eucliddist(X_bs_prob, Y_bs_prob);
				indicesThetaYC=compute_thetayc(X_bs_prob, Y_bs_prob);
				indicesMorHor=compute_morhor(X_bs_prob, Y_bs_prob);
				indicesWRD4=order_dist(X_bs_prob, Y_bs_prob, 4);

				all_indices[b,]=c(indicesSet, indicesChao, indicesProb, indicesEuclid, indicesThetaYC, indicesMorHor, indicesWRD4);
			}


			values=c();
			for(p in 1:ncol(all_indices)){
				ci=get_CI(as.vector(all_indices[,p]));	
				values=c(values,ci);
			}

			outline=paste(sample_names[i], sample_names[j], paste(values, collapse=","), sep=",");
			write(outline,file=fc);

		}else{
			indicesSet=compute_setbased(X, Y);
			indicesChao=compute_chaobased(X, Y);
			indicesProb=compute_probbased(X_prob, Y_prob);
			indicesEuclid=compute_eucliddist(X_prob, Y_prob);
			indicesThetaYC=compute_thetayc(X_prob, Y_prob);
			indicesMorHor=compute_morhor(X_prob, Y_prob);
			indicesWRD4=order_dist(X_prob, Y_prob, 4);

			values=c(indicesSet, indicesChao, indicesProb, indicesEuclid, indicesThetaYC, indicesMorHor, indicesWRD4);
			outline=paste(sample_names[i], sample_names[j], 
				paste(values, collapse=","), sep=",");
			write(outline,file=fc);
		}
	}
	close(fc);
}

writeLines("Done.\n")

warnings();
q(status=0)
