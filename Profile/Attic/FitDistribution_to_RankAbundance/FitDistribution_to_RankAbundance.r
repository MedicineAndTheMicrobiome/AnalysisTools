#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
        "input_file", "i", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage=paste (
		"\nUsage:\n\t", script_name, "\n\t\t-i <Input summary_table.xls FileName>\n\n",
		"This script will try to fit the following distributions to the rank abudance curve.\n",
		"	1.) Pareto 2 (Lomax)\n",
		"	2.) Lognormal\n",
		"	3.) Gamma\n",
		"	4.) Exponential\n",
		"\n");

if(!length(opt$input_file)){
	cat(usage);
	quit(status=-1);
}

InputFileName=opt$input_file;

###############################################################################
# Main program loop

RAPlot= paste(InputFileName, ".fitted_rank_abund.pdf", sep="")
EPFile= paste(InputFileName, ".fitted_rank_abund.param.csv", sep="")

cat("\n")
cat("                 Input File Name: ", InputFileName, "\n")
cat(" Fitted Rank Abundancy Plots PDF: ", RAPlot, "\n")
cat("            Estimated Parameters: ", EPFile, "\n")

###############################################################################
###############################################################################

# Load data
A<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, row.names=1));
#cat("Original Matrix:\n")
#print(A);

# Extract out only the counts (ignore totals)
count_mat=A[,2:ncol(A)];
num_categories=ncol(count_mat);
num_samples=nrow(count_mat);

cat("Num Samples: ", num_samples, "\n");
cat("Num Taxa: ", num_categories, "\n");

###############################################################################

# Compute shortened names
long_names=colnames(count_mat);
short_names=character(num_categories);
for(i in 1:num_categories){
	taxonomy=unlist(strsplit(long_names[i], " "));
	short_names[i]=taxonomy[length(taxonomy)];
}

# Assign shortened names to count matrix
colnames(count_mat)=short_names;
#print(count_mat);

sample_names=rownames(count_mat);

###############################################################################

# Compute normalized
total=apply(count_mat, 1, sum);
#print(total);
norm_mat=count_mat/total;
#print(norm_mat);

###############################################################################

# rank by abundance
ranked_mat=matrix(0, ncol=num_categories, nrow=num_samples);
sortix_mat=matrix(0, ncol=num_categories, nrow=num_samples);
for(i in 1:num_samples){
	if(total[i]>0){
		sinfo=sort(norm_mat[i,], decreasing=TRUE, index.return=TRUE, method="shell");
		ranked_mat[i,]=sinfo$x;
		sortix_mat[i,]=sinfo$ix;
	}
}
#print(ranked_mat);
#print(sortix_mat);

###############################################################################

pdf(RAPlot, width=11,height=8.5)

library(stats4);

# Log normal
loglikfun_lognorm=function(param, data){
        mean=param[1];
        sd=param[2];
        logl=sum(dlnorm(x=data, mean, sd, log=TRUE));
        return(-logl);
}

if(!library(actuar, logical.return=T)){
	dpareto=function(x, shape, scale, log=FALSE){
		a=shape;
		xm=scale;

		y=rep(0, length(x));
		y[x<xm]=0
		y[x>=xm]=(a*(xm^a))/((x+xm)^(a+1));

		if(log){
			return(log(y));
		}else{
			return(y);
		}
	}
}

#overide exp
dhcauchy=function(x, shape, log=FALSE){
	y=2*dcauchy(abs(x),0,shape);
	if(log){
		return(log(y));
	}else{
		return(y);
	}
}

# Pareto
loglikfun_pareto=function(param, data){
	a=param[1]; # Shape
	s=param[2]; # Scale
	logl=sum(dpareto(x=data, a, s, log=TRUE));
	return(-logl);
}

# Gamma 
loglikfun_gamma=function(param,data){
	shape=param[1];
	scale=param[2];
	logl=sum(dgamma(x=data, shape, scale, log=TRUE));
	return(-logl);
}

# Exp
loglikfun_exp=function(param,data){
	lambda=param[1];
	logl=sum(dexp(x=data, lambda, log=TRUE));
	return(-logl);
}

# Half cauchy 
loglikfun_hcauchy=function(param,data){
	scale=param[1];
	logl=sum(dhcauchy(x=data, scale, log=TRUE));
	return(-logl);
}

chisq=function(exp_prob, obs){
	exp_prob=exp_prob/sum(exp_prob);
	res=chisq.test(x=obs, p=exp_prob);
	#print(res$expected):
	#print(res$observed);
	#print("--------------------------");
	return(res$p.value);
}


resample_size=2000;
lognorm_param=matrix(0, nrow=num_samples, ncol=2);
pareto_param =matrix(0, nrow=num_samples, ncol=2);
exp_param    =matrix(0, nrow=num_samples, ncol=2);
gamma_param  =matrix(0, nrow=num_samples, ncol=2);
hcauchy_param =matrix(0, nrow=num_samples, ncol=2);

lnorm_pvalue=rep(0, num_samples);
pareto_pvalue=rep(0, num_samples);
gamma_pvalue=rep(0, num_samples);
exp_pvalue=rep(0, num_samples);
hcauchy_pvalue=rep(0, num_samples);

for(i in 1:num_samples){

	cat("Working on sample: ", i, " of ", num_samples, "\n");

	if(total[i]>0){
		obs=sample(0:(num_categories-1), resample_size, replace=TRUE, prob=ranked_mat[i,]);
		obs=obs+0.5;

		range=(0:(num_categories-1))+.5;

		sorted_counts=count_mat[i,sortix_mat[i,]];

		seed_param=c(1, 1);
		fit=optim(seed_param, loglikfun_lognorm, data=obs); # x in (0, inf)
		lognorm_param[i,]=fit$par;
		lnorm_pvalue[i]=chisq(dlnorm(range, fit$par[1], fit$par[2]), sorted_counts);

		seed_param=c(1, 1);
		fit=optim(seed_param, loglikfun_pareto, data=obs); # x in [Xm, inf), but Xm > 0
		pareto_param[i,]=fit$par;
		pareto_pvalue[i]=chisq(dpareto(range, fit$par[1], fit$par[2]), sorted_counts);

		seed_param=c(1, 2);
		fit=optim(seed_param, loglikfun_gamma, data=obs); # x in [0, inf)
		gamma_param[i,]=fit$par;
		gamma_pvalue[i]=chisq(dgamma(range, fit$par[1], fit$par[2]), sorted_counts);

		seed_param=c(1);
		fit=optim(seed_param, loglikfun_exp, data=obs); # x in [0, inf)
		exp_param[i,1]=fit$par;
		exp_pvalue[i]=chisq(dexp(range, fit$par[1]), sorted_counts);

		seed_param=c(1);
		fit=optim(seed_param, loglikfun_hcauchy, data=obs); # x in [0, inf)
		hcauchy_param[i,1]=fit$par;
		hcauchy_pvalue[i]=chisq(dhcauchy(range, fit$par[1]), sorted_counts);

	}
}

warnings();
###############################################################################

par(oma=c(5.1,.5, 1,.5));


rmsd=function(x,y){
	n=length(x);
	return(sqrt(sum((x-y)^2)/n));
}


for(i in 1:num_samples){
	
	cat("Plotting: ", i, " of ", num_samples, "\n");

	if(total[i]>0){
		centers=barplot(ranked_mat[i,], plot=FALSE);
		
		nonzero_idx=sum(ranked_mat[i,]>0);

		obs=0:(num_categories-1);
		obs=obs+.5;
		estimated_lnorm=dlnorm(obs, lognorm_param[i,1], lognorm_param[i,2]);
		estimated_pareto=dpareto(obs, pareto_param[i,1], pareto_param[i,2]);
		estimated_gamma=dgamma(obs, gamma_param[i,1], gamma_param[i,2]);
		estimated_exp=dexp(obs, exp_param[i,1]);
		estimated_hcauchy=dexp(obs, hcauchy_param[i,1]);

		rmsd_lnorm=rmsd(estimated_lnorm, ranked_mat[i,]);
		rmsd_pareto=rmsd(estimated_pareto, ranked_mat[i,]);
		rmsd_gamma=rmsd(estimated_gamma, ranked_mat[i,]);
		rmsd_exp=rmsd(estimated_exp, ranked_mat[i,]);
		rmsd_hcauchy=rmsd(estimated_hcauchy, ranked_mat[i,]);
		

		xlimit=c(0,centers[max(nonzero_idx)+1]);
		ylimit=c(0, max(c(ranked_mat[i,1], estimated_lnorm, estimated_pareto, estimated_gamma, estimated_exp, estimated_hcauchy))*1.2);

		centers=barplot(
			ranked_mat[i,], 
			ylim=ylimit, xlim=xlimit, 
			main = sample_names[i],
			names.arg=short_names[sortix_mat[i,]],
			las=2
		);

		lines(centers, estimated_lnorm, type="l", pch=1, lwd=3, col="green");
		lines(centers, estimated_pareto, type="l", pch=1, lwd=3, col="blue");
		lines(centers, estimated_gamma, type="l", pch=1, lwd=3, col="red");
		lines(centers, estimated_exp, type="l", pch=1, lwd=3, col="purple");
		lines(centers, estimated_hcauchy, type="l", pch=1, lwd=3, col="orange");

		legend_info=c(
			sprintf("lognormal (mean=%3.3f, variance=%3.3f) p-val=%g rmsd=%g",  
				lognorm_param[i,1], lognorm_param[i,2], lnorm_pvalue[i], rmsd_lnorm),	
			sprintf("pareto 2 (Lomax) (shape=%3.3f, scale=%3.3f) p-val=%g rmsd=%g",       
				pareto_param[i,1], pareto_param[i,2], pareto_pvalue[i], rmsd_pareto),	
			sprintf("gamma (shape=%3.3f, scale=%3.3f) p-val=%g rmsd=%g",        
				gamma_param[i,1], gamma_param[i,2], gamma_pvalue[i], rmsd_gamma),
			sprintf("exponential (lambda=%3.3f) p-val=%g rmsd=%g",              
				exp_param[i,1], exp_pvalue[i], rmsd_exp),
			sprintf("half-cauchy (scale=%3.3f) p-val=%g rmsd=%g",
				hcauchy_param[i,1], hcauchy_pvalue[i], rmsd_hcauchy)
		);

		legend(mean(xlimit)*.75, mean(ylimit), legend=legend_info, fill=c("green", "blue", "red", "purple", "orange"), cex=.8); 
	}
}

dev.off();

###############################################################################

# Output indices

fc=file(EPFile, "w")

header=paste(
	"sample_id",
	"lnorm-mean","lnorm-var",
	"pareto-shape", "pareto-scale",
	"gamma-shape", "gamma-scale",
	"expon-lambda",
	"halfcauchy-scale",
	sep=","
);

cat(header, "\n", sep="", file=fc);

for(i in 1:num_samples){

	dataline=paste(sample_names[i], 
		lognorm_param[i,1], lognorm_param[i,2], 
		pareto_param[i,1], pareto_param[i,2], 
		gamma_param[i,1], gamma_param[i,2],
		exp_param[i,1],
		hcauchy_param[i,1],
		sep=","
	);

	cat(dataline, "\n", sep="", file=fc);

}


close(fc);

###############################################################################

writeLines("Done.\n")

q(status=0)
