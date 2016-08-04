#!/usr/bin/env Rscript

###############################################################################

EXTRAP_LIMIT=3000;

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 2, "character",
	"max_samples", "s", 2, "numeric",
	"extrapolation", "e", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n   ", script_name, "\n",
	"	-i <input rarefaction median file>\n",
	"	[-o <output filename root>\n",
	"	[-s <Num samples to use, ie max, default=use all>]\n",
	"	[-e <How many samples to extrapolate to, default=", EXTRAP_LIMIT, ">]\n",
	"\n",
	"This script will estimate the number of taxa left to be sampled given a rarefaction curve.\n",
	"\n", sep="");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
SampleSubset=opt$max_samples;

if(length(opt$output_file)){
	OutputFileNameRoot=opt$output_file;
}else{
	OutputFileNameRoot=opt$input_file;
}

if(length(opt$extrapolation)){
	EXTRAP_LIMIT=as.numeric(opt$extrapolation);
}


###############################################################################
# Load counts from file

cat("Working on ", InputFileName, "\n", sep="");

# Load summary_table.xls
mat<-as.matrix(read.table(InputFileName, sep="\t", check.names=FALSE, row.names=1, fill=TRUE))

print(mat);
numSamples=nrow(mat);
numCols=ncol(mat);
sample_names=rownames(mat);

# Clean up sample names
for(i in 1:numSamples){
	split_res=strsplit(sample_names[i], "/");
	sample_names[i]=split_res[[1]][length(split_res[[1]])];
	sample_names[i]=gsub("\\.summary_table\\.xls", "", sample_names[i]);
}

cat("Input Sizes:\n");
cat("Num Samples: ", numSamples, "\n", sep="");
cat("Num Columns: ", numCols, "\n", sep="");


# Add 0 to EDF to make a formal CDF
EDF=cbind(0, mat);
colnames(EDF)=0:numCols;
rownames(EDF)=substr(rownames(EDF), 1, 10);

# Take sample subset if requested
if(length(SampleSubset)>0){
	cat("Asked for ", SampleSubset, " samples.\n", sep="");
	if(SampleSubset+1 > numCols){
		SampleSubset = numCols;
	}
	cat("Only using data from 1 to ", SampleSubset, "\n", sep="");
	if(nrow(EDF)==1){
		EDF=matrix(EDF[,1:(SampleSubset+1)], nrow=1);
	}else{
		EDF=EDF[,1:(SampleSubset+1)];
	}
}else{
	cat("Utilizing full dataset.\n");
}
cat("EDF:\n");
print(EDF);

rmsd=function(obs, exp, norm_fact){
	n=length(obs);
	obs=obs*norm_fact;
	exp=exp*norm_fact;
	return(sqrt(sum(((obs-exp))^2)/n));	
}

# We need this for the pareto function
library(actuar); 

####################################################################

fit_gamma_to_raref=function(param,data){
	norm_fact=param[1];
	shape=param[2];
	scale=param[3];
	
	# Compute observed
	norm_data=data/norm_fact;
	
	# Compute expected
	num_obs=length(data);
	x=0:(num_obs-1);
	expected=pgamma(x, shape, scale);
	
	return(rmsd(norm_data, expected, norm_fact));
}

fit_lnorm_to_raref=function(param,data){
	norm_fact=param[1];
	mean=param[2];
	variance=param[3];
	
	# Compute observed
	norm_data=data/norm_fact;
	
	# Compute expected
	num_obs=length(data);
	x=0:(num_obs-1);
	expected=plnorm(x, mean, variance);
	
	return(rmsd(norm_data, expected, norm_fact));
}

fit_pareto_to_raref=function(param,data){
	norm_fact=param[1];
	shape=param[2];
	scale=param[3];
	
	# Compute observed
	norm_data=data/norm_fact;
	
	# Compute expected
	num_obs=length(data);
	x=0:(num_obs-1);
	expected=ppareto(x, shape, scale);
	
	return(rmsd(norm_data, expected, norm_fact));
}

pfrechet=function(x, shape, scale, location=0){
	#y=exp(-exp(-(x-location)/scale)); # Gumbel
	#y=exp(-(-(x-location)/scale)^shape); # Reversed Weibull
	y=exp(-((x-location)/scale)^-shape); # Frechet	

	return(y);
}

dfrechet=function(x, shape, scale, location=0){
	a=shape;
	b=scale;
	g=location;
	bxg=b/(x-g);
	y=(a/b)*((bxg^(a+1))*exp(-(bxg)^a));
	return(y);
}

qfrechet=function(x, shape, scale, location=0){
	a=shape;
	m=location;
	s=scale;
	y=(s*(-log(x))^(-1/a))+m;
	return(y);
}



fit_frechet_to_raref=function(param,data){
	norm_fact=param[1];
	shape=param[2];
	scale=param[3];
	#location=param[4];
	location=0;
	
	# Compute observed
	norm_data=data/norm_fact;
	
	# Compute expected
	num_obs=length(data);
	x=0:(num_obs-1);
	expected=pfrechet(x, shape, scale, location);
	
	return(rmsd(norm_data, expected, norm_fact));
}

fit_exp_to_raref=function(param,data){
	norm_fact=param[1];
	lambda=param[2];
	
	# Compute observed
	norm_data=data/norm_fact;
	
	# Compute expected
	num_obs=length(data);
	x=0:(num_obs-1);
	expected=pexp(x, lambda);
	
	return(rmsd(norm_data, expected, norm_fact));
}

####################################################################

pdf(paste(OutputFileNameRoot, ".fitted.pdf", sep=""), height=8.5, width=11);

out_labels=c(
	"SampleName", "NumDonors", "NumTaxaSampled", "DistType", "RMSD", "EstMaxTaxa", "PercTaxaSeen",
	"PercIncOfTotalTaxa", "NumIncOfTaxaExp",
	"75p", "80p", "85p", "90p", "95p", "97.5p", "99p", 
	"p1name", "p1", "p2name", "p2");

NUM_DIST=4;
out_estimates=matrix("", nrow=numSamples*NUM_DIST, ncol=length(out_labels));

best_optim=function(seed, fit_fun, edf){
	fit=optim(seed, fit_fun, data=edf, method="Nelder-Mead");

	best_fit=fit;
	best_score=fit$value;
	cat("*************************************************\n");
	cat("Initial found parameters: ", paste(best_fit$par, sep=","), "\n");
	cat("Initial score: ", best_score, "\n");

	search_limits=6;
	p1=seq(fit$par[1]/search_limits, fit$par[1]*search_limits, length.out=4);
	p2=seq(fit$par[2]/search_limits, fit$par[2]*search_limits, length.out=4);

	if(length(seed)==3){
		p3=seq(fit$par[3]/search_limits, fit$par[3]*search_limits, length.out=4);
	}else{
		p3=c();
	}	

	for(i in 1:length(p1)){
		for(j in 1:length(p2)){
			for(k in 1:length(p3)){

				if(length(seed)==3){
					seed=c(p1[i], p2[j], p3[k]);
				}else{
					seed=c(p1[i], p2[j]);
				}

				fit_ijk=optim(seed, fit_fun, data=edf, method="Nelder-Mead");

				if(best_score>fit_ijk$value){
					cat("\t*** Found better score! ***\n");
					cat("Best Seed:", paste(seed, collapse=","), "\n", sep="");
					cat("Best Parameters:", paste(fit_ijk$par, collapse=","), "\n", sep="");
					cat("Best Score: ", fit_ijk$value, "\n", sep="");
					best_fit=fit_ijk;
					best_score=fit_ijk$value;
				}else{
					#cat("\tScore not as good.\n");
				}
			}
		}
	}

	return(best_fit);
}

exterp=10;
layout_matrix=matrix(c(1,1,1,2,3,3),nrow=2, byrow=T);
print(layout_matrix);
layout(layout_matrix);

dist_colors=c("blue", "red", "green", "purple", "orange");

for(i in 1:numSamples){

	cat("\n\nWorking on: ", sample_names[i], "\n");	

	# Estimate max taxa as max on curve
	edf=EDF[i,is.finite(EDF[i,])];
	max=3*max(edf);
	numDonors=length(edf)-1;
	numTaxaSampled=max(edf);
	cat("Num Taxa Sampled: ", numTaxaSampled, "\n");
	
	# Perform search on parameters and max num taxa	
	cat("Fitting log normal...\n");
	seed=c(max, 1, 1);
	fit_lnorm=best_optim(seed, fit_lnorm_to_raref, edf);

	cat("Fitting gamma...\n");
	seed=c(max, 1, 1);
	fit_gamma=best_optim(seed, fit_gamma_to_raref, edf);

	cat("Fitting pareto...\n");
	seed=c(max, 1, 1);
	fit_pareto=best_optim(seed, fit_pareto_to_raref, edf);

	cat("Fitting frechet...\n");
	seed=c(max, .5, .5);
	fit_frechet=best_optim(seed, fit_frechet_to_raref, edf);

	#seed=c(max, 1);
	#fit_exp=best_optim(seed, fit_exp_to_raref, edf);

	#-----------------------------------------------------------------------------------------
	# Plot the Empirical Distribution Function
	disp_range=0:(numDonors+exterp);
	x=disp_range;

	# Compute estimated based on discovered parameters/max num taxa
	lnorm_estim=plnorm(x,fit_lnorm$par[2], fit_lnorm$par[3])*fit_lnorm$par[1];
	gamma_estim=pgamma(x,fit_gamma$par[2], fit_gamma$par[3])*fit_gamma$par[1];
	pareto_estim=ppareto(x,fit_pareto$par[2], fit_pareto$par[3])*fit_pareto$par[1];
	frechet_estim=pfrechet(x,fit_frechet$par[2], fit_frechet$par[3], 0)*fit_frechet$par[1];
	#exp_estim=pexp(x,fit_exp$par[2])*fit_exp$par[1];

	xlimit=range(disp_range);
	ylimit=c(0, max(edf,gamma_estim,lnorm_estim,pareto_estim,frechet_estim));
	#ylimit=c(0, max(edf,gamma_estim,lnorm_estim,pareto_estim,frechet_estim, exp_estim));

	plot(0:numDonors, edf, 
		main=sample_names[i],
		ylim=ylimit,
		xlim=xlimit,
		xlab="Num Samples",
		ylab="Num Discovered Taxa"
	);
	mtext(sprintf("Num Samples: %i", numDonors), side=3, line=0, cex=.6);

	# Plot the 3 estimates
	lines(disp_range, lnorm_estim, col=dist_colors[1]);
	lines(disp_range, gamma_estim, col=dist_colors[2]);
	lines(disp_range, pareto_estim, col=dist_colors[3]);
	lines(disp_range, frechet_estim, col=dist_colors[4]);
	#lines(disp_range, exp_estim, col="orange");
	
	# Plot the legend
	legend_info=c(
		sprintf("lognormal (mean=%3.3g, var=%3.3g), rmsd=%3.4g, estim taxa: %i",
			fit_lnorm$par[2], fit_lnorm$par[3], fit_lnorm$value, as.integer(fit_lnorm$par[1])),
		sprintf("gamma (shape=%3.3g, scale=%3.3g), rmsd=%3.4g, estim taxa: %i", 
			fit_gamma$par[2], fit_gamma$par[3], fit_gamma$value, as.integer(fit_gamma$par[1])),
		sprintf("pareto (shape=%3.3g, scale=%3.3g), rmsd=%3.4g, estim taxa: %i",
			fit_pareto$par[2], fit_pareto$par[3], fit_pareto$value, as.integer(fit_pareto$par[1])),
		sprintf("frechet (shape=%3.3g, scale=%3.3g), rmsd=%3.4g, estim taxa: %i",
			fit_frechet$par[2], fit_frechet$par[3], fit_frechet$value, as.integer(fit_frechet$par[1]))
		#sprintf("exp (lambda=%3.3g), rmsd=%3.4g, estim taxa: %i",
		#	fit_exp$par[2], fit_exp$value, as.integer(fit_exp$par[1]))

	);
	legend(median(xlimit), median(ylimit), legend=legend_info, fill=dist_colors, cex=1.2);

	#-----------------------------------------------------------------------------------------
	# Plot the PDF
	x=1:as.integer(numDonors*1.15);
	dlnorm=dlnorm(x, fit_lnorm$par[2], fit_lnorm$par[3]);
	dgamma=dgamma(x, fit_gamma$par[2], fit_gamma$par[3]);
	dpareto=dpareto(x, fit_pareto$par[2], fit_pareto$par[3]);
	dfrechet_v=dfrechet(x, fit_frechet$par[2], fit_frechet$par[3]);
	#dexp=dexp(x, fit_exp$par[2]);

	pdf_values=c(dlnorm,dgamma,dpareto,dfrechet_v);
	#pdf_values=c(dlnorm,dgamma,dpareto,dfrechet_v,dexp);
	pdf_values=pdf_values[!is.nan(pdf_values)];
	ylimit=max(pdf_values)*1.15;

	plot(0, 0,
		main="PDF",
		xlim=c(0,numDonors+1),
		ylim=c(0,ylimit),
		xlab="Num Samples",
		ylab="Increase in Porportion of Est Taxa",
		type="n"
	);

	lines(x, dlnorm, col=dist_colors[1]);
	lines(x, dgamma, col=dist_colors[2]);
	lines(x, dpareto, col=dist_colors[3]);
	lines(x, dfrechet_v, col=dist_colors[4]);
	#lines(x, dexp, col=dist_colors[5]);
	abline(v=numDonors, col="grey");

	#-----------------------------------------------------------------------------------------
	# Plot CDF to EXTRAP_LIMIT

	# Plot the Empirical Distribution Function
	disp_range=0:EXTRAP_LIMIT;
	x=disp_range;

	# Compute estimated based on discovered parameters/max num taxa
	lnorm_estim=plnorm(x,fit_lnorm$par[2], fit_lnorm$par[3])*fit_lnorm$par[1];
	gamma_estim=pgamma(x,fit_gamma$par[2], fit_gamma$par[3])*fit_gamma$par[1];
	pareto_estim=ppareto(x,fit_pareto$par[2], fit_pareto$par[3])*fit_pareto$par[1];
	frechet_estim=pfrechet(x, fit_frechet$par[2], fit_frechet$par[3])*fit_frechet$par[1];
	#exp_estim=pexp(x, fit_exp$par[2])*fit_exp$par[1];

	xlimit=c(0, max(disp_range)*1.05);
	ylimit=c(0, max(edf,gamma_estim,lnorm_estim,pareto_estim,frechet_estim));
	#ylimit=c(0, max(edf,gamma_estim,lnorm_estim,pareto_estim,frechet_estim,exp_estim));
	
	plot(0:numDonors, edf, 
		main=sprintf("Rarefaction Extrapolated to %i Samples", EXTRAP_LIMIT),
		ylim=ylimit,
		xlim=xlimit,
		xlab="Num Samples",
		ylab="Num Discovered Taxa"
	);

	# Plot the estimates
	lines(disp_range, lnorm_estim, col=dist_colors[1]);
	lines(disp_range, gamma_estim, col=dist_colors[2]);
	lines(disp_range, pareto_estim, col=dist_colors[3]);
	lines(disp_range, frechet_estim, col=dist_colors[4]);
	#lines(disp_range, exp_estim, col=dist_colors[5]);

	perc_taxa_lnorm=plnorm(EXTRAP_LIMIT, fit_lnorm$par[2], fit_lnorm$par[3]); 
	perc_taxa_gamma=pgamma(EXTRAP_LIMIT, fit_gamma$par[2], fit_gamma$par[3]); 
	perc_taxa_pareto=ppareto(EXTRAP_LIMIT, fit_pareto$par[2], fit_pareto$par[3]); 
	perc_taxa_frechet=pfrechet(EXTRAP_LIMIT, fit_frechet$par[2], fit_frechet$par[3]); 
	#perc_taxa_exp=pexp(EXTRAP_LIMIT, fit_exp$par[2]); 

	text(EXTRAP_LIMIT, perc_taxa_lnorm*fit_lnorm$par[1], sprintf("%3.1f%%", perc_taxa_lnorm*100), pos=4, col="grey33", cex=.9);
	text(EXTRAP_LIMIT, perc_taxa_gamma*fit_gamma$par[1], sprintf("%3.1f%%", perc_taxa_gamma*100), pos=4, col="grey33", cex=.9);
	text(EXTRAP_LIMIT, perc_taxa_pareto*fit_pareto$par[1], sprintf("%3.1f%%", perc_taxa_pareto*100), pos=4, col="grey33", cex=.9);
	text(EXTRAP_LIMIT, perc_taxa_frechet*fit_frechet$par[1], sprintf("%3.1f%%", perc_taxa_frechet*100), pos=4, col="grey33", cex=.9);
	#text(EXTRAP_LIMIT, perc_taxa_exp*fit_exp$par[1], sprintf("%3.1f%%", perc_taxa_exp*100), pos=4, col="grey33", cex=.9);

	approachmax_lnorm=(fit_lnorm$par[1]-1)/fit_lnorm$par[1];
	approachmax_gamma=(fit_gamma$par[1]-1)/fit_gamma$par[1];
	approachmax_pareto=(fit_pareto$par[1]-1)/fit_pareto$par[1];
	approachmax_frechet=(fit_frechet$par[1]-1)/fit_frechet$par[1];
	#approachmax_exp=(fit_exp$par[1]-1)/fit_exp$par[1];

	legend_info=c(
		sprintf("at %3.2f%%, NumSampl=%3.2e ", approachmax_lnorm*100, qlnorm(approachmax_lnorm, fit_lnorm$par[2], fit_lnorm$par[3])),
		sprintf("at %3.2f%%, NumSampl=%3.2e ", approachmax_gamma*100, qgamma(approachmax_gamma, fit_gamma$par[2], fit_gamma$par[3])),
		sprintf("at %3.2f%%, NumSampl=%3.2e ", approachmax_pareto*100,qpareto(approachmax_pareto, fit_pareto$par[2], fit_pareto$par[3])),
		sprintf("at %3.2f%%, NumSampl=%3.2e ", approachmax_frechet*100,qfrechet(approachmax_frechet, fit_frechet$par[2], fit_frechet$par[3]))
		#sprintf("at %3.2f%%, NumSampl=%3.2e ", approachmax_exp*100,qexp(approachmax_exp, fit_exp$par[2]))
	);

	legend(median(xlimit)/2, median(ylimit)/2, legend=legend_info, fill=dist_colors, cex=1.2);

	#-----------------------------------------------------------------------------------------
	# Output parameters to file

	percentiles=c(.75, .80, .85, .90, .95, .975, .99);

	# Save log normal values
	k=((i-1)*NUM_DIST)+1;
	out_estimates[k,1]=sample_names[i];
	out_estimates[k,2]=numDonors;
	out_estimates[k,3]=numTaxaSampled;
	out_estimates[k,4]="lnorm";
	out_estimates[k,5]=sprintf("%5.5g", fit_lnorm$value);
	out_estimates[k,6]=floor(fit_lnorm$par[1]);
	out_estimates[k,7]=sprintf("%3.2f", 100*numTaxaSampled/floor(fit_lnorm$par[1]));
	out_estimates[k,8]=dlnorm(numDonors+1, fit_lnorm$par[2], fit_lnorm$par[3])*100;
	out_estimates[k,9]=(as.numeric(out_estimates[k,8])/100*fit_lnorm$par[1]);
	max_estim=floor(qlnorm(approachmax_lnorm, fit_lnorm$par[2], fit_lnorm$par[3]));
	for(p in 1:length(percentiles)){
		perc=percentiles[p];
		estim=floor(qlnorm(perc, fit_lnorm$par[2], fit_lnorm$par[3]));
		if(estim<max_estim){
			out_estimates[k, 9+p]=estim;
		}else{
			out_estimates[k, 9+p]=max_estim;
		}
	}
	out_estimates[k,17]="mean";
	out_estimates[k,18]=sprintf("%5.5g",fit_lnorm$par[2]);
	out_estimates[k,19]="variance";
	out_estimates[k,20]=sprintf("%5.5g",fit_lnorm$par[3]);

	# Save gamma values
	k=k+1;
	out_estimates[k,1]=sample_names[i];
	out_estimates[k,2]=numDonors;
	out_estimates[k,3]=numTaxaSampled;
	out_estimates[k,4]="gamma";
	out_estimates[k,5]=sprintf("%5.5g", fit_gamma$value);
	out_estimates[k,6]=floor(fit_gamma$par[1]);
	out_estimates[k,7]=sprintf("%3.2f", 100*numTaxaSampled/floor(fit_gamma$par[1]));
	out_estimates[k,8]=dgamma(numDonors+1, fit_gamma$par[2], fit_gamma$par[3])*100;
	out_estimates[k,9]=(as.numeric(out_estimates[k,8])/100*fit_gamma$par[1]);
	max_estim=floor(qgamma(approachmax_gamma, fit_gamma$par[2], fit_gamma$par[3]));
	for(p in 1:length(percentiles)){
		perc=percentiles[p];
		estim=floor(qgamma(perc, fit_gamma$par[2], fit_gamma$par[3]));
		if(estim<max_estim){
			out_estimates[k, 9+p]=estim;
		}else{
			out_estimates[k, 9+p]=max_estim;
		}
	}
	out_estimates[k,17]="shape";
	out_estimates[k,18]=sprintf("%5.5g",fit_gamma$par[2]);
	out_estimates[k,19]="scale";
	out_estimates[k,20]=sprintf("%5.5g",fit_gamma$par[3]);

	# Save pareto values
	k=k+1;
	out_estimates[k,1]=sample_names[i];
	out_estimates[k,2]=numDonors;
	out_estimates[k,3]=numTaxaSampled;
	out_estimates[k,4]="pareto";
	out_estimates[k,5]=sprintf("%5.5g", fit_pareto$value);
	out_estimates[k,6]=floor(fit_pareto$par[1]);
	out_estimates[k,7]=sprintf("%3.2f", 100*numTaxaSampled/floor(fit_pareto$par[1]));
	out_estimates[k,8]=dpareto(numDonors+1, fit_pareto$par[2], fit_pareto$par[3])*100;
	out_estimates[k,9]=(as.numeric(out_estimates[k,8])/100*fit_pareto$par[1]);
	max_estim=floor(qpareto(approachmax_pareto, fit_pareto$par[2], fit_pareto$par[3]));
	for(p in 1:length(percentiles)){
		perc=percentiles[p];
		estim=floor(qpareto(perc, fit_pareto$par[2], fit_pareto$par[3]));
		if(estim<max_estim){
			out_estimates[k, 9+p]=estim;
		}else{
			out_estimates[k, 9+p]=max_estim;
		}
	}
	out_estimates[k,17]="shape";
	out_estimates[k,18]=sprintf("%5.5g",fit_pareto$par[2]);
	out_estimates[k,19]="scale";
	out_estimates[k,20]=sprintf("%5.5g",fit_pareto$par[3]);

	# Save frechet values
	k=k+1;
	out_estimates[k,1]=sample_names[i];
	out_estimates[k,2]=numDonors;
	out_estimates[k,3]=numTaxaSampled;
	out_estimates[k,4]="frechet";
	out_estimates[k,5]=sprintf("%5.5g", fit_frechet$value);
	out_estimates[k,6]=floor(fit_frechet$par[1]);
	out_estimates[k,7]=sprintf("%3.2f", 100*numTaxaSampled/floor(fit_frechet$par[1]));
	out_estimates[k,8]=dfrechet(numDonors+1, fit_frechet$par[2], fit_frechet$par[3])*100;
	out_estimates[k,9]=(as.numeric(out_estimates[k,8])/100*fit_frechet$par[1]);
	max_estim=floor(qfrechet(approachmax_frechet, fit_frechet$par[2], fit_frechet$par[3]));
	for(p in 1:length(percentiles)){
		perc=percentiles[p];
		estim=floor(qfrechet(perc, fit_frechet$par[2], fit_frechet$par[3]));
		if(estim<max_estim){
			out_estimates[k, 9+p]=estim;
		}else{
			out_estimates[k, 9+p]=max_estim;
		}
	}
	out_estimates[k,17]="shape";
	out_estimates[k,18]=sprintf("%5.5g",fit_frechet$par[2]);
	out_estimates[k,19]="scale";
	out_estimates[k,20]=sprintf("%5.5g",fit_frechet$par[3]);
	
}

fit_out_fh=file(paste(OutputFileNameRoot, ".fit_and_estimates.csv", sep=""), "wt");
cat(file=fit_out_fh, c(out_labels, "\n"), sep=",");

for(i in 1:(NUM_DIST*numSamples)){
	cat(file=fit_out_fh, c(out_estimates[i,], "\n"), sep=",");
}
close(fit_out_fh);

cat("Done.\n");
warnings();
q(status=0);
