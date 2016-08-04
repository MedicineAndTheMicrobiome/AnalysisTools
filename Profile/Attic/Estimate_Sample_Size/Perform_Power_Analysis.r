#!/usr/bin/env Rscript

###############################################################################

library('getopt');
library('vegan');
#source('WeightedRankDifference.r');

ALPHA=0.05;
GAMMA=0.80;
NUM_BOOTSTRAPS=1000;
DEF_DIST="horn";
SAMPLE_SIZE_STR="5,10,15,20,25,30,40,50,70,100";

params=c(
	"input_file_A", "A", 1, "character",
	"input_file_B", "B", 1, "character",
	"output_filename_root", "o", 2, "character",
	"alpha", "a", 2, "numeric",
	"gamma", "g", 2, "numeric",
	"num_bootstraps", "b", 2, "numeric",
	"distance_type", "d", 2, "character",
	"min_effect_size", "e", 2, "numeric",
	"sample_sizes_str", "n", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"	-A <Input file A>\n",
	"	-B <Input file B>\n",
	"	[-o <Output filename root>]\n",
	"	[-a <alpha, default=", ALPHA, ">]\n",
	"	[-g <gamma, default=", GAMMA, ">]\n",
	"	[-b <num bootstraps, default=", NUM_BOOTSTRAPS, ">]\n",
	"	[-d <distance type, default=", DEF_DIST, ">]\n",
	"	[-e <induce a minimum effect size (eta^2), default=0>]\n",
	"	[-n <sample sizes string (", SAMPLE_SIZE_STR, ")>\n",
	"\n",
	"Given two cohorts that are putatively different, this\n",
	"script will calculate the number of samples necessary\n",
	"to differentiate them at various levels of critical cutoffs.\n",
	"\n",
	"Valid distances are those in vegdist, for example:\n",
	"	euclidean, manhattan, bray, kulczynski, horn, jaccard, etc.\n",
	"\n",
	"Alpha (significance level) defines the critical cutoff, which defines\n",	
	"the Gamma (power).\n",
	"\n", sep="");

if(!length(opt$input_file_A) || !length(opt$input_file_B)){
	cat(usage);
	q(status=-1);
}

InputFileNameA=opt$input_file_A;
InputFileNameB=opt$input_file_B;
FilenameRootA=gsub(".summary_table.xls", "", InputFileNameA);
FilenameRootB=gsub(".summary_table.xls", "", InputFileNameB);

if(length(opt$output_filename_root)){
	OutputRoot=opt$output_filename_root;
}else{
	OutputRoot=paste(FilenameRootA, "_vs_", FilenameRootB, sep="");
}

if(length(opt$alpha)){
	Alpha=opt$alpha;
}else{
	Alpha=ALPHA;
}

if(length(opt$gamma)){
	Gamma=opt$gamma;
}else{
	Gamma=GAMMA;
}

if(length(opt$num_bootstraps)){
	NumBootstraps=opt$num_bootstraps;
}else{
	NumBootstraps=NUM_BOOTSTRAPS;
}

if(length(opt$distance_type)){
	DistanceType=opt$distance_type;
}else{
	DistanceType=DEF_DIST;
}

if(length(opt$min_effect_size)){
	MinimumEffectSize=opt$min_effect_size;
}else{
	MinimumEffectSize=0;
}	

if(length(opt$sample_sizes_str)){
	sample_size_str=opt$sample_sizes_str;
}else{
	sample_size_str=SAMPLE_SIZE_STR;
}
sample_sizes=sort(as.numeric(strsplit(sample_size_str, ",")[[1]]));


################################################################################

cat("Filename A: ", InputFileNameA, "\n");
cat("Filename B: ", InputFileNameB, "\n");
cat("Output Root: ", OutputRoot, "\n");
cat("Alpha: ", Alpha, "\n");
cat("Gamma: ", Gamma, "\n");
cat("Num Bootstrap: ", NumBootstraps, "\n");
cat("Distance Type: ", DistanceType, "\n");
cat("Minimum Effect Size Induced: ", MinimumEffectSize, "\n");
cat("Samples Sizes to Compute: ", paste(sample_sizes, collapse=","), "\n");
cat("\n");

################################################################################

load_summary_table=function(filename){
	st=as.matrix(read.table(filename, header=TRUE, sep="\t", row.names=1, check.names=F));
	return(st[,2:ncol(st)]);
}

#-------------------------------------------------------------------------------

normalize=function(st){
	sums=apply(st, 1, sum);
	n=matrix(0,nrow=nrow(st), ncol=ncol(st));
	rownames(n)=rownames(st);
	colnames(n)=colnames(st);
	for(i in 1:nrow(st)){
		n[i,]=st[i,]/sums[i];
	}
	return(n);
}

#-------------------------------------------------------------------------------

join_summary_table=function(st1, st2){
	cnames1=colnames(st1);
	rnames1=rownames(st1);
	cnames2=colnames(st2);
	rnames2=rownames(st2);

	cnames_combined=unique(c(cnames1, cnames2));
	nrows_combined=nrow(st1)+nrow(st2);

	new_st=matrix(0, nrow=nrows_combined, ncol=length(cnames_combined));
	colnames(new_st)=cnames_combined;
	rownames(new_st)=c(rnames1, rnames2);
	
	# Copy counts over for 1
	for(i in 1:nrow(st1)){
		for(j in 1:ncol(st1)){
			new_st[rnames1[i], cnames1[j]]=st1[i,j];
		}
	}

	# Copy counts over for 2
	for(i in 1:nrow(st2)){
		for(j in 1:ncol(st2)){
			new_st[rnames2[i], cnames2[j]]=st2[i,j];
		}
	}

	return(new_st);
	
}

#-------------------------------------------------------------------------------

sum_of_squares_from_centroid=function(dist_mat, group){

	# SS to centroid == (SS among points / num points)
	# print(dist_mat);
	
	ss=0;
	glength=length(group);
	for(i in 1:glength){
		for(j in 1:i){
			ss=ss+dist_mat[group[i],group[j]]^2;
		}
	}
	return(ss/glength);

}

#-------------------------------------------------------------------------------

matrix_anova=function(dist, groupA, groupB){
	#print(groupA);
	#print(groupB);

	dist_mat=as.matrix(dist);

	# Sum of squares within A
	ssa=sum_of_squares_from_centroid(dist_mat, groupA);

	# Sum of squares within B
	ssb=sum_of_squares_from_centroid(dist_mat, groupB);

	# Sum of squares Total
	ssT=sum_of_squares_from_centroid(dist_mat, c(groupA, groupB));

	# Compute SSB and SSW
	ssW=ssa+ssb;
	ssB=ssT-ssW;

	# Compute variables for degrees of freedom
	treatments=2;
	num_samples=length(groupA)+length(groupB);

	# Compute cohen's f^2
	R_sqrd=ssB/ssT;
	cohens_fsqrd=R_sqrd/(1-R_sqrd);

	# Return all results
	result=list();
	result$SST=ssT;
	result$SSB=ssB;
	result$SSW=ssW;
	result$dfB=treatments-1;
	result$dfW=num_samples-treatments;
	result$F=(ssB/result$dfB)/(ssW/result$dfW);
	result$cohens_fsqrd=cohens_fsqrd;
	result$eta_sqrd=R_sqrd;

	return(result);

}

#-------------------------------------------------------------------------------

combined_hist=function(dataA, dataB, mainA="A", mainB="B", alpha, beta, crit_val, color, range=NULL){
	dataAB=c(dataA, dataB);
	num_breaks=ceiling(log2(length(dataAB))+2)*2
	#cat("Num Breaks: ", num_breaks, "\n");
	hist_comb=hist(c(dataA, dataB), breaks=num_breaks, plot=FALSE);
	breaks=hist_comb$breaks;
	histA=hist(dataA, breaks=breaks, plot=FALSE);
	histB=hist(dataB, breaks=breaks, plot=FALSE);

	xrange=range(c(histA$breaks, histB$breaks));
	ymax=max(c(histA$counts, histB$counts));
	#par(mfrow=c(2,1));
	#cat("Bins: \n");
	#print(breaks);
	#cat("Max Counts: ", ymax, "\n");

	if(!is.null(range)){
		xrange=range;
	}
	
	hist(dataA, breaks=breaks, xlim=xrange, ylim=c(0, ymax*1.2), main=mainA, xlab="", col=color);
	abline(v=crit_val, col="blue");
	text(crit_val, ymax, pos=4, sprintf("alpha=%3.3f", alpha), cex=.8);

	hist(dataB, breaks=breaks, xlim=xrange, ylim=c(0, ymax*1.2), main=mainB, xlab="", col=color);
	abline(v=crit_val, col="blue");
	text(crit_val, ymax, pos=2, sprintf("beta=%3.3f", beta), cex=.8);
	text(crit_val, ymax, pos=4, sprintf("power=%3.3f", 1-beta), cex=.8);
	
}

#-------------------------------------------------------------------------------

find_critical_value=function(data, alpha=.05){
	sorted_data=sort(data);
        n=length(data);
        ba=1-(n-1)/n;
        #cat("Best alpha = ", ba, "\n");
        if(ba <= (alpha+.0000000000000001)){
                crit=sorted_data[ceiling(n*(1-(alpha)))];
                return(crit);
        }else{
                return(NA);
        }
}

#-------------------------------------------------------------------------------

find_beta=function(data, critical_value){
	sorted_data=sort(data);
	num_samples=length(data);
	num_beta=sum(sorted_data<=critical_value);
	beta=num_beta/num_samples;
	return(beta);
}

#-------------------------------------------------------------------------------

compute_roc=function(null_dist, alternate_dist){

	null_len=length(null_dist);
	alt_len=length(alternate_dist);
	
	range=range(c(null_dist, alternate_dist));
	steps=seq(range[1], range[2], length.out=200);

	roc=matrix(0, nrow=length(steps), ncol=2);
	idx=1;
	for(crit_cutoff in steps){
		roc[idx, 1]=sum(null_dist>crit_cutoff);
		roc[idx, 2]=sum(alternate_dist>crit_cutoff);
		idx=idx+1;			
	}

	roc[,1]=roc[,1]/null_len; # False Positive (Alpha), X-axis
	roc[,2]=roc[,2]/alt_len;  # True Positive (1-Beta), Y-Axis
	
	return(roc);
}

#-------------------------------------------------------------------------------

add_to_SSB=function(dist, groupA, groupB, value){
	dist_mat=as.matrix(dist);
	mat_dim=nrow(dist_mat);
	new_dist_mat=as.matrix(dist)
	for(ga in groupA){
		for(gb in groupB){
			# Which one is saved in as.dist?
			new_dist_mat[ga, gb] = sqrt(dist_mat[ga, gb]^2+value);
			new_dist_mat[gb, ga] = sqrt(dist_mat[gb, ga]^2+value);
		}
	}
	return(as.dist(new_dist_mat));
}

adjust_etasqrd=function(dist, groupA, groupB, targeted_etasqrd){

	# Interpretation of eta^2 effect sizes
	# eta^2: .1  small
	# eta^2: .25 medium
	# eta^2: .4  large

	# Compute eta^2 based on original distance matrix
	anova_result=matrix_anova(dist, groupA, groupB);

	cat("Targeted Eta^2: ", targeted_etasqrd, "\n");
	cat("  Original Eta^2: ", anova_result$eta_sqrd, "\n");

	# Don't perform adjustment to SSB if there is already an effect greater than targeted
	if(targeted_etasqrd > anova_result$eta_sqrd){

		# Define function to apply adjustment, and test for it's effect
		eta_sqrd_diff=function(adj){
			adj_dist=add_to_SSB(dist, groupA, groupB, adj);
			new_result=matrix_anova(adj_dist, groupA, groupB);
			diff=abs(targeted_etasqrd - new_result$eta_sqrd);
			return(diff);
		}

		# Run 1-dimensional optimization.  
		optim_res=optim(.5, eta_sqrd_diff, method="Brent", lower=0, upper=1000000);

		cat("  Optimized Adjustment to SSB: ", optim_res$par, "\n");

		# Apply final result to distance
		final_adj_dist=add_to_SSB(dist, groupA, groupB, optim_res$par);

		# Compute final eta^2.  This is only used to sanity check.
		new_result=matrix_anova(final_adj_dist, groupA, groupB);
		cat("  Adjusted Eta^2: ", new_result$eta_sqrd, "\n");
		dist=final_adj_dist;

	}else{
		cat("--------------------------------------------\n");
		cat("Actual Eta^2 is already greater than target.\nReturning unadjusted distance matrix.\n");
		cat("--------------------------------------------\n");
	}
	return(dist);

}

################################################################################

library(MASS);
plot_mds=function(dist, groupA, groupB, title){

	nonzero_dist=dist>0;
	mindist=min(dist[nonzero_dist]);

	for(i in 1:length(dist)){
		if(dist[i]==0){
			#dist[i]=1e-323;
			dist[i]=abs(rnorm(1, 0, mindist));
		}
	}

	#print(dist);

	iso=isoMDS(dist);
	num_samples=nrow(iso$points);
	nA=length(groupA);
	nB=length(groupB);

	colors=rep("green", num_samples);
	shapes=rep(2, num_samples);

	point_names=rownames(iso$points);

	for(i in 1:num_samples){
		if(any(point_names[i]==groupA)){
			colors[i]="blue";
			shapes[i]=6;
		}
	} 

	plot(iso$points[,1], iso$points[,2], col=colors, pch=shapes, xlab="Dimension 1", ylab="Dimension 2", main=title);
	
	x_legend=min(iso$points[,1]);
	y_legend=max(iso$points[,2]);
	legend(x_legend, y_legend, legend=c("Group A", "Group B"), col=c("blue", "green"), pch=c(6,2));
}

################################################################################

pdf(paste(OutputRoot, ".pdf", sep=""), height=8.5, width=11);

# Load summary table and get basic info
cat("Loading Summary Table: ", InputFileNameA, "\n");
stA=load_summary_table(InputFileNameA);
rownames(stA)=paste("A:", rownames(stA), sep="");

cat("Loading Summary Table: ", InputFileNameB, "\n");
stB=load_summary_table(InputFileNameB);
rownames(stB)=paste("B:", rownames(stB), sep="");

# Join summary table
cat("Merging counts.\n");
combined_st=join_summary_table(stA, stB);
num_categories=ncol(combined_st);
cat("  Num combined unique categories: ", num_categories, "\n");

# Normalize counts
cat("Normalizing counts.\n");
nst=normalize(combined_st);
#print(nst);

# Compute Distance
cat("Compute pair-wise distance matrix.\n");
dist=vegdist(nst, DistanceType);
min_dist=min(dist);
max_dist=max(dist);
cat("  Min dist: ", min_dist, "\n");
cat("  Max dist: ", max_dist, "\n");

# Compute pseudo-F-stat
groupA=rownames(stA);
groupB=rownames(stB);
#anova_res=matrix_anova(dist, groupA, groupB);

# Adjust effect size if asked for
if(MinimumEffectSize>0){
	par(mfrow=c(2,1));
	plot_mds(dist, groupA, groupB, title="Original Effect");
	dist=adjust_etasqrd(dist, groupA, groupB, TargetEffectSize);
	plot_mds(dist, groupA, groupB, title="Adjusted Effect");
}

# Set up layout for graphs
par(oma=c(.1, .1, 3, .1));
layout_matrix=matrix(c(
	1,3,5,
	2,4,5
	), nrow=2, byrow=T);
layout(layout_matrix);


num_sample_sizes=length(sample_sizes);
colors=rainbow(num_sample_sizes);

result_fields=c("nA", "nB", "alpha", "beta", "power", "critical_value", "gamma", "cohens_fsqrd", "eta_sqrd");
results_matrix=matrix(0, nrow=num_sample_sizes, ncol=length(result_fields));
colnames(results_matrix)=result_fields;

groupsAB=c(groupA, groupB);
roc_list=list();
for(i in 1:num_sample_sizes){

	n=sample_sizes[i];
	cat("\nComputing on sample size: ", n, "\n");

	null_amova_res=list();
	alternate_amova_res=list();

	null_dist=numeric();
	alternate_dist=numeric();
	
	cohens_fsqrd_null=numeric();
	cohens_fsqrd_alt=numeric();

	bs_sampleA_history=list(NumBootstraps);
	bs_sampleB_history=list(NumBootstraps);

	for(b in 1:NumBootstraps){

		# Compute Null F-distribution
		bs_groupA=sample(groupsAB, size=n, replace=TRUE);
		bs_groupB=sample(groupsAB, size=n, replace=TRUE);

		null_amova_res[[b]]=matrix_anova(dist, bs_groupA, bs_groupB);
		null_dist[b]=null_amova_res[[b]]$F;
		cohens_fsqrd_null[b]=null_amova_res[[b]]$cohens_fsqrd;

		# Compute non-central F-distribution
		bs_groupA=sample(groupA, size=n, replace=TRUE);
		bs_groupB=sample(groupB, size=n, replace=TRUE);
		bs_sampleA_history[[b]]=bs_groupA;
		bs_sampleB_history[[b]]=bs_groupB;
		alternate_amova_res[[b]]=matrix_anova(dist, bs_groupA, bs_groupB);
		alternate_dist[b]=alternate_amova_res[[b]]$F;
		cohens_fsqrd_alt[b]=alternate_amova_res[[b]]$cohens_fsqrd;

	}

	# Compute critical cutoff and beta, based on fixed alpha
	crit_val=find_critical_value(null_dist, Alpha);
	cohens_fsqrd_at_alpha=find_critical_value(cohens_fsqrd_null, Alpha);
	beta=find_beta(alternate_dist, crit_val);
	power=1-beta;	

	# Compute detectable effect size, based on fixed gamma
	crit_at_gamma=find_critical_value(alternate_dist, Gamma);
	cohens_fsqrd_at_gamma=find_critical_value(cohens_fsqrd_alt, Gamma);
	eta_sqrd_at_gamma=cohens_fsqrd_at_gamma/(1+cohens_fsqrd_at_gamma);

	# Plot combined histogram based on F-statistic
	combined_hist(null_dist, alternate_dist, 
		paste("Null Distribution\nn=",n,sep=""),
		paste("Alternate Distribution\nn=",n,sep=""),
		Alpha, beta, crit_val,
		col=colors[i]
	);

	# Plot combined histogram based on Cohen's f^2
	combined_hist(cohens_fsqrd_null, cohens_fsqrd_alt, 
		paste("Null (Cohen's f^2)\nn=",n,sep=""),
		paste("Alternate Distribution (Cohen's f^2)\nn=",n,sep=""),
		Alpha, beta, cohens_fsqrd_at_alpha,
		col=colors[i],
		c(0,1)
	);

	# Label the graph
	if((i%%ncol(layout_matrix))==1){
		mtext(paste(FilenameRootA, " vs. ", FilenameRootB, sep=""), outer=TRUE);
	}

	# Plot the MDS of the median sample
	median_alternate_Fstat=median(alternate_dist);
	med_alt_idx=min(which(median_alternate_Fstat>=median_alternate_Fstat)); # If even number of samples, median will be average of the two values closest to center
	bs_groupA=bs_sampleA_history[[med_alt_idx]];
	bs_groupB=bs_sampleB_history[[med_alt_idx]];
	bs_combined=c(bs_groupA, bs_groupB);
	median_dist=as.dist(as.matrix(dist)[bs_combined, bs_combined]);
	plot_mds(median_dist, bs_groupA, bs_groupB, "Median Alternate Configuration");

	# Compute ROC Curve
	cat("Computing for n =", n, "\n");
	n_str=sprintf("%02i", n);
	roc_list[[n_str]]=compute_roc(null_dist, alternate_dist);

	cat("  Alpha = ", Alpha, "\n");
	cat("  Critical Value = ", crit_val, "\n");
	cat("  Beta = ", beta, "\n");
	cat("  Power = ", power, "\n");
	cat("\n");
	cat("  Gamma = ", Gamma, "\n");
	cat("  Cohen's f^2 = ", cohens_fsqrd_at_gamma, "\n");
	cat("  Eta^2 = ", eta_sqrd_at_gamma, "\n");

	# Store key results in matrix
	results_matrix[i,]=c(n, n, Alpha, beta, power, crit_val, Gamma, cohens_fsqrd_at_gamma, eta_sqrd_at_gamma);

}

################################################################################
# Set up a black ROC plot 
par(mfrow=c(1,1));
plot(
	c(0,1), c(0,1),
 	type="l", lty=2, lwd=2, col="grey",
	xlim=c(0,1), ylim=c(0,1),
	main="Receiver Operating Characteristic",
	xlab="False Positive Rate", ylab="True Positive Rate"
);

title(xlab=expression(alpha), line=2, cex=.7);
title(ylab=expression((1-beta)==pi), line=2, cex=.7);

# Cycle through all the sample sizes and draw the ROC curves 
idx=1;
sizes=names(roc_list);
for(n in sizes){
	x=roc_list[[n]][,1];
	y=roc_list[[n]][,2];
	points(x, y, type="l", col=colors[idx]);
	idx=idx+1;
}
legend(.75,.7, legend=sizes, fill=colors, title="Samples Sizes: ", cex=.7)

################################################################################
# Output core statistics into file

mat_col_names=colnames(results_matrix);

fh=file(paste(OutputRoot, ".power_calc.txt", sep=""), "w");
cat(file=fh, paste("#SampleA,SampleB,distance_type,", paste(mat_col_names, collapse=",") ), "\n");
for(i in 1:num_sample_sizes){
	cat(file=fh, paste(c(FilenameRootA, FilenameRootB, DistanceType, results_matrix[i,]), collapse=","), "\n");
}
close(fh);

################################################################################

cat("Done.\n")

q(status=0)
