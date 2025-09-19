#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

options(useFancyQuotes=F);

params=c(
	"summary_file", "s", 1, "character",
	"prop_extr", "p", 1, "numeric",
	"output_file_root", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

DEF_EXTRACT_PROP=0.95;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table for taxa/function>\n",
	"	-p <proportion of categorical variables to extract, default=", DEF_EXTRACT_PROP, ">\n",
	"	[-o <output filename root>]\n",
	"\n",
	"Compute ALR transformed values and export to tsv file.\n",
	"\n",
	"The -p variable, e.g. .95, tell us to extract 95% of the categories, leaving 5% for the\n",
	"denominator.  The 95% references to the average abundance across all samples.\n",
	"\n", sep="");

if(!length(opt$summary_file)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$output_file_root)){
	OutputRoot=gsub(".summary_table.xls", "", opt$summary_file);
	OutputRoot=gsub(".summary_table.tsv", "", OutputRoot);
}else{
	OutputRoot=opt$output_file_root;
}

if(!length(opt$denom_prop)){
	ExtractProp=DEF_EXTRACT_PROP;
}else{
	ExtractProp=opt$prop_extr;
}

SummaryFile=opt$summary_file;
OutputRoot=paste(OutputRoot, ".alr.", ExtractProp*100, sep="");

cat("\n");
cat("Summary File: ", SummaryFile, "\n", sep="");
cat("Output File: ", OutputRoot, "\n", sep="");
cat("Extraction Proportion: ", ExtractProp, "\n", sep="");
cat("\n");

##############################################################################

load_summary_file=function(fname){
	inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, 
		check.names=FALSE, comment.char="", quote="", row.names=1))

	counts_mat=inmat[,2:(ncol(inmat))];

	# Clean category names a little
	cat_names=colnames(counts_mat);
	cat_names=gsub("-", "_", cat_names);
	colnames(counts_mat)=cat_names;
	
	cat("Num Categories in Summary Table: ", ncol(counts_mat), "\n", sep="");
	return(counts_mat);
}

normalize=function(counts){
	totals=apply(counts, 1, sum);
	num_samples=nrow(counts);
	normalized=matrix(0, nrow=nrow(counts), ncol=ncol(counts));

	for(i in 1:num_samples){
		normalized[i,]=counts[i,]/totals[i];
	}
	
	colnames(normalized)=colnames(counts);
	rownames(normalized)=rownames(counts);	
	return(normalized);
}

plot_text=function(strings){

	orig.par=par(no.readonly=T);
	par(family="Courier");
	par(oma=rep(.1,4));
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

	cat(strings, sep="\n");

	par(orig.par);
}

##############################################################################
##############################################################################

# Load summary file table counts 
cat("Loading summary table...\n");
counts=load_summary_file(SummaryFile);

# Remove zero count samples
cat("Looking for zero-count samples to remove.\n");
tot=apply(counts, 1, sum);
nonzero=tot>0;
if(!(all(nonzero))){
	cat("WARNING: Zero count samples found:\n");
	samp_names=rownames(counts);
	print(samp_names[!nonzero]);
	cat("\n");
	counts=counts[nonzero,,drop=F];
}else{
	cat("\tGood.  None found.\n");
}

num_categories=ncol(counts);
num_samples=nrow(counts);

# Normalize
cat("Normalizing counts to proportions.\n");
pure_counts=counts;
pure_normalized=normalize(pure_counts);
pure_totals=apply(counts, 1, sum);

pure_mean_abund=apply(pure_normalized, 2, mean);
ix=order(pure_mean_abund, decreasing=F);

pure_counts=pure_counts[,ix];
pure_normalized=pure_normalized[,ix];
pure_mean_abund=pure_mean_abund[ix];
pure_cum_sum=cumsum(pure_mean_abund);

sorted_taxa_names=colnames(pure_normalized);

cat("Top 10 categories:\n");
print(pure_mean_abund[1:10]);

##############################################################################

pdf_fname=paste(OutputRoot, ".alr_summary.pdf", sep="");
pdf(pdf_fname, height=11, width=9.5);

par(mfrow=c(1,1));
plot_text(c(
	paste("Summary File: ", SummaryFile, "\n", sep=""),
	paste("Output File: ", OutputRoot, "\n", sep=""),
	paste("Extraction Proportion: ", ExtractProp, "\n", sep="")
));

##############################################################################
#------------------------------------------------------------------------------
# Sampling Depth Histograms

par(mfrow=c(2,1));
hist(pure_totals, main="Sample Depths", 
	xlab="Read Counts", ylab="Frequency of Samples", breaks=100);
hist(log10(pure_totals), main="Sample Depths (Log Scale)", 
	xlab="Log10(Read Counts)", ylab="Frequency of Samples", breaks=100);

##############################################################################
# Compute categories counts necessary to achieve landmarks

samp_prop_landmarks=c(.05, .1, .15, .2, .5, 1-ExtractProp);
num_samp_prop_landmarks=length(samp_prop_landmarks);
samp_prop_landmark_category_counts=numeric(num_samp_prop_landmarks);
names(samp_prop_landmark_category_counts)=samp_prop_landmarks;

landmark_widths=rep(1, num_samp_prop_landmarks);
landmark_widths[num_samp_prop_landmarks]=2;
landmark_colors=rep("blue", num_samp_prop_landmarks);
landmark_colors[num_samp_prop_landmarks]="red";

for(i in 1:num_samp_prop_landmarks){
	samp_prop_landmark_category_counts[i]=sum(pure_cum_sum < samp_prop_landmarks[i]);
}

cat("Number of Categories to achieve landmark:\n");
print(samp_prop_landmark_category_counts);

##############################################################################
# Plot cumulative abundances and proportion of samples with nonzero counts as abundance increases

par(mfrow=c(3,1));

plot(1:num_categories, pure_cum_sum, main="Cumulative Abundance",
	xlab="Number of Categories Included (Increasing Abundance)", 
	ylab="Proportion of Total Abundance",
	cex=.5
	);
abline(h=samp_prop_landmarks, col=landmark_colors, lty="dashed", lwd=landmark_widths);
abline(v=samp_prop_landmark_category_counts, col=landmark_colors, lty="dotted", lwd=landmark_widths);


zeros=apply(pure_counts, 2, function(x){ sum(x>0)});
prop_zeros=zeros/num_samples;

plot(1:num_categories, prop_zeros, main="Proportion of Samples with Non-zero Category Counts",
	xlab="Category Identifier (Increasing Abundance)", ylab="Proportion of Samples",
	cex=.5
	);
abline(h=samp_prop_landmarks, col=landmark_colors, lty="dashed", lwd=landmark_widths);
abline(v=samp_prop_landmark_category_counts, col=landmark_colors, lty="dotted", lwd=landmark_widths);

#..............................................................................

calc_resample_prob=function(norm, samp_tot){

	num_samples=nrow(norm);
	num_categories=ncol(norm);
	
	prob_nonzero_mat=matrix(NA, nrow=num_samples, ncol=num_categories);
	colnames(prob_nonzero_mat)=colnames(norm);
	rownames(prob_nonzero_mat)=rownames(norm);

	print(samp_tot)
	for(i in 1:num_samples){
		for(j in 1:num_categories){
			prob_nonzero_mat[i,j]=1-pbinom(0, samp_tot[i], norm[i,j]);
		}
	}

	return(prob_nonzero_mat);
}

#..............................................................................

p_nz_mat=calc_resample_prob(pure_normalized, pure_totals);

nonz_resamp=apply(p_nz_mat, 2, function(x){sum(x>.95)});
prop_nonz_resamp=nonz_resamp/num_samples;

plot(1:num_categories, prop_nonz_resamp, ylim=c(0,1), 
	main="Proportion of Samples with Non-zero Category Counts (95% of time, if resampled)",
	xlab="Category Identifier (Increasing Abundance)", ylab="Proportion of Samples",
	cex=.5
	);
abline(h=samp_prop_landmarks, col=landmark_colors, lty="dashed", lwd=landmark_widths);
abline(v=samp_prop_landmark_category_counts, col=landmark_colors, lty="dotted", lwd=landmark_widths);

##############################################################################
# Non-zero resampling vs. observed

par(mfrow=c(1,1));
par(mar=c(4,4,4,1));

plot(0, type="n", xlim=c(0,1), ylim=c(0,1),
	main="Proportion of Samples Remaining Non-zero if Resampled",
	xlab="Observed Proportion NonZero", ylab="Predicted Proportion NonZero (95% of time, if resampled)");
abline(a=0,b=1, col="blue");
points(prop_zeros, prop_nonz_resamp, cex=.5);

##############################################################################
# Mean Abundance vs. Prop Non zero

plot(log10(pure_mean_abund), prop_zeros, 
	main="Observed Proportion of Samples Non-Zero vs. Mean Abundance (Log Scale)",
	xlab="Log10(Mean Abundance)", ylab="Proportion of Samples NonZero",
	cex=.5
	);

##############################################################################
# Accumulation of denomator plots

cumul_denominator=numeric(num_samples);
mean_denom=numeric(num_categories);
sd_denom=numeric(num_categories);
cov_denom=numeric(num_categories);

for(cat_ix in 1:num_categories){
	cumul_denominator=cumul_denominator+pure_normalized[,cat_ix];
	mean_denom[cat_ix]=mean(cumul_denominator);
	sd_denom[cat_ix]=sd(cumul_denominator);
	cov_denom[cat_ix]=sd_denom[cat_ix]/mean_denom[cat_ix];
}

par(mfrow=c(4,1));
par(oma=c(3,0,3,1));

# Accumulation means
plot(1:num_categories, mean_denom, ylab="Mean", xlab="", main="Mean");
abline(h=samp_prop_landmarks, col=landmark_colors, lty="dashed", lwd=landmark_widths);
abline(v=samp_prop_landmark_category_counts, col=landmark_colors, lty="dotted", lwd=landmark_widths);

# Accumulation Std Dev.
plot(1:num_categories, sd_denom, ylab="St. Dev.", xlab="", main="Standard Deviation");
abline(v=samp_prop_landmark_category_counts, col=landmark_colors, lty="dotted", lwd=landmark_widths);

# Accumulation of coeff of variance
plot(1:num_categories, cov_denom, ylab="StDev/Mean", xlab="", 
	main="Coefficient of Variance (lower is better)");
abline(v=samp_prop_landmark_category_counts, col=landmark_colors, lty="dotted", lwd=landmark_widths);

mtext("Accumulation of 'Denominator' (Increasing Abundance)\n", 
	side=3, line=-1, outer=T, font=2, cex=1.3);
mtext("Number of Categories Included", 
	side=1, line=1, outer=T, cex=1);

# Non-cumulative Mean abundances
plot(1:num_categories, log10(pure_mean_abund),
	main="Mean Abundances (Not Cumulative)",
	xlab="Category Identifer (Increasing Abundance)",
	ylab="Log10(Abundance)", cex=.5);
abline(v=samp_prop_landmark_category_counts, col=landmark_colors, lty="dotted", lwd=landmark_widths);

##############################################################################
# Calculate ALR 
par(mfrow=c(2,1));

denom_ix=pure_cum_sum<=(1-ExtractProp);
denominator_mat=pure_normalized[,denom_ix];
numerators_mat=pure_normalized[,!denom_ix];

num_var_in_denom=sum(denom_ix);
num_vars_for_alr=sum(!denom_ix);

summed_denom=apply(denominator_mat, 1, sum);
acquired_denom_prop=mean(summed_denom);

# 
hist(summed_denom, breaks=100, main="Acquired Denominators (Normalized)",
	xlab="Denominators", ylab="Frequency of Samples"
);


# Calculate the ALR
# Add .5 to the counts of the numerator and denominator

alr_matrix=matrix(NA, nrow=num_samples, ncol=num_vars_for_alr);
rownames(alr_matrix)=rownames(numerators_mat);
colnames(alr_matrix)=colnames(numerators_mat);

samples_with_zeros=numeric(num_samples);

denominators_counts_mat=pure_counts[,denom_ix];
numerators_counts_mat=pure_counts[,!denom_ix];
summed_denom_counts=apply(denominators_counts_mat, 1, sum);

hist(log10(summed_denom_counts+.5), breaks=100, main="Acquired Denominators (Counts)",
	xlab="Log10(Denominators+.5)", ylab="Frequency of Samples"
);

for(i in 1:num_samples){
	samp_counts=numerators_counts_mat[i,];
	samples_with_zeros[i]=sum(samp_counts==0);
	numerator_plus=samp_counts+.5;
	denom_plus=summed_denom_counts[i]+.5;
	logratio=log(numerator_plus/denom_plus);
	alr_matrix[i,]=logratio;
}

hist(log10(samples_with_zeros), breaks=100, 
	main="Distribution of Zero-Abundances Categories (Numerators) Before ALR Transformation",
	xlab="Log10(Number of Zero Count Categories in each Sample)",
	ylab="Frequency of Samples"
);

##############################################################################
# Plot some histograms for the top categories

var_names=colnames(alr_matrix);
par(mfrow=c(3,3));

plot_histogram=function(ix, values, vname, col){

	mean_alr=mean(values);
	sd_alr=sd(values);
	hist(values, breaks=25, 
		main=paste(ix, ".) ", vname, sep=""),
		xlab="ALR", ylab="Frequency",
		col=col
	);
	mtext(paste(
		"Mean: ", round(mean_alr, 4), "\n",
		"St Dev: ", round(sd_alr, 4), "\n",
		"Coef Var: ", round(sd_alr/mean_alr, 4), 
		sep=""), cex=.75, line=-2, side=3);

};

# Show top and bottom 18 alr histograms
num_alr_val=ncol(alr_matrix);

par(mfrow=c(3,3));
 
# Bottom
for(i in 1:(9*2)){
	val=alr_matrix[,i];
	names=var_names[i];
	plot_histogram(i, val, names, col="skyblue");
}

# Top
for(i in num_alr_val-(1:(9*2))){
	val=alr_matrix[,i];
	names=var_names[i];
	plot_histogram(i, val, names, col="pink");
}


##############################################################################
par(mfrow=c(1,1));

plot_text(c(
	"Extraction Information:",
	paste("Target Proportion: ", ExtractProp, sep=""), 
	paste("Denominator Proportion: ", 1-ExtractProp, sep=""),
	paste("Acquired Denominator Proportion: ", acquired_denom_prop, sep=""),
	"",
	paste("Number of Variables in Denominator: ", num_var_in_denom, sep=""),
	paste("Number of Variables available for ALR: ", num_vars_for_alr, sep=""),	
	""
));

##############################################################################
##############################################################################

# Write ALR to tsv file 
tsv_fname=paste(OutputRoot, ".tsv", sep="");
cat("Writing TSV of ALR values to file: ", tsv_fname, "\n", sep="");
fh=file(tsv_fname, "w");
cat(file=fh, "SampleID");
close(fh);
write.table(alr_matrix, file=tsv_fname, quote=F, sep="\t", row.names=T, col.names=NA, append=T);

##############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
