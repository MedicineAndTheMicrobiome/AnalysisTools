#!/usr/bin/env Rscript

###############################################################################
#                                                                             #
#       Copyright (c) 2013 J. Craig Venter Institute.                         #
#       All rights reserved.                                                  #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                             #
###############################################################################
###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"num_categories", "p", 2, "numeric",
	"output_fname_root", "o", 2, "character",
	"testing_flag", "t", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

NumTargetedCategories=100;

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	[-p <number top categories to output, default = ", NumTargetedCategories, ">\n",
	"	[-o <output filename root>]\n",
	"	[-t (Testing flag)]\n",
	"\n",	
	"Reads in a summary table, and then based on the categories\n",
	"with the most variation and abundance, will generate a new summary table\n",
	"\n",
	"First, each abundance is transformed into a log(odds), i.e. logit\n",
	"to improve the normalness of the distribution of abundances within\n",
	"each category.  The mean and variance of the logit for each category is then\n",
	"calculated and a scatter plot is generated.\n",
	"\n",
	"Each category is then assigned two quantile scores by sorting all the\n",
	"categories with respect to their abundance and also with respect to\n",
	"their variance.  A final score is computed based on the euclidean distance from\n",
	"the worst score (low abundance and low variance), so that the combined score is:\n",
	"\n",
	"	score = sqrt(abundance_quantile^2 + variance_quantile^2)\n",
	"\n",
	"Then the top p categories (with the best combination of high variance and high\n",
	"abundance) based on the combined score are selected.\n",
	"The remaining categories are collapsed into a single \"remaining\" \n",
	"category, and a new summary table is generated.\n",	
	"\n",
	"The outputs that are generated are:\n",
	"	1.) A new summary table with the selected high variance/abundance categories.\n",
	"	2.) A scatter plot identifying which categories were selected.\n",
	"	3.) A table of log odd means, variances, quantiles, scores, etc...\n",
	"	4.) A list of categories that were kept.\n",
	"\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$output_fname_root)){
    	OutputFnameRoot=gsub(".summary_table.xls", "", opt$input_file);
    	OutputFnameRoot=gsub(".summary_table.tsv", "", OutputFnameRoot);
}else{
    	OutputFnameRoot=opt$output_fname_root;
}

if(length(opt$num_categories)){
	NumTargetedCategories=opt$num_categories;
}

Testing_Flag=F;
if(length(opt$testing_flag)){
	Testing_Flag=T;
}

###############################################################################

InputFileName=opt$input_file;

cat("\n");
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name Root: ", OutputFnameRoot, "\n");
cat("Number of Most Ideal Categories to Keep: ", NumTargetedCategories, "\n");
cat("\n");

###############################################################################
###############################################################################

# Load data
inmat=as.matrix(read.delim(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="", quote="", row.names=1))
#cat("Original Matrix:\n")
#print(inmat);

# Grab columns we need into a vector, ignore totals, we won't trust it.
counts_mat=inmat[,2:(ncol(inmat)), drop=F];
#print(counts_mat);

num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);

category_names=colnames(counts_mat);
sample_names=rownames(counts_mat);

num_categories=ncol(counts_mat);
num_samples=nrow(counts_mat);
cat("Number of Categories: ", num_categories, "\n");
cat("Number of Samples: ", num_samples, "\n");

if(num_categories<=NumTargetedCategories){
	cat("\n\n");
	cat("******************************************************************************\n");
	cat("*                                                                            *\n");
	cat("*  WARNING: You want to keep ", NumTargetedCategories, 
		" categories, but there are only ", num_categories, "!!!\n", sep="");
	cat("*                                                                            *\n");
	cat("******************************************************************************\n");
	cat("\n\n");
	NumTargetedCategories=num_categories;
}

if(Testing_Flag){
	rnd=sprintf(".%03g", ceiling(runif(1, 0, 1000)));
}else{
	rnd="";
}
pdf(paste(OutputFnameRoot, rnd, ".var_clps.pdf", sep=""), height=8.5, width=11);

###############################################################################
# Normalize

# Sum sample totals
sample_totals=apply(counts_mat, 1, sum);

# normalize, to compute probabilities
normalized=matrix(0, nrow=num_samples, ncol=num_categories);
for(i in 1:num_samples){
	normalized[i,]=counts_mat[i,]/sample_totals[i];
}
#print(normalized);

# Set minimum abundance to be 1/10 of min
min=min(normalized[normalized!=0]);
max=max(normalized[normalized!=1]);
cat("Min non-zero abundance:", min, "\n");

zerosubstitute=min/10;
onesubstitute=1-(1-max)/10;

cat("\nTo avoid infinities...\n");
cat("Setting zeros to: ", zerosubstitute, "\n");
cat("Setting ones  to: ", onesubstitute, "\n");
cat("\n");


norm_no_zeros=matrix(NA, ncol=num_categories, nrow=num_samples);
for(j in 1:num_categories){
	for(i in 1:num_samples){
		if(normalized[i,j]==0){
			norm_no_zeros[i,j]=zerosubstitute;
		}else if(normalized[i,j]==1){
			norm_no_zeros[i,j]=onesubstitute;
		}else{
			norm_no_zeros[i,j]=normalized[i,j];
		}
	}
}

###############################################################################
# Compute log(odds) for abundance

logodds=function(x){
	return(log10(x/(1-x)));
}

logodds_mat=logodds(norm_no_zeros);
colnames(logodds_mat)=category_names;

logodds_range=range(logodds_mat);
cat("Range Log10(Odds): (", logodds_range[1], ",  ", logodds_range[2], ") \n", sep="");

avg_logodds=apply(logodds_mat, 2, mean);
names(avg_logodds)=category_names;

###############################################################################

# Compute variances on log(odds)
variances=apply(logodds_mat, 2, var);
names(variances)=category_names;

###############################################################################

# Get quantiles so we can mark them on the graph
quant_intvl=seq(0,1,.1);
var_quantiles=quantile(variances, quant_intvl);
abnd_quantiles=quantile(avg_logodds, quant_intvl);

# Plot histograms on both mean and variance of logit
hist(avg_logodds, main="Distribution of Means", xlab="Mean(Logit)",  breaks=nclass.Sturges(avg_logodds)*5);
hist(variances, main="Distribution of Variances", xlab="Var(Logit)",  breaks=nclass.Sturges(variances)*5);

# Compute quantile of each category based on variance and abundance
var_pt_qntl=(rank(variances)/num_categories);
abn_pt_qntl=(rank(avg_logodds)/num_categories);

# Use the euclidean distance from (0,0) (worst case) as a score 
score_pt=sqrt(var_pt_qntl^2 + abn_pt_qntl^2);
top_score_pt_ix=order(score_pt, decreasing=T)[1:NumTargetedCategories];
selected=top_score_pt_ix;

# Color selected points purple, and unselected points black.
colors=rep("black", num_categories);
colors[selected]="purple";

# Scatter plot relationship between variance and abundance
par(oma=c(2.1, 2.1, 2.1, 2.1));
plot(avg_logodds, variances, main="", cex=.4, col=colors, xlab="Mean(Logit)", ylab="Var(Logit)"); 
mtext(sprintf("Selected: %g", NumTargetedCategories), line=-1, side=3);
abline(h=var_quantiles, col="grey", lty=2);
abline(v=abnd_quantiles, col="grey", lty=2);
axis(side=4, at=var_quantiles, labels=sprintf("%i%%", 100*quant_intvl), cex.axis=.5, las=2);
axis(side=3, at=abnd_quantiles, labels=sprintf("%i%%", 100*quant_intvl), cex.axis=.5, las=2);
title(main="Variance vs. Abundance", line=3);

###############################################################################

selected_counts_mat=counts_mat[, selected, drop=F];

print(selected);

if(length(selected) !=0){
	remaining_counts_mat=counts_mat[, -selected, drop=F];
}else{
	cat("No remaining unselected counrs left.\n");
	remaining_counts_mat=NULL;
}
sum_of_selected=apply(selected_counts_mat, 1, sum);

if(!is.null(remaining_counts_mat)){
	sum_of_remaining=apply(remaining_counts_mat, 1, sum);
	ratios_of_remaining=sum_of_selected/sum_of_remaining;
}else{
	sum_of_remaining=0;
	ratios_of_remaining="No ratio.";
}


cat("\n");
cat("Sum of Selected:\n");
print(sum_of_selected);
cat("\n");
cat("Sum of Remaining:\n");
print(sum_of_remaining);
cat("\n");
cat("Ratios of Remaining:\n");
print(ratios_of_remaining);
cat("\n");

output_summary_table=cbind(selected_counts_mat, sum_of_remaining);
colnames(output_summary_table)=c(colnames(selected_counts_mat), "Remaining");

res=dev.off();

###############################################################################

# Write out summary file table

write_summary_table=function(counts_mat, fname){

        stfname=paste(fname, ".summary_table.tsv", sep="");
        cat("Writing summary file table: ", stfname, "\n", sep="");

        st_fh=file(stfname, "w");

        num_samples=nrow(counts_mat);
        sample_names=rownames(counts_mat);
        category_names=colnames(counts_mat);

        # Output Header
        cat(file=st_fh, paste("sample_id", "total", paste(category_names, collapse="\t"), sep="\t"));
        cat(file=st_fh, "\n");

        # Output rows
        for(samp_idx in 1:num_samples){
                total=sum(counts_mat[samp_idx,]);
                outline=paste(sample_names[samp_idx],total,paste(counts_mat[samp_idx,], collapse="\t"), sep="\t");
                cat(file=st_fh, outline, "\n", sep="");
        }

        close(st_fh);
        cat("Done writing.\n");
}

write_summary_table(output_summary_table, paste(OutputFnameRoot, ".var_clps", sep=""));

###############################################################################

# Output all stats

cat("Outputing statististics...\n");
fh=file(paste(OutputFnameRoot, ".var_clps.stats.tsv", sep=""), "w");

out_hdr=paste(c(
	"Category", 
	"E(LogOdds)",
	"Var(LogOdds)",
	"Quantile(E[LogOdds])",
	"Quantile(Var[LogOdds])",
	"Score",
	"Kept"
), collapse="\t");
cat(file=fh, out_hdr, "\n", sep="");

for(i in 1:num_categories){
	out_data=c(
		category_names[i],
		avg_logodds[i],
		variances[i],
		abn_pt_qntl[i],
		var_pt_qntl[i],
		score_pt[i],
		ifelse(any(i==selected), 1, 0)
	);

	cat(file=fh, paste(out_data, collapse="\t"), "\n", sep="");
}
close(fh);

###############################################################################

# Output list of kept categories

cat("Outputing list of kept categories...\n");
fh=file(paste(OutputFnameRoot, ".var_clps.kept.txt", sep=""), "w");

for(i in selected){
	cat(file=fh, category_names[i], "\n", sep="");
}

close(fh);

###############################################################################

cat("Done.\n")
print(warnings());
q(status=0)
