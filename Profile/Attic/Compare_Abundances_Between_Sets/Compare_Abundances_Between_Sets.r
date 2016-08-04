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
	"sets_file", "s", 1, "character",
	"output_fname_root", "o", 2, "character",
	"test_run", "t", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table file>\n",
	"	-s <sets input tsv file>\n",
	"	[-o <output filename root>]\n",
	"	[-t (test run flag)]\n",
	"\n",
	"Reads in the summary table and sets file and computes\n",
	"log(odds2/odds1), i.e. log odd ratios, between the pair of sets\n",
	"\n",
	"Bootstrapping is used to compute the confidence intervals\n",
	"between the two sets, by resampling with replacement both\n",
	"the samples and the counts in each category.\n",
	"\n",
	"The first column of the sets input file should be the sample\n",
	"id.  The subsequence columns should consist of 1's and 2's.\n",
	"0's or blanks are considered not in the set.  1's are considered\n",
	"the reference (denominator of the ratio) and 2's will be in the\n",
	"subject (numerator of the ratio).\n",
	"\n",
	"The outputs are:\n",
	"	1.) PDF with each category sorted from negative LogOddsRatio\n",
	"		to positive.  The median is marked with a black +,\n",
	"		the 95% CI are parked with green |'s.  The p-value\n",
	"		for testing if the LogOddsRatio is non-zero is provided.\n",
	"		The FDR column contains the FDR adjusted p-values.\n",
	"	2.) A table which contains the same data as the PDF, but\n",
	"		which can be read into a spreadsheet and sorted by the user\n",
	"		arbitrarily.\n",
	"\n",
	"For each valid subsets column in the sets input file the described\n",
	"outputs will be generated.\n",
	"\n");

if(!length(opt$input_file) || !length(opt$sets_file)){
	cat(usage);
	q(status=-1);
}

if(length(opt$correlation_cutoff)){
	CorrelCutoff=opt$correlation_cutoff;
}

if(!length(opt$output_fname_root)){
    	OutputFnameRoot=gsub(".summary_table.xls", "", opt$input_file);
    	OutputFnameRoot=gsub(".summary_table.tsv", "", OutputFnameRoot);
}else{
    	OutputFnameRoot=opt$output_fname_root;
}

TestRun=F;
if(length(opt$test_run)){
	TestRun=T;
}

###############################################################################

InputFileName=opt$input_file;
SetsFileName=opt$sets_file;

cat("\n");
cat("Input File Name: ", InputFileName, "\n");
cat("Sets File Name: ", SetsFileName, "\n");
cat("Output File Name Root: ", OutputFnameRoot, "\n");
cat("\n");
cat("\n");

if(TestRun){
	cat("WARNING: Test Run.\n");
}

###############################################################################

# Load data
inmat=as.matrix(read.delim(InputFileName, sep="\t", header=TRUE, check.names=FALSE, 
	comment.char="", quote="", row.names=1))
#cat("Original Matrix:\n")
#print(inmat);

# Grab columns we need into a vector, ignore totals, we won't trust it.
counts_mat=inmat[,2:(ncol(inmat)), drop=F];
#print(counts_mat);

num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);

category_names=colnames(counts_mat);
sample_names=rownames(counts_mat);

# Report data summary
num_categories=ncol(counts_mat);
num_samples=nrow(counts_mat);
cat("Number of Categories: ", num_categories, "\n");
cat("Number of Samples: ", num_samples, "\n");
cat("\n");

if(length(unique(category_names))!=num_categories){
	cat("Error:  Category names are not unique.\n");
	quit(status=-1);
}else{
	cat("Category names are indeed unique.\n");
}

###############################################################################
# Load Sets

sets_table=read.table(SetsFileName, sep="\t", header=T, check.names=F,
	comment.char="", quote="", row.names=1);
orig_tabl_ncol=ncol(sets_table);
orig_tabl_colnames=colnames(sets_table);
orig_tabl_rownames=rownames(sets_table);

clean_sets_table=numeric();
clean_sets_colname=character();

for(i in 1:orig_tabl_ncol){
	valid_col=T;
	set_val=unique(sets_table[,i]);

	extra_char=setdiff(set_val, c(0, 1, 2, NA))
	if(length(extra_char)){
		cat("  ", orig_tabl_colnames[i], "has invalid values:", extra_char, "\n");
		valid_col=F;
	}

	missing_char=setdiff(c(1,2), set_val);
	if(length(missing_char)){
		cat("  ", orig_tabl_colnames[i], "is missing values:", missing_char, "\n");
		valid_col=F;
	}

	if(valid_col){
		clean_sets_table=cbind(clean_sets_table, sets_table[,i]);
		clean_sets_colname=c(clean_sets_colname, orig_tabl_colnames[i]);
	}
}

colnames(clean_sets_table)=clean_sets_colname;
rownames(clean_sets_table)=orig_tabl_rownames;
num_sets=length(clean_sets_colname);

cat("\n");
#print(clean_sets_table);

###############################################################################

sets_table_samples=orig_tabl_rownames;
counts_table_samples=sample_names;

shared_samples=sort(intersect(sets_table_samples, counts_table_samples));
cat("Samples shared between summary table and sets file:\n");
print(shared_samples);
cat("\n");

shared_clean_sets_table=clean_sets_table[shared_samples,];
shared_counts_mat=counts_mat[shared_samples,];
shared_num_samples=length(shared_samples);
shared_sample_sizes=apply(shared_counts_mat, 1, sum);

#print(shared_clean_sets_table)
#print(shared_counts_mat);

###############################################################################
# Normalize

# Sum sample totals
sample_totals=apply(shared_counts_mat, 1, sum);

# normalize, to compute probabilities
normalized=matrix(0, nrow=shared_num_samples, ncol=num_categories);
colnames(normalized)=category_names;
rownames(normalized)=rownames(shared_counts_mat);
for(i in 1:shared_num_samples){
	normalized[i,]=shared_counts_mat[i,]/sample_totals[i];
}
#print(normalized);

# Set minimum abundance to be 1/10 of min
min=min(normalized[normalized!=0]);
cat("Min non-zero abundance:", min, "\n");

zerosubstitute=min/10;
cat("Setting zeros to: ", zerosubstitute, "\n\n");

# Generate abundances assuming 0 are 1/10 minimum
norm_no_zeros=matrix(NA, ncol=num_categories, nrow=shared_num_samples);
for(j in 1:num_categories){
	for(i in 1:shared_num_samples){
		if(normalized[i,j]==0){
			norm_no_zeros[i,j]=zerosubstitute;
		}else{
			norm_no_zeros[i,j]=normalized[i,j];
		}
	}
}

# Compute average abundance across all samples
avg_abundance=apply(normalized, 2, mean);
names(avg_abundance)=colnames(normalized);

#print(avg_abundance);
#print(normalized);

#print(apply(normalized, 1, sum));
#print(apply(normalized, 2, mean));

###############################################################################
# Bootstrap the log odds

resample=function(normalized, counts){
	num_samp=nrow(normalized);
	seq=sample(1:num_samp, replace=T);
	new_counts=numeric();
	for(i in seq){
		new_counts=cbind(new_counts, rmultinom(1, counts[i], prob=normalized[i,]));
	}
	colnames(new_counts)=sprintf("%s.%i", rownames(normalized)[seq], 1:num_samp);
	return(t(new_counts));
}

zero_sub=function(x){
	is.zero=(x==0);
	x[is.zero]=zerosubstitute;
	return(x);
}

compute_mean_odds=function(counts_table){
	norm=t(apply(counts_table, 1, function(x){x/sum(x)}));
	norm=zero_sub(norm);
	odds=norm/(1-norm);
	mean=apply(odds, 2, mean);
	return(mean);
}

plot_log_odd_ratios=function(lor_table, pval, cor_pval, title=""){
	num_cat=ncol(lor_table);
	num_samp=nrow(lor_table);	
	cat_names=colnames(lor_table);
	#cat_names=1:num_cat;

	plots_per_page=50;

	lor_val=as.vector(lor_table);
	finite_lor=lor_val[is.finite(lor_val)];
	range=range(finite_lor);
	magnitude=max(abs(range));
	brks=seq(-magnitude*1.1, magnitude*1.1, length.out=100);


	new_page=function(){
		par(mar=c(5,20,1,3));
		plot(0,0, type="n",
			ylim=c(1, plots_per_page),
			xlim=c(-magnitude, magnitude+3),
			xlab="Log10(odds(Comparator)/odds(Reference))",
			ylab="",
			yaxt="n",
			xaxt="n",
			bty="n",
			main=title
		);
		axis(side=1, at=seq(-magnitude, magnitude, length.out=5), 
			labels=sprintf("%3.1f", seq(-magnitude, magnitude, length.out=5)));
		abline(v=0, col="red", lty=2);
		abline(h=0:plots_per_page, col="gray75");
		text(magnitude+.5, plots_per_page+.5, labels="p-value        FDR", cex=.7, font=2);
	}

	for(i in 1:num_cat){
		cat_lor=as.vector(lor_table[,i]);
		finite_cat_lor=cat_lor[is.finite(cat_lor)];
		#hist(finite_cat_lor, xlim=c(-magnitude, magnitude), breaks=brks, freq=F);
		ci95=quantile(finite_cat_lor, c(.025, .5, .975));
		smoothed=density(finite_cat_lor, bw=.2);

		y_basis=plots_per_page-((i-1)%%plots_per_page) -1;
		y=(.8*smoothed$y/max(smoothed$y))+ y_basis;

		if(((i-1)%%plots_per_page)==0){
			new_page();
			axis(side=2, at=(plots_per_page:1)-.5, 
				labels=cat_names[(1:plots_per_page)+((i-1)%/%plots_per_page)*plots_per_page],
				las=2);
		}

		points(smoothed$x, y, type="l", col="blue");
		#text(-5, y_basis, labels=i);

		points(ci95, c(y_basis, y_basis , y_basis), 
			pch=c("|", "+", "|"), 
			cex=c(1,.7,1),
			col=c("green", "black", "green")
		);
		text(magnitude+.5, y_basis+.5, labels=sprintf("%8.4f   %8.4f", pval[i], cor_pval[i]), cex=.7);
	}


}

if(TestRun){
	rnd=sprintf(".%i", ceiling(runif(1, 1,1000)));
	B=40;
}else{
	rnd="";
	B=2000;
}


compute_pval=function(x){
	len=length(x);
	gt_pval=1-(sum(x>0)/len);
	lt_pval=1-(sum(x<0)/len);
	pval=min(1, min(c(lt_pval, gt_pval))*2);
	return(pval);
}

for(set_ix in 1:num_sets){

	cur_setname=clean_sets_colname[set_ix];
	cat("Analyzing: ", cur_setname, "\n");

	
	# Extract out members in the set
	ref_ix=shared_clean_sets_table[,set_ix]==1;
	ref_ix[is.na(ref_ix)]=F;
	comp_ix=shared_clean_sets_table[,set_ix]==2;
	comp_ix[is.na(comp_ix)]=F;

	#print(ref_ix);
	#print(comp_ix);

	# Store bootstrap results
	bs_lor=matrix(NA, ncol=num_categories, nrow=B);
	colnames(bs_lor)=category_names;

	# Run bootstrap
	for(bix in 1:B){

		# Resample
		ref_counts=resample(normalized[ref_ix,], shared_sample_sizes[ref_ix]);	
		comp_counts=resample(normalized[comp_ix,], shared_sample_sizes[comp_ix]);	

		# Compute mean odds within each group
		ref_odds=compute_mean_odds(ref_counts);
		comp_odds=compute_mean_odds(comp_counts);

		# Compute log odds ration of the means
		log_odds_ratio=log(comp_odds/ref_odds);

		# Store the results of this bootstrap iteration
		bs_lor[bix, ]=log_odds_ratio;

	}

	# Sort by mean, then by median
	mean_logodds=apply(bs_lor, 2, mean);
	mean_sortix=order(mean_logodds);
	bs_lor_mean_sorted=bs_lor[, mean_sortix];
	med_logodds=apply(bs_lor_mean_sorted, 2, median);
	median_sortix=order(med_logodds);
	bs_lor_median_sorted=bs_lor_mean_sorted[, median_sortix];

	# Compute pvalues for each category
	pvalues=apply(bs_lor_median_sorted, 2, compute_pval);
	corr_pval=p.adjust(pvalues, method="fdr");

	# Plot odds and p-values
	pdf(paste(OutputFnameRoot, ".", cur_setname, rnd, ".comp_sets.pdf", sep=""), width=8.5, height=14);
	plot_log_odd_ratios(bs_lor_median_sorted, pvalues, corr_pval, title=cur_setname);
	dev.off();

	# Output statistics in a table
	fh=file(paste(OutputFnameRoot, ".", cur_setname, rnd, ".comp_sets.tsv", sep=""), "w");
	qtl=apply(bs_lor_median_sorted, 2, function(x){ quantile(x, c(.5, .025, .975))});
	median_sorted_colnames=colnames(bs_lor_median_sorted);
	headers=c("Category", "median", "2.5%", "97.5%", "uncorrected p-value", "FDR p-value");
	cat(file=fh, paste(headers, collapse="\t"), "\n");
	
	# Output lines
	for(cat_ix in 1:num_categories){
		dataline=paste(c(
			median_sorted_colnames[cat_ix],
			qtl[1,cat_ix],
			qtl[2,cat_ix],
			qtl[3,cat_ix],
			pvalues[cat_ix],
			corr_pval[cat_ix]
			),	
			collapse="\t");
		cat(file=fh, dataline, "\n");
	}

	close(fh);
	cat("\n");

}






###############################################################################

cat("Done.\n")
if(length(warnings())){
	print(warnings());
}
q(status=0)
