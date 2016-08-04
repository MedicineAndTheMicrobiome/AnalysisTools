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

NUM_BS=2000;

params=c(
	"input_A_file", "a", 1, "character",
	"input_B_file", "b", 1, "character",
	"output_file", "o", 2, "character",
	"num_bootstraps", "i", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n",
	"\n",
	script_name, "\n",
	"\n",
	"	-a <input summary_table.xls file A>\n",
	"	-b <input summary_table.xls file B>\n",
	"	[-o <output file name>]\n",
	"	[-i <num bootstraps, default=", NUM_BS, ">]\n",
	"\n",
	"This script will compute the AWKS statistic\n",
	"by bootstrapping the AWK null distribution assuming\n",
	"that A and B are not different.\n",
	"\n",
	"This means creating a A' and B', by sampling from\n",
	"the actual members of A and B.  Then computing\n",
	"the AWKS statistics BS times, to form a the null distribution.\n",
	"Then the observed AWKS statistics between A and B is\n",
	"compared against the null distribution to see if it is\n",
	"significant.\n",
	"\n",
	"\n", sep="");

if(!length(opt$input_A_file) || !length(opt$input_B_file)){
	cat(usage);
	q(status=-1);
}

InputFileAName=opt$input_A_file;
InputFileANameClean=gsub("\\.summary_table\\.xls$", "", InputFileAName);
InputFileANameClean=gsub("\\.summary_table\\.tsv$", "", InputFileANameClean);

InputFileBName=opt$input_B_file;
InputFileBNameClean=gsub("\\.summary_table\\.xls$", "", InputFileBName);
InputFileBNameClean=gsub("\\.summary_table\\.tsv$", "", InputFileBNameClean);

RootA=tail(strsplit(InputFileANameClean, "/")[[1]], 1);
RootB=tail(strsplit(InputFileBNameClean, "/")[[1]], 1);

OutputRoot=paste(RootA, "_vs_", RootB, sep="");

if(length(opt$output_file)){
	OutputRoot=opt$output_file;
}

num_bs=NUM_BS;
if(length(opt$num_bootstraps)){
	num_bs=opt$num_bootstraps;
}

cat("Input File A: ", InputFileAName, "\n");
cat("Input File B: ", InputFileBName, "\n");
cat("Output File Root: ", OutputRoot, "\n");
cat("Number of Bootstraps: ", num_bs, "\n");

################################################################################

load_summary_table=function(fname){

        inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", quote="", row.names=1))
        counts_mat=inmat[,2:(ncol(inmat))];
        cat("Num Categories in Summary Table: ", ncol(counts_mat), "\n", sep="");
        return(counts_mat);
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

combine_st=function(stA, stB){
	
	A_cat_names=colnames(stA);
	B_cat_names=colnames(stB);
	combined_cat_names=sort(unique(c(A_cat_names, B_cat_names)));
	
	A_num_samples=nrow(stA);
	B_num_samples=nrow(stB);
	combined_num_categories=length(combined_cat_names);
	combined_num_samples=A_num_samples+B_num_samples;
	A_sample_names=rownames(stA);
	B_sample_names=rownames(stB);

	cat("Combined Matrix: ", combined_num_samples, " x ", combined_num_categories, "\n");
	combined_matrix=matrix(0, nrow=combined_num_samples, ncol=combined_num_categories);
	colnames(combined_matrix)=combined_cat_names;
	rownames(combined_matrix)=c(A_sample_names, B_sample_names);

	combined_idx=1;
	for(i in 1:A_num_samples){
		combined_matrix[combined_idx, A_cat_names]=stA[i, A_cat_names];
		combined_idx=combined_idx+1;
	}
	for(i in 1:B_num_samples){
		combined_matrix[combined_idx, B_cat_names]=stB[i, B_cat_names];
		combined_idx=combined_idx+1;
	}

	return(combined_matrix);

}

#-------------------------------------------------------------------------------

avg_abund=function(norm_st){
	avg=apply(norm_st, 2, mean);
	return(avg);
}

category_ks_stat=function(a_norm_abd, b_norm_abd){

	# Only test abundances that are in the dataset
	# instead of all of them
	unique_abund=sort(unique(c(a_norm_abd, b_norm_abd)));
	num_unique_abund=length(unique_abund);

	num_a=length(a_norm_abd);
	num_b=length(b_norm_abd);

	max_diff=0;
	for(i in 1:num_unique_abund){
		num_a_above=sum(a_norm_abd>unique_abund[i])/num_a;
		num_b_above=sum(b_norm_abd>unique_abund[i])/num_b;
		diff=abs(num_a_above-num_b_above);
		if(diff>max_diff){
			max_diff=diff;
		}
	}
	return(max_diff);

}

compute_AWKS=function(st_a, st_b){

	norm_a=normalize(st_a);
	norm_b=normalize(st_b);

	# Compute average abundances in each group, by treat each sample in group equally
	avg_abd_a=avg_abund(norm_a);
	avg_abd_b=avg_abund(norm_b);

	# Compute average abundances, treat samples in A and B equally.
	comb_avg=(avg_abd_a+avg_abd_b)/2;
	num_categories=ncol(st_a);

	# Compute the KS values for each taxon
	ks_val=numeric(num_categories);
	for(cat_ix in 1:num_categories){
		ks_val[cat_ix]=category_ks_stat(norm_a[,cat_ix], norm_b[,cat_ix]);
	}
	
	# Weighted average of KS values, using average abundance of each taxon as weight
	awks=sum(ks_val*comb_avg); 
	return(awks);
		
}

compute_AWKS_Null=function(norm_combined, totals_combined, num_samples_A, num_samples_B, num_bootstraps){

	num_comb_samples=nrow(norm_combined);
	num_categories=ncol(norm_combined);

	st_a=matrix(0, nrow=num_samples_A, ncol=num_categories);
	st_b=matrix(0, nrow=num_samples_B, ncol=num_categories);
	null_awks=numeric(num_bootstraps);

	for(bs in 1:num_bootstraps){

		# Randomly pick samples
		resamp_ix_a=sample(1:num_comb_samples, num_samples_A, replace=T);
		resamp_ix_b=sample(1:num_comb_samples, num_samples_B, replace=T);

		#cat("Resamples:\n");
		#print(resamp_ix_a);
		#print(resamp_ix_b);

		for(i in 1:num_samples_A){
			# Randomly pick reads
			rs_ix=resamp_ix_a[i];
			st_a[i,]=rmultinom(1, size=totals_combined[rs_ix], prob=norm_combined[rs_ix,]);
		}

		for(i in 1:num_samples_B){
			# Randomly pick reads
			rs_ix=resamp_ix_b[i];
			st_b[i,]=rmultinom(1, size=totals_combined[rs_ix], prob=norm_combined[rs_ix,]);
		}

		# Compute difference between samples
		null_awks[bs]=compute_AWKS(st_a, st_b);

		if(!(bs %% 100)){
			cat(".");
		}
	}

	cat("\n");

	return(null_awks);
}

################################################################################

st_A=load_summary_table(InputFileAName);
st_B=load_summary_table(InputFileBName);

cat("Input Summary Tables:\n");
#print(st_A);
#print(st_B);

norm_A=normalize(st_A);
norm_B=normalize(st_B);

#print(norm_A);
#print(norm_B);

num_categories_A=ncol(st_A);
num_categories_B=ncol(st_B);

cat("Num Taxa: \n");
cat("	File A: ", num_categories_A, "\n");
cat("	File B: ", num_categories_B, "\n");

num_samples_A=nrow(st_A);
num_samples_B=nrow(st_B);

cat("Num Samples: \n");
cat("	File A: ", num_samples_A, "\n");
cat("	File B: ", num_samples_B, "\n");

# Combined samples
combined_st=combine_st(st_A, st_B);
norm_combined_st=normalize(combined_st);
totals_combined=apply(combined_st, 1, sum);

# Use the zero padded samples from the combined table
st_A=combined_st[1:num_samples_A,];
st_B=combined_st[(num_samples_A+1):(num_samples_A+num_samples_B),];

# Compute AWKS between observed
obs_awks=compute_AWKS(st_A, st_B);
cat("Observed AWKS: ", obs_awks, "\n");

# Compuute AWKS null distribution
cat("Computing null distribution for AWKS...\n");
awks_null_val=compute_AWKS_Null(norm_combined_st, totals_combined, num_samples_A, num_samples_B, num_bs);

# Get quantiles for diagnostics
marker_quants=c(.25,.5,.75,.9,.95,.975);
quantiles=quantile(awks_null_val, p=marker_quants);
null_range=range(awks_null_val);

cat("ok...\n");

# Compute p-value
num_gt_obs=(awks_null_val>=obs_awks);
pvalue=sum(num_gt_obs)/num_bs;
cat("p-value: ", pvalue, "\n");

#------------------------------------------------------------------------------

# Plot histogram
max_awks=max(c(awks_null_val, obs_awks));
pdf(paste(OutputRoot, ".awks_test.pdf", sep=""), height=11, width=8.5);
par(mar=c(5.1,4.1,5.1,2.1));
hist_rec=hist(awks_null_val, freq=F, xlim=c(0, max_awks*1.1),  xlab="AWKS values", main="Null Distribution of AWKS Statistic");

# Draw lines
abline(v=obs_awks, col="blue", lwd=3);
q95=quantile(awks_null_val, .95);
abline(v=q95, col="red", lty="dotted");

# Label lines
max_y=max(hist_rec$density);
text(x=obs_awks, y=max_y*.97, pos=4, labels="Observed", col="blue", font=2);
text(x=q95,      y=max_y*.95, pos=4, labels="One-Tailed 95%", col="red", font=3);

# Label pvalue and obs'd
mtext(sprintf("Observed AWKS = %g", obs_awks), side=3, line=0);
mtext(sprintf("p-value = %g", pvalue), side=3, line=1);

#------------------------------------------------------------------------------

fh=file(paste(OutputRoot, ".awks_test.tsv", sep=""), "w");

quantile_str=paste("Null Q(", marker_quants*100, "%)",sep="");

out_str_hdr=paste(c("FileA","FileB","ObsAWKS","P-value","NumBootstraps", "Null Min", quantile_str, "Null Max"), collapse="\t");
out_str_data=paste(c(RootA, RootB, obs_awks, pvalue, num_bs, null_range[1], quantiles, null_range[2]), collapse="\t");
cat(file=fh, out_str_hdr, "\n", sep="");
cat(file=fh, out_str_data, "\n", sep="");

################################################################################

cat("Done.\n")
print(warnings());

q(status=0)
