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
	"num_cat", "p", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	[-p <number of top categories to search through, default is all>]\n",
	"\n",	
	"This script assumes that each sample has a special\n",
	"characteristic in some set of categories, but\n",
	"most of the categories are not special.\n",
	"\n",
	"For each sample it will try to determine which categories\n",
	"are over or under represented, relative to the remaining\n",
	"samples.  Essentially, the non-scrutinized samples are\n",
	"used as a null distribution.\n",
	"\n",
	"The comparisons are made in log(a/(1-a)) space.\n",
	"\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

###############################################################################

InputFileName=opt$input_file;

NumSearchCat=Inf;
if(length(opt$num_cat)){
	NumSearchCat=opt$num_cat;
}

cat("\n");
cat("Parameters:\n");
cat("  Input file name: ", InputFileName, "\n");
cat("  Number of categories targeted: ", NumSearchCat, "\n");
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

cat("Summary Table: \n");
cat("  Num Categories: ", num_categories, "\n");
cat("  Num Samples:    ", num_samples, "\n");
cat("\n");

if(NumSearchCat>num_categories){
	NumSearchCat=num_categories;
} 

cat("Num Categories to Analyze: ", NumSearchCat, "\n");
cat("\n");

###############################################################################

# Normalize
cat("Normalizing individual samples...\n");
sample_totals=apply(counts_mat, 1, sum);
normalized=matrix(0, nrow=num_samples, ncol=num_categories);
colnames(normalized)=category_names;
rownames(normalized)=sample_names;
for(i in 1:num_samples){
	normalized[i,]=counts_mat[i,]/sample_totals[i];
}
#print(normalized);

overall_avg=apply(normalized, 2, mean);
#print(overall_avg);

# Order from greatest to smallest abundance
oa_sort=order(overall_avg, decreasing=TRUE);
#print(overall_avg[oa_sort]);

# Reorder the original data
category_names=category_names[oa_sort];
normalized=normalized[,oa_sort];
counts_mat=counts_mat[,oa_sort];
#print(normalized);

# Move Remainder category to end of table
uc_category_names=toupper(category_names);
which_isRemaining="REMAINING"==uc_category_names;
if(any(which_isRemaining)){
	if(sum(which_isRemaining)>1){
		cat("Error:  More than one Remaining category.\n");
		quit(-1);
	}else{
		cat("Remaining category identified.\n");
		rem_ix=which(which_isRemaining);
		new_order=c((1:num_categories)[-rem_ix],rem_ix); 
		category_names=category_names[new_order];
		normalized=normalized[,new_order];
		counts_mat=counts_mat[,new_order];
	}
}
#print(normalized);

###############################################################################

logodd_mat=log(normalized/(1-normalized));

###############################################################################

for(i in 1:(NumSearchCat)){
	cat_name=category_names[i];
	cat(i, ".) Category: ", cat_name, "\n", sep="");
	val=logodd_mat[,i];
	print(val);
	h=hist(val, main=paste(i, ".) ", cat_name, sep=""), breaks=nclass.Sturges(val)*3);

	breaks=h$breaks;
	mids=h$mids;
	for(bix in 0:((length(breaks)-1))){
		inbin=(val>breaks[bix]) & (val<breaks[bix+1]);
		if(sum(inbin)==1){
			print(sample_names[inbin]);
			text(mids[bix], 1, sample_names[inbin], adj=c(-.5,0), srt=90, cex=.75);
		}
	}

	cat("\n");
}

###############################################################################

cat("Done.\n")
print(warnings());

q(status=0)
