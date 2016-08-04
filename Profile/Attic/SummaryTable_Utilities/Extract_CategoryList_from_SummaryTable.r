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
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	[-o <output list name>]\n",
	"\n",	
	"This script will place all the columns/category names in your\n",
	"summary table into a list.\n",
	"\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

###############################################################################

InputFileName=opt$input_file;
OutputFileName=opt$output_file;

if(length(OutputFileName)==0){
	OutputFileName=InputFileName;
	OutputFileName=gsub("\\.summary_table\\.tsv$", "", OutputFileName);
	OutputFileName=gsub("\\.summary_table\\.xls$", "", OutputFileName);
	OutputFileName=paste(OutputFileName, ".txt", sep="");
}

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");

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

###############################################################################
# Normalize

# Sum sample totals
sample_totals=apply(counts_mat, 1, sum);
zero_category_totals=sum(sample_totals==0);
#print(sample_totals);

# normalize, to compute probabilities
normalized=matrix(0, nrow=num_samples, ncol=num_categories);
for(i in 1:num_samples){
	normalized[i,]=counts_mat[i,]/sample_totals[i];
}
#print(normalized);

normalized_category_average=apply(normalized, 2, mean);
#print(normalized_category_average);

nca_sort=sort(normalized_category_average, decreasing=TRUE, index.return=TRUE);

cat("Num Categories: ", num_categories, "\n");

fh=file(OutputFileName, "w");
for(i in 1:num_categories){
	cat(file=fh, category_names[nca_sort$ix[i]], "\n", sep="");
} 

###############################################################################

cat("Done.\n")
print(warnings());

q(status=0)
