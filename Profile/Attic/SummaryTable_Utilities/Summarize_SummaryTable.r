#!/usr/bin/env Rscript

###############################################################################
#                                                                             #
#       Copyright (c) 2013 J. Craig Venter Institute.                         #
#       Copyright (c) 2018 University of Pittsburgh                           #
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
	"input_file", "i", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary_table.tsv>\n",
	"\n",	
	"This script will report a quick summary of the contents of your summary_table.tsv file.\n",
	"\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

###############################################################################

InputFileName=opt$input_file;

cat("\n")
cat("Input File Name: ", InputFileName, "\n\n");

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
smin=min(counts_mat);
smax=max(counts_mat);
smed=median(counts_mat);
smean=mean(counts_mat);
ssd=sd(as.vector(counts_mat));

cat("Num Categories: ", num_categories, "\n");
cat("Top Categories:\n");
for(i in 1:(min(num_categories, 150))){
	cat("\t[", sprintf("%0.4f",nca_sort$x[i]), "]\t", category_names[nca_sort$ix[i]], "\n", sep="");
} 
cat("\n");
cat("Data Depth Min: ", smin, "\n", sep="");
cat("           Max: ", smax, "\n", sep="");
cat("           Med: ", smed, "\n", sep="");
cat("          Mean: ", smean, "\n", sep="");
cat("        St Dev: ", ssd, "\n", sep="");
cat("\n");
cat("Num Zero Count Categories: ", zero_category_totals, "\n");
cat("\n");

st_sort=sort(sample_totals, decreasing=TRUE, index.return=TRUE);
min=min(sample_totals);
max=max(sample_totals);
med=median(sample_totals);
mean=mean(sample_totals);
sd=sd(sample_totals);

cat("Num Samples: ", num_samples, "\n");
cat("Deepest Samples:\n");
for(i in 1:(min(num_samples, 5))){
	cat("\t[", st_sort$x[i], "]\t", sample_names[st_sort$ix[i]], "\n", sep="");
} 

cat("\n");
cat("Sample Depth Min: ", min, "\n", sep="");
cat("             Max: ", max, "\n", sep="");
cat("             Med: ", med, "\n", sep="");
cat("            Mean: ", mean, "\n", sep="");
cat("          St Dev: ", sd, "\n", sep="");





###############################################################################


###############################################################################

cat("Done.\n")
print(warnings());

q(status=0)
