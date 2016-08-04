#!/usr/bin/env Rscript

###############################################################################
#                                                                             # 
#       Copyright (c) 2009 J. Craig Venter Institute.                         #     
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
library(MASS);
library(gdata);

params=c(
		"distance_matrix", "d", 1, "character",
		"output_file", "o", 2, "character",
		"antigen_list", "s", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage <- paste (
		"\nUsage:\n", script_name, "\n",
		"	-d <distance matrix>\n",
		"	[-o <output file name>]\n",
		"	[-s <antigen list>]\n",
		"\n",
		"This script will read in a distance matrix (eg. from sequence comparisons)\n",
		"build a Ward's minimum variance clustering, then utilize the accomplice groups\n",
		"to identify which sequence groups may be similar to each other based on the\n",
		"accomplices position to each other.\n",
		"\n");


## Parameters Parsing and Variable Assignment
if(!length(opt$distance_matrix)){
	cat(usage);
	quit(status=0);	
}

DistanceMatrix=opt$distance_matrix;

OutputFileName=gsub("\\.r_distmat", "", DistanceMatrix);
if(length(opt$output_file)){
	OutputFileName=opt$output_file;
}

AntigenList="";
if(length(opt$antigen_list)){
	AntigenList=opt$antigen_list;
}

#------------------------------------------------------------------------------

cat("\n");
cat("Distance Matrix: ", DistanceMatrix, "\n", sep="");
cat("Output Filename Root: ", OutputFileName, "\n", sep="");
cat("Antigen List: ", AntigenList, "\n", sep="");
cat("\n");

## Functions
###############################################################################

load_list=function(filename){
# Loads a list
	members=as.matrix(read.table(filename, sep="\t", header=FALSE));
	return(members[,1]);
}

# Output Confidence Intervals
ci=function(x, alpha){
	n=length(x);

	med=median(x);

	ba=1-(n-2)/n;
	#cat("Best alpha = ", ba, "\n");
	if(ba <= (alpha+.0000000000000001)){
		sorted=sort(x);
		lb=sorted[floor(n*(alpha/2))+1];
		ub=sorted[ceiling(n*(1-(alpha/2)))];
		return(c(med,lb,ub));
	}else{
		return(c(med,NA,NA))
	}
}

######################################################################################################

## MAIN

# Load data
distance_matrix=as.matrix(read.table(DistanceMatrix, header=FALSE, check.names=FALSE));

num_samples=nrow(distance_matrix);

min_d = round(min(as.numeric(distance_matrix[,3])))-.5;
max_d = round(max(as.numeric(distance_matrix[,3])))+.5;

range = max_d - min_d;

b = seq(min_d,max_d,.5);

# Get the sample names
sample_names=c(distance_matrix[,1:2]);

# Load special sample list
if(AntigenList!=""){
	antigens=load_list(AntigenList);
}else{
	antigens=sort(unique(sample_names));
}

cat("Num samples in Distance Matrix:", num_samples, "\n");

######################################################################################################
# Output dendrograms in pdf

dendro_pdf=paste(OutputFileName, ".antigen_analysis.pdf", sep="")
pdf(dendro_pdf, height=8.5, width=11);
fh=file(paste(OutputFileName, ".top_antigens.csv", sep=""), "w");

par(mfrow=c(2,2));
perc=c(0,.25,.50,.75,1);

length(antigens);

samps=c();

for (j in 1:length(antigens)) {
	dist_dat=c();
	for (i in 1:num_samples) {
		if (distance_matrix[i,1]==antigens[j] || distance_matrix[i,2]==antigens[j]) {
			dist_dat <- as.vector(c(dist_dat, distance_matrix[i,3]));
		}
		dist_dat = as.numeric(dist_dat);
	}
	conf_vals = c(ci(dist_dat,.05));
	samps=rbind(samps, c(antigens[j], conf_vals, length(dist_dat)));
	num_assays=length(dist_dat);
	if (!(is.na(conf_vals[2])) && !(is.na(conf_vals[3]))) {
		h = hist(dist_dat, main=paste("Frequency of distances for ", antigens[j]), breaks=b, xlab="Antigenic Distance", axes=FALSE);
		axis(side=1,at=b,labels=sprintf("%3.1f", b), tick=TRUE, cex=.5, las=2);
		#axis(side=3,at=counts[j],labels=sprintf("%3.1f", counts[j]/num_samples*100), tick=TRUE, cex=.5, las=2);
		f_max = max(h$counts);
		ax=seq(0,f_max, f_max/4);
		axis(side=2,at=ax,labels=sprintf("%3.0f", ax), tick=TRUE, cex=.5, las=2);
		#axis(side=2,at=counts[j],labels=sprintf("%3.1f", counts[j]), tick=TRUE, cex=.5, las=2);
		abline(v=conf_vals[1], col="red");
		abline(v=conf_vals[2], col="purple");
		abline(v=conf_vals[3], col="purple");
		cat(file=fh, c(antigens[j], conf_vals, length(dist_dat)), "\n", sep='\t');
	}
}

dev.off()
cat("Done.\n");

