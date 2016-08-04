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

params=c(
		"summary_matrix", "d", 1, "character",
		"output_filename_root", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];
bin_dir=dirname(script_name);

usage <- paste (
		"\nUsage:\n", script_name, "\n",
		"\n",
		"	-d <summary matrix>\n",
		"	[-o <output file name root>]\n",
		"\n",
		"\n");

if(!length(opt$summary_matrix)){
	cat(usage);
	quit(status=-1);	
}

SummaryMatrix=opt$summary_matrix;

OutputFileName=gsub("\\.r_distmat", "", SummaryMatrix);
if(length(opt$output_file)){
	OutputFileName=opt$output_filename_root;
}

#------------------------------------------------------------------------------

cat("\n");
cat("Summary Matrix: ", SummaryMatrix, "\n", sep="");
cat("Output Filename Root: ", OutputFileName, "\n", sep="");
cat("\n");

###############################################################################

cat("Sourcing libraries from: ", bin_dir, "\n");
source(paste(bin_dir, "antigenic_dist_lib.r", sep="/"));
source(paste(bin_dir, "dendrogram_lib.r", sep="/"));
source(paste(bin_dir, "sequence_dist_lib.r", sep="/"));
library(MASS)

cat("\n\n");

###############################################################################

# Set up output file
dendro_pdf_name=paste(OutputFileName, ".antigen_summary.pdf", sep="")
pdf(dendro_pdf_name, height=14, width=9);

# Keep track of all In's
matrix = as.matrix(read.table(SummaryMatrix, sep=",", quote=""));
colnames(matrix)=rep("",ncol(matrix));
candidate_cov_list=list();
coverage_confidence_intervals=list();

for (i in 1:nrow(matrix)) {
	antigen=matrix[i,1];
	coverage_confidence_intervals$antigen=antigen;
	coverage_confidence_intervals$perc_in=c(as.numeric(matrix[i,8]),as.numeric(matrix[i,9]),as.numeric(matrix[i,10]));
	coverage_confidence_intervals$perc_out=c(as.numeric(matrix[i,11]),as.numeric(matrix[i,12]),as.numeric(matrix[i,13]));
	coverage_confidence_intervals$perc_unk=c(as.numeric(matrix[i,14]),as.numeric(matrix[i,15]),as.numeric(matrix[i,16]));
	candidate_cov_list[[antigen]]=coverage_confidence_intervals;
}

par(mfrow=c(1,1));
plot_all_relative_values(candidate_cov_list);


#cov_stat_fname=paste(OutputFileName, ".antigenic_coverage.csv", sep="")
#output_coverage_statistics(cov_stat_fname, candidate_cov_list);

######################################################################################################

dev.off()
cat("Done.\n");
