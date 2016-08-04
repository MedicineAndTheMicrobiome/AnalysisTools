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
library(MASS);

params=c(
	"distance_list", "d", 1, "character",
	"distance_cutoff", "c", 1, "numeric",
	"reference", "l", 1, "character",
	"candidate_fn", "a", 2, "character",
	"output_root", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n    ", script_name, "\n",
	"	-d <input pairwise distance list>\n",
	"	-c <distance cutoff>\n",
	"	-l <leader/reference identifier>\n",
	"	[-a <candidate identifiers>]\n",
	"	[-o <output filename root>]\n",
	"\n",
	"Reads in a pairwise distance list, then based on the leader/reference\n",
	"identifier, computes the quantile of the distance between each candidate\n",
	"and the leader, and also the quantile (proportion) of the antigens\n",
	"that are within the distance cutoff to the leader/reference\n",
	"\n",
	"If the -a option is specified, a separate file with only the candidate antigens\n",
	"specified will be output.\n",
	"\n");

if(!length(opt$distance_list) || !length(opt$distance_cutoff) || !length(opt$reference)){
	cat(usage);
	q(status=-1);
}

InputDistanceList=opt$distance_list;
DistanceCutoff=opt$distance_cutoff;
Reference=opt$reference;
CandidatesFn=opt$candidate_fn;

if(length(opt$output_root)){
	OutputFileNameRoot=opt$output_root;
}else{
	clean_ref_name=gsub("/", "-", Reference);
	clean_input_dist_list=gsub(".dist_flat", "", InputDistanceList);
	OutputFileNameRoot=paste(
		clean_input_dist_list, ".", 
		DistanceCutoff, ".", 
		clean_ref_name, ".",
		"quantiles",
		sep="");
}

###############################################################################

cat("\n");
cat("Input Distance List: ", InputDistanceList, "\n");
cat("Input Distance Cutoff: ", DistanceCutoff, "\n");
cat("Leader/Reference Name: ", Reference, "\n");
cat("Candidate File Name: ", CandidatesFn, "\n");
cat("Output Filename Root:", OutputFileNameRoot, "\n");

###############################################################################
# Load data
cat("\nLoading list...\n");
list=as.matrix(read.table(InputDistanceList, header=FALSE, check.names=FALSE));

# Extract relevent distances which are associated with the reference
colA=which(list[,1]==Reference);
colB=which(list[,2]==Reference);

distances=as.numeric(list[c(colA, colB),3]);
num_distances=length(distances);
antigens=c(list[colA,2], list[colB,1]);

# Collapse/Average down antigen distances with multiple measurements
unique_antigens=unique(antigens);
num_unique_antigens=length(unique_antigens);
averaged_dist=numeric(num_unique_antigens);
names(averaged_dist)=unique_antigens;

# All Raw distances
cat("Raw distances:\n");
for(i in 1:num_distances){
	cat(antigens[i], "\t", distances[i], "\n", sep="");
}
cat("\n");

# Collapse by averaging
for(ag in unique_antigens){
	all_ag_spec_dist=distances[ag==antigens];
	averaged_dist[ag]=mean(all_ag_spec_dist);
}

# Averaged distances
cat("Averaged distances\n");
print(t(t(averaged_dist)));
cat("\n");

# Order everything by averaged distance
sortix=order(averaged_dist);
averaged_dist=averaged_dist[sortix];
unique_antigens=unique_antigens[sortix];

# Compute the quantile for each taxa
quantiles=numeric(num_unique_antigens);
names(quantiles)=unique_antigens;

for(i in 1:num_unique_antigens){
	cur_dist=averaged_dist[i];
	quantiles[i]=100*(1-sum(cur_dist<=averaged_dist)/num_unique_antigens);
}
cat("Quantiles (%):\n");
print(quantiles);

# Compute quantile at the specified cutoff
cutoff_quantile=100*(1-sum(DistanceCutoff<=averaged_dist)/num_unique_antigens);
cat("\n");
cat("Quantile at cutoff: ", cutoff_quantile, "%\n", sep="");

###############################################################################
# Output results

fh=file(paste(OutputFileNameRoot, ".summary.tsv", sep=""), "w");

cat(file=fh, paste("Reference:", Reference, sep="\t"), "\n"); 
cat(file=fh, paste("Total Distances to Reference:", num_distances, sep="\t"), "\n");
cat(file=fh, paste("Total Antigens to Reference:", num_unique_antigens, sep="\t"), "\n");
cat(file=fh, paste("Quantile at Cutoff (", DistanceCutoff, "):\t", cutoff_quantile, sep=""), "\n");
cat(file=fh, "\n");
cat(file=fh, paste("Antigen", "Avg_Dist", "Quantile", sep="\t"), "\n");
for(i in 1:num_unique_antigens){
	cat(file=fh, paste(i, unique_antigens[i], averaged_dist[i], quantiles[i], sep="\t"), "\n");
}

close(fh);

################################################################################

if(length(CandidatesFn)){

	# Load candidates file
	candidates=scan(CandidatesFn, "character");
	print(candidates);

	# Output candidates quantiles only
	fh=file(paste(OutputFileNameRoot, ".candidates.tsv", sep=""), "w");
	for(ag in candidates){
		cat(file=fh, paste(ag, quantiles[ag], sep="\t"), "\n");
	}

	close(fh);
}

################################################################################

cat("done.\n\n");
