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

params=c(
	"distance_list", "d", 1, "character",
	"distance_cutoff", "c", 1, "numeric",
	"leader_name", "l", 1, "character",
	"output_root", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n    ", script_name, "\n",
	"	-d <input pairwise distance list>\n",
	"	-c <distance cutoff>\n",
	"	-l <leader sample name>\n",
	"	[-o <output filename root>]\n",
	"\n",
	"Reads in a list of pairwise distances and generates a grouping of in or out\n",
	"If there are redundant pair names, distances will be averaged.\n",
	"\n");

if(!length(opt$distance_list) || !length(opt$distance_cutoff) || !length(opt$leader_name)){
	cat(usage);
	q(status=-1);
}

InputDistanceList=opt$distance_list;
DistanceCutoff=opt$distance_cutoff;
LeaderName=opt$leader_name;

if(length(opt$output_root)){
	OutputFileNameRoot=opt$output_root;
}else{
	clean_leader_name=gsub("/", "-", LeaderName);
	OutputFileNameRoot=paste(clean_leader_name, "_at_", sprintf("%2.1f", DistanceCutoff), sep="");
}

###############################################################################

cat("Input Distance List: ", InputDistanceList, "\n");
cat("Input Distance Cutoff: ", DistanceCutoff, "\n");
cat("Leader Name: ", LeaderName, "\n");
cat("Output Filename Root:", OutputFileNameRoot, "\n");

###############################################################################
# Load data
cat("Loading list...\n");
list=as.matrix(read.table(InputDistanceList, header=FALSE, check.names=FALSE));
list=rbind(list, c(LeaderName,LeaderName,"0"));
numSamples=nrow(list);

distances=as.numeric(list[,3]);


subset_c1=(list[,1]==LeaderName);
subset_c2=(list[,2]==LeaderName);

subset_cb_idx=which(subset_c1 | subset_c2);
num_rel_dist=length(subset_cb_idx);

relevant_dist=distances[subset_cb_idx];
relevant_names=list[subset_cb_idx,1:2];
non_leader_name=character();

for(i in 1:num_rel_dist){
	if(relevant_names[i,1]==LeaderName){
		non_leader_name[i]=relevant_names[i,2];
	}else{
		non_leader_name[i]=relevant_names[i,1];
	}
}

###############################################################################
# Average duplicate data

cat("Averaging duplicated pairwise distances...\n");

#print(non_leader_name);
#print(relevant_dist);

unique_names=unique(non_leader_name);
num_unique=length(unique_names);

avg_non_leader_name=character();
avg_dist=numeric();
for(i in 1:num_unique){
	
	cur_name=unique_names[i];
	dup_idx=(non_leader_name==cur_name);
	avg=mean(relevant_dist[dup_idx]);

	avg_dist[i]=avg;
	avg_non_leader_name[i]=cur_name;
}

#print(avg_non_leader_name);
#print(avg_dist);

###############################################################################
# Determine In/Out nes

cat("Determining In or Outness...\n");

in_grp=avg_dist<=DistanceCutoff;
out_grp=avg_dist>DistanceCutoff;

# Assign to group
group_label=character();
group_label[in_grp]=sprintf("In_%3.1f", DistanceCutoff);
group_label[out_grp]=sprintf("Out_%3.1f", DistanceCutoff);

# Put into easy to output data structure
outmatrix=matrix("", nrow=num_unique, ncol=2);
outmatrix[,1]=group_label;
outmatrix[,2]=avg_non_leader_name;

print(outmatrix);

#print(list);

###############################################################################
# Output results

fh=file(paste(OutputFileNameRoot, ".groups", sep=""), "w");
for(i in 1:num_unique){
	cat(file=fh, paste(outmatrix[i,1], outmatrix[i,2], sep="\t"), "\n", sep="");
}
close(fh);

#-----------------------------------------------------------------------------

fh=file(paste(OutputFileNameRoot, ".oneway_dist", sep=""), "w");
for(i in 1:num_unique){
	cat(file=fh, paste(outmatrix[i,2], avg_dist[i], sep="\t"), "\n", sep="");
}
close(fh);

################################################################################

cat("done.\n\n");
