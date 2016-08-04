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
		"sequence_distance_matrix", "d", 1, "character",
		"sequence_weighting", "w", 2, "character",
		"accession_to_strain_name_map", "m", 2, "character",
		"antigenic_distances", "a", 2, "character",
		"antigenic_distance_cutoff", "c", 2, "numeric",
		"num_bootstraps", "b", 2, "numeric",
		"output_filename_root", "o", 2, "character",
		"antigen_test_list", "t", 2, "character",
		"suppress_summary", "s", 2, "character",
		"accession_clade_map", "l", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];
bin_dir=dirname(script_name);


NUM_BOOTSTRAPS=1000;
AG_DIST_CUTOFF=2;

usage <- paste (
		"\nUsage:\n", script_name, "\n",
		"\n",
		"	sequence options:\n",
		"	-d <sequence-based distance matrix>\n",
		"	[-w <sequence weighting file>]\n",
		"	[-m <accession to strain name mapping>]\n",
		"\n",
		"	antigen options:\n",
		"	-a <antigenic distances>\n",
		"	[-c <coverage cutoff, default = ", AG_DIST_CUTOFF, ">]\n",
		"\n",
		"	output options:\n",
		"	[-b <number of bootstrap, default = ", NUM_BOOTSTRAPS, ">]\n",
		"	[-o <output file name root>]\n",
		"	[-l <accession to clade/group map>]\n",
		"\n",
		"	target options:\n",
		"	[-t <list antigens to test>]\n",
		"	[-s (suppress summary page)]\n",
		"\n",
		"\n",
		"This script will:\n",
		"	1.) Use the sequenced-based distance matrix to build a\n",
		"		Ward's minimum variance hierarchical cluster.\n",
		"	2.) Use the pair-wise list of antigenic distances to\n",
		"		estimate which strains are covered\n",
		"	3.) For each antigen, bootstrap it's sequence coverage\n",
		"		based on bootstrapping which HI information is\n",
		"		involved, and the standard deviation of the distances.\n",
		"\n",
		"The option, -l, is a name of a file which contains an accession \n",
		"	to color/group mapping.\n",
		"\n");

if(!length(opt$sequence_distance_matrix) || !length(opt$antigenic_distances)){
	cat(usage);
	quit(status=-1);	
}

start_time_str=date();

SequenceDistanceMatrix=opt$sequence_distance_matrix;
AntigenicDistances=opt$antigenic_distances;

SequenceWeightingFile="";
if(length(opt$sequence_weighting)){
	SequenceWeightingFile=opt$sequence_weighting;
}

AntigenicDistanceCutoff=AG_DIST_CUTOFF;
if(length(opt$antigenic_distance_cutoff)){
	AntigenicDistanceCutoff=opt$antigenic_distance_cutoff;
}

NumBootstraps=NUM_BOOTSTRAPS;
if(length(opt$num_bootstraps)){
	NumBootstraps=opt$num_bootstraps;
}

OutputFileName=gsub("\\.r_distmat", "", SequenceDistanceMatrix);
if(length(opt$output_file)){
	OutputFileName=opt$output_filename_root;
}

AccessionMapName="";
if(length(opt$accession_to_strain_name_map)){
	AccessionMapName=opt$accession_to_strain_name_map;
}

TestAntigenList="";
if(length(opt$antigen_test_list)){
	TestAntigenList=opt$antigen_test_list;
}

SuppressSummary=FALSE;
if(length(opt$suppress_summary)){
	SuppressSummary=TRUE;
}

AccessionCladeMapName="";
if(length(opt$accession_clade_map)){
	AccessionCladeMapName=opt$accession_clade_map;	
}


#------------------------------------------------------------------------------

cat("\n");
cat("Sequence Distance Matrix: ", SequenceDistanceMatrix, "\n", sep="");
cat("Sequence Weighting File: ", SequenceWeightingFile, "\n", sep="");
cat("Accession Map: ", AccessionMapName, "\n", sep="");
cat("Accession to Clade Map: ", AccessionCladeMapName, "\n", sep="");
cat("\n");
cat("Antigenic Distances: ", AntigenicDistances, "\n", sep="");
cat("Antigenic Distance Cutoff: ", AntigenicDistanceCutoff, "\n", sep="");
cat("\n");
cat("Number of Bootstraps: ", NumBootstraps, "\n", sep="");
cat("Output Filename Root: ", OutputFileName, "\n", sep="");
cat("\n");
cat("Antigen Test List: ", TestAntigenList, "\n", sep="");
cat("Suppress Summary Page: ", SuppressSummary, "\n", sep="");

###############################################################################

cat("Sourcing libraries from: ", bin_dir, "\n");
source(paste(bin_dir, "antigenic_dist_lib.r", sep="/"));
source(paste(bin_dir, "dendrogram_lib.r", sep="/"));
source(paste(bin_dir, "sequence_dist_lib.r", sep="/"));
library(MASS)

cat("\n\n");

###############################################################################

# Sequence-related data load
seq_distmat=load_distance_matrix(SequenceDistanceMatrix);

accession_map=load_map(AccessionMapName, numeric=F);
weighting_map=load_map(SequenceWeightingFile, numeric=T);
clade_map=load_map(AccessionCladeMapName, numeric=F);
#print(weighting_map);

cat("Renaming Distance Matrix:\n");
seq_distmat=rename_distance_matrix(seq_distmat, accession_map);
cat("done.\n");

if(length(weighting_map)>0){
	cat("Renaming Weighting Map...\n");
	weighting_map=rename_weighting_map(weighting_map, accession_map);
	cat("done.\n");
}

if(length(clade_map)>0){
	cat("Renaming Group Map...\n");
	clade_map=rename_weighting_map(clade_map, accession_map);
	cat("done.\n");
}

# HI-related data loads
antigenic_distances_records=load_antigenic_distances(AntigenicDistances);
antigenic_distance_range=range(antigenic_distances_records$Distance);

###############################################################################

# Build the hierarchical clustering for the sequences
cat("Building hierarchical clustering for sequences...\n");
hcl=hclust(as.dist(seq_distmat), method="ward");
sample_names_by_position=hcl$labels[hcl$order];
max_clust_height=max(hcl$height);

# Convert hierarchical cluster to dendrogram
seq_dend=as.dendrogram(hcl);
seq_dend_height=attributes(seq_dend)$height;

###############################################################################

# Rescale labels in dendrogram
# Compute the scale factor of the labels
num_samples=nrow(seq_distmat);

# Plot overall dendrogram
color_map=list();
color_map$IN="green";
color_map$OUT="red";
color_map$UNKNOWN="grey";
color_map$CONFLICT="black";

###############################################################################

# Get list of antigen candidates
if(TestAntigenList==""){
	candidate_antigens=get_list_of_all_antigens(antigenic_distances_records);
}else{
	candidate_antigens=load_list(TestAntigenList);
}
cat("\nCandidate Antigens to process: \n");
print(candidate_antigens);
cat("\n");

num_candidate_antigens=length(candidate_antigens);
cat("Num antigen candidates: ", num_candidate_antigens, "\n");

# Set up output file
dendro_pdf_name=paste(OutputFileName, ".antigenic_coverage.pdf", sep="")
pdf(dendro_pdf_name, height=8.5, width=11);

# Coverage stats
coverage_records=list();
coverage_counts=list();

# Keep track of all In's
candidate_cov_list=list();

global_counts=c(); #=matrix(nrow=num_samples);

for(ag_idx in 1:num_candidate_antigens){

	cur_antigen_candidate=candidate_antigens[ag_idx];
	cat("Estimate coverage for: ", cur_antigen_candidate, "\n");

	# Extract distances by candidate
	distances_from_candidate=extract_distances_by_antigen(cur_antigen_candidate, antigenic_distances_records);
	num_distances=length(distances_from_candidate$Distance);
	inout=list();
	inout$dist_in=sum(distances_from_candidate$Distance<=AntigenicDistanceCutoff);
	inout$dist_out=sum(distances_from_candidate$Distance>AntigenicDistanceCutoff);

	# Compute observed ocverage
	group_map=get_antigens_within_cutoff(AntigenicDistanceCutoff, distances_from_candidate);
	covered_sequences_list=infer_leaf_coverage(group_map, seq_dend, cur_antigen_candidate);
        coverage_counts=accumulate_coverage(covered_sequences_list, coverage_counts);
        obs_perc_coverage=compute_percent_coverage(weighting_map, covered_sequences_list);

	test=vector();
	c_names=vector();
	for (i in 1:length(coverage_counts)) {
		test[i] = unlist(coverage_counts[[i]])[1];
		c_names[i]=strsplit(names(coverage_counts)[i],":")[[1]][2];
		#global_counts=cbind(global_counts,test);
	}
		
	#print(coverage_counts);
	#cat("test vector: ",test, "\n");
	global_counts=rbind(global_counts,test);
	colnames(global_counts) = c_names;
			
	# Bootstrap coverage
	coverage=list();
	coverage$perc_in =numeric(NumBootstraps);
	coverage$perc_out=numeric(NumBootstraps);
	coverage$perc_unk=numeric(NumBootstraps);
	coverage_counts=NULL;
	if(NumBootstraps>=1){
		for(bs in 1:NumBootstraps){

			#cat("Bootstrap: ", bs, "\n");

			# Generate IN/OUT antigen groups
			perturbed_distances=perturb_distances(distances_from_candidate);
			group_map=get_antigens_within_cutoff(AntigenicDistanceCutoff, perturbed_distances);

			# Apply antigens to dendrogram to get sequence coverage
			covered_sequences_list=infer_leaf_coverage(group_map, seq_dend, cur_antigen_candidate);
			coverage_counts=accumulate_coverage(covered_sequences_list, coverage_counts);

			# Apply weights to sequences to get proportion covered
			perc_coverage=compute_percent_coverage(weighting_map, covered_sequences_list);

			# Save each bs result
			coverage$perc_in[bs] =perc_coverage$perc_in;
			coverage$perc_out[bs]=perc_coverage$perc_out;
			coverage$perc_unk[bs]=perc_coverage$perc_unk;
		
			cat(".");
			
		}
	}
	cat("\n");

	# Store the confidence intervals and other
	coverage_records$antigen=cur_antigen_candidate;
	coverage_records$num_distances=num_distances;
	coverage_records$inout=inout;
	coverage_records$observed_coverages=obs_perc_coverage;
	coverage_records$perc_in=get_CI(coverage$perc_in, .05);
	coverage_records$perc_out=get_CI(coverage$perc_out, .05);
	coverage_records$perc_unk=get_CI(coverage$perc_unk, .05);
	candidate_cov_list[[cur_antigen_candidate]]=coverage_records;

	print(coverage_records);

	# Plot the observed (not a bootstrap instance) dendrogram with HI
	group_map=get_antigens_within_cutoff(AntigenicDistanceCutoff, distances_from_candidate);
	observed_seq_dend=dendrapply(seq_dend, mark_leaf_nodes_with_group_id, group_map, cur_antigen_candidate);
	plot_summary_dendrogram(observed_seq_dend, coverage_counts, cur_antigen_candidate, 
		coverage_records, obs_perc_coverage, weighting_map, clade_map, label_scale=-1);


	if(NumBootstraps>=1){
		# Generate plots for coverage statistics
		plot_coverage_statistics(coverage, distances_from_candidate, AntigenicDistanceCutoff, antigenic_distance_range, cur_antigen_candidate);
	}

}


rownames(global_counts) = c(candidate_antigens);
#colnames(global_counts) = c(colnames(seq_distmat));
#colnames(global_counts) = colnames(global_counts, do.NULL=FALSE, prefix="C");
cat("global_counts matrix: ", "\n");
global_counts;
cat("rows: ", nrow(global_counts), "cols: ", ncol(global_counts), "names: ", rownames(global_counts), "\n");
write.table(global_counts, "tmp.csv",row.names=TRUE,col.names=TRUE,sep=",")

if(!SuppressSummary){
	par(mfrow=c(1,1));
	plot_all_relative_values(candidate_cov_list);
}

cov_stat_fname=paste(OutputFileName, ".antigenic_coverage.csv", sep="")
output_coverage_statistics(cov_stat_fname, candidate_cov_list);

# Output timing 
end_time_str=date();
timing=proc.time();
time_str=paste(sprintf("%s: %f", names(timing), timing), collapse="\n");
status_fname=paste(OutputFileName, ".completion", sep="");
stat_fh=file(status_fname, "w");
cat(file=stat_fh, script_name, "\n");
cat(file=stat_fh, "\n");
cat(file=stat_fh, "Start Time: ", start_time_str, "\n", sep="");
cat(file=stat_fh, "  End Time: ", end_time_str, "\n", sep="");
cat(file=stat_fh, "\n");
cat(file=stat_fh, time_str, "\n", sep="");
close(stat_fh);

######################################################################################################

dev.off()
cat("Done.\n");
