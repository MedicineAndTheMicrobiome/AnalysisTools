#!/usr/bin/env Rscript

###############################################################################

library('getopt');
library(MASS);
options(width=127);
source("wilcoxon_cluster_difference.r");

params=c(
	"distance_matrix", "d", 1, "character",
	"treatments", "t", 1, "character",
	"actions", "a", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"\t\t-d <Distance Matrix>\n",
	"\t\t-t <Treatments>\n",
	"\t\t-a <list of A,B,C,H,I separated by commas>\n",
	"\n",
	"This script will read in the treatments and distance matrix and determine if there\n",
	"is a statistically significant difference between treatments.\n",
	"\n",
	"The first column should always be the sample name\n",
	"\n",
	"A = Average samples together, e.g., for replicates\n",
	"B = Break into independent subgroups.\n",
	"C = Compare matching identifiers, i.e., compare this these to each other. \n",
	"H = Treat as homogenous, i.e., \n",
	"I = Ignore, e.g., descriptions or other free text\n",
	"\n",
	"1.) Using B,C and H build composite key identifiers for each sample.\n",
	"2.) For each composite key, pull up samples.  Those keys with multiple samples will be treated as replicates.\n",
	"3.) Break composite keys by B categories.\n",
	"4.) For each treatment to compare (C), do a pairwise comparison and report stats.\n",
	"\n",
	"Note: That the p-values calculated have not been corrected for multiple comparisions.\n",
	"      You need to adjust your alpha if you want to test the alternative hypothesis that a treatment has had an effect.\n",
	"\n",
	"Also note: That H's and I's can be used multiple times, but A, B, and C can only be used once.\n",
	"\n");

if(!length(opt$distance_matrix) || !length(opt$treatments) || !length(opt$actions)){
	cat(usage);
	q(status=-1);
}

DistanceMatrixFile=opt$distance_matrix;
TreatmentsFile=opt$treatments;
Actions=opt$actions;
OutputFilenameRoot=gsub("\\.dist_mat$", "", DistanceMatrixFile);

cat("Distance Matrix : ", DistanceMatrixFile, "\n", sep="");
cat("Treatments      : ", TreatmentsFile, "\n", sep="");
cat("Actions         : ", Actions, "\n", sep="");
cat("Output File Root: ", OutputFilenameRoot, "\n", sep="");

###############################################################################
# Load distance matrix
cat("\nLoading Distance Matrix...\n");
dist_mat=as.matrix(read.table(DistanceMatrixFile, sep=" ", header=TRUE, check.names=FALSE, row.names=1))

#print(dist_mat[1:10,1:10]);

# Get and report num samples
num_samples=nrow(dist_mat);
cat("Num samples: ", num_samples, "\n", sep="");
if(num_samples != ncol(dist_mat)){
	cat("Num columns: ", ncol(dist_mat), "\n");
	cat("Error, distance matrix is not squared.\n");
	quit(status=-1);
}

###############################################################################
# Load treatment table
cat("\nLoading treatments...\n");
treatments=as.matrix(read.table(TreatmentsFile, sep="\t", header=TRUE, check.names=FALSE, row.names=1));
#print(treatments[1:10,]);
treatment_names=colnames(treatments);
num_treatments=length(treatment_names);
#print(num_treatments);
sample_names=rownames(dist_mat);

###############################################################################
# Load actions
actions=strsplit(Actions, ",")[[1]];
num_actions=length(actions);
#print(actions);

cat("Num Treatments: ", num_treatments, "\n");
cat("Num Actions: ", num_actions, "\n");
if(num_actions != num_treatments){
	cat("Error: Number of specified actions does not match number of treatments.\n");
	cat("\t", num_actions, "!=", num_treatments, "\n");
	quit(status=-1);
}

###############################################################################

# If there is no Break (B) action, then treat all samples as part of All group
treatment_action=cbind(treatment_names, actions);
cat("\nTreatments:\n");
print(treatment_action);

if(!any(actions=="B")){
	treatment_action=rbind(treatment_action, c("All", "B"));
	treatments=cbind(treatments, "All");
	actions=c(actions, "B");
}
print(treatment_action);
print(treatments);

###############################################################################

cat("\nBuilding composite keys (out of attributes marked with action B, C and H):\n");
ms_idx=which(actions=="B" | actions=="C" | actions=="H");
#print(ms_idx);
composite_keys=matrix(character(), nrow=num_samples, ncol=1);
colnames(composite_keys)="CompositeKey";
rownames(composite_keys)=sample_names;
for(i in 1:num_samples){
	composite_keys[i,1]=paste(treatments[i,ms_idx], collapse=".");
}
#cat("Sample ID -> Target Composite Key:\n");
#print(composite_keys);

unique_composite_keys=unique(composite_keys[,1]);
num_composite_keys=length(unique_composite_keys);
cat("\n");
cat("Num composite keys: ", num_composite_keys, "\n");
cat("\n");

treatments=cbind(treatments,composite_keys);
composite_key_col=ncol(treatments);

#cat("Composite Keys:\n");
#print(unique_composite_keys);
#cat("\n");

###############################################################################

#cat("Compressing/Averaging replicates\n");
#print(dist_mat);

sample_names=rownames(dist_mat);
#print(sample_names);

# Compute centroids for each row
for(i in 1:num_composite_keys){
	cur_comp_key=unique_composite_keys[i];
	sample_members=sample_names[which(composite_keys==cur_comp_key)];
	num_members=length(sample_members);
	sample_member_rows=dist_mat[sample_members,];
	if(num_members>1){
		combined=apply(sample_member_rows, 2, mean);
	}else{
		combined=sample_member_rows;
	}
	dist_mat=rbind(dist_mat, combined);
}

# Compute centroids for each column
for(i in 1:num_composite_keys){
	cur_comp_key=unique_composite_keys[i];
	sample_members=sample_names[which(composite_keys==cur_comp_key)];
	num_members=length(sample_members);
	sample_member_col=dist_mat[,sample_members];
	if(num_members>1){
		combined=apply(sample_member_col, 1, mean);
	}else{
		combined=sample_member_col;
	}
	dist_mat=cbind(dist_mat, combined);
}

rownames(dist_mat)=c(sample_names, unique_composite_keys);
colnames(dist_mat)=c(sample_names, unique_composite_keys);
#cat("After:\n");
#print(dist_mat);

###############################################################################

# Compute composite key range
composite_key_range=num_samples+(1:num_composite_keys);

# Set non self-self 0's to a small number
for(i in 1:nrow(dist_mat)){
	for(j in 1:ncol(dist_mat)){
		if(i!=j && dist_mat[i,j]==0){
			dist_mat[i,j]=1e-323;
		}
	}
}

###############################################################################
# Open output files

output_pdf_fname=paste(OutputFilenameRoot, ".pdf", sep="");
pdf(output_pdf_fname, height=8.5, width=11);

output_txt_fname=paste(OutputFilenameRoot, ".results", sep="");
fh=file(output_txt_fname, "w");

###############################################################################

# Extract out distance matrix for the centroids
avgd_distmat=dist_mat[composite_key_range,composite_key_range];

# Plot all samples with and averages
mds=isoMDS(dist_mat);

xlimits=range(mds$points[,1])
ylimits=range(mds$points[,2])
xlimits=c(xlimits[1]*1.1, xlimits[2]*1.1);
ylimits=c(ylimits[1]*1.1, ylimits[2]*1.1);
text_size=42/num_samples;

# Plot all together
plot(mds$points[,1], mds$points[,2], xlim=xlimits, ylim=ylimits, type="n", main="All Samples", xlab="", ylab="");
text(mds$points[,1], mds$points[,2], colnames(dist_mat), col=c(rep(1,num_samples),rep(2,num_composite_keys)), cex=c(rep(text_size*.9,num_samples),rep(text_size*1.2,num_composite_keys)));

# Plot only collapsed replicates
plot(mds$points[composite_key_range,1], mds$points[composite_key_range,2], xlim=xlimits, ylim=ylimits, type="n", main="Collapsed Replicates", xlab="", ylab="");
text(mds$points[composite_key_range,1], mds$points[composite_key_range,2], colnames(dist_mat)[composite_key_range],cex=text_size*1.2, col=2);

#plot(hclust(as.dist(dist_mat)), cex=text_size);

###############################################################################

# Given a list of vectors of names, returns a list of vectors of indices
extract_index_from_distmat_by_sample_name=function(distmat, list_of_sample_lists){
	matrix_names=colnames(distmat);	
	group_names=names(list_of_sample_lists);
	index_list=list();

	# For each group, grab look up it's index by comparing column names
	for(group in group_names){
		sample_list=list_of_sample_lists[[group]];
		index=numeric();
		i=1;
		for(sample_name in sample_list){
			index[i]=which(sample_name==matrix_names);
			i=i+1;
		}
		index_list[[group]]=index;
	}
	return(index_list);
}

# Given a list of indices and a distmat, computes the wilcoxon difference between mean intra and inter distances
compare_groups_pairwise=function(group_indices, avgd_distmat){
	# group_indices: contains a list of vectors, where each vector represents the members in each group
	groups=names(group_indices);
	num_groups=length(groups);

	all_results=list();
	num_comparisons=0;
	
	for(i in 1:num_groups){
		for(j in 1:num_groups){
			if(i<j){
				compare_str=paste(groups[i], " vs ", groups[j], sep="");
				cat("Comparing: ", compare_str, "\n");

				i_vec=as.vector(group_indices[[groups[i]]]);
				j_vec=as.vector(group_indices[[groups[j]]]);

				if(length(i_vec)==0 || length(j_vec)==0){
					cat("Could not compare ", groups[i], " vs ", groups[j], "\n");
				}else{

					# Compute the significance of the two proposed groups
					result=wilcoxon_cluster_difference(i_vec, j_vec, avgd_distmat);

					# Generate some plots
					plot_mds(i_vec, j_vec, avgd_distmat, title=compare_str);
					pval_str=sprintf("Wilcoxon Rank Sum Test, p-value: %3.4f", result$wilcoxon$p.value);
					plot_combo_histogram(result$intra_distances, result$inter_distances, notes=c(pval_str));

					# Index the number of comparisons
					num_comparisons=num_comparisons+1;

					# Build output record for this comparison
					cur_result=list();
					cur_result$comparison=c(groups[i], groups[j]);
					cur_result$p_value=result$wilcoxon$p.value;
					cur_result$a_count=length(i_vec);
					cur_result$b_count=length(j_vec);
					cur_result$median_intra=median(result$intra_distances);
					cur_result$median_inter=median(result$inter_distances);

					all_results[[num_comparisons]]=cur_result;
				}
			}
		}

	}
	cat("Num pairwise comparisons made: ", num_comparisons, "\n");
		
	#print(all_results);
	return(all_results);
}

###############################################################################

output_comparision_results=function(fh, comparison_results){
	num_records=length(comparison_results);
	cat("Outputing ", num_records, " record(s)\n");

	cat(file=fh,
		"A", "B", "p-value", "ACount", "BCount", "MedIntra", "MedInter", "\n", sep=",");

	for(i in 1:num_records){
		rec=comparison_results[[i]];

		cat(file=fh,
			rec$comparison,
			rec$p_value,
			rec$a_count,
			rec$b_count,
			rec$median_intra,
			rec$median_inter,
		"\n", sep=",");
	}
	
}

###############################################################################

a_column_set=(actions=="A");
b_column_set=(actions=="B");
c_column_set=(actions=="C");
h_column_set=(actions=="H");

a_column_idc=which(a_column_set);
b_column_idc=which(b_column_set);
c_column_idc=which(c_column_set);
h_column_idc=which(h_column_set);

a_treatments=unique(treatments[,a_column_idc]);
b_treatments=unique(treatments[,b_column_idc]);
c_treatments=unique(treatments[,c_column_idc]);
h_treatments=unique(treatments[,h_column_idc]);


cat("\n");
cat("Breaking by: ", paste(b_treatments, collapse=", "), "\n");

cat("Composite Key Column: ", composite_key_col, "\n", sep="");
for(break_trt in b_treatments){

	# Output page separator for the treatment that is being broken out
	plot(0,0, xlim=c(-5,5), ylim=c(-5,5), type="n", xlab="", ylab="", xaxt="n", yaxt="n");
	text(0,0, break_trt, cex=3, font=2);

	# Break out samples according to the break out treatment
	breakout_rows=which(treatments[,b_column_idc]==break_trt);
	breakout_samples=treatments[breakout_rows,];
	num_breakout_samples=nrow(breakout_samples);
	cat("\n\tNum (", break_trt, ") samples (replicates not removed): ", num_breakout_samples, "\n", sep="");
	#print(breakout_samples);
	
	# Separate the samples according to the separate treatment
	cat("\n\t\tComparing: ", paste(c_treatments, collapse=","), "\n", sep="");
	compare_sample_list=list();
	for(compare_trt in c_treatments){	
		comparison_rows=which(breakout_samples[,c_column_idc]==compare_trt);
		comparison_samples=breakout_samples[comparison_rows,];
		num_compare_samples=nrow(comparison_samples);
		cat("\t\t\tNum (", compare_trt, ") samples (replicates not removed): ", num_compare_samples, "\n", sep="");
		compare_sample_list[[compare_trt]]=unique(comparison_samples[,composite_key_col]);
	}

	# Extract out relevant columns from distance matrix based on composite key
	#print(compare_sample_list); # Before translation
	group_indices=extract_index_from_distmat_by_sample_name(avgd_distmat, compare_sample_list);
	#print(group_indices); # After translation to index

	# Compare composite keys pairwise
	comparison_results=compare_groups_pairwise(group_indices, avgd_distmat);

	# Output results
	cat(file=fh, "\n", break_trt, "\n");
	output_comparision_results(fh, comparison_results);
}

###############################################################################
cat("\n\n");
print(warnings());
dev.off();
cat("Done.\n");
q(status=0);
