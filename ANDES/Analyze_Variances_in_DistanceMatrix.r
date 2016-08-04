#!/usr/bin/env Rscript

###############################################################################
#                                                                             # 
#       Copyright (c) 2011 J. Craig Venter Institute.                         #     
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
		"input_file", "i", 1, "character",
		"groups_file", "l", 1, "character",
		"output_file", "o", 2, "character",
		"quantile_cutoff", "q", 2, "numeric",
		"accession_to_sample_name_map", "m", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage <- paste (
		"\nUsage:\n", script_name, "\n",
		"	-i <input distance matrix>\n",
		"	-l <group labels file>\n",
		"	[-o <output file name root>]\n",
		"	[-q <quantile cutoff, default >90%>]\n",
		"	[-m <accssion to sample name map>]\n",
		"\n",
		"Analyzes the variances within the supplied distance matrix.\n",
		"	1.) Computes the variance in each group\n",
		"	2.) For each member, computes the proportion of it's distance from the centroid\n",
		"		to the expected variance of the group.\n",
		"\n",
		"The format of the group labels file is 2-column tab-deliminated file:\n",
		"	<group id>\\t<member id>\\n\n",
		"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileRoot=InputFileName;
OutputFileRoot=gsub("\\.r_distmat$","", OutputFileRoot);
OutputFileRoot=gsub("\\.distmat$","", OutputFileRoot);

if(length(opt$output_file)){
	OutputFileRoot=opt$output_file;
}

GroupFilename=character();
GroupProvided=FALSE;
if(length(opt$groups_file)){
	GroupFilename=opt$groups_file;
	GroupProvided=TRUE;
}

cat("Input File: ", InputFileName, "\n");
if(GroupProvided){
	cat("Label File: ", GroupFilename, "\n");
}else{
	cat("No Label File specified.\n");
}

QuantileCutoff=.9;
if(length(opt$quantile_cutoff)){
	QuantileCutoff=opt$quantile_cutoff/100;
}

AccessionMapName="";
if(length(opt$accession_to_sample_name_map)){
        AccessionMapName=opt$accession_to_sample_name_map;
}

OutputFilePDF=paste(OutputFileRoot, ".amova.pdf", sep="");
OutputFilesIntraClusterDeviations=paste(OutputFileRoot, ".intra_cluster_stdev.txt", sep="");
OutputFilesSortedClusterMembers=paste(OutputFileRoot, ".cluster_members_by_dev.txt", sep="");
OutputFilesClosestToCentroid=paste(OutputFileRoot, ".closest_to_centroid.csv", sep="");

###############################################################################

load_map=function(filename){
# Loads a two column matrix so that the first column can be translated into the second column
        members=as.matrix(read.table(filename, sep="\t", header=FALSE));
        num_members=nrow(members);
        cat("Num members in map: ", num_members, "\n");
        map=list();
        for(i in 1:num_members){
                map[[members[i,1]]]=members[i,2];
        }
        return(map);
}

###############################################################################
options(width=200);

# Load the accession to name mapping, if necessary
if(AccessionMapName!=""){
        accession_map=load_map(AccessionMapName);
        map_length=length(accession_map);
        cat("Accession map length: ", map_length, "\n");
}else{
        cat("No accession map file specified.\n");
}

# Load data
dist_mat=as.matrix(read.table(InputFileName, header=TRUE, check.names=FALSE));
num_samples=nrow(dist_mat);

# If accession map provide, rename the distance matrix now.
if(AccessionMapName!=""){
	dup_name=list();
	rname=rownames(dist_mat);
	# First identify duplicate names;
	for(i in 1:num_samples){
		new_name=accession_map[[rname[i]]];
		if(length(new_name)>0){
			if(is.null(dup_name[[new_name]])){
				dup_name[[new_name]]=1;
			}else{
				dup_name[[new_name]]=dup_name[[new_name]]+1;
			}
		}else{
			cat("Could not find name for: '", rname[i], "'\n", sep="");
		}
	}
	# Second, put accession in name if new name is duplicated.
	for(i in 1:num_samples){	
		new_name=accession_map[[rname[i]]];
		if(length(new_name)>0){
			if(dup_name[[new_name]]>1){
				rname[i]=paste(new_name, ":", rname[i], "", sep="");
			}else{
				rname[i]=new_name;
			}
		}
	}	
	rownames(dist_mat)=rname;
	colnames(dist_mat)=rname;
}

# Compute label scales based on number of samples
sample_names=rownames(dist_mat);
label_scale=63/num_samples;
if(label_scale > 2){
	label_scale=2;
}

group_to_sample_map=list();
members_per_group=numeric();
# Assign colors to groups by sample id
if (GroupProvided) {
	# Load cluster assignments
	cat("Loading: ", GroupFilename, "\n");
	cluster_assignments=as.matrix(read.table(GroupFilename, header=FALSE, check.names=FALSE, sep="\t"));
	rownames(cluster_assignments)=NULL;
	colnames(cluster_assignments)=NULL;
	num_assignments=nrow(cluster_assignments);

	# First column is group, second column is sample id
	group_names=unique(cluster_assignments[,1])
	num_groups=length(group_names);

	# Assign members to group
	sample_to_group_map=list()
	for(i in 1:num_assignments) {
		sample_to_group_map[cluster_assignments[i,2]]=cluster_assignments[i,1];
	}
	
	# Assign group members
	for(i in 1:num_assignments){
		group_id=cluster_assignments[i,1];
		if(is.null(group_to_sample_map[[group_id]])){
			group_to_sample_map[[group_id]]=character();
		}
		group_to_sample_map[[group_id]]=c(group_to_sample_map[[group_id]], cluster_assignments[i,2]);
	}

	# Compute num members per group
	for(i in 1:num_groups){
		members_per_group[i]=length(group_to_sample_map[[i]]);
	}
}else{
	# If no groups specified, then assume it is a group of 1.
	num_groups=1;
	sample_to_group_map=list();
	group_names="Default";
	sample_names=colnames(dist_mat);
	for(i in 1:num_samples){
		sample_to_group_map[sample_names]="Default";
		group_to_sample_map[["Default"]]=sample_names;
	}
	members_per_group[1]=num_samples;
}

######################################################################################################

names(members_per_group)=group_names;

cat("Num Groups: ", num_groups, "\n");
print(group_names);
cat("\n");

cat("Members per group: \n");
print(members_per_group);
cat("\n");

######################################################################################################

if(num_groups>1){
	# Allocate colors
	rainbow_natural=rainbow(num_groups)
	rainbow_lowval=rainbow(num_groups,v=.5)
	rainbow_lowsat=rainbow(num_groups,s=.75)
	odd=seq(1,num_groups,2);
	even=seq(2,num_groups,2);
	rainbow_colors=character();
	rainbow_colors[odd]=rainbow_natural[odd];
	rainbow_colors[even]=rainbow_lowsat[even];
	yll=round((.5/6)*(num_groups))+1;
	yul=round((1.5/6)*(num_groups))+1;
	yellows=seq(yll,yul,1);
	rainbow_colors[yellows]=rainbow_lowval[yellows]

	# Shuffle colors
	matrix_dim=ceiling(sqrt(num_groups));
	shuf_vect=rep(0,matrix_dim^2);
	shuf_vect[1:num_groups]=1:num_groups;
	shuf_vect=as.vector(t(matrix(shuf_vect,ncol=matrix_dim)));
	shuf_vect=shuf_vect[shuf_vect>0];
	rainbow_colors=rainbow_colors[shuf_vect];

	# Assign colors to groups
	group_to_color_map=list();
	for(i in 1:num_groups) {
		group_to_color_map[[group_names[i]]]=rainbow_colors[i];
	}
}else{
	rainbow_colors="black";
}

######################################################################################################
# Compute the centroids for each group

centroid_matrix=matrix(0, ncol=num_samples, nrow=num_groups);
rownames(centroid_matrix)=group_names;
colnames(centroid_matrix)=sample_names;
for(i in 1:num_groups){
	cur_group=group_names[i];
	#cat("Working on ", cur_group, "\n", sep="");	
	samples_in_group=group_to_sample_map[[cur_group]];
	num_members=length(samples_in_group);

	# The centroid is just the mean of each of the distances
	sum=numeric(num_samples);
	for(sample_idx in 1:num_members){
		cur_sample=samples_in_group[sample_idx];
		sum=sum+dist_mat[cur_sample,];
	}
	centroid_matrix[i,]=sum/num_members;
}
cat("\n");
#print(centroid_matrix);

######################################################################################################
# Compute the variance for each group.  Make sure we only include the members that are in the group

group_variance=numeric(); # per group: sum of square distances from cnetroid, divided by size of group
per_group_sqrDev=list();  # per group: list of squares distances from centroid
sample_centroid_sqr_dev=list(); # For each sample: distance from it's group centroid

for(i in 1:num_groups){
	cur_group=group_names[i];
        #cat("Working on ", cur_group, "\n", sep="");
	samples_in_group=group_to_sample_map[[cur_group]];
        num_members=length(samples_in_group);

	sum_sqr=0;
	per_group_sqrDev[[cur_group]]=numeric(); # Per group, sum of squares Model/Within
	for(sample_id in samples_in_group){
		sqr_dev_from_centroid=centroid_matrix[cur_group, sample_id]^2;
		# For variance calculation:
		sum_sqr=sum_sqr+sqr_dev_from_centroid;
		# For per group analysis:
		per_group_sqrDev[[cur_group]]=c(per_group_sqrDev[[cur_group]], sqr_dev_from_centroid);
		# For global analysis:
		sample_centroid_sqr_dev[[sample_id]]=sqr_dev_from_centroid;
	}
	names(per_group_sqrDev[[cur_group]])=samples_in_group;
	group_variance[cur_group]=sum_sqr/num_members;
}

#cat("\n");
#print(group_variance);
group_stddev=numeric();
group_stddev=sqrt(group_variance);
#cat("\n");
#cat("Group standard deviations:\n");
#print(group_stddev);

######################################################################################################
# Compute the SSE for this dataset, and use that as our NULL distribution.

# Number of distances to generate
bs_num_samples_per_group=80;
num_distances=num_groups*bs_num_samples_per_group;
cat("Num distances to bootstrap: ", num_distances, "\n");
null_distrib=c();
bs_idx=1;
for(grp_id in 1:num_groups){
	null_distrib=rbind(null_distrib,
			sample(per_group_sqrDev[[grp_id]],bs_num_samples_per_group,replace=TRUE)
			);
}
null_sorted=sort(null_distrib);
quantiles=c(.50,.75,.90,.95);
num_null_samples=length(null_distrib);
quantile_offsets=quantiles*num_null_samples;
quantile_cutoff_value=null_sorted[QuantileCutoff*num_null_samples];

#print(null_distrib);

######################################################################################################
# Compute quantile for each sample based on null distribution

sorted_null_distrib=sort(null_distrib, decreasing=FALSE);

sample_quantile=list();
for(sample_id in sample_names){
	sample_sqr_dev=sample_centroid_sqr_dev[[sample_id]];
	sample_quantile[[sample_id]]=sum(sorted_null_distrib<sample_sqr_dev)/num_distances;
}
#print(sample_quantile);

######################################################################################################
# Compute the MDS per group
mds_list_result=list();
pca_list_result=list();

mds_point_limits=list();
mds_point_limits$min_x=Inf;
mds_point_limits$min_y=Inf;
mds_point_limits$max_x=-Inf;
mds_point_limits$max_y=-Inf;

for(grp_id in 1:num_groups){

	cur_group=group_names[grp_id];
	cat("Preparing group specific distance matrix: ", cur_group, "\n");

	# Extract out group specific distance matrix
	samples_in_group=group_to_sample_map[[cur_group]];
	num_samples_in_group=length(samples_in_group);

	if(num_samples_in_group<2){
		cat("Num samples in group: ", num_samples_in_group, " too small to perform analysis.\n");
		break;		
	}

	grp_dist_mat=matrix(0, ncol=num_samples_in_group, nrow=num_samples_in_group);
	colnames(grp_dist_mat)=samples_in_group;
	rownames(grp_dist_mat)=samples_in_group;
	#grp_dist_mat[grp_dist_mat<0]<-10^-300;
	for(a_idx in 2:num_samples_in_group){
		for(b_idx in 1:(a_idx-1)){
			dist_val=dist_mat[samples_in_group[a_idx], samples_in_group[b_idx]];
			if(dist_val==0){
				dist_val=1e-315;
			}
			grp_dist_mat[a_idx,b_idx]=dist_val;
		}
	}
	
	grp_dist_mat=as.dist(grp_dist_mat);
	grp_dist_mat[grp_dist_mat<0]<-1e-300;
	# Compute MDS
	if(num_samples_in_group>3){

		if(min(grp_dist_mat)==max(grp_dist_mat) && max(grp_dist_mat)==1e-315){
			mds_result=NULL;
		}else{
			# Do a prelim mds scale to determine if there is enough dimensionality to the data
			test_mds=cmdscale(grp_dist_mat,k=3);
			if(ncol(test_mds)<3 || all(is.nan(test_mds[,3]))){
				# Not enough complexity to do 3D, do 2D, and fill last dimension with 0's
				mds_result=isoMDS(grp_dist_mat,k=2);
				mds_result$points=cbind(mds_result$points, 0);
			}else{
				#print(grp_dist_mat);
				#grp_dist_mat[grp_dist_mat<0]<-10^-300;
				mds_result=isoMDS(grp_dist_mat,k=3);
			}
		}

		cat("Setting group id: ", grp_id, "\n");
		mds_list_result[[as.character(grp_id)]]=mds_result;

		mds_point_limits$min_x=min(mds_point_limits$min_x, mds_result$points[,1])
		mds_point_limits$min_y=min(mds_point_limits$min_y, mds_result$points[,2])
		mds_point_limits$max_x=max(mds_point_limits$max_x, mds_result$points[,1])
		mds_point_limits$max_y=max(mds_point_limits$max_y, mds_result$points[,2])
	
	}

}
#print( mds_point_limits);


######################################################################################################

transform_range=function(x, new_min, new_max){

	x_min=min(x);
	x_max=max(x);
	x_toZero=x-x_min;

	x_range=x_max-x_min;
	x_norm=x_toZero/x_range;

	new_range=new_max-new_min;

	new_x=x_norm*new_range;
	new_x=new_x+new_min;

	if(all(is.nan(new_x))){
		new_x=rep(1,length(new_x));
	}

	return(new_x);
}

######################################################################################################
######################################################################################################

pdf(OutputFilePDF, width=11, height=8.5);

#----------------------------------------------------------------------------------------------------

# Plot bar chart of standard deviations
par(mfrow=c(1,1));
barpos=barplot(group_stddev,
	col=rainbow_colors, 
	ylim=c(0,max(group_stddev)*1.2),
	main=OutputFileRoot, ylab="Percent Difference among Sequences", xlab="Group ID");
mtext("Intracluster Standard Deviations");
text(barpos,group_stddev, labels=sprintf("%3.4f", group_stddev), pos=3, cex=.8);

#----------------------------------------------------------------------------------------------------

# Plot histograms for the squared distances from the centroid for each group
par(mfrow=c(1,3));
par(oma=c(.5,.5,2,.5));
max_variance=max(null_distrib);
hist_breaks=seq(0,max_variance*1.5,length.out=20);
for(grp_id in 1:num_groups){

	group_name=group_names[grp_id];
	variances=per_group_sqrDev[[group_name]];
	num_members=length(variances);

	
	cat("BREAKS:",hist_breaks,"\nVAR:",max_variance,"\n");
	# Plot Histograph
	par(mar=c(3,3,3,.5));
	hist_rec=hist(variances, 
		xlim=c(0,max_variance), breaks=hist_breaks, 
		main=group_names[grp_id], xlab="Squared Distance from Group Centroid", col=rainbow_colors[grp_id]);
	mtext(sprintf("Num Members = %i", num_members), line=-1, cex=.8);
	abline(v=quantile_cutoff_value, col="red", lty=3);
	text(quantile_cutoff_value, max(hist_rec$counts), sprintf("%3.1f%%", QuantileCutoff*100),col="red", cex=.8, pos=4);

	# Plot MDS in "3D"
	point_min=.5;
	point_max=3;
	target_range=range(c(point_min, point_max));
	par(mar=c(3,1,3,0));
	plot(0,0,
		xlim=c(mds_point_limits$min_x, mds_point_limits$max_x),
		ylim=c(mds_point_limits$min_y, mds_point_limits$max_y),
		xlab="Dimension 1",
		ylab="Dimension 2", type="n",
		main=OutputFileRoot
	);

	grp_id=as.character(grp_id);
	if(!is.null(mds_list_result[[grp_id]])){
		pts=mds_list_result[[grp_id]]$points;
		num_pts=nrow(pts);
		labels=rownames(pts);

		# Plot points, using cex to represent Z axis
		point_sizes=transform_range(c(pts[,3],0), point_min, point_max);
		points(0,0, pch=16, col="grey", cex=tail(point_sizes,1));
		points(pts[,1],pts[,2], cex=point_sizes);

		for(pt_idx in 1:num_pts){
			#text(pts[pt_idx,1],pts[pt_idx,2], labels=labels[pt_idx], pos=1, cex=.5);

			text(pts[pt_idx,1],pts[pt_idx,2], 
				labels=sprintf("%3.1f",100*sample_quantile[[labels[pt_idx]]]), cex=.8, pos=3);

			if(sample_quantile[[labels[pt_idx]]]>=QuantileCutoff){
				points(pts[pt_idx,1],pts[pt_idx,2], pch=4, cex=point_sizes[pt_idx], col="red");
			}
		}
	}else{
		text(0,0, "No MDS plot for <4 points.\n");
		cat("Null MDS calculation for: ", group_names[grp_id], "\n");
	}


	# Plot list of names
	par(mar=c(3,0,3,.5));
	plot(0,0, xlim=c(-5,100), ylim=c(-100, 0), ylab="", xlab="", type="n", xaxt="n", yaxt="n", bty="n");
	par(family="Courier");
	text(0, 0, " n   Qntl   SampleName", pos=4, font=2);
	sorted_variances=sort(variances, decreasing=FALSE);
	for(pt_idx in 1:num_members){
		sample_id=names(sorted_variances)[pt_idx];
		quantile=sprintf("%3.1f",100*sample_quantile[[sample_id]])
		index=sprintf("%2i", pt_idx);
		line_str=paste(index, " ", quantile, " ", sample_id);
		text(0, -(pt_idx+1)*1.5, line_str, pos=4);
	}
	par(family="");

}

#----------------------------------------------------------------------------------------------------
# Plot the overall distribution
par(mfrow=c(1,1));
par(oma=c(1,1,2,.5));
par(mar=c(4,4,2,.5));
cat("Plotting the bootstrapped squared distances distribution.\n");
hist_rec=hist((as.vector(null_distrib)), 
	xlim=c(0,(max_variance)),
	main=OutputFileRoot, breaks=seq(0, (max_variance), length.out=20), xlab="Square Distances from Each Group Centroid");
mtext("Global Intracluster Distances", line=-.5, cex=.75);
abline(v=null_sorted[quantile_offsets],col="grey");
text(null_sorted[quantile_offsets], rep(max(hist_rec$counts),length(quantiles)), sprintf("%i%%",quantiles*100), 
	cex=.8, col="grey", pos=4);
abline(v=quantile_cutoff_value, col="red", lty=3);
text(quantile_cutoff_value, max(hist_rec$counts), sprintf("%i%%", QuantileCutoff*100),col="red", cex=.8, pos=4);

######################################################################################################
# Output standard deviations for each group.

cat("Writing intra cluster standard deviations...\n");
fh=file(OutputFilesIntraClusterDeviations, "w");
for(grp_id in 1:num_groups){
	cat(file=fh, group_names[grp_id], group_stddev[grp_id], sep="\t");
	cat(file=fh, "\n");
}
close(fh);

#----------------------------------------------------------------------------------------------------
# Output members by variance

cat("Writing cluster members...\n");
fh=file(OutputFilesSortedClusterMembers, "w");
cat(file=fh,"Group_ID", "\t", "Sample_ID", "\t", "Sqrd_Dev", "\t", "Contrib_To_Variance", "\t", "Global_Quantile", "\n");
for(grp_id in 1:num_groups){
	group_name=group_names[grp_id];

	sqr_deviations=per_group_sqrDev[[group_name]];
	num_members=length(sqr_deviations);
	sorted_sqr_dev=sort(sqr_deviations, decreasing=TRUE);
	sorted_member_names=names(sorted_sqr_dev);
	sum_sqr_dev=sum(sqr_deviations);
	contrib_to_var=sorted_sqr_dev/sum_sqr_dev;

	cat(file=fh, group_name, "\n", sep="");
	for(member_idx in 1:num_members){
		sample_name=sorted_member_names[member_idx];
		global_quantile=sample_quantile[[sample_name]]*100;
		cat(file=fh, "\t", sample_name, "\t", sorted_sqr_dev[member_idx], "\t", contrib_to_var[member_idx], "\t", sprintf("%3.1f", global_quantile), "\n", sep="");
	}	
	total_variance=sum(variances);
}
close(fh);

#----------------------------------------------------------------------------------------------------
# Output member closest to centroid

cat("Outputing closest member to centroid of each group...\n");
fh=file(OutputFilesClosestToCentroid, "w");
cat(file=fh, paste(c("#Group_ID", "Member_Name", "Distance_to_Centroid", "Global_Quantile"), collapse=","));
cat(file=fh, "\n");
for(grp_id in 1:num_groups){
	group_name=group_names[grp_id];
	sqr_deviations=per_group_sqrDev[[group_name]];
	sorted_sqrdeviations=sort(sqr_deviations, decreasing=FALSE);
	closest_name=names(sorted_sqrdeviations)[1];	
	closest_distance=sqrt(sorted_sqrdeviations[1]);
	global_quantile=sprintf("%3.1f", sample_quantile[[closest_name]]*100);
	cat(file=fh, 
		paste(c(group_name, closest_name, sprintf("%3.4f", closest_distance[1]), global_quantile), collapse=","),
		"\n");
}
close(fh);

######################################################################################################

dev.off()
cat("Printing warnings.\n");
print(warnings());
cat("Done.\n");
