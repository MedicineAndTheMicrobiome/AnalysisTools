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
# Get parameters from the user/command line

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"group_regex", "r", 1, "character",
	"output_file", "o", 1, "character",
	"outlier_threshold", "t", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
        "\nUsage:\n   ", script_name, "\n",
        "       -i <input distance matrix>\n",
        "       -o <output file name root>\n",
	"	-r <group regular expression>\n",
	"	[-t <outlier threshold, eg. 95 is 95%]\n",
        "\n",
	"Reads in a distance matrix and generates an MDS and dendrogram.",
	" The regular expression could look something like: -r '\\.z..(.*$)'\n",
	"\n",
        "\n");

if(!length(opt$input_file) || !length(opt$output_file) || !length(opt$group_regex)){
        cat(usage);
        q(status=-1);
}
InputFileName=opt$input_file;
GroupRegex=opt$group_regex;
OutputFileName=opt$output_file;
OutlierThreshold=opt$outlier_threshold;

if(length(OutlierThreshold)==0){
	OutlierThreshold=95;	
}

###############################################################################

library(MASS)

cat("Input File:", InputFileName, "\n");
cat("Outlier Threshold:", OutlierThreshold, "\n");

# Load data
A=as.matrix(read.table(InputFileName, sep=",", header=TRUE, check.names=FALSE, row.names=1))
#print(A);

#regex="\\.(z.+)$"; # sites
#regex="\\.(z.)"; # regions
#regex="\\.(\\d+)"; # donors

names=names(A[1,]);
num_samples=nrow(A);
cat("Num samples: ", num_samples, "\n", sep="");

###############################################################################
# Determine group memberships

# Try to parse out the regular expression in order to assign each sample to a group
group=character(num_samples);
for(i in 1:num_samples){
	pos=regexpr(GroupRegex, names[i], perl=TRUE);
	#print(pos);
	group[i]=substring(names[i], first=pos+1, last=pos+attr(pos,"match.length")-1);
	#cat(names[i], " -> ", group[i], "\n", sep="");
}

# Determine which groups were detected
unique_groups=unique(group);
cat("Detected groups (based on your regex: ", GroupRegex, ")\n", sep="");
print(unique_groups);
num_unique_groups=length(unique_groups);

# Keep track of indices into the matrix for each group
group_indices=list();
for(i in 1:num_unique_groups){
	group_indices[[i]]=numeric();
	for(j in 1:num_samples){
		if(group[j]==unique_groups[i]){
			group_indices[[i]]=c(group_indices[[i]], j);
		}
	}
}
#print(group_indices);

###############################################################################
# Compute centroids for each group
centr_dist=matrix(0, nrow=num_unique_groups, ncol=num_samples);
for(centr_id in 1:num_unique_groups){
	for(sample_id in 1:num_samples){
		#cat("Centroid ID: ", centr_id, "\n", sep="");
		#cat("Sample ID: ", sample_id, "\n", sep="");
		centr_dist[centr_id,sample_id]=mean(A[group_indices[[centr_id]],sample_id]);
	}
}
#print(centr_dist);

# Compute centroid-centroid distances
centr_centr_dist=matrix(0, nrow=num_unique_groups, ncol=num_unique_groups);
for(centr_id1 in 1:num_unique_groups){
	for(centr_id2 in 1:num_unique_groups){
		centr_centr_dist[centr_id1, centr_id2]=mean(centr_dist[centr_id1, group_indices[[centr_id2]]]);
	}
}

# Insert centroid info into distance matrix

# Insert bottom
A=rbind(A,centr_dist); 	

# Insert right columns, but don't worry about setting distances, we will only use bottom left anyway
A=cbind(A,matrix(0,nrow=nrow(A), ncol=num_unique_groups));
centroid_matrix_pos=numeric(0);
for(centr_id1 in 1:num_unique_groups){
	for(centr_id2 in 1:num_unique_groups){
		A[num_samples+centr_id1, num_samples+centr_id2]=centr_centr_dist[centr_id1,centr_id2];
	}
}
rownames(A)=c(names, unique_groups);
names=rownames(A);

total_rows_w_centroids=nrow(A);
total_cols_w_centroids=ncol(A);

cat("Rows w/ centroids: ", total_rows_w_centroids, "\n", sep="");
cat("Cols w/ centroids: ", total_cols_w_centroids, "\n", sep="");

###############################################################################

centr_idx=(num_samples+1):(num_samples+num_unique_groups);
sample_idx=seq(1,num_samples);

cat("Sample indices: ", min(sample_idx), "-", max(sample_idx), "\n", sep="");
cat("Centroid indices: ", paste(centr_idx, collapse=","), "\n", sep="");

###############################################################################
# Assign colors to groups and samples

# Create color palette
predef_colors=c("red", "blue", "green", "orange", "purple", "brown", "black", "pink", "gray");
if(num_unique_groups<=length(predef_colors)){
	colors=predef_colors;
}else{
	colors=rainbow(num_unique_groups);
}

# Assign colors to each group by creating a group->color map
color_map=list(0);
for(i in 1:num_unique_groups){
	color_map[[unique_groups[i]]]=colors[i];
}

# Assign colors to each sample by going from group(sample)->color
assigned_colors=character(0);
for(i in 1:num_samples){
	assigned_colors[i]=color_map[[group[i]]];
}

# Assign black to centroids entries
assigned_colors=c(assigned_colors, rep("black",num_unique_groups));
#print(assigned_colors);

# Compute label scales based on number of samples
numSamples=nrow(A);
label_scale=35/numSamples;
if(label_scale > 1){
       label_scale=1;
 }
#label_scale=.3;

##################################################################################
# Detect outliers

cutoff=OutlierThreshold/100.0;
#cutoff=.8;

outliers_in_matrix_idx=numeric(0);
outliers_by_group=list();
for(group_id in 1:num_unique_groups){
	cat("Group: ", unique_groups[group_id], "\n");
	#print(group_indices[[group_id]]);
	dist_to_centroid=A[centr_idx[group_id], group_indices[[group_id]]];
	#print(dist_to_centroid);
	sort_info=sort(dist_to_centroid, index.return=TRUE);
	num_members=length(group_indices[[group_id]]);

	lb=1+as.integer(ceiling(cutoff*num_members));
	ub=num_members;
	if(lb>ub){
		lb=ub;
	}
	outliers=lb:ub;

	#print(outliers);
	num_outliers=length(outliers);
	cat("Num outliers: ", num_outliers, " / Num members: ", num_members, "\n", sep="");

	cat("outliers: ", outliers, "\n");
	print(sort_info$x);
	outlier_matrix_idx=group_indices[[group_id]][sort_info$ix[outliers]];
	print(outlier_matrix_idx);
	#print(outliers);
	#print(outlier_matrix_idx);
	outliers_in_matrix_idx=c(outliers_in_matrix_idx, outlier_matrix_idx);
	outliers_by_group[[group_id]]=outlier_matrix_idx;
	cat("\n\n\n");
}


##################################################################################
pdf(paste(OutputFileName, ".pdf", sep=""),width=11,height=8.5)

for(i in 1:nrow(A)){
        for(j in 1:ncol(A)){
                if(A[i,j]==0){
                        A[i,j]=1e-99;
                }
        }
}

dist=as.dist(A);
#print(dist);

##################################################################################
# Compute MDS
mds<-isoMDS(dist)


# Plot the overall MDS
# Compute the padding inside of the graph so the text labels don't fall off the screen
min<-min(mds$points[,1])
max<-max(mds$points[,1])
width=(max-min)
margin=width*(0.1)

# Create empty plot so that the ranges are good to go 
plot(mds$points,type="n", main="All", xlim=c(min-margin, max+margin), xlab="", ylab="")

# Plot X's under outliers
points(mds$points[outliers_in_matrix_idx,], col="grey20", pch=13, cex=label_scale*13);

# Plot individual samples
text(mds$points[sample_idx,],labels=names[sample_idx], cex=label_scale*3, col=assigned_colors[sample_idx])

# Plot centroids
points(mds$points[centr_idx,],cex=2, col=colors, pch=19);
text(mds$points[centr_idx,],labels=names[centr_idx], cex=label_scale*5.5);

#print(mds$points);
#print(nrow(mds$points));
#print(length(names(A[1,])));
#print(length(assigned_colors));

# Plot the MDS expanded by group
for(group_idx in 1:num_unique_groups){

	cur_group=unique_groups[[group_idx]];
	cat("Working on group: ", cur_group, "\n", sep="");

	x_coord=numeric(0);
	y_coord=numeric(0);
	ingroup=numeric(0);
	outgroup=numeric(0);
	for(sample_idx in 1:num_samples){
		if(group[sample_idx] == cur_group){
			x_coord=c(x_coord, mds$points[sample_idx,1]);
			y_coord=c(y_coord, mds$points[sample_idx,2]);
			ingroup=c(ingroup,sample_idx);
		}else{
			outgroup=c(outgroup,sample_idx);
		}
	}

	min=min(x_coord);
	max=max(x_coord);
	width=(max-min);
	margin=width*(0.1);

	plot(mds$points,type="n", main=cur_group, xlim=c(min-margin, max+margin), ylim=c(min(y_coord), max(y_coord)), xlab="", ylab="")
#	text(mds$points,labels=names(A[1,]), cex=label_scale*3, col=assigned_colors)

	# Under X the outliers
	ingroup_outliers=intersect(outliers_in_matrix_idx, ingroup);
	outgroup_outliers=intersect(outliers_in_matrix_idx, outgroup);

	points(matrix(mds$points[ingroup_outliers,], ncol=2), col="grey30", pch=13, cex=label_scale*25, font=2);
	points(matrix(mds$points[outgroup_outliers,], ncol=2), col=assigned_colors[outgroup_outliers], pch=13, cex=label_scale*15);

	# Plot samples
	text(matrix(mds$points[outgroup,], ncol=2),labels=names[outgroup], cex=label_scale*3, col=assigned_colors[outgroup])
	text(matrix(mds$points[ingroup,], ncol=2),labels=names[ingroup], cex=label_scale*7, col="black")

	# Plot centroids
	points(matrix(mds$points[centr_idx,], ncol=2),cex=2, col=colors,pch=19);
	text(matrix(mds$points[centr_idx,], ncol=2),labels=names[centr_idx], cex=label_scale*5.5);

	
}


##################################################################################
# Plot dendrogram

cluster=hclust(dist);
dendrogram=as.dendrogram(cluster);

colorFun = function(node){
	if(is.leaf(node)){
		node_label=(attr(node, "label"));
		pos=regexpr(GroupRegex, node_label, perl=TRUE);
		group=substring(node_label, first=pos+1, last=pos+attr(pos,"match.length")-1);
		node_attr=attributes(node);
		attr(node, "nodePar") = c(node_attr$nodePar, list(lab.col=color_map[[group]], lab.cex=label_scale*1.5, cex=0.1));
	}
	node;
}

colored_dendrogram=dendrapply(dendrogram, colorFun);

plot(colored_dendrogram, main=InputFileName);

##################################################################################
# Write out outliers

outlier_filename=paste(OutputFileName, ".outliers", sep="");
fc=file(outlier_filename, "w");
for(group_id in 1:num_unique_groups){
	outlier_idx=outliers_by_group[[group_id]];
	#cat("outlier idx:", outlier_idx, "\n");
	for(i in 1:length(outlier_idx)){
	#	cat("i: ", i, " group: ", group_id, "\n");
	#	cat("Outlier index: ", outlier_idx[i], "\n");
	#	cat("centroid index: ", centr_idx[group_id], "\n");
	#	print(names[outlier_idx[i]])
	#	print(A[centr_idx[group_id], outlier_idx[i]]);
	#	cat("\n");
		outline=paste(unique_groups[group_id], ",", names[outlier_idx[i]], ",", A[ centr_idx[group_id], outlier_idx[i]], sep="");
		write(outline,file=fc);
	}
}











