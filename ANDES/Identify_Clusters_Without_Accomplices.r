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
		"distance_matrix", "d", 1, "character",
		"accomplices", "a", 1, "character",
		"num_target_clusters", "n", 1, "numeric",
		"clustering_method", "c", 2, "character",
		"accession_to_sample_name_map", "m", "2", "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage <- paste (
		"\nUsage:\n", script_name, "\n",
		"	-d <distance matrix>\n",
		"	-a <accomplice list>\n",
		"	-n <num clusters to identify>\n",
		"	[-c <clustering method, default wards>]\n",
		"	[-m <accession to sample name map>]\n",
		"\n",
		"This script will read in a distance matrix, cluster it\n",
		"using Ward's minimum variance, and then identify which\n",
		"clusters do not contain accomplices.\n",
		"\n");

if(!length(opt$distance_matrix) || !length(opt$accomplices)){
	cat(usage);
	quit(status=0);	
}

DistanceMatrix=opt$distance_matrix;
Accomplices=opt$accomplices;
NumTargetClusters=opt$num_target_clusters;
OutputFileName=gsub("\\.r_distmat", "", DistanceMatrix);

ClusteringMethod="ward";

if(length(opt$clustering_method)){
	ClusteringMethod=opt$clustering_method;
}

AccessionMapName="";
if(length(opt$accession_to_sample_name_map)){
	AccessionMapName=opt$accession_to_sample_name_map;
}

#------------------------------------------------------------------------------

cat("\n");
cat("Distance Matrix: ", DistanceMatrix, "\n", sep="");
cat("Accomplices: ", Accomplices, "\n", sep="");
cat("Num Target Clusters: ", NumTargetClusters, "\n", sep="");
cat("Clustering Method: ", ClusteringMethod, "\n", sep="");
cat("Output Filename Root: ", OutputFileName, "\n", sep="");
cat("Accession Map: ", AccessionMapName, "\n", sep="");
cat("\n");

###############################################################################

load_list=function(filename){
# Loads a list
	members=as.matrix(read.table(filename, sep="\t", header=FALSE));
	return(members[,1]);
}

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

library(MASS)

# Load the list of accomplices to avoid in clusters
accomplice_list=load_list(Accomplices);
num_accomplices=length(accomplice_list);
cat("Num accomplices read: ", num_accomplices, "\n");

# Load the accession to name mapping, if necessary
if(AccessionMapName!=""){
	accession_map=load_map(AccessionMapName);
	map_length=length(accession_map);
	cat("Accession map length: ", map_length, "\n");
}else{
	cat("No accession map file specified.\n");
}


###############################################################################

# Load data
distance_matrix=as.matrix(read.table(DistanceMatrix, header=TRUE, check.names=FALSE));

num_samples=nrow(distance_matrix);

# If accession map provide, rename the distance matrix now.
#if(AccessionMapName!=""){
#	rname=rownames(distance_matrix);
#	for(i in 1:num_samples){
#		if(length(accession_map[[rname[i]]])>0){
#			rname[i]=accession_map[[rname[i]]];
#		}else{
#			cat("Could not find name for: '", rname[i], "'\n", sep="");
#		}
#	}	
#	rownames(distance_matrix)=rname;
#}

if(AccessionMapName!=""){
	dup_name=list();
	rname=rownames(distance_matrix);
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
	rownames(distance_matrix)=rname;
}

# Get the sample names
sample_names=rownames(distance_matrix);

cat("Num samples in Distance Matrix:", num_samples, "\n");

# Compute the scale factor of the labels
label_scale=63/num_samples;
if(label_scale > 2){
	label_scale=2;
}

###############################################################################

# Define dendrogram apply function to colors notes according to group a sample belongs to
color_denfun_byGroup=function(n){
	if(is.leaf(n)){
		leaf_attr=attributes(n);
		leaf_name=leaf_attr$label;
		if(sum(samples_assigned_to_group==leaf_name)==1){
			group_color=group_to_color_map[[sample_to_group_map[[leaf_name]]]];
			group_id=group_to_group_id_map[[sample_to_group_map[[leaf_name]]]];

			attr(n, "nodePar") = c(leaf_attr$nodePar, 
						list(lab.col=group_color, group_id=group_id));
		}else{
			attr(n, "nodePar") = c(leaf_attr$nodePar, 
						list(lab.col="black", group_id=0));
		}
	}
	return(n);
}

color_denfun_byIndividual=function(n, sample_to_color_map){
	if(is.leaf(n)){
		leaf_attr=attributes(n);
		leaf_name=leaf_attr$label;
		ind_color=sample_to_color_map[[leaf_name]];
		if(is.null(ind_color)){
			ind_color="grey50";
		}

		attr(n, "nodePar") = c(leaf_attr$nodePar, 
						list(lab.col=ind_color));
	}
	return(n);
}

underline_sample_name_denfun_byIndividual=function(n, sample_list){
	if(is.leaf(n)){
		leaf_attr=attributes(n);
		#leaf_name=leaf_attr$label;
		leaf_name=strsplit(leaf_attr$label, ":")[[1]][1];	# Remove the accession
		if(any(leaf_name==sample_list)){
			att=attr(n, "nodePar");
			att$cex=.7;
			att$pch=18;
			att$col="red";
			att$lab.col="black";
			attr(n, "nodePar")=att;
		}
	}
	return(n);
}

get_color=function(den){
	if(is.leaf(den)){
		leaf_color=(attributes(den)$nodePar[["lab.col"]]);
		return(leaf_color);
	}else{
		num_nodes=length(den);
		colors=character();
		for(i in 1:num_nodes){
			colors[i]=get_color(den[[i]]);
		}
		unique_color=unique(colors);
		if(length(unique_color)==1){
			return(unique_color);
		}else{
			return("black");
		}
	}
}

get_group_id=function(den){
	if(is.leaf(den)){
		leaf_group_id=(attributes(den)$nodePar[["group_id"]]);
		return(leaf_group_id);
	}else{
		num_nodes=length(den);
		group_ids=numeric();
		for(i in 1:num_nodes){
			group_ids[i]=get_group_id(den[[i]]);
		}
		mean_group_id=mean(group_ids);
		return(mean_group_id);
	}

}

#-----------------------------------------------------------------------------------------------------

color_edges=function(den){
	consensus_color=get_color(den);
	edge_attr=attributes(den);
	attr(den, "edgePar") = c(edge_attr$edgePar, list(col=consensus_color));		
	return(den);
}

#-----------------------------------------------------------------------------------------------------

text_scale_denfun=function(n, label_scale){
	if(is.leaf(n)){
		leaf_attr=attributes(n);
		leaf_name=leaf_attr$label;
		#attr(n, "nodePar") = c(leaf_attr$nodePar, 
		#			cex=0,
		#			lab.cex=label_scale);
	        att=attr(n, "nodePar");
		att$lab.cex=label_scale;
		att$cex=0;
		attr(n, "nodePar")=att;
	}
	return(n);
}

######################################################################################################

allocate_colors=function(num_groups){
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

	return(rainbow_colors);
}

######################################################################################################

reorder_cluster_ids_by_size=function(clusters){
# Reassigns the cluster IDs by the cluster size

	tab=table(clusters);
	num_clusters=length(tab);
	ordered_tab=sort(tab, decreasing=TRUE);
	sorted_clusters=rep(0, length(clusters));
	for(i in 1:num_clusters){
		sorted_clusters[clusters==as.integer(names(ordered_tab[i]))]=i;
	}
	names(sorted_clusters)=names(clusters);
	return(sorted_clusters);
}

######################################################################################################

identify_acc_free_clusters=function(hclust, num_targets, accomplices){
# Keeps cutting until the number of clusters without members from the accomplice list equals the number of targets

	num_samples=length(hclust$labels);
	num_accomplice_free_clusters=0;

	# Cut until we find num_targets clusters without any accomplices in them
	k=0;
	while(k<num_samples && num_accomplice_free_clusters<num_targets){

		k=k+1;
		cat("Trying k =", k, "\n");
		clusters=cutree(hclust, k=k);
		clusters=reorder_cluster_ids_by_size(clusters);

		accomplices_found=numeric(k);
		for(i in 1:k){
			member_names=names(which(clusters==i));
			new_list =  c();
			for (c in 1:length(member_names)) {
				leaf_name=strsplit(member_names[c], ":")[[1]][1];	# Remove the accession
				new_list[c] = leaf_name;
			}
			if(length(intersect(new_list, accomplices))>0){
				accomplices_found[i]=1;
			}else{
				accomplices_found[i]=0;
			}
		}
		num_accomplice_free_clusters=k-sum(accomplices_found);
		cat("Num clusters without acc: ", num_accomplice_free_clusters, "\n");
	}

	# Extract out list for each accomplice free cluster
	free_clusters=list();
	cl_id=1;
	for(i in 1:k){
		if(accomplices_found[i]==0){
			free_clusters[[cl_id]]=names(which(clusters==i));
			cl_id=cl_id+1;
		}
	}

	free_clusters$k=k;
	free_clusters$free=num_accomplice_free_clusters;

	return(free_clusters);
}

######################################################################################################

find_height_at_k=function(hclust, k){
# Computes the height on the dendrogram for a particular k

	heights=hclust$height;
	num_heights=length(heights);
	num_clust=numeric(num_heights);
	for(i in 1:num_heights){
		num_clust[i]=length(unique(cutree(hclust, h=heights[i])));
	}
	height_idx=which(num_clust==k);
	midpoint=(heights[height_idx+1]+heights[height_idx])/2;
	return(midpoint);
}

######################################################################################################

# Computer cluster
hcl=hclust(as.dist(distance_matrix), method=ClusteringMethod);
hcl$height=sort(hcl$height);
dend=as.dendrogram(hcl);
max_clust_height=max(hcl$height);

# Rescale labels
dend=dendrapply(dend, text_scale_denfun, label_scale)

# Compute free clusters
free_clusters=identify_acc_free_clusters(hcl, NumTargetClusters, accomplice_list);

# Find height of k
cut_height=find_height_at_k(hcl, free_clusters$k);

# Allocate colors
colors=allocate_colors(free_clusters$free);

# Assign colors to samples
sample_to_color_map=list();
cat("Num free clusters: ", free_clusters$free, "\n", sep="");
for(i in 1:free_clusters$free){
	members=free_clusters[[i]];
	for(j in 1:length(members)){
		#cat(members[j], " -> ", i , "\n");
		sample_to_color_map[[members[j]]]=colors[i];	
	}
}

# Color and decorate dendrogram
dend=dendrapply(dend, color_denfun_byIndividual, sample_to_color_map);
dend=dendrapply(dend, color_edges);
dend=dendrapply(dend, underline_sample_name_denfun_byIndividual, accomplice_list);

######################################################################################################
# Output dendrograms in pdf

dendro_pdf=paste(OutputFileName, ".accomplice_free_clusters.pdf", sep="")
pdf(dendro_pdf, height=8.5, width=11);

# Plot overall dendrogram
plot(dend, main=OutputFileName);
abline(h=cut_height, col="red", lty=2, lwd=.75);
text(0, cut_height, label=sprintf("k=%g", free_clusters$k), pos=3, cex=.8);
legend(num_samples*.75, max_clust_height, legend=1:free_clusters$free, fill=colors);

# Plot zoomed in
par(mfrow=c(2,1));
par(mar=c(10,0,1,0))
par(oma=c(3,0,1,0))
dend=dendrapply(dend, text_scale_denfun, label_scale*2);
dend=dendrapply(dend, underline_sample_name_denfun_byIndividual, accomplice_list);
plot(dend, ylim=c(0,attributes(dend)$height), xlim=c(0,(num_samples/2)*1.05), main=OutputFileName);
abline(h=cut_height, col="red", lty=2, lwd=.75);
plot(dend, ylim=c(0,attributes(dend)$height), xlim=c((num_samples/2)/1.05,num_samples));
abline(h=cut_height, col="red", lty=2, lwd=.75);

######################################################################################################
# Output cluster members list

members_file=paste(OutputFileName, ".accomplice_free_clusters.tsv", sep="");
fh=file(members_file, "w");

for(i in 1:free_clusters$free){
	cluster_members=free_clusters[[i]];
	num_members=length(cluster_members);
	for(j in 1:num_members){
		cat(file=fh, paste(i, cluster_members[j], sep="\t"), "\n", sep="");
	}
}

close(fh);

## output cluster count k
clust_file=file("cluster_count", "w");
cat(file=clust_file, free_clusters$k, "\n", sep="");
close(clust_file);
######################################################################################################

dev.off()
cat("Done.\n");
