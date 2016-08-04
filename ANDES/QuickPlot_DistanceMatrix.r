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
		"input_file", "i", 1, "character",
		"groups_file", "l", 2, "character",
		"color_map_file", "c", 2, "character",
		"legend", "L", 2, "logical",
		"minflip", "f", 2, "logical",
		"description_file", "d", 2, "character",
		"paper_dimension", "p", 2, "character"
);

PAPERDIM="8.5x11";

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage <- paste (
		"\nUsage:\n", script_name, "\n",
		"	-i <input file>\n",
		"	[-l <group labels file>]\n",
		"	[-L (plot a legend for the labels)]\n",
		"	[-c <individual colormap file>]\n",
		"	[-f (turn on branch flipping to push lower groups to left)]\n",
		"	[-d <description map>]\n",
		"	[-p <paper dimensions, default = ", PAPERDIM, ">]\n",
		"\n",
		"Plots a distance matrix using R's isoMDS and dendrogram.\n",
		"This serves as a quick way to visualize your distance matrix.\n",
		"It's not publication quality by any means...\n",
		"\n",
		"The format of the group labels file is 2-column tab-deliminated file:\n",
		"	<group id>\\t<member id>\\n\n",
		"\n",
		"The format of the individual color map file is:\n",
		"	<sample name>\\t<#RGB color>\\n\n",
		"\n",
		"Note that the descript map will be applied after coloring, so the colormap\n",
		"must be associated with the accessions in the distance matrix, not the descriptions\n",
		"in the description map.\n",
		"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;

GroupFilename=character();
GroupProvided=FALSE;
if(length(opt$groups_file)){
	GroupFilename=opt$groups_file;
	GroupProvided=TRUE;
}

ColorMapFilename=character();
ColorMapProvided=FALSE;
if(length(opt$color_map_file)){
	ColorMapFilename=opt$color_map_file;
	ColorMapProvided=TRUE;
}


if(ColorMapProvided && GroupProvided){
	cat("Error:  You can't provide both a group coloring and individual sample coloring.\n");
	q(status=-1);
}

PlotLegend=FALSE;
if(length(opt$legend)){
	PlotLegend=TRUE;
}

MinFlip=FALSE;
if(length(opt$minflip)){
	MinFlip=TRUE;
}

DescriptionFile="";
if(length(opt$description_file)){
	DescriptionFile=opt$description_file;
}

PaperDim=PAPERDIM;
if(length(opt$paper_dimension)){
	PaperDim=opt$paper_dimension;
}

hw=as.numeric(strsplit(PaperDim, "x")[[1]]);
cat("Paper Dimensions: ", hw[1], " by ", hw[2], "\n");


###############################################################################

library(MASS)

# Get input file name, so we can automatically set an output file name
cat("Input File: ", InputFileName, "\n");
if(GroupProvided){
	cat("Label File: ", GroupFilename, "\n");
}
MdsPDF=paste(InputFileName, ".mds.pdf", sep="")
DendrogramPDF=paste(InputFileName, ".dendrogram.pdf", sep="")

# Load data
A<-as.matrix(read.table(InputFileName, header=TRUE, check.names=FALSE));

# Remove 0 distances, by setting them to a very very small number instead
noZeroA=A;
for(x in 1:nrow(A)){
	for(y in 1:ncol(A)){
		if(A[x,y]==0){
			noZeroA[x,y]=0.0000000000001;
		}
	}
}

# Compute label scales based on number of samples
numSamples=nrow(A);
label_scale=50/numSamples;
if(label_scale > 2){
	label_scale=2;
}

# Assign colors to groups by sample id
if (GroupProvided) {
	
	# Load cluster assignments
	cat("Loading: ", GroupFilename, "\n");
	cluster_assignments=as.matrix(read.table(GroupFilename, header=FALSE, check.names=FALSE, sep="\t"));

	# First column is group, second column is sample id
	num_samples_in_group=length(cluster_assignments[,2])
	samples_assigned_to_group=cluster_assignments[,2];
	group_names=unique(cluster_assignments[,1])
	num_groups=length(group_names);

	# Assign members to group
	sample_to_group_map=list()
	for(i in 1:num_samples_in_group) {
		sample_to_group_map[cluster_assignments[i,2]]=cluster_assignments[i,1]
	}

	# Allocate colors
	rainbow_natural=rainbow(num_groups)
	rainbow_lowval=rainbow(num_groups,v=.5)
	rainbow_lowsat=rainbow(num_groups,s=.75)
	odd=seq(1,num_groups,2);
	even=seq(2,max(num_groups,2),2);
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

	# Assign group number to group
	group_to_group_id_map=list();
	for(i in 1:num_groups) {
		group_to_group_id_map[[group_names[i]]]=i;
	}

}

if(ColorMapProvided){

	# Load colormap
	cat("Loading: ", ColorMapFilename, "\n");
	color_assignments=as.matrix(read.table(ColorMapFilename, header=FALSE, check.names=FALSE, comment.char="", sep="\t"));
	num_color_assignments=nrow(color_assignments);

	# Convert matrix into hash
	sample_to_color_map=list();
	for(i in 1:num_color_assignments){
		sample_to_color_map[[color_assignments[i,1]]]=color_assignments[i,2];
	}

	#print(sample_to_color_map);
}


if(DescriptionFile!=""){
	
	# Load accession to description mapping file
	cat("Loading: ", DescriptionFile, "\n");

	description_mat=as.matrix(read.table(DescriptionFile, header=FALSE, check.names=FALSE, comment.char="", sep="\t"));
	num_descriptions=nrow(description_mat);

	# Convert matrix into hash
	sample_to_description_map=list();
	for(i in 1:num_descriptions){
		sample_to_description_map[[description_mat[i,1]]]=description_mat[i,2];
	}
		
	cat("Num descriptions loaded: ", num_descriptions, "\n");
	#print(sample_to_description_map);
}

######################################################################################################

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

color_denfun_byIndividual=function(n){
	if(is.leaf(n)){
		leaf_attr=attributes(n);
		leaf_name=leaf_attr$label;
		ind_color=sample_to_color_map[[leaf_name]];
		if(is.null(ind_color)){
			ind_color="black";
		}

		attr(n, "nodePar") = c(leaf_attr$nodePar, 
						list(lab.col=ind_color));
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

assign_edge_group_id=function(den){
	mean_group_id=get_group_id(den);
	edge_attr=attributes(den);
	attr(den, "edgePar") = c(edge_attr$edgePar, list(group_id=mean_group_id));		
	return(den);
}

#-----------------------------------------------------------------------------------------------------

text_scale_denfun=function(n){
	if(is.leaf(n)){
		leaf_attr=attributes(n);
		leaf_name=leaf_attr$label;
		attr(n, "nodePar") = c(leaf_attr$nodePar, 
					cex=0,
					lab.cex=label_scale);
	}
	return(n);
}

#---------------------------------------------------------------------------------------------------------

minflip_den=function(den){

	if(!is.leaf(den)){
		
		num_edges=length(den);
		for(i in 1:num_edges){
			den[[i]]=minflip_den(den[[i]]);
		}

		temp_edges=list();
		group_ids=numeric();
		for(i in 1:num_edges){
			group_ids[i]=attributes(den[[i]])$edgePar$group_id;
			temp_edges[[i]]=den[[i]];
		}
		sort_info=sort(group_ids, index.return=TRUE, decreasing=FALSE);
		for(i in 1:num_edges){
			den[[i]]=temp_edges[[sort_info$ix[i]]];
		}
	}
	return(den)
}

#---------------------------------------------------------------------------------------------------------

rename_labels_with_descriptions=function(den){
	if(is.leaf(den)){
		accession=attr(den, "label");
		description=sample_to_description_map[[accession]];
		if(length(description)==0 || is.null(description)){
			description=accession;
		}
		attr(den, "label")=description;
	}
	return(den);
}


#---------------------------------------------------------------------------------------------------------

calc_entropy=function(x){
	xnorm=x/sum(x);
	return(-sum(xnorm*log(xnorm)));
}

#---------------------------------------------------------------------------------------------------------

plot_cluster_curve=function(hclust_out, type){
	height_min=min(hclust_out$height);
	height_max=max(hclust_out$height);
	precision=(height_max-height_min)/200;

	increments=seq(from=height_max, to=height_min, by=-precision);
	num_increments=length(increments);

	entropy=numeric(num_increments);
	evenness=numeric(num_increments);
	singletons=numeric(num_increments);
	num_clusters=numeric(num_increments);

	# hack for known hclust bizarreness
	hclust_out$height=sort(hclust_out$height);
	
	i=1;
	for(height in increments){
		clusters=cutree(hclust_out, h=height);
		frequencies=table(clusters);
		num_clusters[i]=length(frequencies);
		entropy[i]=calc_entropy(frequencies);
		evenness[i]=entropy[i]/log(num_clusters[i]);
		singletons[i]=sum(frequencies==1);		
		i=i+1;
	}

	par(mfrow=c(2,2));
	rev_increments=max(increments)-increments;
	label_ix=floor(seq(1,length(increments),length.out=5));

	plot(rev_increments,num_clusters,xaxt="n", xlab="Height", ylab="Num Clusters", main=type);
	axis(side=1, at=rev_increments[label_ix], labels=sprintf("%3.2f",rev(rev_increments[label_ix])));
	plot(rev_increments,entropy,xaxt="n", xlab="Height", ylab="Entropy", main=type);
	axis(side=1, at=rev_increments[label_ix], labels=sprintf("%3.2f",rev(rev_increments[label_ix])));
	plot(rev_increments,evenness,xaxt="n", xlab="Height", ylab="Evenness", main=type);
	axis(side=1, at=rev_increments[label_ix], labels=sprintf("%3.2f",rev(rev_increments[label_ix])));
	plot(rev_increments,singletons,xaxt="n", xlab="Height", ylab="Num Singletons", main=type);
	axis(side=1, at=rev_increments[label_ix], labels=sprintf("%3.2f",rev(rev_increments[label_ix])));
	par(mfrow=c(1,1));
}

######################################################################################################

AasDist=as.dist(A);
sorteddist=sort(AasDist);
cat("10 Most Distant\n");
tail(sorteddist,10);
cat("10 Most Similar\n");
head(sorteddist,10);

######################################################################################################

pdf(DendrogramPDF, height=hw[1], width=hw[2]);

clust_methods=c("ward","single","average","complete");
clust_descr=c("Ward's Minimum Variance", "Single Linkage (Nearest Neighbor)", "Average Linkage", "Complete Linkage (Farthest Neighbor)");

#for(method_idx in 1:1){
for(method_idx in 1:length(clust_methods)){

	cat("Working on: ", clust_descr[method_idx], "\n", sep="");
	# Compute dendrogram
	hclust_out=hclust(AasDist, method=clust_methods[method_idx]);
	dendro=as.dendrogram(hclust_out);

	# Apply colors if group specified
	if(GroupProvided){
		dendro=dendrapply(dendro, color_denfun_byGroup);
		dendro=dendrapply(dendro, color_edges);

		dendro=dendrapply(dendro, assign_edge_group_id);

		if(MinFlip){
			dendro=minflip_den(dendro);
		}

	}
	if(ColorMapProvided){
		dendro=dendrapply(dendro, color_denfun_byIndividual);
	}

	# Change accessions to description
	if(DescriptionFile!=""){
		dendro=dendrapply(dendro, rename_labels_with_descriptions);
	}

	dendro=dendrapply(dendro, text_scale_denfun);


	# Color edges
	# Plot dendrogram
	par(mar=c(10,4,4,4))
	plot(dendro, main=InputFileName, ylab=clust_descr[method_idx]);

	if(PlotLegend){
	    	legend(numSamples*.75, attributes(dendro)$height, fill=rainbow_colors, legend=group_names);
	}

	# Compute structural statistics
	par(mar=c(3,3,3,3))
	plot_cluster_curve(hclust_out, clust_descr[method_idx]);

	cat("\tcompleted.\n");

}

######################################################################################################

# Compute MDS
mds<-isoMDS(noZeroA)

# Compute the padding inside of the graph so the text labels don't fall off the screen
min<-min(mds$points[,1])
max<-max(mds$points[,1])
width=(max-min)
margin=width*(0.1)

# Assign colors if group is specified
mds_colors=rep("black",numSamples);
if(GroupProvided){
	sample_names=colnames(A);
	for(i in 1:numSamples){
		if(sum(samples_assigned_to_group==sample_names[i])==1){
			mds_colors[i]=group_to_color_map[[sample_to_group_map[[sample_names[i]]]]];
		}else{
			mds_colors[i]="black";
		}
	}
}
if(ColorMapProvided){
	sample_names=colnames(A);
	for(i in 1:numSamples){
		ind_color=sample_to_color_map[[sample_names[i]]];
		if(!is.null(ind_color)){
			mds_colors[i]=ind_color;	
		}
	}
}

pdf(MdsPDF,width=11,height=8.5)

# Plot with labels
plot(mds$points,type="n", main=InputFileName, xlim=c(min-margin, max+margin),xlab="",ylab="")

mdslabels=names(A[1,]);
if(DescriptionFile!=""){
	num_labels=length(mdslabels);
	for(i in 1:num_labels){
		desc=sample_to_description_map[[mdslabels[i]]];
		if(!is.null(desc)){
			mdslabels[i]=desc;
		}
	}
}

text(mds$points,labels=mdslabels, cex=label_scale, col=mds_colors)

# Plot with symbols
plot(mds$points, main=InputFileName, xlim=c(min-margin, max+margin),xlab="",ylab="", col=mds_colors)

######################################################################################################

dev.off()
cat("Done.\n");
