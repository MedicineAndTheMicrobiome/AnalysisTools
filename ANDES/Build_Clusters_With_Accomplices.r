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
		"accomplice_groups", "a", 1, "character",
		"output_file", "o", 2, "character",
		"accession_to_sample_name_map", "m", 2, "character",
		"special_list", "s", 2, "character",
		"target_k", "k", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage <- paste (
		"\nUsage:\n", script_name, "\n",
		"	-d <distance matrix>\n",
		"	-a <accomplice groups>\n",
		"	[-o <output file name>]\n",
		"	[-m <accession to sample name map>]\n",
		"	[-k <target number of clusters, k>]\n",
		"	[-s <special sample list>]\n",
		"\n",
		"This script will read in a distance matrix (eg. from sequence comparisons)\n",
		"build a Ward's minimum variance clustering, then utilize the accomplice groups\n",
		"to identify which sequence groups may be similar to each other based on the\n",
		"accomplices position to each other.\n",
		"\n",
		"The -k option is only for annotation.  It does not affect the calculation.\n",
		"The -s option is for annotating special samples.\n",
		"\n");


## Parameters Parsing and Variable Assignment
if(!length(opt$distance_matrix) || !length(opt$accomplice_groups)){
	cat(usage);
	quit(status=0);	
}

DistanceMatrix=opt$distance_matrix;
AccompliceGroups=opt$accomplice_groups;

OutputFileName=gsub("\\.r_distmat", "", DistanceMatrix);
if(length(opt$output_file)){
	OutputFileName=opt$output_file;
}

TargetK=-1;
if(length(opt$target_k>0)){
	TargetK=opt$target_k;
}

AccessionMapName="";
if(length(opt$accession_to_sample_name_map)){
	AccessionMapName=opt$accession_to_sample_name_map;
}

SpecialSampleList="";
if(length(opt$special_list)){
	SpecialSampleList=opt$special_list;
}

#------------------------------------------------------------------------------

cat("\n");
cat("Distance Matrix: ", DistanceMatrix, "\n", sep="");
cat("Accomplice Groups: ", AccompliceGroups, "\n", sep="");
cat("Output Filename Root: ", OutputFileName, "\n", sep="");
cat("Accession Map: ", AccessionMapName, "\n", sep="");
cat("Target K: ", TargetK, "\n", sep="");
cat("Special Sample List: ", SpecialSampleList, "\n", sep="");
cat("\n");

## Functions
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

load_groups=function(filename){
# Loads a two column matrix so that the first column can be translated into the second column
	members=as.matrix(read.table(filename, sep="\t", header=FALSE));
	num_members=nrow(members);
	cat("Num members in map: ", num_members, "\n");
	map=list();
	for(i in 1:num_members){
		map[[members[i,2]]]=members[i,1];	
	}
	return(map);
}

# Define dendrogram apply function to colors notes according to group a sample belongs to
color_denfun_byGroup=function(n, sample_to_group_map){
	if(is.leaf(n)){
		leaf_attr=attributes(n);
		leaf_name=strsplit(leaf_attr$label, ":")[[1]][1];	# Remove the accession

		if(length(sample_to_group_map[[leaf_name]])){
			att=attr(n, "nodePar");
			att$lab.col=sample_to_group_map[[leaf_name]];
			att$cex=.5;
			att$pch=18;
			att$col="red";
			att$lab.font=4; #2=bold, 3=italic 4=bold_italic
			att$lab.cex=att$lab.cex*1.2
			attr(n, "nodePar")=att;
		}
	}
	return(n);
}

get_color=function(den){
	if(is.leaf(den)){
		leaf_color=(attributes(den)$nodePar[["lab.col"]]);
		if(length(leaf_color)==0){
			return("grey");
		}else{
			return(leaf_color);
		}
	}else{
		num_nodes=length(den);
		colors=character();
		for(i in 1:num_nodes){
			colors[i]=get_color(den[[i]]);
		}
		unique_color=unique(setdiff(colors, "grey"));
		if(length(unique_color)==0){
			return("grey");
		}else if(length(unique_color)==1){
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

color_all=function(den, color){
	if(is.leaf(den)){
		att=attr(den, "nodePar");
		att$lab.col=color;
		attr(den, "nodePar")=att;
	}
	attributes(den)$edgePar[["col"]]=color;
	return(den);
}

color_edges=function(den){
	consensus_color=get_color(den);
	edge_attr=attributes(den);
	attr(den, "edgePar") = c(edge_attr$edgePar, list(col=consensus_color));		
	return(den);
}

text_scale_denfun=function(n, label_scale){
	if(is.leaf(n)){
		leaf_attr=attributes(n);
		leaf_name=leaf_attr$label;
		    att=attr(n, "nodePar");
		    att$lab.cex=label_scale;
		    att$cex=0;
		attr(n, "nodePar")=att;
	}
	return(n);
}

color_children_of_homogeneous_parents=function(n){
	if(is.leaf(n)){
		# Do nothing
		return(n);
	}else{
		num_children=length(n);

		# If a edge, if both children are the same color, set their children to the same color.
		left_attr=attr(n[[1]], "edgePar");
		right_attr=attr(n[[2]], "edgePar");
	
		if((left_attr$col==right_attr$col) &&
			(left_attr$col!="grey") &&
			(left_attr$col!="black")){

			n[[1]]=dendrapply(n[[1]], color_all, left_attr$col);
			n[[2]]=dendrapply(n[[2]], color_all, left_attr$col);
		}else{
			n[[1]]=color_children_of_homogeneous_parents(n[[1]]);
			n[[2]]=color_children_of_homogeneous_parents(n[[2]]);
		}

		return(n);
	}
}

get_color_of_all_leafs=function(n){
	mapping=list();
	save_name_to_color_mapping=function(den, mapping){
		if(is.leaf(den)){
			attr=attributes(den);
			#cat(attr$label, " ", attr$nodePar$lab.col, "\n");
			if(length(attr$nodePar$lab.col)==0){
				attr$nodePar$lab.col=-1;
			}
			mapping[[attr$label]]=as.numeric(attr$nodePar$lab.col);
		}else{
			num_children=length(den)
			for(i in 1:num_children){
				mapping=save_name_to_color_mapping(den[[i]], mapping);
			}
		}
		return(mapping);
	}
	mapping=save_name_to_color_mapping(n, mapping);
	return(mapping);
}

allocate_colors=function(num_groups){
	if (num_groups == 2 || num_groups == 1) {
		colors=c("#0e8406", "#8a231f");
		rainbow_colors=colors[1:2];
	} else {
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
	}
	return(rainbow_colors);
}

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

# Load the accession to name mapping, if necessary
if(AccessionMapName!=""){
	accession_map=load_map(AccessionMapName);
	map_length=length(accession_map);
	cat("Accession map length: ", map_length, "\n");
}else{
	cat("No accession map file specified.\n");
}

accomplice_groups=load_groups(AccompliceGroups);
accomplice_list=names(accomplice_groups);
accomplice_group_list=sort(unique(paste(accomplice_groups)));
num_accomplice_groups=length(accomplice_group_list);

cat("Accomplice Groups: \n");
print(accomplice_group_list);

# Map group names to color by just counting off from 1.
accomplice_group_color_map=list();
for(i in 1:length(accomplice_list)){
	key=accomplice_list[i]
	accomplice_group_color_map[[key]]=which(accomplice_groups[[key]]==accomplice_group_list);
}

###############################################################################

# Load data
distance_matrix=as.matrix(read.table(DistanceMatrix, header=TRUE, check.names=FALSE));

num_samples=nrow(distance_matrix);

# If accession map provide, rename the distance matrix now.
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

# Load special sample list
if(SpecialSampleList!=""){
	special_samples=load_list(SpecialSampleList);
}else{
	special_samples=c();
}

# Get the sample names
sample_names=rownames(distance_matrix);

cat("Num samples in Distance Matrix:", num_samples, "\n");

# Compute the scale factor of the labels
label_scale=50/num_samples;
if(label_scale > 2){
	label_scale=2;
}

# Compute cluster
hcl=hclust(as.dist(distance_matrix), method="ward");
sample_names_by_position=hcl$labels[hcl$order];
dend=as.dendrogram(hcl);
max_clust_height=max(hcl$height);

# Find height at K
if(TargetK>0){
	cut_height=find_height_at_k(hcl, TargetK);	
}

# Rescale labels
dend=dendrapply(dend, text_scale_denfun, label_scale)
dend=dendrapply(dend, color_denfun_byGroup, accomplice_group_color_map);
dend=dendrapply(dend, color_edges);

dend=color_children_of_homogeneous_parents(dend);

# Find position of special samples
special_sample_pos=which(sample_names_by_position==special_samples);

# Extract interpolated clusters
interpolated_leaf_colors=get_color_of_all_leafs(dend);
counts=numeric(num_accomplice_groups+1);
names(counts)=c(accomplice_group_list, "Unknown");
j=1;
for(i in c(1:num_accomplice_groups, -1)){
	counts[j]=sum(interpolated_leaf_colors==i);
	j=j+1;
}
group_perc=counts/sum(counts);
print(counts);
if(sum(counts)!=num_samples){
	cat("Error:  Sum of interpolated clusters, does not match total input samples.\n");
}

cat("Resampling...\n");

re_num = 40;
resampcounts=matrix(nrow=re_num,ncol=num_accomplice_groups+1);
colnames(resampcounts)=c(accomplice_group_list, "Unknown");
boot=file(paste(OutputFileName, ".boot", sep=""), "w");
#len=length(counts);
#cats=c(accomplice_group_list, "Unknown");

for (s in 1:re_num) {
	cat(".");
	rdend=as.dendrogram(hcl);
	
	re_acc_group_map = resample(accomplice_group_color_map, length(accomplice_group_color_map), replace=TRUE)
	
	# Rescale labels
	rdend=dendrapply(rdend, text_scale_denfun, label_scale)
	rdend=dendrapply(rdend, color_denfun_byGroup, re_acc_group_map);
	rdend=dendrapply(rdend, color_edges);
	
	rdend=color_children_of_homogeneous_parents(rdend);
	
# Find position of special samples
	special_sample_pos=which(sample_names_by_position==special_samples);
	
# Extract interpolated clusters
	interpolated_leaf_colors=get_color_of_all_leafs(rdend);

	j=1;
	for(i in c(1:num_accomplice_groups, -1)){
		#if (j == 1) {
			cat(file=boot, sprintf("%s\t%02i\n", accomplice_group_list[j], sum(interpolated_leaf_colors==i)));
		#}
		resampcounts[s,j]=sum(interpolated_leaf_colors==i);
		j=j+1;
	}
}
cat("done.\n");
print(resampcounts);

med=numeric(num_accomplice_groups+1);
lower=numeric(num_accomplice_groups+1);
upper=numeric(num_accomplice_groups+1);
j=1;
alpha=.05;
for(i in c(1:num_accomplice_groups, -1)){
	vals = ci(resampcounts[,j],alpha);
	med[j] = vals[1];
	lower[j] = vals[2];
	upper[j] = vals[3];
	j=j+1
}
######################################################################################################
# Output dendrograms in pdf

dend_height=attributes(dend)$height;

dendro_pdf=paste(OutputFileName, ".clusters_wAccomplices.pdf", sep="")
pdf(dendro_pdf, height=8.5, width=11);

# Plot overall dendrogram
palette(allocate_colors(num_accomplice_groups));
plot(dend, main=OutputFileName, ylim=c(dend_height*-.2, dend_height));
par(mar=c(6,1,1,1));
legend(
	x=6*num_samples/10,
	y=dend_height,
	legend=sprintf("%s: %i (%3.1f%%) median %3.1f%% CI: %3.1f%% - %3.1f%%", c(accomplice_group_list, "Unknown"), counts, group_perc*100, med/sum(counts)*100, lower/sum(counts)*100, upper/sum(counts)*100),
	fill=c(1:num_accomplice_groups, "grey"),
	cex=.8
);

## make coverage file
covfile=file(paste(OutputFileName, ".cov", sep=""), "w");
len=length(counts);
cats=c(accomplice_group_list, "Unknown");
for (i in 1:len) {
	cat(file=covfile, sprintf("%s\t%02i\n", cats[i], counts[i]));
}

# Highlight samples of interest
if(length(special_sample_pos)>0){
	segments(
		x0=special_sample_pos, x1=special_sample_pos,
		y0=-10, y1=0,
		col="#FFFF0048", lwd=3);
}

# If target K specified, draw line
if(TargetK>0){
	abline(h=cut_height, col="red", lty=2, lwd=.75);
	text(0, cut_height, sprintf("k=%g",TargetK), cex=.75, pos=3);
}

# Plot zoomed in
par(mfrow=c(2,1));
par(mar=c(5,0,3,0))
par(oma=c(1,0,1,0))
dend=dendrapply(dend, text_scale_denfun, label_scale*2);
dend=dendrapply(dend, color_denfun_byGroup, accomplice_group_color_map);

# Top
plot(dend, ylim=c(dend_height*-5,dend_height), xlim=c(0,(num_samples/2)*1.05), main=OutputFileName);
if(TargetK>0){
	abline(h=cut_height, col="red", lty=2, lwd=.75);
}
if(length(special_sample_pos)>0){
	segments(
		x0=special_sample_pos, x1=special_sample_pos,
		y0=-10, y1=0,
		col="#FFFF0048", lwd=3);
}

# Bottom
plot(dend, ylim=c(dend_height*-5,dend_height), xlim=c((num_samples/2)/1.05,num_samples));
if(TargetK>0){
	abline(h=cut_height, col="red", lty=2, lwd=.75);
}
if(length(special_sample_pos)>0){
	segments(
		x0=special_sample_pos, x1=special_sample_pos,
		y0=-10, y1=0,
		col="#FFFF0048", lwd=3);
}


## stats graphs
perc=c(0,.25,.50,.75,1);

par(mfrow=c(2,2),oma=c(2,2,0,2));
par(mar=c(4.8, 4, 7.5, 2.1));
boxplot(resampcounts, main="Spread of Leaves Covered During Resampling by Group", yaxt="n", ylab="Coverage", xlab="Accomplice Group");
axis(side=4,at=(num_samples*perc),labels=sprintf("%3.1f%%", perc*100),tick=TRUE, cex=.5);
axis(side=2,at=(num_samples*perc),labels=sprintf("%2.0f", num_samples*perc), tick=TRUE, cex=.5);
j=1;
for(i in c(1:num_accomplice_groups, -1)){
	abline(h=counts[j], col="red");
	text(j,counts[j], sprintf("%3.1f%%", counts[j]/num_samples*100), col="red", adj = c(0, -.1))
	j=j+1;
}
j=1;
for(i in c(1:num_accomplice_groups, -1)){
	hist(resampcounts[,j], main=paste("Frequency of Coverage During Resampling for", colnames(resampcounts)[j]), ylab="Frequency", xlab=colnames(resampcounts)[j], xlim=c(1,num_samples), xaxt="n", breaks=c(seq(0, num_samples+10, 10)), ylim=c(0,re_num));
	axis(side=3,at=num_samples*perc,labels=sprintf("%3.1f", perc*100), tick=TRUE, cex=.5, las=2);
	axis(side=3,at=counts[j],labels=sprintf("%3.1f", counts[j]/num_samples*100), tick=TRUE, cex=.5, las=2);
	axis(side=1,at=num_samples*perc,labels=sprintf("%2.0f", num_samples*perc), tick=TRUE, cex=.5, las=2);
	axis(side=1,at=counts[j],labels=sprintf("%2.0f", counts[j]), tick=TRUE, cex=.5, las=2);
	abline(v=counts[j], col="red");
	j=j+1;
}
######################################################################################################

dev.off()
cat("Done.\n");
