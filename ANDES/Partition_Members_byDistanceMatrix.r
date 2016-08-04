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
	"cluster_method", "c", 2, "character",
	"height_cutoff", "h", 2, "numeric",
	"num_clusters", "k", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n    ", script_name, "\n",
	"	-d <input distance matrix>\n",
	"	[-c <cluster method name, default 'complete'>]\n",
	"	[-h <cluster height cutoff, default automatic>]\n",
	"	[-k <num clusters, default automatic>]\n",
	"\n",
	"Reads in a distance matrix and generates a list of ids for each partition\n",
	"that belong to the same cluster at the specified or automatic cutoff.\n",
	"\n",
	"The possible values for cluster method are:\n",
	"\tward: Ward's minimum variance\n",
	"\tsingle: Single linkage / Nearest Neighbor\n",
	"\tcomplete: Complete linkage / Farthest Neightbor\n",
	"\taverage: Average linkage\n",
	"\n",
	"You may also use mcquitty, median or centroid.\n",
	"\n");

if(!length(opt$distance_matrix)){
	cat(usage);
	q(status=-1);
}

InputDistanceMatrix=opt$distance_matrix;

ClusterMethod="complete";
if(length(opt$cluster_method)){
	ClusterMethod=opt$cluster_method;
}

HeightCutoff=0;
if(length(opt$height_cutoff)){
	HeightCutoff=opt$height_cutoff;
}

KCutoff=0;
if(length(opt$num_clusters)){
	KCutoff=opt$num_clusters;
}

###############################################################################

cat("Input Distance Matrix: ", InputDistanceMatrix, "\n");
cat("Clustering Method: ", ClusterMethod, "\n");

if(HeightCutoff!=0){
	cat("Height Cutoff: ", HeightCutoff, "\n");
}

if(KCutoff!=0){
	cat("Num Cluster Defined: ", KCutoff, "\n");
}

###############################################################################

plot_colored_dendrogram=function(dendro_in, names, clusters, colors=NULL){

        num_clusters=length(unique(clusters));
        # Allocate mapping variables
        name_to_cluster_mapping=list();

        # Get all the lengths
        names_length=length(names);
        clusters_length=length(clusters);
        num_clusters=length(unique(clusters));
        if(is.null(colors)){
                colors=rainbow(num_clusters);
        }
        num_colors=length(colors);

        # Confirm name lengths match number of cluster assignments
        if(names_length!=clusters_length){
                cat("Error: names length don't match cluster assignment lengths.\n");
                return(NULL);
                #q(status=-1);
        }

        # Create cluster mapping/hash
        for(i in 1:names_length){
                name_to_cluster_mapping[[names[i]]]=clusters[i];
        }

        # Assign colors to clusters;
        if(num_colors<num_clusters){
                cat("Not enough colors.");
                return(NULL);
                #q(status=-1);
        }

        # Resize labels based on number of samples
        label_scale=63/names_length;
        if(label_scale>2){
                label_scale=2;
        }

        # Insert color attribute into node
        color_denfun=function(n){
                if(is.leaf(n)){
                        leaf_attr=attributes(n);
                        leaf_name=leaf_attr$label;

                        color_assignment=colors[name_to_cluster_mapping[[leaf_name]]];

                        attr(n, "nodePar") = c(leaf_attr$nodePar,
                                                list(
                                                        lab.col=color_assignment,
                                                        lab.cex=label_scale,
                                                        cex=0
                                                ));
                }
                return(n);
        }

        # Execute coloring function on dendrogram
        dendro_out=dendrapply(dendro_in, color_denfun);

        # Plot
        plot(dendro_out, main=paste("k = ", num_clusters, sep=""));
}

###############################################################################
# Plots a tree with a line marking where the tree was cut, calls code to color nodes

plot_marked_tree=function(hclust, sample_names, cluster, colors, height=NULL){

        dendro=as.dendrogram(hclust);
        plot_colored_dendrogram(dendro, sample_names, cluster, colors);
        dendro_attrib=attributes(dendro);
        k=length(unique(cluster));

	if(height==0){
		height_a=0;
		height_b=0;
		hclust$height=sort(hclust$height);
		for(i in seq(0, dendro_attrib$height, length.out=1000)){
			clusters=cutree(hclust, h=i);
			num_clusters=length(unique(clusters));
			if(num_clusters>(k-1)){
				height_a=i;
			}
			if(num_clusters>k){
				height_b=i;
			}
			if(num_clusters<k){
				break;
			}
		}

		abline(h=(height_a+height_b)/2, col="red", lty=4);
	}else{
		abline(h=height, col="red", lty=4);
	}

}


###############################################################################

library(MASS);

# Get input file name, so we can automatically set an output file name
OutputFileNamePDF=paste(InputDistanceMatrix, ".cutmarked.pdf", sep="");

# Load data
cat("Loading matrix.\n");
A=as.matrix(read.table(InputDistanceMatrix, header=TRUE, check.names=FALSE));
numSamples=nrow(A);
sample_names=colnames(A);

# Compute Clusters
cat("Performing clustering.\n");
tree = hclust(as.dist(A), method=ClusterMethod);
tree$height=sort(tree$height);

#print(tree);
#print(tree$height);

entropy=function(x){
	x=x/sum(x);
	return(-sum(x*log(x)));
}

if(HeightCutoff==0 && KCutoff==0){
	max_height=max(tree$height);
	cat("Max Height: ", max_height, "\n");

	# Use a constant precision, independent of tree height
	length=100;
	increments=rev(seq(max_height/length,max_height,length.out=length));

	# Compute the entropy across the clusters for each cutoff
	entropies=numeric(length);
	idx=1;
	for(i in increments){
		clusters=cutree(tree, h=i);
		cluster_dist=table(clusters);
		entropies[idx]=entropy(cluster_dist);
		idx=idx+1;
	}

	# Remove zero/Nan entropies
	max_entropy=max(entropies);
	nonzero_ent_idx=entropies>0;
	nz_entropy=entropies[nonzero_ent_idx];
	nz_heights=increments[nonzero_ent_idx];

	# Compute change in entropy and index of max entropy
	delta_entropy=diff(nz_entropy/max_entropy)
	max_delta=which(delta_entropy==max(delta_entropy));

	# Computes k based on max delta entropy
	cutoff_idx=max_delta+1;
	cluster_cutoff=nz_heights[cutoff_idx];
	cat("Suggested Cluster Cutoff: ", cluster_cutoff, "\n");
	clusters=cutree(tree, h=cluster_cutoff);
}else{
	if(HeightCutoff){
		clusters=cutree(tree, h=HeightCutoff);
	}else if(KCutoff){
		clusters=cutree(tree, k=KCutoff);
	}
}

################################################################################

pdf(OutputFileNamePDF, width=11, height=8.5);

if(HeightCutoff==0 && KCutoff==0){
	# Plot entropy over cutoffs
	par(mfrow=c(2,1));
	plot(increments,entropies, main="Entropy", xlab="Height", ylab="Entropy");

	# Plot delta entropy over cutoffs
	colors=rep("black",length);
	colors[max_delta]="red";
	symbols=rep(1,length);
	symbols[max_delta]=16;
	plot(head(nz_heights, length(nz_heights)-1), delta_entropy, 
		main="Delta Entropy", xlab="Height", ylab="Percent Increase", col=colors, pch=symbols);

	par(mfrow=c(1,1));
}

# Plot colored and marked trees
max_k=max(clusters);
root=ceiling(sqrt(max_k));
colors=as.vector(t(matrix(rainbow(root^2), nrow=root)));
plot_marked_tree(tree, sample_names, clusters, colors, HeightCutoff);

dev.off();

################################################################################
# Output clusters

num_clusters=max(clusters);
for(cl in 1:num_clusters){
	outputfilename=paste(InputDistanceMatrix, ".", sprintf("%02i",cl), sep="");
	fc=file(outputfilename, "w");
	members=names(clusters[clusters==cl]);
	write(members,file=fc);
	close(fc);
}

cm=file(paste("color", ".map", sep=""), "w");
for(cl in 1:num_clusters){
	color = substr(colors[cl], 1, 7)
	cat(file=cm, sprintf("%i\t%s\n", cl, color));	
}
close(cm);

#------------------------------------------------------------------------------
# Output counts

fc=file(paste(InputDistanceMatrix, ".counts", sep=""), "w");
for(cl in 1:num_clusters){
	len=length(clusters[clusters==cl]);
	cat(file=fc, sprintf("%02i\t%i\n", cl, len));	
}
close(fc);

#------------------------------------------------------------------------------
# Output cluster members and sample ids

fc=file(paste(InputDistanceMatrix, ".groups", sep=""), "w");
width=floor(log(num_clusters, 10))+1;
fmt_str=paste("%0", width, "i", sep="");
for(cl in 1:num_clusters){
	members=names(clusters[clusters==cl]);
	num_members=length(members);
	for(m in 1:num_members){
		cat(file=fc, sprintf(fmt_str,cl) , "\t", members[m], "\n", sep="");
	}
}
close(fc);


################################################################################

cat("done.\n\n");
