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
library(seqinr);

ALPHA=0.05;
NUM_CLUST=6;
CL_METH="ward";

params=c(
                "clustal_aln", "a", 1, "character",
                "distmat", "d", 1, "character",
		"num_clusters", "k", 2, "numeric",
		"clustering_method", "m", 2, "character",
		"accession_to_sample_name_map", "c", 2, "character",
		"alpha", "p", 2, "numeric",
		"hide_identifiers", "h", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage <- paste (
                "\nUsage:\n", script_name, "\n",
		"\t-a <clustal aln file>\n",
		"\t-d <distance matrix>\n",
		"\t[-k <num clusters to compute on, default =", NUM_CLUST, ">]\n",
		"\t[-m <clustering method, default = ", CL_METH, ">]\n",
		"\t[-c <accession to name map>]\n",
		"\t[-p <p-value cutoff, default alpha = ", ALPHA, ">]\n",
		"\t[-h (hide nuemrical branch point identifiers)]\n",
		"\n",
		"Computes the differences between two clusters at the branch points.\n",
		"\n",
		"Choices for clustering method are:  'ward', 'single', 'complete', 'average', \n",
		"					'mcquitty', 'median' or 'centroid'\n",
                "\n");

if(!length(opt$clustal_aln) || !length(opt$distmat)){
        cat(usage);
        q(status=-1);
}

InputAlignment=opt$clustal_aln;
DistanceMatrix=opt$distmat;
NumClusters=NUM_CLUST;
AccessionNameMap="";

if(length(opt$accession_to_sample_name_map)){
	AccessionNameMap=opt$accession_to_sample_name_map;
}

OutputFilenameRoot=opt$distmat;
OutputFilenameRoot=gsub("\\.r_distmat$", "", OutputFilenameRoot);
OutputFilenameRoot=gsub("\\.distmat$", "", OutputFilenameRoot);

if(length(opt$num_clusters)){
	NumClusters=opt$num_clusters;
}

ClusteringMethod=CL_METH;
if(length(opt$clustering_method)){
	ClusteringMethod=opt$clustering_method;
}


AlphaCutoff=ALPHA;
if(length(opt$alpha)){
	AlphaCutoff=opt$alpha;
}

HideBranchPointIdentifiers=F;
if(length(opt$hide_identifiers)){
	HideBranchPointIdentifiers=opt$hide_identifiers;
}


cat("Input Alignment File: ", InputAlignment, "\n");
cat("Distance Matrix: ", DistanceMatrix, "\n");
cat("OutputFilenameRoot: ", OutputFilenameRoot, "\n");
cat("Alpha Cutoff: ", AlphaCutoff, "\n");

###############################################################################

# Figure out where we are running
path_comp=strsplit(script_name, "/")[[1]];
bin_comp_idx=length(path_comp);
bin_path=paste(path_comp[-bin_comp_idx], collapse="/", sep="");
if(nchar(bin_path)==0){
        bin_path=".";
}
cat("Binary path: '", bin_path, "'\n", sep="");

source(paste(bin_path, "SequenceCompare.r", sep="/"));

###############################################################################

text_scale_denfun=function(n, label_scale){
        if(is.leaf(n)){
                leaf_attr=attributes(n);
                leaf_name=leaf_attr$label;
                attr(n, "nodePar") = c(leaf_attr$nodePar,
                                        cex=0,
                                        lab.cex=label_scale);
        }
        return(n);
}

#------------------------------------------------------------------------------

get_leaves_names=function(dend){
	num_branches=length(dend);
	names=c();
	attr=attributes(dend);
	if(length(attr$leaf)){
		return(attr$label);
	}else{
		for(i in 1:num_branches){
			names=c(names, get_leaves_names(dend[[i]]));
		}	
		return(names);
	}
}

#------------------------------------------------------------------------------

get_all_heights=function(dend){
	num_branches=length(dend);
	heights=numeric();
	attr=attributes(dend);
	if(length(attr$leaf)){
		return();
	}else{
		for(i in 1:num_branches){
			heights=c(heights, get_all_heights(dend[[i]]));
		}	
		heights=c(heights, attr$height);
		return(heights);
	}
}

#------------------------------------------------------------------------------

join_list=function(list1, list2){
	len1=length(list1);
	len2=length(list2);
	if(len2>=1){
		for(i in 1:len2){
			list1[[i+len1]]=list2[[i]];
		}
	}
	return(list1);
}

get_dendrograms_above_cutoff=function(dend, cutoff){
	num_branches=length(dend);
	dend_list=list();
	attr=attributes(dend);

	if(attr$height<cutoff){
		return();
	}else{
		cat("Height: ", attr$height, "\n");
		for(i in 1:num_branches){
			dend_list=join_list(dend_list, get_dendrograms_above_cutoff(dend[[i]], cutoff));	
		}
		dend_list=join_list(dend_list, list(dend));
		return(dend_list);
	}

}

#------------------------------------------------------------------------------

make_color_map=function(cluster_assignments, names){
	num_names=length(names);
	map=list(num_names);
	for(i in 1:num_names){
		map[[names[i]]]=cluster_assignments[i];
	}
	return(map);
}

color_denfun_byIndividual=function(n, color_map){
        if(is.leaf(n)){
                leaf_attr=attributes(n);
                leaf_name=leaf_attr$label;
                ind_color=color_map[[leaf_name]];
                if(is.null(ind_color)){
                        ind_color="black";
                }

                attr(n, "nodePar") = c(leaf_attr$nodePar,
                                                list(lab.col=ind_color));
        }
        return(n);
}

#------------------------------------------------------------------------------

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

color_edges=function(den){
        consensus_color=get_color(den);
        edge_attr=attributes(den);
        attr(den, "edgePar") = c(edge_attr$edgePar, list(col=consensus_color));
        return(den);
}


###############################################################################

label_centers=function(dend, left_offset, height_cutoff){
	num_branches=length(dend);
	attr=attributes(dend);
	if(length(attr$leaf)>0 && attr$leaf==TRUE){
		return(1);
	}else{
		# label left child
		left_count=label_centers(dend[[1]], left_offset, height_cutoff);
	
		# label right child
		right_count=label_centers(dend[[2]], left_offset+left_count, height_cutoff);

		# Label self
		#text(left_offset+attr$midpoint+1, attr$height, labels=sprintf("%2.1f",attr$midpoint));
		if(attr$height>=(height_cutoff*.95)){
			divisor=10^floor(log(height_cutoff,10));
			divisor=divisor/1000;
			text(left_offset+attr$midpoint+1, attr$height, labels=sprintf("%g",floor(attr$height/divisor)), pos=3, adj=c(1,1));
		}

		return(left_count+right_count);
	}
}

###############################################################################


alignment_info=read.alignment(InputAlignment, format="clustal"); 
	# Composed of $com, $nam=seq names, $nb=num sequences, $seq

if(AccessionNameMap!=""){
	accession_map=load_accession_to_strainname_map(AccessionNameMap);
}


#------------------------------------------------------------------------------

matrix=as.matrix(read.table(DistanceMatrix, sep=" ", header=TRUE, check.names=FALSE, row.names=1));
num_samples=ncol(matrix);

# If accession map provide, rename the distance matrix now.
if(AccessionNameMap!=""){
	dup_name=list();
	rname=rownames(matrix);
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
	num_alignment_samples=length(alignment_info$nam);
	# rename the aligment file too
	for(i in 1:num_alignment_samples){
		new_name=accession_map[[alignment_info$nam[i]]];
		if (length(new_name)>0) {
			if(dup_name[[new_name]]>1){
				alignment_info$nam[i]=paste(new_name, ":", alignment_info$nam[i], "", sep="");
			}else{
				alignment_info$nam[i]=accession_map[[alignment_info$nam[i]]];
			}
		}
	}
	rownames(matrix)=rname;
}

align_matrix=seqinr_to_matrix(alignment_info);
align_matrix=clean_leadtrail_gaps(align_matrix);

sample_names=rownames(matrix);

distance=as.dist(matrix);
hclust=hclust(distance, method=ClusteringMethod);
hclust$height=sort(hclust$height);

dendro=as.dendrogram(hclust);

#names=get_leaves_names(dendro);
#print(names);

heights=get_all_heights(dendro);
sorted_heights=sort(heights, decreasing=TRUE);

sorted_heights=sorted_heights[sorted_heights>0];
if(NumClusters>length(sorted_heights)){
	cat("Number of clusters requested is less than the number available.\n");
}
cutoff=min(head(sorted_heights, NumClusters));

cat("Height Cutoff: ", cutoff, "\n");

clusters=cutree(hclust, h=cutoff*(1-1e-10));

num_clusters=length(unique(clusters))-1;
color_map=make_color_map(clusters, sample_names);

#color_map;

# Color dendrograms
dendro=dendrapply(dendro, color_denfun_byIndividual, color_map);
dendro=dendrapply(dendro, color_edges);

cut_dendros=get_dendrograms_above_cutoff(dendro, cutoff);

# Open output files
pdf(paste(OutputFilenameRoot, ".res_change_annot.pdf", sep=""), height=8.5, width=11);
fh=file(paste(OutputFilenameRoot, ".res_change_annot.csv", sep=""), "w");

crayola16=c("#1F75FE", "#B4674D", "#1CAC78", "#FF7538", "#EE204D", "#926EAE", "#7366BD", "#FFAACC", "#0D98BA", "#C0448F", "#FF5349");

palette(crayola16);

align_length=ncol(align_matrix);
mtc_alpha=AlphaCutoff/(align_length*num_clusters); # Multiple Testing Corrected Alpha
cat("Multiple Testing (Bonferroni) Corrected Alpha: ", mtc_alpha, "\n");

layout_matrix=t(matrix(c(
	1,1,1,1,2,2,
	1,1,1,1,2,2
), nrow=6));
layout(layout_matrix);
#print(layout_matrix);
for(i in num_clusters:1){
	
	cat("Working on cluster: ", i, "\n");

	# Get subtree to process
	cur_den=cut_dendros[[i]];

	# Rescale names
	attr=attributes(cur_den);
	text_scale=min(1, 40/attr$members);
	cur_den=dendrapply(cur_den, text_scale_denfun, text_scale);

	# Plot this subtree
	par(mar=c(11,1,3,.5));
	par(oma=c(.5,1,3,0));
	par(family="");
	plot(cur_den, main=OutputFilenameRoot);
	abline(h=cutoff, col="grey", lty=2);
	par(family="");

	points(attributes(cur_den)$midpoint+1, attributes(cur_den)$height, pch=17, cex=2);

	if(!HideBranchPointIdentifiers){
		label_centers(cur_den,0, cutoff);
	}

	
	# Get members from both sides of tree
	left_leaves=get_leaves_names(cur_den[[1]]);
	right_leaves=get_leaves_names(cur_den[[2]]);

	cat("-------------------------------------------------------------------\n");
	cat("i: ", i, "\n");
	#cat("left:\n");
	#print(left_leaves);
	#cat("right:\n");
	#print(right_leaves);
	
	mut_count=1;
	mut_info=c();
	par(mar=c(0,0,.5,.5));
	plot(0,0, xlim=c(-10,100), ylim=c(-100, 0), ylab="", xlab="", type="n", xaxt="n", yaxt="n", bty="n");
	par(family="Courier");
	
	text(0,0, "  Pos  LPerc  Res  RPerc  Res  p-val", pos=4, font=2);
	for(pos in 1:align_length){
		ct=build_contingency_table(align_matrix, left_leaves, right_leaves, pos);
		if(nrow(ct)>1){
			ft=fisher.test(ct);
			pval=ft$p.val;
			if(pval<mtc_alpha){

				ct_summary=summarize_contingency_table(ct);

				g1_prop=ct_summary[1];
				g1_nuc=ct_summary[2];
				g2_prop=ct_summary[3];
				g2_nuc=ct_summary[4];

				ct_summary_str=paste(
					g1_prop, "% ", 
					g1_nuc, " -> ", 
					g2_prop, "% ", 
					g2_nuc, sep="");
			
				mut_info=paste(sprintf("%5i", pos), "  ", ct_summary_str, "   ", 
					sprintf("%10.3e", pval), "\n", sep="");

				cat(file=fh, paste(
					attributes(cur_den)$height,
					pos,
					g1_prop, g1_nuc, g2_prop, g2_nuc, 
					sprintf("%5.4g",pval) , sep=","), "\n");

				text(0,-(mut_count+1)*1.5, mut_info, pos=4)
				mut_count=mut_count+1;
			}
		}
	}
	cat(file=fh, "\n");
}

close(fh);

###############################################################################
###############################################################################

dev.off();
cat("Done.\n");
