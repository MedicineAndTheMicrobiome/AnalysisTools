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

params=c(
                "matrix_a", "A", 1, "character",
                "matrix_b", "B", 1, "character",
		"output_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage <- paste (
                "\nUsage:\n", script_name, "\n",
                "       -A <input distance matrix file A>\n",
                "       -B <input distance matrix file B>\n",
		"	-o <output filename root>\n",
		"\n",
		"Reads in two distance matrices and compares them with Mantel's Test.\n",
                "\n",
                "\n");

if(
	!length(opt$matrix_a) ||
	!length(opt$matrix_b) ||
	!length(opt$output_root)
){
        cat(usage);
        q(status=-1);
}

MatrixA=opt$matrix_a;
MatrixB=opt$matrix_b;
OutputFileRoot=opt$output_root;
OutputPDF=paste(OutputFileRoot, ".dist_comp.pdf", sep="");
OutputTXT=paste(OutputFileRoot, ".dist_comp.txt", sep="");

cat("Matrix A: ", MatrixA, "\n");
cat("Matrix B: ", MatrixB, "\n");
cat("Output Root: ", OutputFileRoot, "\n");

MatrixAShortName=tail(strsplit(MatrixA,"/")[[1]], 1);
MatrixBShortName=tail(strsplit(MatrixB,"/")[[1]], 1);

###############################################################################

library(MASS);
library(vegan);

###############################################################################
# Plots a dendrogram that has been colored by cluster assignment
plot_colored_dendrogram=function(dendro_in, name_to_color_mapping, title, k){

        # Resize labels based on number of samples
	names_length=length(name_to_color_mapping);
        label_scale=33/names_length;
        if(label_scale>2){
                label_scale=2;
        }

        # Insert color attribute into node
        color_denfun=function(n){
                if(is.leaf(n)){
                        leaf_attr=attributes(n);
                        leaf_name=leaf_attr$label;

                        color_assignment=name_to_color_mapping[[leaf_name]];

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
        plot(dendro_out, main=title, yaxt="n");
	mtext(paste("k = ", k), side=2)
}

###############################################################################

# Load data
dist_mat_A=as.matrix(read.table(MatrixA, header=TRUE, check.names=FALSE));
dist_mat_B=as.matrix(read.table(MatrixB, header=TRUE, check.names=FALSE));

n_samples_A=nrow(dist_mat_A);
n_samples_B=nrow(dist_mat_B);

cat("Num samples in A: ", n_samples_A, "\n");
cat("Num samples in B: ", n_samples_B, "\n");

names_A=colnames(dist_mat_A);
names_B=colnames(dist_mat_B);

#print(names_A);
#print(names_B);

# Compute common samples
both_names=sort(intersect(names_A,names_B));
print(both_names);
n_samples_common=length(both_names);
cat("Num samples common to A&B: ", n_samples_common, "\n");

# Generate common distance matrix
common_dist_mat_A=dist_mat_A[both_names, both_names];
common_dist_mat_B=dist_mat_B[both_names, both_names];

# Compute hierarchical clusters
hclust_A=hclust(as.dist(common_dist_mat_A), method="ward");
hclust_B=hclust(as.dist(common_dist_mat_B), method="ward");
dendro_A=as.dendrogram(hclust_A);
dendro_B=as.dendrogram(hclust_B);

# Compute correlation coefficient between identical samples of the different matrices
corr_arr=numeric(n_samples_common);
for(i in 1:n_samples_common){
	corr_arr[i]=cor(common_dist_mat_A[i,], common_dist_mat_B[i,], method="spearman");
}
sorted_corr=sort(corr_arr, decreasing=TRUE, index.return=TRUE);

# Compute mantel statistic and significance
man_out=mantel(common_dist_mat_A, common_dist_mat_B, method="spearman");
print(man_out);
significance=man_out$signif;
statistic=man_out$statistic;

###############################################################################
# Generate plots

heatcol=(grey(1:100/100));

a_stdev=sd(common_dist_mat_A);
a_norm=(common_dist_mat_A)/a_stdev;

b_stdev=sd(common_dist_mat_B);
b_norm=(common_dist_mat_B)/b_stdev;

pdf(OutputPDF, height=11, width=8.5);

par(mfrow=c(3,1));
hclust_order=hclust_A$order;
par(mar=c(1,3,1,3))
image(a_norm[hclust_order, hclust_order], col=heatcol, main=MatrixAShortName, xaxt="n", yaxt="n");
image(b_norm[hclust_order, hclust_order], col=heatcol, main=MatrixBShortName, xaxt="n", yaxt="n");
par(mar=c(14,3,5,3))
image(matrix(corr_arr[hclust_order], ncol=1), col=heatcol, xaxt="n", yaxt="n");
axis(side=1, labels=both_names[hclust_order], at=seq(0,1,length.out=n_samples_common), las=2);
mtext(side=3, paste("Mantel Statistic:", statistic), line=4);
mtext(side=3, paste("Mantel Significance:", significance), line=3);

alpha=.05;
if(significance<alpha){
		mtext(side=3, paste("At alpha=", alpha, ", the two distance matrices are correlated.", sep=""), line=1);
}else{
		mtext(side=3, paste("At alpha=", alpha, ", the two distance matrices are NOT correlated.", sep=""), line=1);
}

#------------------------------------------------------------------------------

par(oma=c(1,2,0,0));
layout_matrix=matrix(c(1,2,3,3), byrow=T, nrow=2);
layout(layout_matrix);

#------------------------------------------------------------------------------
plot_color_correlation=function(correlation, names, sort_order, color_mapping){
	num_samples=length(names);
	color=character(num_samples);
	for(i in 1:num_samples){
		color[i]=color_mapping[[names[sort_order[i]]]]
	}
	barplot(correlation[sort_order], names.arg=names[sort_order], 
		ylim=c(min(correlation)*1.1,1.1), col=color, las=2,
		main="Correlation between Distances"
		);
}

#------------------------------------------------------------------------------
# Plot A to B Colors
for(k in 2:5){
	clusters=cutree(hclust_A, k=k);

	color_mapping=list(n_samples_common);
	color_palette=rainbow(k);
	for(i in 1:n_samples_common){
		color_mapping[[both_names[i]]]=color_palette[clusters[i]];
	}

	par(mar=c(7,0,1,0));
	plot_colored_dendrogram(dendro_A, color_mapping, MatrixAShortName, k);
	plot_colored_dendrogram(dendro_B, color_mapping, MatrixBShortName, k);
	par(mar=c(13,1,1,0));
	plot_color_correlation(corr_arr, both_names, sorted_corr$ix, color_mapping);
}

#------------------------------------------------------------------------------
# Plot B to A Colors
for(k in 2:5){
	clusters=cutree(hclust_B, k=k);

	color_mapping=list(n_samples_common);
	color_palette=rainbow(k);
	for(i in 1:n_samples_common){
		color_mapping[[both_names[i]]]=color_palette[clusters[i]];
	}

	par(mar=c(7,0,1,0));
	plot_colored_dendrogram(dendro_A, color_mapping, MatrixAShortName, k);
	plot_colored_dendrogram(dendro_B, color_mapping, MatrixBShortName, k);
	par(mar=c(13,1,1,0));
	plot_color_correlation(corr_arr, both_names, sorted_corr$ix, color_mapping);
}

#------------------------------------------------------------------------------

fh=file(OutputTXT, "w");

cat(file=fh, paste("#Keyword", "MatAName", "MatBName", "MantelStatistic", "Significance", "NumCommonSamples", sep=","), "\n");
cat(file=fh, paste("MantelStat", MatrixAShortName, MatrixBShortName, statistic, significance, n_samples_common, sep=","), "\n");

close(fh);

###############################################################################

cat("Done.\n");


