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
                "input_file", "i", 1, "character",
                "output_file", "o", 2, "character",
		"SSB_alpha_threshold", "t", 2, "numeric",
		"max_clusters", "m", 2, "numeric",
		"num_bootstraps", "b", 2, "numeric",
		"output_partition_members", "p", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage <- paste (
                "\nUsage:\n", script_name, "\n",
                "       -i <input distance matrix file>\n",
		"	[-t <SSB threshold, default alpha=.10>]\n",
		"	[-o <output file root name, default: input file name>]\n",
		"	[-m <max clusters to analyze to, default: 35>]\n",
		"	[-b <num bootstraps to perform, default: 200>]\n",
		"	[-p (output partition members into files flag)]\n",
		"\n",
		"Reads in a generic square distance matrix file and\n",
		"estimates the number of clusters that should be formed\n",
		"based on an analysis of extreme SS between samples.\n",
                "\n",
                "\n");

if(!length(opt$input_file)){
        cat(usage);
        q(status=-1);
}

InputFilename=opt$input_file;

OutputFilenameRoot=opt$input_file;
if(length(opt$output_file)>0){
	OutputFilenameRoot=opt$output_file;
}

OutputPDF=paste(OutputFilenameRoot, ".SSAllocClusters.pdf", sep="");
OutputTXT=paste(OutputFilenameRoot, ".SSAllocClusters.txt", sep="");

if(length(opt$max_clusters)>0){
	MaxClusters=opt$max_clusters;
}else{
	MaxClusters=35;
}

if(length(opt$num_bootstraps)>0){
	NumBootstraps=opt$num_bootstraps;
}else{
	NumBootstraps=200;
}

if(length(opt$SSB_alpha_threshold)>0){
	SSB_alpha_threshold=opt$SSB_alpha_threshold;
}else{
	SSB_alpha_threshold=.10;
}

if(length(opt$output_partition_members)){
	OutputPartitionMembers=TRUE;
}else{
	OutputPartitionMembers=FALSE;
}


cat("Input File Name: ", InputFilename, "\n");
cat("Output Filename Root: ", OutputFilenameRoot, "\n");
cat("Output PDF: ", OutputPDF, "\n");
cat("\n");
cat("Max Clusters: ", MaxClusters, "\n");
cat("Num Bootstraps: ", NumBootstraps, "\n");
cat("SSB threshold (alpha): ", SSB_alpha_threshold, "\n");
cat("Outputing Parition Members: ", OutputPartitionMembers, "\n");

###############################################################################

library(MASS);

###############################################################################
# Computes the sum of squared distances (According to Calinski and Harabasz)

sum_sqrs=function(dist_mat){
	if(nrow(dist_mat)==1){
		return(0);
	}else{
		n=ncol(dist_mat);
		ss=sum(as.dist(dist_mat)^2);
		mean=ss/n;
		return(mean);
	}
}

###############################################################################
# Computes the mean squared distance for a given full matrix

mean_sqrd_dist=function(dist_mat){
	if(nrow(dist_mat)==1){
		return(-1);
	}else{
		half_mat=as.dist(dist_mat);
		n=length(half_mat);
		ssd=sum(half_mat^2);
		mean=ssd/n;
		return(mean);
	}
}

###############################################################################
# Makes non diagonal distances non zero, so isoMDS won't complain

make_nonzero=function(in_mat){
	for(i in 1:nrow(in_mat)){
		for(j in 1:ncol(in_mat)){
			if(i!=j){
				if(in_mat[i,j]==0){
					in_mat[i,j]=.0000000000000000001;
				}
			}
		}
	}
	return(in_mat);
}

###############################################################################
# Plots a dendrogram that has been colored by cluster assignment

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

plot_marked_tree=function(hclust, sample_names, cluster, colors){

	dendro=as.dendrogram(hclust);
	plot_colored_dendrogram(dendro, sample_names, cluster, colors);
	dendro_attrib=attributes(dendro);
	k=length(unique(cluster));

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

}

###############################################################################
# Plots a MDS plot

plot_colored_mds=function(mds_result, cluster, colors){
	plot(mds_result$points[,1], mds_result$points[,2], col=colors[cluster], xlab="Dimension 1", ylab="Dimension 2");
}

###############################################################################
# Computes the fstat for a set of clusters and a distance matrix (According to Calinski and Harabasz)

compute_fstat=function(dist_mat, clusters){
        NumSamples=nrow(dist_mat);
        num_clusters=max(clusters);

        #cat("Num Samples: ", NumSamples, "\n");
        #cat("Num Clusters: ", num_clusters, "\n");
        #print(clusters);

        #-----------------------------------------------------------------------------
        # Compute SSW: Within SS
        gss=numeric(num_clusters);
        for(cl in 1:num_clusters){
                #cat("cluster: ", cl, "\n");
                grp_idx=which(clusters==cl);
                group_dist_mat=as.matrix(dist_mat[grp_idx, grp_idx], nrow=length(grp_idx));
                gss[cl]=sum_sqrs(group_dist_mat);
        }
        wgss=sum(gss);
        #cat("WGSS: ", wgss, "\n");

        #-----------------------------------------------------------------------------
        # Compute SST: Total SS
        tss=sum_sqrs(dist_mat);
        #cat("TSS: ", tss, "\n");

        #-----------------------------------------------------------------------------
        # Compute F-stat, Total SS = Within SS + Between SS
        bgss=tss-wgss;
        k=num_clusters;
        n=NumSamples;
        fstat=(bgss/(k-1))/(wgss/(n-k));
        #print(fstat);

        #-----------------------------------------------------------------------------
        # Output F-stat and components
        out=list();
        out$SSB=bgss;
        out$SSW=wgss;
        out$SST=tss;
        out$Fstat=fstat;
        return(out);
}

###############################################################################
# Computes the threshold for SSB/SST by looking at the SS in the top alpha of all distances

compute_SSB_threshold=function(sqrd_dist_mat, alpha_cutoff){
	if((alpha_cutoff>=1)||(alpha_cutoff<=0)){
		cat("Error, sample cutoff should be between 0 and 1.\n");
		cat("In particular, it should be around .05 or .1.\n");
		return(NaN);
	}
	half_mat=as.dist(sqrd_dist_mat);
	sorted_distances=sort(half_mat);
	num_distances=length(half_mat);

	#cat("Num distances (half matrix): ", num_distances, "\n");
	
	bottom_n=ceiling((1-alpha_cutoff)*num_distances);
	#cat("Num samples in bottom ", sprintf("%3.2f%%", (1-alpha_cutoff)*100), ": ", bottom_n, "\n", sep="");

	bottom_ssd=sum(sorted_distances[1:bottom_n]);
	total_ssd=sum(sorted_distances);
	portion_ssd_in_bottom=bottom_ssd/total_ssd;
	SSB_threshold=1-portion_ssd_in_bottom;

	#cat("Bottom SSD: ", bottom_ssd, "\n");
	#cat("Total SSD: ", total_ssd, "\n");
	#cat("Perc of Bottom: ", sprintf("%3.2f%%", portion_ssd_in_bottom*100), "\n");
	#cat("SSB Threshold: ",  sprintf("%3.2f%%", SSB_threshold*100), "\n");

	return(SSB_threshold);

}

###############################################################################

get_CI = function(x, alpha){
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

###############################################################################
###############################################################################

# Load data
dist_mat=as.matrix(read.table(InputFilename, header=TRUE, check.names=FALSE));

nrow=nrow(dist_mat);
ncol=ncol(dist_mat);

if(nrow!=ncol){
	cat("Error: The distance matrix is not square.\n");
}

NumSamples=nrow;
cat("Num Samples: ", NumSamples, "\n", sep="");

# Determine how the max range of num clusters to look at
if(MaxClusters>NumSamples){
	MaxClusters=NumSamples;
}

sample_names=colnames(dist_mat);

# Eliminate 0 distance between non self-self
dist_mat=make_nonzero(dist_mat);

# Cluster distances
hcluster=hclust(as.dist(dist_mat),"ward");

# Generate clusters for each K
cluster_assgn=list();
for(k in 1:MaxClusters){
	cluster_assgn[[k]]=cutree(hcluster, k=k);
}

# Allocate arrays for each k for each of the measures
ssb_over_sst=rep(0, MaxClusters);

for(k in 2:MaxClusters){
	cat(".");

	# Grab the cluster assignments
	cl_asgn=cluster_assgn[[k]];	

	# Compute SS statistics
	ss=compute_fstat(dist_mat, cl_asgn);

	# Compute SSB/SST, percent of variance assigned to intercluster
	ssb_over_sst[k]=ss$SSB/ss$SST;

}
cat("\n");

###############################################################################

plot_hist_with_alpha=function(x, alpha){
	cat("Plotting histogram with alpha = ", alpha, "\n");
	h=hist(x, breaks=200, plot=FALSE);
	one_minus_alpha=1-alpha;
	num_bins=length(h$counts);
	cumsum=cumsum(h$counts)
	total=sum(h$counts);
	norm_cum=cumsum/total;
	col=rep("black", sum(norm_cum<one_minus_alpha));
	col=c(col, rep("red", sum(norm_cum>=one_minus_alpha)));
	barplot(h$counts, 
		names.arg=h$mids,
		#xlab="Distances", las=2,
		ylab="Counts",
		col=col, 
		border=col,
		space=0,
		main=sprintf("Histogram of Distances: \nRed is Alpha = %3.2f", alpha));
}

###############################################################################

plot_sorted_SD_with_alpha=function(d, alpha){
	cat("Plotting sorted SD with alpha = ", alpha, "\n");
	
	num_distances=length(d);
	count_cutoff=num_distances*(1-alpha);
	d_sorted=sort(d, decreasing=FALSE);
	d_sqrd=d_sorted^2;
	#cumsum_d_sqrd=cumsum(d_sqrd);
	#total_ssd=sum(d_sqrd);
	colors=rep("red", num_distances);
	colors[1:count_cutoff]="black";
	subidx=seq(1,num_distances, length.out=500);
	barplot(d_sqrd[subidx], xlab="Sorted Distances", ylab="Distances Squared", 
		main="Sorted Squared Distances", 
		col=colors[subidx], border=colors[subidx], space=0);
	#plot(d_sqrd[subidx], type="p", xlab="Sorted d^2", ylab="d^2", col=colors[subidx] );
	#abline(v=count_cutoff, col="red");
		
}

###############################################################################

recommend_k=function(ssb_threshold, ssb_over_sst){
	num_k=length(ssb_over_sst);
	rec_k=min(which(ssb_over_sst>=ssb_threshold));
	if(is.nan(rec_k)){
		cat("SSB threshold is greater than last calculated SSB/SST for k's (", num_k, ") calculated.\n");
	}
	return(rec_k);
}


###############################################################################

pdf(OutputPDF, height=18.5, width=11.5);

# Compute MDS layout
mds=isoMDS(dist_mat);

sqrd_dist=dist_mat^2;
half_sqrd_dist=as.dist(sqrd_dist);
half_dist=as.dist(dist_mat);

SSB_threshold=compute_SSB_threshold(sqrd_dist, SSB_alpha_threshold);

rec_k=recommend_k(SSB_threshold, ssb_over_sst);
effective_ssb=ssb_over_sst[rec_k];

#-------------------------------------------------------------------------------
par(mfrow=c(3,1));

# Plot histogram of distances with alpha
plot_hist_with_alpha(half_dist, SSB_alpha_threshold);

# Plot CDF of SSs with alpha
plot_sorted_SD_with_alpha(half_dist, SSB_alpha_threshold);

# Plot SSW versus SSB ratio
SS_prop=c(1-SSB_threshold, SSB_threshold);
pos=barplot(SS_prop, col=c("black", "red"), names.arg=c("SSW", "SSB"), 
	main="SSW and SSB Proportions",
	ylab="Proportion", ylim=c(0, max(SS_prop)*1.3));
text(pos, SS_prop, labels=sprintf("%3.2f%%", 100*SS_prop), cex=1.1, pos=3);

#-------------------------------------------------------------------------------
par(mfrow=c(3,1));

cat("SSB Threshold:", sprintf("%3.2f%%", SSB_threshold*100), "\n");

plot_all=function(attribute, title, hcluster, mds, sample_names, threshold=NULL){
	par(mfrow=c(3,1));
	k=length(attribute);

	attribute[is.nan(attribute)]=0;
	if(is.null(threshold)){
		# Find max attribute and k for it
		max_att=max(attribute);
		max_k=min(which(max_att==attribute));
	}else{
		max_k=min(which(attribute>threshold));
		max_att=attribute[max_k];
	}
	
	# Plot k versus attibute
	plot(2:k, attribute[2:k], main=title, xlab="k", ylab="SSB/SST");
	text(2:k, attribute[2:k], 2:k, pos=3);			# Label each point by num clusters, k
	points(max_k, max_att, col="red", pch=16, cex=2); 	# Highlight best point

	# Draw tree/mds colored by best cluster recommendation
	cluster=cutree(hcluster, k=max_k);
	root=ceiling(sqrt(max_k));
	colors=as.vector(t(matrix(rainbow(root^2), nrow=root)));
	plot_marked_tree(hcluster, sample_names, cluster, colors);
	plot_colored_mds(mds, cluster, colors);
}

plot_all(ssb_over_sst, sprintf("alpha: %1.2f  SSB/SST > %2.2f%%", SSB_alpha_threshold, 100*SSB_threshold), hcluster, mds, sample_names, SSB_threshold);


###############################################################################

cat("\nPerforming bootstrapping to generate confidence intervals around k:\n");
bs_rec_k=numeric(0);
num_samples=nrow(dist_mat);
for(bs in 1:NumBootstraps){
	#cat("\nBootstrap Iteration: ", bs, "\n");
	cat(".");
	bs_resample=sample(1:num_samples, num_samples, replace=TRUE);
	bs_distmat=dist_mat[bs_resample, bs_resample];
	bs_hcluster=hclust(as.dist(bs_distmat),"ward");
	bs_sqrd_dist=bs_distmat^2;
	bs_SSB_threshold=compute_SSB_threshold(bs_sqrd_dist, SSB_alpha_threshold);
	for(bs_k in 2:num_samples){
		bs_clusters=cutree(bs_hcluster, bs_k);
		bs_ss=compute_fstat(bs_distmat, bs_clusters);
		bs_ss_alloc=bs_ss$SSB/bs_ss$SST;
	
		#cat("Allocation : ", bs_ss_alloc, "\n");
		if(bs_ss_alloc>=bs_SSB_threshold){
			bs_rec_k[bs]=bs_k;
			break;
		}
	}
}
cat("\n");

###############################################################################
# Generate histogram on bootstapping

hist(bs_rec_k, breaks=.5+(0:(max(bs_rec_k))), main="Bootstrapped Recommended K Values", xlab="Recommended K");
intervals=get_CI(bs_rec_k, .05);
cat("Median: ", intervals[1], "\n");
cat("LowerBound: ", intervals[2], "\n");
cat("UpperBound: ", intervals[3], "\n");
mtext(sprintf("Median: %f", intervals[1]), side=3, line=-1);
mtext(sprintf("LB: %f", intervals[2]), side=3, line=-2);
mtext(sprintf("UB: %f", intervals[3]), side=3, line=-3);

###############################################################################

cluster_sizes=table(cluster_assgn[[rec_k]]);
cluster_size_str=paste(cluster_sizes, collapse=";");
percent_matching_input=sum(bs_rec_k==rec_k)/NumBootstraps;

fh=file(OutputTXT, "w");
cat(file=fh, paste("#InputFileName", "alpha", "Rec_k", "Cluster_Sizes", "Median_k", "95%LB_k", "95%UB_k", "PercOfRecK", "Eff_SSB", "Thres_SSB", sep=","), "\n", sep="");
cat(file=fh, paste(InputFilename, SSB_alpha_threshold, rec_k, cluster_size_str, intervals[1], intervals[2], intervals[3], percent_matching_input, effective_ssb, SSB_threshold, sep=","), "\n", sep="");
close(fh);

###############################################################################

if(OutputPartitionMembers){
	cat("Outputing partition members and counts...\n");

	members_fname_root=OutputFilenameRoot;
	counts_fname=paste(OutputFilenameRoot, ".counts", sep="");
	counts_fh=file(counts_fname, "w");

	rec_clus=cluster_assgn[[rec_k]];

	string_width=ceiling(log(rec_k, 10));

	for(i in 1:rec_k){
		cat("Working on cluster: ", i, "\n");
		istring=sprintf(paste("%0", string_width, "i", sep=""),i);
		members_fname=paste(members_fname_root, ".", istring, sep="");
		members_fh=file(members_fname, "w");
		members_idx=which(rec_clus==i);
		num_members=length(members_idx);
		cat("  Num Members: ", num_members, "\n");
		for(j in 1:num_members){
			cat(file=members_fh, sample_names[members_idx[j]], "\n", sep="");
		}

		cat(file=counts_fh, istring, "\t", num_members, "\n");
	}
	
}

dev.off();

cat("Done.\n");
