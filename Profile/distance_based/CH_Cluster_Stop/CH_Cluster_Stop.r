#!/usr/bin/env Rscript

###############################################################################

cat("\n\n");

library(MASS)
library('getopt');
library('vegan');
library('plotrix');
library('parallel');


DEF_DISTTYPE="man";
DEF_NUM_TOP_CAT=35;
DEF_NUM_CLUS=-1;
DEF_SPLIT_CHAR=";";

DEF_NUM_BS=320;
DEF_NUM_PROC=25;

params=c(
	"input_summary_table", "i", 1, "character",
	"output_filename_root", "o", 2, "character",
	"dist_type", "d", 2, "character",
	"max_clusters", "K", 2, "numeric",
	"num_sub_samp", "s", 2, "numeric",
	"sample_inclusion_fname", "c", 2, "character",
	"num_bootstraps", "b", 2, "numeric",
	"tag_name", "t", 2, "character",
	"num_proc", "p", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\n\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	[-o <output file root name, default is input file base name>]\n",
	"	[-d <euc/wrd/man/bray/horn/bin/gow/tyc/minkp5/minkp3, default =", DEF_DISTTYPE, ">]\n",
	"\n",
	"	[-K <maximum number of clusters to cut, default=ceiling(log2(num_samples))>]\n",
	"	[-s <sub sample size, default=all>]\n",
	"	[-c <sample inClusion filename>]\n",
	"	[-b <num bootstraps, default=", DEF_NUM_BS, ">]\n",
	"	[-t <tag name>]\n",
	"\n",
	"	[-p <num processors/threads to use for bootstrapping, default=", DEF_NUM_PROC, ">]\n",
	"\n",
	"This script identifies the recommended number of clusters to split the data\n",
	"into based on the Calinski-Harabasz statistic.  It is a Pseudo-F statistic\n",
	"in the sense that is based on the F-statistic, but doesn't enforce any of\n",
	"the standard rules of F-statistic applicability.  In any case, the F statistic\n",
	"tries to measure the signal-to-noise ratio (SNR) for a particular number of clusters\n",
	"where the separation between clusters is the Signal (Sum of Squares Between) and\n",
	"the noise is the separation within clusters (Sum of Squares Within), while\n",
	"adjusting for the number of clusters (degrees of freedom). \n",
	"\n",
	"The basical idea is the number of clusters that is idea, is the one that maximizes\n",
	"the SNR, or Pseudo-F stat across all the number of clusters.\n",
	"\n",
	"This script will:\n",
	"	1.) Read in a summary table and compute a full/complete distance matrix.\n",
	"	2.) Hierarchically cluster with Wards, then use cuttree to generate clustering.\n",
	"	3.) Generate num_bootstrap resamples of the samples.\n",
	"	3.) For each cluster, compute CH statistic, for each bootstrap\n",
	"	4.) For each number of clusters, k, plot 95% CI for each CH stat\n",
	"	5.) Across each of the bootstraps, find max CH stat, and plot histogram\n",
	"\n",
	"For the distance types:\n",
	" minkp5 is the minkowski with p=.5, i.e. sum((x_i-y_i)^1/2)^2\n",
	" minkp3 is the minkowski with p=.3, i.e. sum((x_i-y_i)^1/3)^3\n",
	"\n",
	"The sample inclusion filename -c option specifies a file that contains all the\n",
	"sample names to include.\n",
	"\n");

if(!length(opt$input_summary_table)){
	cat(usage);
	q(status=-1);
}
InputFileName=opt$input_summary_table;

if(length(opt$output_file)>0){
	OutputFileRoot=opt$output_file;
}else{
	OutputFileRoot=InputFileName;
	OutputFileRoot=gsub(".summary_table.xls$", "", OutputFileRoot);
	OutputFileRoot=gsub(".summary_table.tsv$", "", OutputFileRoot);
	cat("No output file root specified.  Using input file name as root.\n");
}

dist_type=DEF_DISTTYPE;
if(length(opt$dist_type)){
	dist_type=opt$dist_type;
}

num_bs=DEF_NUM_BS;
if(length(opt$num_bootstraps)){
	num_bs=opt$num_bootstraps;
}

sub_sample=0;
if(length(opt$num_sub_samp)){
	sub_sample=opt$num_sub_samp;
}

if(!any(dist_type == c("wrd","man","bray","horn","bin","gow","euc","tyc","minkp3","minkp5"))){
	cat("Error: Specified distance type: ", dist_type, " not recognized.\n");
	quit(status=-1);
}

SampleIncFname="";
if(length(opt$sample_inclusion_fname)){
	SampleIncFname=opt$sample_inclusion_fname;
}

MaxClusters=-1;
if(length(opt$max_clusters)){
	MaxClusters=opt$max_clusters;
}

NumProc=DEF_NUM_PROC;
if(length(opt$num_proc)){
	NumProc=opt$num_proc;
}

if(length(opt$tag_name)){
        TagName=opt$tag_name;
        cat("Setting TagName Hook: ", TagName, "\n");
        setHook("plot.new",
                function(){
                        cat("Hook called.\n");
                        if(par()$page==T){
                                oma_orig=par()$oma;
                                exp_oma=oma_orig;
                                exp_oma[1]=max(exp_oma[1], 1);
                                par(oma=exp_oma);
                                mtext(paste("[", TagName, "]", sep=""), side=1, line=exp_oma[1]-1,
                                        outer=T, col="steelblue4", font=2, cex=.8, adj=.97);
                                par(oma=oma_orig);
                        }
                }, "append");

}else{
        TagName="";
}

###############################################################################
# See http://www.mothur.org/wiki/Thetayc for formula

tyc_fun=function(v1, v2){
	sum_intersect=sum(v1*v2);
	sum_sqrd_diff=sum((v1-v2)^2);
	denominator=sum_sqrd_diff + sum_intersect;
	tyc=1-(sum_intersect/denominator);
	return(tyc);
}

thetaYC=function(matrix){
	
	nsamples=nrow(matrix);
	ycdist=matrix(0, ncol=nsamples, nrow=nsamples);
	for(i in 1:nsamples){
		for(j in 1:nsamples){
			if(i<j){
				ycdist[i,j]=tyc_fun(matrix[i,], matrix[j,]);
			}else{
				ycdist[i,j]=ycdist[j,i];			
			}
		}
	}
	
	as.dist(return(ycdist));
}

weight_rank_dist_opt=function(M, deg){
        NumSamples=nrow(M);
        order_matrix=matrix(0, nrow=nrow(M), ncol=ncol(M));
        for(i in 1:NumSamples){
                order_matrix[i,]=rank(M[i,], ties.method="average");
        }

        dist_mat=matrix(0, nrow=NumSamples, ncol=NumSamples);
        colnames(dist_mat)=rownames(M);
        rownames(dist_mat)=rownames(M);
        for(i in 1:NumSamples){
                for(j in 1:i){
                        dist_mat[i,j]=
                                sqrt(sum((
                                        (order_matrix[i,]-order_matrix[j,])^2)*
                                        (((M[i,]+M[j,])/2)^deg)
                                        )
                                );
                }
        }
        return(as.dist(dist_mat));

}


###############################################################################

load_summary_table=function(st_fname){
	inmat=as.matrix(read.delim(st_fname, sep="\t", header=TRUE, row.names=1, check.names=FALSE, comment.char="", quote=""))

	num_categories=ncol(inmat)-1;
	num_samples=nrow(inmat);

	cat("Loaded Summary Table: ", st_fname, "\n", sep="");
	cat("  Num Categories: ", num_categories, "\n", sep="");
	cat("  Num Samples: ", num_samples, "\n", sep="");

	countsmat=inmat[,2:(num_categories+1)];

	return(countsmat);
}

#------------------------------------------------------------------------------

normalize=function(st){
	num_samples=nrow(st);
	num_categories=ncol(st);

	normalized=matrix(0, nrow=num_samples, ncol=num_categories);
	colnames(normalized)=colnames(st);
	rownames(normalized)=rownames(st);

	sample_counts=apply(st, 1, sum);
	for(i in 1:num_samples){
		normalized[i,]=st[i,]/sample_counts[i];
	}
	return(normalized);
}

#------------------------------------------------------------------------------

plot_text=function(strings){
        orig_par=par(no.readonly=T);

        par(mfrow=c(1,1));
        par(family="Courier");
        par(oma=rep(.5,4));
        par(mar=rep(0,4));

        num_lines=length(strings);

        top=max(as.integer(num_lines), 40);

        plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                );
        for(i in 1:num_lines){
                #cat(strings[i], "\n", sep="");
                text(0, top-i, strings[i], pos=4, cex=.8);
        }

        par(orig_par);
}

#------------------------------------------------------------------------------

compute_dist=function(norm_st, type){

	if(type=="euc"){
		dist_mat=dist(norm_st);
	}else if (type=="wrd"){
		dist_mat=weight_rank_dist_opt(norm_st, deg=4);
	}else if (type=="man"){
		dist_mat=vegdist(norm_st, method="manhattan");
	}else if (type=="bray"){
		dist_mat=vegdist(norm_st, method="bray");
	}else if (type=="horn"){
		dist_mat=vegdist(norm_st, method="horn");
	}else if (type=="bin"){
		dist_mat=vegdist(norm_st, method="bin");
	}else if (type=="gow"){
		dist_mat=vegdist(norm_st, method="gower");
	}else if (type=="tyc"){
		dist_mat=thetaYC(norm_st);
	}else if (type=="minkp3"){
		dist_mat=dist(norm_st, method="minkowski", p=1/3);
	}else if (type=="minkp5"){
		dist_mat=dist(norm_st, method="minkowski", p=1/2);
	}

	dist_mat[dist_mat==0]=1e-323;

	return(dist_mat);
}

#------------------------------------------------------------------------------

# This is the single-threaded base function
compute_pseudoF=function(dist_mat, memberships, resample_sequence){

	dist_mat=as.matrix(dist_mat);

	num_samples=length(memberships);
	num_clusters=length(unique(memberships));

	n=num_samples;
	k=num_clusters;

	sample_list=names(memberships);

	num_bs=nrow(resample_sequence);
	pseudo_f=numeric(num_bs);

	for(bs_ix in 1:num_bs){
		SSB=0;
		SSW=0;

		# Resample
		resamp_ix=resample_sequence[bs_ix,]; #sample(num_samples, replace=T);

		# Compute SSB and SSW
		for(i in 1:num_samples){
			samp_i=sample_list[resamp_ix[i]];	

			for(j in 1:num_samples){

				if(i<j){

					samp_j=sample_list[resamp_ix[j]];

					i_mem=memberships[samp_i];
					j_mem=memberships[samp_j];

					#cat(i_mem, "/", j_mem, "\n");
					if(i_mem==j_mem){
						SSW=SSW+dist_mat[samp_i, samp_j];
						#cat("SSW: ", SSW, "\n");
					}else{
						SSB=SSB+dist_mat[samp_i, samp_j];
						#cat("SSB: ", SSB, "\n");
					}
				}
			}
		}

		# Compute/Store this instance of bootstrap 
		pseudo_f[bs_ix]=SSB/(k-1)*(n-k)/SSW;
	}

	return(pseudo_f);
}

# This is the multi-threaded function
compute_pseudoF_parallel=function(dist_mat, memberships, resample_sequence, nproc=20){

	# wrap compute function down to a single parameter / list item
	comp_pF=function(params){
		l_dist_mat=params[["dist_mat"]];	
		l_memberships=params[["memberships"]];	
		l_res_seq=params[["resample_sequence"]];	
		l_batch_id=params[["id"]];
		cat("Starting new process for: ", l_batch_id, "\n");
		psfs=compute_pseudoF(l_dist_mat, l_memberships, l_res_seq);		
		cat("Done (", l_batch_id, ")\n");
		return(psfs);
	}

	num_bs=nrow(resample_sequence);

	cat("Number of Bootstraps Requested: ", num_bs, "\n");
	cat("Number of Processors Supplied: ", nproc, "\n");

	splits=ceiling(seq(0, num_bs, length.out=(nproc+1)));
	starts=splits[1:nproc]+1;
	ends=splits[2:(nproc+1)];
	split_mat=matrix(c(starts, ends), ncol=2, byrow=F);
	colnames(split_mat)=c("start", "end");

	cat("Bootstrap Splits:\n");
	print(split_mat);
	
	# Place the 3 parameter needed for each thread into a list
	# The resample sequence is what is split across threads, other parameters are
	#   the same.
	input_list=list();
	for(batch_ix in 1:nproc){
		bat_start=split_mat[batch_ix, "start"];
		bat_end=split_mat[batch_ix, "end"];

		params=list();
		params[["id"]]=batch_ix;
		params[["dist_mat"]]=dist_mat;
		params[["memberships"]]=memberships;
		params[["resample_sequence"]]=resample_sequence[bat_start:bat_end,,drop=F];

		input_list[[batch_ix]]=params
	}

	#print(input_list);

	cat("'Making' Compute Cluster...\n");
	cl=makeCluster(nproc);

	clusterExport(cl, c("compute_pseudoF"));

	cat("Running parLapply...\n");
	results=parLapply(cl, input_list, comp_pF);

	cat("Stopping Cluster...\n");
	stopCluster(cl);

	cat("Joining results across splits...\n");
	pseudoF_array=do.call(c, results);

	cat("Accumulated length of results: ", length(pseudoF_array), "\n");

	return(pseudoF_array);

}


###############################################################################

color_denfun_bySample=function(n){
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

###############################################################################

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

get_clstrd_leaf_names=function(den){
# Get a list of the leaf names, from left to right
        den_info=attributes(den);
        if(!is.null(den_info$leaf) && den_info$leaf==T){
                return(den_info$label);
        }else{
                lf_names=character();
                for(i in 1:2){
                        lf_names=c(lf_names, get_clstrd_leaf_names(den[[i]]));
                }
                return(lf_names);
        }
}

get_middle_of_groups=function(clustered_leaf_names, group_asgn){
# Finds middle of each group in the plot
        num_leaves=length(group_asgn);
        groups=sort(unique(group_asgn));
        num_groups=length(groups);

        reord_grps=numeric(num_leaves);
        names(reord_grps)=clustered_leaf_names;
        reord_grps[clustered_leaf_names]=group_asgn[clustered_leaf_names];

        mids=numeric(num_groups);
        names(mids)=1:num_groups;
        for(i in 1:num_groups){
                mids[i]=mean(range(which(reord_grps==i)));
        }
        return(mids);

}

color_denfun_bySample=function(n){
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

text_scale_denfun=function(n){
        if(is.leaf(n)){
                leaf_attr=attributes(n);
                leaf_name=leaf_attr$label;
                attr(n, "nodePar") = c(leaf_attr$nodePar,
					cex=0,
                                        lab.cex=denfun.label_scale);
        }
        return(n);
}

color_leaves=function(n){
        if(is.leaf(n)){
                leaf_attr=attributes(n);
                leaf_name=leaf_attr$label;
                attr(n, "nodePar") = c(leaf_attr$nodePar,
                                        col=as.character(denfun.label_col_map[leaf_name]),
                                        cex=denfun.label_scale*2.75,
                                        pch=15
                                        );
        }
        return(n);
}

color_edges=function(den){
        edge_attr=attributes(den);
        attr(den, "edgePar") = c(edge_attr$edgePar, list(col=consensus_color));
        return(den);
}

###############################################################################

output_fname_root = paste(OutputFileRoot, ".", dist_type, sep="");

pdf(paste(output_fname_root, ".ch_stop.pdf", sep=""), height=8.5, width=14);

parameter_summary_text=capture.output({
	cat("\n");
	cat("Input Summary Table Name: ", InputFileName, "\n", sep="");
	cat("Output Filename Root: ", output_fname_root, "\n", sep="");
	cat("Distance Type: ", dist_type, "\n", sep="");
	cat("Num Bootstraps: ", num_bs, "\n", sep="");
	cat("Sample Inc List: ", SampleIncFname, "\n", sep="");
	cat("Subsample Size: ", sub_sample, "\n", sep="");
	cat("Max Clusters: ", MaxClusters, "\n", sep="");
	cat("Num Processors: ", NumProc, "\n", sep="");
	cat("Tag Name:", TagName, "\n", sep="");
	cat("\n");
});

print(parameter_summary_text, quote=F);
plot_text(parameter_summary_text);

cat("Loading summary table...\n");
counts_mat=load_summary_table(InputFileName);

if(SampleIncFname!=""){
	cat("Loading sample inclusion list...\n");
	samp_incl_list=scan(SampleIncFname, "character");
	sumtab_samp_ids=rownames(counts_mat);
	shared_samp_ids=intersect(samp_incl_list, sumtab_samp_ids);
	num_shared=length(shared_samp_ids);
	cat("Num samples overlapping summary table and inclusion list: ", num_shared, "\n");
	counts_mat=counts_mat[shared_samp_ids,,drop=F];
}

if(sub_sample>0){
	cat("Subsampling counts table to: ", sub_sample, "\n", sep="");
	counts_mat=counts_mat[sample(nrow(counts_mat),sub_sample),];
}

# Normalize counts
cat("Normalizing counts...\n");
norm_mat=normalize(counts_mat);
#print(norm_mat);

# Compute Avg Abund
avg_abd=apply(norm_mat, 2, mean);
rank=order(avg_abd, decreasing=T);

# Reorder norm matrix
norm_mat=norm_mat[, rank, drop=F];

category_names=colnames(norm_mat);
sample_names=rownames(norm_mat);
num_samples=nrow(norm_mat);
num_categories=ncol(norm_mat);

###############################################################################

# Compute full distances
cat("Computing distances...\n");
full_dist_mat=compute_dist(norm_mat, dist_type);

###############################################################################

hcl=hclust(full_dist_mat, method="ward.D2");

# Find height where cuts are made
max_clusters=ifelse(MaxClusters==-1, ceiling(log(num_samples, 2)), MaxClusters);
cat("Max Clusters to compute: ", max_clusters, "\n");
cut_midpoints=numeric(max_clusters);
for(k in 2:max_clusters){
        cut_midpoints[k]=find_height_at_k(hcl, k);
}

orig_dendr=as.dendrogram(hcl);

# Extract names of leaves from left to right
lf_names=get_clstrd_leaf_names(orig_dendr);

palette_col=c("red", "green", "blue", "cyan", "magenta", "orange", "gray", "pink", "black", "purple", "brown", "aquamarine");
palette(palette_col);

# Compute ISO and classical MDS
nonparm_mds_res=matrix(0,nrow=num_samples,ncol=2);
classic_mds_res=matrix(0,nrow=num_samples,ncol=2);

cat("Computing classic and non parametric MDS...\n");
imds_res=isoMDS(full_dist_mat);
nonparm_mds_res[,1]=imds_res$points[,1]
nonparm_mds_res[,2]=imds_res$points[,2]
classic_mds_res=cmdscale(full_dist_mat);

###############################################################################

denfun.label_scale=min(2,55/num_samples);

###############################################################################

dendro_layout=matrix(c(
	1,1,1,1,1,1,1,1,1,1,1,1,2
), byrow=T, ncol=13);


mds_layout=matrix(c(
	1,1,1,1,2,2,2,2,3,
	1,1,1,1,2,2,2,2,3,
	1,1,1,1,2,2,2,2,3,
	1,1,1,1,2,2,2,2,3
), byrow=T, ncol=9);


#------------------------------------------------------------------------------

# Pre-randomize samples across all cluster sizes so finding max will be consistent
resample_sequence=matrix(0, nrow=num_bs, ncol=num_samples);
# Rows are each BS instance
# Columns are the selected samples
for(i in 1:num_bs){
	resample_sequence[i,]=sample(num_samples, replace=T);
}
#print(resample_sequence);

pseudoF_mat=matrix(0, nrow=num_bs, ncol=max_clusters);

membership_matrix=matrix(1, nrow=num_samples, ncol=max_clusters);
rownames(membership_matrix)=lf_names;
colnames(membership_matrix)=1:max_clusters;
print(membership_matrix);

#------------------------------------------------------------------------------

for(num_cl in 2:max_clusters){

        # Cut tree at target number of clusters
        cat("Cutting for ", num_cl, " clusters...\n", sep="");
        memberships=cutree(hcl, k=num_cl);

	# Compute CH stats
	cat("  Computing PseudoF Stats...\n");
	#pseudoF_res=compute_pseudoF(full_dist_mat, memberships, resample_sequence);
	pseudoF_res=compute_pseudoF_parallel(full_dist_mat, memberships, resample_sequence, NumProc);
	pseudoF_mat[,num_cl]=pseudoF_res;

        grp_mids=get_middle_of_groups(lf_names, memberships);

        # Reorder cluster assignments to match dendrogram left/right
        plot_order=order(grp_mids);
        mem_tmp=numeric(num_samples);
        for(gr_ix in 1:num_cl){
                old_id=(memberships==plot_order[gr_ix]);
                mem_tmp[old_id]=gr_ix;
        }
        names(mem_tmp)=names(memberships);
        memberships=mem_tmp;
        grp_mids=grp_mids[plot_order];

	# Store memberships for later export
	membership_matrix[lf_names, num_cl]=memberships[lf_names];

        # Plot Dendrogram
        par(oma=c(4,0,5,0));
        par(mar=c(5.1,4.1,4.1,2.1));
	layout(dendro_layout);
        sample_to_color_map=as.list(memberships);

        tweaked_dendro=dendrapply(orig_dendr, color_denfun_bySample);
        tweaked_dendro=dendrapply(tweaked_dendro, text_scale_denfun);

	cat("  Plotting dendrogram...\n");
        plot(tweaked_dendro, horiz=F, lwd=2);
        for(cl_ix in 1:num_cl){
                lab_size=3/ceiling(log10(cl_ix+1));
                axis(side=1, outer=T, at=grp_mids[cl_ix], labels=cl_ix, cex.axis=lab_size, col.ticks=cl_ix,
                        lend=1, lwd=10, padj=1, line=-3);
        }
        abline(h=cut_midpoints[num_cl], col="red", lty=2);

        par(mar=c(0,0,0,0));
        plot(0, type="n", xlab="", ylab="", main="", bty="n", xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1));
        legend(0,1, fill=1:num_cl, legend=c(as.character(1:num_cl)), bty="n", cex=2);

        mtext(paste("Distance Type: ", dist_type), side=3, line=0, outer=T);
        mtext(paste("Num Clusters: ", num_cl), side=3, line=1, outer=T);

	layout(mds_layout);
	
        # Generate MDS plots
        par(oma=c(0,0,4,0));
        par(mar=c(5.1,4.1,4.1,2.1));
	cat("  Plotting MDS Plots...\n");
        plot(nonparm_mds_res, col=memberships, xlab="Dim 1", ylab="Dim 2", main="non-metric MDS", cex=2);
        plot(classic_mds_res, type="n", col=memberships, xlab="Dim 1", ylab="Dim 2", main="classical MDS");
        points(classic_mds_res, col=memberships, cex=2);

        # MDS Legend
        par(mar=c(0,0,0,0));
        plot(0, type="n", xlab="", ylab="", main="", bty="n", xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1));
        legend(0,1, fill=1:num_cl, legend=c(as.character(1:num_cl)), bty="n", cex=2);

        mtext(paste("Distance Type: ", dist_type, sep=""), side=3, outer=T, line=.5);
        mtext(paste("Num Clusters: ", num_cl, sep=""), side=3, outer=T, line=2);

}


cat("\n");
cat("PseudoF Matrix:\n");
print(pseudoF_mat);
cat("\n");

# Remove infinite values, probably due to low sample size
pseudoF_mat=apply(pseudoF_mat, 1:2, function(x){ if(is.infinite(x)){return(NA);}else{return(x);}});

med_pseudoF=apply(pseudoF_mat, 2, function(x){quantile(x, .5, na.rm=T)});
lb_pseudoF=apply(pseudoF_mat, 2, function(x){quantile(x, .025, na.rm=T)});
ub_pseudoF=apply(pseudoF_mat, 2, function(x){quantile(x, .975, na.rm=T)});
max_pseudoF=apply(pseudoF_mat, 1, function(x){min(which(max(x,na.rm=T)==x))});

max_med_ix=which(max(med_pseudoF, na.rm=T)==med_pseudoF);

###############################################################################

# Plot CH over clusters
par(mfrow=c(1,1));
par(mar=c(8,8,3,3));

plotCI(x=2:max_clusters, y=med_pseudoF[2:max_clusters], li=lb_pseudoF[2:max_clusters], ui=ub_pseudoF[2:max_clusters],
	sfrac=0.01,
	xaxt="n",
	main="Cluster Separation w/ 95% CI",
	xlab="Number of Clusters (k)", ylab="Calinski-Harabasz Pseudo-F Statistic");
axis(side=1, at=2:max_clusters, labels=2:max_clusters);
mtext(paste("Num Bootstraps per cluster cut: ", nrow(pseudoF_mat), sep=""), side=3, line=0, cex=.75);
abline(h=med_pseudoF[max_med_ix], lty=2, col="red");

# Plot Histogram for best
#print(max_pseudoF);
quants=quantile(max_pseudoF, c(.025, .5, .975));
hist(max_pseudoF, breaks=c(1, (1:max_clusters)+.5), 
	xlab="Num Clusters at Max(CH PseudoF)",
	main="Recommended Number of Clusters", xaxt="n");
mtext(paste("Median = ", quants[2], sep=""), side=3, line=-1, cex=.75);
mtext(paste("95% CI = (", quants[1], ", ", quants[3], ")", sep=""), side=3, line=-2, cex=.75);
axis(side=1, at=1:max_clusters, labels=1:max_clusters);
abline(v=quants[2], lwd=2, col="blue4");
abline(v=quants[c(1,3)], lty=2, col="blue");

###############################################################################

membership_matrix[lf_names, num_cl]=memberships[lf_names];

colnames(membership_matrix)=sprintf("k%02i", as.numeric(colnames(membership_matrix)));
membership_matrix=apply(membership_matrix, c(1,2), function(x){sprintf("cl%02i", as.numeric(x));});

print(membership_matrix);
mat_rnames=rownames(membership_matrix);
out_mat=cbind(mat_rnames, membership_matrix);
colnames(out_mat)=c("SampleID", colnames(membership_matrix));

write.table(out_mat, file=paste(OutputFileRoot, ".cls_mem.csv", sep=""),
	quote=F, sep=",", row.names=F, col.names=T);

write.table(out_mat, file=paste(OutputFileRoot, ".cls_mem.tsv", sep=""),
	quote=F, sep="\t", row.names=F, col.names=T);

###############################################################################

cat("\nDone.\n")
dev.off();
warn=warnings();
if(length(warn)){
        print(warn);
}
q(status=0);
