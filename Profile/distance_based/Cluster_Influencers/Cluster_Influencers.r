#!/usr/bin/env Rscript

###############################################################################

cat("\n\n");

library(MASS)
library('getopt');
library('vegan');
library('plotrix');


DEF_DISTTYPE="man";
DEF_NUM_TOP_CAT=35;
DEF_NUM_CLUS=-1;
DEF_SPLIT_CHAR=";";

params=c(
	"input_summary_table", "i", 1, "character",
	"output_filename_root", "o", 2, "character",
	"sample_id_list", "l", 2, "character",
	"dist_type", "d", 2, "character",
	"num_top_cat", "p", 2, "numeric",
	"split_char", "s", 2, "character",
	"max_cuts", "k", 2, "numeric",
	"only_at_k", "K", 2, "numeric",
	"factor_filename", "f", 2, "character",
	"factor_names_list", "n", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\n\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	[-o <output file root name, default is input file base name>]\n",
	"	[-l <sample id list, otherwise use all samples in summary table>]\n",
	"	[-d <euc/wrd/man/bray/horn/bin/gow/tyc/minkp5/minkp3, default =", DEF_DISTTYPE, ">]\n",
	"	[-p <num of top categories to probe, default =", DEF_NUM_TOP_CAT, " >\n",
	"	[-s <split char for long category names, default =", DEF_SPLIT_CHAR, ">]\n",
	"\n",
	"	Clustering-Based (Hierarchical Clustering) Options:\n",
	"	[-k <max number of clusters to split into, default =", DEF_NUM_CLUS, ">\n",
	"	[-K <only compute at K cuts>\n",
	"\n",
	"	Metadata-Based (User-defined Factors) Options:\n",
	"	[-f <factor/metadata file]\n",
	"	[-n <file containing column names to analyze in factor file>]\n",
	"\n",
	"This script will:\n",
	"	1.) Read in a summary table and compute a full/complete distance matrix.\n",
	"	2.) Precompute reduced distance matrices for each of the top categories to probe.\n",
	"	3.) Cluster full distance matrix hierarchically.\n",
	"\n",
	"	If Clustering-based:\n",
	"	4.) Iteratively cut tree until max_c clusters.\n",
	"	If Metadata-based:\n",
	"	4.) Step through each of the groupings specified in the metadata.\n",
	"\n",
	"	5.) For each pair of clusters, compute:\n",
	"		a.) compute complete and reduced R^2=SSB/SST\n",
	"		b.) compute log10 Ratio between reduced(R^2)/complete(R^2)\n",
	"		c.) compute which cluster has greater amount of that category\n",
	"		d.) Accumulate log(R^2 ratios) across clusters to estimate cluster differentiators\n",
	"\n",
	"For the distance types:\n",
	" minkp5 is the minkowski with p=.5, i.e. sum((x_i-y_i)^1/2)^2\n",
	" minkp3 is the minkowski with p=.3, i.e. sum((x_i-y_i)^1/3)^3\n",
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

if(length(opt$sample_id_list)){
	SampleIDListFname=opt$sample_id_list;
}else{
	SampleIDListFname="";
}

dist_type=DEF_DISTTYPE;
if(length(opt$dist_type)){
	dist_type=opt$dist_type;
}

if(!any(dist_type == c("wrd","man","bray","horn","bin","gow","euc","tyc","minkp3","minkp5"))){
	cat("Error: Specified distance type: ", dist_type, " not recognized.\n");
	quit(status=-1);
}

max_clusters=DEF_NUM_CLUS;
if(length(opt$max_cuts)){
	max_clusters=opt$max_cuts;
}

if(length(opt$only_at_k)){
	OnlyAtK=opt$only_at_k;
}else{
	OnlyAtK=0;
}	


num_top_cat=DEF_NUM_TOP_CAT;
if(length(opt$num_top_cat)){
	num_top_cat=opt$num_top_cat;
}

SplitChar=DEF_SPLIT_CHAR;
if(length(opt$split_char)){
	SplitChar=opt$split_char;
}

FactorFilename="";
if(length(opt$factor_filename)){
	FactorFilename=opt$factor_filename;
}

FactorListFile="";
if(length(opt$factor_names_list)){
	FactorListFile=opt$factor_names_list;
}

useMetadata=T;
w_meta_ext="";
if(FactorFilename==""){
	useMetadata=F;
}else{
	w_meta_ext=".meta";
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

load_list=function(list_fname){
	list=read.delim(list_fname, sep="\t", header=F, row.names=NULL, as.is=T, check.names=F, comment.char="#", quote="");
	return(list[,1]);
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

compute_anova=function(dist_mat, members_a, members_b){

	# Compute SST (SS Total)
	# Compute SSW (SS Error)
	# Compute SSB (SS Treatment)

	# Fstat = 

	dist_mat=as.matrix(dist_mat);

	SSB=0;
	SSW=0;

	between_cts=0;
	within_cts=0;

	for(memb1 in members_a){
		for(memb2 in members_a){
			if(memb1!=memb2){
				SSW=SSW+dist_mat[memb1, memb2];
				within_cts=within_cts+1;
			}
		}
	}

	for(memb1 in members_b){
		for(memb2 in members_b){
			if(memb1!=memb2){
				SSW=SSW+dist_mat[memb1, memb2];
				within_cts=within_cts+1;
			}
		}
	}

	for(memb1 in members_a){
		for(memb2 in members_b){
			SSB=SSB+dist_mat[memb1, memb2];
			between_cts=between_cts+1;
		}
	}

	# Correct double adding of members (that were on the other side of diagonal)
	SSW=SSW/2;
	within_cts=within_cts/2;

	SST=SSW+SSB;

	MSW=SSW/within_cts;
	MSB=SSB/between_cts;
	
	pseudoF=MSB/MSW;
	Rsqrd=SSB/SST;

	results=list();
	results$SSW=SSW;
	results$SSB=SSB;
	results$SST=SST;
	results$MSW=MSW;
	results$MSB=MSB;
	results$pseudoF=pseudoF;
	results$Rsqrd=Rsqrd;
	
	return(results);

}

###############################################################################

analyze_cluster_pairs=function(st, members_a, members_b, num_cat_to_probe, dist_mat_list){

	num_mem_a=length(members_a);
	num_mem_b=length(members_b);
	cat_names=colnames(st);

	cat("Num members in clusters: ", num_mem_a, " vs. ", num_mem_b, "\n");

	cat("\nMembers A:\n");
	print(members_a);
	cat("\nMembers B:\n");
	print(members_b);
	cat("\n\n");

	complete_dist=dist_mat_list[["full"]];
	complete_anova=compute_anova(complete_dist, members_a, members_b);
	cat("Complete ANOVA:\n");
	print(complete_anova);

	ratios=matrix(0, nrow=num_cat_to_probe, ncol=1);
	rownames(ratios)=cat_names[1:num_cat_to_probe];
	colnames(ratios)=c("rsqrd");

	for(i in 1:num_cat_to_probe){
		cat("Leaving out: ", cat_names[i], "\n", sep="");
		part_dist=dist_mat_list[[i]];
		part_anova=compute_anova(part_dist, members_a, members_b);
		cat("Partial ANOVA:\n");
		#print(part_anova);
		ratios[i,"rsqrd"]=part_anova$Rsqrd/complete_anova$Rsqrd;
		#print(ratios);
	}

	results=list();
	results[["partial_anova"]]=part_anova;
	results[["ratios"]]=ratios;
	return(results);

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

load_factor_file=function(fn){
        inmat=read.delim(fn, sep="\t", header=TRUE, row.names=1, check.names=F, comment.char="", quote="");

        # Changes spaces to underscore
        var_names=colnames(inmat);
        var_names=gsub(" ", "_", var_names);
        colnames(inmat)=var_names;

        cat("  Num Factors: ", ncol(inmat), "\n", sep="");
        cat("  Num Samples: ", nrow(inmat), "\n", sep="");
        return(inmat);
}

reconcile_factors_and_summary_table=function(fact, st){

	cat("\nReconciling Factor and Summary Tables...\n");
	num_fact_samp=nrow(fact);	
	num_st_samp=nrow(st);	

	fact_samp_ids=rownames(fact);
	st_samp_ids=rownames(st);

	shared_samp_ids=intersect(fact_samp_ids, st_samp_ids);
	fact_only_sample_ids=setdiff(fact_samp_ids, shared_samp_ids);
	st_only_sample_ids=setdiff(st_samp_ids, shared_samp_ids);
	num_shared_samp=length(shared_samp_ids);

	results=list();
	results[["summary_table"]]=st[shared_samp_ids,, drop=F];
	results[["factors_table"]]=fact[shared_samp_ids,, drop=F];
	results[["summary"]]=c(
		"Original Number of Samples:",
		paste("  Summary Table: ", num_st_samp, sep=""),
		paste("  Factors Table: ", num_fact_samp, sep=""),
		"",
		paste("Shared Number of Samples: ",  num_shared_samp, sep=""),
		"",
		"Samples exclusive to Summary Table:",
		capture.output(print(st_only_sample_ids)),
		"Samples exclusive to Factors Table:",
		capture.output(print(fact_only_sample_ids)),
		""
	);
	return(results);
}

###############################################################################

if(OnlyAtK){
	kext=paste(".k", OnlyAtK, sep="");
}else{
	kext="";
}

output_fname_root = paste(OutputFileRoot, ".", dist_type, w_meta_ext, kext, sep="");
cat("\n");
cat("Input Summary Table Name: ", InputFileName, "\n", sep="");
cat("Output Filename Root: ", output_fname_root, "\n", sep="");
cat("Distance Type: ", dist_type, "\n", sep="");
cat("Num Top categories to analyze: ", num_top_cat, "\n", sep="");
cat("Targeted Variables: ", FactorListFile, "\n", sep="");

if(OnlyAtK){
	cat("Only computing at: ", OnlyAtK, "\n", sep="");
}
cat("\n");

if(useMetadata){
	cat("Loading Factor Table...\n");
	factors_matrix=load_factor_file(FactorFilename);
	
	if(FactorListFile!=""){
		cat("Loading Targeted Factor List...\n");
		target_factors=load_list(FactorListFile);
		cat("Extracting factors of interest:\n");
		print(target_factors);
		factors_matrix=factors_matrix[,target_factors, drop=F];
	}
	target_factors=colnames(factors_matrix);
	num_target_factors=ncol(factors_matrix);
}else{
	cat("Max Num clusters: ", max_clusters, "\n", sep="");
}

###############################################################################

cat("Loading summary table...\n");
counts_mat=load_summary_table(InputFileName);
#print(counts_mat);

if(SampleIDListFname!=""){
	sample_keep_list=load_list(SampleIDListFname);
	st_sample_ids=rownames(counts_mat);
	keep_list=intersect(st_sample_ids, sample_keep_list);
	counts_mat=counts_mat[keep_list,,drop=F];
}

# Reconcile
if(useMetadata){
	reconciliation_res=reconcile_factors_and_summary_table(factors_matrix, counts_mat);
	counts_mat=reconciliation_res$summary_table;
	factors_matrix=reconciliation_res$factors_table;
	cat(reconciliation_res$summary, sep="\n");
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

cat("\nTop Categories: \n");
print(category_names[1:num_top_cat]);
cat("\n");

# Shorten category names
short_cat_names=character();
for(i in 1:num_categories){
	short_cat_names[i]=tail(strsplit(category_names[i], SplitChar)[[1]],1);
	short_cat_names[i]=gsub("_unclassified$", "_uncl", short_cat_names[i]);
	short_cat_names[i]=gsub("_group", "_grp", short_cat_names[i]);
}
colnames(norm_mat)=short_cat_names;
cat("Shorted Top Categories: \n");
print(short_cat_names[1:num_top_cat]);
cat("\n");

###############################################################################
###############################################################################

# Compute full distances
cat("Computing distances...\n");
full_dist_mat=compute_dist(norm_mat, dist_type);

# Compute partial distances
cat("Precomputing distances without top categories...\n");
dist_mat_list=list();
for(target_cat_ix in 1:num_top_cat){
	cat("\t", target_cat_ix, ".) ", short_cat_names[target_cat_ix], "\n", sep="");
	dist_mat_list[[target_cat_ix]]=compute_dist(norm_mat[,-target_cat_ix], dist_type);
}

dist_mat_list[["full"]]=full_dist_mat;

###############################################################################

hcl=hclust(full_dist_mat, method="ward.D2");

# Find height where cuts are made
if(max_clusters==-1){
	max_clusters=ceiling(log(num_samples,2));
}
cut_midpoints=numeric(max_clusters);
for(k in 2:max_clusters){
        cut_midpoints[k]=find_height_at_k(hcl, k);
}

orig_dendr=as.dendrogram(hcl);

# Exactract names of leaves from left to right
lf_names=get_clstrd_leaf_names(orig_dendr);

pdf(paste(output_fname_root, ".cl_inf.pdf", sep=""), height=8.5, width=14);
palette_col=c("red", "green", "blue", "cyan", "magenta", "orange", "gray", "pink", "black", "purple", "brown", "aquamarine");
#palette_col=c("blue", "red", "green", "cyan", "magenta", "orange", "gray", "pink", "black", "purple", "brown", "aquamarine");
palette(palette_col);

# Compute ISO and classical MDS
nonparm_mds_res=matrix(0,nrow=num_samples,ncol=2);
classic_mds_res=matrix(0,nrow=num_samples,ncol=2);

imds_res=isoMDS(full_dist_mat);
nonparm_mds_res[,1]=imds_res$points[,1]
nonparm_mds_res[,2]=imds_res$points[,2]

classic_mds_res=cmdscale(full_dist_mat);

if(useMetadata){
	# Give more room for factor level names
	mds_layout=matrix(c(
		1,1,1,1,2,2,2,2,3,3), byrow=T, nrow=1, ncol=10);
}else{
	mds_layout=matrix(c(
		1,1,1,1,2,2,2,2,3), byrow=T, nrow=1, ncol=9);
}

barplot_layout=matrix(c(
	1,3,5,7,9,11,
	1,3,5,7,9,11,
	1,3,5,7,9,11,
	2,4,6,8,10,12), byrow=T, nrow=4, ncol=6);
barplots_per_page=6;


# Compute necessary radii of grid lines based on range of MDS points
mds_range=range(classic_mds_res);
num_grid_lines=6;
grid_lines=(1:num_grid_lines)*(mds_range[2]-mds_range[1])/((num_grid_lines-2)*2);
cat("Grid Line radii:\n");
print(grid_lines);

# For a run with metadata, we'll look across the number of factors.
# For a run without metadata, we'll cut the tree until the max clusters.

if(useMetadata){
	num_iterations=num_target_factors;
	target_factors=colnames(factors_matrix);
	num_target_factors=ncol(factors_matrix);
}else{
	num_iterations=max_clusters-1;
}

if(OnlyAtK!=0){
	iterations=OnlyAtK;
}else{
	iterations=1:num_iterations;
}

cat("Using Metadata:", useMetadata, "\n");

# Begin pair-wise cluster analyses
for(ix in iterations){

	if(useMetadata){
		fact_val=factors_matrix[,ix];
		fact_name=colnames(factors_matrix)[ix];
		num_samps=length(fact_val);

		memberships=rep(NA, num_samples);
		names(memberships)=rownames(factors_matrix);
		levels=sort(unique(fact_val));

		num_levels=length(levels);
		if(num_levels>10){
			cat("Factor has ", num_levels, " levels.\n", sep="");
			print(levels);
			cat("Error: Too many factor levels!\n");
			next;
		}else if(num_levels==1){
			cat("Factor has only 1 level.  Skipping...\n");
			next;
		}

		for(lix in 1:num_levels){
			memberships[levels[lix]==fact_val]=lix;
		}

		memberships=memberships[!is.na(memberships)];
		num_cl=num_levels;
		grp_mids=NA;

		legend_labels=paste(1:num_cl, ": ", levels, sep="");
	}else{
		#for(num_cl in 2:max_clusters):
		if(OnlyAtK>0){
			num_cl=OnlyAtK;
		}else{
			num_cl=ix+1;
		}

		cat("Cutting for ", num_cl, " clusters...\n", sep="");
		memberships=cutree(hcl, k=num_cl);
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

		legend_labels=as.character(1:num_cl);
	}

	# Plot Dendrogram
	cat("Prepping dendrogram...\n");
	max_sample_name_length=max(nchar(sample_names));
	cat("Max sample name length: ", max_sample_name_length, "\n");
	label_scale=min(2,20/max_sample_name_length);
	par(oma=c(0,0,0,0));
	par(mar=c(max_sample_name_length/2,4.1,4.1,2.1));
	par(mfrow=c(1,1));
	sample_to_color_map=as.list(memberships);
	tweaked_dendro=dendrapply(orig_dendr, color_denfun_bySample);
	tweaked_dendro=dendrapply(tweaked_dendro, text_scale_denfun);
	plot(tweaked_dendro);
	
	if(!useMetadata){
		# Label cluster numbers underneath the samples
		for(cl_ix in 1:num_cl){
			lab_size=3/ceiling(log10(cl_ix+1));
			axis(side=1, outer=T, at=grp_mids[cl_ix], labels=cl_ix, cex.axis=lab_size, col.ticks=cl_ix,
				lend=1, lwd=10, padj=1, line=-3);
		}
		abline(h=cut_midpoints[num_cl], col="red", lty=2);
	}

	# Legend at top left of dendrogram
	ranges=par()$usr;
	legend(ranges[1], ranges[4], fill=1:num_cl, legend=legend_labels, bty="n");

	# Label page with cut/factor specific info
	mtext(paste("Distance Type: ", dist_type), side=3, line=1, outer=T);
	mtext(paste("Num Clusters: ", num_cl), side=3, line=0, outer=T);
	if(useMetadata){
		mtext(paste("Factor Name: ", fact_name), side=3, line=-1, outer=T);
	}

	# Generate MDS plots
	par(oma=c(0,0,2,0));
	par(mar=c(5.1,4.1,4.1,2.1));
	layout(mds_layout);

	label_centroids=function(points, memberships){
                members=sort(unique(memberships));
                for(i in members){
                        cur_grp=(i==memberships);
                        midx=mean(points[cur_grp,1]);
                        midy=mean(points[cur_grp,2]);
                        text(midx, midy, label=i, col="black", cex=3/ceiling(log10(i+1)));
                }

        }

	plot(nonparm_mds_res, col=memberships, xlab="Dim 1", ylab="Dim 2", main="non-metric MDS");
	label_centroids(nonparm_mds_res, memberships);
	plot(classic_mds_res, type="n", col=memberships, xlab="Dim 1", ylab="Dim 2", main="classical MDS");
	for(grid_radius in grid_lines){
		draw.ellipse(0,0, a=grid_radius, b=grid_radius, border="grey90");
	}
	points(classic_mds_res, col=memberships);
	label_centroids(classic_mds_res, memberships);

	# Plot legend in far right
	par(mar=c(0,0,0,0));
	plot(0, type="n", xlab="", ylab="", main="", bty="n", xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1));

	legend(0,1, fill=1:num_cl, legend=legend_labels, bty="n", cex=2);
	mtext(paste("Num Clusters: ", num_cl), side=3, outer=T);
	
	# Compute R^2 pairwise between clusters
	ratios_list=list();
	log_ratios_ranges=numeric();
	for(i in 1:num_cl){

		members_i=names(memberships[memberships==i]);

		for(j in 1:num_cl){

			if(i>=j){
				next;
			}

			cat("Working on cluster[", i, "] vs cluster[", j, "]\n");

			members_j=names(memberships[memberships==j]);

			sub_mat=norm_mat[c(members_i, members_j),,drop=F];
			results=analyze_cluster_pairs(sub_mat, members_i, members_j, num_top_cat, dist_mat_list);

			lograt=log10(results$ratios[,"rsqrd"]);
			ratios_list[[ paste(i, "#", j, sep="")]]=lograt;
			log_ratios_ranges=range(log_ratios_ranges, range(lograt));

		}

	}

	# Generate R^2 Ratio Plots
	layout(barplot_layout);
	par(oma=c(.5,13,3,1));
	plot_count=0;
	for(i in 1:num_cl){

		members_i=names(memberships[memberships==i]);

		for(j in 1:num_cl){

			if(i>=j){
				next;
			}

			members_j=names(memberships[memberships==j]);
			both_mem=c(members_i, members_j);
		
			# Only plot the left most category labels
			if(!(plot_count %% barplots_per_page)){
				plot_cat_names=short_cat_names[1:num_top_cat];
			}else{
				plot_cat_names=rep("",num_top_cat);
			}

			# Pull up the previously calculated ratios
			ratios=ratios_list[[ paste(i, "#", j, sep="")]];
			
			neg_ratios_ix=which(ratios<0);
			pos_contrib=rep("grey", num_top_cat);
			pos_contrib[neg_ratios_ix]="red";

			# Compute differences in abundance.
			i_means=apply(norm_mat[members_i,1:num_top_cat, drop=F], 2, mean);
			j_means=apply(norm_mat[members_j,1:num_top_cat, drop=F], 2, mean);
			greater_thans=(i_means > j_means);

			greater_than_col=rep(j, num_top_cat);
			greater_than_col[greater_thans]=i;

			# Plot bars
			# Extend the plot range so we have room to label the cluster with greater abundance
			plot_range=numeric(2);
			plot_range[1]=log_ratios_ranges[1]-diff(log_ratios_ranges)/10;
			plot_range[2]=log_ratios_ranges[2];

			par(mar=c(5,1,1,1));
			barpos=barplot(ratios, names.arg=plot_cat_names, 
				las=2, horiz=T, xlab="", col=greater_than_col,
				xlim=plot_range
			);
			midlines=diff(barpos)/2+barpos[-num_top_cat];
			abline(h=midlines, col="grey90");

			# Label greater thans for categories that contribute to a greater R^2
			for(cat in neg_ratios_ix){
				if(greater_thans[cat]){
					text(ratios[cat], barpos[cat], i, font=2, col=i, pos=2);
					text(ratios[cat], barpos[cat], i, font=1, col="black", pos=2);
				}else{
					text(ratios[cat], barpos[cat], j, font=2, col=j, pos=2);
					text(ratios[cat], barpos[cat], j, font=1, col="black", pos=2);
				}
			}

			# Label the 0 point
			abline(v=0);


			# Plot thumbnail MDS
			par(mar=c(3,1,0,1));
			plot(classic_mds_res[both_mem,], type="n", xaxt="n", yaxt="n", main="");
			for(grid_radius in grid_lines){
				draw.ellipse(0,0, a=grid_radius, b=grid_radius, border="grey90");
			}
			points(classic_mds_res[both_mem,], col=memberships[both_mem]);
			x_plot_range=par()$usr;
			axis(side=1, at=mean(x_plot_range[c(1,2)]), labels=sprintf("\n%i vs %i", i, j), tick=F)

			# Label number of clusters in the margins
			if(!(plot_count %% barplots_per_page)){
				mtext(paste("Num Clusters: ", num_cl), side=3, line=1.5, outer=T);
				mtext(paste("Cluster Differentiators:  Log10(R^2 Ratio) Plot"), side=3, line=0, outer=T);
			}

			plot_count=plot_count+1;

		}

	}

	# Compute cluster unifiers
	layout(barplot_layout);
	par(oma=c(.5,13,3,1));

	plot_count=0;
	pairs_names=names(ratios_list);
	for(i in 1:num_cl){

		ratios_by_cluster=matrix(0, nrow=num_cl-1, ncol=num_top_cat)
		greater_thans=matrix(F, nrow=num_cl-1, ncol=num_top_cat)

		colnames(ratios_by_cluster)=short_cat_names[1:num_top_cat];	
		colnames(greater_thans)=short_cat_names[1:num_top_cat];	

		mem_j=numeric(num_cl-1);

		ix=1;
		
		# Pull out the matching cluster comparisons
		for(j in 1:num_cl){

			if(i==j){
				next;
			}

			combin_name=paste(i, "#", j, sep="");

			if(any(combin_name==pairs_names)){
				ratios=ratios_list[[combin_name]];
			}else{
				combin_name=paste(j, "#", i, sep="");
				ratios=ratios_list[[combin_name]];
			}

			cat("Pulled: ", combin_name, "\n", sep="");
			ratios_by_cluster[ix,]=ratios;

			# Compute differences in abundance.
			members_i=names(memberships[memberships==i]);
			members_j=names(memberships[memberships==j]);

			i_means=apply(norm_mat[members_i,1:num_top_cat, drop=F], 2, mean);
			j_means=apply(norm_mat[members_j,1:num_top_cat, drop=F], 2, mean);
			greater_thans[ix,]=(i_means > j_means);

			mem_j[ix]=j;

			ix=ix+1;

		}

		cat("Focusing on: ", i, "\n");

		#cat("Greater thans:\n");
		#print(greater_thans);
		#cat("All Ratios:\n");
		#print(ratios_by_cluster);

		# Compute the mean ratios only cross the categories that have a greater abundance
		cluster_unifiers_matrix=matrix(0, nrow=num_cl-1, ncol=num_top_cat, 
			dimnames=list(1:(num_cl-1),short_cat_names[1:num_top_cat]));	

		for(cl_ix in 1:(num_cl-1)){
			for(cat_ix in 1:num_top_cat){
				if(ratios_by_cluster[cl_ix, cat_ix]<0 && greater_thans[cl_ix, cat_ix]){
					val=ratios_by_cluster[cl_ix, cat_ix];
				}else{
					val=NA;
				}
				cluster_unifiers_matrix[cl_ix, cat_ix]=val;
			}
		}

		#cat("Unifiers Matrix:\n");
		#print(cluster_unifiers_matrix);

		mean_wo_na=function(x){
			return(mean(x, na.rm=T));
		}
		mean_unifiers=apply(cluster_unifiers_matrix, 2, mean_wo_na);
		#cat("Kept:");
		#print(mean_unifiers);
		
		# Only plot the left most category labels
		if(!(plot_count %% barplots_per_page)){
			plot_cat_names=short_cat_names[1:num_top_cat];
		}else{
			plot_cat_names=rep("",num_top_cat);
		}

		par(mar=c(5,1,1,1));
		barpos=barplot(mean_unifiers, names.arg=plot_cat_names, 
			las=2, horiz=T, xlab="", col=i,
			xlim=plot_range
		);
		for(cat_ix in 1:num_top_cat){
			data_points=cluster_unifiers_matrix[, cat_ix];
			for(dp_ix in 1:length(data_points)){
				points(data_points[dp_ix], barpos[cat_ix], pch=23, bg=i, col="black");
			}
		}
		midlines=diff(barpos)/2+barpos[-num_top_cat];
		abline(h=midlines, col="grey90");

		# Label the 0 point
		abline(v=0);

		# Plot thumbnail MDS
		par(mar=c(3,1,0,1));
		cols=rep("grey50", num_samples);
		names(cols)=rownames(classic_mds_res);
		cols[members_i]=i;
		plot(classic_mds_res, type="n", xaxt="n", yaxt="n", main="");
		for(grid_radius in grid_lines){
			draw.ellipse(0,0, a=grid_radius, b=grid_radius, border="grey90");
		}
		points(classic_mds_res, col=cols);
		x_plot_range=par()$usr;
		axis(side=1, at=mean(x_plot_range[c(1,2)]), labels=paste("Cluster ", i, sep=""), tick=F)

		# Label number of clusters in the margins
		if(!(plot_count %% barplots_per_page)){
			mtext(paste("Num Clusters: ", num_cl), side=3, line=1.5, outer=T);
			mtext(paste("Cluster Unifiers:  Mean (Log10(R^2 Ratio) w/ Greater Abundances) Plot"), side=3, line=0, outer=T);
		}

		plot_count=plot_count+1;

		
	}
	

}


###############################################################################

cat("\nDone.\n")
dev.off();
warn=warnings();
if(length(warn)){
        print(warn);
}
q(status=0);
