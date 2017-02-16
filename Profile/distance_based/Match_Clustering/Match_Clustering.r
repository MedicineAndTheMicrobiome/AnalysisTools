#!/usr/bin/env Rscript

###############################################################################

cat("\n\n");

library(MASS)
library('getopt');
library('vegan');
library('plotrix');


DEF_DISTTYPE="euc";
DEF_NUM_TOP_CAT=35;
DEF_NUM_CLUS=10;
DEF_SPLIT_CHAR=";";
DEF_COLUMN=1;

params=c(
	"input_summary_table", "i", 1, "character",
	"input_grouping_file", "r", 1, "character",
	"output_filename_root", "o", 2, "character",
	"grouping_column", "c", 2, "numeric",
	"dist_type", "d", 2, "character",
	"num_clus", "k", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\n\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	-r <input grouping/metadata .tsv file>\n",
	"	[-o <output file root name, default is input file base name>]\n",
	"	[-c <column with grouping information (column headers expected), default=", DEF_COLUMN, ">]\n",
	"	[-d <euc/wrd/man/bray/horn/bin/gow/tyc/minkp5/minkp3, default =", DEF_DISTTYPE, ">]\n",
	"	[-k <max num of clusters to split into, default =", DEF_NUM_CLUS, ">\n",
	"\n",
	"This script will:\n",
	"	1.) Load groups (metadata).\n",
	"	2.) Read in a summary table and compute a distance matrix.\n",
	"	3.) Cluster hierarchically.\n",
	"	4.) Iteratively cut tree until max_c clusters.\n",
	"	5.) For each group of clusters, compute:\n",
	"		a.) pseudo F-stat\n",
	"		b.) Bootstrap Null distribution\n",
	"		c.) Bootstrap based on group information\n",
	"		d.) Compute p-value, where null hypothesis is clusters are compatible\n",
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
InputGroupingFile=opt$input_grouping_file;

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

if(!any(dist_type == c("wrd","man","bray","horn","bin","gow","euc","tyc","minkp3","minkp5"))){
	cat("Error: Specified distance type: ", dist_type, " not recognized.\n");
	quit(status=-1);
}

max_clusters=DEF_NUM_CLUS;
if(length(opt$num_clus)){
	max_clusters=opt$num_clus;
}

GroupingColumn=DEF_COLUMN;
if(length(opt$grouping_column)){
	GroupingColumn=opt$grouping_column;
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

load_mapping=function(fn, target_col=1){
	inmat=as.matrix(read.delim(fn, sep="\t", header=TRUE, row.names=1, check.names=F, comment.char="", quote=""));
	map_val=inmat[,target_col, drop=F];
	map_val=map_val[!is.na(map_val),, drop=F];
	return(map_val);
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

shannon_entropy=function(counts, num_bs=1000){
	tot=sum(counts);
	prob=counts/tot;

	resamples=t(rmultinom(num_bs, tot, prob));

	ent=function(cts){
		cts=cts[cts>0];
		tot=sum(cts);
		pr=cts/tot;
		return(-sum(pr*log2(pr)));
	}

	se_arr=apply(resamples, 1, ent);

	ci=quantile(se_arr, c(.5, .025, .975));
	return(ci);
	
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

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, counts=F){
        num_row=nrow(mat);
        num_col=ncol(mat);

        cat("Num Rows: ", num_row, "\n");
        cat("Num Cols: ", num_col, "\n");

        mat=mat[rev(1:num_row),, drop=F];

        num_colors=50;
        color_arr=rainbow(num_colors, start=0, end=4/6);
        if(high_is_hot){
                color_arr=rev(color_arr);
        }

        remap=function(in_val, in_range, out_range){
                in_prop=(in_val-in_range[1])/(in_range[2]-in_range[1])
                out_val=in_prop*(out_range[2]-out_range[1])+out_range[1];
                return(out_val);
        }

        if(is.na(plot_min)){
                plot_min=min(mat);
        }
        if(is.na(plot_max)){
                plot_max=max(mat);
        }
        cat("Plot min/max: ", plot_min, "/", plot_max, "\n");

        plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), xaxt="n", yaxt="n", bty="n", xlab="", ylab="", main=title);

        # x-axis
        axis(side=1, at=seq(.5, num_col-.5, 1), labels=colnames(mat), las=2);
        axis(side=2, at=seq(.5, num_row-.5, 1), labels=rownames(mat), las=2);

        if(log_col){
                plot_min=log10(plot_min+.0125);
                plot_max=log10(plot_max+.0125);
        }

        for(x in 1:num_col){
                for(y in 1:num_row){

                        if(log_col){
                                col_val=log10(mat[y,x]+.0125);
                        }else{
                                col_val=mat[y,x];
                        }

                        remap_val=remap(col_val, c(plot_min, plot_max), c(1, num_colors));
                        col_ix=ceiling(remap_val);

                        rect(x-1, y-1, (x-1)+1, (y-1)+1, border=NA, col=color_arr[col_ix]);

			if(counts){
				text_lab=sprintf("%i", mat[y,x]);
			}else{
				text_lab=sprintf("%0.4f", mat[y,x]);
			}
			text(x-.5, y-.5, text_lab, srt=45, cex=2, font=2);
                }
        }

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

###############################################################################

output_fname_root = paste(OutputFileRoot, ".", dist_type, sep="");
cat("\n");
cat("Input Summary Table Name: ", InputFileName, "\n", sep="");
cat("Input Grouping File: ", InputGroupingFile, "\n", sep="");
cat("  Grouping Column: ", GroupingColumn, "\n", sep="");
cat("Output Filename Root: ", output_fname_root, "\n", sep="");
cat("Distance Type: ", dist_type, "\n", sep="");
cat("Max Num clusters: ", max_clusters, "\n", sep="");
cat("\n");

cat("Loading summary table...\n");
counts_mat=load_summary_table(InputFileName);
#print(counts_mat);

# Normalize counts
cat("Normalizing counts...\n");
norm_mat=normalize(counts_mat);
#print(norm_mat);

category_names=colnames(norm_mat);
sample_names=rownames(norm_mat);
num_samples=nrow(norm_mat);
num_categories=ncol(norm_mat);

# Group Information
grouping_info=load_mapping(InputGroupingFile, GroupingColumn);
group_name=colnames(grouping_info);
grouping_map=grouping_info[,1];
cat("Group Name: ", group_name, "\n");
print(grouping_map);
num_group_samples=length(grouping_map);
groups=sort(unique(grouping_map));
group_samples=names(grouping_map);
num_groups=length(groups);
cat("Groups:\n");
print(groups);
cat("Num Groups: ", num_groups, "\n");
group_id=1:num_groups;
names(group_id)=groups;
cat("Group IDs:\n");
print(group_id);
group_color_palette=c("aquamarine", "coral", "cornflowerblue", "gold", "deeppink", "limegreen", "chocolate", "cornsilk4", "blueviolet");
group_color_palette=group_color_palette[1:num_groups];
group_colors=character(num_group_samples);

for(i in group_id){
	group_colors[grouping_map==groups[i]]=group_color_palette[i];
}
names(group_colors)=group_samples;

# Reconcile samples between groupings and summary table
shared_samples=intersect(sample_names, group_samples);
num_shared_samples=length(shared_samples);
cat("Num Shared Samples between Groupings/Summary Table: ", num_shared_samples, "\n");

norm_mat=norm_mat[shared_samples,, drop=F];
grouping_map=grouping_map[shared_samples];

excl_st_samples=setdiff(sample_names, shared_samples);
excl_gr_samples=setdiff(group_samples, shared_samples);

if(length(excl_st_samples)){
	cat("Samples exclusive to Summary Table:\n");
	print(excl_st_samples);
}

if(length(excl_gr_samples)){
	cat("Samples exclusive to Groupings:\n");
	print(excl_gr_samples);
}

###############################################################################

# Compute full distances
cat("Computing distances...\n");
full_dist_mat=compute_dist(norm_mat, dist_type);

###############################################################################

# Generate hierarchical clustering
hcl=hclust(full_dist_mat, method="ward.D2");

# Find height where cuts are made
cut_midpoints=numeric(max_clusters);
for(k in 2:max_clusters){
	cut_midpoints[k]=find_height_at_k(hcl, k);
}

# Convert to dendrogram
orig_dendr=as.dendrogram(hcl);

# Exactract names of leaves from left to right
lf_names=get_clstrd_leaf_names(orig_dendr);

# Open output PDF file
pdf(paste(output_fname_root, ".", group_name, ".cl_match.pdf", sep=""), height=8.5, width=14);

# Assign cluster colors
palette_col=c("red", "green", "blue", "cyan", "magenta", "orange", "gray", "pink", "black", "purple", "brown", "aquamarine");
palette(palette_col);

# Compute ISO and classical MDS
nonparm_mds_res=matrix(0,nrow=num_shared_samples,ncol=2);
classic_mds_res=matrix(0,nrow=num_shared_samples,ncol=2);

# -- ISO MDS
imds_res=isoMDS(full_dist_mat);
#print(imds_res$points);
nonparm_mds_res[,1]=imds_res$points[,1]
nonparm_mds_res[,2]=imds_res$points[,2]

# -- Classical MDS
classic_mds_res=cmdscale(full_dist_mat);

###############################################################################
# Define layout for multiple plot pages

mds_layout=matrix(c(
	1,1,1,1,2,2,2,2,3), byrow=T, nrow=1, ncol=9);

cont_tab_layout=matrix(c(
	1,1,1,1,1,1,2,
	1,1,1,1,1,1,2,
	1,1,1,1,1,1,2,
	3,3,3,3,3,3,4),
	byrow=T, nrow=4, ncol=7);

indiv_pval_layout=matrix(c(
	1,1,1,1,1,1,2,
	1,1,1,1,1,1,2,
	1,1,1,1,1,1,2,
	1,1,1,1,1,1,2,
	3,3,3,3,3,3,4),
	byrow=T, nrow=5, ncol=7);

# Define static variables used by dendro call back functions
denfun.label_scale=min(2,55/num_shared_samples);
denfun.label_col_map=group_colors;

# Begin per cluster count analyses
pvalues=numeric(max_clusters-1);
i=1;
for(num_cl in 2:max_clusters){

	# Cut tree at target number of clusters
	cat("Cutting for ", num_cl, " clusters...\n", sep="");
	memberships=cutree(hcl, k=num_cl);
	grp_mids=get_middle_of_groups(lf_names, memberships);

	# Reorder cluster assignments to match dendrogram left/right
	plot_order=order(grp_mids);
	mem_tmp=numeric(num_shared_samples);
	for(gr_ix in 1:num_cl){
		old_id=(memberships==plot_order[gr_ix]);
		mem_tmp[old_id]=gr_ix;
	}
	names(mem_tmp)=names(memberships);
	memberships=mem_tmp;
	grp_mids=grp_mids[plot_order];

	# Plot Dendrogram
	par(oma=c(4,0,5,0));
	par(mar=c(5.1,4.1,4.1,2.1));
	par(mfrow=c(1,1));
	sample_to_color_map=as.list(memberships);

	tweaked_dendro=dendrapply(orig_dendr, color_denfun_bySample);
	tweaked_dendro=dendrapply(tweaked_dendro, text_scale_denfun);
	tweaked_dendro=dendrapply(tweaked_dendro, color_leaves);

	plot(tweaked_dendro, horiz=F);
	for(cl_ix in 1:num_cl){
		lab_size=3/ceiling(log10(cl_ix+1));
		axis(side=1, outer=T, at=grp_mids[cl_ix], labels=cl_ix, cex.axis=lab_size, col.ticks=cl_ix, 
			lend=1, lwd=10, padj=1, line=-3);
	}

	abline(h=cut_midpoints[num_cl], col="red", lty=2);
	ranges=par()$usr;
	legend(ranges[1]+(ranges[2]-ranges[1])/100, ranges[4], 
		fill=c("white", 1:num_cl, "white", "white", group_color_palette), 
		legend=c("Cluster IDs:", as.character(1:num_cl), "", "Groups:", groups), 
		border=c("white", rep("black", num_cl), "white", "white", rep("black", num_groups)),
		text.font=c(2, rep(1, num_cl), 1, 2, rep(1, num_groups)),
		box.col="white", bg="white");
	mtext(paste("Distance Type: ", dist_type), side=3, line=0, outer=T);
	mtext(paste("Num Clusters: ", num_cl), side=3, line=1, outer=T);

	# Generate MDS plots
	par(oma=c(0,0,4,0));
	par(mar=c(5.1,4.1,4.1,2.1));
	layout(mds_layout);
	plot(nonparm_mds_res, col=memberships, xlab="Dim 1", ylab="Dim 2", main="non-metric MDS");
	plot(classic_mds_res, type="n", col=memberships, xlab="Dim 1", ylab="Dim 2", main="classical MDS");
	points(classic_mds_res, col=memberships);
	
	# MDS Legend
	par(mar=c(0,0,0,0));
	plot(0, type="n", xlab="", ylab="", main="", bty="n", xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1));
	legend(0,1, fill=1:num_cl, legend=c(as.character(1:num_cl)), bty="n", cex=2);

	mtext(paste("Distance Type: ", dist_type, sep=""), side=3, outer=T, line=.5);
	mtext(paste("Num Clusters: ", num_cl, sep=""), side=3, outer=T, line=2);

	# Generate contigency table
	cont_tab=matrix(0, nrow=num_groups, ncol=num_cl, dimnames=list(groups, 1:num_cl));
	for(sample_id in shared_samples){
		#cat(paste(sample_id, memberships[sample_id], grouping_map[sample_id], sep="/"), "\n");
		cell_cont=cont_tab[grouping_map[sample_id], memberships[sample_id]];
		cont_tab[grouping_map[sample_id], memberships[sample_id]]=cell_cont+1;
	}
	cat("Contingency Table:\n");
	print(cont_tab);
	
	#-----------------------------------------------------------------------------
	# Compute fisher exact for the ngrps x nclusters table
	#ft=fisher.test(cont_tab, workspace=200000*10000);
	ft=fisher.test(cont_tab, simulate.p.value=T, B=20000*max_clusters)
	pvalues[i]=ft$p.value;

	# Plot contigency table 
	par(oma=c(0,0,5,0));
	layout(cont_tab_layout);
	par(mar=c(2,15,1,2));
	paint_matrix(cont_tab, title=paste("Contingency Table: p-value = ", sprintf("%.3g", pvalues[i]), sep=""), plot_min=0, counts=T);
	title(ylab=group_name, line=8, cex.lab=3);
	par(mar=c(2,0,1,2));
	paint_matrix(matrix(apply(cont_tab, 1, sum), ncol=1, dimnames=list(dimnames(cont_tab)[[1]],"total")), plot_min=0, counts=T);
	par(mar=c(5,15,1,2));
	paint_matrix(matrix(apply(cont_tab, 2, sum), nrow=1, dimnames=list("total",dimnames(cont_tab)[[2]])), plot_min=0, counts=T);
	title(xlab="Cluster Number")

	mtext(paste("Distance Type: ", dist_type, sep=""), side=3, outer=T, line=.5);
	mtext(paste("Num Clusters: ", num_cl, sep=""), side=3, outer=T, line=2);

	# Compute p-values for each cluster
	cl_pval=matrix(0, nrow=1, ncol=num_cl, dimnames=list("Fisher Exact 2-way p-values", 1:num_cl));
	for(cl_ix in 1:num_cl){
		# Extract 2x2 from complete contigency table to test single cluster
		two_by_two=matrix(0, nrow=num_groups, ncol=2);
		two_by_two[,1]=cont_tab[,cl_ix, drop=F];
		two_by_two[,2]=apply(cont_tab[,-cl_ix, drop=F], 1, sum);
		#res=fisher.test(two_by_two);
		res=fisher.test(two_by_two, simulate.p.value=T, B=20000*max_clusters);
		cl_pval[1,cl_ix]=res$p.val;
	}
	
	#-----------------------------------------------------------------------------
	# Compute probabilities across each group
	group_sums=apply(cont_tab, 1, sum);
	norm_tab=matrix(0, nrow=num_groups, ncol=num_cl, dimnames=list(groups, 1:num_cl));
	shan_ent_mat=matrix(0, nrow=num_groups, ncol=3, dimnames=list(groups, c("Median", "95 LB", "95 UB")));
	for(grp_ix in 1:num_groups){
		norm_tab[grp_ix,]=cont_tab[grp_ix,]/group_sums[grp_ix];	
		shan_ent_mat[grp_ix,]=shannon_entropy(cont_tab[grp_ix,]);
	}

	# Plot probabilities across each group
	par(oma=c(0,0,5,0));
	par(mar=c(5,15,1,2));
	layout(indiv_pval_layout);

	paint_matrix(norm_tab, title="Normalized by Group Size", plot_min=0, counts=F);
	title(ylab=group_name, line=8, cex.lab=3);

	par(mar=c(5,0,1,2));
	paint_matrix(shan_ent_mat, title="Shan. Ent.", counts=F);

	par(mar=c(5,15,1,2));
	paint_matrix(cl_pval, plot_min=0, plot_max=1, counts=F, log_col=F, high_is_hot=F);
	title(xlab="Cluster Number");
	mtext(paste("Distance Type: ", dist_type, sep=""), side=3, outer=T, line=.5);
	mtext(paste("Num Clusters: ", num_cl, sep=""), side=3, outer=T, line=2);

	i=i+1;	
	
}

# Generate p-value vs. cluster count plot
log_pval=log10(pvalues);
min_pval_ix=which(min(log_pval)==log_pval);

pval_ticks=c(.1, .05, .025, .0125);
log_pval_ticks=log10(pval_ticks);

pval_range=range(c(log_pval, log_pval_ticks));

par(oma=c(0,0,4,0));
par(mar=c(5,5,5,1));
plot(cbind(2:max_clusters, log_pval), type="b", 
	ylim=pval_range,
	xaxt="n",
	xlab="Num Clusters", ylab="Log10(p-value)", main="Log10(p-value) vs. Num Clusters");

# Draw bulls eye for min pvalue
abline(h=log_pval_ticks, col="blue", lty=3);
axis(side=4, at=log_pval_ticks, label=pval_ticks);
points(min_pval_ix+1, log_pval[min_pval_ix], pch=1, col="red", cex=1, lwd=2);
points(min_pval_ix+1, log_pval[min_pval_ix], pch=1, col="red", cex=2.5, lwd=2);
points(min_pval_ix+1, log_pval[min_pval_ix], pch=1, col="red", cex=4, lwd=2);
points(min_pval_ix+1, log_pval[min_pval_ix], pch=1, col="red", cex=5.5, lwd=2);
mtext(paste("Distance Type: ", dist_type, sep=""), side=3, outer=T, line=1);

axis(side=1, at=2:max_clusters, labels=2:max_clusters);

###############################################################################

cat("\nDone.\n")
dev.off();
warn=warnings();
if(length(warn)){
        print(warn);
}
q(status=0);
