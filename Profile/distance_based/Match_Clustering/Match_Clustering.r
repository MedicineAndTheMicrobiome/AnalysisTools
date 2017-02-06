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
	return(inmat[,target_col]);
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

        mat=mat[rev(1:num_row),];

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

        #par(oma=c(12,12,6,6));
        par(mar=c(12, 12, 8, 6));
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
			text(x-.5, y-.5, text_lab, srt=45);
                }
        }

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
grouping_map=load_mapping(InputGroupingFile, GroupingColumn);
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
group_color_palette=c("aquamarine", "coral", "cornflowerblue", "gold", "deeppink", "limegreen");
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

hcl=hclust(full_dist_mat, method="ward.D2");

print(hcl);

orig_dendr=as.dendrogram(hcl);

pdf(paste(output_fname_root, ".cl_match.pdf", sep=""), height=8.5, width=14);

palette_col=c("red", "green", "blue", "cyan", "magenta", "orange", "gray", "pink", "black", "purple");
palette(palette_col);

# Compute ISO and classical MDS
nonparm_mds_res=matrix(0,nrow=num_samples,ncol=2);
classic_mds_res=matrix(0,nrow=num_samples,ncol=2);

imds_res=isoMDS(full_dist_mat);
nonparm_mds_res[,1]=imds_res$points[,1]
nonparm_mds_res[,2]=imds_res$points[,2]

classic_mds_res=cmdscale(full_dist_mat);

mds_layout=matrix(c(
	1,1,1,1,2,2,2,2,3), byrow=T, nrow=1, ncol=9);

denfun.label_scale=min(2,55/num_samples);
denfun.label_col_map=group_colors;

# Begin per cluster count analyses
pvalues=numeric(max_clusters-1);
i=1;
for(num_cl in 2:max_clusters){

	cat("Cutting for ", num_cl, " clusters...\n", sep="");
	memberships=cutree(hcl, k=num_cl);

	# Plot Dendrogram
	par(oma=c(0,0,2,0));
	par(mar=c(5.1,4.1,4.1,2.1));
	par(mfrow=c(1,1));
	sample_to_color_map=as.list(memberships);

	tweaked_dendro=dendrapply(orig_dendr, color_denfun_bySample);
	tweaked_dendro=dendrapply(tweaked_dendro, text_scale_denfun);
	tweaked_dendro=dendrapply(tweaked_dendro, color_leaves);

	plot(tweaked_dendro, horiz=F);
	ranges=par()$usr;
	legend(ranges[1], ranges[4], 
		fill=c("white", 1:num_cl, "white", "white", group_color_palette), 
		legend=c("Cluster IDs:", as.character(1:num_cl), "", "Groups:", groups), 
		border=c("white", rep("black", num_cl), "white", "white", rep("black", num_groups)),
		text.font=c(2, rep(1, num_cl), 1, 2, rep(1, num_groups)),
		bty="n");
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

	mtext(paste("Distance Type: ", dist_type, sep=""), side=3, outer=T, line=0);
	mtext(paste("Num Clusters: ", num_cl, sep=""), side=3, outer=T, line=1);

	# Generate contigency table
	cont_tab=matrix(0, nrow=num_groups, ncol=num_cl, dimnames=list(groups, 1:num_cl));
	for(sample_id in shared_samples){
		#cat(paste(sample_id, memberships[sample_id], grouping_map[sample_id], sep="/"), "\n");
		cell_cont=cont_tab[grouping_map[sample_id], memberships[sample_id]];
		cont_tab[grouping_map[sample_id], memberships[sample_id]]=cell_cont+1;
	}
	
	print(cont_tab);
	print(cont_tab/num_shared_samples);
	ft=fisher.test(cont_tab, workspace=200000*10000);
	pvalues[i]=ft$p.value;

	par(oma=c(0,0,1,0));
	par(mfrow=c(1,1));
	paint_matrix(cont_tab, title=paste("Contingency Table: p-value = ", sprintf("%.3g", pvalues[i]), sep=""), plot_min=0, counts=T);
	title(xlab="Cluster Number", ylab="Grouping");
	mtext(paste("Distance Type: ", dist_type, sep=""), side=3, outer=T, 1);
	mtext(paste("Num Clusters: ", num_cl, sep=""), side=3, outer=T, line=2);

	group_sums=apply(cont_tab, 1, sum);
	norm_tab=matrix(0, nrow=num_groups, ncol=num_cl, dimnames=list(groups, 1:num_cl));
	for(grp_ix in 1:num_groups){
		norm_tab[grp_ix,]=cont_tab[grp_ix,]/group_sums[grp_ix];	
	}

	par(oma=c(0,0,1,0));
	paint_matrix(norm_tab, title="Normalized by Group Size", plot_min=0, counts=F);
	title(xlab="Cluster Number", ylab="Grouping");
	mtext(paste("Distance Type: ", dist_type, sep=""), side=3, outer=T, 1);
	mtext(paste("Num Clusters: ", num_cl, sep=""), side=3, outer=T, line=2);

	i=i+1;	
	
}

par(oma=c(2,2,1,1));
par(mar=c(4,4,4,4));
print(2:max_clusters);
print(pvalues);

log_pval=log(pvalues);
min_pval_ix=which(min(log_pval)==log_pval);
par(oma=c(0,0,1,0));
plot(cbind(2:max_clusters, log_pval), type="b", 
	xaxt="n",
	xlab="Num Clusters", ylab="Log10(p-value)", main="Log10(p-value) vs. Num Clusters");
points(min_pval_ix+1, log_pval[min_pval_ix], pch=1, col="red", cex=1, lwd=2);
points(min_pval_ix+1, log_pval[min_pval_ix], pch=1, col="red", cex=2.5, lwd=2);
points(min_pval_ix+1, log_pval[min_pval_ix], pch=1, col="red", cex=4, lwd=2);
points(min_pval_ix+1, log_pval[min_pval_ix], pch=1, col="red", cex=5.5, lwd=2);
mtext(paste("Distance Type: ", dist_type, sep=""), side=3, outer=T, 1);

axis(side=1, at=2:max_clusters, labels=2:max_clusters);




###############################################################################

cat("\nDone.\n")
dev.off();
warn=warnings();
if(length(warn)){
        print(warn);
}
q(status=0);
