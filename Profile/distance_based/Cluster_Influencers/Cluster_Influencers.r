#!/usr/bin/env Rscript

###############################################################################

cat("\n\n");

library(MASS)
library('getopt');
library('vegan');


DEF_DISTTYPE="euc";
DEF_NUM_TOP_CAT=10;
DEF_NUM_CLUS=8;
DEF_SPLIT_CHAR=";";

params=c(
	"input_summary_table", "i", 1, "character",
	"output_filename_root", "o", 2, "character",
	"dist_type", "d", 2, "character",
	"num_top_cat", "p", 2, "numeric",
	"num_clus", "k", 2, "numeric",
	"split_char", "s", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\n\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	[-o <output file root name, default is input file base name>]\n",
	"	[-d <euc/wrd/man/bray/horn/bin/gow/tyc, default =", DEF_DISTTYPE, ">]\n",
	"	[-p <num of top categories to probe, default =", DEF_NUM_TOP_CAT, " >\n",
	"	[-k <num of clusters to split into, default =", DEF_NUM_CLUS, ">\n",
	"	[-s <split char for long category names, default =", DEF_SPLIT_CHAR, "\n",
	"\n",
	"This script will:\n",
	"	1.) Read in a summary table and compute a distance matrix.\n",
	"	2.) Cluster hierarchically.\n",
	"	3.) Iteratively cut tree until max_c clusters.\n",
	"	4.) For each pair of clusters, compute:\n",
	"		a.) computer SSB w/ all categories\n",
	"		b.) for each category, compute SSB w/o category of interest\n",
	"		c.) Compute SSB[i]/SSB[all] for top p taxa.\n",
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

if(!any(dist_type == c("wrd","man","bray","horn","bin","gow","euc","tyc"))){
	cat("Error: Specified distance type: ", dist_type, " not recognized.\n");
	quit(status=-1);
}

max_clusters=DEF_NUM_CLUS;
if(length(opt$num_clus)){
	max_clusters=opt$num_clus;
}

num_top_cat=DEF_NUM_TOP_CAT;
if(length(opt$num_top_cat)){
	num_top_cat=opt$num_top_cat;
}

SplitChar=DEF_SPLIT_CHAR;
if(length(opt$split_char)){
	SplitChar=opt$split_char;
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
                                        cex=0,
                                        lab.cex=label_scale);
        }
        return(n);
}

###############################################################################

output_fname_root = paste(OutputFileRoot, ".", dist_type, sep="");
cat("\n");
cat("Input Summary Table Name: ", InputFileName, "\n", sep="");
cat("Output Filename Root: ", output_fname_root, "\n", sep="");
cat("Distance Type: ", dist_type, "\n", sep="");
cat("Max Num clusters: ", max_clusters, "\n", sep="");
cat("Num Top categories to analyze: ", num_top_cat, "\n", sep="");
cat("\n");

cat("Loading summary table...\n");
counts_mat=load_summary_table(InputFileName);
#print(counts_mat);

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
}
colnames(norm_mat)=short_cat_names;
cat("Shorted Top Categories: \n");
print(short_cat_names[1:num_top_cat]);
cat("\n");

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

print(hcl);

orig_dendr=as.dendrogram(hcl);

pdf(paste(output_fname_root, ".cl_inf.pdf", sep=""), height=8.5, width=14);
palette_col=c("black", "red", "green", "blue", "cyan", "magenta", "orange", "gray");
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

barplot_layout=matrix(c(
	1,3,5,7,9,11,
	1,3,5,7,9,11,
	1,3,5,7,9,11,
	2,4,6,8,10,12), byrow=T, nrow=4, ncol=6);
barplots_per_page=6;

label_scale=min(2,50/num_samples);

# Begin pair-wise cluster analyses
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
	plot(tweaked_dendro);
	ranges=par()$usr;
	legend(ranges[1], ranges[4], fill=1:num_cl, legend=c(as.character(1:num_cl)), bty="n");
	mtext(paste("Num Clusters: ", num_cl), side=3, outer=T);

	# Generate MDS plots
	par(oma=c(0,0,2,0));
	par(mar=c(5.1,4.1,4.1,2.1));
	layout(mds_layout);
	plot(nonparm_mds_res, col=memberships, xlab="Dim 1", ylab="Dim 2", main="non-metric MDS");
	plot(classic_mds_res, col=memberships, xlab="Dim 1", ylab="Dim 2", main="classical MDS");
	par(mar=c(0,0,0,0));
	plot(0, type="n", xlab="", ylab="", main="", bty="n", xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1));
	legend(0,1, fill=1:num_cl, legend=c(as.character(1:num_cl)), bty="n", cex=2);
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
			sub_mat=norm_mat[c(members_i, members_j),];
			results=analyze_cluster_pairs(sub_mat, members_i, members_j, num_top_cat, dist_mat_list);

			lograt=log10(results$ratios[,"rsqrd"]);
			ratios_list[[ paste(i, "#", j, sep="")]]=lograt;
			log_ratios_ranges=range(log_ratios_ranges, range(lograt));

		}

	}

	# Generate Plots
	layout(barplot_layout);
	par(oma=c(.5,10,2,1));
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
			i_means=apply(norm_mat[members_i,1:num_top_cat], 2, mean);
			j_means=apply(norm_mat[members_j,1:num_top_cat], 2, mean);
			greater_thans=(i_means > j_means);

			# Plot bars
			# Extend the plot range so we have room to label the cluster with greater abundance
			plot_range=numeric(2);
			plot_range[1]=log_ratios_ranges[1]-diff(log_ratios_ranges)/10;
			plot_range[2]=log_ratios_ranges[2];

			par(mar=c(5,1,1,1));
			barpos=barplot(ratios, names.arg=plot_cat_names, 
				las=2, horiz=T, xlab="", col=pos_contrib,
				xlim=plot_range
			);

			# Label greater thans for categories that contribute to a greater R^2
			for(cat in neg_ratios_ix){
				if(greater_thans[cat]){
					text(ratios[cat], barpos[cat], i, font=2, col=i, pos=2);
				}else{
					text(ratios[cat], barpos[cat], j, font=2, col=j, pos=2);
				}
			}

			# Label the 0 point
			abline(v=0);


			# Plot thumbnail MDS
			par(mar=c(3,1,0,1));
			plot(classic_mds_res[both_mem,], col=memberships[both_mem], xaxt="n", yaxt="n", 
				xlab="",
				ylab="", main="");
			x_plot_range=par()$usr;
			axis(side=1, at=mean(x_plot_range[c(1,2)]), labels=sprintf("\n%i vs %i", i, j), tick=F)

			# Label number of clusters in the margins
			if(!(plot_count %% barplots_per_page)){
				mtext(paste("Num Clusters: ", num_cl), side=3, outer=T);
			}

			plot_count=plot_count+1;

		}

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
