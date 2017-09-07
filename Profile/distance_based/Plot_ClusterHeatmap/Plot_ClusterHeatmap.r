#!/usr/bin/env Rscript

###############################################################################

cat("\n\n");

library(MASS)
library('getopt');
library('vegan');
library('plotrix');
library('nnet');


DEF_DISTTYPE="euc";
DEF_NUM_TOP_CAT=35;
DEF_NUM_CLUS=10;
DEF_SPLIT_CHAR=";";
DEF_COLUMN=1;

params=c(
	"input_summary_table", "i", 1, "character",
	"input_factor_file", "f", 1, "character",
	"model_string", "m", 2, "character",
	"model_filename", "M", 2, "character",
	"output_filename_root", "o", 2, "character",
	"dist_type", "d", 2, "character",
	"num_clus", "k", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\n\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	-f <input factor file>\n",
	"	[-M <model variables list file>\n",
	"\n",
	"	[-o <output file root name, default is input file base name>]\n",
	"	[-d <euc/wrd/man/bray/horn/bin/gow/tyc/minkp5/minkp3, default =", DEF_DISTTYPE, ">]\n",
	"	[-k <max num of clusters to split into, default =", DEF_NUM_CLUS, ">\n",
	"\n",
	"This script will:\n",
	"	1.) Load metadata.\n",
	"	2.) Read in a summary table and compute a distance matrix.\n",
	"	3.) Cluster hierarchically.\n",
	"	4.) Generate a cluster heat map for each variable specified.\n",
	"\n");

if(!length(opt$input_summary_table)){
	cat(usage);
	q(status=-1);
}
InputFileName=opt$input_summary_table;
InputFactorFile=opt$input_factor_file;

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

#color_leaves=function(n){
#        if(is.leaf(n)){
#                leaf_attr=attributes(n);
#                leaf_name=leaf_attr$label;
#                attr(n, "nodePar") = c(leaf_attr$nodePar,
#					col=as.character(denfun.label_col_map[leaf_name]),
#					cex=denfun.label_scale*2.75,
#					pch=15
#                                        );
#        }
#        return(n);
#}

color_edges=function(den){
        edge_attr=attributes(den);
        attr(den, "edgePar") = c(edge_attr$edgePar, list(col=consensus_color));
        return(den);
}

plot_text=function(strings){

	orig_par=par(no.readonly=T);	

        par(family="Courier");
        par(oma=rep(.5,4));
        par(mar=rep(0,4));

        num_lines=length(strings);

        top=max(as.integer(num_lines), 52);

        plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                );

        text_size=max(.01, min(.8, .8 - .003*(num_lines-52)));
        #print(text_size);

        for(i in 1:num_lines){
                #cat(strings[i], "\n", sep="");
                strings[i]=gsub("\t", "", strings[i]);
                text(0, top-i, strings[i], pos=4, cex=text_size);
        }

	par(orig_par);

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

plot_heatmap=function(sample_names, factors){

	orig_par=par(no.readonly=T);

	# Create dummy variables out of factors
	model_string=paste("~ ", paste(colnames(factors), collapse="+"), " -1", sep="");
	print(model_string);
	factors=model.matrix(as.formula(model_string), data=factors);
	
        num_samples=nrow(factors);
        num_variables=ncol(factors);
	variable_names=colnames(factors);

        cat("Num Samples: ", num_samples, "\n");
        cat("Num Variables: ", num_variables, "\n");

	# Define number of different colors
        num_colors=50;
        color_arr=rev(rainbow(num_colors, start=0, end=4/6));

	# Compute ranges
	var_ranges=matrix(NA, nrow=num_variables, ncol=3, dimnames=list(variable_names, c("Min", "Mean", "Max")));
	for(i in 1:num_variables){
		val=factors[,i];
		var_ranges[i,"Min"]=min(val, na.rm=T);
		var_ranges[i,"Mean"]=mean(val, na.rm=T);
		var_ranges[i,"Max"]=max(val, na.rm=T);
	}

	cat("Variable Ranges:\n");
	print(var_ranges);

        remap=function(in_val, in_range, out_range){
                in_prop=(in_val-in_range[1])/(in_range[2]-in_range[1])
                out_val=in_prop*(out_range[2]-out_range[1])+out_range[1];
                return(out_val);
        }

	variable_names=colnames(factors);
	sample_names=rownames(factors);

	color_mat=matrix(0, nrow=num_samples, ncol=num_factors);
	for(i in 1:num_factors){
		color_mat[,i]=color_arr[round(remap(factors[,i], c(var_ranges[i, "Min"], var_ranges[i, "Max"]), c(1, num_colors)))];
	}
	print(color_mat);

        plot(0,0, type="n", xlim=c(0,num_samples), ylim=c(0,num_factors), xaxt="n", yaxt="n", bty="n", xlab="", ylab="");

        axis(side=2, at=seq(.5, num_variables-.5, 1), labels=variable_names, las=2);

        for(x in 1:num_samples){
                for(y in 1:num_factors){

                        rect(x-1, y-1, (x-1)+1, (y-1)+1, border=NA, col=color_mat[x, y]);

			#text_lab=sprintf("%1.1g", color_mat[x,y]);
			text_lab=signif(factors[x,y], 2);
			text(x-.5, y-.5, text_lab, srt=45, cex=.5, font=2);
		}
        }


}

###############################################################################

output_fname_root = paste(OutputFileRoot, ".", dist_type, sep="");
cat("\n");
cat("Input Summary Table Name: ", InputFileName, "\n", sep="");
cat("Input Factor File: ", InputFactorFile, "\n", sep="");
cat("Output Filename Root: ", output_fname_root, "\n", sep="");
cat("Distance Type: ", dist_type, "\n", sep="");
cat("Max Num clusters: ", max_clusters, "\n", sep="");
cat("\n");

cat("Loading summary table...\n");
counts_mat=load_summary_table(InputFileName);

# Normalize counts
cat("Normalizing counts...\n");
norm_mat=normalize(counts_mat);

category_names=colnames(norm_mat);
sample_names=rownames(norm_mat);
num_samples=nrow(norm_mat);
num_categories=ncol(norm_mat);

# Factor Files
cat("Loading factor file...\n");
all_factors=load_factor_file(InputFactorFile);
factor_samples=rownames(all_factors);
factor_names=colnames(all_factors);

###############################################################################

# Reconcile samples between groupings and summary table
shared_samples=sort(intersect(sample_names, factor_samples));
num_shared_samples=length(shared_samples);
cat("Num Shared Samples between Groupings/Summary Table: ", num_shared_samples, "\n");

norm_mat=norm_mat[shared_samples,, drop=F];
shared_factors=all_factors[shared_samples,, drop=F];
num_factors=ncol(shared_factors);

excl_st_samples=setdiff(sample_names, shared_samples);
excl_gr_samples=setdiff(factor_samples, shared_samples);

if(length(excl_st_samples)){
	cat("Samples exclusive to Summary Table:\n");
	print(excl_st_samples);
}

if(length(excl_gr_samples)){
	cat("Samples exclusive to Factor File:\n");
	print(excl_gr_samples);
}

###############################################################################

fact_per_inch=15;
samp_per_inch=10;

pdf(paste(OutputFileRoot, ".clhm.pdf", sep=""), height=5+num_factors/fact_per_inch, width=1+num_shared_samples/samp_per_inch);
par(oma=c(1,1,1,1));

###############################################################################

# Output run info
plot_text(c(
	paste("Input Summary Table: ", InputFileName, sep=""),
	paste("Input Factor File: ", InputFactorFile, sep=""),
	"",
	paste("Num Shared Samples: ", num_shared_samples, sep=""),
	"",
	paste("Distance Type: ", dist_type, sep="")
));

###############################################################################
# Precompute distances and perform clustering

# Compute full distances
cat("Computing distances...\n");
full_dist_mat=compute_dist(norm_mat, dist_type);

# Generate hierarchical clustering
cat("Clustering...\n");
hcl=hclust(full_dist_mat, method="ward.D2");

# Find height where cuts are made
cat("Finding Midpoints...\n");
cut_midpoints=numeric(max_clusters);
for(k in 2:max_clusters){
	cut_midpoints[k]=find_height_at_k(hcl, k);
}

# Convert to dendrogram
orig_dendr=as.dendrogram(hcl);

# Exactract names of leaves from left to right
lf_names=get_clstrd_leaf_names(orig_dendr);

# Assign cluster colors
palette_col=c("red", "green", "blue", "cyan", "magenta", "orange", "gray", "pink", "black", "purple", "brown", "aquamarine");
palette(palette_col);

###############################################################################
dend_height=5;
dend_width=10; #num_shared_samples;
hmap_height=15; #num_factors;
key_width=3;

layout_mat=matrix(c(
	rep(
		c(rep(1, dend_width), rep(2, key_width))
	, dend_height),
	
	rep(
		c(rep(3, dend_width), rep(4, key_width))
	, hmap_height)),
	
	byrow=T, ncol=dend_width+key_width);

print(layout_mat);
layout(layout_mat);

###############################################################################

# Define static variables used by dendro call back functions
denfun.label_scale=min(2,50/num_shared_samples);


par(oma=c(1,10,4,1));
for(num_cl in 1:1){
#for(num_cl in 2:max_clusters){

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
	sample_to_color_map=as.list(memberships);

	tweaked_dendro=dendrapply(orig_dendr, color_denfun_bySample);
	tweaked_dendro=dendrapply(tweaked_dendro, text_scale_denfun);

	par(mar=c(4,0,0,0));
	plot(tweaked_dendro, horiz=F);
	for(cl_ix in 1:num_cl){
		lab_size=2/ceiling(log10(cl_ix+1));
		axis(side=1, outer=F, at=grp_mids[cl_ix], labels=cl_ix, cex.axis=lab_size, col.ticks=cl_ix, 
			lend=1, lwd=10, padj=1, line=2);
	}

	abline(h=cut_midpoints[num_cl], col="red", lty=2);
	ranges=par()$usr;

	par(mar=c(0,0,0,0));
        plot(0,0, xlim=c(0,10), ylim=c(0,10), type="n",  xaxt="n", yaxt="n", bty="n", xlab="", ylab="");
	legend(0, 10, 
		fill=c("white", 1:num_cl), 
		legend=c("Cluster IDs:", as.character(1:num_cl)), 
		border=c("white", rep("black", num_cl)),
		text.font=c(2, rep(1, num_cl)),
		box.col="white", bg="white");
	mtext(paste("Distance Type: ", dist_type), side=3, line=0, outer=T);
	mtext(paste("Num Clusters: ", num_cl), side=3, line=1, outer=T);

	# Plot heatmap and key
	plot_heatmap(lf_names, shared_factors[lf_names,, drop=F]);
	#plot(0,0, type="n"); text(0,0, "1", cex=5);
	plot(0,0, type="n"); text(0,0, "2", cex=5);


}

###############################################################################

cat("\nDone.\n")
dev.off();
warn=warnings();
if(length(warn)){
        print(warn);
}
q(status=0);
