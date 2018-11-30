#!/usr/bin/env Rscript

###############################################################################

cat("\n\n");

library(MASS)
library('getopt');
library('vegan');
library('plotrix');


DEF_DISTTYPE="euc";
DEF_NUM_TOP_CAT=35;
DEF_NUM_CLUS=-1;
DEF_SPLIT_CHAR=";";


params=c(
	"input_summary_table_A", "a", 1, "character",
	"input_summary_table_B", "b", 1, "character",
	"mapping_file", "m", 1, "character",
	"output_filename_root", "o", 1, "character",
	"dist_type", "d", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\n\nUsage:\n", script_name, "\n",
	"	-a <input summary_table.tsv file A>\n",
	"	-b <input summary_table.tsv file B>\n",
	"	-m <mapping file, from A to B Identifiers>\n",
	"	-o <output file root name>\n",
	"\n",
	"	Options:\n",
	"	[-d <euc/wrd/man/bray/horn/bin/gow/tyc/minkp5/minkp3, default =", DEF_DISTTYPE, ">]\n",
	"\n",
	"This script will compute distance matrices for the two summary tables independently, then\n",
	"compare them using Mantel's statistic.  The summary tables do not have to have matching\n",
	"categories.\n",
	"\n\n",
	"The mapping file must be specified so that the sample IDs can be paired up.\n",
	"The first row of each column in the mapping file will be used in the analyses.\n",
	"Eg.:\n",
	"\n",
	"<sample group A name>\\t<sample group B name>\n",
	"<id.A.1>\\t<id.B.1>\n",
	"<id.A.2>\\t<id.B.2>\n",
	"<id.A.3>\\t<id.B.3>\n",
	"...\n",
	"<id.A.N>\\t<id.B.N>\n",
	
	"\n\n",
	"For the distance types:\n",
	" minkp5 is the minkowski with p=.5, i.e. sum((x_i-y_i)^1/2)^2\n",
	" minkp3 is the minkowski with p=.3, i.e. sum((x_i-y_i)^1/3)^3\n",
	"\n");

if(
	!length(opt$input_summary_table_A) || 
	!length(opt$input_summary_table_B) || 
	!length(opt$mapping_file) || 
	!length(opt$output_filename_root) 
){
	cat(usage);
	q(status=-1);
}

InputSumTabA=opt$input_summary_table_A;
InputSumTabB=opt$input_summary_table_B;
MappingFile=opt$mapping_file;
OutputFileRoot=opt$output_filename_root;

DistType=DEF_DISTTYPE;
if(length(opt$dist_type)){
	DistType=opt$dist_type;
}


if(!any(DistType== c("wrd","man","bray","horn","bin","gow","euc","tyc","minkp3","minkp5"))){
	cat("Error: Specified distance type: ", DistType, " not recognized.\n");
	quit(status=-1);
}

###############################################################################

cat("Input Summary Table A:", InputSumTabA, "\n");
cat("Input Summary Table B:", InputSumTabB, "\n");
cat("Mapping File         :", MappingFile, "\n");
cat("Output File          :", OutputFileRoot, "\n");
cat("Distance Type        :", DistType, "\n");



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

load_mapping_file=function(mp_fname, keep_a_ids, keep_b_ids){

	num_keep_a=length(keep_a_ids);
	num_keep_b=length(keep_b_ids);
	cat("Num A's IDs to keep: ", num_keep_a, "\n");
	cat("Num B's IDs to keep: ", num_keep_b, "\n");

	inmat=as.matrix(read.delim(mp_fname, sep="\t", header=TRUE, check.names=F, comment.char="", quote=""));

	# Keep Entry if record is in both lists
	keep_ix=c();
	orig_mat_rows=nrow(inmat);
	cat("Number of Mapping Entries Read: ", orig_mat_rows, "\n");
	for(i in 1:orig_mat_rows){
		if(any(inmat[i,1]==keep_a_ids) && any(inmat[i,2]==keep_b_ids)){
			keep_ix=c(keep_ix, i);
		}
	}
	inmat=inmat[keep_ix,];
	num_kept_matrows=nrow(inmat);
	cat("Number of Mapping Entries Kept: ", num_kept_matrows, "\n");

	mapping=as.list(x=inmat[,1]);
	names(mapping)=inmat[,2];

	coln=colnames(inmat);
	map_info=list();
	map_info[["map"]]=mapping;
	map_info[["a"]]=coln[1];
	map_info[["b"]]=coln[2];
	map_info[["a_id"]]=inmat[,1];
	map_info[["b_id"]]=inmat[,2];

	return(map_info);	
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

###############################################################################

output_fname_root = paste(OutputFileRoot, ".", DistType, sep="");

cat("\n");
cat("Loading summary table A:", InputSumTabA, "\n");
counts_mat_A=load_summary_table(InputSumTabA);
cat("Loading summary table B:", InputSumTabB, "\n");
counts_mat_B=load_summary_table(InputSumTabB);

# Reconcile summary table IDs through mapping file
samples_stA=rownames(counts_mat_A);
samples_stB=rownames(counts_mat_B);

cat("Loading Mapping file:", MappingFile, "\n");
map_info=load_mapping_file(MappingFile, samples_stA, samples_stB);

cat("Removing samples without complete mappings...\n");
counts_mat_A=counts_mat_A[map_info[["a_id"]],];
counts_mat_B=counts_mat_B[map_info[["b_id"]],];

###############################################################################

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
max_clusters=ceiling(log(num_samples, 2));
cat("Max Clusters to compute: ", max_clusters, "\n");
cut_midpoints=numeric(max_clusters);
for(k in 2:max_clusters){
        cut_midpoints[k]=find_height_at_k(hcl, k);
}

orig_dendr=as.dendrogram(hcl);

# Extract names of leaves from left to right
lf_names=get_clstrd_leaf_names(orig_dendr);

pdf(paste(output_fname_root, ".ch_stop.pdf", sep=""), height=8.5, width=14);

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


###############################################################################

cat("\nDone.\n")
dev.off();
warn=warnings();
if(length(warn)){
        print(warn);
}
q(status=0);
