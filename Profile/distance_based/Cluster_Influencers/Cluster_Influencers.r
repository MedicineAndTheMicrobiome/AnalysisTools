#!/usr/bin/env Rscript

###############################################################################

cat("\n\n");

library(MASS)
library('getopt');
library('vegan');


DEF_DISTTYPE="euc";
DEF_NUM_TOP_CAT=10;
DEF_NUM_CLUS=8;

params=c(
	"input_summary_table", "i", 1, "character",
	"output_filename_root", "o", 2, "character",
	"dist_type", "d", 2, "character",
	"num_top_cat", "p", 2, "numeric",
	"num_clus", "k", 2, "numeric"
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

compute_sum_of_squares=function(dist_mat, members_a, members_b){

	# Compute SST (SS Total)
	# Compute SSW (SS Error)
	# Compute SSB (SS Treatment)

	# Fstat = 

}

###############################################################################

analyze_cluster_pairs=function(st, members_a, members_b, num_cat_to_probe, dist_type){

	num_mem_a=length(members_a);
	num_mem_b=length(members_b);

	cat("Num members in clusters: ", num_mem_a, " vs. ", num_mem_b, "\n");

	cat("\nMembers A:\n");
	print(members_a);
	cat("\nMembers B:\n");
	print(members_b);
	cat("\n\n");

	complete_dist=compute_dist(st, dist_type);
	complete_ss=compute_sum_of_squares(complete_dist, members_a, members_b);

	for(i in 1:num_cat_to_probe){
		red_st=st[, -i, drop=F];
		red_dist=compute_dist(red_st, dist_type);
		red_ss=compute_sum_of_squares(red_dist, members_a, members_b);
	}

	cat("\n");
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

cat("Top Categories: \n");
print(category_names[1:30]);
cat("\n\n");


# Compute distances
cat("Computing distances...\n");
all_samp_dist_mat=compute_dist(norm_mat, dist_type);

###############################################################################

hcl=hclust(all_samp_dist_mat, method="ward.D2");

print(hcl);

dendr=as.dendrogram(hcl);

pdf(paste(output_fname_root, ".cl_inf.pdf", sep=""), height=8.5, width=14);
plot(dendr);

nonparm_mds_res=matrix(0,nrow=num_samples,ncol=2);
classic_mds_res=matrix(0,nrow=num_samples,ncol=2);

imds_res=isoMDS(all_samp_dist_mat);
nonparm_mds_res[,1]=imds_res$points[,1]
nonparm_mds_res[,2]=imds_res$points[,2]
#plot(mds_res[,1], mds_res[,2]);

classic_mds_res=cmdscale(all_samp_dist_mat);
#plot(mds_res[,1], mds_res[,2]);

# Output distance matrix
#asFull=as.matrix(dist_mat);
#print(asFull);

for(num_cl in 2:max_clusters){

	cat("Cutting for ", num_cl, " clusters...\n", sep="");
	memberships=cutree(hcl, k=num_cl);
	plot(nonparm_mds_res, col=memberships, main=paste("Num Clusters: ", num_cl));
	plot(classic_mds_res, col=memberships, main=paste("Num Clusters: ", num_cl));
	
	for(i in 1:num_cl){

		members_i=names(memberships[memberships==i]);

		for(j in 1:num_cl){

			if(i>=j){
				next;
			}

			cat("Working on cluster[", i, "] vs cluster[", j, "]\n");

			members_j=names(memberships[memberships==j]);
			sub_mat=norm_mat[c(members_i, members_j),];
			analyze_cluster_pairs(sub_mat, members_i, members_j, num_top_cat, dist_type);
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
