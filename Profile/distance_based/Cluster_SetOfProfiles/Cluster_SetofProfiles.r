#!/usr/bin/env Rscript

###############################################################################

cat("\n\n");

library(MASS)
library('getopt');
library('vegan');
library('plotrix');

source('~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r');
source('~/git/AnalysisTools/Metadata/OutputFileLibrary/OutputFileLibrary.r');


DEF_DISTTYPE="man";
DEF_SUBSAMPLE_SIZE=200;
DEF_COR_CUT=.75;

params=c(
	"input_sumtab_dir", "r", 1, "character",
	"output_filename_root", "o", 1, "character",
	"dist_type", "d", 2, "character",
	"subsamp_size", "s", 2, "numeric",
	"abscorcut", "c", 2, "numeric"
);

options(width=200);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\n\nUsage:\n", script_name, "\n",
	"	-r <input summary_table.tsv directory>\n",
	"	-o <output filename root>\n",
	"\n",
	"	Options:\n",
	"	[-d <euc/wrd/man/bray/horn/bin/gow/tyc/minkp5/minkp3, default =", DEF_DISTTYPE, ">]\n",
	"	[-s <subsample size, default = ", DEF_SUBSAMPLE_SIZE, ">]\n",
	"	[-p <cluster abs correlation cutoff, default = ", DEF_COR_CUT, ">]\n",
	"\n",
	"\n");

if(
	!length(opt$input_sumtab_dir) || 
	!length(opt$output_filename_root) 
){
	cat(usage);
	q(status=-1);
}

InputSumTabDir=opt$input_sumtab_dir;
OutputFileRoot=opt$output_filename_root;

DistType=DEF_DISTTYPE;
if(length(opt$dist_type)){
	DistType=opt$dist_type;
}

SubsampSize=DEF_SUBSAMPLE_SIZE;
if(length(opt$subsamp_size)){
	SubsampSize=opt$subsamp_size;
}

ClusterAbsCorCutoff=DEF_COR_CUT;
if(length(opt$abscorcut)){
	ClusterAbsCorCutoff=opt$abscorcut;
}

if(!any(DistType== c("wrd","man","bray","horn","bin","gow","euc","tyc","minkp3","minkp5"))){
	cat("Error: Specified distance type: ", DistType, " not recognized.\n");
	quit(status=-1);
}


###############################################################################

cat("Input Summary Table Directory:", InputSumTabDir, "\n");
cat("Output Filename Root:", OutputFileRoot, "\n");
cat("Distance Type        :", DistType, "\n");
cat("Subsample Size       :", SubsampSize, "\n");

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

mean_nonrem=function(sumtab){
	means=apply(sumtab, 2, mean);
	sum_means=sum(means);
	return(sum_means);
}

print_stl=function(stl, field=""){
	nm_arr=names(stl);
	for(nm in nm_arr){
		cat(nm, "\n");
		if(field!=""){
			cat(field, ": ", stl[[nm]][[field]], "\n");
		}else{
			rec_names=names(stl[[nm]]);

			rec_names=setdiff(rec_names, c("cnts_mat", "norm_mat", "dist"));

			for(rn in rec_names){

				if(rn=="dist"){
					cat(rn, ":", "\n");
					print(stl[[nm]][["dist"]]);
				}else{
					cat(rn, ":", stl[[nm]][[rn]], "\n");
				}
			}
		}
		cat("\n");
	}
}


###############################################################################

# Get list of summary tables from directory
# Read in summary tables into list
# Sort by abundance of "non remaining"
# Subsample BS subjects (200, for 99%)
# Compute distance matrix for all summary tables
# Starting with most abundant, compute Mantel stat, then 

output_fname_root = paste(OutputFileRoot, ".", DistType, sep="");

##############################################################################
# Find summary tables

get_summary_tables_from_directory=function(dir){
	cat("Looking for summary_table.tsv's in : ", dir, "\n");
	files=system(paste("find -L ", dir), intern=T);
	matches=grep("\\.summary_table\\.tsv$", files, value=T);
	return(matches);
}

sumtab_arr=get_summary_tables_from_directory(InputSumTabDir);
print(sumtab_arr);

num_sumtabs=length(sumtab_arr);
cat("Number of Summary Tables Found: ", num_sumtabs, "\n");

##############################################################################
# Load summary tables

calc_even=function(vals){
	tots=sum(vals);	
	prop=vals/tots;
	shan=-sum(prop*log(prop));
	smax=log(length(vals));
	even=shan/smax;
	return(even);
}

sumtab_list=list();

intersected_samples=c();
for(st in sumtab_arr){
	cnts_mat=load_summary_file(st);
	norm_mat=normalize(cnts_mat);
	num_samples=nrow(norm_mat);

	if(length(intersected_samples)==0){
		intersected_samples=rownames(cnts_mat);
	}else{
		intersected_samples=intersect(intersected_samples, rownames(cnts_mat));
	}

	cat_names=colnames(cnts_mat);
	no_rem_cat=setdiff(cat_names, "Remaining");
	norm_no_remaining=norm_mat[,no_rem_cat,drop=F];
	
	mean_abd=apply(norm_no_remaining, 2, mean);
	print(mean_abd);
	evenness=calc_even(mean_abd);
	mean_cat_values=sum(apply(norm_no_remaining, 2, mean));

	prop_nonzero=sum(apply(norm_no_remaining, 1, function(x){ sum(x)>0;}))/num_samples;

	simp_name=tail(strsplit(st, "/")[[1]],1);
	simp_name=gsub(".summary_table.tsv$", "", simp_name);

	rec=list();
	rec[["simple_name"]]=simp_name;
	rec[["cnts_mat"]]=cnts_mat;
	rec[["norm_mat"]]=norm_mat;
	rec[["mean_non_remaining"]]=signif(mean_cat_values, 4);
	rec[["prop_nonzero"]]=signif(prop_nonzero, 4);
	rec[["num_cat"]]=length(no_rem_cat);
	rec[["evenness"]]=evenness;
	rec[["num_samples"]]=nrow(cnts_mat);
	sumtab_list[[st]]=rec;
}

sort_by_MNR=function(stl){
	mnr_val=sapply(stl, function(x){return(x[["mean_non_remaining"]])});
	sorted_stl=stl[order(mnr_val, decreasing=T)];
	return(sorted_stl);
}

cat("\n\n");
cat("Sorting by Mean Non-Remaining:\n");
sumtab_list=sort_by_MNR(sumtab_list);

num_intersected_samples=length(intersected_samples);

cat("Num Interesected Samples: ", num_intersected_samples, "\n");


#SubsampSize=5;
cat("Subsample Size for Distance Mat Calc: ", SubsampSize, "\n");
subsamp_ix=sample(num_intersected_samples, SubsampSize);
subsamp_id=intersected_samples[subsamp_ix];

print(subsamp_id);


##############################################################################
# Compute distance matrices

cat("Calculating Distance Matrices (on subsample)\n");

for(st in sumtab_arr){
	norm_mat=sumtab_list[[st]][["norm_mat"]];
	norm_mat_subsmp=norm_mat[subsamp_id,];
	sumtab_list[[st]][["dist"]]=compute_dist(norm_mat_subsmp, DistType);
}

print_stl(sumtab_list);

##############################################################################

calc_cor=function(dist_mat_a, dist_mat_b){
	matcor=cor(dist_mat_a, dist_mat_b);
	return(matcor);
}

cluster_sumtab=function(st_list, cutoff){

	cat("Cutoff: ", cutoff, "\n");
	num_st=length(st_list);
	st_names=names(st_list);

	# Overal Cluster Matrix
	# Cluster ID, Name, abscor_from_rep, non_rem_abund
	header_names=c("ClusterID", "CorFromRep", "SimpleName", 
		"NonRemAbund", "PropNonzero", "NumChildren", "ChildEvenness");
	memb_info_mat=matrix(NA, nrow=num_st, ncol=length(header_names));
	colnames(memb_info_mat)=header_names;
	rownames(memb_info_mat)=st_names;

	init_row=function(cl_id, abs_cor, name){
		arr=c(cl_id, abs_cor, 
			st_list[[name]][["simple_name"]],
			st_list[[name]][["mean_non_remaining"]],
			round(st_list[[name]][["prop_nonzero"]], 4),
			st_list[[name]][["num_cat"]],
			round(st_list[[name]][["evenness"]], 4)
		);
		return(arr);	
	}

	cluster_list=list(); #list of arrays, where first member of array is it's main rep
	
	first_rep_name=st_names[1];
	cluster_list[[1]]=first_rep_name;
	memb_info_mat[st_names[1],]=init_row(1, 1, first_rep_name);
	st_names=st_names[2:length(st_names)];

	# Loop for each of the members to cluster
	st_count=2;
	for(nm in st_names){

		cat("Working on: ", nm, ": ", st_count, "/", num_st, "\n");
		# Loop for each of the cluster to look through
		cluster_found=F;
		num_clusters=length(cluster_list);
		cor_arr=numeric(num_clusters);
		for(i in 1:length(cluster_list)){
			# compare
			cluster_rep=cluster_list[[i]][1];
			dist_cor=calc_cor(st_list[[cluster_rep]][["dist"]], st_list[[nm]][["dist"]]);
			cor_arr[i]=dist_cor;
		}

		# Find closest cluster
		#cat("Correlations: \n");
		#print(abs_cor_arr);
		max_cor=max(cor_arr);
		cat("Max Cor: ", max_cor, "\n");
		if(is.na(max_cor)){
			memb_info_mat[nm,]=init_row(NA, NA, nm);
		}else if(max_cor>cutoff){
			clus_id=which(max_cor==cor_arr);
			# |cor| is > cutoff put it in the cluster
			cluster_list[[clus_id]]=c(cluster_list[[clus_id]], nm);
			memb_info_mat[nm,]=
				init_row(clus_id, 
					round(max_cor, 4), 
					nm);
		}else{	
			# add member to end of cluster list
			num_cluster=length(cluster_list);
			cluster_list[[num_cluster+1]]=nm;
			memb_info_mat[nm,]=
				init_row(num_cluster+1, 1, nm);
		}

		st_count=st_count+1;
	}

	# Sort by cluster ID
	by_id=order(as.numeric(memb_info_mat[,"ClusterID"]), method="shell");
	memb_info_mat=memb_info_mat[by_id,,drop=F];
	return(memb_info_mat);

}

#------------------------------------------------------------------------------

export_cluster_table=function(cltab, fname_root, cutoff){
	cnames=colnames(cltab);
	out_tab=cbind(rownames(cltab), cltab);
	colnames(out_tab)=c("Name", cnames);
	outtabfn=paste(fname_root, ".", round(cutoff*100,0), ".cluster_table.tsv", sep="");
	write.table(out_tab, outtabfn, quote=F, sep="\t", row.names=F);
}

export_cluster_representatives_table=function(cltab, fname_root, cutoff){

	print(cltab[,"ClusterID"]);
	cluster_numbers=as.numeric(cltab[,"ClusterID"]);
	cluster_numbers=cluster_numbers[!is.na(cluster_numbers)];
	num_clusters=max(cluster_numbers, na.rm=T);
	cat("Num Clusters: ", num_clusters, "\n");

	num_rows=length(cluster_numbers);
	prev_cl=0;
	clust_rep_ix=numeric(num_clusters);
	crp_ix=1;
	for(i in 1:num_rows){
		if(cluster_numbers[i]!=prev_cl){
			clust_rep_ix[crp_ix]=i;
			prev_cl=cluster_numbers[i];
			crp_ix=crp_ix+1;
		};
	}

	cltab=cltab[clust_rep_ix,,drop=F];

	cnames=colnames(cltab);
	out_tab=cbind(rownames(cltab), cltab);
	colnames(out_tab)=c("Name", cnames);
	outtabfn=paste(fname_root, ".", round(cutoff*100,0), ".cluster_repres.tsv", sep="");
	write.table(out_tab, outtabfn, quote=F, sep="\t", row.names=F);
}

#------------------------------------------------------------------------------

for(cut in ClusterAbsCorCutoff+c(-.1, 0, .1)){
	cluster_table=cluster_sumtab(sumtab_list, cut);
	export_cluster_table(cluster_table, OutputFileRoot, cut);
	export_cluster_representatives_table(cluster_table, OutputFileRoot, cut);
}

##############################################################################

cat("\nDone.\n")
dev.off();
warn=warnings();
if(length(warn)){
        print(warn);
}
q(status=0);
