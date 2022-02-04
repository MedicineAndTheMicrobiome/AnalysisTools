#!/usr/bin/env Rscript

###############################################################################

cat("\n\n");

library(MASS)
library('getopt');
library('vegan');
library('plotrix');


DEF_DISTTYPE="man";

params=c(
	"input_summary_table_A", "a", 1, "character",
	"input_summary_table_B", "b", 1, "character",
	"mapping_file", "m", 1, "character",
	"mapping_file_Aname", "A", 2, "character",
	"mapping_file_Bname", "B", 2, "character",
	
	"include_list_A", "i", 2, "character",
	"include_list_B", "j", 2, "character",
	"output_filename_root", "o", 1, "character",
	"dist_type", "d", 2, "character",
	"tag_name", "t", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\n\nUsage:\n", script_name, "\n",
	"	-a <input summary_table.tsv file A>\n",
	"	-b <input summary_table.tsv file B>\n",
	"	-m <mapping file, from A to B Identifiers>\n",
	"	-A <for mapping file, A's column name>\n",
	"	-B <for mapping file, B's column name>\n",
	"	-o <output file root name>\n",
	"\n",
	"	Options:\n",
	"	[-d <euc/wrd/man/bray/horn/bin/gow/tyc/minkp5/minkp3, default =", DEF_DISTTYPE, ">]\n",
	"	[-i <include list of IDs for summary_table.tsv A, may be a factor file>]\n",
	"	[-j <include list of IDs for summary_table.tsv B, may be a factor file>]\n",
	"\n",
	"	[-t <tag name>]\n",
	"\n",
	"\n",
	"This script will compute distance matrices for the two summary tables independently, then\n",
	"use various methods to compare their clustering.  Since the summary tables do not have\n",
	"to share the same categories/taxa/function, the method is useful for integrating\n",
	"'omics results with different biomarkers.\n",
	"\n",
	"The first half of comparisons compare hierarchical clusters at different cuts, so that\n",
	"contingency tables can be calculated.  Table-wise and cell-wise chi-squared tests are\n",
	"performed.  Table-wise and cell-wise conditional probabilities are also calculated.\n",
	"Dendrogram samples/leaves are colored by the group to see how clusters members in the\n",
	"first set of clusters distribute in the second set of clusters.\n",
	"\n",
	"Jumping down to around 3/4 though the output, a table/heatmap identifies which pair of\n",
	"cut combinations leads to the best predictable (least random) relationship between the\n",
	"sample types.  Summary pvalues and conditional probabilities summary tables are provide\n",
	"and the best combination of clusters cuts is reported based on minimum table-wise p-value\n",
	"\n",
	"The second set of analyses focuses on using the Mantel statistic (correlation of distances)\n",
	"between 2 distance matrices to understand overall clustering structure.  Since it's\n",
	"possible for some clusters in a hierarchical cluster to share structure, while other\n",
	"clusters do not, the Mantel statistic is calculated iteratively across cuts using the\n",
	"clustering of either set of samples in two passes.  Summary plots for each of the\n", 
	"cuts show how correlation between sample types change as they are subdivided.\n",
	"\n",
	"The format of the sample mapping file:\n",
	"<sample group A name>\\t<sample group B name>\n",
	"<id.A.1>\\t<id.B.1>\n",
	"<id.A.2>\\t<id.B.2>\n",
	"<id.A.3>\\t<id.B.3>\n",
	"...\n",
	"<id.A.N>\\t<id.B.N>\n",
	"\n",
	"\n",
	"For the distance types:\n",
	" minkp5 is the minkowski with p=.5, i.e. sum((x_i-y_i)^1/2)^2\n",
	" minkp3 is the minkowski with p=.3, i.e. sum((x_i-y_i)^1/3)^3\n",
	"\n");

if(
	!length(opt$input_summary_table_A) || 
	!length(opt$input_summary_table_B) || 
	!length(opt$mapping_file) || 
	!length(opt$mapping_file_Aname) || 
	!length(opt$mapping_file_Bname) || 
	!length(opt$output_filename_root) 
){
	cat(usage);
	q(status=-1);
}

InputSumTabA=opt$input_summary_table_A;
InputSumTabB=opt$input_summary_table_B;
MappingFile=opt$mapping_file;
MappingFileAname=opt$mapping_file_Aname;
MappingFileBname=opt$mapping_file_Bname;
OutputFileRoot=opt$output_filename_root;

DistType=DEF_DISTTYPE;
if(length(opt$dist_type)){
	DistType=opt$dist_type;
}


if(!any(DistType== c("wrd","man","bray","horn","bin","gow","euc","tyc","minkp3","minkp5"))){
	cat("Error: Specified distance type: ", DistType, " not recognized.\n");
	quit(status=-1);
}

OutputFileRoot_nodist=OutputFileRoot;
OutputFileRoot=paste(OutputFileRoot, ".", DistType, sep="");

Include_A_list="";
Include_B_list="";
if(length(opt$include_list_A)){
	Include_A_list=opt$include_list_A;
}
if(length(opt$include_list_B)){
	Include_B_list=opt$include_list_B;
}

if(length(opt$tag_name)){
        TagName=opt$tag_name;
        cat("Setting TagName Hook: ", TagName, "\n");
        setHook("plot.new",
                function(){
                        #cat("Hook called.\n");
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

cat("Input Summary Table A:", InputSumTabA, "\n");
cat("Input Summary Table B:", InputSumTabB, "\n");
cat("Mapping File         :", MappingFile, "\n");
cat("Mapping File Aname   :", MappingFileAname, "\n");
cat("Mapping File Bname   :", MappingFileBname, "\n");
cat("Output File          :", OutputFileRoot, "\n");
cat("Distance Type        :", DistType, "\n");
cat("\n");
if(Include_A_list!=""){
	cat("Include List for A:", Include_A_list, "\n");
}
if(Include_B_list!=""){
	cat("Include List for B:", Include_B_list, "\n");
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
	inmat=as.matrix(read.delim(st_fname, sep="\t", header=TRUE, row.names=1, 
		check.names=FALSE, comment.char="", quote=""))

	num_categories=ncol(inmat)-1;
	num_samples=nrow(inmat);

	cat("Loaded Summary Table: ", st_fname, "\n", sep="");
	cat("  Num Categories: ", num_categories, "\n", sep="");
	cat("  Num Samples: ", num_samples, "\n", sep="");

	countsmat=inmat[,2:(num_categories+1)];

	return(countsmat);
}

#------------------------------------------------------------------------------

load_mapping_file=function(mp_fname, keep_a_ids, keep_b_ids, colname_a, colname_b){

	num_keep_a=length(keep_a_ids);
	num_keep_b=length(keep_b_ids);
	cat("Num A's IDs to keep: ", num_keep_a, "\n");
	cat("Num B's IDs to keep: ", num_keep_b, "\n");

	inmat=as.matrix(read.delim(mp_fname, sep="\t", header=TRUE, check.names=F, comment.char="", quote=""));

	inmat=inmat[,c(colname_a, colname_b)];

	# Remove unpaired
	keep_ix=apply(inmat, 1, function(x){ all(!is.na(x))});
	inmat=inmat[keep_ix,,drop=F];

	# Keep Entry if record is in both lists
	keep_ix=c();
	orig_mat_rows=nrow(inmat);
	cat("Number of Mapping Entries Read: ", orig_mat_rows, "\n");
	for(i in 1:orig_mat_rows){
		if(any(inmat[i,1]==keep_a_ids) &&  any(inmat[i,2]==keep_b_ids)){
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

plot_text=function(strings){

        orig_par=par(no.readonly=T);
	options(width=100);

        par(family="Courier");
        par(mar=rep(0,4));

        num_lines=length(strings);

        top=max(as.integer(num_lines), 50);

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
	
	if(orig_par$fig[1]<0){
		orig_par$fig[1]=0;
	}

	#print(orig_par);
        par(orig_par);

}

load_list=function(fn){
	cat("Loading ", fn, " as list...\n");
	data=read.delim(fn, header=F, sep="\t", quote="", as.is=T);
	res=data[,1];
	cat(" Number of ID's loaded: ", length(res), "\n");
	return(res);
}


###############################################################################

output_fname_root = paste(OutputFileRoot, ".", DistType, sep="");

cat("\n");
cat("Loading summary table A:", InputSumTabA, "\n");
counts_mat_A=load_summary_table(InputSumTabA);
cat("Loading summary table B:", InputSumTabB, "\n");
counts_mat_B=load_summary_table(InputSumTabB);

if(Include_A_list!=""){
	cat("Include List for A:", Include_A_list, "\n");
	include_a_arr=load_list(Include_A_list);
	nsamp_before=nrow(counts_mat_A);
	counts_mat_A=counts_mat_A[intersect(include_a_arr, rownames(counts_mat_A)),,drop=F];
	nsamp_after=nrow(counts_mat_A);
	cat("Num Samples kept: ", nsamp_after, "/", nsamp_before, "\n", sep="");
}
if(Include_B_list!=""){
	cat("Include List for B:", Include_B_list, "\n");
	include_b_arr=load_list(Include_B_list);
	nsamp_before=nrow(counts_mat_A);
	counts_mat_B=counts_mat_B[intersect(include_b_arr, rownames(counts_mat_B)),,drop=F];
	nsamp_after=nrow(counts_mat_A);
	cat("Num Samples kept: ", nsamp_after, "/", nsamp_before, "\n", sep="");
}
cat("\n\n");

# Reconcile summary table IDs through mapping file
samples_stA=rownames(counts_mat_A);
samples_stB=rownames(counts_mat_B);

cat("Loading Mapping file:", MappingFile, "\n");
map_info=load_mapping_file(MappingFile, samples_stA, samples_stB, MappingFileAname, MappingFileBname);

print(map_info);

cat("Removing samples without complete mappings...\n");
counts_mat_A=counts_mat_A[map_info[["a_id"]],];
counts_mat_B=counts_mat_B[map_info[["b_id"]],];

num_samples=nrow(counts_mat_A);
cat("Num usable samples: ", num_samples, "\n");

###############################################################################

# Normalize counts
cat("Normalizing counts...\n");
norm_mat_A=normalize(counts_mat_A);
norm_mat_B=normalize(counts_mat_B);

###############################################################################

# Compute full distances
cat("Computing distances...\n");
dist_mat_A=compute_dist(norm_mat_A, DistType);
dist_mat_B=compute_dist(norm_mat_B, DistType);

###############################################################################

pdf(paste(OutputFileRoot, ".cmp_dist.pdf", sep=""), height=8.5, width=8);

param_summary=capture.output({
	cat("Input Summary Table A:", InputSumTabA, "\n");
	cat("Input Summary Table B:", InputSumTabB, "\n");
	cat("Mapping File         :", MappingFile, "\n");
	cat("Output File          :", OutputFileRoot, "\n");
	cat("Distance Type        :", DistType, "\n");
	cat("\n");
	
});
plot_text(param_summary);

###############################################################################

test=F;

max_cuts=log2(num_samples);

palette_col=c("red", "green", "blue", "cyan", "magenta", "orange", "gray", 
	"pink", "black", "purple", "brown", "aquamarine");
num_pref_col=length(palette_col);
num_clus=length(unique(max_cuts));	

if(max_cuts>num_pref_col){
	palette_col=rainbow(n=max_cuts, start=0, end=4/6);
}

palette(palette_col);

###############################################################################

find_extremes=function(pts){

	num_pts=nrow(pts);

	pts_sorted_by_x=sort(pts[,1]);
	pts_sorted_by_y=sort(pts[,2]);

	# Find distances from center
	centroid=apply(pts, 2, median);
	print(centroid);
	dc=dist_from_centr=apply(pts, 1, function(x){
		sqrt((x[1]-centroid[1])^2+(x[2]-centroid[2])^2)})
	dc_sort=sort(dc, decreasing=T);

	sample_ids=unique(c(
		names(pts_sorted_by_x[1]),
		names(pts_sorted_by_y[1]),
		names(pts_sorted_by_x[num_pts]),
		names(pts_sorted_by_y[num_pts]),
		names(dc_sort[1:4])
	));

	return(sample_ids);

}

expand_range=function(x, fact=.1){
	mn=min(x);
	mx=max(x);
	diff=(mx-mn);
	return(c(mn-fact*diff, mx+fact*diff));
	
}

compare_mds=function(apts, bpts, type, aclus, bclus, aname, bname){

	num_samp=nrow(apts);
	
	aout=find_extremes(apts);
	bout=find_extremes(bpts);

	ax=1:num_samp;
	bx=1:num_samp;

	names(ax)=rownames(apts);
	names(bx)=rownames(bpts);

	ax=ax[aout];
	bx=bx[bout];

	bothx=unique(c(ax, bx));

	cat("Extremes/Outliers Labeled:\n");
	aoutnames=rownames(apts)[bothx];
	boutnames=rownames(bpts)[bothx];
	print(aoutnames);
	print(boutnames);

	par(mfrow=c(2,2));

	plot(apts, main=paste(type, ": ", aname, sep=""), xlab="Dim 1", ylab="Dim 2", 
		col=aclus, xlim=expand_range(apts));
	text(apts[bothx,], aoutnames, cex=.3);
	mtext(paste("Colored by ", aname, sep=""));
	plot(bpts, main=paste(type, ": ", bname, sep=""), xlab="Dim 1", ylab="Dim 2", 
		col=aclus, xlim=expand_range(bpts));
	text(bpts[bothx,], boutnames, cex=.3);

	plot(apts, main=paste(type, ": ", aname, sep=""), xlab="Dim 1", ylab="Dim 2", 
		col=bclus, xlim=expand_range(apts));
	text(apts[bothx,], aoutnames, cex=.3);
	plot(bpts, main=paste(type, ": ", bname, sep=""), xlab="Dim 1", ylab="Dim 2", 
		col=bclus, xlim=expand_range(bpts));
	text(bpts[bothx,], boutnames, cex=.3);
	mtext(paste("Colored by ", bname, sep=""));

}

compute_pseudof=function(clsmem, distmat){

	# Compute SS Between
	#centroid_distmat=compute_centroids(clsmem, distmat);

	# Compute SS Within
	num_grps=length(unique(clsmem));
	num_samples=length(clsmem);
	distmat2d=as.matrix(distmat);

	tsw=0;
	tsb=0;
	nsb=0;
	nsw=0;

	for(i in 1:num_grps){
		grp_ix=(clsmem==i);
		for(j in 1:num_grps){
			grp_jx=(clsmem==j);

			subdist=distmat2d[grp_ix, grp_jx];

			if(i==j){
				
				# cat("SSW:\n");
				# print(dim(subdist));
				ssd=sum(subdist^2)/2;

				rows=nrow(subdist);
				num_dist=rows*(rows-1)/2;	# Number of distances in half matrix	
				
				tsw=tsw+ssd;
				nsw=nsw+num_dist;
			}else{

				# cat("SSB:\n");
				# print(dim(subdist));
				ssd=sum(subdist^2);

				num_dist=nrow(subdist)*ncol(subdist);

				tsb=tsb+ssd;
				nsb=nsb+num_dist;
			}

		}
	}

	# Compute sum of SS's based on Mean Squared
	ssb=(tsb/nsb)*num_grps;
	ssw=(tsw/nsw)*num_samples;

	# Compute Variance by adjusting for degree of freedom
	b_var=ssb/(num_grps-1);
	w_var=ssw/(num_samples-num_grps);

	# Return pseudo F-stat
	pseudof=b_var/w_var;
	if(!length(pseudof)){
		pseudof=0;
	}
	return(pseudof);

}

compare_pseudof=function(dista, distb, grp_hcla, grp_hclb, max_k, namea, nameb){

	amsd_byA=numeric(max_k-1);
	bmsd_byA=numeric(max_k-1);
	amsd_byB=numeric(max_k-1);
	bmsd_byB=numeric(max_k-1);

	for(clix in 2:max_k){	
		cat("Cutting tree to k=", clix, "\n");
		mem_byA=cutree(grp_hcla, k=clix);
		mem_byB=cutree(grp_hclb, k=clix);
		
		cat("Computing SSD:\n");
		amsd_byA[clix-1]=compute_pseudof(mem_byA, dista);
		bmsd_byA[clix-1]=compute_pseudof(mem_byA, distb);
		amsd_byB[clix-1]=compute_pseudof(mem_byB, dista);
		bmsd_byB[clix-1]=compute_pseudof(mem_byB, distb);
	}

	max_asmd=max(amsd_byA, amsd_byB, na.rm=T);
	max_bsmd=max(bmsd_byA, bmsd_byB, na.rm=T);
	min_asmd=min(amsd_byA, amsd_byB, na.rm=T);
	min_bsmd=min(bmsd_byA, bmsd_byB, na.rm=T);

	asmd_range=c(min_asmd, max_asmd);
	bsmd_range=c(min_bsmd, max_bsmd);

	par(mfrow=c(2,3));


	# Clustering from A's perspective
	plot(2:max_k, amsd_byA, xlab="Num Clusters, k", ylab="Pseudo F-Stat", 
		main=paste("Pseudo F-Stat: ", namea, sep=""), type="b", ylim=asmd_range);
	mtext(paste("Clustered by Optimal ", namea, " Groupings", sep=""), cex=.6);

	plot(2:max_k, amsd_byB, xlab="Num Clusters, k", ylab="Pseudo F-Stat", 
		main=paste("Pseudo F-Stat: ", namea, sep=""), type="b", ylim=asmd_range);
	mtext(paste("Clustered by Optimal ", nameb, " Groupings", sep=""), cex=.6);

	lograt=log(amsd_byB/amsd_byA);
	lograt[!is.finite(lograt)]=0;
	lims=c(-1,1)*max(abs(lograt), rm.na=T);
	plot(2:max_k, lograt, xlab="Num Clusters, k", 
		ylab=paste("Pseudo F-Stat LogRatio(", namea, ")", sep=""),  
		main=paste("Pseudo F-stat Ratio"), type="b", ylim=lims);
	abline(h=0, col="blue");
	mtext(paste(namea, ": By ", nameb, "/", namea, " Groupings", sep=""), cex=.6);


	# Clustering from B's perspective
	plot(2:max_k, bmsd_byB, xlab="Num Clusters, k", ylab="Pseudo F-Stat", 
		main=paste("Pseudo F-stat: ", nameb, sep=""), type="b", ylim=bsmd_range);
	mtext(paste("Clustered by Optimal ", nameb, " Groupings", sep=""), cex=.6);

	plot(2:max_k, bmsd_byA, xlab="Num Clusters, k", ylab="Pseudo F-Stat", 
		main=paste("Pseudo F-stat: ", nameb, sep=""), type="b", ylim=bsmd_range);
	mtext(paste("Clustered by Optimal ", namea, " Groupings", sep=""), cex=.6);

	lograt=log(bmsd_byA/bmsd_byB);
	lograt[!is.finite(lograt)]=0;
	lims=c(-1,1)*max(abs(lograt), rm.na=T);
	plot(2:max_k, lograt, xlab="Num Clusters, k", 
		ylab=paste("Pseudo F-Stat LogRatio(", nameb, ")", sep=""),  
		main=paste("Pseudo F-stat Ratio"), type="b", ylim=lims);
	abline(h=0, col="blue");
	mtext(paste(nameb, ": By ", namea, "/", nameb, " Groupings", sep=""), cex=.6);


}

compare_dendrograms=function(hclA, hclB, num_cuts, namea, nameb, idsb){

	color_denfun_bySample=function(n){
		if(is.leaf(n)){
			leaf_attr=attributes(n);
			leaf_name=leaf_attr$label;
			ind_color=sample_to_color_map[leaf_name];
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

	orig_par=par(no.readonly=T);

	par(mar=c(10,2,3,.5));
	par(mfrow=c(2,1));

	dendra=as.dendrogram(hclA);
	dendrb=as.dendrogram(hclB);

	label_scale=.2;
	dendra=dendrapply(dendra, text_scale_denfun);
	dendrb=dendrapply(dendrb, text_scale_denfun);

	mem_byA=cutree(hclA, k=num_cuts);
	sample_to_color_map=mem_byA;
	dendra=dendrapply(dendra, color_denfun_bySample);
	names(sample_to_color_map)=idsb;
	dendrb=dendrapply(dendrb, color_denfun_bySample);
	

	plot(dendra, main=paste("Cut by ", namea, "'s clustering", sep=""));
	plot(dendrb, main=nameb);

	# Plot shared statistics
	
	shared_mat=matrix(0, nrow=num_cuts, ncol=num_cuts);
	samples_a=names(mem_byA);

	mem_byB=cutree(hclB, k=num_cuts);
	names(mem_byB)=names(mem_byA);

	num_samples=length(mem_byA);
	cat("Num Samples: ", num_samples, "\n", sep="");

	for(r in 1:num_cuts){

		samples_in_r=samples_a[(mem_byA==r)];

		for(c in 1:num_cuts){
		
			total=sum(mem_byB==c);
			samples_in_c=(mem_byB[samples_in_r]==c);
			overlapping=sum(mem_byB[samples_in_r]==c);

			shared_mat[r, c]=overlapping/num_samples;
		}
	}
	
	#print(shared_mat);

	par(orig_par);



}

##############################################################################


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

reorder_member_ids=function(members_cut, dendr_names){

	grp_mids=get_middle_of_groups(dendr_names, members_cut);

        # Reorder cluster assignments to match dendrogram left/right
        plot_order=order(grp_mids);
        mem_tmp=numeric(num_samples);
	num_cl=length(unique(members_cut));
        for(gr_ix in 1:num_cl){
		old_id=(members_cut==plot_order[gr_ix]);
		mem_tmp[old_id]=gr_ix;
        }
        names(mem_tmp)=names(members_cut);
        members_cut=mem_tmp;
	return(members_cut);
} 

remap_coord=function(x, sbeg, send, dbeg, dend){
	srang=send-sbeg;
	norm=(x-sbeg)/srang;
	drang=dend-dbeg;
	return(norm*drang+dbeg);
}
##############################################################################

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, deci_pts=4,
        label_zeros=T, counts=F, value.cex=1,
        plot_col_dendr=F,
        plot_row_dendr=F
){

        num_row=nrow(mat);
        num_col=ncol(mat);

        row_names=rownames(mat);
        col_names=colnames(mat);

        orig.par=par(no.readonly=T);

        cat("Num Rows: ", num_row, "\n");
        cat("Num Cols: ", num_col, "\n");

        # Flips the rows, so becuase origin is bottom left
        mat=mat[rev(1:num_row),, drop=F];

        # Generate a column scheme
        num_colors=50;
        color_arr=rainbow(num_colors, start=0, end=4/6);
        if(high_is_hot){
                color_arr=rev(color_arr);
        }

        # Provide a means to map values to an (color) index
        remap=function(in_val, in_range, out_range){
                in_prop=(in_val-in_range[1])/(in_range[2]-in_range[1])
                out_val=in_prop*(out_range[2]-out_range[1])+out_range[1];
                return(out_val);
        }

        # If range is not specified, find it based on the data
        if(is.na(plot_min)){
                plot_min=min(mat, na.rm=T);
        }
        if(is.na(plot_max)){
                plot_max=max(mat, na.rm=T);
        }

        if(plot_min>=-1 && plot_max<=1){
                fractions_only=T;
        }else{
                fractions_only=F;
        }
        cat("Plot min/max: ", plot_min, "/", plot_max, "\n");

        # Get Label lengths
        row_max_nchar=max(nchar(row_names));
        col_max_nchar=max(nchar(col_names));
        cat("Max Row Names Length: ", row_max_nchar, "\n");
        cat("Max Col Names Length: ", col_max_nchar, "\n");

        ##################################################################################################

        get_dendrogram=function(in_mat, type){
                if(type=="row"){
                        dendist=dist(in_mat);
                }else{
                        dendist=dist(t(in_mat));
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


                hcl=hclust(dendist, method="ward.D2");
                dend=list();
                dend[["tree"]]=as.dendrogram(hcl);
                dend[["names"]]=get_clstrd_leaf_names(dend[["tree"]]);
                return(dend);
        }


        ##################################################################################################
        # Comput Layouts
        col_dend_height=ceiling(num_row*.1);
        row_dend_width=ceiling(num_col*.2);

        heatmap_height=num_row;
        heatmap_width=num_col;

        if(plot_col_dendr && plot_row_dendr){
                layoutmat=matrix(
                        c(
                        rep(c(rep(4, row_dend_width), rep(3, heatmap_width)), col_dend_height),
                        rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height)
                        ), byrow=T, ncol=row_dend_width+heatmap_width);

                col_dendr=get_dendrogram(mat, type="col");
                row_dendr=get_dendrogram(mat, type="row");

                mat=mat[row_dendr[["names"]], col_dendr[["names"]]];

        }else if(plot_col_dendr){
                layoutmat=matrix(
                        c(
                        rep(rep(2, heatmap_width), col_dend_height),
                        rep(rep(1, heatmap_width), heatmap_height)
                        ), byrow=T, ncol=heatmap_width);

                col_dendr=get_dendrogram(mat, type="col");
                mat=mat[, col_dendr[["names"]]];

        }else if(plot_row_dendr){
                layoutmat=matrix(
                        rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height),
                        byrow=T, ncol=row_dend_width+heatmap_width);

                row_dendr=get_dendrogram(mat, type="row");
                mat=mat[row_dendr[["names"]],];
        }else{
                layoutmat=matrix(
                        rep(1, heatmap_height*heatmap_width),
                        byrow=T, ncol=heatmap_width);
        }

        #print(layoutmat);
        layout(layoutmat);

        ##################################################################################################

        par(oma=c(col_max_nchar*.60, 0, 3, row_max_nchar*.60));
        par(mar=c(0,0,0,0));
        plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), xaxt="n", yaxt="n", bty="n", xlab="", ylab="");
        mtext(title, side=3, line=0, outer=T, font=2);

        # x-axis
        axis(side=1, at=seq(.5, num_col-.5, 1), labels=colnames(mat), las=2, line=-1.75);
        axis(side=4, at=seq(.5, num_row-.5, 1), labels=rownames(mat), las=2, line=-1.75);

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

                        if(is.na(mat[y,x]) || mat[y,x]!=0 || label_zeros){
                                if(counts){
                                        text_lab=sprintf("%i", mat[y,x]);
                                }else{
                                        text_lab=sprintf(paste("%0.", deci_pts, "f", sep=""), mat[y,x]);
                                        if(fractions_only){
                                                if(!is.na(mat[y,x])){
                                                        if(mat[y,x]==-1 || mat[y,x]==1){
                                                                text_lab=as.integer(mat[y,x]);
                                                        }else{
                                                                text_lab=gsub("0\\.","\\.", text_lab);
                                                        }
                                                }
                                        }
                                }
                                text(x-.5, y-.5, text_lab, srt=atan(num_col/num_row)/pi*180, cex=value.cex, font=2);
                        }
                }
        }

        ##################################################################################################

        par(mar=c(0, 0, 0, 0));

        if(plot_row_dendr && plot_col_dendr){
                rdh=attributes(row_dendr[["tree"]])$height;
                cdh=attributes(col_dendr[["tree"]])$height;
                plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", bty="n", xlim=c(rdh, 0));
                plot(col_dendr[["tree"]], leaflab="none",xaxt="n", yaxt="n", bty="n", ylim=c(0, cdh));
                plot(0,0, type="n", bty="n", xaxt="n", yaxt="n");
                #text(0,0, "Placeholder");
        }else if(plot_row_dendr){
                rdh=attributes(row_dendr[["tree"]])$height;
                plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", bty="n", xlim=c(rdh, 0));
                #text(0,0, "Row Dendrogram");
        }else if(plot_col_dendr){
                cdh=attributes(col_dendr[["tree"]])$height;
                plot(col_dendr[["tree"]], leaflab="none", xaxt="n", yaxt="n", bty="n", ylim=c(0, cdh));
                #text(0,0, "Col Dendrogram");
        }

        par(orig.par);

}


##############################################################################

plot_dendro_contigency=function(hclA, hclB, acuts, bcuts, namea, nameb, idsb){

	color_denfun_bySample=function(n){
		if(is.leaf(n)){
			leaf_attr=attributes(n);
			leaf_name=leaf_attr$label;
			ind_color=sample_to_color_map[leaf_name];
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

	# Compute
	cat("Working on: ", nameb, ": ", bcuts, " x ", namea, ": ", acuts, "\n", sep="");

	dendra=as.dendrogram(hclA);
	dendrb=as.dendrogram(hclB);

	memb_byA=cutree(hclA, k=acuts);
	memb_byB=cutree(hclB, k=bcuts);

	dendr_names_a=get_clstrd_leaf_names(dendra);
	dendr_names_b=get_clstrd_leaf_names(dendrb);

	memb_byA=reorder_member_ids(memb_byA, dendr_names_a);
	memb_byB=reorder_member_ids(memb_byB, dendr_names_b);

	dend_mids_a=get_middle_of_groups(dendr_names_a, memb_byA);
	dend_mids_b=get_middle_of_groups(dendr_names_b, memb_byB);

	num_members=length(memb_byA);

	grp_cnts_a=(table(memb_byA)[1:acuts]);
	grp_cnts_b=(table(memb_byB)[1:bcuts]);;
	
	grp_prop_a=grp_cnts_a/num_members;
	grp_prop_b=grp_cnts_b/num_members;;

	cat("Group Proportions A:\n");
	print(grp_prop_a);

	cat("Group Proportions B:\n");
	print(grp_prop_b);

	# Count up observed
	ab_cnts_obs_mat=matrix(0, nrow=bcuts, ncol=acuts);
	for(i in 1:num_members){
		ma=memb_byA[i];
		mb=memb_byB[i];
		ab_cnts_obs_mat[mb, ma]=ab_cnts_obs_mat[mb, ma]+1;
	}

	ab_prop_obs_mat=ab_cnts_obs_mat/num_members;
	
	# Calculate expected
	ab_prop_exp_mat=grp_prop_b %*% t(grp_prop_a);
	ab_cnts_exp_mat=ab_prop_exp_mat * num_members;

	cat("Observed Counts\n");
	print(ab_cnts_obs_mat);

	cat("Observed Proportions\n");
	print(ab_prop_obs_mat);

	cat("Expected Counts:\n");
	print(ab_cnts_exp_mat);

	cat("Expected Proportions:\n");
	print(ab_prop_exp_mat);
	
	# Compute pvalue for contingency table
	cst=chisq.test(ab_cnts_obs_mat);
	ct_cst_pval=cst$p.value;

	#print(cst);
	#print(cst$observed);
	#print(cst$expected);

	fish_exact_mat=matrix(0, nrow=bcuts, ncol=acuts);

	for(rix in 1:bcuts){
		for(cix in 1:acuts){

			twobytwo=matrix(0, nrow=2, ncol=2);
			twobytwo[1,1]=ab_cnts_obs_mat[rix, cix];
			twobytwo[1,2]=sum(ab_cnts_obs_mat[-rix, cix]);
			twobytwo[2,1]=sum(ab_cnts_obs_mat[rix, -cix]);
			twobytwo[2,2]=sum(ab_cnts_obs_mat[-rix, -cix]);
			ind_cst_res=fisher.test(twobytwo);
			fish_exact_mat[rix, cix]=ind_cst_res$p.value;
		}
	}
	#print(fish_exact_mat);

	# Calculate conditional probabilities
	cnd_pr_agb=matrix(0, nrow=bcuts, ncol=acuts);
	cnd_pr_bga=matrix(0, nrow=bcuts, ncol=acuts);

	for(rix in 1:bcuts){
		cnd_pr_agb[rix,]=ab_prop_obs_mat[rix,]/grp_prop_b[rix];
	}

	for(cix in 1:acuts){
		cnd_pr_bga[,cix]=ab_prop_obs_mat[,cix]/grp_prop_a[cix];
	}

	# Calculate cumulative probability for highest probability mapping
	cum_pr_bga=sum(apply(ab_prop_obs_mat, 2, max));
	cum_pr_agb=sum(apply(ab_prop_obs_mat, 1, max));
	

	##########################################
	# Plot

	orig_par=par(no.readonly=T);

	par(oma=c(1,1,1,1));

	table_sp=5;
	layout_mat=matrix(c(
		1,rep(2, table_sp),
		rep(c(3,rep(4, table_sp)), table_sp)),
		nrow=table_sp+1, byrow=T);
	#print(layout_mat);
	layout(layout_mat);

	# plot top/left spacer
	par(mar=c(0,0,0,0));
	plot(0,0,type="n", bty="n", xlab="", ylab="", main="", xaxt="n", yaxt="n");
	text(0,0, paste
		(nameb, ": ", bcuts, "\n x \n", namea, ": ", acuts, 
		"\n\nX^2 Test p-value:\n", sprintf("%1.3g", ct_cst_pval),
		"\n\nCumulative Top:\nPr(", nameb, "|", namea, ")=\n", round(cum_pr_bga, 3), 
		"\nPr(", namea, "|", nameb, ")=\n", round(cum_pr_agb, 3),
		 sep=""), cex=.8, font=2);



	# Scale leaf sample IDs
	label_scale=.2;
	dendra=dendrapply(dendra, text_scale_denfun);
	dendrb=dendrapply(dendrb, text_scale_denfun);
	
	# Color both dendrograms by A clustering
	sample_to_color_map=memb_byA;
	dendra=dendrapply(dendra, color_denfun_bySample);
	names(sample_to_color_map)=idsb;
	dendrb=dendrapply(dendrb, color_denfun_bySample);
	
	# Find height where clusters separate
	acutheight=find_height_at_k(hclA, acuts);
	bcutheight=find_height_at_k(hclB, bcuts);

	top_label_spc=4;
	left_label_spc=4;
	title_spc=2;

	# Plot A Dendrogram
	par(mar=c(5,left_label_spc,title_spc,0));
	plot(dendra, main=namea, horiz=F, yaxt="n", xaxt="n", xlab="", ylab="", xlim=c(-1,num_members+1));
	abline(h=acutheight, col="blue", lty=2, lwd=.7);
	abline(v=c(0,cumsum(grp_cnts_a)+.5), col="grey75", lwd=.5);
	trans_dend_mids_a=remap_coord(dend_mids_a, 0, num_members, 0, 1);

	# Plot B Dendrogram
	par(mar=c(0,title_spc,top_label_spc,5));
	plot(dendrb, main=nameb, horiz=T, xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(-1,num_members+1));
	abline(v=bcutheight, col="blue", lty=2, lwd=.7);
	abline(h=c(0, cumsum(grp_cnts_b)+.5), col="grey75", lwd=.5);
	trans_dend_mids_b=remap_coord(dend_mids_b, 0, num_members, 0, 1);

	# Plot shared statistics
	par(mar=c(0,left_label_spc,top_label_spc,0));
	plot(0,0, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n");
	#points(c(0,0,1,1), c(0,1,0,1));
	axis(3, at=trans_dend_mids_a, 1:acuts, tick=F, line=NA, font=2, cex.axis=2);
	axis(3, at=trans_dend_mids_a, grp_cnts_a, tick=F, line=-1, font=2, cex.axis=1);
	axis(2, at=trans_dend_mids_b, 1:bcuts, tick=F, line=NA, font=2, cex.axis=2);
	axis(2, at=trans_dend_mids_b, grp_cnts_b, tick=F, line=-1, font=2, cex.axis=1);

	cellab_size=min(1, 4/sqrt(acuts^2+bcuts^2));

	minpval=min(fish_exact_mat);
	cell_bounds_x=c(0,cumsum(grp_cnts_a)+.5)/num_members;
	cell_bounds_y=c(0,cumsum(grp_cnts_b)+.5)/num_members;

	for(colx in 1:acuts){
		for(rowx in 1:bcuts){

			cur_pval=fish_exact_mat[rowx, colx];

			cell_info=paste(
				"ob ct: ", ab_cnts_obs_mat[rowx, colx], "\n",
				"ob pr: ", round(ab_prop_obs_mat[rowx, colx], 3), "\n",
				"ex ct: ", round(ab_cnts_exp_mat[rowx, colx], 1), "\n",
				"ex pr: ", round(ab_prop_exp_mat[rowx, colx], 3), "\n",
				"fe pv: ", sprintf("%3.3g", cur_pval), "\n",
				"\n",
				"Pr(", namea, "|", nameb, ")=\n", sprintf("%3.3g", cnd_pr_agb[rowx, colx]), "\n",
				"Pr(", nameb, "|", namea, ")=\n", sprintf("%3.3g", cnd_pr_bga[rowx, colx]), "\n",
				sep="");

			
			# Highlight significant cells labels
			col="grey";
			font=1;
			if(cur_pval<.05){
				col="blue";
			}
			if(cur_pval<=minpval){
				font=2;
				points(c(cell_bounds_x[colx], cell_bounds_x[colx], cell_bounds_x[colx+1], 
					cell_bounds_x[colx+1], cell_bounds_x[colx]),
					c(cell_bounds_y[rowx], cell_bounds_y[rowx+1], cell_bounds_y[rowx+1], 
					cell_bounds_y[rowx], cell_bounds_y[rowx]), 
					col="cornflowerblue", type="l");
			}
				
			text(trans_dend_mids_a[colx], trans_dend_mids_b[rowx], cell_info, cex=cellab_size,
				font=font, col=col);
		}
	}

	par(orig_par);

	result=list();
	result[["chisqr_test_pval"]]=ct_cst_pval;
	result[["cumulative_cond_prob_agb"]]=cum_pr_agb;
	result[["cumulative_cond_prob_bga"]]=cum_pr_bga;

	return(result);

}

##############################################################################

plot_dendro_group_compare=function(hclA, hclB, acuts, bcuts, namea, nameb, idsb){

	color_denfun_bySample=function(n){
		if(is.leaf(n)){
			leaf_attr=attributes(n);
			leaf_name=leaf_attr$label;
			ind_color=sample_to_color_map[leaf_name];
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

	# Compute
	cat("Working on: ", nameb, ": ", bcuts, " x ", namea, ": ", acuts, "\n", sep="");

	dendra=as.dendrogram(hclA);
	dendrb=as.dendrogram(hclB);

	memb_byA=cutree(hclA, k=acuts);
	memb_byB=cutree(hclB, k=bcuts);

	dendr_names_a=get_clstrd_leaf_names(dendra);
	dendr_names_b=get_clstrd_leaf_names(dendrb);

	memb_byA=reorder_member_ids(memb_byA, dendr_names_a);
	memb_byB=reorder_member_ids(memb_byB, dendr_names_b);

	dend_mids_a=get_middle_of_groups(dendr_names_a, memb_byA);
	dend_mids_b=get_middle_of_groups(dendr_names_b, memb_byB);

	num_members=length(memb_byA);

	grp_cnts_a=(table(memb_byA)[1:acuts]);
	grp_cnts_b=(table(memb_byB)[1:bcuts]);;
	
	grp_prop_a=grp_cnts_a/num_members;
	grp_prop_b=grp_cnts_b/num_members;;

	cat("Group Proportions A:\n");
	print(grp_prop_a);

	cat("Group Proportions B:\n");
	print(grp_prop_b);

	##########################################

	orig_par=par(no.readonly=T);
	#par(mfrow=c(4,1));
	plsz=3;
	layout_mat=matrix(c(
		rep(1, plsz),
		rep(2, plsz),
		3,
		rep(4, plsz),
		rep(5, plsz)
		), byrow=T, ncol=1);
	layout(layout_mat);

	# Find height where clusters separate
	acutheight=find_height_at_k(hclA, acuts);
	bcutheight=find_height_at_k(hclB, bcuts);

	top_label_spc=4;
	left_label_spc=1;
	title_spc=2;

	#-----------------------------------------------------------------------------

	# Scale leaf sample IDs
	label_scale=.2;
	dendraA=dendrapply(dendra, text_scale_denfun);
	dendrbA=dendrapply(dendrb, text_scale_denfun);
	
	# Color both dendrograms by A clustering
	sample_to_color_map=memb_byA;
	dendraA=dendrapply(dendraA, color_denfun_bySample);
	idsa=names(sample_to_color_map);
	names(sample_to_color_map)=idsb;
	dendrbA=dendrapply(dendrbA, color_denfun_bySample);

	# Plot A Dendrogram
	par(mar=c(5,left_label_spc,title_spc,0));
	plot(dendraA, main=paste(namea, ": ", acuts, " cuts", sep=""), 
		horiz=F, yaxt="n", xaxt="n", xlab="", ylab="", xlim=c(-1,num_members+1));
	abline(h=acutheight, col="blue", lty=2, lwd=.7);
	abline(v=c(0,cumsum(grp_cnts_a)+.5), col="grey75", lwd=.5);

	# Plot B Dendrogram
	par(mar=c(5,left_label_spc,title_spc,0));
	plot(dendrbA, main=paste(nameb, ": colored by ", namea, sep=""), 
		horiz=F, xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(-1,num_members+1));
	abline(h=bcutheight, col="blue", lty=2, lwd=.7);
	abline(v=c(0, cumsum(grp_cnts_b)+.5), col="grey75", lwd=.5);

	#-----------------------------------------------------------------------------
	par(mar=c(0,0,0,0));
	plot(0,0, type="n", main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n");
	abline(h=0, col="grey", lwd=5);

	#-----------------------------------------------------------------------------

	# Color both dendrograms by B clustering
	label_scale=.2;
	dendraB=dendrapply(dendra, text_scale_denfun);
	dendrbB=dendrapply(dendrb, text_scale_denfun);

	sample_to_color_map=memb_byB;
	dendrbB=dendrapply(dendrbB, color_denfun_bySample);
	names(sample_to_color_map)=idsa;
	dendraB=dendrapply(dendraB, color_denfun_bySample);

	# Plot B Dendrogram
	par(mar=c(5,left_label_spc,title_spc,0));
	plot(dendrbB, main=paste(nameb, ": ", bcuts, " cuts", sep=""),
		horiz=F, xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(-1,num_members+1));
	abline(h=bcutheight, col="blue", lty=2, lwd=.7);
	abline(v=c(0, cumsum(grp_cnts_b)+.5), col="grey75", lwd=.5);

	# Plot A Dendrogram
	par(mar=c(5,left_label_spc,title_spc,0));
	plot(dendraB, main=paste(namea, ": colored by ", nameb, sep=""),
		horiz=F, yaxt="n", xaxt="n", xlab="", ylab="", xlim=c(-1,num_members+1));
	abline(h=acutheight, col="blue", lty=2, lwd=.7);
	abline(v=c(0,cumsum(grp_cnts_a)+.5), col="grey75", lwd=.5);

	par(orig_par);


}

##############################################################################

classic_mds_pts_A=matrix(0, nrow=num_samples, ncol=2); 
nonparm_mds_pts_A=matrix(0, nrow=num_samples, ncol=2); 
classic_mds_pts_B=matrix(0, nrow=num_samples, ncol=2); 
nonparm_mds_pts_B=matrix(0, nrow=num_samples, ncol=2); 

class_mdsA_res=cmdscale(dist_mat_A, k=2);
class_mdsB_res=cmdscale(dist_mat_B, k=2);
classic_mds_pts_A=class_mdsA_res;
classic_mds_pts_B=class_mdsB_res;

nonpr_mdsA_res=isoMDS(dist_mat_A, k=2);
nonpr_mdsB_res=isoMDS(dist_mat_B, k=2);
nonparm_mds_pts_A=nonpr_mdsA_res$points;
nonparm_mds_pts_B=nonpr_mdsB_res$points;

hcl_A=hclust(dist_mat_A, method="ward.D2");
hcl_B=hclust(dist_mat_B, method="ward.D2");

clus4_A=cutree(hcl_A, k=max_cuts);
clus4_B=cutree(hcl_B, k=max_cuts);

cat("Comparing MDS plots:\n");
par(mfrow=c(2,2));
compare_mds(classic_mds_pts_A, classic_mds_pts_B, "Classical MDS", 
	clus4_A, clus4_B, map_info[["a"]], map_info[["b"]]);
compare_mds(nonparm_mds_pts_A, nonparm_mds_pts_B, "NonMetric MDS", 
	clus4_A, clus4_B, map_info[["a"]], map_info[["b"]]);

##############################################################################

analyze_dendro_cont=T;

if(analyze_dendro_cont){

	probmat=matrix(0, nrow=max_cuts, ncol=max_cuts);
	rownames(probmat)=c(paste(map_info[["b"]], 1:max_cuts));
	colnames(probmat)=c(paste(map_info[["a"]], 1:max_cuts));

	pval_mat=probmat;
	cum_agb_mat=probmat;
	cum_bga_mat=probmat;

	for(acuts in 2:max_cuts){
		for(bcuts in 2:max_cuts){

			cont_res=plot_dendro_contigency(hcl_A, hcl_B, acuts, bcuts, 
				map_info[["a"]], map_info[["b"]], map_info[["b_id"]]);

			plot_dendro_group_compare(hcl_A, hcl_B, acuts, bcuts,
                                map_info[["a"]], map_info[["b"]], map_info[["b_id"]]);

			pval_mat[bcuts, acuts]=cont_res[["chisqr_test_pval"]];
			cum_agb_mat[bcuts, acuts]=cont_res[["cumulative_cond_prob_agb"]];
			cum_bga_mat[bcuts, acuts]=cont_res[["cumulative_cond_prob_bga"]];
		}
	}

	pval_mat=pval_mat[2:max_cuts, 2:max_cuts, drop=F];
	cum_agb_mat=cum_agb_mat[2:max_cuts, 2:max_cuts, drop=F];
	cum_bga_mat=cum_bga_mat[2:max_cuts, 2:max_cuts, drop=F];

	logodds=log2(cum_agb_mat/cum_bga_mat);


	print(pval_mat);
	min_cont_pval=min(pval_mat);
	min_idx=which(pval_mat==min_cont_pval, arr.ind=T);

	anames=colnames(pval_mat);
	bnames=rownames(pval_mat);

	paint_matrix(-log10(pval_mat), title="Contigency Table Dimension Log10(P-Values)");

	par(mfrow=c(1,1));
	plot_text(c(
		"Contingency Table Chi-Squared Tests by Num Clusters, p-value:",
		"",
		capture.output(print(signif(pval_mat, 3))),
		"",
		"",
		"",
		"Contingency Table Chi-Squared Tests by Num Clusters, -log10(p-value):",
		"",
		capture.output(print(-log10(pval_mat))),
		"",
		"",
		paste("Min P-Value: ", sprintf("%3.3g", min_cont_pval), 
			" at (", bnames[min_idx[1]], ", ", anames[min_idx[2]], ")", sep="")
	));

	plot_text(c(
		paste("Given a sample from ", map_info[["b"]], 
			" what's the probability classifying it in ",  map_info[["a"]], "?", sep=""),
		paste("Cumulative Pr(", map_info[["a"]], "|", map_info[["b"]], "):", sep=""),
		"",
		capture.output(print(round(cum_agb_mat,3))),
		"",
		"",
		paste("Given a sample from ", map_info[["a"]], 
			" what's the probability classifying it in ",  map_info[["b"]], "?", sep=""),
		paste("Cumulative Pr(", map_info[["b"]], "|", map_info[["a"]], "):", sep=""),
		"",
		capture.output(print(round(cum_bga_mat,3))),
		"",
		"",
		paste("Log(Pr(", map_info[["a"]], "|", map_info[["b"]], ")/Pr(", 
			map_info[["b"]], "|", map_info[["a"]], ")):", sep=""),
		paste("  Positive Log Prob Ratio implies ", map_info[["b"]], 
			" predicts ", map_info[["a"]], " better than vice versa.", sep=""),
		"",
		capture.output(print(round(logodds, 2)))
	));

	plot_dendro_contigency(hcl_A, hcl_B, min_idx[2]+1, min_idx[1]+1, 
		map_info[["a"]], map_info[["b"]], map_info[["b_id"]]);

	plot_dendro_group_compare(hcl_A, hcl_B, min_idx[2]+1, min_idx[1]+1,
                map_info[["a"]], map_info[["b"]], map_info[["b_id"]]);

	# Write optimal cuts to file
	cnt_fh=file(paste(OutputFileRoot, ".", map_info[["a"]], ".cuts", sep=""), "w");
	cat(file=cnt_fh, min_idx[2]+1, "\n", sep="");
	close(cnt_fh);

	cnt_fh=file(paste(OutputFileRoot, ".", map_info[["b"]], ".cuts", sep=""), "w");
	cat(file=cnt_fh, min_idx[1]+1, "\n", sep="");
	close(cnt_fh);
}


##############################################################################

cat("Comparing cluster separation:\n");
compare_pseudof(dist_mat_A, dist_mat_B, hcl_A, hcl_B,  max_k=max_cuts, map_info[["a"]], map_info[["b"]]);


##############################################################################

cat("Comparing distance distributions:\n");

par(mfrow=c(2,2));

hist(dist_mat_A, xlab="Distances", main=map_info[["a"]]);
hist(dist_mat_B, xlab="Distances", main=map_info[["b"]]);

if(!test){
	cat("Computing 2D histogram/Heatmap...\n");
	k=kde2d(dist_mat_A, dist_mat_B, n=500);
	image(k, col=rev(rainbow(100, start=0, end=4/6)), xlab=map_info[["a"]], ylab= map_info[["b"]], main=sprintf("Correlation: %3.3f", cor(dist_mat_A, dist_mat_B)));
}

max_dist_A=max(dist_mat_A);
max_dist_B=max(dist_mat_B);

plot(dist_mat_A, dist_mat_B, main="All Distances", xlab=map_info[["a"]], ylab= map_info[["b"]], cex=.5, xlim=c(0, max_dist_A), ylim=c(0, max_dist_B));

# Subsample so we can differentiate all the points in the scatter plot
num_dist=length(dist_mat_A);

##############################################################################

cat("Analyze sub-cluster statistics\n");

plot_subcluster_cuts=function(hclA, hclB, distmatA, distmatB, num_cuts, namea, nameb, idsb){

	color_denfun_bySample=function(n){
		if(is.leaf(n)){
			leaf_attr=attributes(n);
			leaf_name=leaf_attr$label;
			ind_color=sample_to_color_map[leaf_name];
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

	# Mean sum of squared distances 
	msd=function(dist){
		# Square matrix
		num_samp=ncol(dist);
		num_dist=(num_samp*(num_samp-1));
		#cat("Num Dist: ", num_dist, "\n");
		msd=sum(dist^2)/num_dist;
		return(msd);
	}

	# Spearman correlation
	distcorel=function(dista, distb){
		return(cor(dista, distb, method="spearman"));	
	}

	# Compute
	cat("Working on: ", namea, ": ", num_cuts, "\n", sep="");

	dendra=as.dendrogram(hclA);
	memb_byA=cutree(hclA, k=num_cuts);
	dendr_names_a=get_clstrd_leaf_names(dendra);
	memb_byA=reorder_member_ids(memb_byA, dendr_names_a);
	dend_mids_a=get_middle_of_groups(dendr_names_a, memb_byA);
	num_members=length(memb_byA);


	results=list();
	#msd_mat_a[clx, 1:num_cuts]=res$msdA;
	#msd_mat_b[clx, 1:num_cuts]=res$msdB;
	#cor_mat[clx, 1:num_cuts]=res$cor;

	memnamesa=names(memb_byA);
	memnamesb=idsb;

	msdA=numeric();
	msdB=numeric();
	man_corr=numeric();
	man_pval=numeric();

	##########################################
	# Plot

	orig_par=par(no.readonly=T);

	par(oma=c(1,1,1,1));

	den_spacing=2;
	hst_spacing=1;
	txt_spacing=1;
	cor_spacing=2;

	seqs=seq(0, num_cuts-1);

	layout_mat=matrix(c(
		rep(rep(1, num_cuts),den_spacing),
		rep((seqs*4)+2, hst_spacing),
		rep((seqs*4)+3, hst_spacing),
		rep((seqs*4)+4, txt_spacing),
		rep((seqs*4)+5, cor_spacing)
		),
		ncol=num_cuts, 
		byrow=T);
	#print(layout_mat);
	print(layout_mat);
	layout(layout_mat);

	# Scale leaf sample IDs
	label_scale=.2;
	dendra=dendrapply(dendra, text_scale_denfun);
	
	# Color both dendrograms by A clustering
	sample_to_color_map=memb_byA;
	dendra=dendrapply(dendra, color_denfun_bySample);
	
	# Find height where clusters separate
	acutheight=find_height_at_k(hclA, num_cuts);

	top_label_spc=4;
	left_label_spc=4;
	title_spc=2;

	# Plot A Dendrogram
	par(mar=c(5,left_label_spc,title_spc,0));
	plot(dendra, main=paste(namea, ", k=", num_cuts, sep=""), horiz=F, 
		yaxt="n", xaxt="n", xlab="", ylab="", xlim=c(-1,num_members+1));
	axis(side=1, at=dend_mids_a, labels=1:num_cuts, font=2, tick=F);
	abline(h=acutheight, col="blue", lty=2, lwd=.7);

	dist_max_a=max(distmatA);
	dist_max_b=max(distmatB);

	abreaks=seq(0, dist_max_a, length.out=15);
	bbreaks=seq(0, dist_max_b, length.out=15);

	for(clx in 1:num_cuts){
		membersA=memnamesa[memb_byA==clx];
		membersB=memnamesb[memb_byA==clx];

		num_members=length(membersA);
		cat("Cluster Cuts: ", clx, " of ", num_cuts, "\n");
		cat("Num A Members: ", length(membersA), "\n");
		cat("Num B Members: ", length(membersB), "\n");

		dista=distmatA[membersA, membersA, drop=F];
		distb=distmatB[membersB, membersB, drop=F];

		dist_arr_a=as.dist(dista);
		dist_arr_b=as.dist(distb);

		# Calcuate Ststs
		msdA[clx]=msd(dista);
		msdB[clx]=msd(distb);

		if(num_members>1){
			mantel_res=mantel(dista, distb, permutations=10000);
			man_corr[clx]=mantel_res[["statistic"]];
			man_pval[clx]=mantel_res[["signif"]];
			if(is.na(man_pval[clx])){
				man_pval[clx]=1;
			}	
		}else{
			man_corr[clx]=1;
			man_pval[clx]=1;
		}


		par(mar=c(.5,1,.5,0));
		hist(dist_arr_a, xlim=c(0, dist_max_a), xlab="", ylab="", main="", xaxt="n", yaxt="n", breaks=abreaks);
		if(clx==1){
			title(ylab=namea, line=0, cex.lab=1);
		}
		hist(dist_arr_b, xlim=c(0, dist_max_b), xlab="", ylab="", main="", xaxt="n", yaxt="n", breaks=bbreaks);
		if(clx==1){
			title(ylab=nameb, line=0, cex.lab=1);
		}

		par(mar=c(0,0,0,0));

		if(man_pval[clx]<.05){
			statcol="blue";
			statfon=2;
		}else{
			statcol="grey";	
			statfon=1;
		}
		plot(0,0, type="n", xlab="", ylab="", main="", xaxt="n", yaxt="n", bty="n");
		text(0,0, paste(
			c(
			paste("N=", num_members, sep=""),
			"",
			"Mantel:",
			sprintf("Cor: %2.3f", man_corr[clx]),
			sprintf("Pvl: %1.3g", man_pval[clx]),
			"",
			"MSD:",
			sprintf("%s: %1.2f", namea, msdA[clx]),
			sprintf("%s: %1.2f", nameb, msdB[clx])
		), collapse="\n"), cex=.7, col=statcol, font=statfon);

		par(mar=c(3,1,.5,0));

		# Draw regresson line for correlation
		plot(dist_arr_a, dist_arr_b, xaxt="n", yaxt="n", ylab="", main="", 
			xlim=c(0, dist_max_a), ylim=c(0, dist_max_b), col=clx);
		if(num_members>1){
			fit=lm(dist_arr_b~dist_arr_a);
			if(!is.na((fit$coefficients[2]))){
				abline(fit, col="grey", lty=2);
			}
		}
		title(xlab=clx, line=1.5, font.lab=2, cex.lab=2);
		title(xlab=namea, line=0, cex.lab=1);

		if(clx==1){
			title(ylab=nameb, line=0, cex.lab=1);
		}



	}

	results[["msdA"]]=msdA;
	results[["msdB"]]=msdB;
        results[["mantel_corr"]]=man_corr;
        results[["mantel_pval"]]=man_pval;

	par(orig_par);

	return(results);

}

analyze_subcluster_distances=function(hclA, hclB, distmatA, distmatB, num_cuts, namea, nameb, idsb){


	dist2da=as.matrix(distmatA);
	dist2db=as.matrix(distmatB);

	msd_mat_a=matrix(NA, nrow=num_cuts, ncol=num_cuts);
	msd_mat_b=matrix(NA, nrow=num_cuts, ncol=num_cuts);
	cor_mat=matrix(NA, nrow=num_cuts, ncol=num_cuts);
	cpv_mat=matrix(NA, nrow=num_cuts, ncol=num_cuts);

	for(clx in 1:num_cuts){
		cat("Analyzing ", clx, " cuts to: ", namea, "\n");

		res=plot_subcluster_cuts(hclA, clx, dist2da, dist2db, clx, namea, nameb, idsb);	

		msd_mat_a[clx, 1:clx]=res$msdA;
		msd_mat_b[clx, 1:clx]=res$msdB;
		cor_mat[clx, 1:clx]=res$mantel_corr;
		cpv_mat[clx, 1:clx]=res$mantel_pval;

		print(msd_mat_a);
		print(msd_mat_b);
		print(cor_mat);

	}

	# Plot optimal cluster cuts
	par(mfrow=c(2,1));
	maxmsd=max(msd_mat_a, na.rm=T);
	plot(0,0, type="n", xlim=c(1, num_cuts), ylim=c(0, maxmsd), 
		xlab="Cluster Cuts", ylab="MSD", main=paste(namea, " clustered by ", namea, sep=""), xaxt="n");
	axis(side=1, at=1:num_cuts, labels=T);
	abline(h=msd_mat_a[1,1], col="blue", lty=2);
	for(clx in 1:num_cuts){
		points(rep(clx, clx), msd_mat_a[clx, 1:clx], col=1:clx);	
		text(rep(clx, clx), msd_mat_a[clx, 1:clx], labels=1:clx, col="black", pos=4, cex=.7);
	}

	# Plot alternative cluster cuts
	maxmsd=max(msd_mat_b, na.rm=T);
	plot(0,0, type="n", xlim=c(1, num_cuts), ylim=c(0, maxmsd), 
		xlab="Cluster Cuts", ylab="MSD", main=paste(nameb, " clustered by ", namea, sep=""), xaxt="n");
	axis(side=1, at=1:num_cuts, labels=T);
	abline(h=msd_mat_b[1,1], col="blue", lty=2);
	for(clx in 1:num_cuts){
		points(rep(clx, clx), msd_mat_b[clx, 1:clx], col=1:clx);	
		text(rep(clx, clx), msd_mat_b[clx, 1:clx], labels=1:clx, col="black", pos=4, cex=.7);
	}

	# Plot cluster correlation across cuts
	cor_range=c(-1,1)*max(abs(cor_mat), na.rm=T);
	plot(0,0, type="n", xlim=c(1, num_cuts), ylim=cor_range, 
		xlab="Cluster Cuts", ylab="Distance Correlation", 
		main=paste("Distance Correlation: ", nameb, " clustered by ", namea, sep=""), xaxt="n");
	axis(side=1, at=1:num_cuts, labels=T);
	abline(h=cor_mat[1,1], col="blue", lty=2);
	abline(h=0, col="grey");

	for(clx in 1:num_cuts){
		points(rep(clx, clx), cor_mat[clx, 1:clx], col=1:clx);	

		sigchar=rep("", clx);
		sigchar[cpv_mat[clx,1:clx]<=0.05]="*";
		cl_label=paste(as.character(1:clx), sigchar, sep="");

		text(rep(clx, clx), cor_mat[clx, 1:clx], labels=cl_label, col="black", pos=4, cex=.7);
	}

	# Plot cluster correlation pvalues
	pval_ladder= c(1,.1,.05,.01,.001);
	nlog_cpv_mat=-log10(cpv_mat);
	nlog_pval_range=c(0, max(max(nlog_cpv_mat, na.rm=T), -log10(pval_ladder)));
	plot(0,0, type="n", xlim=c(1, num_cuts), ylim=nlog_pval_range, 
		xlab="Cluster Cuts", ylab="-log(pval)", 
		main=paste("Mantel P-Value: ", nameb, " clustered by ", namea, sep=""), xaxt="n");
	axis(side=1, at=1:num_cuts, labels=T);
	axis(side=4, at=-log10(pval_ladder), labels=pval_ladder, cex.axis=.7, tick=F, line=-.5, las=2);
	abline(h=-log10(pval_ladder), col="grey", lwd=.5, lty=2);
	for(clx in 1:num_cuts){
		points(rep(clx, clx), nlog_cpv_mat[clx, 1:clx], col=1:clx);
		text(rep(clx, clx), nlog_cpv_mat[clx, 1:clx], labels=1:clx, col="black", pos=4, cex=.7);
	}


}


analyze_subcluster_distances(hcl_A, hcl_B, dist_mat_A, dist_mat_B, max_cuts, map_info[["a"]],  map_info[["b"]], map_info[["b_id"]]);
analyze_subcluster_distances(hcl_B, hcl_A, dist_mat_B, dist_mat_A, max_cuts, map_info[["b"]],  map_info[["a"]], map_info[["a_id"]]);


##############################################################################

cat("Exporting shared summary tables...\n");

write_list=function(data, fname){
	fh=file(fname, "w");
	cat(file=fh, data, sep="\n");
	close(fh);
}

write_list(map_info[["a_id"]], paste(OutputFileRoot_nodist, ".", map_info[["a"]], ".shrd_w", map_info[["b"]], ".lst", sep=""));
write_list(map_info[["b_id"]], paste(OutputFileRoot_nodist, ".", map_info[["b"]], ".shrd_w", map_info[["a"]], ".lst", sep=""));

##############################################################################

cat("\nDone.\n")
dev.off();
warn=warnings();
if(length(warn)){
        print(warn);
}
q(status=0);
