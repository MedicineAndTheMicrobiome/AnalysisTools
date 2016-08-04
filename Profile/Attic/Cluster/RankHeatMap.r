#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"max_k", "k", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"\t\t-i <Input summary_table.xls FileName>\n\n",
	"\t\t[-k <Maximum k, default=10>]\n\n",
	"\n",
	"This program will read in a summary file table and generate a heatmap with ranks labeled on it.\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileNameRoot=gsub(".summary_table.xls", "", opt$input_file);

MaximumK=opt$max_k;
if(length(opt$max_k)==0){
	MaximumK=10;
}

# Figure out where we are running, this script so we can find the "WeightedRankDifference.r" code.
path_comp=strsplit(script_name, "/")[[1]];
bin_comp_idx=length(path_comp);
bin_path=paste(path_comp[-bin_comp_idx], collapse="/", sep="");
cat("Binary path: ", bin_path, "\n", sep="");
source(paste(bin_path, "WeightedRankDifference.r", sep="/"));

###############################################################################
# Load counts from file

cat("Working on ", InputFileName, "\n", sep="");

# Load summary_table.xls
mat<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, row.names=1))

# Exclude total counts column
count_mat=mat[,2:ncol(mat)]
num_cat=ncol(count_mat);
num_samples=nrow(count_mat);

if(num_samples<MaximumK){
	MaximumK=(num_samples-1);
}
cat("Maximum K set to: ", MaximumK, "\n", sep="");

#cat("Input matrix:\n");
#print(count_mat);

#------------------------------------------------------------------------------
# Get sample names and assign integer ids
sample_names=rownames(mat);
rownames(count_mat)=sprintf("%i",1:num_samples);

# Get column/category names
categories=as.vector(colnames(count_mat));
num_categories=length(categories);

#cat("Original category names:\n");
#print(categories);

# Compute shorted names
short_names=character(num_categories);
for(i in 1:num_categories){
	taxonomy=unlist(strsplit(categories[i], " "));
	short_names[i]=paste(taxonomy[1],taxonomy[length(taxonomy)],sep="_");
#	short_names[i]=taxonomy[length(taxonomy)];
}

#cat("Shortened category names:\n");
#print(short_names);

#------------------------------------------------------------------------------
# Identify 0 count columns
col_sum=apply(count_mat, 2, sum);
col_idx_iszero=(1:num_cat)[col_sum==0];
if(length(col_idx_iszero)>0){
	cat("Zero count categories:\n");
	print(col_idx_iszero);
}else{
	cat("No zero count categories detected.\n");
}

if(length(col_idx_iszero)>0){
	cat("Zero count names.\n");
	print(short_names[col_idx_iszero]);
}

# Remove all zero columns
if(length(col_idx_iszero)>0){
	cat("Removing zero count categories.\n");
	short_names=short_names[-col_idx_iszero];
	count_mat=count_mat[,-col_idx_iszero];
}

#------------------------------------------------------------------------------
# Show num samples/categories to be used

NumSamples=nrow(count_mat);
NumCategories=ncol(count_mat);

cat("\n");
cat("Num Samples: ", NumSamples, "\n");
cat("Num Categories: ", NumCategories, "\n");
cat("\n");

#------------------------------------------------------------------------------
# Normalize
sample_totals=numeric(NumSamples);
prob_mat=matrix(nrow=NumSamples,ncol=NumCategories);
colnames(prob_mat)=short_names;
rownames(prob_mat)=rownames(count_mat);
for(i in 1:NumSamples){
    	sample_totals[i]=sum(count_mat[i,]);
	prob_mat[i,]=count_mat[i,]/sample_totals[i];
}
#cat("Normalized matrix:\n");
#print(prob_mat);

#------------------------------------------------------------------------------
# Compute overall averages and sort by them
cat_averages=apply(prob_mat, 2, mean);
#print(cat_averages);
average_order=sort(cat_averages, decreasing=TRUE, index.return=TRUE);
prob_mat=prob_mat[,average_order$ix];
short_names_avgorderd=short_names[average_order$ix];

strA=as.integer(charToRaw("A"));
strZ=as.integer(charToRaw("Z"));
stra=as.integer(charToRaw("a"));
strz=as.integer(charToRaw("z"));
str0=as.integer(charToRaw("0"));
str9=as.integer(charToRaw("9"));

abbreviated=rawToChar(
	as.raw(
		c(
		    seq(strA,strZ,1),
		    seq(stra,strz,1),
		    seq(str0,str9,1)
		)
	),
	multiple=T
);

###############################################################################
# Output abbreviation to full name mapping

output_abbrev_mapping=function(abbr, full, file){
	num_names=length(full);
	for(i in 1:num_names){
		cat(abbr[i], "\t", full[i], "\n", file=file);
	}
}

abbrev_map_fh=file(paste(OutputFileNameRoot, ".abbrev_mapping", sep=""), "w");
output_abbrev_mapping(abbreviated, short_names_avgorderd, abbrev_map_fh);
close(abbrev_map_fh);

###############################################################################
###############################################################################

library(MASS);
library(stats);

###############################################################################

#source("WeightedRankDifference.r");

compute_topshared=function(hclust_res, ranks, k, topconsidered=10){
	cat("\tWorking on k = ", k, "\n", sep="");
	cuts=cutree(hclust_res, k=k);

	top_shared=list();
	for(i in 1:k){
		cluster_members=cuts[cuts==i];
		cluster_member_ids=as.numeric(names(cluster_members));
		cluster_ranks=ranks[cluster_member_ids,]<=topconsidered;
		sums=apply(as.matrix(cluster_ranks), 2, sum);
		num_members=length(cluster_members);
		shared_top=which(sums==num_members);

		struct=list();
		struct[["Samples"]]=cluster_member_ids;
		struct[["Categories"]]=shared_top;
		cs=paste("c",i,sep="");
		top_shared[[cs]]=struct;
	}

	return(top_shared);
}

output_cluster_drivers=function(top_shared, sample_names, category_names, file){
	keys=names(top_shared);
	num_clusters=length(keys);

	for(cl_id in 1:num_clusters){
		cat(file=file, "\tCluster(", cl_id, ")\n", sep="");
		#sample_str=paste(top_shared[[cl_id]][["Samples"]], collapse=",");
		#categories_str=paste(top_shared[[cl_id]][["Categories"]], collapse=",");

		sample_str=paste(sample_names[top_shared[[cl_id]][["Samples"]]], collapse=",");
		categories_str=paste(category_names[top_shared[[cl_id]][["Categories"]]], collapse=",");

		cat(file=file, "\t\t", sample_str, "\n");
		cat(file=file, "\t\t", categories_str, "\n");
	}
	cat(file=file, "\n");
}

output_drivers_rename_map=function(top_shared, sample_names, abbrev, file){
	keys=names(top_shared);
	num_clusters=length(keys);

	#sample_str=paste(top_shared[[cl_id]][["Samples"]], collapse=",");
	#categories_str=paste(top_shared[[cl_id]][["Categories"]], collapse=",");

	for(i in 1:num_clusters){
		categories_str=paste(abbrev[top_shared[[i]][["Categories"]]], collapse="");
		samples=sample_names[top_shared[[i]][["Samples"]]];

		for(samp_idx in 1:length(samples)){
			if(nchar(categories_str)>0){
				rename=paste(samples[samp_idx], "_", categories_str, sep="");
			}else{
				rename=samples[samp_idx];
			}
			cat(file=file, samples[samp_idx], "\t", rename, "\n", sep="");
		}

	}

}

###############################################################################


pdf(paste(OutputFileNameRoot, ".ranked_heatmap.pdf", sep=""), height=8.5, width=11);
drivers_full_fh=file(paste(OutputFileNameRoot, ".cluster_drivers_full", sep=""), "w");
drivers_abbrev_fh=file(paste(OutputFileNameRoot, ".cluster_drivers_abbrev", sep=""), "w");
rename_map_fh=file(paste(OutputFileNameRoot, ".sample_rename_map", sep=""), "w");

dist_list=list();
dist_list[[1]]=weight_rank_dist(prob_mat, 4);

dist_name=vector();
dist_name[1]="Weighted Rank Difference, power=4";

num_distances=length(dist_list);

inten_col=grey(rev(seq(0,1,length.out=100)));

for(dist_idx in 1:num_distances){

	cat("Working on: ", dist_name[dist_idx], "\n\n", sep="");

	samp_dist=dist_list[[dist_idx]];
	max_dist=max(samp_dist);

	# Remove 0 distances so isoMDS doesn't freakout
	for(i in 1:length(samp_dist)){
		if(samp_dist[i]==0){
			samp_dist[i]=1e-323;
		}
	}

	# Compute dendrogram
	hclust_res=hclust(samp_dist, method="ward");
	dendro=as.dendrogram(hclust_res);
	dend_order=(order.dendrogram(dendro));

	layout_mat=matrix(c(
		2,1,1,1,1,
		2,1,1,1,1,
		2,1,1,1,1,
		2,1,1,1,1,
		2,1,1,1,1
	), nrow=5, ncol=5, byrow=TRUE);
		
	layout(layout_mat);

	reordered_matrix=prob_mat[dend_order,];
	transposed_matrix=t(reordered_matrix);

	par(mar=c(5,0,2,0));
	par(oma=c(5,0,3,10));

	# Compute Rankings
	ranks=matrix(0, nrow=NumSamples, ncol=NumCategories);
	for(sampl_idx in 1:NumSamples){
		ranks[sampl_idx,]=(NumCategories-rank(reordered_matrix[sampl_idx,], ties.method="min"))+1;
	}

	# Plot the rank/heat map along different top rankings
	for(top in c(5, 10, 20, 30)){
	
		# Compute display range
		num_disp_cat=min(NumCategories, top);
		range=1:num_disp_cat;

		# Draw heatmap
		image(transposed_matrix[range,], col=inten_col, xaxt="n", yaxt="n", main=paste(dist_name[dist_idx], ", Top ", top, sep=""));

		# Label the rankings
		samp_pos=seq(0,1,length.out=NumSamples);
		cat_pos=seq(0,1,length.out=num_disp_cat);

		scale=35/max(NumSamples, num_disp_cat);
		for(sampl_idx in 1:NumSamples){
			for(cat_idx in 1:num_disp_cat){

				rank_val=ranks[sampl_idx, cat_idx];

				if(rank_val==1){
					color="red";
				}else if(rank_val==2){
					color="orange";
				}else if(rank_val==3){
					color="green";
				}else if(rank_val==4){
					color="blue";
				}else if(rank_val==5){
					color="purple";
				}else{
					color="black";
				}

				text(cat_pos[cat_idx], samp_pos[sampl_idx], label=rank_val, cex=scale, col=color);
			}
		}

		# Label sample names on right
		sample_name_scale=70/length(sample_names);
		if(sample_name_scale>1){
			sample_name_scale=1;
		}
		cat("Sample Name Scale: ", sample_name_scale, "\n", sep="");
		axis(side=4, at=samp_pos, sample_names[dend_order], las=2, cex.axis=sample_name_scale);

		# Label taxon on bottom
		category_name_scale=70/num_disp_cat;
		if(category_name_scale>1){
			category_name_scale=1;
		}
		axis(side=1, at=cat_pos, (short_names[average_order$ix])[range], las=2, cex.axis=category_name_scale);

		# Draw dendrogram on left
		# Try to estimate resizing the dendrogram so it matches the heatmap
		#plot(dendro, horiz=TRUE, leaflab="none", ylim=c(.75,NumSamples+.25)); # 7
		#plot(dendro, horiz=TRUE, leaflab="none", ylim=c(1.5,NumSamples-.5)); # 25
		m=(1.5-.75)/(25-7);
		m2=(-.5-.25)/(25-7);
		plot(dendro, horiz=TRUE, leaflab="none", xaxt="n", yaxt="n", ylim=c(m*NumSamples+.45833 ,NumSamples+(m2*NumSamples+.5416)));

	}

	# Output members and drivers of each cluster
	unorder=order(dend_order);
	cat(file=drivers_full_fh, dist_name[dist_idx], "\n\n", sep="");
	cat(file=drivers_abbrev_fh, dist_name[dist_idx], "\n\n", sep="");
	for(k in 1:MaximumK){
		top_shared=compute_topshared(hclust_res, ranks[unorder,], k);
		cat(file=drivers_full_fh, "*******************************************************************\n", sep="");
		cat(file=drivers_full_fh, "k = ", k, "\n", sep="");
		output_cluster_drivers(top_shared, sample_names, short_names_avgorderd, drivers_full_fh);

		cat(file=drivers_abbrev_fh, "*******************************************************************\n", sep="");
		cat(file=drivers_abbrev_fh, "k = ", k, "\n", sep="");
		output_cluster_drivers(top_shared, sample_names, abbreviated, drivers_abbrev_fh);

		if(k==MaximumK){
			output_drivers_rename_map(top_shared, sample_names, abbreviated, rename_map_fh);
		}
	}
}





cat("Done.\n");
print(warnings());
q(status=0);
