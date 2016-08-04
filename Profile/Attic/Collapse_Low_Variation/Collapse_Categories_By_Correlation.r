#!/usr/bin/env Rscript

###############################################################################
#                                                                             #
#       Copyright (c) 2013 J. Craig Venter Institute.                         #
#       All rights reserved.                                                  #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                             #
###############################################################################
###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"correlation_cutoff", "r", 2, "numeric",
	"output_fname_root", "o", 2, "character",
	"testing_flag", "t", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

CorrelCutoff=.8;

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	[-r <correlation cutoff, default =", CorrelCutoff, ">]\n",
	"	[-o <output filename root>]\n",
	"	[-t (Testing flag)\n",
	"\n",
	"Reads in the summary table and computes correlations between\n",
	"all categories across the samples in the summary table.\n",
	"Abundances are first transformed into log(odds) to make their distribution\n",
	"more normal before computing correlations.\n",
	"The correlations are then converted into distances, d = 1-|cor|, so\n",
	"that the greater the magnitude of the correlation, the closer\n",
	"two categories are.   Hierarchical clustering is then performed\n",
	"using complete linkage, so that for any cluster, the minimum required\n",
	"correlation between two categories will be enforced when cutree is\n",
	"performed.  For each cluster, a reference category is selected based\n",
	"on the member with the larger abundance.  If the cluster contains\n",
	"positively and negatively correlated members, it is split into two\n",
	"clusters.  The new cluster is then named after the reference cluster\n",
	"prefixed with a P if the members are positively correlated with the\n",
	"reference, or else N if the members are negatively correlated.\n",
	"\n",
	"The following output is generated:\n",
	"	1.) Summary table with collapsed members.\n",
	"	2.) Dendrogram containing the categories, with cutoff marked.\n",
	"	3.) A list of cluster and their members, including abundances and correlations.\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

if(length(opt$correlation_cutoff)){
	CorrelCutoff=opt$correlation_cutoff;
}

if(!length(opt$output_fname_root)){
    	OutputFnameRoot=gsub(".summary_table.xls", "", opt$input_file);
    	OutputFnameRoot=gsub(".summary_table.tsv", "", OutputFnameRoot);
}else{
    	OutputFnameRoot=opt$output_fname_root;
}

OutputFnameRoot=paste(OutputFnameRoot, sprintf(".%g", CorrelCutoff*100), sep="");

TestFlag=F;
if(length(opt$testing_flag)){
	TestFlag=T;
}

###############################################################################

InputFileName=opt$input_file;

cat("\n");
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name Root: ", OutputFnameRoot, "\n");
cat("\n");
cat("Correlation cutoff: ", CorrelCutoff, "\n");
cat("\n");

###############################################################################
###############################################################################

compute_correlation_matrix=function(data){

	nvar=ncol(data);
	names=colnames(data);
	cat("Num categories to compare: ", nvar, "\n");

	# Matrices to store pvalues and correlation values
	cor_mat=matrix(1, nrow=nvar, ncol=nvar);
	cor_pval_mat=matrix(NA, nrow=nvar, ncol=nvar);
	adj_cor_pval_mat=matrix(NA, nrow=nvar, ncol=nvar);
	triangle=numeric(nvar*(nvar-1)/2)

	colnames(cor_mat)=names;
	rownames(cor_mat)=names;
	colnames(cor_pval_mat)=names;
	rownames(cor_pval_mat)=names;
	colnames(adj_cor_pval_mat)=names;
	rownames(adj_cor_pval_mat)=names;

	# Compute correlation and p-values
	k=1;
	for(i in 1:nvar){
		for(j in 1:i){
			cor_res=cor.test(data[,i], data[,j]);
			cor_pval_mat[i,j]=cor_res$p.value;
			cor_mat[i,j]=cor_res$estimate;
			triangle[k]=cor_res$p.value;
			k=k+1;

			# Populate other side
			cor_pval_mat[j,i]=cor_pval_mat[i,j];
			cor_mat[j,i]=cor_mat[i,j]
		}
	}

	# Adjust for multiple testing
	triangle_adjusted=p.adjust(triangle, method="fdr");

	# Move adjust pvalues back into matrix 
	k=1;
	for(i in 1:nvar){
		for(j in 1:i){
			adj_cor_pval_mat[i,j]=triangle_adjusted[k];
			adj_cor_pval_mat[j,i]=triangle_adjusted[k];
			k=k+1;
		}
	}

	# Save results in structure
	res=list();
	res$cor=cor_mat;
	res$unadj.pvalues=cor_pval_mat;
	res$adj.pvalues=adj_cor_pval_mat;

	return(res);
}

###############################################################################
###############################################################################

###############################################################################
# Load data
inmat=as.matrix(read.delim(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="", quote="", row.names=1))
#cat("Original Matrix:\n")
#print(inmat);

# Grab columns we need into a vector, ignore totals, we won't trust it.
counts_mat=inmat[,2:(ncol(inmat)), drop=F];
#print(counts_mat);

num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);

category_names=colnames(counts_mat);
sample_names=rownames(counts_mat);

# Report data summary
num_categories=ncol(counts_mat);
num_samples=nrow(counts_mat);
cat("Number of Categories: ", num_categories, "\n");
cat("Number of Samples: ", num_samples, "\n");
cat("\n");

if(length(unique(category_names))!=num_categories){
	cat("Error:  Category names are not unique.\n");
	quit(status=-1);
}else{
	cat("Great, categories names are indeed unique.\n");
}

###############################################################################
# Normalize

# Sum sample totals
sample_totals=apply(counts_mat, 1, sum);

# normalize, to compute probabilities
normalized=matrix(0, nrow=num_samples, ncol=num_categories);
colnames(normalized)=category_names;
for(i in 1:num_samples){
	normalized[i,]=counts_mat[i,]/sample_totals[i];
}
#print(normalized);

# Set minimum abundance to be 1/10 of min
min=min(normalized[normalized!=0]);
cat("Min non-zero abundance:", min, "\n");

zerosubstitute=min/10;
cat("Setting zeros to: ", zerosubstitute, "\n");

# Generate abundances assuming 0 are 1/10 minimum
norm_no_zeros=matrix(NA, ncol=num_categories, nrow=num_samples);
for(j in 1:num_categories){
	for(i in 1:num_samples){
		if(normalized[i,j]==0){
			norm_no_zeros[i,j]=zerosubstitute;
		}else{
			norm_no_zeros[i,j]=normalized[i,j];
		}
	}
}

# Compute average abundance across all samples
avg_abundance=apply(normalized, 2, mean);
names(avg_abundance)=colnames(normalized);

###############################################################################
# Log(Odds)

# Compute the logodds so transformed abundances look more normal
logodds=function(x){
	return(log10(x/(1-x)));
}

# Generate the log(odds) matrix
logodds_mat=logodds(norm_no_zeros);
colnames(logodds_mat)=category_names;

# Compute range of the log(odds)
logodds_range=range(logodds_mat);
cat("Range Log10(Odds): (", logodds_range[1], ",  ", logodds_range[2], ") \n", sep="");

logodds_width=abs(diff(logodds_range));

###############################################################################
# Compute correlation on Log(odds)

cor_mat=cor(logodds_mat);
cor_mat[is.na(cor_mat)]=0;
colnames(cor_mat)=category_names;
rownames(cor_mat)=category_names;
num_entries=ncol(cor_mat);

cat("Correlation Matrix Dimensions: ", num_entries, "\n");

# Treat high pos or negative correlation as similar
dist=1-abs(cor_mat);
hcl=hclust(as.dist(dist));
den=as.dendrogram(hcl);

#-------------------------------------------------------------------------------

SkipDendrogram=F;
if(!SkipDendrogram){

	if(TestFlag){
		rnd=sprintf(".%03g", ceiling(runif(1, 0, 1000)));
	}else{
		rnd="";
	}
	pdf(paste(OutputFnameRoot, ".dendr",  rnd, ".pdf", sep=""), height=11*17, width=8.5);
	par(mar=c(5.1, 4.1, 4.1, 10));

	den_resize_labels=function(n, size){
		if(is.leaf(n)){
			leaf_attr=attributes(n);
			leaf_name=leaf_attr$label;
			attr(n, "nodePar") = c(leaf_attr$nodePar,
						cex=0,
						lab.cex=size);
		}
		return(n);
	}

	label_scale=min(1, 800/num_entries);
	cat("Dendrogram label scale: ", label_scale, "\n");
	den=dendrapply(den, den_resize_labels, label_scale);

	#-------------------------------------------------------------------------------

	ylimit=c(1-5/num_entries, num_entries+5/num_entries);

	tryCatch(
	{
		plot(den, horiz=T, ylim=ylimit, yaxt="n", ylab="abs(correlation)");
		axis(side=1, at=seq(0,1,.05), labels=sprintf("%3.2f", 1-seq(0,1,.05)), las=2);
		axis(side=3, at=seq(0,1,.05), labels=sprintf("%3.2f", 1-seq(0,1,.05)), las=2);
		abline(v=.05, col="grey", lwd=.25);
		abline(v=.1, col="grey", lwd=.25);
		abline(v=.15, col="grey", lwd=.25);
		abline(v=.2, col="grey", lwd=.25);
		abline(v=(1-CorrelCutoff), col="blue", lty=2, lwd=1.5);
	}, error=function(e){
		cat("*******************************************\n");
		cat("*                                         *\n");
		cat("*  Error plotting dendrogram.  Skipping.  *\n");
		cat("*  (Maybe out of memory?)                 *\n");
		cat("*                                         *\n");
		cat("*******************************************\n");
	}
	);

	d=dev.off();
}

###############################################################################
# Cut tree

clusters=cutree(hcl, h=(1-CorrelCutoff));
cluster_counts=table(clusters);
clusters_gt1_idx=which(cluster_counts>1);
num_clusters_gt1=length(clusters_gt1_idx);
if(num_clusters_gt1==0){
	cat("No clusters with greater than a single member.  Try a lower correlation cutoff.\n");
	quit(status=-1);
}

singleton_idx=which(cluster_counts==1);
num_singletons=length(singleton_idx);

cat("Number of clusters with more than 1 member: ", num_clusters_gt1, "\n");
cat("Number of singleton clusters: ", num_singletons, "\n");

###############################################################################
# Collapse counts and Output cluster members

rp_fh=file(paste(OutputFnameRoot, ".cor_col.cluster_reps.tsv", sep=""), "w");
cat(file=rp_fh, paste(c("", "Category", "Avg Abundance", "Correlation"), collapse="\t"), "\n", sep="");

new_summary_table=matrix(NA, nrow=num_samples, ncol=num_categories);
new_st_colnames=character(num_categories);

num_new_categories=0;
prefix=c("P", "N");
which_cor=list();

cluster_list=list();
for(cl_ix in clusters_gt1_idx){
	cur_cl=clusters[clusters==cl_ix];
	clnames_arr=names(cur_cl);

	# Select member with greatest abundance as representative
	cl_abun=avg_abundance[clnames_arr];
	rep_cat=which(max(cl_abun)==cl_abun)[1];
	rep_name=clnames_arr[rep_cat];

	# Extract out correlations specific to this cluster
	cl_cor=cor_mat[rep_name, clnames_arr];

	# Split correlations into pos and neg
	which_cor[[1]]=which(cl_cor>0);
	which_cor[[2]]=which(cl_cor<0)

	for(negpos_ix in c(1,2)){
		np=which_cor[[negpos_ix]];
		if(length(np)>0){

			# Generate new representative name
			output_rep_name=paste(prefix[negpos_ix], ";", rep_name, sep="");

			# Get category names for all cluster members
			member_names=clnames_arr[np];

			# Output member information
			cat(file=rp_fh, "Representative:\t", output_rep_name, "\n", sep="");
			for(cat_name in member_names){

				flag="";
				if(cat_name==rep_name){
					flag="*";
				}
				outstr=paste(
					flag, cat_name, avg_abundance[cat_name], cor_mat[rep_name, cat_name],
					sep="\t");
				cat(file=rp_fh, outstr, "\n", sep="");
			}
			cat(file=rp_fh, "\n");

			# Keep track of number of collapsed/new categories
			num_new_categories=num_new_categories+1;

			# Sum up the counts for this cluster
			sums=apply(counts_mat[, member_names, drop=F], 1, sum);
			new_summary_table[, num_new_categories]=sums;
			new_st_colnames[num_new_categories]=output_rep_name;

		}
	}


}

# Include singletons back into new table
for(cl_ix in singleton_idx){
	cl_name=names(clusters[clusters==cl_ix]);
	num_new_categories=num_new_categories+1;
	new_summary_table[, num_new_categories]=counts_mat[, cl_name];
	new_st_colnames[num_new_categories]=cl_name;

	# Output singleton is cluster rep file
	outstr=paste("Representative:\t", cl_name, sep="");
	cat(file=rp_fh, outstr, "\n", sep="");
	outstr=paste(c("*", cl_name, avg_abundance[cl_name], cor_mat[cl_name, cl_name]), collapse="\t");
	cat(file=rp_fh, outstr, "\n\n", sep="");
}

# Remove unused columns (because they were collapsed)
output_summary_table=new_summary_table[, 1:num_new_categories];
colnames(output_summary_table)=new_st_colnames[1:num_new_categories];
rownames(output_summary_table)=sample_names;

# Confirm that counts still add up
new_sample_counts=apply(output_summary_table, 1, sum);

almost_same=function(x, y){
	diff=x-y;
	return(abs(diff)<1e-10);
}

if(!all(almost_same(new_sample_counts, sample_totals))){
	cat("Error:  Collapsed sample counts don't equal original sample counts.\n");
	quit(status==-1);
}

###############################################################################

# Write out summary file table

write_summary_table=function(counts_mat, fname){

	stfname=paste(fname, ".summary_table.tsv", sep="");
	cat("Writing summary file table: ", stfname, "\n", sep="");

	st_fh=file(stfname, "w");

	num_samples=nrow(counts_mat);
	sample_names=rownames(counts_mat);
	category_names=colnames(counts_mat);

	# Output Header
	cat(file=st_fh, paste("sample_id", "total", paste(category_names, collapse="\t"), sep="\t"));
	cat(file=st_fh, "\n");

	# Output rows
	for(samp_idx in 1:num_samples){
		total=sum(counts_mat[samp_idx,]);
		outline=paste(sample_names[samp_idx],total,paste(counts_mat[samp_idx,], collapse="\t"), sep="\t");
		cat(file=st_fh, outline, "\n", sep="");
	}

	close(st_fh);
	cat("Done writing.\n");
}

write_summary_table(output_summary_table, paste(OutputFnameRoot, ".cor_clps", sep=""));


###############################################################################

cat("Done.\n")
print(warnings());
q(status=0)
