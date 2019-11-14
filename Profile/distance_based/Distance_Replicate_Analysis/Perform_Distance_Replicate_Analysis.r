#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
library('vegan');

options(useFancyQuotes=F);
options(width=200);

DEF_DISTTYPE="man";

source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');

params=c(
	"summary_file", "s", 1, "character",
	"factor_file", "f", 1, "character",
	"sample_coln", "S", 1, "character",
	"replicate_coln", "R", 1, "character",
	"outputroot", "o", 1, "character",
	"distance_type", "d", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table>\n",
	"	-f <factors (replicate map)>\n",
	"	-S <Sample Name column name>\n",
	"	-R <Replicate Name column name>\n",
	"	-o <output filename root>\n",
	" 	[-d <distance metric, default=", DEF_DISTTYPE, ">]\n",
	"\n",
	"This script will perform an analyses across a set of\n",
	"of replicates using bootstrapping to look at the effects\n",
	"of sequencing depth.\n",
	"\n",
	"Only a few columns are used in the factor file:\n",
	"	1.) Sample ID (to look up in the summary table)\n",
	"	2.) Sample Name (Underlying sample being replicated)\n",
	"	3.) Replicate Name (e.g. Run ID)\n",
	"\n");

if(
	!length(opt$summary_file) || 
	!length(opt$factor_file) || 
	!length(opt$sample_coln) || 
	!length(opt$replicate_coln) ||
	!length(opt$outputroot)
){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputRoot=gsub(".summary_table.xls", "", opt$summary_file);
	OutputRoot=gsub(".summary_table.tsv", "", OutputRoot);
}else{
	OutputRoot=opt$outputroot;
}


if(!length(opt$distance_type)){
	DistanceType=DEF_DISTTYPE;
}else{
	DistanceType=opt$distance_type;
}

SummaryFile=opt$summary_file;
FactorsFile=opt$factor_file;
OutputRoot=opt$outputroot;
Replicate_ColumnName=opt$replicate_coln;
Sample_ColumnName=opt$sample_coln;


cat("\n");
cat("Summary File : ", SummaryFile, "\n", sep="");
cat("Factors File: ", FactorsFile, "\n", sep="");
cat("Output File: ", OutputRoot, "\n", sep="");
cat("Replicate ID Column Name: ", Replicate_ColumnName, "\n", sep="");
cat("Sample Name Column Name: ", Sample_ColumnName, "\n", sep="");
cat("Distance Type: ", DistanceType, "\n", sep="");
cat("\n");


##############################################################################

load_factors=function(fname){
	cat("Loading Factors: ", fname, "\n");
	factors=data.frame(read.table(fname, sep="\t",  header=TRUE, row.names=1, 
		check.names=FALSE, comment.char="", quote=""));
	factor_names=colnames(factors);

	ignore_idx=grep("^IGNORE\\.", factor_names);

	if(length(ignore_idx)!=0){
		return(factors[-ignore_idx]);
	}else{
		return(factors);
	}
}

load_summary_file=function(fname){
	cat("Loading Summary Table: ", fname, "\n");
	inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, 
		comment.char="", row.names=1))
	counts_mat=inmat[,2:(ncol(inmat))];
	return(counts_mat);
}

normalize=function(counts){
	totals=apply(counts, 1, sum);
	num_samples=nrow(counts);
	normalized=matrix(0, nrow=nrow(counts), ncol=ncol(counts));

	for(i in 1:num_samples){
		normalized[i,]=counts[i,]/totals[i];
	}
	
	colnames(normalized)=colnames(counts);
	rownames(normalized)=rownames(counts);	
	return(normalized);
}

plot_text=function(strings){
	par(mfrow=c(1,1));
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
}

##############################################################################

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


##############################################################################

# Load matrix
counts=load_summary_file(SummaryFile);
num_taxa=ncol(counts);
num_samples=nrow(counts);
#print(counts);

##############################################################################

# Load factors
factors=load_factors(FactorsFile);
factor_names=colnames(factors);
num_factors=ncol(factors);
factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);

cat("\n");
cat(num_factors, " Factor(s) Loaded:\n", sep="");
#print(factor_names);
cat("\n");

##############################################################################

# Reconcile factors with samples
factor_sample_ids=rownames(factors);
counts_sample_ids=rownames(counts);

shared_sample_ids=intersect(factor_sample_ids, counts_sample_ids);
num_shared_sample_ids=length(shared_sample_ids);
num_factor_sample_ids=length(factor_sample_ids);
num_counts_sample_ids=length(counts_sample_ids);

cat("Num counts sample IDs: ", num_counts_sample_ids, "\n");
cat("Num factor sample IDs: ", num_factor_sample_ids, "\n");
cat("Num shared sample IDs: ", num_shared_sample_ids, "\n");
cat("\n");

cat("Samples missing from count information:\n");
print(setdiff(factor_sample_ids, counts_sample_ids));
cat("\n");
cat("Samples missing from factor information:\n");
print(setdiff(counts_sample_ids, factor_sample_ids));
cat("\n");
cat("Total samples shared: ", num_shared_sample_ids, "\n");

shared_sample_ids=sort(shared_sample_ids);

##############################################################################

# Reorder data by sample id
counts=counts[shared_sample_ids,];
num_samples=nrow(counts);
factors=factors[shared_sample_ids,, drop=F];

normalized=normalize(counts);

##############################################################################

pdf(paste(OutputRoot, ".rep_analys.dist.pdf", sep=""), height=14, width=8.5);

plot_text(c(
	paste("Summary File : ", SummaryFile, sep=""),
	paste("Factors File: ", FactorsFile, sep=""),
	paste("Output File: ", OutputRoot, sep=""),
	paste("Replicate ID Column Name: ", Replicate_ColumnName, "\n", sep=""),
	paste("Sample Name Column Name: ", Sample_ColumnName, "\n", sep=""),
	paste("Distance Type: ", DistanceType, "\n", sep="")
));

##############################################################################
# For each sample name

replicate_group=factors[, Replicate_ColumnName];
sample_names=factors[, Sample_ColumnName];
sample_ids=rownames(factors);

uniq_samp_names=unique(sample_names);
uniq_repl_names=unique(replicate_group);

num_uniq_samp=length(uniq_samp_names);
num_uniq_repl=length(uniq_repl_names);

cat("Number of Unique Sample Names: ", num_uniq_samp, "\n");
cat("\n");
cat("Number of Unique Replicate Group: ", num_uniq_repl, "\n");
print(uniq_repl_names);

##############################################################################

signf_char=function(x){
	if(x<.001){
		return("***");
	}else if(x<.01){
		return("**");
	}else if(x<.05){
		return("*");
	}else if(x<.1){
		return(".");
	}else{
		return("");
	}
}

##############################################################################

bootstrap_profiles=function(depths_arr, norm_arr, num_bs){

	num_samps=length(depths_arr);
	samp_names=names(depths_arr);

	total_bs=num_bs*num_samps;
	bs_dm=matrix(0, nrow=total_bs, ncol=ncol(norm_arr));
	bs_dm_rownames=character(total_bs);

	for(samp_ix in 1:num_samps){
		row_ix=(1:num_bs)+num_bs*(samp_ix-1);

		bs_dm[row_ix,]=
			t(rmultinom(num_bs, depths_arr[samp_ix], norm_arr[samp_ix,]))/depths_arr[samp_ix];

		bs_dm_rownames[row_ix]=sprintf("%s.%03g", samp_names[samp_ix], 1:num_bs);
		
		samp_ix=samp_ix+1;
	}
	
	print(bs_dm_rownames);
	rownames(bs_dm)=bs_dm_rownames;

	return(bs_dm);

} 




diversity_95ci_rarefaction=function(depth_range, distribution, div_fun, num_bs=80){

        num_depths=length(depth_range);
        ci95_matrix=matrix(NA, ncol=3, nrow=num_depths);
        colnames(ci95_matrix)=c("lb95", "median", "ub95");

        dpix=1;
        bs_div_vals=numeric(num_bs);
        for(dep in depth_range){

                bs_div_vals=diversity_bs(dep, distribution, div_fun, num_bs);
                ci95_matrix[dpix,]=quantile(bs_div_vals, c(.025, .5, .975));

                dpix=dpix+1;
        }

        raref_results=list();
        raref_results[["ci95_matrix"]]=ci95_matrix;
        raref_results[["div_at_max"]]=bs_div_vals;
        # Store bs values at greatest depth

        return(raref_results);

}



##############################################################################

nbs=80;
iter=0;

for(cur_sample_name in uniq_samp_names){

	cat("Analyzing Sample Group: ", cur_sample_name, "\n");
	
	samp_ix=(sample_names==cur_sample_name);
	samp_ids=sample_ids[samp_ix];
	repl_ids=replicate_group[samp_ix];
	names(repl_ids)=samp_ids;

	num_samp=length(samp_ids);
	if(num_samp<2){
		next;
	}

	colors=rainbow(num_samp, start=0, end=.8);
	names(colors)=samp_ids;

	cat("Samples ID in Sample Group:\n");
	print(samp_ids);
	cur_counts=counts[samp_ids,, drop=F];
	cur_normal=normalized[samp_ids,, drop=F];
	cur_centroid=apply(cur_normal, 2, mean);
	combined_norm=rbind(cur_normal, cur_centroid);

	read_depths=apply(cur_counts, 1, sum);
	print(read_depths);
	max_depth=max(read_depths);

	cat("Bootstraping Profile Distribution...\n");
	bs_profs=bootstrap_profiles(read_depths, cur_normal, num_bs=nbs);

	cat("Calculating Distributions...\n");
	distmat=compute_dist(bs_profs, DistanceType);

	cat("Performing MDS...\n");
	clmds=cmdscale(distmat);


	par(oma=c(1,1,4,1));
        par(mfrow=c(2,1));
        par(mar=c(5,20,5,1));
        barplot(read_depths, horiz=T, las=1, col=colors, xlab="Reads/Sample", main="Replicate Sequencing Depth");


	par(mar=c(10,4,10,4));
	rang=range(clmds);

	bscol=as.vector(matrix(rep(colors,nbs), byrow=T, nrow=nbs));
	plot(clmds[,1], clmds[,2], xlim=rang, ylim=rang, col=bscol);

	cat("ok.\n");


	#gap=150;
	#if(max_depth<gap){
#		gap=max_depth/10;
	#}
	#depth_range=sort(unique(c(seq(gap, max_depth, gap), read_depths)));
	#cat("Rarefying at: \n");
	#print(depth_range);

	#cat("\n");

	iter=iter+1;
	if(iter==20){
		break;
	}
}


##############################################################################

cat("Done.\n");
print(warnings());
q(status=0);
