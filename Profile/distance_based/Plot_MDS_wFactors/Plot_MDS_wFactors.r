#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"factor_file", "f", 1, "character",
	"factor_subset", "M", 2, "character",
	"output_file", "o", 2, "character",
	"distance_type", "d", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

DEF_DISTANCE="euc";
TOP_CATEGORIES=20;
NUM_ROWS=2;
NUM_COLS=3;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	-f <factor file>\n",
	"	[-M <subset of factors/variables to group by>]\n",
	"	[-o <output file root name>]\n",
	"	[-d <distance, default=", DEF_DISTANCE, ".]\n",
	"\n",
	"	This script will read in the summary table\n",
	"	and the factor file.\n",
	"\n",
	"	For each factor in the factor/metadata file, an MDS plot will be generated\n",
	"	with the samples grouped by color.\n",
	"\n",
	"	Distance Types: euc,wrd,man,bray,horn,bin,gow,tyc,minkp5,minkp3\n",
	"\n",
	"\n", sep="");

if(!length(opt$input_file) || !length(opt$factor_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
FactorFileName=opt$factor_file;

FactorSubset=NULL;
if(length(opt$factor_subset)){
	FactorSubset=opt$factor_subset;
}

DistanceType=DEF_DISTANCE;
if(length(opt$distance_type)){
	DistanceType=opt$distance_type;
}

if(length(opt$output_file)>0){
	OutputFileRoot=opt$output_file;
}else{
	OutputFileRoot=InputFileName;
	OutputFileRoot=gsub("\\.summary_table\\.tsv$", "", OutputFileRoot);
	OutputFileRoot=gsub("\\.summary_table\\.xls$", "", OutputFileRoot);
	cat("No output file root specified.  Using input file name as root.\n");
}

###############################################################################

OutputFileRoot=paste(OutputFileRoot, ".", substr(DistanceType, 1, 4), sep="");
OutputPDF = paste(OutputFileRoot, ".mds.pdf", sep="");
cat("Output PDF file name: ", OutputPDF, "\n", sep="");

inch_p_plot=3;
pdf(OutputPDF,width=inch_p_plot*(3+.5), height=inch_p_plot*3)

###############################################################################

load_factors=function(fname){
        factors=data.frame(read.table(fname,  header=TRUE, check.names=FALSE, 
		row.names=1, comment.char="", quote="", sep="\t"));
        dimen=dim(factors);
        cat("Rows Loaded: ", dimen[1], "\n");
        cat("Cols Loaded: ", dimen[2], "\n");
        return(factors);
}

load_summary_file=function(fname){
        cat("Loading Summary Table: ", fname, "\n");
        inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", row.names=1))
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

	orig.par=par(no.readonly=T);
        par(family="Courier");
        par(oma=rep(.1,4));
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

	par(orig.par);
}

###############################################################################

tail_statistic=function(x){
        sorted=sort(x, decreasing=TRUE);
        norm=sorted/sum(x);
        n=length(norm);
        tail=0;
        for(i in 1:n){
                tail=tail + norm[i]*((i-1)^2);
        }
        return(sqrt(tail));
}

###############################################################################

orig_factors_mat=load_factors(FactorFileName);
#print(factors_mat);

###############################################################################

orig_counts_mat=load_summary_file(InputFileName);
#print(counts_mat);

###############################################################################

orig_factors_samples=rownames(orig_factors_mat);
orig_counts_samples=rownames(orig_counts_mat);
shared=intersect(orig_factors_samples, orig_counts_samples);

cat("\n\n");
cat("Samples not represented in summary table file:\n");
excl_to_st=setdiff(orig_counts_samples, shared);
print(excl_to_st);
cat("Samples not represented in offsets file:\n");
excl_to_fct=setdiff(orig_factors_samples, shared);
print(excl_to_fct);
cat("\n\n");

num_shared=length(shared);
cat("Number of Shared Samples: ", num_shared, "\n");

factors_mat=orig_factors_mat[shared,];
counts_mat=orig_counts_mat[shared,];

###############################################################################

normalized_mat=normalize(counts_mat);

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

###############################################################################

plot_text(c(
	paste("Summary Table File: ", InputFileName),
	paste("Factor File: ", FactorFileName),
	paste("Output File Root: ", OutputFileRoot),
	"",
	paste("Distance Metric:", DistanceType),
	"",
	"Summary Table:",
	paste("      Num Samples:", nrow(orig_counts_mat)),
	paste("   Num Categories:", ncol(orig_counts_mat)),
	"",
	"Factor Table:",
	paste("      Num Samples:", nrow(orig_factors_mat)),
	paste("      Num Factors:", ncol(orig_factors_mat)),
	"",
	paste("Shared Samples:", num_shared),
	"",
	"Samples exclusive to Summary Table:",
	capture.output(print(excl_to_st)),
	"",
	"Samples exclusive to Factor Table:",
	capture.output(print(excl_to_fct))
));

###############################################################################

map_val_to_grp=function(fact_mat){
	# This function will convert a factor matrix, into
	# a grouping matrix to reduce the number of continous values

	num_factors=ncol(factors_mat);
	num_values=nrow(factors_mat);
	fact_names=colnames(factors_mat);
	map_mat=as.data.frame(fact_mat);

	for(fidx in 1:num_factors){
		fact_name=fact_names[fidx];
		cat("\nMapping on: ", fact_name, "\n");

		fact_val=factors_mat[,fidx];
		print(fact_val);
		non_na_ix=!is.na(fact_val);
		fact_val=fact_val[non_na_ix];
		num_fact_val=length(fact_val);

		if(is.factor(fact_val)){
			cat(fact_name, " is a factor.\n", sep="");
			fact_lev=levels(fact_val);
			print(fact_lev);
		}else{
			unique_val=unique(fact_val);
			num_unique=length(unique_val);

			if(num_unique<=2){
				cat(fact_name, ": few enough unique values, NOT grouping\n", sep="");
				map_mat[,fidx]=as.character(map_mat[,fidx]);
			}else{
				cat(fact_name, ": too many unique values, grouping...\n", sep="");
				hist_res=hist(fact_val,breaks=nclass.Sturges(fact_val), plot=F);
				cat("Values\n");
				print(fact_val);
				num_grps=length(hist_res$breaks);
				cat("Num Groups: ", num_grps, "\n");

				# Attach an order to the grp level to keep sorting correct
				order_pad=sprintf(paste("%0",floor(log10(num_grps))+1, "i", sep=""), 
					1:(num_grps-1));
				print(order_pad);
			
				grp_levels=paste(
					paste("(", order_pad, ") ", sep=""),
					hist_res$breaks[1:(num_grps-1)], 
					"-", hist_res$breaks[2:num_grps], sep="");
				cat("Group:\n");
				print(grp_levels);

				grp_asn=character(num_fact_val);
				lowerbounds=hist_res$breaks[1:(num_grps-1)];
				for(i in 1:num_fact_val){
					grp_asn[i]=grp_levels[max(which(fact_val[i]>=lowerbounds))];
					#cat(fact_val[i], "->", grp_levels[grp_asn[i]],"\n");	
				}
				cat("Assigned Groups:\n");
				print(grp_asn);

				# Convert strings to factors
				grp_as_factor=factor(grp_asn, levels=grp_levels, ordered=F);
				# Initialize an array
				tmp=rep(grp_as_factor[1],num_values);
				# Copy values over
				tmp[non_na_ix]=grp_asn;
				# Replace NAs
				tmp[setdiff(1:num_values, which(non_na_ix))]=NA;
				map_mat[,fidx]=tmp;
			}

			cat("Unique Val:", unique_val, "\n");
		
		}	
	}
	return(map_mat);
}

###############################################################################

if(!is.null(FactorSubset)){
	fact_subset_arr=scan(FactorSubset, character(), comment.char="#");
	cat("Focusing on subset of factors:\n");
	print(fact_subset_arr);
	factors_mat=factors_mat[,fact_subset_arr, drop=F];
}

grp_mat=map_val_to_grp(factors_mat);
print(grp_mat);

sample_names=rownames(grp_mat);
grp_names=colnames(grp_mat);

###############################################################################

cat("Computing distance matrix...\n");
distmat=compute_dist(normalized_mat, DistanceType);

for(i in 1:length(distmat)){
        if(distmat[i]==0){
                distmat[i]=1e-323;
        }
}

cat("Calculating Ordination...\n");

# PCA
cor_mat=cor(as.matrix(distmat));
eigen_out=eigen(cor_mat);

PC_contributions=eigen_out$values/sum(eigen_out$values);
pc1=eigen_out$vectors[,1];
pc2=eigen_out$vectors[,2];
names(pc1)=rownames(normalized_mat);
names(pc2)=rownames(normalized_mat);

# Non metric
mds=isoMDS(distmat);
mds1_nonmet=mds$points[,1];
mds2_nonmet=mds$points[,2];

# Metric
mds=cmdscale(distmat);
mds1_pcoa=mds[,1];
mds2_pcoa=mds[,2];

plot_mds=function(x, y, samp_grp_map, title, lab=F, cntrd=F){

	cat("Plotting MDS:\n");
	samp_ids=names(x);	
	print(samp_grp_map);

	xrange=range(x);
	yrange=range(y);

	xspan=diff(xrange);
	yspan=diff(yrange);
	
	xlim=c(xrange[1]-xspan*.15, xrange[2]+xspan*.15);
	ylim=c(yrange[1]-yspan*.1, yrange[2]+yspan*.1);

	grps=unique(samp_grp_map[!is.na(samp_grp_map)]);
	ngrps=length(grps)
	centroids=matrix(NA, nrow=ngrps, ncol=2);

	if(cntrd){
		for(i in 1:ngrps){
			gr_ix=which(samp_grp_map==grps[i]);
			centroids[i,]=c(mean(x[gr_ix]), mean(y[gr_ix]));
		}
	}

	plot(0, type="n", xlim=xlim, ylim=ylim, 
		main=title,
		xlab="", ylab="");

	if(cntrd){
		for(i in 1:ngrps){
			points(centroids[i,1], centroids[i,2], cex=3, pch=21, bg=grps[i], col="black");
		}
	}


	if(lab==F){
		points(x, y, col=samp_grp_map[samp_ids], cex=1.1, pch=1, bg="white");
		#text(x, y, samp_grp_map[samp_ids]);
	}else{
		for(cur_samp in samp_ids){
			color=samp_grp_map[cur_samp];
			if(!is.na(color)){
				text(x[cur_samp], y[cur_samp], cur_samp, col=color, cex=.8);
			}
		}
	}

}

plot_legend=function(levels, counts){
	plot(0,0, type="n", bty="n", 
		xaxt="n", yaxt="n",
		xlab="", ylab="", main="", 
		xlim=c(0,1), ylim=c(0,1));
	
	leg_text=paste(levels, " [n=", counts, "]", sep="");

	legend(0,1, legend=leg_text, fill=1:length(levels), bty="n");
}


#------------------------------------------------------------------------------

simple_colors=c(
        "blue", "red", "green", "orange", "violet", "pink", "deepskyblue", "black");
num_simple_colors=length(simple_colors);
palette(simple_colors);

par(oma=c(0,0,2,0));

layout_mat=matrix(c(
	1,1,2,2,3,3,10,
	4,4,5,5,6,6,10,
	7,7,8,8,9,9,10), nrow=3, byrow=T);
layout(layout_mat);

for(i in 1:ncol(grp_mat)){
		
	values=grp_mat[,i];

	if(is.factor(values)){
		all_levels=levels(values);
	}else{
		all_levels=sort(as.numeric(unique(values)));
	}
	colmap=as.numeric(values);
	if(min(colmap, na.rm=T)==0){
		colmap=colmap+1;
	}
	names(colmap)=rownames(grp_mat);
	
	groups=sort(unique(values[!is.na(values)]));
	grp_name=grp_names[i];
	num_grps=length(groups);

	cat("Plotting: ", grp_name, "\n");
	cat("Available Groups: \n");
	print(groups);

	if(num_grps>num_simple_colors){
		cat("Num groups (", num_grps, ") exceeded max (", num_simple_colors, ")\n\n", sep="");
		next;
	}

	group_counts=table(values);

	cat("Num Available Groups: ", num_grps, "\n");
	print(group_counts);

	if(num_grps==0){
		cat("No informations... skipping...\n");
		next;
	}

	par(mar=c(3,3,4,1));
	plot_mds(pc1, pc2, colmap, title="PCA (Principal Comp Analysis)", lab=F, cntrd=F);
	plot_mds(mds1_nonmet, mds2_nonmet, colmap, title="Non-metric MDS", lab=F, cntrd=F);
	plot_mds(mds1_pcoa, mds2_pcoa, colmap, title="PCoA (Classical MDS)", lab=F, cntrd=F);

	par(mar=c(3,3,1,1));
	plot_mds(pc1, pc2, colmap, title="", lab=F, cntrd=T);
	plot_mds(mds1_nonmet, mds2_nonmet, colmap, title="", lab=F, cntrd=T);
	plot_mds(mds1_pcoa, mds2_pcoa, colmap, title="", lab=F, cntrd=T);

	par(mar=c(3,3,1,1));
	plot_mds(pc1, pc2, colmap, title="", lab=T, cntrd=F);
	plot_mds(mds1_nonmet, mds2_nonmet, colmap, title="", lab=T, cntrd=F);
	plot_mds(mds1_pcoa, mds2_pcoa, colmap, title="", lab=T, cntrd=F);

	mtext(grp_name, side=3, line=0, outer=T, font=2);
	
	par(mar=c(3,0,4,0));
	plot_legend(all_levels, group_counts);

	cat("\n");

}

dev.off();


###############################################################################

cat("Done.\n")
warn=warnings();
if(length(warn)){
	print(warn);
}
q(status=0)
