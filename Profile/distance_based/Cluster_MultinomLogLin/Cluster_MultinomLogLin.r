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
DEF_SPLIT_CHAR=";";
DEF_NUM_CLUS=-1;
DEF_COLUMN=1;
RM_NA_TRIALS=10000*64;

params=c(
	"input_summary_table", "i", 1, "character",
	"input_factor_file", "f", 1, "character",
	"model_string", "m", 2, "character",
	"model_filename", "M", 2, "character",
	"subset_samp_ids_fname", "l", 2, "character",
	"output_filename_root", "o", 2, "character",
	"dist_type", "d", 2, "character",
	"num_clus", "k", 2, "numeric",
	"only_at_k", "K", 2, "numeric",
	"rm_na_trials", "N", 2, "numeric",
	"required", "q", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\n\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	-f <input factor file>\n",
	"	[-M <model variables list file>\n",
	"	[-m <\"model string\">]\n",
	"\n",
	"	[-l <list of sample IDs to focus on>\n",
	"	[-o <output file root name, default is input file base name>]\n",
	"	[-d <euc/wrd/man/bray/horn/bin/gow/tyc/minkp5/minkp3, default =", DEF_DISTTYPE, ">]\n",
	"	[-k <max num of clusters to split into, default = ceiling(log2(num_samples))>\n",
	"	[-r <reference file>]\n",
	"	[-q <required variables>]\n",
	"	[-K <only compute at K cuts>]\n",
	"\n",
	"	[-N <remova NA trials, trials=", RM_NA_TRIALS, "\n",
	"\n",
	"This script will:\n",
	"	1.) Load metadata.\n",
	"	2.) Read in a summary table and compute a distance matrix.\n",
	"	3.) Cluster hierarchically.\n",
	"	4.) Iteratively cut tree until max_c clusters.\n",
	"	5.) For each group of clusters:\n",
	"		a.) Fit multinomial log-linear model\n",
	"		b.) Plot coefficients and p-values\n",
	"\n",
	"For the distance types:\n",
	" minkp5 is the minkowski with p=.5, i.e. sum((x_i-y_i)^1/2)^2\n",
	" minkp3 is the minkowski with p=.3, i.e. sum((x_i-y_i)^1/3)^3\n",
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

ModelString="";
if(length(opt$model_string)){
	ModelString=opt$model_string;
}

ModelFilename="";
if(length(opt$model_filename)){
	ModelFilename=opt$model_filename;
}

Num_Remove_NA_Trials=RM_NA_TRIALS;
if(length(opt$rm_na_trials)){
	Num_Remove_NA_Trials=opt$rm_na_trials;
}

if(length(opt$required)){
        RequiredFile=opt$required;
}else{
        RequiredFile="";
}

if(length(opt$subset_samp_ids_fname)){
	SubsetSampIDFname=opt$subset_samp_ids_fname;
}else{
	SubsetSampIDFname="";
}

if(length(opt$only_at_k)){
	OnlyAtK=opt$only_at_k;
}else{
	OnlyAtK=0;
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

load_list=function(list_fname){
	list=read.delim(list_fname, sep="\t", header=F, row.names=NULL, as.is=T, check.names=F, comment.char="#", quote="");
        return(list[,1]);
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


###############################################################################

get_main_effects_from_formula_string=function(formula_string){
	no_space=gsub(" ", "", formula_string);
	no_oper=gsub("[:\\+\\*]", " ", no_space);
	main_effects=strsplit(no_oper, " ")[[1]];
	return(main_effects);
}

###############################################################################
draw_mean_ses=function(mean_matrix, stderr_matrix, title="Mean/SEs Plot", grp_col=NULL, include_zero_ref=F){

	cat("Plotting Mean and Std Errors:\n");
	#print(mean_matrix);
	#print(stderr_matrix);

	num_groups=nrow(mean_matrix);
	num_stats=ncol(mean_matrix);

	stat_names=colnames(mean_matrix);

	par(mar=c(15,4,1,2));

	if(is.null(grp_col)){
		grp_col=1:num_groups;
	}

	extra_pad=.2;
	num_dens_points=200;
	dens_plot_scale=.5; # to keep curves within plot space

	plot(0, xlim=c(0, num_stats), ylim=c(-extra_pad,1+extra_pad), xaxt="n", yaxt="n", xlab="", ylab="", type="n",
		bty="n");
	title(main=title, outer=T);
        axis(side=1, at=1:num_stats, labels=stat_names, las=2, lwd=0);
        axis(side=2, at=c(0,1), labels=c("Min", "Max"), las=2, lwd=0);
	abline(h=c(0,1), lwd=.5, col="grey", lty=2);

	pl_range=par()$usr; # x1, x2, y1, y2
	points(c(0, num_stats), c(0-extra_pad,0-extra_pad), type="l", lwd=2); # Horizontal
	#points(c(0,0), c(0,1), type="l", lwd=2); # Verticle

	for(stt in 1:num_stats){

		if(include_zero_ref){
			zero_ref=0;
		}else{
			zero_ref=c();
		}
		stat_min=min(c(mean_matrix[,stt], zero_ref));
		stat_max=max(c(mean_matrix[,stt], zero_ref));
		stat_rng=stat_max-stat_min;

		orig_val=mean_matrix[,stt];
		norm_val=(mean_matrix[,stt]-stat_min)/stat_rng;

		dens_mat=matrix(0, nrow=num_groups, ncol=num_dens_points);
		for(grp in 1:num_groups){
			dens_x=seq(stat_min-extra_pad*stat_rng, stat_max+extra_pad*stat_rng, length.out=num_dens_points);
			dens_mat[grp,]=dnorm(dens_x, mean_matrix[grp,stt], stderr_matrix[grp,stt]);
		}
		max_dens=max(dens_mat);
		dens_mat=dens_mat/max_dens;

		if(include_zero_ref){
			zero_pos=-stat_min/stat_rng;
			# Label values on axis
			text(stt, zero_pos, 0, cex=.65, pos=4, col="grey");
			# Draw ticks for values
			points(stt, zero_pos, col="grey", pch=3, cex=1.5);
		}

		for(grp in 1:num_groups){

			# Translate/Scale curves to new location and scale
			dens_curv=stt-(dens_mat[grp,])*dens_plot_scale;

			# Label values on axis
			text(stt, norm_val[grp], sprintf("%+-0.4f", orig_val[grp]), cex=.65, pos=4);

			# Draw ticks for values
			points(stt, norm_val[grp], col="black", pch=3, cex=1.5);

			# Draw graph separators
			abline(v=stt, col="black", lwd=2);

			# Draw density curves and label group at peak of curve
			# If stderr is 0, because N=1, set it to zero
			if(is.na(stderr_matrix[grp,stt])){
				stderr_matrix[grp,stt]=0;
				line_type=3;
			}else{
				line_type=1;
			}
			
			if(stderr_matrix[grp,stt]==0){
				points(x=c(stt, stt-1*dens_plot_scale), y=c(norm_val[grp], norm_val[grp]), 
					col=grp_col[grp], lwd=2, type="l", lty=line_type);
				text(stt-1*dens_plot_scale, norm_val[grp], 
					grp, cex=.75, adj=c(1.5,.5), col=grp_col[grp], font=2);
			}else{
				points(x=dens_curv, y=seq(0-extra_pad,1+extra_pad,length.out=num_dens_points), 
					type="l", col=grp_col[grp], lwd=1.75, lty=line_type);
				text(min(dens_curv), norm_val[grp],	
					grp, cex=.75, adj=c(1.5,.5), col=grp_col[grp], font=2);
			}

		}
	}
	
	cat("Done.\n");

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


paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, label_zeros=T, counts=F){

	orig_par=par(no.readonly=T);

        num_row=nrow(mat);
        num_col=ncol(mat);

        cat("Num Rows: ", num_row, "\n");
        cat("Num Cols: ", num_col, "\n");

        mat=mat[rev(1:num_row),, drop=F];

        num_colors=50;
        color_arr=rainbow(num_colors, start=0, end=4/6);
        if(high_is_hot){
                color_arr=rev(color_arr);
        }

        remap=function(in_val, in_range, out_range){
                in_prop=(in_val-in_range[1])/(in_range[2]-in_range[1])
                out_val=in_prop*(out_range[2]-out_range[1])+out_range[1];
                return(out_val);
        }

        if(is.na(plot_min)){
                plot_min=min(mat);
        }
        if(is.na(plot_max)){
                plot_max=max(mat);
        }
        cat("Plot min/max: ", plot_min, "/", plot_max, "\n");

        plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), xaxt="n", yaxt="n", bty="n", xlab="", ylab="", main=title);

        # x-axis
        axis(side=1, at=seq(.5, num_col-.5, 1), labels=colnames(mat), las=2);
        axis(side=2, at=seq(.5, num_row-.5, 1), labels=rownames(mat), las=2);

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

			if(mat[y,x]!=0 || label_zeros){
				if(counts){
					text_lab=sprintf("%i", mat[y,x]);
				}else{
					text_lab=sprintf("%0.4f", mat[y,x]);
				}
				text(x-.5, y-.5, text_lab, srt=45, cex=1, font=2);
			}
        	        }
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

if(OnlyAtK>0){
	ext_k=paste(".k", OnlyAtK, sep="");
}else{
	ext_k="";
}

output_fname_root = paste(OutputFileRoot, ".", dist_type, ext_k, sep="");
cat("\n");
cat("Input Summary Table Name: ", InputFileName, "\n", sep="");
cat("Input Factor File: ", InputFactorFile, "\n", sep="");
cat("Output Filename Root: ", output_fname_root, "\n", sep="");
cat("Distance Type: ", dist_type, "\n", sep="");
cat("Max Num clusters: ", max_clusters, "\n", sep="");
cat("Required Variables File: ", RequiredFile, "\n", sep="");
cat("\n");

if(OnlyAtK>0){
	cat("Only computing at k=", OnlyAtK, "\n");
}

cat("Loading summary table...\n");
counts_mat=load_summary_table(InputFileName);

if(SubsetSampIDFname!=""){
        sample_keep_list=load_list(SubsetSampIDFname);
        st_sample_ids=rownames(counts_mat);
        keep_list=intersect(st_sample_ids, sample_keep_list);
        counts_mat=counts_mat[keep_list,,drop=F];
}

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

required_arr=NULL;
if(RequiredFile != ""){
	required_arr=load_list(RequiredFile);
	cat("Required Variables:\n");
	print(required_arr);
	cat("\n");
}

###############################################################################

if(ModelString==""){
	if(ModelFilename==""){
		ModelString=paste(factor_names, collapse="+");
		main_effects=factor_names;
	}else{
		filelist_var=scan(ModelFilename, "character", comment.char="#");
		filelist_var=gsub("^\\s*","", filelist_var);
		filelist_var=gsub("\\s*$","", filelist_var);
		ModelString=paste(filelist_var, collapse="+");
		main_effects=intersect(factor_names, filelist_var);
	}
}else{
	main_effects=get_main_effects_from_formula_string(ModelString);
}

cat("\n");
cat("Formula String: ", ModelString, "\n");
num_main_effects=length(main_effects);
cat("Number of Main Effects: ", num_main_effects, "\n");
cat("Main Effects:\n");
print(main_effects);
cat("\n");

main_factors=all_factors[,main_effects, drop=F];

###############################################################################

# Reconcile samples between groupings and summary table
shared_samples=sort(intersect(sample_names, factor_samples));
num_shared_samples=length(shared_samples);
cat("Num Shared Samples between Groupings/Summary Table: ", num_shared_samples, "\n");

norm_mat=norm_mat[shared_samples,, drop=F];

factors=main_factors[shared_samples,, drop=F];
num_factors=ncol(factors);

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

na_info=c();
if(any(is.na(main_factors))){
	cat("NAs's found in factors...\n");
	
	script_path=paste(head(strsplit(script_name, "/")[[1]], -1), collapse="/");
	source(paste(script_path, "/../../../Metadata/RemoveNAs/Remove_NAs.r", sep=""));

	before_num_factors=ncol(factors);
	before_num_samples=nrow(factors);

	orig_factor_names=colnames(factors);
	noNA_result=remove_sample_or_factors_wNA_parallel(
			factors, required=required_arr, num_trials=Num_Remove_NA_Trials, 
			num_cores=64, outfile=OutputFileRoot);
	factors=noNA_result$factors;

	after_num_factors=ncol(factors);
	after_num_samples=nrow(factors);

	#print(rownames(factors));
	#print(rownames(norm_mat));
	shared_samples=sort(rownames(factors));

	num_factors=ncol(factors);
	factor_names=colnames(factors);
	factors=factors[shared_samples,, drop=F];
	norm_mat=norm_mat[shared_samples,, drop=F];
	num_shared_samples=nrow(norm_mat);

	removed_factors=setdiff(orig_factor_names, colnames(factors));

	ModelString=rem_missing_var_from_modelstring(ModelString, factor_names);
	main_effects=factor_names;

	na_info=c(
		"NAs were found and removed:",
		"Original:",
		paste("  Num Samples: ", before_num_samples),	
		paste("  Num Factors: ", before_num_factors),	
		"After NA Removal:",
		paste("  Num Samples: ", after_num_samples),	
		paste("  Num Factors: ", after_num_factors),
		"",
		"Factors removed:",
		capture.output(print(removed_factors))
	);

}else{
	na_info=c("No relevant NAs found in metadata.");
}

# Open output PDF file
paper_width=max(28, 2+num_factors/1.75);
pdf(paste(output_fname_root, ".cl_mll.pdf", sep=""), height=8.5, width=paper_width);

# Output run info
plot_text(c(
	paste("Input Summary Table: ", InputFileName, sep=""),
	paste("Input Factor File: ", InputFactorFile, sep=""),
	"",
	paste("Num Shared Samples: ", num_shared_samples, sep=""),
	"",
	paste("Distance Type: ", dist_type, sep=""),
	"",
	na_info
));

###############################################################################

# Output list of used samples
cat("Writing used/shared samples list to file:\n");
fh=file(paste(output_fname_root, ".cl_mll.used_samp.lst", sep=""), "w");
for(samp_id in shared_samples){
	cat(file=fh, samp_id, "\n", sep="");
}
close(fh);
cat("done.\n");

###############################################################################
# Precompute distances and perform clustering

# Compute full distances
cat("Computing distances...\n");
full_dist_mat=compute_dist(norm_mat, dist_type);
for(i in 1:length(full_dist_mat)){
        if(full_dist_mat[i]==0){
                full_dist_mat[i]=1e-323;
        }
}

# Generate hierarchical clustering
hcl=hclust(full_dist_mat, method="ward.D2");

# Find height where cuts are made
if(max_clusters==-1){
	max_clusters=ceiling(log(num_shared_samples, 2));
}
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

# Compute ISO and classical MDS
nonparm_mds_res=matrix(0,nrow=num_shared_samples,ncol=2);
classic_mds_res=matrix(0,nrow=num_shared_samples,ncol=2);

# -- ISO MDS
imds_res=isoMDS(full_dist_mat);
#print(imds_res$points);
nonparm_mds_res[,1]=imds_res$points[,1]
nonparm_mds_res[,2]=imds_res$points[,2]

# -- Classical MDS
classic_mds_res=cmdscale(full_dist_mat);

###############################################################################
# Define layout for multiple plot pages

mds_layout=matrix(c(
	1,1,1,1,2,2,2,2,3), byrow=T, nrow=1, ncol=9);

cont_tab_layout=matrix(c(
	1,1,1,1,1,1,2,
	1,1,1,1,1,1,2,
	1,1,1,1,1,1,2,
	3,3,3,3,3,3,4),
	byrow=T, nrow=4, ncol=7);

indiv_pval_layout=matrix(c(
	1,1,1,1,1,1,2,
	1,1,1,1,1,1,2,
	1,1,1,1,1,1,2,
	1,1,1,1,1,1,2,
	3,3,3,3,3,3,4),
	byrow=T, nrow=5, ncol=7);

# Define static variables used by dendro call back functions
denfun.label_scale=min(2,55/num_shared_samples);


mean_summary=function(x){
	if(is.character(x)){
		if(length(levels(x))==2){
			return(mean(as.numeric(x)));
		}
	}else{
		return(mean(x));
	}
}

se_summary=function(x){
	if(is.factor(x)){
		if(length(levels(x))==2){
			return(sd(as.numeric(x))/sqrt(n));
		}
	}else{
		n=length(x);
		return(sd(x)/sqrt(n));
	}
}

# Begin per cluster count analyses
aics=numeric(max_clusters-1);
i=1;

min_pval_matrix=numeric();
geom_pval_matrix=numeric();
cluster_labels=numeric();


# Store all coefficients and pvalues for building phenotype tree
pvalue_mat_list=list();
coeff_mat_list=list();

geomean=function(x){
	return(exp(mean(log(x))));
}

if(OnlyAtK==0){
	cl_cuts=2:max_clusters;
}else{
	cl_cuts=OnlyAtK;
}

for(num_cl in cl_cuts){

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
	par(oma=c(4,0,5,0));
	par(mar=c(10.1,4.1,4.1,2.1));
	par(mfrow=c(1,1));
	sample_to_color_map=as.list(memberships);

	tweaked_dendro=dendrapply(orig_dendr, color_denfun_bySample);
	tweaked_dendro=dendrapply(tweaked_dendro, text_scale_denfun);

	plot(tweaked_dendro, horiz=F);
	for(cl_ix in 1:num_cl){
		lab_size=3/ceiling(log10(cl_ix+1));
		axis(side=1, outer=T, at=grp_mids[cl_ix], labels=cl_ix, cex.axis=lab_size, col.ticks=cl_ix, 
			lend=1, lwd=10, padj=1, line=-3);
	}

	abline(h=cut_midpoints[num_cl], col="red", lty=2);
	ranges=par()$usr;
	legend(ranges[1]+(ranges[2]-ranges[1])/100, ranges[4], 
		fill=c("white", 1:num_cl), 
		legend=c("Cluster IDs:", as.character(1:num_cl)), 
		border=c("white", rep("black", num_cl)),
		text.font=c(2, rep(1, num_cl)),
		box.col="white", bg="white");
	mtext(paste("Distance Type: ", dist_type), side=3, line=0, outer=T);
	mtext(paste("Num Clusters: ", num_cl), side=3, line=1, outer=T);

	# Generate MDS plots
	par(oma=c(0,0,4,0));
	par(mar=c(5.1,4.1,4.1,2.1));
	layout(mds_layout);

        label_centroids=function(points, memberships){
                members=sort(unique(memberships));
                for(i in members){
                        cur_grp=(i==memberships);
                        midx=mean(points[cur_grp,1]);
                        midy=mean(points[cur_grp,2]);
                        text(midx, midy, label=i, col="black", cex=3/ceiling(log10(i+1)));
                }

        }

	plot(nonparm_mds_res, col=memberships, xlab="Dim 1", ylab="Dim 2", main="non-metric MDS");
	label_centroids(nonparm_mds_res, memberships);
	plot(classic_mds_res, col=memberships, xlab="Dim 1", ylab="Dim 2", main="classical MDS");
	label_centroids(classic_mds_res, memberships);
	
	# MDS Legend
	par(mar=c(0,0,0,0));
	plot(0, type="n", xlab="", ylab="", main="", bty="n", xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1));
	legend(0,1, fill=1:num_cl, legend=c(as.character(1:num_cl)), bty="n", cex=2);

	mtext(paste("Distance Type: ", dist_type, sep=""), side=3, outer=T, line=.5);
	mtext(paste("Num Clusters: ", num_cl, sep=""), side=3, outer=T, line=2);

	#----------------------------------------------------------------------
	# Prepare groups and metadata for multinom function
	
	memberships=memberships[shared_samples];
	memberships_asfactors=as.factor(memberships);

	#print(memberships);
	
	coeff_matrix=numeric();
	coef_se_matrix=numeric();
	pval_matrix=numeric();
	mean_matrix=numeric();
	se_matrix=numeric();

	#pvalue_mat_list[[num_cl]]=matrix(NA, ncol=num_cl, nrow=num_coeff);
	#coeff_mat_list[[num_cl]]=matrix(NA, ncol=num_cl, nrow=num_coeff);

	pvalue_tmp=c();
	coeff_tmp=c();


	aic=NULL;
	for(cl_ix in 1:num_cl){

		cat("Computing Logistic Regression for: ", cl_ix, " of ", num_cl, "\n");
	
		in_group_id = names(memberships[(memberships==cl_ix)]);
		out_group_id = names(memberships[(memberships!=cl_ix)]);

		response=rep(1, num_shared_samples);
		names(response)=shared_samples;
		response[out_group_id]=0;
	
		#print(response);
		mll_formula=paste("response~", ModelString, sep="");
		mm=model.matrix(as.formula(paste("~", ModelString, sep="")), data=factors);
		in_group_id=sort(intersect(in_group_id, rownames(mm)));

		#print(response);
		#print(factors);

		fit=glm(mll_formula, family=binomial, data=factors, control=list(trace=T, maxit=100));

		sum_fit=summary(fit);
		print(sum_fit);

		# Check p-values for degeneracy
		pvals=sum_fit$coefficients[,"Pr(>|z|)"];

		if(all(pvals==0)){
			num_coeff=length(pvals);
			cat("WARNING: glm.fit: fitted probabilities numerically 0 or 1 occurred\n");
			pval_matrix=rbind(pval_matrix, rep(1, num_coeff));
			coeff_matrix=rbind(coeff_matrix, rep(0, num_coeff));
			coef_se_matrix=rbind(coef_se_matrix, rep(0, num_coeff));
		}else{
			pval_matrix=rbind(pval_matrix, pvals);
			coeff_matrix=rbind(coeff_matrix, sum_fit$coefficients[,"Estimate"]);
			coef_se_matrix=rbind(coef_se_matrix, sum_fit$coefficients[,"Std. Error"]);
		}

		# Compute Group Means
		means=apply(mm[in_group_id,,drop=F], 2, mean_summary);
		mean_matrix=rbind(mean_matrix, means);

		# Compute Group SE
		ses=apply(mm[in_group_id,,drop=F], 2, se_summary);
		se_matrix=rbind(se_matrix, ses);
	
		#if(cl_ix!=1){
		#	aic=aic+fit$aic;	
		#}
		if(is.null(aic)){
			aic=fit$aic;
		}else{
			aic=min(c(aic, fit$aic));
		}

		if(length(pvalue_tmp)==0){
			num_coeff=length(sum_fit$coefficients[,"Pr(>|z|)"]);
			pvalue_tmp=matrix(0, nrow=num_coeff, ncol=num_cl);
			coeff_tmp=matrix(0, nrow=num_coeff, ncol=num_cl);
			
			var_names=rownames(sum_fit$coefficients);
			rownames(pvalue_tmp)=var_names;
			rownames(coeff_tmp)=var_names;
			colnames(pvalue_tmp)=1:num_cl;
			colnames(coeff_tmp)=1:num_cl;

			pvalue_tmp[var_names,cl_ix]=sum_fit$coefficients[var_names,"Pr(>|z|)"];
			coeff_tmp[var_names,cl_ix]=sum_fit$coefficients[var_names,"Estimate"];
		}else{
			var_names=rownames(sum_fit$coefficients);
			pvalue_tmp[var_names, cl_ix]=sum_fit$coefficients[var_names,"Pr(>|z|)"];
			coeff_tmp[var_names, cl_ix] =sum_fit$coefficients[var_names,"Estimate"];
		}

	}

	print(pvalue_tmp);
	print(coeff_tmp);
	pvalue_mat_list[[num_cl]]=pvalue_tmp;
	coeff_mat_list[[num_cl]]=coeff_tmp;

	aics[i]=aic;

	# Remove Intercept
	coeff_matrix=coeff_matrix[, -1, drop=F];
	coef_se_matrix=coef_se_matrix[, -1, drop=F];
	pval_matrix=pval_matrix[, -1, drop=F];
	mean_matrix=mean_matrix[, -1, drop=F];
	se_matrix=se_matrix[, -1, drop=F];
	
	rownames(coeff_matrix)=1:num_cl;
	rownames(coef_se_matrix)=1:num_cl;
	rownames(pval_matrix)=1:num_cl;
	rownames(mean_matrix)=1:num_cl;
	rownames(se_matrix)=1:num_cl;
	
	if(1){
		cat("Coefficients:\n");
		print(coeff_matrix);
		cat("\n");
		cat("P-values:\n");
		print(pval_matrix);
		cat("\n");
		cat("Means:\n");
		print(mean_matrix);
		cat("\n");
		cat("SEs:\n");
		print(se_matrix);
		cat("\n");
	}

	signif_coeff_matrix=coeff_matrix * (pval_matrix <= .05);

	par(mfrow=c(1,1));
	par(mar=c(10,10,1,1));
	paint_matrix(mean_matrix, title=paste(num_cl, " Clusters: Means", sep=""));
	draw_mean_ses(mean_matrix, se_matrix, 
		title=paste(num_cl, " Clusters: Means and Standard Errors", sep=""), grp_col=palette_col);
	#paint_matrix(coeff_matrix, title=paste(num_cl, " Clusters, Coefficients", sep=""));
	paint_matrix(coeff_matrix, title=paste(num_cl, " Clusters: Logistic Regression Coefficients", sep=""));
	draw_mean_ses(coeff_matrix, coef_se_matrix, 
		title=paste(num_cl, " Clusters: Logistic Regression Coefficient Standard Errors", sep=""), 
		grp_col=palette_col, include_zero_ref=T);
	paint_matrix(pval_matrix, title=paste(num_cl, " Clusters: P-values", sep=""), 
		high_is_hot=F, plot_max=1, plot_min=0);
	paint_matrix(signif_coeff_matrix,
		title=paste(num_cl, " Clusters: Significient Coefficients (Uncorrected p-value < 0.05)", 
			sep=""), label_zeros=F);

	# Keep min value for each coefficient
	min_pvalues=apply(pval_matrix, 2, min);
	geomean_pvalues=apply(pval_matrix, 2, geomean);
	min_pval_matrix=rbind(min_pval_matrix, min_pvalues);
	geom_pval_matrix=rbind(geom_pval_matrix, geomean_pvalues);
	cluster_labels=cbind(cluster_labels, num_cl);

	cat("*********************************************************************************************************\n");
	i=i+1;	
	
}


plot_coeff_pvalues=function(pval_matrix, line_col, title){

	cat("Plotting: ", title, "\n", sep="");
	matrix_dim=dim(pval_matrix);
	
	num_coefficients=ncol(pval_matrix);
	clus_pts=as.numeric(rownames(pval_matrix));
	max_clusters=max(clus_pts);
	cat("Max Clusters:", max_clusters, "\n");

	# Plot P-values across clusters
	log_min_pval_matrix=log10(pval_matrix);
	print(log_min_pval_matrix);
	log_references=log10(c(0.1, 0.05, 0.025));
	min_log_pval=min(c(log_min_pval_matrix, log_references));
	max_log_pval=max(c(log_min_pval_matrix, log_references));
	cat("Min/Max P-value: ", min_log_pval, "/", max_log_pval, "\n");

	cat("Plotting Number of Clusters vs Log10(P-values)\n");
	if(min_log_pval==-Inf){
		min_log_pval=min(log_min_pval_matrix[log_min_pval_matrix!=-Inf])-1;
	}
	par(mar=c(4, 4, 2, 2));
	plot(0, type="n", xlim=c(2, max_clusters+3), ylim=c(min_log_pval, max_log_pval),
		xlab="Number of Clusters", ylab="Log10(P-values)",
		xaxt="n", bty="l"
	);

	mtext(title, side=3, outer=T, line=2, cex=2, font=2);
	mtext(paste("Num Coefficients Displayed: ", num_coefficients, sep=""), side=3, outer=T, line=0);

	# Reference lines
	points(c(1, max_clusters), rep(log_references[1],2), col="grey", type="l", lwd=.5, lty=2);
	points(c(1, max_clusters), rep(log_references[2],2), col="grey", type="l", lwd=.5, lty=2);
	points(c(1, max_clusters), rep(log_references[3],2), col="grey", type="l", lwd=.5, lty=2);

	if(matrix_dim[2]>0){
		for(i in 1:num_coefficients){
			cur_pvals=log_min_pval_matrix[,i];
			points(clus_pts, cur_pvals, col=line_col[i], lwd=4, type="l", pch=16);
			
			first_lowest=min(which(min(cur_pvals)==cur_pvals));
			points(first_lowest+1, cur_pvals[first_lowest], col=line_col[i], pch=19, cex=3);
			points(first_lowest+1, cur_pvals[first_lowest], col="white",     pch=19, cex=2);
			points(first_lowest+1, cur_pvals[first_lowest], col=line_col[i], pch=19, cex=1);
		}
		axis(side=1, at=clus_pts, labels=clus_pts, cex.axis=2);

		print(log_min_pval_matrix);
		max_clust_pvals=log_min_pval_matrix[as.character(max_clusters),,drop=F];
		order_ix=order(max_clust_pvals);
		coeff_names=colnames(log_min_pval_matrix);

		# Calculate space for labels
		param=par();
		plot_range=param$usr;
		plot_dim=c(plot_range[2]-plot_range[1], plot_range[4]-plot_range[3]);
		char_size=param$cxy;
		max_lines_per_plot=plot_dim[2]/char_size[2];

		# Labels
		padding=abs(min_log_pval-max_log_pval)*.1
		label_y_pos=seq(min_log_pval+padding, max_log_pval-padding, length.out=num_coefficients);
		text(rep(max_clusters+.5, num_coefficients), label_y_pos, 
			pos=4, coeff_names[order_ix], col=line_col[order_ix],
			cex=(max_lines_per_plot/4)/num_coefficients
		);

		# Draw lines from end of plot to label
		for(i in 1:num_coefficients){
			points(c(max_clusters+.1, max_clusters+.5-.1), 
				c(max_clust_pvals[order_ix[i]], label_y_pos[i]), 
				type="l", col=line_col[order_ix[i]], lty=3, lwd=1);
		}
	}
}

num_coeff=ncol(min_pval_matrix);
coef_colors=1:num_coeff;

rownames(min_pval_matrix)=cluster_labels;
rownames(geom_pval_matrix)=cluster_labels;

remove_low_variability_coeff=function(pvmat, delta=.5){
	log_pval=log10(pvmat);
	min_val=apply(log_pval, 2, min);
	max_val=apply(log_pval, 2, max);
	range=(max_val-min_val);
	keep=range>=delta;
	return(keep);
}

remove_coeff_never_below_thres=function(pvmat, thres=.05){
	min_val=apply(pvmat, 2, min);
	keep=min_val<=thres;
	return(keep);
}


# Min P-values across clusters
plot_coeff_pvalues(min_pval_matrix, coef_colors, "All Coefficient Min P-Values");
keep=remove_low_variability_coeff(min_pval_matrix);
plot_coeff_pvalues(min_pval_matrix[,keep, drop=F], coef_colors[keep], "Most Dynamic Coefficient Min P-Values");
keep=remove_coeff_never_below_thres(min_pval_matrix);
plot_coeff_pvalues(min_pval_matrix[,keep, drop=F], coef_colors[keep], "Coefficient P-Values Sometime Significant");

# Geometric Mean P-values across clusters
plot_coeff_pvalues(geom_pval_matrix, coef_colors, "All Coefficient Geometric Mean P-Values");
keep=remove_low_variability_coeff(geom_pval_matrix);
plot_coeff_pvalues(geom_pval_matrix[,keep, drop=F], coef_colors[keep], "Most Dynamic Geometric Mean P-Values");
keep=remove_coeff_never_below_thres(geom_pval_matrix);
plot_coeff_pvalues(geom_pval_matrix[,keep, drop=F], coef_colors[keep], "Coefficient Geometric Mean P-Values Sometime Significant");

# Extract min pvalues
best_cl_cut=as.data.frame(matrix(NA, nrow=ncol(min_pval_matrix), ncol=3));
colnames(best_cl_cut)=c("cluster_cut", "p-value", "signf_char");
rownames(best_cl_cut)=colnames(min_pval_matrix);
cl_cuts=rownames(min_pval_matrix);

sig_char=function(val){
	if(val == 0){ return("***");}
	if(val <= .001){ return("** ");}
	if(val <= .01){ return("*  ");}
	if(val <= .05){ return(".  ");}
	return(" ");
}

for(i in 1:ncol(min_pval_matrix)){
	min_pv_ix=min(which(min_pval_matrix[,i]==min(min_pval_matrix[,i])));
	best_cl_cut[i, "p-value"]=min_pval_matrix[min_pv_ix, i];
	best_cl_cut[i, "cluster_cut"]=as.numeric(cl_cuts[min_pv_ix]);
	best_cl_cut[i, "signf_char"]=sig_char(min_pval_matrix[min_pv_ix, i]);
}

pval_ord_ix=order(best_cl_cut[,"p-value"]);
best_cl_cut=best_cl_cut[pval_ord_ix,];
print(best_cl_cut);

vnam_ord_ix=order(rownames(best_cl_cut));
clct_ord_ix=order(best_cl_cut[, "cluster_cut"]);


options(width=200);
plot_text(c(
	"Ordered By: Variable Name",
	capture.output(print(best_cl_cut[vnam_ord_ix,], quote=F))
));

plot_text(c(
	"Ordered By: Min P-Value",
	capture.output(print(best_cl_cut, quote=F))
));

plot_text(c(
	"Ordered By: Cluster Cuts",
	capture.output(print(best_cl_cut[clct_ord_ix,], quote=F))
));

###############################################################################

dev.off();

plot_tree_phenotypes=function(hcl, coef_mat_list, pval_mat_list, alpha=.10){

	print(coeff_mat_list);
	print(pvalue_mat_list);

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
	
	#-----------------------------------------------------------------------------

	num_cuts=length(coeff_mat_list);
	num_variables=nrow(coeff_mat_list[[2]]);
	
	cat("Num Cuts:", num_cuts, "\n");
	cat("Num Variables:", num_variables, "\n");

	best_cut_mat=matrix(NA, nrow=num_variables, ncol=4);
	rownames(best_cut_mat)=rownames(coeff_mat_list[[2]]);
	colnames(best_cut_mat)=c("pvalue", "coeff", "cut", "cluster");
	best_cut_mat[,"pvalue"]=1;
	best_cut_mat[,"coeff"]=0;
	
	# Find cluster cut with lowest pvalue for each variable
	for(cut_ix in 2:num_cuts){
		cur_cut_mat=pval_mat_list[[cut_ix]];
		for(var_ix in 1:num_variables){
			for(cl_ix in 1:cut_ix){
				if(cur_cut_mat[var_ix, cl_ix] < best_cut_mat[var_ix, "pvalue"]){
					best_cut_mat[var_ix, "pvalue"]=cur_cut_mat[var_ix, cl_ix];
					best_cut_mat[var_ix, "coeff"]=coef_mat_list[[cut_ix]][var_ix, cl_ix];
					best_cut_mat[var_ix, "cut"]=cut_ix;
					best_cut_mat[var_ix, "cluster"]=cl_ix;
				}
			}	
		}		
	}
	
	cat("\n\nAll Minimums:\n");
	print(best_cut_mat);

	# Remove non-significant
	signf_ix=best_cut_mat[,"pvalue"]<=alpha;
	best_cut_mat=best_cut_mat[signf_ix,,drop=F];	
	cat("\n\nSignificant Minimums:\n");
	print(best_cut_mat);

	if(nrow(best_cut_mat)==0){
		plot(0,0, type="n", xlim=c(0,1), ylim=c(0,1), 
			bty="n", xaxt="n", yaxt="n", main="", xlab="", ylab="");	
		msg=paste("No variables with signficant cuts at alpha<=", alpha, sep="");
		text(0,0, msg, cex=2, font=2);
		cat(msg, "\n");
		return();
	}

	# Sort by cut
	cluster_order=order(best_cut_mat[,"cluster"]);
	best_cut_mat_byClust=best_cut_mat[cluster_order,,drop=F];

	cut_order=order(best_cut_mat_byClust[,"cut"]);
	best_cut_mat=best_cut_mat_byClust[cut_order,,drop=F];

	cat("\n\nSorted by Cut:\n");
	print(best_cut_mat);
	plot_text(capture.output(best_cut_mat));
	num_kept_signf_var=nrow(best_cut_mat);	
	cat("Num Significant Variables:", num_kept_signf_var, "\n");

	# Look for addition clusters with significant associations in same cut
	additional_cuts_list=list();
	num_signf_var=nrow(best_cut_mat);
	signf_var_names=rownames(best_cut_mat);
	add_ix=1;
	for(var_ix in 1:num_signf_var){
		cur_cut=best_cut_mat[var_ix,"cut"];
		cur_var_name=signf_var_names[var_ix];
		all_cuts_pval=pval_mat_list[[cur_cut]][cur_var_name,,drop=F];

		signf_cuts=all_cuts_pval<=alpha;

		add_rec=list();
		add_rec[["cut"]]=cur_cut;
		add_rec[["cluster"]]=which(signf_cuts);
		add_rec[["pvalue"]]=pval_mat_list[[cur_cut]][cur_var_name,signf_cuts];
		add_rec[["coeff"]]=coef_mat_list[[cur_cut]][cur_var_name,signf_cuts];
		additional_cuts_list[[cur_var_name]]=add_rec;
	}
	print(additional_cuts_list);

	#-----------------------------------------------------------------------------

	samp_dendro=as.dendrogram(hcl);
	lf_names=get_clstrd_leaf_names(samp_dendro);
	num_dendro_samples=length(lf_names);
	
	max_cuts=max(best_cut_mat[,"cut"]);
	cat("Max Cuts to Plot:", max_cuts, "\n");

	layout_mat=matrix(c(
		1,
		rep(2, 4),
		3), ncol=1);
	layout(layout_mat);

	left_mar=15;
	right_mar=5;

	membership_matrix=numeric();

	clst_sizes=list();

	for(cut_ix in 2:max_cuts){

		# Cut tree at target number of clusters
		cat("Cutting for ", cut_ix, " clusters...\n", sep="");
		memberships=cutree(hcl, k=cut_ix);
		grp_mids=get_middle_of_groups(lf_names, memberships);

		# Reorder cluster assignments to match dendrogram left/right
		plot_order=order(grp_mids);
		mem_tmp=numeric(num_dendro_samples);
		for(gr_ix in 1:cut_ix){
			old_id=(memberships==plot_order[gr_ix]);
			mem_tmp[old_id]=gr_ix;
		}
		names(mem_tmp)=names(memberships);
		memberships=mem_tmp;
		grp_mids=grp_mids[plot_order];

		membership_matrix=cbind(membership_matrix, memberships);

		if(cut_ix==max_cuts){

			# Plot Dendrogram
			par(mar=c(8,left_mar,2,right_mar));
			sample_to_color_map=as.list(memberships);

			tweaked_dendroT=dendrapply(samp_dendro, color_denfun_bySample);
			tweaked_dendroT=dendrapply(tweaked_dendroT, text_scale_denfun);
			plot(tweaked_dendroT, horiz=F);
			for(cl_ix in 1:cut_ix){
				lab_size=3/ceiling(log10(cl_ix+1));
				axis(side=1, at=grp_mids[cl_ix], labels=cl_ix, 
					cex.axis=lab_size, col.ticks=cl_ix, 
					lend=1, lwd=10, padj=1, line=4);
			}

			#points((1:num_dendro_samples), rep(2, num_dendro_samples));
		}

		clst_sizes[[cut_ix]]=table(memberships);

	}

	# Compute cluster evolution
	memb_map=unique(apply(membership_matrix, 1, function(x){paste(x, collapse=",")}));
	membership_matrix=c();
	for(i in 1:length(memb_map)){
		membership_matrix=rbind(membership_matrix, as.numeric(strsplit(memb_map[i], ",")[[1]]));
	}
	colnames(membership_matrix)=paste("k=", 2:(ncol(membership_matrix)+1), sep="");
	sort_ix=order(membership_matrix[, ncol(membership_matrix)]);
	membership_matrix=membership_matrix[sort_ix,];
	print(membership_matrix);

	#-----------------------------------------------------------------------------

	par(mar=c(1,left_mar,1,right_mar));
	plot(0,0, type="n", 
		xlab="", ylab="", bty="n", xaxt="n", yaxt="n",
		xlim=c(0,num_dendro_samples), ylim=c(0, num_kept_signf_var));

	var_y_pos=rev(1:num_kept_signf_var);
	axis(side=2, at=var_y_pos, labels=rownames(best_cut_mat),
		las=2);

	axis(side=4, at=var_y_pos, labels=paste("k=",best_cut_mat[,"cut"]), las=2);


	get_signf_char= function(x){
		if(x<=.001){return("***");}
		if(x<=.01){return("**");}
		if(x<=.05){return("*");}
		return("");
	};


	abline(h=var_y_pos, col="grey80");

	for(var_ix in 1:num_kept_signf_var){

		cur_var_rec=additional_cuts_list[[var_ix]];
		cut=cur_var_rec[["cut"]];
		num_cl=length(cur_var_rec$cluster);

		cur_clst_size=clst_sizes[[cut]];
		stst_pos=c(0, cumsum(cur_clst_size));
		points(stst_pos, rep(var_y_pos[var_ix], cut+1), cex=2);

		for(cl in 1:num_cl){
	
			#pval=best_cut_mat[var_ix, "pvalue"];	
			#coef=best_cut_mat[var_ix, "coeff"];
			#cut=best_cut_mat[var_ix, "cut"];	
			#clus=best_cut_mat[var_ix, "cluster"];	

			pval=cur_var_rec[["pvalue"]][cl];
			coef=cur_var_rec[["coeff"]][cl];
			clus=cur_var_rec[["cluster"]][cl];

			cat("\n");
			cat("P-val: ", pval, "\n");
			cat("Coeff: ", coef, "\n");
			cat("Cuts : ", cut, "\n");
			cat("Clust: ", clus, "\n");

		
			bar_pos=c(stst_pos[clus], stst_pos[clus+1]);

			if(coef>0){
				bar_col="red";
			}else{
				bar_col="blue";
			}

			# Draw bar
			points(bar_pos, rep(var_y_pos[var_ix],2), type="l", col=bar_col, lend="butt", lwd=5);

			# Draw bar line ends
			points(rep(bar_pos[1],2), var_y_pos[var_ix]+c(.25, -.25), 
				col="black", type="l", lwd=2); 
			points(rep(bar_pos[2],2), var_y_pos[var_ix]+c(.25, -.25), 
				col="black", type="l", lwd=2); 

			# Draw asterisks over center of bar
			text(mean(bar_pos), var_y_pos[var_ix], get_signf_char(pval), font=2, adj=c(.5,.3), cex=2);
		
		}

	}

	par(mar=c(0,0,0,0));

	orig_fam=par()$family;
	par(family="mono");

	plot(0,0, type="n", xlim=c(0,1), ylim=c(0,1), bty="n", xaxt="n", yaxt="n", main="", xlab="", ylab="");
	legend(0, 1, title="Associations:", bty="n", cex=1.7,
		legend=c("Positive", "Negative"),
		fill=c("red", "blue"));

	legend(.25, 1, title="Significance:", bty="n", cex=1.7,
		legend=c(
			"***, p < 0.001", 
			" **, p < 0.01", 
			"  *, p < 0.05", 
			"   , p < 0.1"));	
	par(family=orig_fam);

	#-----------------------------------------------------------------------------
	# Output table:
	# variable name, pval, coeff, cut, clusters_inv, propo
	
	print(membership_matrix);
	summary_matrix=character();
	varnames=names(additional_cuts_list);
	for(var_ix in 1:num_kept_signf_var){
		
		cur_var_rec=additional_cuts_list[[var_ix]];
		cut=cur_var_rec[["cut"]];
		num_cl=length(cur_var_rec$cluster);
		varname=names(cur_var_rec);

		cur_clst_size=clst_sizes[[cut]];

		for(cl in 1:num_cl){

			clus=cur_var_rec[["cluster"]][cl];

			# Look up subclusters
			cutrow=membership_matrix[, paste("k=",cut,sep="")];
			subcl=membership_matrix[cutrow==clus, paste("k=", max_cuts, sep="")];
			subcl_str=paste(subcl, collapse=",");

			summary_matrix=rbind(summary_matrix,
				c(
					varnames[var_ix],
					sprintf("%8.3f", cur_var_rec[["coeff"]][cl]),
					sprintf("%5.3f", cur_var_rec[["pvalue"]][cl]),
					get_signf_char(cur_var_rec[["pvalue"]][cl]),
					cut,
					subcl_str,
					cur_clst_size[clus],
					round(cur_clst_size[clus]/num_dendro_samples, 3)
			));

		
		}
	}

	colnames(summary_matrix)=c(
		"Variable",
		"Coeff",
		"P-Val",
		"Signf",
		"k",
		"SubclstInc",
		"NumSamp",
		"PropSamp");

	rownames(summary_matrix)=1:nrow(summary_matrix);

	# Truncate long variable names
	summary_matrix[,"Variable"]=substr(summary_matrix[,"Variable"],1,35);

	table=capture.output(print(summary_matrix, quote=F));
	par(mfrow=c(1,1));
	plot_text(table);
}

pdf(paste(output_fname_root, ".cl_mll.phenotree.pdf", sep=""), height=11, width=8.5);

plot_tree_phenotypes(hcl=hcl, coef_mat_list=coeff_mat_list, pval_mat_list=pvalue_mat_list);

dev.off();

###############################################################################

cat("\nDone.\n")
warn=warnings();
if(length(warn)){
        print(warn);
}
q(status=0);
