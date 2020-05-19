#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');

source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');
source('~/git/AnalysisTools/Longitudinal/Longitudinal.r');

options(useFancyQuotes=F);

params=c(
	"summary_file", "s", 1, "character",

	"num_top_pred", "v", 2, "numeric",
	"contains_remaining", "R", 2, "logical",
	"shorten_category_names", "x", 2, "character",

	"factor_file", "f", 1, "character",
	"offset_col", "t", 1, "character",
	"subject_id_col", "i", 1, "character",

	"output_root", "o", 1, "character",

	"alpha", "a", 2, "numeric",

	"model_file", "m", 2, "character",
	"group_col", "g", 2, "character",

	"dont_reset_offsets", "n", 2, "logical",
	"begin_offset", "b", 2, "numeric",
	"end_offset", "e", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

NUM_TOP_PRED_CAT=20;
ALPHA=0.1;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table>\n",
	"\n",
	"	[-v <number of top predictor (as ALR) categories to include, default=", NUM_TOP_PRED_CAT, ">]\n",
	"	[-R (pay attention to 'remaining' category)]\n",
	"	[-x <shorten category names, with separator in double quotes (default=\"\")>]\n",
	"\n",
	"	-f <factor file name>\n",
	"	-t <offset column name>\n",
	"	-i <subject identifier column name>\n",

	"	-o <output root>\n",
	"\n",
	"	[-a <p-value cutoff for non parametric longitudinal statistics, default=", ALPHA, ">]\n",
	"\n",
	"	[-m <model variable list>]\n",
	"	[-g <group/cohort variable column name>]\n",
	"\n",
	"	[-n (do not reset earliest offsets to 0 to line up time points, default=reset offsets)]\n",
	"	[-b <begin offset, default=-Inf>]\n",
	"	[-e <end offset, default=Inf>]\n",
	"\n",	
	"\n", sep="");

if(
	!length(opt$summary_file) || 
	!length(opt$factor_file) || 
	!length(opt$offset_col) || 
	!length(opt$subject_id_col) || 
	!length(opt$output_root)
	){
	cat(usage);
	q(status=-1);
}

# Required
SummaryFile=opt$summary_file;
FactorFile=opt$factor_file;
OutputRoot=opt$output_root;
SubjectIDCol=opt$subject_id_col;
TimeOffsetCol=opt$offset_col;

# Optional, i.e. with defaults
NumALRPredictors=NUM_TOP_PRED_CAT;
UseRemaining=F;
ShortenCategoryNames="";
Alpha=ALPHA;
ModelFile="";
ResetOffsets=T;
BeginOffset=-Inf;
EndOffset=Inf;
GroupCol="";

if(length(opt$num_top_pred)){
	NumALRPredictors=opt$num_top_pred;
}

if(length(opt$contains_remaining)){
	UseRemaining=T;
}

if(length(opt$shorten_category_names)){
	ShortenCategoryNames=opt$shorten_category_names;
}

if(length(opt$alpha)){
	Alpha=opt$alpha;
}

if(length(opt$model_file)){
	ModelFile=opt$model_file;
}

if(length(opt$dont_reset_offsets)){
	ResetOffsets=F;
}

if(length(opt$begin_offset)){
	BeginOffset=opt$begin_offset;
}

if(length(opt$end_offset)){
	EndOffset=opt$end_offset;
}

if(length(opt$group_col)){
	GroupCol=opt$group_col;
}

###############################################################################

input_param=capture.output({
	cat("\n");
	cat("Summary File: ", SummaryFile, "\n", sep="");
	cat("  Num Top ALR Predictors: ", NumALRPredictors, "\n", sep="");
	cat("  Contains 'Remaining': ", UseRemaining, "\n", sep="");
	cat("  Shorten Categories: ", ShortenCategoryNames, "\n", sep="");
	cat("\n");
	cat("Factor File: ", FactorFile, "\n", sep="");
	cat("  Subject ID Column Name: ", SubjectIDCol, "\n", sep="");
	cat("  Time Offsets Column Name: ", TimeOffsetCol, "\n", sep="");
	cat("  Group Column Name: ", GroupCol, "\n", sep="");
	cat("  Reset Offsets? ", ResetOffsets, "\n", sep="");
	cat("\n");
	cat("Output File Root: ", OutputRoot, "\n", sep="");
	cat("\n");
});

cat(paste(input_param, collapse="\n"));


if(ShortenCategoryNames==TRUE){
	cat("Error:  You need to specify a delimitor to split the category names.\n");
	cat("        i.e., this -x option is not a flag, it requires a parameter.\n");
	quit(status=-1);
}

options(width=100);
cat("Text Line Width: ", options()$width, "\n", sep="");

##############################################################################
##############################################################################

# General functions

load_summary_file=function(fname){
	inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, 
		comment.char="", quote="", row.names=1))
	counts_mat=inmat[,2:(ncol(inmat))];

	# Clean category names a little
	cat_names=colnames(counts_mat);
	cat_names=gsub("-", "_", cat_names);
	colnames(counts_mat)=cat_names;
	
	return(counts_mat);
}

load_factors=function(fname){
        factors=data.frame(read.table(fname,  sep="\t", header=TRUE, row.names=1,
                check.names=FALSE, comment.char=""));
        factor_names=colnames(factors);

        ignore_idx=grep("^IGNORE\\.", factor_names);

        if(length(ignore_idx)!=0){
                return(factors[-ignore_idx]);
        }else{
                return(factors);
        }
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

load_list=function(filename){
	val=scan(filename, what=character(), comment.char="#");
	return(val);
}

extract_top_categories=function(ordered_normalized, top){

	num_samples=nrow(ordered_normalized);
	num_categories=ncol(ordered_normalized);

	cat("Samples: ", num_samples, "\n");
	cat("Categories: ", num_categories, "\n");
	
	num_saved=min(c(num_categories, top+1));

	cat("Top Requested to Extract: ", top, "\n");
	cat("Columns to Extract: ", num_saved, "\n");

	top_cat=matrix(0, nrow=num_samples, ncol=num_saved);
	top=num_saved-1;

	# Extract top categories requested
	top_cat[,1:top]=ordered_normalized[,1:top];

	# Included remaineder as sum of remaining categories
	top_cat[,(top+1)]=apply(
		ordered_normalized[,(top+1):num_categories, drop=F],
		1, sum);

	rownames(top_cat)=rownames(ordered_normalized);
	colnames(top_cat)=c(colnames(ordered_normalized)[1:top], "Remaining");

	return(top_cat);
			
}

additive_log_rato=function(ordered_matrix){
# Assumes last column will be the denominator

	num_cat=ncol(ordered_matrix);
	num_samp=nrow(ordered_matrix);

	denominator=ordered_matrix[,num_cat];
	alr_mat=matrix(0, nrow=num_samp, ncol=(num_cat-1));
	
	for(i in 1:num_samp){
		alr_mat[i,]=log(ordered_matrix[i,1:(num_cat-1)]/denominator[i]);
		#print(alr_mat[i,]);
	}

	rownames(alr_mat)=rownames(ordered_matrix)
	colnames(alr_mat)=head(colnames(ordered_matrix), num_cat-1);

	alr_struct=list();
	alr_struct[["transformed"]]=alr_mat;
	alr_struct[["denominator"]]=denominator;

	return(alr_struct);
}

title_page=function(string){

	plot(0,0, xlim=c(-1,1), ylim=c(-1,1), type="n",  xaxt="n", yaxt="n",
		xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
		);

	text(0,0, string, font=2, cex=2);

}

plot_text=function(strings, maxlpp=100){

	nlines=length(strings);
	if(nlines>maxlpp){
		plot_text(strings[1:maxlpp]);
		plot_text(strings[(maxlpp+1):nlines]);
	}else{

		par(family="Courier");
		par(mar=rep(0,4));

		num_lines=length(strings);
		cat("Num Plot Text lines:", num_lines, "\n");
		
		top=maxlpp;

		plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
			xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
			);

		text_size=max(.01, min(.8, .8 - .003*(maxlpp-52)));
		#print(text_size);

		for(i in 1:num_lines){
			#cat(strings[i], "\n", sep="");
			strings[i]=gsub("\t", "", strings[i]);
			text(0, top-i, strings[i], pos=4, cex=text_size); 
		}
	}
}

mask_matrix=function(val_mat, mask_mat, mask_thres, mask_val){
        masked_matrix=val_mat;
        masked_matrix[mask_mat>mask_thres]=mask_val;
        return(masked_matrix);
}


paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, deci_pts=4,
        label_zeros=T, counts=F, value.cex=1,
        plot_col_dendr=F,
        plot_row_dendr=F
){

        cat("Working on: ", title, "\n");

        num_row=nrow(mat);
        num_col=ncol(mat);

        if(num_row==0 || num_col==0){
                cat("Nothing to plot.\n");
                return();
        }

        any_nas=any(is.na(mat));

        if(num_row==1 || any_nas){
                plot_row_dendr=F;
        }
        if(num_col==1 || any_nas){
                plot_col_dendr=F;
        }

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

                mat=mat[row_dendr[["names"]], col_dendr[["names"]], drop=F];

        }else if(plot_col_dendr){
                layoutmat=matrix(
                        c(
                        rep(rep(2, heatmap_width), col_dend_height),
                        rep(rep(1, heatmap_width), heatmap_height)
                        ), byrow=T, ncol=heatmap_width);

                col_dendr=get_dendrogram(mat, type="col");
                mat=mat[, col_dendr[["names"]], drop=F];

        }else if(plot_row_dendr){
                layoutmat=matrix(
                        rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height),
                        byrow=T, ncol=row_dend_width+heatmap_width);

                row_dendr=get_dendrogram(mat, type="row");
                mat=mat[row_dendr[["names"]],, drop=F];
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
        axis(side=1, at=seq(.5, num_col-.5, 1), labels=colnames(mat), las=2, line=-1.75, cex.axis=value.cex);
        axis(side=4, at=seq(.5, num_row-.5, 1), labels=rownames(mat), las=2, line=-1.75, cex.axis=value.cex);

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

get_colors=function(num_col, alpha=1){
        colors=hsv(seq(0,1,length.out=num_col+1), c(1,.5), c(1,.75,.5), alpha=alpha);
        color_mat_dim=ceiling(sqrt(num_col));
        color_pad=rep("grey", color_mat_dim^2);
        color_pad[1:num_col]=colors[1:num_col];
        color_mat=matrix(color_pad, nrow=color_mat_dim, ncol=color_mat_dim);
        colors=as.vector(t(color_mat));
        colors=colors[colors!="grey"];
}

#------------------------------------------------------------------------------
# Analysis specific functions

plot_alr_time_indv=function(tar_cat, tar_subj, offsets_rec, alr_categories_val, 
	offset_range, alr_range, alr_med, col="black"){

	indv_offsets=offsets_rec[["OffsetsBySubject"]][[tar_subj]];
	samp_ids=rownames(indv_offsets);

	cat("Working on: ", tar_cat, " for ", tar_subj, "\n");

	x_val=indv_offsets[,"Offsets"];
	y_val=alr_categories_val[samp_ids, tar_cat];

	plot(0, type="n", ylim=alr_range, xlim=offset_range, ylab=tar_subj);

	if(min(offset_range)<0){
		abline(v=0, col="blue", lty=3, lwd=3, lend="butt");
	}

	abline(h=alr_med, col="grey", lty="dotdash");
	points(x_val, y_val, type="l", lwd=5, col=col);
	points(x_val, y_val, type="l", lwd=.5, col="black");
	points(x_val, y_val, type="p", cex=3, col="black");


}

plot_alr_time_grpd=function(tar_cat, subj_arr, grouping, grouping_name, offsets_rec, alr_categories_val, 
	offset_range, alr_range, alr_med, subj_cols){

	plot(0, type="n", ylim=alr_range, xlim=offset_range, ylab=paste(grouping_name, ": ", grouping));
	abline(h=alr_med, col="grey", lty="dotdash");

	if(min(offset_range)<0){
		abline(v=0, col="blue", lty=3, lwd=3, lend="butt");
	}

	for(tar_subj in subj_arr){

		cat("Working on: ", tar_cat, " for ", tar_subj, "\n");

		indv_offsets=offsets_rec[["OffsetsBySubject"]][[tar_subj]];
		samp_ids=rownames(indv_offsets);

		x_val=indv_offsets[,"Offsets"];
		y_val=alr_categories_val[samp_ids, tar_cat];

		points(x_val, y_val, type="l", lwd=5, col=subj_cols[tar_subj]);
		points(x_val, y_val, type="l", lwd=.5, col="black");
		points(x_val, y_val, type="p", cex=1, col="black");
	}

}

compute_and_plot_loess=function(tar_cat, subj_arr, grouping, grouping_name, offsets_rec, alr_categories_val, 
	offset_range, alr_range, alr_med, subj_cols){

	plot(0, type="n", ylim=alr_range, xlim=offset_range, ylab=paste(grouping_name, ": ", grouping));
	abline(h=alr_med, col="grey", lty="dotdash");

	if(min(offset_range)<0){
		abline(v=0, col="blue", lty=3, lwd=3, lend="butt");
	}

	all_x=numeric();
	all_y=numeric();
	for(tar_subj in subj_arr){

		cat("Working on: ", tar_cat, " for ", tar_subj, "\n");

		indv_offsets=offsets_rec[["OffsetsBySubject"]][[tar_subj]];
		samp_ids=rownames(indv_offsets);

		x_val=indv_offsets[,"Offsets"];
		y_val=alr_categories_val[samp_ids, tar_cat];

		all_x=c(all_x, x_val);
		all_y=c(all_y, y_val);

		points(x_val, y_val, type="p", cex=1.5, col=subj_cols[tar_subj]);
	}

	order_x=order(all_x);
	all_x=all_x[order_x];
	all_y=all_y[order_x];

	loess_res=loess(all_y~all_x);

	grp_loess=cbind(loess_res[["x"]], loess_res[["fitted"]]);
	colnames(grp_loess)=c("x", "y");

	points(grp_loess[,"x"], grp_loess[,"y"], type="l", col="blue");

	return(grp_loess);

}

plot_barplot_wsignf_annot=function(title, stat, grps, alpha=0.1, samp_gly=T){
        # Generate a barplot based on stats and groupings
        # Annotat barplot with signficance

        cat("Making Barplot with Significance annotated...\n");
        cat("  Alpha", alpha, "\n");
        group_names=names(grps);
        num_grps=length(group_names);

        # Convert matrix into array, if necessary
        if(!is.null(dim(stat))){
                stat_name=colnames(stat);
                stat=stat[,1];
        }else{
                stat_name="value";
        }

        # Remove NAs
        na_ix=is.na(stat);
        subj=names(stat);
        stat=stat[!na_ix];
        na_subj=names(stat);
        for(grnm in group_names){
                grps[[grnm]]=intersect(grps[[grnm]], na_subj);
                print(stat[grps[[grnm]]]);
        }
        print(grps);

        # Precompute pairwise wilcoxon pvalues
        cat("\n  Precomputing group pairwise p-values...\n");
        pval_mat=matrix(1, nrow=num_grps, ncol=num_grps);
        rownames(pval_mat)=group_names;
        colnames(pval_mat)=group_names;
        signf=numeric();
        for(grp_ix_A in 1:num_grps){
                for(grp_ix_B in 1:num_grps){
                        if(grp_ix_A<grp_ix_B){

                                grpAnm=group_names[grp_ix_A];
                                grpBnm=group_names[grp_ix_B];

                                res=wilcox.test(stat[grps[[grpAnm]]], stat[grps[[grpBnm]]]);

				if(is.na(res$p.value)){
					pval=1;
				}else{
					pval=res$p.value;
				}
                                pval_mat[grpAnm, grpBnm]=pval;
                                if(pval<=alpha){
					amean=mean(stat[grps[[grpAnm]]]);
					bmean=mean(stat[grps[[grpBnm]]]);
                                        signf=rbind(signf, c(
						grpAnm, signif(amean, 3), 
						grpBnm, signif(bmean, 3), 
						signif(res$p.value,3)));
                                }
                        }
                }
        }

        cat("p-value matrix:\n");
        print(pval_mat);

        # Count how many rows have significant pairings
        num_signf=nrow(signf);
        cat("  Num Significant: ", num_signf, "\n");
        signf_by_row=apply(pval_mat, 1, function(x){sum(x<alpha)});
        cat("  Num Significant by Row:\n");
        print(signf_by_row);

        num_signf_rows=sum(signf_by_row>0);
        cat("  Num Rows to plot:", num_signf_rows, "\n");

        #signf_mat=apply(pval_mat, 1:2,
        #       function(x){
        #               if(x<.001){return("***")}
        #               if(x<.01){return("**")}
        #               if(x<.05){return("*")}
        #               else{return("")}
        #       }
        #);

        #print(signf_mat, quote=F);

        # Compute 95% CI around mean
        cat("\n  Precomputing group means and 95% CI...\n");
        num_bs=320;

        grp_means=numeric(num_grps);
        names(grp_means)=group_names;

        ci95=matrix(NA, nrow=num_grps, ncol=2);
        rownames(ci95)=group_names;
        colnames(ci95)=c("LB", "UB");
        samp_size=numeric(num_grps);
        for(grp_ix in 1:num_grps){

                grpnm=group_names[grp_ix];
                grp_means[grpnm]=mean(stat[grps[[grpnm]]]);
                num_samp=length(grps[[grpnm]]);

                if(num_samp>=40){
                        meds=numeric(num_bs);
                        for(i in 1:num_bs){
                                meds[i]=mean(sample(stat[grps[[grpnm]]], replace=T));

                        }
                        ci95[grp_ix,]=quantile(meds, c(.025, .975));
                }else{
                        ci95[grp_ix,]=rep(mean(stat[grps[[grpnm]]]),2);
                }

                samp_size[grp_ix]=num_samp;
        }

        cat("Group Means:\n");
        print(grp_means);
        print(length(grp_means));
        cat("Group Median 95% CI:\n");
        print(ci95);

        # Estimate spacing for annotations
        annot_line_prop=1/5; # proportion of pl
        min_95ci=min(c(ci95[,1], stat), na.rm=T);
        max_95ci=max(c(ci95[,2], stat), na.rm=T);
        minmax_span=max_95ci-min_95ci;
        plotdatamax=max_95ci+minmax_span*0.3;
        plotdatamin=min_95ci-minmax_span*0.3;;
        space_for_annotations=minmax_span*annot_line_prop*(num_signf_rows+2);
        horiz_spacing=annot_line_prop*plotdatamax;
        # Start plot
        par(mar=c(8,8,4,3));
        cat("  Plot Limits: (", plotdatamin, ", ", plotdatamax, ")\n");
        plot(0, type="n",
                ylim=c(plotdatamin, plotdatamax+space_for_annotations),
                xlim=c(0, num_grps+1),
                yaxt="n", xaxt="n", xlab="", ylab="", bty="n");
        for(grp_ix in 1:num_grps){
                points(c(grp_ix-.25, grp_ix+.25), rep(grp_means[grp_ix],2), type="l", lwd=3);
        }
        mids=1:num_grps;

        yticks=seq(min_95ci, max_95ci, length.out=5);
        cat("Y ticks:\n");
        print(yticks);
        signf_digits=max(ceiling(abs(log10(abs(yticks)))));
        yticks=signif(yticks, signf_digits);

        axis(side=2, at=yticks, labels=sprintf("%3.2f", yticks), cex.axis=.5, las=2);
        title(ylab=paste("Mean ", stat_name, "\nwith Bootstrapped\n95% CI", sep=""));
        title(main=title, cex.main=1.5);
        title(main="with Wilcoxon rank sum test (difference between group means) p-values",
                line=.25, cex.main=.7, font.main=3);

        bar_width=mean(diff(mids));
	if(is.na(bar_width)){
		bar_width=1;
	}
	
        qbw=bar_width/6;

	plot_bottom=par()$usr[3];

        # Label x-axis
        text(mids-par()$cxy[1]/2, rep(plot_bottom+6*-par()$cxy[2]/2, num_grps),
                group_names, srt=-45, xpd=T, pos=4, font=2,
                cex=min(c(1,.7*bar_width/par()$cxy[1])));

        # Scatter
        if(samp_gly){
                for(grp_ix in 1:num_grps){
                        grpnm=group_names[grp_ix];
                        pts=stat[grps[[grpnm]]];
                        numpts=length(pts);
                        points(
                                #rep(mids[grp_ix], numpts),
                                mids[grp_ix]+rnorm(numpts, 0, bar_width/10),
                                pts, col="darkblue", cex=.5, type="p");
                }
        }

        # label CI's
        for(grp_ix in 1:num_grps){
                if(samp_size[grp_ix]>=40){
                        points(
                                c(mids[grp_ix]-qbw, mids[grp_ix]+qbw),
                                rep(ci95[grp_ix, 2],2), type="l", col="blue");
                        points(
                                c(mids[grp_ix]-qbw, mids[grp_ix]+qbw),
                                rep(ci95[grp_ix, 1],2), type="l", col="blue");
                        points(
                                rep(mids[grp_ix],2),
                                c(ci95[grp_ix, 1], ci95[grp_ix,2]), type="l", col="blue");
                }
        }

        # label sample size
        for(grp_ix in 1:num_grps){
                text(mids[grp_ix], plot_bottom+3*-par()$cxy[2]/2, paste("mean =", round(grp_means[grp_ix], 2)),
                        cex=.95, xpd=T, font=3, adj=c(.5,-1));

                text(mids[grp_ix], plot_bottom+4*-par()$cxy[2]/2, paste("n =",samp_size[grp_ix]),
                        cex=.85, xpd=T, font=3, adj=c(.5,-1));
        }

        connect_significant=function(A, B, ypos, pval){
                abline(h=ypos);
        }

        sigchar=function(x){
                if(x<=.0001){
                        return("***");
                }else if(x<=.001){
                        return("**");
                }else if(x<=.01){
                        return("*");
                }else{
                        return("");
                }
        }

        row_ix=1;
        for(i in 1:(num_grps-1)){

                pvalrow=pval_mat[i,];
                #print(pvalrow);

                signf_pairs=(pvalrow<alpha);
                if(any(signf_pairs)){
                        signf_grps=which(signf_pairs);
                        cat("Pairs: ", i, " to:\n");
                        print(signf_grps);

                        y_offset=plotdatamax+horiz_spacing*row_ix;

                        # Draw line between left/reference to each paired signf grp
                        points(c(
                                mids[i], mids[max(signf_grps)]),
                                rep(y_offset,2),
                                type="l", lend="square"
                        );

                        # Mark left/ref group
                        points(
                                rep(mids[i],2),
                                c(y_offset,y_offset-horiz_spacing/4),
                                type="l", lwd=3, lend="butt");

                        # Mark each signf paired reference group
                        for(pair_ix in signf_grps){
                                points(
                                        rep(mids[pair_ix],2),
                                        c(y_offset,y_offset-horiz_spacing/4),
                                        type="l", lwd=1, lend="butt");


                                # label pvalue
                                paird_pval=sprintf("%5.4f", pvalrow[pair_ix]);
                                text(mids[pair_ix], y_offset, paird_pval,
                                        adj=c(.5, -1), cex=.7);
                                text(mids[pair_ix], y_offset, sigchar(pvalrow[pair_ix]),
                                        adj=c(.5, -1.25), cex=1);
                        }

                        row_ix=row_ix+1;

                }
	}

	return(signf);
}

calc_longitudinal_stats=function(offset_rec, alr_cat_val){

	cat("Calculating Longitudinal Statistics...\n");

	l_min=function(x, y){
		res=min(y);
		return(res);
	}
	l_max=function(x, y){
		res=max(y);
		return(res);
	}
	l_median=function(x, y){
		res=median(y);
		return(res);
	}
	l_mean=function(x, y){
		res=mean(y);
		return(res);
	}
	l_stdev=function(x, y){
		res=sd(y);
		return(res);
	}
	l_range=function(x, y){
		r=range(y);
		span=abs(r[1]-r[2]);
		return(span);
	}
	l_N=function(x, y){
		return(length(x));
	}
	l_last_time=function(x, y){
		return(max(x));
	}

	l_volatility=function(x, y){
		if(length(x)>1){
			fit=lm(y~x);
			sumfit=summary(fit);
			vol=sd(sumfit$residuals);
			return(vol);
		}else{
			return(NA);
		}
	}
	l_slope=function(x, y){
		if(length(x)>1){
			fit=lm(y~x);
			slope=fit$coefficients["x"];
			return(slope);
		}else{
			return(NA);
		}
	}

	l_time_wght_avg=function(x, y){
		npts=length(x);
		if(npts>1){

			tot_avg=0;
			
			for(i in 1:(npts-1)){
				avg_val=(y[i]+y[i+1])/2;
				duration=(x[i+1]-x[i]);
				tot_avg=tot_avg+(avg_val*duration);
			}

			norm=tot_avg/(x[npts]-x[1]);	
			return(norm);
		}else{
			return(NA);
		}	
	}

	l_time_at_max=function(x, y){
		max_val=max(y);
		ix=min(which(y==max_val));
		return(x[ix]);
	}

	l_time_at_min=function(x, y){
		min_val=min(y);
		ix=min(which(y==min_val));
		return(x[ix]);
	}

	l_time_closest_to_t0=function(x, y){
		starty=y[1];
		y=y[-1];
		dist=abs(y-starty);
		min_dist=min(dist);
		min_ix=min(which(min_dist==dist));
		return(x[min_ix+1]);
	}

	l_time_furthest_fr_t0=function(x, y){
		starty=y[1];
		y=y[-1];
		dist=abs(y-starty);
		max_dist=max(dist);
		max_ix=min(which(max_dist==dist));
		return(x[max_ix+1]);
	}

	l_start_end_diff=function(x,y){
		start=y[1];
		end=tail(y,1);
		return(end-start);
	}

	l_convexcave=function(x, y){
		num_pts=length(x);
		# y=mx+b
		# b=y-mx

		m=(y[num_pts]-y[1])/(x[num_pts]-x[1]);
		b=y[1]-m*x[1];

		lvl_y=numeric(num_pts);
		for(i in 1:num_pts){
			lvl_y[i]=y[i]-(m*x[i]+b);
		}

		cum_sum=0;
		for(i in 1:(num_pts-1)){
			dx=x[i+1]-x[i];
			avgy=(lvl_y[i+1]+lvl_y[i])/2;
			avg=dx*avgy/2;
			cum_sum=cum_sum+avg;
		}
		vexcav=cum_sum/(x[num_pts]-x[1]);
			
		return(vexcav);
	}

	# statistic:
	#    ALR:
	#       individual:	 

	stat_name=c(
		"min", "max", "range",
		"volatility", "slope", "time_wght_avg",
		"time_at_max", "time_at_min",
		"time_closest_to_t0", "time_furthest_fr_t0",
		"start_end_diff",
		"convexcave"
	);

	results=list();

	alrcat=colnames(alr_cat_val);
	individuals=as.character(offset_rec[["SubjectIDs"]]);

	cat("\n");
	cat("Individuals:\n");
	print(individuals);
	
	cat("\n");
	cat("Categories:\n");
	print(alrcat);

	num_cat=ncol(alr_cat_val);
	num_ind=length(individuals);


	for(stat_ix in stat_name){

		results[[stat_ix]]=list();

		tmp_mat=matrix(NA, nrow=num_ind, ncol=num_cat);
		rownames(tmp_mat)=individuals;
		colnames(tmp_mat)=alrcat;

		for(cat_ix in alrcat){

			for(ind_ix in individuals){
					
				indv_offsets=offset_rec[["OffsetsBySubject"]][[ind_ix]];
                		samp_ids=rownames(indv_offsets);

				time=indv_offsets[,"Offsets"];
				val=alr_categories_val[samp_ids, cat_ix];

				funct_name=paste("l_", stat_ix, sep="");

				call_res=do.call(funct_name, list(x=time, y=val));

				tmp_mat[ind_ix, cat_ix]=call_res;

			}
			
		}

		results[[stat_ix]]=tmp_mat;

	}

	return(results);
}

regress_longitudinal_stats=function(longit_stats, variable_list, factors){

	subj_ids=rownames(factors);
	stat_names=names(longit_stats);
	num_variables=length(variable_list);

	all_regr_res=list();
	df.factors=as.data.frame(factors);

	# Get coefficient names from model matrix
	resp=rep(0, length(subj_ids));
	form_str=paste("resp ~ ", paste(variable_list, collapse=" + "), sep="");
	mm=model.matrix(as.formula(form_str), data=df.factors);
	coef_names=setdiff(colnames(mm), "(Intercept)");
	cat("Coefficient Names:\n");
	print(coef_names);

	for(stat_ix in stat_names){
		cat("Working on: ", stat_ix, "\n");
		
		stat_mat=longit_stats[[stat_ix]];
		cat_names=colnames(stat_mat);
		num_cat=length(cat_names);

		pval_mat=matrix(NA, nrow=length(coef_names), ncol=num_cat);
		rownames(pval_mat)=coef_names;
		colnames(pval_mat)=cat_names;
		coef_mat=pval_mat;


		for(cat_ix  in cat_names){
		
			cat("Regression Analsis for: ", stat_ix, " of ", cat_ix, "\n", sep="");	

			resp=stat_mat[subj_ids, cat_ix, drop=F];
			
			form_str=paste("resp ~ ", paste(variable_list, collapse=" + "), sep="");

			cat("Model: ", form_str, "\n");
			fit=lm(as.formula(form_str), data=df.factors);
			sumfit=summary(fit);
			
			kept=intersect(coef_names, rownames(sumfit$coefficients));
			pval_mat[kept, cat_ix]=sumfit$coefficients[kept, "Pr(>|t|)"];
			coef_mat[kept, cat_ix]=sumfit$coefficients[kept, "Estimate"];
	
		}

		reg_res=list();
		reg_res[["pval"]]=pval_mat;
		reg_res[["coef"]]=coef_mat;

		all_regr_res[[stat_ix]]=reg_res;
	}

	return(all_regr_res);
}

create_GrpToSbj_map=function(subjects_arr, groups_arr){
	map=list();
	if(length(subjects_arr) != length(groups_arr)){
		cat("Error: Subj/Grp array lengths do no match.\n");
		cat("Subjects:\n");
		print(subjects_arr);
		cat("\nGroups:\n");
		print(groups_arr);
		quit(-1);
	}
	uniq_grps=sort(unique(groups_arr));
	for(grp_ix in uniq_grps){
		map[[grp_ix]]=unique(sort(as.character(subjects_arr[groups_arr==grp_ix])));
	}
	return(map);
}

##############################################################################
##############################################################################

# Load Factor File
cat("Loading Factor File:\n");
factor_info=load_factors(FactorFile);	
factor_samp_ids=rownames(factor_info);

# Load summary file table counts 
cat("Loading summary table...\n");
counts=load_summary_file(SummaryFile);

# Remove zero count samples
cat("Looking for zero count samples...\n");
tot=apply(counts, 1, sum);
nonzero=tot>0;
if(!(all(nonzero))){
	cat("WARNING: Zero count samples found:\n");
	samp_names=rownames(counts);
	print(samp_names[!nonzero]);
	cat("\n");
	cat("Removing...\n");
	counts=counts[nonzero,,drop=F];
}

num_st_categories=ncol(counts);
num_st_samples=nrow(counts);

cat("Num Categories: ", num_st_categories, "\n");
cat("   Num Samples: ", num_st_samples, "\n");

##############################################################################

# Shorten cateogry names
if(ShortenCategoryNames!=""){
	full_names=colnames(counts);
	splits=strsplit(full_names, ShortenCategoryNames);
	short_names=character();
	for(i in 1:length(full_names)){
		short_names[i]=tail(splits[[i]], 1);
		short_names[i]=gsub("_unclassified$", "_uc", short_names[i]);
		short_names[i]=gsub("_group", "_gr", short_names[i]);
		short_names[i]=gsub("\\[", "", short_names[i]);
		short_names[i]=gsub("\\]", "", short_names[i]);
	}
	colnames(counts)=short_names;
	cat("Names have been shortened.\n");
}else{
	cat("Keeping original category names...\n");
}

##############################################################################

if(NumALRPredictors >= num_st_categories){
	NumALRPredictors= (num_st_categories-1);
	cat("Number of taxa to work on was changed to: ", NumALRPredictors, "\n");
}

##############################################################################

# Reconciling Sample IDs
cat("Reconciling Summary Table samples with Factor File.\n");
sumtab_samp_ids=rownames(counts);
shared_samp_ids=sort(intersect(factor_samp_ids, sumtab_samp_ids));

# Remove exclusive samples
counts=counts[shared_samp_ids,,drop=F];
factor_info=factor_info[shared_samp_ids,,drop=F];

if(GroupCol==""){
	group_names=rep("All", nrow(factor_info));
	GroupCol="Group";
}else{
	group_names=factor_info[,GroupCol,drop=F];
}
grp_to_sbj=create_GrpToSbj_map(factor_info[,SubjectIDCol], group_names);
unique_group_names=sort(unique(group_names));
print(grp_to_sbj);

# Extract offset info
offset_rec=extract_offset(factor_info, SubjectIDCol, TimeOffsetCol);
print(offset_rec);

# Normalize
cat("Normalizing counts...\n");
counts=counts+.5;
normalized=normalize(counts);

cat("Reordering normalized...\n");
mean_norm=apply(normalized, 2, mean);
ord_ix=order(mean_norm, decreasing=T);
normalized=normalized[,ord_ix, drop=F];
counts=counts[,ord_ix, drop=F];

if(UseRemaining){
	category_names=colnames(counts);	
	uc_cat_names=toupper(category_names);
	remaining_ix=which(uc_cat_names=="REMAINDER" | uc_cat_names=="REMAINING");
	if(length(remaining_ix)!=1){
		cat("*******************************************************\n");
		cat("*  WARNING:  Could not identify remaining column.     *\n");
		cat("*******************************************************\n");
		UseRemaining=F;
	}else{
		cat("Remaining original column: ", remaining_ix, "\n");
		# Take out "remaining" column so it doesn't end up as a top column
		normalized_remaining_col_dat=normalized[,remaining_ix, drop=F];
		normalized=normalized[,-remaining_ix];
	}
}else{
	cat("Assuming no categories called 'remainder' or 'remaining'\n");
}

##############################################################################

# Open main output file
OutputRoot=paste(OutputRoot, ".", offset_rec[["RangeTag"]], sep="");
pdf(paste(OutputRoot, ".alr_ts.pdf", sep=""), height=14, width=8.5);

##############################################################################

cat("\n");
cat("Extracting Top categories: ", NumALRPredictors, " from amongst ", ncol(normalized), "\n", sep="");

cat_abundances=extract_top_categories(normalized, NumALRPredictors);
resp_alr_struct=additive_log_rato(cat_abundances);
alr_categories_val=resp_alr_struct$transformed;
alr_cat_names=colnames(alr_categories_val);

plot_text(c(
	paste("ALR Categories (Top ", NumALRPredictors, ")", sep=""),
	capture.output(print(alr_cat_names))
));

##############################################################################

alr_min=apply(alr_categories_val, 2, min);
alr_max=apply(alr_categories_val, 2, max);
alr_med=apply(alr_categories_val, 2, median);

cat("ALR Minimums:\n");
print(alr_min);

cat("ALR Maximums:\n");
print(alr_max);

cat("ALR Medians:\n");
print(alr_med);

cat("ALR Categories:\n");
print(alr_cat_names);

offset_ranges=c(offset_rec[["Earliest"]], offset_rec[["Latest"]]);
cat("Offset Range:\n");
print(offset_ranges);
	
subj_ids=offset_rec[["SubjectIDs"]];
num_subj=length(subj_ids);
subj_col=get_colors(num_subj);
names(subj_col)=subj_ids;

plots_per_page=6;
num_groups=length(unique_group_names);

par(mar=c(2,4,1,1));
par(oma=c(.5,.5,5,.5));

cat_loess=list();

for(cat_ix in alr_cat_names){
	cat("Plotting: ", cat_ix, "\n");
	alr_range=c(alr_min[cat_ix], alr_max[cat_ix]);

	# Plotting individuals by group
	cat("Plotting Individuals by Group:\n");
	for(grp_ix in unique_group_names){
		cat("For group: ", grp_ix, "\n");
	
		
		par(mfrow=c(plots_per_page, 1));	
		plot_ix=0;
		for(subj_ix in grp_to_sbj[[grp_ix]]){
			plot_alr_time_indv(cat_ix, subj_ix, offset_rec, alr_categories_val, 
				offset_ranges, alr_range, alr_med[cat_ix], subj_col[subj_ix]);
			plot_ix=plot_ix+1;
			if(plot_ix==plots_per_page){
				mtext(cat_ix, side=3, outer=T, font=2, cex=2, line=2);
				mtext(paste(GroupCol, ": ", grp_ix, sep=""), 
					side=3, outer=T, font=1, cex=1, line=.25);
				plot_ix=0;
			}
		}
		if(plot_ix<plots_per_page){
			mtext(cat_ix, side=3, outer=T, font=2, cex=2, line=2);
			mtext(paste(GroupCol, ": ", grp_ix, sep=""), 
				side=3, outer=T, font=1, cex=1, line=.25);
		}
	}

	# Plot individuals in single plot
	cat("Plotting Individuals in Single Plot:\n");
	par(mfrow=c(num_groups, 1));
	for(grp_ix in unique_group_names){
		cat("For group: ", grp_ix, "\n");

		subj_in_grp=grp_to_sbj[[grp_ix]];
		plot_alr_time_grpd(cat_ix, subj_in_grp, grp_ix, GroupCol,  offset_rec, alr_categories_val,
			offset_ranges, alr_range, alr_med[cat_ix], subj_col);
	}
	mtext(cat_ix, side=3, outer=T, font=2, cex=2, line=2);

	# plot loess
	cat("Plotting Loess:\n");
	par(mfrow=c(num_groups, 1));
	grp_loess=list();
	for(grp_ix in unique_group_names){
		cat("For group: ", grp_ix, "\n");

		subj_in_grp=grp_to_sbj[[grp_ix]];
		grp_loess[[grp_ix]]=compute_and_plot_loess(cat_ix, subj_in_grp, grp_ix,
			GroupCol, offset_rec, alr_categories_val,
			offset_ranges, alr_range, alr_med[cat_ix], subj_col);
	}
	#plot_grp_loess(grp_loess);
	mtext(cat_ix, side=3, outer=T, font=2, cex=2, line=2);

	cat_loess[[cat_ix]]=grp_loess;

}

#print(cat_loess);

# Plot loess thumbnails

num_rows_pp=5;
par(oma=c(6,6,4,1));
par(mar=c(.5,.5,.5,.5));
par(mfrow=c(num_rows_pp,num_groups));

rownum=1;
last_cat=tail(alr_cat_names,1);
for(cat_ix in alr_cat_names){
	alr_range=c(alr_min[cat_ix], alr_max[cat_ix]);

	colnum=1;
	for(grp_ix in unique_group_names){

		plot(0, type="n", ylim=alr_range, xlim=offset_ranges, ylab=cat_ix, xaxt="n", yaxt="n");

		if(min(offset_ranges)<0){
			abline(v=0, col="blue", lty=3, lwd=3, lend="butt");
		}

		abline(h=alr_med[cat_ix], col="grey", lty="dotdash");
		grp_loess=cat_loess[[cat_ix]][[grp_ix]];
		x=grp_loess[,"x"];
		y=grp_loess[,"y"];
		points(x,y, type="l", col="blue");

		# Label bottom
		if(rownum==num_rows_pp || cat_ix==last_cat){
			axis(side=1);
		}

		# Label top row with group IDs
		if(rownum==1){
			if(nchar(grp_ix)>10){
				cex_adj=20/nchar(grp_ix);
			}else{
				cex_adj=2;
			}
			axis(side=3, at=mean(offset_ranges), labels=grp_ix, cex.axis=cex_adj, 
				line=1, outer=T, tick=F);
		}

		# Label left side with category and alr scale
		if(colnum==1){
			if(nchar(cat_ix)>14){
				cex_adj=2*14/nchar(cat_ix);
			}else{
				cex_adj=2;
			}
			axis(side=2);
			axis(side=2, at=mean(alr_range), labels=cat_ix, cex.axis=cex_adj, line=2, outer=T, tick=F);
		}

		colnum=colnum+1;
	}

	if(rownum==num_rows_pp){
		rownum=0;
	}
	rownum=rownum+1;

}

###############################################################################


long_stats=calc_longitudinal_stats(offset_rec, alr_categories_val);
print(long_stats);
stat_names=names(long_stats);

for(stat_ix in stat_names){
	mat=long_stats[[stat_ix]];
	paint_matrix(mat, stat_ix);
}

stat_table=matrix(NA, nrow=0, ncol=7);
colnames(stat_table)=c(
	"Statistic", 
	"Category", 
	"Gr1",
	"mean(Gr1)",
	"Gr2",
	"mean(Gr2)",
	"p-value");

for(stat_ix in stat_names){

	num_grps=length(unique_group_names);
	num_alr_cat=ncol(alr_categories_val);

	grp_mat=matrix(NA, nrow=num_grps, ncol=num_alr_cat);
	rownames(grp_mat)=unique_group_names;
	colnames(grp_mat)=colnames(alr_categories_val);

	for(grp_ix in unique_group_names){
		grp_members=grp_to_sbj[[grp_ix]];
		grp_mat[grp_ix,]=apply(
			long_stats[[stat_ix]][grp_members,,drop=F], 2, 
			function(x){mean(x, na.rm=T)});
	}

	paint_matrix(grp_mat, paste("mean(", stat_ix, ") for Grouping: ", GroupCol, sep=""));

	plots_pp=5;
	par(mfrow=c(plots_pp,1));
	par(oma=c(0,0,2,0));
	label_oma=F
	plot_ix=0;
	for(cat_ix in alr_cat_names){
		signf=plot_barplot_wsignf_annot(
			title=cat_ix,
			long_stats[[stat_ix]][,cat_ix],
			grp_to_sbj,
			alpha=Alpha
		);
		
		if(length(signf)){
			num_sig_rows=nrow(signf);
			newrows=cbind(rep(stat_ix, num_sig_rows), rep(cat_ix, num_sig_rows), signf);
			stat_table=rbind(stat_table, newrows);
		}

		plot_ix=plot_ix+1;
		if(plot_ix==plots_pp){
			mtext(stat_ix, outer=T, cex=1.5, col="blue", font=2);
			plot_ix=0;
		}
	}
	if(plot_ix!=plots_pp){
		mtext(stat_ix, outer=T, cex=1.5, col="blue", font=2);
	}
}

par(mfrow=c(1,1));
plot_text(c(
	"Explanation of Longitudinal Statistics:",
	"",
	"These stats are calculated independent of the time component:",
	"  min: minimum value", 
	"  max: maximum value", 
	"  median: median value",
	"  mean: average value", 
	"  stdev: standard deviation of values", 
	"  range: (max-min) of values",
	"  N: number of samples",
	"  last_time: last recorded sample",
	"",
	"These stats are based on fitting a linear model (value = intercept + time + error):",
	"  volatility:  standard deviation of residuals", 
	"  slope: slope of line", 
	"",
	"These stats take time into account:",
	"  time_weighted_average: average between two time points weighted by time between points",
	"  time_at_max: time at value's maximum", 
	"  time_at_min: time at value's minimum",
	"  time_closest_to_t0: time when value is closest to first time value", 
	"  time_furthest_fr_t0: time when value is furthest from first time value",
	"  start_end_diff: difference between first and last time points (end-start)",
	"  convexcave: if curve is convex or conave, stat will be <0 or >0, respectively",
	"  	A concave curve looks like a hill, a convex curve looks like a pit" 
	
));

###############################################################################

sigchar=function(x){
	if(x<=.0001){
		return("***");
	}else if(x<=.001){
		return("**");
	}else if(x<=.01){
		return("*");
	}else if(x<=.05){
		return(".");
	}else{
		return("");
	}
}

if(nrow(stat_table)>0){

	pvals=as.numeric(stat_table[,"p-value"]);
	signf=sapply(pvals, sigchar);
	stat_table=cbind(stat_table, signf);
	row_idx_str=paste(1:nrow(stat_table), ".", sep="");

	options(width=1000);

	#------------------------------------------------
	cat("Ordering by P-value...\n");
	# First order by pvalue, then do stable sort
	stat_order_ix=order(pvals);
	stat_table=stat_table[stat_order_ix,];
	print(stat_table);

	#------------------------------------------------
	cat("Ordering by Stat Name...\n");

	stat_order_ix=order(stat_table[,"Statistic"], method="shell");
	out_stat_table=stat_table[stat_order_ix,];
	rownames(out_stat_table)=row_idx_str;
	out=capture.output(print(out_stat_table, quote=F));

	cat("Writing by statistic:\n");
	plot_text(c(
		"By Statistic:",
		"",
		out));

	#------------------------------------------------
	cat("Ordering by ALR Category...\n");

	stat_order_ix=order(stat_table[,"Category"], method="shell");
	out_stat_table=stat_table[stat_order_ix,];
	rownames(out_stat_table)=row_idx_str;
	out=capture.output(print(out_stat_table, quote=F));

	cat("Writing by category:\n");
	plot_text(c(
		"By Category:",
		"",
		out));

	#------------------------------------------------
	cat("Outputing by P-Value...\n");

	out_stat_table=stat_table;
	rownames(out_stat_table)=row_idx_str;
	out=capture.output(print(out_stat_table, quote=F));

	cat("Writing by p-value:\n");
	plot_text(c(
		"By P-value:",
		"",
		out));

	#------------------------------------------------

	cat("Writing stat table...\n");
	write.table(out_stat_table, file=paste(OutputRoot,".lngt_stats.comp.tsv", sep=""), quote=F, sep="\t");
}else{
	plot_text(c(
		"No Significant Differences Identified Between Groups."
	));
	fh=file(paste(OutputRoot,".lngt_stats.comp.tsv", sep=""), "w");
	cat(file=fh, "Sorry.  No significant differences identified between groups.\n");
	close(fh);
}

###############################################################################
# Export stats table


###############################################################################


#cat("--------------------------------------------------------------------------------\n");
#cat("Fitting regression models with longitudinal stats as response:\n");
#cat("Factor File: ", FactorFile, "\n");
#cat("Model File: ", ModelFile, "\n");


if(ModelFile!=""){
	model_var_list=load_list(ModelFile);
	cat("Model variables in: ", ModelFile, "\n");
}else{
	cat("Model File was not specified. Skipping analyses with factors.\n");
	quit(status=-1);
}

# collapse
kept_fact=factor_info[,c(SubjectIDCol,model_var_list)];
uniq_subj=sort(unique(factor_info[,SubjectIDCol]));

num_subj=length(uniq_subj);
num_model_var=length(model_var_list);

# Collapse factor info and only keep those with offset info
colpsd_factors=numeric();
uniq_subj=intersect(uniq_subj, offset_rec[["SubjectIDs"]]);
for(subj_id in uniq_subj){
	ix=(subj_id==factor_info[,SubjectIDCol]);
	first_found=min(which(ix));
	val=factor_info[first_found, model_var_list];
	colpsd_factors=rbind(colpsd_factors, val);
}

rownames(colpsd_factors)=uniq_subj;
colnames(colpsd_factors)=model_var_list;

regres=regress_longitudinal_stats(long_stats, model_var_list, colpsd_factors);
options(width=200);
print(regres);

stat_names=names(regres);

plot_heat_maps=function(coef_mat, pval_mat, title){

	title_page(title);
	
	paint_matrix(
		coef_mat, 
		title=paste(title, ": Coefficients", sep=""));

	min_coef=min(coef_mat);
	max_coef=max(coef_mat);

	paint_matrix(
		pval_mat, 
		title=paste(title, ": P-values", sep=""),
		plot_min=0, plot_max=1, high_is_hot=F);
	
	coef_at_10=mask_matrix(coef_mat, pval_mat, 0.1, 0);
	coef_at_05=mask_matrix(coef_mat, pval_mat, 0.05, 0);
	coef_at_01=mask_matrix(coef_mat, pval_mat, 0.01, 0);

	paint_matrix(
		coef_at_10,
		label_zeros=F, plot_min=min_coef, plot_max=max_coef,
		title=paste(title, ": Coeff with p-val < 0.10", sep=""));
	paint_matrix(
		coef_at_05,
		label_zeros=F, plot_min=min_coef, plot_max=max_coef,
		title=paste(title, ": Coeff with p-val < 0.05", sep=""));
	paint_matrix(
		coef_at_01,
		label_zeros=F, plot_min=min_coef, plot_max=max_coef,
		title=paste(title, ": Coeff with p-val < 0.01", sep=""));
	
}

for(stat_ix in stat_names){
	plot_heat_maps(
		regres[[stat_ix]][["coef"]],
		regres[[stat_ix]][["pval"]],
		stat_ix);
}

###############################################################################

regr_stat_summary=character();

for(stat_ix in stat_names){
	coef_mat=regres[[stat_ix]][["coef"]];
	pval_mat=regres[[stat_ix]][["pval"]];

	predictors=rownames(coef_mat);
	categories=colnames(coef_mat);

	for(pred_ix in predictors){
		for(cat_ix in categories){

			pval=pval_mat[pred_ix, cat_ix];
		
			if(!is.na(pval) && pval<=0.1){
				coef=signif(coef_mat[pred_ix, cat_ix], 3);
				pval=signif(pval_mat[pred_ix, cat_ix], 3);

				row_info=c(
					stat_ix, pred_ix, cat_ix, 
					coef,
					pval,
					sigchar(pval)	
				);

				regr_stat_summary=rbind(regr_stat_summary, row_info);
			}	
		}

	}
}


colnames(regr_stat_summary)=c("Statistic", "Predictor", "Category", "Coefficient", "P-value", "Significance");

pval_numeric=as.numeric(regr_stat_summary[,"P-value"]);
pval_ord_ix=order(pval_numeric);

regr_stat_summary=regr_stat_summary[pval_ord_ix,];
num_sigf_reg_assoc=nrow(regr_stat_summary);

#---------------------------------------------------------------

for(ordering in c("Statistic", "Predictor", "Category", "P-value")){

	if(ordering=="P-value"){
		out_stat_table=regr_stat_summary;
	}else{
		stat_order=order(regr_stat_summary[,ordering], method="shell");
		out_stat_table=regr_stat_summary[stat_order,];
	}
	rownames(out_stat_table)=paste(1:num_sigf_reg_assoc, ".", sep="");

	cat("Outputing by: ", ordering, "\n");

	out=capture.output(print(out_stat_table, quote=F));

	plot_text(c(
		paste("By ", ordering, ":", sep=""), 
		"",
		out));

}

###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
