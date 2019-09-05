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

	"offset_file", "t", 1, "character",
	"output_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

NUM_TOP_PRED_CAT=20;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table>\n",
	"\n",
	"	[-v <number of top predictor (as ALR) categories to include, default=", NUM_TOP_PRED_CAT, ">]\n",
	"	[-R (pay attention to 'remaining' category)]\n",
	"	[-x <shorten category names, with separator in double quotes (default=\"\")>]\n",
	"\n",
	"	-t <offset file>\n",
	"	-o <output root>\n",
	"\n",
	"\n",	
	"\n", sep="");

if(
	!length(opt$summary_file) || 
	!length(opt$offset_file) || 
	!length(opt$output_root)
	){
	cat(usage);
	q(status=-1);
}

# Required
SummaryFile=opt$summary_file;
OffsetFile=opt$offset_file;
OutputRoot=opt$output_root;

# Optional, i.e. with defaults
NumALRPredictors=NUM_TOP_PRED_CAT;
UseRemaining=F;
ShortenCategoryNames="";

if(length(opt$num_top_pred)){
	NumALRPredictors=opt$num_top_pred;
}

if(length(opt$contains_remaining)){
	UseRemaining=T;
}

if(length(opt$shorten_category_names)){
	ShortenCategoryNames=opt$shorten_category_names;
}

###############################################################################

input_param=capture.output({
	cat("\n");
	cat("Summary File: ", SummaryFile, "\n", sep="");
	cat("  Num Top ALR Predictors: ", NumALRPredictors, "\n", sep="");
	cat("  Contains 'Remaining': ", UseRemaining, "\n", sep="");
	cat("  Shorten Categories: ", ShortenCategoryNames, "\n", sep="");
	cat("\n");
	cat("Offset File: ", OffsetFile, "\n", sep="");
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

plot_text=function(strings){
	par(family="Courier");
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

plot_alr_time_indv=function(tar_cat, tar_subj, offsets_rec, alr_categories_val, 
	offset_range, alr_range, alr_med, col="black"){

	indv_offsets=offsets_rec[["OffsetsByIndiv"]][[tar_subj]];
	samp_ids=rownames(indv_offsets);

	cat("Working on: ", tar_cat, " for ", tar_subj, "\n");

	x_val=indv_offsets[,"Offsets"];
	y_val=alr_categories_val[samp_ids, tar_cat];

	plot(0, type="n", ylim=alr_range, xlim=offset_range, ylab=tar_subj);
	abline(h=alr_med, col="grey", lty="dotdash");
	points(x_val, y_val, type="l", lwd=5, col=col);
	points(x_val, y_val, type="l", lwd=.5, col="black");
	points(x_val, y_val, type="p", cex=3, col="black");

}

plot_alr_time_grpd=function(tar_cat, subj_arr, grouping, grouping_name, offsets_rec, alr_categories_val, 
	offset_range, alr_range, alr_med, subj_cols){

	plot(0, type="n", ylim=alr_range, xlim=offset_range, ylab=paste(grouping_name, ": ", grouping));
	abline(h=alr_med, col="grey", lty="dotdash");

	for(tar_subj in subj_arr){

		cat("Working on: ", tar_cat, " for ", tar_subj, "\n");

		indv_offsets=offsets_rec[["OffsetsByIndiv"]][[tar_subj]];
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

	all_x=numeric();
	all_y=numeric();
	for(tar_subj in subj_arr){

		cat("Working on: ", tar_cat, " for ", tar_subj, "\n");

		indv_offsets=offsets_rec[["OffsetsByIndiv"]][[tar_subj]];
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


##############################################################################
##############################################################################

# Open main output file
pdf(paste(OutputRoot, ".alr_ts.pdf", sep=""), height=14, width=8.5);

# Load offset file
offset_raw=load_offset(OffsetFile);
print(offset_raw);
offset_samp_ids=rownames(offset_raw[["matrix"]]);

# Load summary file table counts 
cat("\n");
cat("Loading summary table...\n");
counts=load_summary_file(SummaryFile);

# Remove zero count samples
tot=apply(counts, 1, sum);
nonzero=tot>0;
if(!(all(nonzero))){
	cat("WARNING: Zero count samples found:\n");
	samp_names=rownames(counts);
	print(samp_names[!nonzero]);
	cat("\n");
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
		short_names[i]=gsub("_unclassified$", "_uncl", short_names[i]);
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

sumtab_samp_ids=rownames(counts);
shared_samp_ids=sort(intersect(offset_samp_ids, sumtab_samp_ids));

counts=counts[shared_samp_ids,];
offset_raw[["matrix"]]=offset_raw[["matrix"]][shared_samp_ids,];

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

offset_info=group_offsets(offset_raw);

###############################################################################

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

offset_ranges=range(offset_info[["Offsets"]]);
cat("Offset Range:\n");
print(offset_ranges);
	
print(offset_info);
subj_ids=offset_info[["Individuals"]];
num_subj=length(subj_ids);
subj_col=get_colors(num_subj);
names(subj_col)=subj_ids;

plots_per_page=6;
group_name=offset_raw[["GroupID"]];
num_groups=length(offset_info[["Groups"]]);

par(mar=c(2,4,1,1));
par(oma=c(.5,.5,5,.5));

cat_loess=list();

for(cat_ix in alr_cat_names){
	cat("Plotting: ", cat_ix, "\n");
	alr_range=c(alr_min[cat_ix], alr_max[cat_ix]);

	# Plotting individuals by group
	for(grp_ix in offset_info[["Groups"]]){
		cat("For group: ", grp_ix, "\n");
	
		
		par(mfrow=c(plots_per_page, 1));	
		plot_ix=0;
		for(subj_ix in offset_info[["IndivByGrp"]][[grp_ix]]){
			plot_alr_time_indv(cat_ix, subj_ix, offset_info, alr_categories_val, 
				offset_ranges, alr_range, alr_med[cat_ix], subj_col[subj_ix]);
			plot_ix=plot_ix+1;
			if(plot_ix==plots_per_page){
				mtext(cat_ix, side=3, outer=T, font=2, cex=2, line=2);
				mtext(paste(group_name, ": ", grp_ix, sep=""), 
					side=3, outer=T, font=1, cex=1, line=.25);
				plot_ix=0;
			}
		}
		if(plot_ix<plots_per_page){
			mtext(cat_ix, side=3, outer=T, font=2, cex=2, line=2);
			mtext(paste(group_name, ": ", grp_ix, sep=""), 
				side=3, outer=T, font=1, cex=1, line=.25);
		}
	}

	# Plot individuals in single plot
	par(mfrow=c(num_groups, 1));
	for(grp_ix in offset_info[["Groups"]]){
		cat("For group: ", grp_ix, "\n");

		subj_in_grp=offset_info[["IndivByGrp"]][[grp_ix]];
		plot_alr_time_grpd(cat_ix, subj_in_grp, grp_ix, group_name,  offset_info, alr_categories_val,
			offset_ranges, alr_range, alr_med[cat_ix], subj_col);
	}
	mtext(cat_ix, side=3, outer=T, font=2, cex=2, line=2);

	# plot loess
	par(mfrow=c(num_groups, 1));
	grp_loess=list();
	for(grp_ix in offset_info[["Groups"]]){
		cat("For group: ", grp_ix, "\n");

		subj_in_grp=offset_info[["IndivByGrp"]][[grp_ix]];
		grp_loess[[grp_ix]]=compute_and_plot_loess(cat_ix, subj_in_grp, grp_ix,
			group_name, offset_info, alr_categories_val,
			offset_ranges, alr_range, alr_med[cat_ix], subj_col);
	}
	#plot_grp_loess(grp_loess);
	mtext(cat_ix, side=3, outer=T, font=2, cex=2, line=2);

	cat_loess[[cat_ix]]=grp_loess;

}

print(cat_loess);

num_rows_pp=6;
par(oma=c(6,6,4,1));
par(mar=c(.5,.5,.5,.5));
par(mfrow=c(num_rows_pp,num_groups));

rownum=1;
last_cat=tail(alr_cat_names,1);
for(cat_ix in alr_cat_names){
	alr_range=c(alr_min[cat_ix], alr_max[cat_ix]);

	colnum=1;
	for(grp_ix in offset_info[["Groups"]]){


		plot(0, type="n", ylim=alr_range, xlim=offset_ranges, ylab=cat_ix, xaxt="n", yaxt="n");
		abline(h=alr_med[cat_ix], col="grey", lty="dotdash");
		grp_loess=cat_loess[[cat_ix]][[grp_ix]];
		x=grp_loess[,"x"];
		y=grp_loess[,"y"];
		points(x,y, type="l", col="blue");

		# Label bottom
		if(rownum==num_rows_pp || cat_ix==last_cat){
			axis(side=1);
		}

		# Label to row with group IDs
		if(rownum==1){
			axis(side=3, at=mean(offset_ranges), labels=grp_ix, cex.axis=2, line=1, outer=T, tick=F);
		}

		# Label left side with category and alr scale
		if(colnum==1){
			axis(side=2);
			axis(side=2, at=mean(alr_range), labels=cat_ix, cex.axis=1.5, line=2, outer=T, tick=F);
		}

		colnum=colnum+1;
	}

	if(rownum==num_rows_pp){
		rownum=0;
	}
	rownum=rownum+1;

}

###############################################################################

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

	l_time_weighted_average=function(x, y){
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

	l_time_closest_to_start=function(x, y){
		starty=y[1];
		y=y[-1];
		dist=abs(y-starty);
		min_dist=min(dist);
		min_ix=min(which(min_dist==dist));
		return(x[min_ix+1]);
	}

	l_time_furthest_from_start=function(x, y){
		starty=y[1];
		y=y[-1];
		dist=abs(y-starty);
		max_dist=max(dist);
		max_ix=max(which(max_dist==dist));
		return(x[max_ix+1]);
	}

	# statistic:
	#    ALR:
	#       individual:	 

	stat_name=c(
		"min", "max", "median", "mean", "stdev", "range",
		"volatility", "slope", "time_weighted_average",
		"time_at_max", "time_at_min",
		"time_closest_to_start", "time_furthest_from_start"
	);

	results=list();

	alrcat=colnames(alr_cat_val);
	individuals=as.character(offset_rec[["Individuals"]]);

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
					
				indv_offsets=offset_rec[["OffsetsByIndiv"]][[ind_ix]];
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

long_stats=calc_longitudinal_stats(offset_info, alr_categories_val);
print(long_stats);

stat_names=names(long_stats);

for(stat_ix in stat_names){
	mat=long_stats[[stat_ix]];

	paint_matrix(mat, stat_ix);

}



###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
