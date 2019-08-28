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

get_colors=function(num_col, alpha=1){
        colors=hsv(seq(0,1,length.out=num_col+1), c(1,.5), c(1,.75,.5), alpha=alpha);
        color_mat_dim=ceiling(sqrt(num_col));
        color_pad=rep("grey", color_mat_dim^2);
        color_pad[1:num_col]=colors[1:num_col];
        color_mat=matrix(color_pad, nrow=color_mat_dim, ncol=color_mat_dim);
        colors=as.vector(t(color_mat));
        colors=colors[colors!="grey"];
}

plot_alr_time_indv=function(tar_cat, tar_subj, offsets_rec, alr_categories_val, offset_range, alr_range, col="black"){

	indv_offsets=offsets_rec[["OffsetsByIndiv"]][[tar_subj]];
	samp_ids=rownames(indv_offsets);

	cat("Working on: ", tar_cat, " for ", tar_subj, "\n");

	x_val=indv_offsets[,"Offsets"];
	y_val=alr_categories_val[samp_ids, tar_cat];

	plot(0, type="n", ylim=alr_range, xlim=offset_range, ylab=tar_subj);
	points(x_val, y_val, type="l", lwd=5, col=col);
	points(x_val, y_val, type="l", lwd=.5, col="black");
	points(x_val, y_val, type="p", cex=1, col="black");

}

plot_alr_time_grpd=function(tar_cat, subj_arr, grouping, grouping_name, offsets_rec, alr_categories_val, 
	offset_range, alr_range, subj_cols){

	plot(0, type="n", ylim=alr_range, xlim=offset_range, ylab=paste(grouping_name, ": ", grouping));

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
	offset_range, alr_range, subj_cols){

	plot(0, type="n", ylim=alr_range, xlim=offset_range, ylab=paste(grouping_name, ": ", grouping));

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

offset_raw=load_offset(OffsetFile);
print(offset_raw);

offset_info=group_offsets(offset_raw);

###############################################################################

alr_min=apply(alr_categories_val, 2, min);
alr_max=apply(alr_categories_val, 2, max);

cat("ALR Minimums:\n");
print(alr_min);

cat("ALR Maximums:\n");
print(alr_max);

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
				offset_ranges, alr_range, subj_col[subj_ix]);
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
			offset_ranges, alr_range, subj_col);
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
			offset_ranges, alr_range, subj_col);
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

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
