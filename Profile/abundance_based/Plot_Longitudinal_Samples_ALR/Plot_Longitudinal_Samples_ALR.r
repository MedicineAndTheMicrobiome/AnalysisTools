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
	"end_offset", "e", 2, "numeric",

	"tag_name", "T", 2, "character"
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
	"	[-T <tag name>]\n",
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

if(length(opt$tag_name)){
        TagName=opt$tag_name;
        cat("Setting TagName Hook: ", TagName, "\n");
        setHook("plot.new",
                function(){
                        #cat("Hook called.\n");
                        if(par()$page==T){
                                oma_orig=par()$oma;
                                exp_oma=oma_orig;
                                exp_oma[1]=max(exp_oma[1], 1.5);
                                par(oma=exp_oma);
                                mtext(paste("[", TagName, "]", sep=""), side=1, line=exp_oma[1]-1,
                                        outer=T, col="steelblue4", font=2, cex=.8, adj=.97);
                                par(oma=oma_orig);
                        }
                },
                "append");

}else{
        TagName="";
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

	offset_span=offset_range[2]-offset_range[1];

	plot(0, type="n", ylim=alr_range, xlim=c(offset_range[1], offset_range[2]+0.2*offset_span), 
		ylab=paste(grouping_name, ": ", grouping));
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

		text(tail(x_val, 1), tail(y_val, 1), tar_subj, cex=1.2, adj=c(-0.25, .5));
	}

}

compute_and_plot_loess=function(tar_cat, subj_arr, grouping, grouping_name, offsets_rec, alr_categories_val, 
	offset_range, alr_range, alr_med, subj_cols, grp_cols){

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

	offsets_range=range(offsets_rec[["Offsets"]]);
	num_offsets=offsets_rec[["NumUniqOffsets"]];
	loess_fit_x=seq(offsets_range[1], offsets_range[2], length.out=num_offsets*2);
	loess_fit_y=tryCatch({
			predict(loess_res, loess_fit_x);
		}, error=function(e){
			return(NULL);
		});
	if(is.null(loess_fit_y)){
		loess_fit_x=all_x;
		loess_fit_y=all_y;
	}
	

	grp_loess=cbind(loess_fit_x, loess_fit_y);
	colnames(grp_loess)=c("x", "y");

	points(grp_loess[,"x"], grp_loess[,"y"], type="l", lwd=3, col=grp_cols[grouping]);
	#points(grp_loess[,"x"], grp_loess[,"y"], type="l", lwd=.4, col="black");

	return(grp_loess);

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
	group_names=factor_info[,GroupCol];
}

subject_grouping_rec=create_GrpToSbj_map(factor_info[,SubjectIDCol], group_names);
unique_group_names=as.character(sort(unique(group_names)));

cat("Unique Group Names:\n");
print(unique_group_names);

grp_to_sbj=subject_grouping_rec[["GrpToSbj"]];
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

grp_col=terrain.colors(num_groups+1)[1:num_groups];
names(grp_col)=unique_group_names;

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
			offset_ranges, alr_range, alr_med[cat_ix], subj_col, grp_col);
	}
	#plot_grp_loess(grp_loess);
	mtext(cat_ix, side=3, outer=T, font=2, cex=2, line=2);

	# plot multiple overlapping
	cat("Plotting Overlapping Loess:\n");
	print(grp_loess);
	plot(0, type="n", xlim=offset_ranges, ylim=alr_range, xlab="", ylab="", main=cat_ix);
	abline(h=alr_med[cat_ix], col="grey", lty="dotdash");
	for(grp_ix in unique_group_names){
		points(grp_loess[[grp_ix]][,"x"], grp_loess[[grp_ix]][,"y"], type="l", lwd=3, col=grp_col[grp_ix]);
	}
	
	# plot legend
	plot(0, type="n", bty="n", xlab="", ylab="", main="", xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1));
	legend(0, 1, fill=grp_col, legend=names(grp_col), bty="n");	

	

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
		points(x,y, type="l", col=grp_col[grp_ix], lwd=3);

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

# Plot combined loess

num_rows_pp=5;
par(oma=c(6,6,4,1));
par(mar=c(.5,.5,.5,.5));
par(mfrow=c(num_rows_pp, 1));

rownum=1;
for(cat_ix in alr_cat_names){
	alr_range=c(alr_min[cat_ix], alr_max[cat_ix]);
	plot(0, type="n", ylim=alr_range, xlim=offset_ranges, ylab=cat_ix, xaxt="n", yaxt="n");

	for(grp_ix in unique_group_names){

		abline(h=alr_med[cat_ix], col="grey", lty="dotdash");

		grp_loess=cat_loess[[cat_ix]][[grp_ix]];
		x=grp_loess[,"x"];
		y=grp_loess[,"y"];
		points(x,y, type="l", col=grp_col[grp_ix], lwd=3);

		# Label bottom
		if(rownum==num_rows_pp || cat_ix==last_cat){
			axis(side=1);
		}

		# Label left side with category and alr scale
		if(nchar(cat_ix)>14){
			cex_adj=2*14/nchar(cat_ix);
		}else{
			cex_adj=2;
		}
		axis(side=2);
		axis(side=2, at=mean(alr_range), labels=cat_ix, cex.axis=cex_adj, line=2, outer=T, tick=F);

	}

	if(rownum==num_rows_pp){
		rownum=0;
	}
	rownum=rownum+1;

}

###############################################################################

# Calculate longitudinal stats
long_stats=calc_longitudinal_stats(offset_rec, alr_categories_val);
stat_names=names(long_stats);

# Plot heatmaps
for(stat_ix in stat_names){
	mat=long_stats[[stat_ix]];
	paint_matrix(mat, stat_ix, deci_pts=-1);
}

# Plot and compute pairwise group comparisons
stat_table_grp_cmp=plot_pairwise_grp_comparisons(long_stats, subject_grouping_rec, plots_pp=3);

###############################################################################

# Output significant stats
output_stat_table_alternate_ordering(stat_table_grp_cmp, OutputRoot);

###############################################################################
###############################################################################
# Work on regressing longitudinal stats as response to factors/group info

if(ModelFile!=""){
	model_var_list=load_list(ModelFile);
	model_var_list=c(model_var_list, GroupCol);
	cat("Model variables in: ", ModelFile, "\n");
}else{
	cat("Model File was not specified. Skipping analyses with factors.\n");
	quit();
}


cat("Collapsing Factors...\n");
colpsd_factors=collapse_factors(factor_info, SubjectIDCol, model_var_list);

cat("Regressing Longitudinal Stats...\n");
regres=regress_longitudinal_stats(long_stats, model_var_list, colpsd_factors);

options(width=200);
print(regres);

cat("Generating Heatmaps by Stat...\n");
title_page(paste(
	"Response:\nLongitudinal Stats\n\nPredictors:\n", paste(model_var_list, collapse="\n"), sep=""));

stat_names=names(regres);
for(stat_ix in stat_names){
	plot_heat_maps(
		regres[[stat_ix]][["coef"]],
		regres[[stat_ix]][["pval"]],
		stat_ix);
}

cat("Summarizing Regression Results into Table...\n");
regr_stat_summary=summarize_regression_results(regres, stat_names);
num_sigf_reg_assoc=nrow(regr_stat_summary);

cat("Writing Stats by Alternative Ordering...\n");
output_long_regression_stats_w_alt_ordering(regr_stat_summary);

###############################################################################

# Output descript of stats
par(mfrow=c(1,1));
plot_text(longit_stat_descriptions);

###############################################################################

cat("Done.\n");

dev.off();
print(warnings());
q(status=0);
