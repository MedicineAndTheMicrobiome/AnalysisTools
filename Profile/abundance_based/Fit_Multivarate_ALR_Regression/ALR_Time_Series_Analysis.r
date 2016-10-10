#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
options(useFancyQuotes=F);

NUM_VAR=20;

params=c(
	"summary_file", "s", 1, "character",
	"factors", "f", 1, "character",
	"num_variables", "p", 2, "numeric",
	"treatment", "t", 2, "character",
	"time", "m", 2, "character",
	"subject", "u", 2, "character",
	"outputroot", "o", 2, "character",
	"reset_time_pts", "r", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table>\n",
	"	-f <factors>\n",
	"\n",
	"    Column Names in Factor File:\n",
	"	-t <Treatment variable>\n",
	"	-m <tiMe variable>\n",
	"	-u <sUbject variable>\n",
	"\n",
	"	[-p <number of variables, default=", NUM_VAR, " >]\n",
	"	[-o <output filename root>]\n",
	"	[-r (reset time points so first time is 0)]\n",
	"\n",
	"This script will generate trellis plots grouped by treatment\n",
	"variable across all the subjects across time.\n",
	"\n");

if(!length(opt$summary_file) || !length(opt$factors)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputRoot=gsub(".summary_table.xls", "", opt$summary_file);
	OutputRoot=gsub(".summary_table.tsv", "", OutputRoot);
}else{
	OutputRoot=opt$outputroot;
}

if(!length(opt$num_variables)){
	NumVariables=NUM_VAR;
}else{
	NumVariables=opt$num_variables;
}

if(!length(opt$reset_time_pts)){
	ResetTimePoints=F;
}else{
	ResetTimePoints=T;
}

SummaryFile=opt$summary_file;
FactorsFile=opt$factors;

TreatmentVar=opt$treatment;
TimeVar=opt$time;
SubjectVar=opt$subject;

cat("\n");
cat("Summary File : ", SummaryFile, "\n", sep="");
cat("Factors File: ", FactorsFile, "\n", sep="");
cat("Number of Response Variables: ", NumVariables, "\n", sep="");
cat("Output File: ", OutputRoot, "\n", sep="");
cat("\n");
cat("Treatment Variable: ", TreatmentVar, "\n");
cat("Time Variable: ", TimeVar, "\n");
cat("Subject Variable: ", SubjectVar, "\n");
cat("\n");
cat("Reset Time Points? ", ResetTimePoints, "\n");

options(width=120);
cat("Text Line Width: ", options()$width, "\n", sep="");

##############################################################################

pdf(paste(OutputRoot, ".alr.time_series.pdf", sep=""), height=8.5, width=11);

##############################################################################

load_factors=function(fname){
	factors=data.frame(read.table(fname,  header=TRUE, row.names=1, check.names=FALSE));
	factor_names=colnames(factors);

	ignore_idx=grep("^IGNORE\\.", factor_names);

	if(length(ignore_idx)!=0){
		return(factors[-ignore_idx]);
	}else{
		return(factors);
	}
}

load_summary_file=function(fname){
	inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", row.names=1))
	counts_mat=inmat[,2:(ncol(inmat))];
	return(counts_mat);
}

load_reference_levels_file=function(fname){
        inmat=as.matrix(read.table(fname, sep="\t", header=F, check.names=FALSE, comment.char="#", row.names=1))
        colnames(inmat)=c("ReferenceLevel");
        print(inmat);
        cat("\n");
        if(ncol(inmat)!=1){
                cat("Error reading in reference level file: ", fname, "\n");
                quit(status=-1);
        }
        return(inmat);
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

##############################################################################

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

additive_log_ratio=function(ordered_matrix){
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

# Reset time points to 0.
reset_time_points=function(factors, subj_var, time_var){
	#print(factors);
	groups_arr=factors[,subj_var];
	groups=unique(groups_arr);
	for(grp in groups){
		#cat("resetting: ", grp, "\n");
		times=factors[groups_arr==grp,time_var];
		min_time=min(times);
		factors[groups_arr==grp,time_var]=times-min_time;
	}
	return(factors);	
}

##############################################################################

plot_text=function(strings, plot_only=F){

	if(!plot_only){
		for(line in strings){
			cat(line, "\n", sep="");
		}
		return;
	}

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

plot_heatmap=function(mat, title="", noPrintZeros=F, guideLines=F){

        if(is.null(dim(mat))){
                cat(title, " Matrix is NULL.  No heatmap generated.\n");
                return;
        }

        cat("Plotting: ", title, "\n");

	orig_family=par()$family;
        par(family="Courier");
        par(oma=c(10, 10, 1, .5));
        par(mar=c(5.1, 4.1, .5, .5));

        # Generate colors from red to blue
        colors=(rainbow(2^16, start=0, end=0.65));

        # Remember that rows and columsn are reversed in the image
        image(1:nrow(mat),1:ncol(mat), mat,
                xaxt="n", yaxt="n",
                xlab="", ylab="",
                col=colors
        );

        # Pad strings
        cnames=paste(colnames(mat), " ", sep="");
        rnames=paste(rownames(mat), " ", sep="");

        # Get longest length of each column or row label
        cname_max_len=max(nchar(cnames));
        rname_max_len=max(nchar(rnames));

        # Get the number of rows and columns
        ncols=ncol(mat);
        nrows=nrow(mat);

        cscale=min(c(25/cname_max_len, 25/ncols));
        rscale=min(c(25/rname_max_len, 25/nrows));

        max_width=max(nchar(sprintf("%.2f",mat)));
        cell_cex=(3.5/max_width)*sqrt(min(c(cscale, rscale))^2);

        for(i in 1:nrow(mat)){
                for(j in 1:ncol(mat)){

                        if(!is.na(mat[i,j]) && (noPrintZeros && mat[i,j]==0)){
                                # Skip
                        }else{
                                str=sprintf("%.2f",mat[i,j]);
                                str=gsub("^0\\.",".", str);
                                text(i,j,labels=str, cex=cell_cex, srt=45);
                        }
                }
        }
        # Plot guidelines
        if(guideLines){

                splits=c(2,3,4,5);

                h_remainder=ncols %% splits;
                best_h_split=splits[(max(which(h_remainder==0)))];
                if(ncols>best_h_split){
                        h_line_pos=seq(best_h_split, ncols, best_h_split)+.5;
                        abline(h=h_line_pos, col="black", lty="dashed");
                        abline(h=h_line_pos, col="white", lty="dotted");
                }

                v_remainder=nrows %% splits;
                best_v_split=splits[(max(which(v_remainder==0)))];
                if(nrows>best_v_split){
                        v_line_pos=seq(best_v_split, nrows, best_v_split)+.5;
                        abline(v=v_line_pos, col="black", lty="dashed");
                        abline(v=v_line_pos, col="white", lty="dotted");
                }

        }

        # Plot the labels
        mtext(cnames, at=1:ncols, side=2, las=2, cex=cscale);
        mtext(rnames, at=1:nrows, side=1, las=2, cex=rscale);

        # Plot the title
        mtext(title, line=0, at=nrows*.5, side=3, font=2);

	par(family=orig_family);

}

##############################################################################

# Get color assignments
get_colors=function(num_col, alpha=1){
        colors=hsv(seq(0,1,length.out=num_col), c(1,.5), c(1,.75,.5), alpha);
        color_mat_dim=ceiling(sqrt(num_col));
        color_pad=rep("grey", color_mat_dim^2);
        color_pad[1:num_col]=colors;
        color_mat=matrix(color_pad, nrow=color_mat_dim, ncol=color_mat_dim);
        colors=as.vector(t(color_mat));
        colors=colors[colors!="grey"];
}

##############################################################################
##############################################################################

# Load matrix
counts=load_summary_file(SummaryFile);
num_taxa=ncol(counts);
num_samples=nrow(counts);

cat("Num Taxa Loaded: ", num_taxa, "\n");
cat("Num Samples Loaded: ", num_samples, "\n");

##############################################################################

# Load factors
factors=load_factors(FactorsFile);
factor_names=colnames(factors);
num_factors=ncol(factors);
factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);

if(ResetTimePoints){
	factors=reset_time_points(factors, SubjectVar, TimeVar);
}

cat("\n");
cat(num_factors, " Factor(s) Loaded:\n", sep="");
print(factor_names);
cat("\n");

##############################################################################

# Reconcile factors with samples
factor_sample_ids=rownames(factors);
counts_sample_ids=rownames(counts);

#print(factor_sample_id);
#print(counts_sample_id);

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

counts=counts[shared_sample_ids,,drop=F];
factors=factors[shared_sample_ids,,drop=F];

print(factors);
##############################################################################

# Normalize
normalized=normalize(counts);

# Assign 0's to values smaller than smallest abundance across entire dataset
min_assay=min(normalized[normalized!=0]);
cat("Lowest non-zero value: ", min_assay, "\n\n", sep="");
zero_replacment=min_assay/10;
normalized[normalized==0]=zero_replacment;

##############################################################################

# Reorder categories by abundance
mean_abund=apply(normalized, 2, mean);
ix=order(mean_abund, decreasing=TRUE);
normalized=normalized[,ix];
mean_abund=mean_abund[ix];
sorted_taxa_names=colnames(normalized);

num_top_taxa=NumVariables;
num_top_taxa=min(c(num_top_taxa, num_taxa));
prop_abundance_represented=sum(mean_abund[1:num_top_taxa]);


text_arr=c(
	paste("The top ", num_top_taxa, " categories accounting for ", 
		signif(prop_abundance_represented*100, 4), "% of categories are:", sep=""),
	"\n"
);

for(i in 1:num_top_taxa){
	text_arr=c(text_arr, paste("\t", sorted_taxa_names[i], "\t[", signif(mean_abund[i]*100, 4), "%]", sep=""));
}

plot_text(text_arr);

##############################################################################

top_perc_categ=function(mean_abund, abn_acct=.1){
	ord_ix=order(mean_abund, decreasing=T);
	sorted_abd=mean_abund[ord_ix];
	cum_sum=cumsum(sorted_abd);
	keep=(which(cum_sum<abn_acct))+1;
	if(length(keep)){
		keep_max=max(keep);
	}else{
		keep_max=1;
	}
	keep_names=(names(sorted_abd[1:keep_max]));	
	return(keep_names[!is.na(keep_names)]);
}

bottom_perc_categ=function(mean_abund, abn_acct=.1){
	ord_ix=order(mean_abund, decreasing=F);
	sorted_abd=mean_abund[ord_ix];
	cum_sum=cumsum(sorted_abd);
	keep=(which(cum_sum<abn_acct))+1;
	if(length(keep)){
		keep_max=max(keep);
	}else{
		keep_max=1;
	}
	keep_names=(names(sorted_abd[1:keep_max]));	
	return(keep_names[!is.na(keep_names)]);
}

surround_perc_categ=function(mean_abund, abn_acct=.1, cat_index){
	ord_ix=order(mean_abund, decreasing=T);
	sorted_abd=mean_abund[ord_ix];

	cat_names=names(sorted_abd);
	center_ix=which(cat_index==cat_names);
	num_cat=length(mean_abund);

	lt_cat=c();
	gt_cat=c();

	if(center_ix<num_cat){
		lt_cat=top_perc_categ(sorted_abd[(center_ix+1):num_cat], abn_acct/2);
	}
	if(center_ix>1){
		gt_cat=bottom_perc_categ(sorted_abd[1:(center_ix-1)], abn_acct/2);
	}

	return(c(gt_cat, lt_cat));
}

##############################################################################

# Extract top taxa for ALR
if(num_top_taxa>= num_taxa){
	num_top_taxa = (num_taxa-1);
	cat("Number of taxa to work on was changed to: ", num_top_taxa, "\n");
}

cat("Extracting: ", num_top_taxa, " + 1 (remaining) categories.\n", sep="");

responses=extract_top_categories(normalized, num_top_taxa);
resp_alr_struct=additive_log_ratio(responses);
transformed=resp_alr_struct$transformed;

top_cat_ix=top_perc_categ(mean_abund, abn_acct=.25);
bot_cat_ix=bottom_perc_categ(mean_abund, abn_acct=.25);

text_arr=c(
	"\n",
	paste("Top Categories: (", signif(sum(mean_abund[top_cat_ix])*100, 4), "% of Avg Abn)\n", sep=""),
	capture.output(print(top_cat_ix)),
	"\n"
);
plot_text(text_arr);

text_arr=c(
	paste("Bottom Categories: (", signif(sum(mean_abund[bot_cat_ix])*100, 4), "% of Avg Abn)\n", sep=""),
	capture.output(print(bot_cat_ix)),
	"\n"
);
plot_text(text_arr);

##############################################################################
# Get treatment, time, and subject groups
treatment_levels=levels(as.factor(factors[,TreatmentVar]));
subject_levels=levels(factors[,SubjectVar]);
time_ranges=range(factors[,TimeVar]);

num_treatment_levels=length(treatment_levels);
num_subjects=length(subject_levels);

cat("Num Subjects: ", num_subjects, "\n");
cat("Num Treatments Levels: ", num_treatment_levels, "\n");
cat("Time Levels: ", time_ranges[1], "-", time_ranges[2], "\n");

##############################################################################

transform_names=c("Logit", "ALR Rem", "ALR Top", "ALR Sur", "ALR Bot");
num_transforms=length(transform_names);

##############################################################################

indiv_colors=get_colors(num_subjects);
names(indiv_colors)=subject_levels;
palette(indiv_colors);

trt_colors_transp=get_colors(num_treatment_levels+3, 0.3);
trt_colors_opaque=get_colors(num_treatment_levels+3, 1.0);


#num_top_taxa=1;
for(var_ix in 1:num_top_taxa){

	cat("##########################################################################\n");
	cat("#                                                                        #\n");
	cat("# ", sorted_taxa_names[var_ix], "\n");
	cat("#                                                                        #\n");
	cat("##########################################################################\n");

	par(oma=c(.8,.5,3,3));
	par(mfrow=c(8,5));
	first_plot=T;

	alr_range=range(transformed[,var_ix]);
	abd_range=range(normalized[,var_ix]);
	logit_abd_range=log10(abd_range/(1-abd_range));

	# Get the categories surrounding current target category
	sur_cat_ix=surround_perc_categ(mean_abund, abn_acct=.25, sorted_taxa_names[var_ix]);


	# Compute and organize numbers
	trt_list=list();
	trt_ends=list();
	all_transform_range=numeric();
	for(trt in treatment_levels){
		cat("Working on Treatment: ", trt, "\n");
		trt_samp=(factors[,TreatmentVar]==trt);
		sbj_in_trt=unique(factors[trt_samp,SubjectVar]);
		num_sbj_in_trt=length(sbj_in_trt);
		
		sbj_list=list();

		starts_mat=matrix(0, nrow=num_sbj_in_trt, ncol=num_transforms, dimnames=list(sbj_in_trt, transform_names));
		ends_mat=matrix(0, nrow=num_sbj_in_trt, ncol=num_transforms, dimnames=list(sbj_in_trt, transform_names));

		for(sbj in sbj_in_trt){

			cat("Working on Subject: ", sbj, "\n");

			samps=rownames(factors[factors[,SubjectVar]==sbj,]);
			num_samps=length(samps);
			time_pts=factors[samps,TimeVar];
			ord_ix=order(time_pts);

			samps=samps[ord_ix];
			time_pts=time_pts[ord_ix];

			norm_abd=normalized[samps,var_ix];
			logit_abd=log10(norm_abd/(1-norm_abd));
			alr=transformed[samps,var_ix];
			
			top_den=apply(normalized[samps, top_cat_ix, drop=F], 1, sum);
			sur_den=apply(normalized[samps, sur_cat_ix, drop=F], 1, sum);
			bot_den=apply(normalized[samps, bot_cat_ix, drop=F], 1, sum);

			alr_top=log10(norm_abd/top_den);
			alr_sur=log10(norm_abd/sur_den);
			alr_bot=log10(norm_abd/bot_den);

			sbj_list[[sbj]]=matrix(c(time_pts, logit_abd, alr, alr_top, alr_sur, alr_bot),
						ncol=6, byrow=F, 
						dimnames=list(
							c(samps),
							c("Time", transform_names)
						)
					);
			all_transform_range=range(c(all_transform_range, sbj_list[[sbj]][,2:6]));
			
			for(trans_id in transform_names){
				starts_mat[sbj, trans_id]=sbj_list[[sbj]][1, trans_id];
				ends_mat[sbj, trans_id]=sbj_list[[sbj]][length(time_pts), trans_id];
			}

			
		}

		ends_list=list();
		ends_list[["starts"]]=starts_mat;
		ends_list[["ends"]]=ends_mat;

		trt_list[[trt]]=sbj_list;
		trt_ends[[trt]]=ends_list;
	}

	#print(trt_list);
	#print(trt_ends);

	#------------------------------------------------------------------------------------------------
	# Plot individual time series
	for(trt in treatment_levels){

		trt_samp=(factors[,TreatmentVar]==trt);
		sbj_in_trt=unique(factors[trt_samp,SubjectVar]);
		
		for(sbj in sbj_in_trt){

			# Determine whether to label each column of plots
			if(first_plot){
				first_plot=F;
				title=c("Logit", "ALR Remain", "ALR Top", "ALR Surr", "ALR Bottom");
				top_mar=1;
			}else{
				title=c("","","","","");
				top_mar=0;
			}

			time_pts=trt_list[[trt]][[sbj]][,"Time"];
			logit_abd   =trt_list[[trt]][[sbj]][,"Logit"];
			alr     =trt_list[[trt]][[sbj]][,"ALR Rem"];
			alr_top =trt_list[[trt]][[sbj]][,"ALR Top"];
			alr_sur =trt_list[[trt]][[sbj]][,"ALR Sur"];
			alr_bot =trt_list[[trt]][[sbj]][,"ALR Bot"];

			indi_col=indiv_colors[sbj];
			num_pts=length(time_pts);

			# Compute position of axis labels and guide lines
			axis_ticks=seq(all_transform_range[1], all_transform_range[2], length.out=7);
			axis_ticks=axis_ticks[c(-1,-7)];
			guide_line_col="grey90";
			guide_line_lwd=.3;

			plot_mini=function(x, y, xlim, ylim, main){
				plot(x, y,
					xlim=xlim, ylim=ylim,
					xlab="", ylab="",
					xaxt="n", yaxt="n",
					main=main, type="n"
					);
				last=length(x);
				abline(h=axis_ticks, col=guide_line_col, lwd=guide_line_lwd);
				points(x, y, col=indi_col, type="l", lwd=2);
				points(x, y, col="black", type="b", pch=16, cex=.25, lwd=.25);
				points(x[c(1,1,last)], y[c(1,1,last)], col=indi_col, pch=c(17,1,15), cex=c(1,2,1.25));
			}
			
			par(mar=c(0, 1, top_mar, 0));
			plot_mini(time_pts, logit_abd, time_ranges, all_transform_range, title[1]);
			title(ylab=sbj, line=.05, cex.lab=1);	

			plot_mini(time_pts, alr, time_ranges, all_transform_range, title[2]);
			plot_mini(time_pts, alr_top, time_ranges, all_transform_range, title[3]);
			plot_mini(time_pts, alr_sur, time_ranges, all_transform_range, title[4]);
			plot_mini(time_pts, alr_bot, time_ranges, all_transform_range, title[5]);

			axis(side=4, at=axis_ticks, labels=signif(axis_ticks,2),
				las=2, cex.axis=.75);


			# Label top margin with category name being analyzed
			if(par()$page){
				mtext(paste(
					var_ix, ".) ", sorted_taxa_names[var_ix], " (", signif(mean_abund[var_ix]*100,3), "%) ",
					sep=""),
				side=3, line=1.5, outer=T);
			}
		}
	}

	#------------------------------------------------------------------------------------------------
	# Plot individuals combined per treatment

	par(mfrow=c(num_treatment_levels, 5));
	num_treatment_levels=length(treatment_levels);

	for(trt_ix in 1:num_treatment_levels){
		trt=treatment_levels[trt_ix];

		trt_samp=(factors[,TreatmentVar]==trt);
		sbj_in_trt=unique(factors[trt_samp,SubjectVar]);

		for(i in 1:num_transforms){

			if(trt_ix==1){
				transform_label=transform_names[i];
				top_mar=1;	
			}else{
				transform_label="";
				top_mar=0;
			}	

			par(mar=c(2, 1, top_mar, 0));
			plot(0, xlim=time_ranges, ylim=all_transform_range,
				type="n", xlab="Time", ylab="",
				xaxt="n", yaxt="n"
			);
			title(main=transform_label);

			if(i==1){
				mtext(text=trt, side=2, cex=.75);
			}

			if(i==num_transforms){
				axis(side=4, at=axis_ticks, labels=signif(axis_ticks,2),
					las=2, cex.axis=.75);
			}

			axis_ticks=seq(all_transform_range[1], all_transform_range[2], length.out=7);
			axis_ticks=axis_ticks[c(-1,-7)];
			guide_line_col="grey90";
			guide_line_lwd=.3;
			abline(h=axis_ticks, col=guide_line_col, lwd=guide_line_lwd);

			for(sbj in sbj_in_trt){
				time_pts =trt_list[[trt]][[sbj]][,"Time"];
				transform=trt_list[[trt]][[sbj]][,transform_names[i]];
				last=length(time_pts);
				points(time_pts, transform, col=indiv_colors[sbj], type="l", lwd=2);
				points(time_pts, transform, col="black", type="b", pch=16, cex=.25, lwd=.25);
				points(time_pts[c(1,1,last)], transform[c(1,1,last)], 
					col=indiv_colors[sbj], pch=c(17,1,15), cex=c(1,2,1.25));
			}


			if(trt_ix==num_treatment_levels){
				time_ticks=round(seq(time_ranges[1], time_ranges[2], length.out=10));
				axis(side=1, line=0,
					at=time_ticks,
					labels=time_ticks,
					las=2
				);
			}

		}
	}
	mtext(paste(
		var_ix, ".) ", sorted_taxa_names[var_ix], " (", signif(mean_abund[var_ix]*100,3), "%) ",
		sep=""),
	side=3, line=1.5, outer=T);

	#------------------------------------------------------------------------------------------------
	# Plot before/after scatter plot

	
	par(oma=c(0,0,3,0));
	par(mar=c(1, 1, 1, 1.1));
	par(mfrow=c(num_treatment_levels+1, num_transforms));

	# Plot individual
	for(trt_ix in 1:num_treatment_levels){
		trt=treatment_levels[trt_ix];

		trt_samp=(factors[,TreatmentVar]==trt);
		sbj_in_trt=unique(factors[trt_samp,SubjectVar]);

		for(i in 1:num_transforms){

			left_label=ifelse(i==1, "after", "");
			top_label=ifelse(trt_ix==1, transform_names[i], "");

			starts_val=trt_ends[[trt_ix]]$starts[,i];		
			ends_val=trt_ends[[trt_ix]]$ends[,i];		

			plot(0,0, xlim=all_transform_range, ylim=all_transform_range, 
				col=trt_ix, type="n", xaxt="n", yaxt="n"
			);
			title(ylab=left_label, line=0, cex=.8);
			title(main=top_label, line=.2, cex=.8);

			if(i==num_transforms){
				mtext(text=trt, side=4, line=.1);
			}
			abline(0, 1, col="grey85", lwd=2);

			#points_mat=cbind(starts_val, ends_val);
			#ch_ix=chull(points_mat);
			#polygon(points_mat[ch_ix,], border=trt_ix, col=trt_colors_transparent);	

			abline(v=mean(starts_val), col=trt_colors_transp[trt_ix]);
			abline(h=mean(ends_val), col=trt_colors_transp[trt_ix]);
			points(starts_val, ends_val, col=trt_colors_transp[trt_ix], cex=3, pch=16);
			points(starts_val, ends_val, col=trt_colors_opaque[trt_ix], cex=.25, pch=16);

		}
	}

	# Plot Combined
	for(i in 1:num_transforms){
		plot(0,0, xlim=all_transform_range, ylim=all_transform_range, 
			type="n", xaxt="n", yaxt="n",
			);

		left_label=ifelse(i==1, "after", "");
		title(ylab=left_label, line=0);
		title(xlab="before", line=0);
		if(i==num_transforms){
			mtext(text="Combined", side=4, line=.1, cex=.8);
		}

		abline(0, 1, col="grey85", lwd=2);

		for(trt_ix in 1:num_treatment_levels){
			trt=treatment_levels[trt_ix];

			trt_samp=(factors[,TreatmentVar]==trt);
			sbj_in_trt=unique(factors[trt_samp,SubjectVar]);

			if(trt_ix==1){
				transform_label=transform_names[i];
			}else{
				transform_label="";
			}	

			starts_val=trt_ends[[trt_ix]]$starts[,i];		
			ends_val=trt_ends[[trt_ix]]$ends[,i];		

			abline(v=mean(starts_val), col=trt_colors_transp[trt_ix]);
			abline(h=mean(ends_val), col=trt_colors_transp[trt_ix]);

			points(starts_val, ends_val, col=trt_colors_transp[trt_ix], cex=3, pch=16);
			points(starts_val, ends_val, col=trt_colors_opaque[trt_ix], cex=.25, pch=16);

		}
	}

	# Label plot
	mtext(paste(
		var_ix, ".) ", sorted_taxa_names[var_ix], " (", signif(mean_abund[var_ix]*100,3), "%) ",
		sep=""),
		side=3, line=1.5, outer=T);

	

}

#############################################################################	
dev.off();

cat("Done.\n");
print(warnings());
q(status=0);
