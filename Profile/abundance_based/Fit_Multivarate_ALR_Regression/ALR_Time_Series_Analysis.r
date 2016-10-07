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

plot_text=function(strings){
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

cat("\nThe top ", num_top_taxa, " taxa are:\n", sep="");
for(i in 1:num_top_taxa){
	cat("\t", sorted_taxa_names[i], "\t[", mean_abund[i], "]\n", sep="");
}
cat("\n");

cat("Accounting for ", prop_abundance_represented, " of taxa.\n", sep="");
cat("\n");

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

cat("\n");
cat("Top Categories: (", sum(mean_abund[top_cat_ix])*100, "% of Avg Abn)\n", sep="");
print(top_cat_ix);
cat("\n");
cat("Bottom Categories: (", sum(mean_abund[bot_cat_ix])*100, "% of Avg Abn)\n", sep="");
print(bot_cat_ix);
cat("\n");

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

pdf(paste(OutputRoot, ".alr.time_series.pdf", sep=""), height=11, width=8.5);

par(oma=c(.8,.5,3,.8));

for(var_ix in 1:num_top_taxa){

	cat("##########################################################################\n");
	cat("#                                                                        #\n");
	cat("# ", sorted_taxa_names[var_ix], "\n");
	cat("#                                                                        #\n");
	cat("##########################################################################\n");

	par(mfrow=c(15,5));
	first_plot=T;

	alr_range=range(transformed[,var_ix]);
	abd_range=range(normalized[,var_ix]);
	logit_abd_range=log10(abd_range/(1-abd_range));

	sur_cat_ix=surround_perc_categ(mean_abund, abn_acct=.25, sorted_taxa_names[var_ix]);

	for(trt in treatment_levels){
		cat("Working on Treatment: ", trt, "\n");
		trt_samp=(factors[,TreatmentVar]==trt);

		sbj_in_trt=unique(factors[trt_samp,SubjectVar]);
		
		for(sbj in sbj_in_trt){

			cat("\nWorking on Subject: ", sbj, "\n");

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

			mat=matrix(c(norm_abd, top_den, sur_den, bot_den), nrow=num_samps, ncol=4, byrow=T);
			colnames(mat)=c("Categ Abund", "Top", "Surrounding", "Bottom");
			rownames(mat)=samps;
			print(mat);

			alr_top=log10(norm_abd/top_den);
			alr_sur=log10(norm_abd/sur_den);
			alr_bot=log10(norm_abd/bot_den);

			if(first_plot){
				first_plot=F;
				title=c("Logit", "ALR Remain", "ALR Top", "ALR Surr", "ALR Bottom");
				top_mar=1;
			}else{
				title=c("","","","","");
				top_mar=0;
			}

			# Plot logit(abundance)
			par(mar=c(0, 1, top_mar, 0));
			plot(time_pts, logit_abd,
				xlim=time_ranges, ylim=logit_abd_range,
				type="b",
				xlab="", ylab="",
				xaxt="n", yaxt="n",
				main=title[1]
				);
			title(ylab=sbj, line=.05, cex.lab=.55);	

			par(mar=c(0, 1, top_mar, 0));

			# Plot alr, with denominator = not top N categories
			plot(time_pts, alr,
				xlim=time_ranges, ylim=alr_range,
				type="b",
				xlab="", ylab="",
				xaxt="n", yaxt="n",
				main=title[2]
				);

			# Plot with top cat as denominator
			plot(time_pts, alr_top,
				xlim=time_ranges, ylim=alr_range,
				type="b",
				xlab="", ylab="",
				xaxt="n", yaxt="n",
				main=title[3]
				);

			# Plot with surrounding cat as denominator
			plot(time_pts, alr_sur,
				xlim=time_ranges, ylim=alr_range,
				type="b",
				xlab="", ylab="",
				xaxt="n", yaxt="n",
				main=title[4]
				);

			# Plot with low cat as denominator
			plot(time_pts, alr_bot,
				xlim=time_ranges, ylim=alr_range,
				type="b",
				xlab="", ylab="",
				xaxt="n", yaxt="n",
				main=title[5]
				);


			if(par()$page){
				mtext(paste(
					var_ix, ".) ", sorted_taxa_names[var_ix], " (", signif(mean_abund[var_ix]*100,3), "%) ",
					sep=""),
				side=3, line=1.5, outer=T);
			}
		}
	}


}

#############################################################################	
dev.off();

cat("Done.\n");
print(warnings());
q(status=0);
