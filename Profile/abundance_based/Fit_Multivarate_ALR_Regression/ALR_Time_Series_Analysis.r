#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
library(car);
library(boot);
options(useFancyQuotes=F);

params=c(
	"summary_file", "s", 1, "character",
	"factors", "f", 1, "character",
	"num_variables", "p", 2, "numeric",
	"treatment", "t", 2, "character",
	"time", "m", 2, "character",
	"subject", "u", 2, "character",
	"outputroot", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table>\n",
	"	-f <factors>\n",
	"	[-p <number of variables>]\n",
	"\n",
	"	-t <Treatment variable>\n",
	"	-m <tiMe variable>\n",
	"	-u <sUbject variable>\n",
	"\n",
	"	[-o <output filename root>]\n",
	"\n",
	"\n");

if(!length(opt$summary_file) || !length(opt$factors)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputRoot=gsub(".summary_table.xls", "", opt$summary_file);
}else{
	OutputRoot=opt$outputroot;
}

if(!length(opt$num_variables)){
	NumVariables=20;
}else{
	NumVariables=opt$num_variables;
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

options(width=120);
cat("Text Line Width: ", options()$width, "\n", sep="");

rnd=paste(".", sprintf("%i",sample(1000, 1)), sep="");

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

# Load matrix
counts=load_summary_file(SummaryFile);
num_taxa=ncol(counts);
num_samples=nrow(counts);

# Shorten names:
full_names=colnames(counts);
splits=strsplit(full_names, " ");
short_names=character();
for(i in 1:length(full_names)){
        short_names[i]=tail(splits[[i]], 1);
}
colnames(counts)=short_names;

# Normalize
normalized=normalize(counts);

# Reorder by abundance
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

# Load factors
factors=load_factors(FactorsFile);
factor_names=colnames(factors);
num_factors=ncol(factors);
factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);

cat("\n");
cat(num_factors, " Factor(s) Loaded:\n", sep="");
print(factor_names);
cat("\n");

#print(factors);
#summary(factors);

continuous_factors=factors;
is_continous_factor=logical(num_factors);

for(f in 1:num_factors){
	level_info=levels(factors[,f]);
	is_continous_factor[f]=is.null(level_info);

	if(is_continous_factor[f]){
		# do nothing
	}else if(length(level_info)==2){
		# Convert two level factors to numeric
		is_continous_factor[f]=TRUE;
		continuous_factors[,f]=as.integer(continuous_factors[,f]);
	}else{
		is_continous_factor[f]=FALSE;
	}
}

continuous_factors=continuous_factors[,is_continous_factor, drop=F];
print(continuous_factors);

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

# Reorder data by sample id
normalized=normalized[shared_sample_ids,];
num_samples=nrow(normalized);
factors=factors[shared_sample_ids,,drop=F];

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

additive_log_rato=function(ordered_matrix){
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

##############################################################################

# Assign 0's to values smaller than smallest abundance across entire dataset
min_assay=min(normalized[normalized!=0]);
cat("Lowest non-zero value: ", min_assay, "\n\n", sep="");
zero_replacment=min_assay/10;
normalized[normalized==0]=zero_replacment;

##############################################################################

if(num_top_taxa>= num_taxa){
	num_top_taxa = (num_taxa-1);
	cat("Number of taxa to work on was changed to: ", num_top_taxa, "\n");
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

pdf(paste(OutputRoot, rnd, ".alr.time_series.pdf", sep=""), height=8.5, width=11);

##############################################################################

cat("Performing regression.\n");
cat("Extracting: ", num_top_taxa, " + 1 (remaining) categories.\n", sep="");

responses=extract_top_categories(normalized, num_top_taxa);
resp_alr_struct=additive_log_rato(responses);
transformed=resp_alr_struct$transformed;

#print(factors);
treatment_levels=levels(factors[,TreatmentVar]);
time_levels=levels(factors[,TimeVar]);
subject_levels=levels(factors[,SubjectVar]);

num_treatments=length(treatment_levels);
num_times=length(time_levels);
num_subjects=length(subject_levels);

alr_range=range(transformed);


col_layout=c(
	1,
	1,
	2,
	3
);

mat_layout=col_layout;
for(i in 1:(num_treatments-1)){
	mat_layout=cbind(mat_layout, (col_layout+i*max(col_layout)));
}
print(mat_layout);

colors=c(
	"black",
	"red",
	"green",
	"blue",
	"orange",
	"purple",
	"brown",
	"grey",
	"pink",
	"orange4"
);

palette(colors);

sd_diff_fun=function(data, index){
	sub_data=data[index,];
	grps=(sub_data[,"grp"]==0);
	grpA=sub_data[grps, "values"];
	grpB=sub_data[!grps, "values"];
	stdA=sd(grpA);
	stdB=sd(grpB);
	return(stdA-stdB);
}

# Generate names for comparison matrices
mat_side_dim=(num_times-1)*(num_treatments);
changes_compare_names=character(mat_side_dim);
i=1;
for(trt_ix in 1:num_treatments){
	for(tim_ix in 1:(num_times-1)){
		changes_compare_names[i]=paste(treatment_levels[trt_ix], ":", 
			time_levels[tim_ix], "-",
			time_levels[tim_ix+1], 
			sep="");
		i=i+1;
	}
}

# Generates names for changes comparison matrices
all_comparison_dim=num_times*num_treatments;
all_compar_names=character(all_comparison_dim);
i=1;
for(trt_ix in 1:num_treatments){
	for(tim_ix in 1:num_times){
		all_compar_names[i]=
			paste(treatment_levels[trt_ix], ":", time_levels[tim_ix], sep="");
		i=i+1;
	}
}	

for(var_ix in 1:num_top_taxa){

	grouped_values=array(NA, c(num_subjects, num_times, num_treatments));
	dimnames(grouped_values)=list(subject_levels, time_levels, treatment_levels);

	cat("##########################################################################\n");
	cat("#                                                                        #\n");
	cat("# ", sorted_taxa_names[var_ix], "\n");
	cat("#                                                                        #\n");
	cat("##########################################################################\n");

	alr=transformed[,var_ix];
	
	for(samp_ix in 1:num_samples){
		subj_id=factors[samp_ix, SubjectVar];
		time_id=factors[samp_ix, TimeVar];
		treat_id=factors[samp_ix, TreatmentVar];
		grouped_values[subj_id, time_id, treat_id]=alr[samp_ix];
	}

	# Compute changes between time periods within the same donor.
	changes_matrix=array(NA, c(num_subjects, num_times-1, num_treatments));
	dimnames(changes_matrix)=list(subject_levels, tail(time_levels, num_times-1), treatment_levels);
	for(trt_ix in 1:num_treatments){
		for(subj_ix in 1:num_subjects){
			for(time_ix in 1:(num_times-1)){
				changes_matrix[subj_ix, time_ix, trt_ix]=
					abs(
						grouped_values[subj_ix, time_ix+1, trt_ix]-
						grouped_values[subj_ix, time_ix, trt_ix]
					);
			}
		}
	}
	
	# Create names for changes matrix
	changes_names=character();
	for(tim_ix in 1:(num_times-1)){
		changes_names[tim_ix]=paste(time_levels[tim_ix], "-", time_levels[tim_ix+1], sep="");
	}
	
	# Find max change
	all_changes=as.vector(changes_matrix);
	all_changes=all_changes[!is.na(all_changes)];
	max_change=max(all_changes);
	
	# Generate Speghetti Plots for each treatment
	par(mfrow=c(1, num_treatments));
	par(oma=c(0,0,4,0));
	par(mar=c(5.1, 4.1, 1, 1));

	layout(mat_layout);

	for(trt_ix in 1:num_treatments){

		# Set up empty plot
		plot(1:num_times, time_levels, type="n", 
			xlim=c(1-.5, num_times+.5), ylim=alr_range,
			xlab="Time",
			ylab="ALR Transformed Abundances",
			main=treatment_levels[trt_ix],
			xaxt="n",
			yaxt="n"
		);
		axis(side=1, at=1:num_times, labels=time_levels);
		tickmark=seq(floor(alr_range[1]), ceiling(alr_range[2]), 1);
		axis(side=2, at=tickmark, labels=sprintf("%3.0f", tickmark), las=2, cex=.8);

		# Generate Spaghetti Plot each subject
		plotted=0;
		for(subj_ix in 1:num_subjects){

			if(all(is.na(grouped_values[subj_ix, time_levels, trt_ix]))){
				next;
			}else{
				plotted=plotted+1;

				points(1:num_times, grouped_values[subj_ix, time_levels, trt_ix], type="b",
					col=plotted);

				subj_name=subject_levels[subj_ix];

				# Label left-most point with donor ID
				text(1, grouped_values[subj_ix, 1, trt_ix], 
					pos=2,
					labels=subj_name, col=plotted, cex=.5);

				# Label right-most point with donor ID
				text(num_times, grouped_values[subj_ix, num_times, trt_ix], 
					pos=4,
					labels=subj_name, col=plotted, cex=.5);
			}
		}

		# Compute standard deviation at each point
		sds=apply(grouped_values[, , trt_ix], 2, function(x){sd(x, na.rm=T)});
		barplot(sds, xlab="Standard Deviations", ylim=c(0,5));

		# Plot changes boxplot
		changes_list=list();
		for(tim_ix in 1:(num_times-1)){
			changes_list[[tim_ix]]=changes_matrix[,tim_ix,trt_ix];
		}
		boxplot(changes_list, xlim=c(0, num_times), cex=0,
			xaxt="n",	
			ylim=c(0, max_change), xlab="Changes");
		stripchart(changes_list, method="jitter", add=T, vertical=T, pch=1, col="blue");

		axis(side=1, at=1:(num_times-1), labels=changes_names);
		
		

	}
	mtext(sorted_taxa_names[var_ix], outer=T, side=3, cex=1.5, line=1);
	
	#--------------------------------------------------------------------
	# Compute p-values for differences in medians and standard deviation	

	median_pval_mat=matrix(NA, ncol=all_comparison_dim, nrow=all_comparison_dim);
	rownames(median_pval_mat)=all_compar_names;
	colnames(median_pval_mat)=all_compar_names;

	stdev_pval_mat=matrix(NA, ncol=all_comparison_dim, nrow=all_comparison_dim);
	rownames(stdev_pval_mat)=all_compar_names;
	colnames(stdev_pval_mat)=all_compar_names;
	
	pm_ixA=1;
	for(trt_ixA in 1:num_treatments){
		for(time_ixA in 1:num_times){
			pm_ixB=1;
			for(trt_ixB in 1:num_treatments){
				for(time_ixB in 1:num_times){

					# Medians
					res=wilcox.test(
						grouped_values[, time_ixA, trt_ixA],
						grouped_values[, time_ixB, trt_ixB],
						alternative="less"
					);
					median_pval_mat[pm_ixA, pm_ixB]=res$p.value;
						
					# Std. Dev
					grpA=grouped_values[,time_ixA, trt_ixA];
					grpB=grouped_values[,time_ixB, trt_ixB];

					grpA=grpA[!is.na(grpA)];
					grpB=grpB[!is.na(grpB)];
					values=c(grpA, grpB);
					groups=c(rep(0, length(grpA)), rep(1, length(grpB)));
					values_mat=cbind(values, groups);
					colnames(values_mat)=c("values", "grp");

					# Boot strap differences in standard deviation
					b_res=boot(values_mat, sd_diff_fun, R=1000, strata=groups);
					sd_pval=sum(b_res$t>0)/length(b_res$t);
					stdev_pval_mat[pm_ixA, pm_ixB]=sd_pval;

					pm_ixB=pm_ixB+1;
				}
			}
			pm_ixA=pm_ixA+1;
		}
	}

	#--------------------------------------------------------------------
	# Compute p-values for changes between groups

	changes_pval_mat=matrix(NA, ncol=mat_side_dim, nrow=mat_side_dim);
	colnames(changes_pval_mat)=changes_compare_names;
	rownames(changes_pval_mat)=changes_compare_names;

	pm_ixA=1;
	for(trt_ixA in 1:num_treatments){
		for(time_ixA in 1:(num_times-1)){
			pm_ixB=1;
			for(trt_ixB in 1:num_treatments){
				for(time_ixB in 1:(num_times-1)){
					res=wilcox.test(
						changes_matrix[,time_ixA, trt_ixA],			
						changes_matrix[,time_ixB, trt_ixB],
						alternative="less"
					);
					changes_pval_mat[pm_ixA, pm_ixB]=res$p.value;
					pm_ixB=pm_ixB+1;
				}
			}
			pm_ixA=pm_ixA+1;
		}
	}

	#--------------------------------------------------------------------
	# Plot pvalue matrices

	par(mfrow=c(1,1));
	par(oma=c(0,0,4,0));

	plot_heatmap(median_pval_mat, 
		title=paste(sorted_taxa_names[var_ix], ", Alternate Hypothesis: Row > Col Median", sep=""),
		guideLines=T);
	plot_heatmap(stdev_pval_mat, 
		title=paste(sorted_taxa_names[var_ix], ", Alternate Hypothesis: Row > Col Std Dev", sep=""),
		guideLines=T);
	plot_heatmap(changes_pval_mat, 
		title=paste(sorted_taxa_names[var_ix], ", Alternate Hypothesis: Row > Col Changes", sep=""),
		guideLines=T);



}

#############################################################################	
dev.off();

cat("Done.\n");
print(warnings());
q(status=0);
