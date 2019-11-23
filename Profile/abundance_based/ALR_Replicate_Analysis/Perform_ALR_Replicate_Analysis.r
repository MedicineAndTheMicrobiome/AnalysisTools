#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
library('vegan');

options(useFancyQuotes=F);
options(width=80);

DEF_NUM_TOP_CAT=25;

source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');

params=c(
	"summary_file", "s", 1, "character",
	"factor_file", "f", 1, "character",
	"sample_coln", "S", 1, "character",
	"replicate_coln", "R", 1, "character",
	"outputroot", "o", 1, "character",

	"num_top_alr_cat", "v", 2, "numeric",
	"contains_remaining", "r", 2, "logical",
	"shorten_category_names", "x", 2, "character"
	
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
	"\n",
	"	[-R (pay attention to 'remaining' category)]\n",
	"	[-v number of ALR categories to include, default=", DEF_NUM_TOP_CAT, ">]\n",
	"	[-x <shorten category names, with separator in double quotes (default=\"\">)]\n",
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


NumTopALRCat=DEF_NUM_TOP_CAT;
UseRemaining=F;
ShortenCategoryNames="";

if(length(opt$shorten_category_names)){
	ShortenCategoryNames=opt$shorten_category_names;
}

if(length(opt$contains_remaining)){
	UseRemaining=T
}

if(length(opt$num_top_alr_cat)){
	NumTopALRCat=opt$num_top_alr_cat;
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
cat("\n");
cat("Num ALR Categories: ", NumTopALRCat, "\n");
cat("Use Remaining?: ", UseRemaining, "\n");
cat("Shorten Category Names Seperator: ", ShortenCategoryNames, "\n");
cat("\n");

pdf(paste(OutputRoot, ".rep_analys.alr.pdf", sep=""), height=14, width=8.5);



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

##############################################################################

plot_text(c(
	paste("Summary File : ", SummaryFile, sep=""),
	paste("Factors File: ", FactorsFile, sep=""),
	paste("Output File: ", OutputRoot, sep=""),
	paste("Replicate ID Column Name: ", Replicate_ColumnName, "\n", sep=""),
	paste("Sample Name Column Name: ", Sample_ColumnName, "\n", sep=""),
	"",
	paste("Num ALR Categories: ", NumTopALRCat, "\n"),
	paste("Use Remaining?: ", UseRemaining, "\n"),
	paste("Shorten Category Names Seperator: ", ShortenCategoryNames, "\n")
));

# Load matrix
counts=load_summary_file(SummaryFile);
num_taxa=ncol(counts);
num_samples=nrow(counts);
#print(counts);

# Clean category names
cat_names=colnames(counts);
cat_names=gsub("\\[", "", cat_names);
cat_names=gsub("\\]", "", cat_names);
cat_names=gsub("\\(", "", cat_names);
cat_names=gsub("\\)", "", cat_names);
colnames(counts)=cat_names;

if(ShortenCategoryNames!=""){
        full_names=colnames(counts);
        splits=strsplit(full_names, ShortenCategoryNames);
        short_names=character();
        for(i in 1:length(full_names)){
                short_names[i]=tail(splits[[i]], 1);
                short_names[i]=gsub("_unclassified$", "_uncl", short_names[i]);
                short_names[i]=gsub("_Incertae_Sedis_", "_InctSeds_", short_names[i]);
        }
        colnames(counts)=short_names;
        cat("Names have been shortened.\n");
}else{
        cat("Keeping original category names...\n");
}

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

unadj_counts=counts;
counts=counts+.5;
normalized=normalize(counts);

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

# Reorder by abundance
cat("Reordering summary table categories by abundance...\n");
mean_abund=apply(normalized, 2, mean);
ix=order(mean_abund, decreasing=TRUE);
normalized=normalized[,ix];
counts=counts[,ix];
unadj_counts=counts[,ix];
mean_abund=mean_abund[ix];

if(UseRemaining){
        normalized=cbind(normalized, normalized_remaining_col_dat);
        mean_abund=c(mean_abund, mean(normalized_remaining_col_dat));
}

sorted_taxa_names=colnames(normalized);

num_top_taxa=NumTopALRCat;
num_top_taxa=min(c(num_top_taxa, num_taxa));
prop_abundance_represented=sum(mean_abund[1:num_top_taxa]);

cat("\nThe top ", num_top_taxa, " taxa are:\n", sep="");
for(i in 1:num_top_taxa){
        cat("\t", sorted_taxa_names[i], "\t[", mean_abund[i], "]\n", sep="");
}
cat("\n");

cat("Accounting for ", prop_abundance_represented, " of taxa.\n", sep="");
cat("\n");

cat_abundances=extract_top_categories(normalized, NumTopALRCat);
alr_struct=additive_log_rato(cat_abundances);
alr_categories_val=alr_struct$transformed;
alr_cat_names=colnames(alr_categories_val);

plot_text(c(
        paste("ALR Categories (Top ", NumTopALRCat, ")", sep=""),
        capture.output(print(alr_cat_names))
));

cat("Categories:\n");
print(alr_cat_names);

##############################################################################
##############################################################################
# For each sample name

replicate_group=as.character(factors[, Replicate_ColumnName]);
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
	bs_prof=matrix(0, nrow=total_bs, ncol=ncol(norm_arr));
	bs_prof_rownames=character(total_bs);

	ncats=ncol(norm_arr);

	for(samp_ix in 1:num_samps){
		row_ix=(1:num_bs)+num_bs*(samp_ix-1);

		bs_cts=t(.5+rmultinom(num_bs, depths_arr[samp_ix], norm_arr[samp_ix,]));
		bs_prof[row_ix,]=bs_cts/(depths_arr[samp_ix]+(.5*ncats));

		bs_prof_rownames[row_ix]=sprintf("%s.%03g", samp_names[samp_ix], 1:num_bs);
		
		samp_ix=samp_ix+1;
	}

	colnames(bs_prof)=colnames(norm_arr);	
	rownames(bs_prof)=bs_prof_rownames;

	return(bs_prof);

} 

plot_smp_vs_grp_barplots=function(bootstrapped_distmat, num_bootstraps, sample_ids, replicate_ids, sample_colors){
	sqr_dm=as.matrix(bootstrapped_distmat);
	num_samps=length(sample_ids);

	medians=numeric(num_samps);
	lb=numeric(num_samps);
	ub=numeric(num_samps);

	for(i in 1:num_samps){
		tar_smp=sample_ids[i];
		tar_bs_ix=(i-1)*num_bootstraps+1:num_bootstraps;
		#print(tar_bs_ix);

		sub_mat=apply(sqr_dm[tar_bs_ix, -tar_bs_ix], 1, mean);

		
		#print(rownames(sub_mat));
		#print(colnames(sub_mat));
	
		qres=quantile(sub_mat, c(0.025, .5, 0.975));
		medians[i]=qres[2];
		lb[i]=qres[1];
		ub[i]=qres[3];
	}	

	
	max_dist=max(ub);
	mids=barplot(medians, names.arg=replicate_ids[sample_ids], col=sample_colors[sample_ids], 
		main="Median and 95%CI Distance to Centroid of Others",
		ylab="Distance",
		ylim=c(0, max_dist*1.2),
		las=2);

	ep=.15;
	for(i in 1:num_samps){
		points(c(mids[i]-ep, mids[i]+ep), rep(ub[i],2), type="l");
		points(c(mids[i]-ep, mids[i]+ep), rep(lb[i],2), type="l");
		points(rep(mids[i],2), c(ub[i], lb[i]), type="l");
	}

	result=list();
	result[["medians"]]=medians;
	result[["repl_ids"]]=as.character(replicate_ids);

	return(result);

}

##############################################################################

# Precompute the densities for all the samples are reference
dens_bw_adj=1;
densities=list();
xranges=list();
ymaxes=list();
alr_categories_val=alr_struct$transformed;
for(alr_ix in alr_cat_names){
	densities[[alr_ix]]=density(alr_categories_val[,alr_ix], adjust=dens_bw_adj);
	xranges[[alr_ix]]=range(alr_categories_val[,alr_ix]);
	ymaxes[[alr_ix]]=max(densities[[alr_ix]]$y);
}

# Variables for accumulate differences across all samples
cum_differences=list();
cum_depths=list();
for(alr_ix in alr_cat_names){
	cum_differences[[alr_ix]]=list();
	for(repix in uniq_repl_names){
		cum_differences[[alr_ix]][[repix]]=numeric();
	}
}
cum_depths[[repix]]=numeric();

# Number of points to bootstrap each replicate
nbs=320;
iter=0;
for(cur_sample_name in uniq_samp_names){

	cat("\n", iter, ": Analyzing Sample Group: ", cur_sample_name, "\n", sep="");
	
	samp_ix=(sample_names==cur_sample_name);
	samp_ids=sample_ids[samp_ix];

	cat("Samples ID in Sample Group:\n");
	print(samp_ids);

	repl_ids=replicate_group[samp_ix];
	names(repl_ids)=samp_ids;
	num_samp=length(samp_ids);
	tot_bs=nbs*num_samp;

	# Abort if no replicate
	if(num_samp<2){
		next;
	}

	# Assign colors
	colors=rainbow(num_samp, start=0, end=.8);
	names(colors)=samp_ids;

	# Profile
	cur_counts=unadj_counts[samp_ids,, drop=F];
	cur_normal=normalize(cur_counts);
	focused_prof=cur_normal[,alr_cat_names];
	remainder=apply(focused_prof, 1, function(x){1-sum(x)});
	cur_prof=cbind(focused_prof, remainder);
	obs_alr=additive_log_rato(cur_prof);
	print(obs_alr);

	# Depth
	read_depths=apply(cur_counts, 1, sum);
	print(read_depths);
	max_depth=max(read_depths);

	# Store rep to depth for later global analysis
	for(sid in samp_ids){
		rep=repl_ids[sid];
		cum_depths[[rep]]=c(cum_depths[[rep]], read_depths[sid]);
	}
	

	# Bootstrap
	bs_profs=bootstrap_profiles(read_depths, cur_prof, nbs);
	bs_alr_struct=additive_log_rato(bs_profs);
	bs_alr_values=bs_alr_struct$transformed;

	#-----------------------------------------------------------------------------
	# Read Depth Bar Plot

	par(oma=c(1,1,5,1));
        #par(mfrow=c(4,2));
	layout(
		matrix(c(
			1,1,
			2,3,
			4,5,
			6,7), byrow=T, nrow=4));	

        par(mar=c(5,20,5,1));
        barplot(read_depths, horiz=T, las=1, col=colors, xlab="Reads/Sample", main="Repl Seq Depth");
	#plot(0, type="n", xlab="", ylab="", bty="n", xaxt="n", yaxt="n", main="");


	#-----------------------------------------------------------------------------
	# Plot replicates in context

        par(mar=c(3,3,4,1));
	for(alr_ix in alr_cat_names){

		if(par()$page){
			mtext(cur_sample_name, outer=T, font=2, cex=1.5, line=1.5);
			par(mfrow=c(4,2));
		}

		# Calculate densities and find max
		max_density=0;
		i=1;
		dens_list=list();
		for(samp_ix in samp_ids){
			alrval=bs_alr_values[(1:nbs)+(i-1)*nbs, alr_ix];
			dens_list[[samp_ix]]=density(alrval, adjust=dens_bw_adj);
			max_density=max(c(max_density,dens_list[[samp_ix]]$y));
			i=i+1;
		}

		# Start empty plot
		maxy=max(max_density, ymaxes[[alr_ix]]);
		plot(0, type="n", xlim=xranges[[alr_ix]], ylim=c(0, maxy),
			xlab="ALR", ylab="Probability",
			main=alr_ix);

		# Draw reference 0 line
		abline(h=0, col="black", lty=2);	

		# Plot reference density
		points(densities[[alr_ix]]$x, densities[[alr_ix]]$y, type="l",
			lwd=3, col="grey");

		# Plot densities for each bs sample
		for(samp_ix in samp_ids){
			dens=dens_list[[samp_ix]];
			points(dens$x, dens$y, type="l", col=colors[samp_ix]);
		}

		#-----------------------------------------------------------------------------

		# Calculate difference from others
		diff_matrix=matrix(NA, nrow=num_samp, ncol=3);
		rownames(diff_matrix)=samp_ids;
		colnames(diff_matrix)=c("lb95", "median", "ub95");

		i=1;
		for(samp_ix in samp_ids){
			samp_bs_ix=(1:nbs)+(i-1)*nbs;

			samp_bs_alr=bs_alr_values[samp_bs_ix, alr_ix];
			other_bs_alr=bs_alr_values[-samp_bs_ix, alr_ix];

			num_tar_samp=length(samp_bs_alr);
			smp_tar_alr=sample(samp_bs_alr, nbs);
			smp_oth_alr=sample(other_bs_alr, nbs, replace=T);
			smp_diff=(smp_tar_alr-smp_oth_alr)
			perc=quantile(smp_diff, c(0.025, .5, .975));	
			diff_matrix[samp_ix, ]=perc;

			i=i+1;
		}

		# calculate ranges of differences
		diff_min=min(diff_matrix[,1]);
		diff_max=max(diff_matrix[,3]);
		diff_range=diff_max-diff_min;

		# Plot differences in barplot	
		par(mar=c(5,3,3,1));
		mids=barplot(diff_matrix[,2], col=colors, 
			names.arg=repl_ids[samp_ids],
			main=alr_ix,
			ylim=c(diff_min-diff_range/8, diff_max+diff_range/8));

		# Draw with 95%CI error bars
		pad=1/8;
		names(mids)=samp_ids;
		for(samp_ix in samp_ids){
			# lb bar
			points(mids[samp_ix]+c(-pad,+pad), rep(diff_matrix[samp_ix, 1],2),
				type="l");
			# ub bar
			points(mids[samp_ix]+c(-pad,+pad), rep(diff_matrix[samp_ix, 3],2),
				type="l");
			# connector
			points(rep(mids[samp_ix],2), diff_matrix[samp_ix, c(1,3)],
				type="l");
				
		}

		#-----------------------------------------------------------------------------

		# Store differences cum_differences[[alr_ix]][[repix]]=numeric();
		for(sid in samp_ids){
			rep=repl_ids[sid];
			cum_differences[[alr_ix]][[rep]]=
				c(cum_differences[[alr_ix]][[rep]], diff_matrix[sid, 2]);
		}
	}
	mtext(cur_sample_name, outer=T, font=2, cex=1.5, line=1.5);
	
	if(iter==20){
		#break;
	}
	iter=iter+1;

}

#print(cum_depths);
#print(cum_differences);

# Plot biases in ALR by for each category by run

sd_repcat_cnames=c("Replicate", "Category")
sd_effpv_cnames=c("Effect", "P-value");
sd_repcat_mat=matrix(NA, nrow=0, ncol=length(sd_repcat_cnames));
sd_effpv_mat=matrix(NA, nrow=0, ncol=length(sd_effpv_cnames));
colnames(sd_repcat_mat)=sd_repcat_cnames
colnames(sd_effpv_mat)=sd_effpv_cnames;

for(repid in uniq_repl_names){
	
	par(mfrow=c(4,2));

	for(alrix in alr_cat_names){
		vals=cum_differences[[alrix]][[repid]];
		mags=max(abs(vals));
		medn=median(vals);
		pad=mags/10;
		hist(vals, breaks=seq(-mags-pad, mags+pad, length.out=20), 
			main=alrix, xlim=c(-mags, mags));
		abline(v=0, col="black", lty=2, lwd=2);
		
		tres=t.test(vals, alternative="two.sided");
		if(tres$p.value<.05){
			if(medn>0){
				dircol="red";
			}else{
				dircol="blue";
			}

			sd_repcat_mat=rbind(sd_repcat_mat, 
				c(repid, alrix));
			sd_effpv_mat=rbind(sd_effpv_mat, 
				c(mean(vals), tres$p.value));

		}else{
			dircol="grey";
		}

		

		abline(v=medn, col=dircol, lwd=3);
		if(par()$page){
			mtext(repid, outer=T, font=2, cex=1.5, line=1.5);
		}
	}

	mtext(repid, outer=T, font=2, cex=1.5, line=1.5);
}

###############################################################################
# Output Tables

pval_order=order(sd_effpv_mat[,"P-value"]);
sd_effpv_mat=sd_effpv_mat[pval_order,,drop=F];
sd_repcat_mat=sd_repcat_mat[pval_order,,drop=F];

cat_order=order(sd_repcat_mat[,"Category"]);
repl_order=order(sd_repcat_mat[,"Replicate"]);
effect_order=order(sd_effpv_mat[,"Effect"]);

num_signf=nrow(sd_effpv_mat);
linenum=paste(1:num_signf, ".)", sep="");

outmat_bypv=cbind(sd_repcat_mat, signif(sd_effpv_mat, 4));
colnames(outmat_bypv)=c(sd_repcat_cnames, sd_effpv_cnames);

# By pvalue
rownames(outmat_bypv)=linenum;
plot_text(c(
	capture.output(print(outmat_bypv, quote=F))
	));
mtext("By P-value", outer=T, font=2, cex=1.5, line=1.5);

# By category
outmat=outmat_bypv[cat_order,,drop=F];
rownames(outmat)=linenum;
plot_text(c(
	capture.output(print(outmat, quote=F))
	));
mtext("By Category", outer=T, font=2, cex=1.5, line=1.5);

# By Effect 
outmat=outmat_bypv[effect_order,,drop=F];
rownames(outmat)=linenum;
plot_text(c(
	capture.output(print(outmat, quote=F))
	));
mtext("By Effect", outer=T, font=2, cex=1.5, line=1.5);

# By Replicate 
outmat=outmat_bypv[repl_order,,drop=F];
rownames(outmat)=linenum;
plot_text(c(
	capture.output(print(outmat, quote=F))
	));
mtext("By Replicate", outer=T, font=2, cex=1.5, line=1.5);

###############################################################################
# Plot depth vs Error

par(mar=c(5, 4, 3, 1));
for(repid in uniq_repl_names){
	
	par(mfrow=c(4,1));

	for(alrix in alr_cat_names){
		vals=cum_differences[[alrix]][[repid]];
		mags=max(abs(vals));

		plot(cum_depths[[repid]], vals, main=alrix, ylim=c(-mags*1.1, mags*1.1),
			xlab="reads/sample", ylab="estimated error");

		lws=lowess(cum_depths[[repid]], vals);
		points(lws$x, lws$y, type="l", col="red");
		abline(h=0, lty=2, col="blue");

		if(par()$page){
			mtext(repid, outer=T, font=2, cex=1.5, line=1.5);
		}
	}

	mtext(repid, outer=T, font=2, cex=1.5, line=1.5);

}

##############################################################################

cat("Done.\n");
print(warnings());
q(status=0);
