#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
library(car);
options(useFancyQuotes=F);

source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');

params=c(
	"summary_file", "s", 1, "character",
	"factor_file", "f", 1, "character",
	"sample_coln", "S", 1, "character",
	"replicate_coln", "R", 1, "character",
	"outputroot", "o", 1, "character"
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
	!length(opt$replicate_coln) 
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


##############################################################################

# Load matrix
counts=load_summary_file(SummaryFile);
num_taxa=ncol(counts);
num_samples=nrow(counts);
#print(counts);

# Shorten names:
full_names=colnames(counts);
splits=strsplit(full_names, " ");
short_names=character();
for(i in 1:length(full_names)){
	short_names[i]=tail(splits[[i]], 1);
}
colnames(counts)=short_names;

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

normalized=normalize(counts);

##############################################################################

pdf(paste(OutputRoot, ".rep_analys.div.pdf", sep=""), height=14, width=8.5);

plot_text(c(
	paste("Summary File : ", SummaryFile, sep=""),
	paste("Factors File: ", FactorsFile, sep=""),
	paste("Output File: ", OutputRoot, sep=""),
	paste("Replicate ID Column Name: ", Replicate_ColumnName, "\n", sep=""),
	paste("Sample Name Column Name: ", Sample_ColumnName, "\n", sep="")
));

##############################################################################
# For each sample name

replicate_group=factors[, Replicate_ColumnName];
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

diversity_bs=function(depth, distribution, div_fun, num_bs){

	samps=rmultinom(num_bs, depth, distribution);

	# Columns sum up to the depth
	# Number of columns equals num_bs

	inst_dist=samps/depth;
	
	div_arr=numeric(num_bs);
	for(i in 1:num_bs){
		div_arr[i]=div_fun(inst_dist[,i]);
	}

	return(div_arr);

}

diversity_95ci_rarefaction=function(depth_range, distribution, div_fun, num_bs=80){

	num_depths=length(depth_range);
	ci95_matrix=matrix(NA, ncol=3, nrow=num_depths);
	colnames(ci95_matrix)=c("lb95", "median", "ub95");

	dpix=1;
	bs_div_vals=numeric(num_bs);
	for(dep in depth_range){
		
		bs_div_vals=diversity_bs(dep, distribution, div_fun, num_bs);
		ci95_matrix[dpix,]=quantile(bs_div_vals, c(.025, .5, .975)); 

		dpix=dpix+1;
	}

	raref_results=list();
	raref_results[["ci95_matrix"]]=ci95_matrix;
	raref_results[["div_at_max"]]=bs_div_vals; 
	# Store bs values at greatest depth
	
	return(raref_results);

}

plot_combined_rarefaction=function(depth_ranges, rarefy_list, depth_markers, samp_colors){

	samp_names=names(rarefy_list);

	num_samp=length(samp_names);
	max_div=0;
	min_div=NA;

	max_med=numeric(num_samp);
	names(max_med)=samp_names;
	for(sampnm in samp_names){
		max_div=max(max_div, rarefy_list[[sampnm]]);
		min_div=min(min_div, rarefy_list[[sampnm]], na.rm=T);

		max_med[sampnm]=max(rarefy_list[[sampnm]][,2]);
	}

	samp_order_by_div=order(max_med, decreasing=T);

	max_depth=max(depth_ranges);
	num_depths=length(depth_ranges);
	
	cat("Max Diversity: ", max_div, "\n");
	
	plot(0,0, type="n", xlim=c(0, max_depth), ylim=c(min_div/1.1, max_div*1.1),
		xlab="Depth", ylab="Diversity");
	
	ix=0;

	abline(v=depth_markers[samp_names], col=samp_colors[samp_names]);

	# Draw CI lines first
	for(sampnm in samp_names[samp_order_by_div]){
		div_mat=rarefy_list[[sampnm]];
		cur_col=samp_colors[sampnm];	

		for(d in 1:num_depths){
			points(depth_ranges[d], div_mat[d, 1], 
				col=cur_col, bg=cur_col, type="p", pch=24, cex=1);
			points(depth_ranges[d], div_mat[d, 3], 
				col=cur_col, bg=cur_col, type="p", pch=25, cex=1);

			points(rep(depth_ranges[d], 2), div_mat[d, c(1,3)], 
				col=cur_col, type="l", lwd=num_samp-ix+.25);
		}

		ix=ix+1;
	}

	# Draw medians over CI lines for clarity
	for(sampnm in samp_names[samp_order_by_div]){
		div_mat=rarefy_list[[sampnm]];
		cur_col=samp_colors[sampnm];	

		# Draw median line
		points(depth_ranges, div_mat[,2], col=cur_col, type="l", lwd=5);
		points(depth_ranges, div_mat[,2], col="black", type="l", lwd=.25);

	}
}

##############################################################################

calc_shannon=function(distr){
	zeroFreeNorm=distr[distr>0];
	shan=-sum(zeroFreeNorm*log(zeroFreeNorm));
	return(shan);
}

calc_tail=function(x){
        sorted=sort(x, decreasing=TRUE);
	sumx=sum(x)
	if(sumx==0){return(0)}
        norm=sorted/sumx;
        n=length(norm);
        tail=0;
        for(i in 1:n){
                tail=tail + norm[i]*((i-1)^2);
        }
        return(sqrt(tail));
}

compare_indiv_to_group=function(bs_div_list){
	print(bs_div_list);
	
	ind_names=names(bs_div_list);

	target_diff=list();
	for(target_indv in ind_names){
		
		grp_combined=numeric();
		for(indiv in ind_names){
			if(indiv != target_indv){
				grp_combined=c(grp_combined, bs_div_list[[indiv]]);	
			}
		}

		samp_size=2000;

		target_val=sample(bs_div_list[[target_indv]], samp_size, replace=T);
		group_val=sample(grp_combined, samp_size, replace=T);

		diff=target_val-group_val;
		mean_diff=mean(diff);

		num_gt_zero=sum(diff>0);
		prob_gt_zero=num_gt_zero/samp_size;
		if(prob_gt_zero>.5){
			pval=(1-prob_gt_zero)*2;
		}else{
			pval=2*prob_gt_zero;
		}
	
		cat(target_indv, ":\n");
		cat(" Diff: ", mean_diff, "\n");
		cat(" Pvalue: ", pval, "\n");

		comp_res=list();
		comp_res[["diff"]]=mean_diff;
		comp_res[["pval"]]=pval;

		target_diff[[target_indv]]=comp_res;

	}

	return(target_diff);
}

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

plot_diff_from_group=function(diff_rec, ind_cols){
	print(diff_rec);

	ind_nm=names(diff_rec);
	num_ind=length(ind_nm);

	diff_arr=numeric(num_ind);
	names(diff_arr)=ind_nm;
	pval_arr=diff_arr;

	for(ind in ind_nm){
		diff_arr[ind]=diff_rec[[ind]][["diff"]];
		pval_arr[ind]=diff_rec[[ind]][["pval"]];
	}

	max_diff=max(diff_arr);
	min_diff=min(diff_arr);
	range=max_diff-min_diff
	mids=barplot(diff_arr, col=ind_cols[ind_nm], 
		ylim=c(min_diff-range/10, max_diff+range/10),
		ylab="Difference from Mean\n(Excluding Sample of Interest)");

	pval_text=character(num_ind);
	names(pval_text)=ind_nm;
	for(ind in ind_nm){
		sc=signf_char(pval_arr[ind]);
		pval_text[ind]=paste("p = ", round(pval_arr[ind],3), " ", sc, sep="");
	}	
	
	text(mids, diff_arr, labels=pval_text, cex=1.1, pos=3);
	


}

for(cur_sample_name in uniq_samp_names){

	cat("Analyzing Sample Group: ", cur_sample_name, "\n");
	
	samp_ix=(sample_names==cur_sample_name);
	samp_ids=sample_ids[samp_ix];
	num_samp=length(samp_ids);
	if(num_samp<2){
		next;
	}

	colors=rainbow(num_samp, start=0, end=.8);
	names(colors)=samp_ids;

	cat("Samples ID in Sample Group:\n");
	print(samp_ids);
	cur_counts=counts[samp_ids,, drop=F];
	cur_normal=normalized[samp_ids,, drop=F];

	read_depths=apply(cur_counts, 1, sum);
	print(read_depths);
	max_depth=max(read_depths);

	gap=150;
	if(max_depth<gap){
		gap=max_depth/10;
	}
	depth_range=sort(unique(c(seq(gap, max_depth, gap), read_depths)));
	cat("Rarefying at: \n");
	print(depth_range);

	cat("Calculating Shannon...\n");
	raref_95ci_shannon_list=list();
	bs_div_at_max=list();
	for(cur_samp_ids in samp_ids){
		raref_res=diversity_95ci_rarefaction(
			depth_range, cur_normal[cur_samp_ids,], calc_shannon);
		raref_95ci_shannon_list[[cur_samp_ids]]=raref_res[["ci95_matrix"]];
		bs_div_at_max[[cur_samp_ids]]=raref_res[["div_at_max"]];
	}
	shan_diff=compare_indiv_to_group(bs_div_at_max);

	cat("Calculating Tail...\n");
	raref_95ci_tail_list=list();
	for(cur_samp_ids in samp_ids){
		raref_res=diversity_95ci_rarefaction(
			depth_range, cur_normal[cur_samp_ids,], calc_tail);
		raref_95ci_tail_list[[cur_samp_ids]]=raref_res[["ci95_matrix"]];
		bs_div_at_max[[cur_samp_ids]]=raref_res[["div_at_max"]];
	}
	tail_diff=compare_indiv_to_group(bs_div_at_max);
	
	cat("ok.\n");

	par(oma=c(1,1,4,1));
	par(mfrow=c(5,1));
	par(mar=c(5,20,5,1));
	barplot(read_depths, horiz=T, las=1, col=colors, xlab="Reads/Sample", main="Replicate Sequencing Depth");
	
	par(mar=c(4,4,3,1));
	plot_combined_rarefaction(depth_range, raref_95ci_shannon_list, read_depths, colors);
	title(main="Shannon Bootstrapped Rarefaction Curve");
	plot_diff_from_group(shan_diff, colors);

	plot_combined_rarefaction(depth_range, raref_95ci_tail_list, read_depths, colors);
	title(main="Tail Bootstrapped Rarefaction Curve");
	plot_diff_from_group(tail_diff, colors);

	mtext(cur_sample_name, outer=T, font=2, cex=1.5);

	cat("\n");

}
 


# Bootstrap/Rarefy Samples
# Plot overlapping curves
# Perform pairwise comparison between replicates, relative to smaller depth and relative to larger depth

























quit();


##############################################################################
# Compute diversity indices

div_names=c("Tail", "Shannon", "Simpson", "Evenness", "SimpsonsRecip");
num_div_idx=length(div_names);

div_mat=matrix(0, nrow=num_samples, ncol=num_div_idx);
colnames(div_mat)=div_names;
rownames(div_mat)=rownames(normalized);

cat("Computing diversity indices across samples.\n");
for(i in 1:num_samples){
	curNorm=normalized[i,];
	zeroFreeNorm=curNorm[curNorm>0]
	div_mat[i,"Tail"]=tail_statistic(zeroFreeNorm);
	div_mat[i,"Shannon"]=-sum(zeroFreeNorm*log(zeroFreeNorm));
	div_mat[i,"Simpson"]=1-sum(curNorm^2);
	div_mat[i,"Evenness"]=div_mat[i,"Shannon"]/log(length(zeroFreeNorm));
	div_mat[i,"SimpsonsRecip"]=1/sum(curNorm^2);
}

# Set evenness for NAs to 0.  Assume evenness is very low, but it's degenerate
# because of insufficient sequencing depth.
evenness=div_mat[,"Evenness"];
div_mat[is.na(evenness),"Evenness"]=0;

cat("Plotting histograms of raw diversity indices.\n");
par(mfrow=c(3,2));
par(oma=c(1, 1, 1, 1));
par(mar=c(5,4,4,2));
for(i in 1:num_div_idx){
	hist(div_mat[, div_names[i]], main=div_names[i], xlab=div_names[i], 
		breaks=15);
}
mtext("Sample Diversity Indices Distributions", outer=T);

##############################################################################

##############################################################################
# Perform Box-Cox transformation

lambda=numeric(num_div_idx);
transformed=matrix(0, nrow=num_samples, ncol=num_div_idx);
colnames(transformed)=div_names;

cat("Computing and plotting Box-Cox transformations.\n\n");
BOX_COX_SEARCH_RANGE=2;
SEARCH_FREQUENCY=10;
par(mfrow=c(3,2));
for(i in 1:num_div_idx){

	cat("Computing lambda for: ", div_names[i], "\n", sep="");

	raw=div_mat[, div_names[i]];

	# Look for lambda approximately
	lambda_found=FALSE;
	lambda_start=-BOX_COX_SEARCH_RANGE;
	lambda_end=BOX_COX_SEARCH_RANGE;
	search_trial=0;
	while(!lambda_found){

		search_trial=search_trial+1;

		# Perform trial boxcox search
		cat("[", search_trial, "] Search range: ", lambda_start, " to ", lambda_end, "\n", sep="");
		cat("Model: ", model_string, "\n");


		degenerates=which(is.nan(raw));
		if(length(degenerates)){
			cat("*************************************\n");
			cat("*  Degenerate diversities: \n");
			print(raw[degenerates]);
			cat("*************************************\n");
			raw[degenerates]=0;
		}

		zero_ix=(raw==0);
		if(any(zero_ix)){
			cat("Zeros found in diversity:\n");
			print(raw[zero_ix]);
			min_nonzero=min(raw[!zero_ix]);
			min_subst=min_nonzero/10;
			cat("Adding ", min_nonzero, " / 10 = ", min_subst, " to all responses...\n", sep="");
			raw=raw+min_subst;
			print(sort(raw));
		}
		
		bc=tryCatch({
			boxcox(as.formula(model_string), data=factors, 
				lambda=seq(lambda_start, lambda_end, length.out=SEARCH_FREQUENCY),
				plotit=FALSE);
			}, 
			error=function(cond){
				cat("Error finding Box-Cox transform lambda value.\n");
				print(cond);
				return(NULL);
			});

		if(is.null(bc)){
			approx_lambda=1;
			lambda_found=TRUE;
		}else{

			# Determine where to expand search depending if peak is on left or right
			max_idx=which(bc$y==max(bc$y));
			if(max_idx==1){
				#lambda_end=lambda_start+1;
				lambda_start=lambda_start-(BOX_COX_SEARCH_RANGE^search_trial);
			}else if(max_idx==length(bc$y)){
				#lambda_start=lambda_end-1;
				lambda_end=lambda_end+BOX_COX_SEARCH_RANGE^search_trial;
			}else{
				search_tolerance=(lambda_end-lambda_start)/SEARCH_FREQUENCY;
				approx_lambda=bc$x[max_idx];
				lambda_found=TRUE;
				cat("Lambda found around: ", approx_lambda, "\n");	
				cat("   Search tolerance: ", search_tolerance, "\n");	
			}
		}
	}


	if(approx_lambda!=1){
		# Rerun and plot the lambda search with a smaller increments
		par(mar=c(5,4,4,2));
		cat("Refining search around: ", approx_lambda-search_tolerance, " - ", approx_lambda+search_tolerance, "\n");
		bc=boxcox(as.formula(model_string), data=factors, 
			lambda=seq(approx_lambda-search_tolerance, approx_lambda+search_tolerance, length.out=40));
		title(main=div_names[i]);

		# Store find grain results
		max_idx=which(bc$y==max(bc$y));
		lambda[i]=bc$x[max_idx];
	}else{
		cat("Warning: Box-Cox Lambda value not found.  Going with 1 (i.e. no transform)\n");
		cat("\n");
		cat("If the sample size was close to the number of variables in the model, it is likely\n");
		cat("that the model was overfit, so it was not possible to find lambda that minimizes the residuals.\n");
		cat("\n");
		cat("Num Model Variables: ", num_model_variables, "\n");
		cat("Num Samples: ", num_samples, "\n");
		cat("\n");
		lambda[i]=1;

	}
	
	cat(div_names[i], ": Box-Cox Transformation Lambda = ", lambda[i], "\n\n", sep="");

	# Apply transform to raw data
	if(lambda[i]==0){
		transformed[, div_names[i]]=log(raw);
	}else{
		transformed[, div_names[i]]=((raw^lambda[i])-1)/lambda[i];
	}
}
mtext("Box-Cox Lambda Intervals", outer=T);

cat("Plotting transformed histograms...\n");
par(mfrow=c(3,2));
par(oma=c(1, 1, 1, 1));
par(mar=c(5,4,4,2));
for(i in 1:num_div_idx){
	hist(transformed[, div_names[i]], 
		main=paste(div_names[i], " (lambda=", sprintf("%3.2f", lambda[i]), ")", sep=""),
		xlab=div_names[i], 
		breaks=15);
}
mtext("Transformed Sample Diversity Indices Distributions", outer=T);

##############################################################################
# Output reference factor levels

text=character();
text[1]="Reference factor levels:";
text[2]="";

cat("Outputing factor levels...\n");
for(i in 1:num_factors){
	fact_levels=levels(factors[,i]);
	if(!is.null(fact_levels)){
		fact_info=paste(factor_names[i], ": ", fact_levels[1], sep="");	
	}else{
		fact_info=paste(factor_names[i], ": None (ordered factor)", sep="");	
	}
	text=c(text, fact_info);
}

text=c(text, "");
text=c(text, paste("Number of Samples: ", num_samples, sep=""));
text=c(text, "");
text=c(text, "Description of Factor Levels and Samples:");
text=c(text, capture.output(summary(factors)));

plot_text(text);

##############################################################################

bin_continuous_values=function(values, num_bins=10){
        minv=min(values);
        maxv=max(values);
        range=maxv-minv;
        # Map values between 0 and 1
        prop=(values-minv)/range;
        # Scale value up to bin, and round, to quantize
        closest=round(prop*(num_bins-1),0);
        log10range=log10(range);
        trunc=signif(closest/(num_bins-1)*range+minv, 5)
	bin_range=mean(diff(sort(unique(trunc))))/2;

	lb=trunc-bin_range;
	ub=trunc+bin_range;
        # Remap values to original range and location
        return(paste("(", lb, ", ", ub, ")", sep=""));
}

plot_diversity_with_factors=function(raw, factors, model_string, stat_name, bin_cont=5){

	palette(c(
		"red",
		"green",
		"blue",
		"orange",
		"purple",
		"black",
		"brown",
		"darkgoldenrod3"
	));
	
	# Extract out predictors
	predictor_string=strsplit(model_string, "~")[[1]][2];
	pred_arr=strsplit(predictor_string, "\\+")[[1]];
	
	# Remove interaction terms
	interact_ix=grep(":", pred_arr);
	if(length(interact_ix)){
		pred_arr=pred_arr[-interact_ix]
	}

	pred_arr=gsub(" ", "", pred_arr);
	num_pred=length(pred_arr);

	num_values=length(raw);
	raw_range=range(raw);

	# Sort records by magnitude of raw values
	sort_ix=order(raw);
	raw=raw[sort_ix];
	factors=factors[sort_ix, , drop=F];

	par(mar=c(5,5,5,5));
	zeros=rep(0, num_values);
	
	min_spacing=(raw_range[2]-raw_range[1])/60;
	cat("Min spacing: ", min_spacing, "\n");
	#abline(h=seq(raw_range[1], raw_range[2], length.out=50));

	adj=min_spacing/10;
	cat("Adjustment Increment: ", adj, "\n");
	adj_pos=raw;

	for(adj_ix in 1:10000){
		adjusted=F
		for(forw in 1:(num_values-1)){
			if(abs(adj_pos[forw]-adj_pos[forw+1])<min_spacing){
				adj_pos[forw]=adj_pos[forw]-adj;
				adjusted=T;
				break;
			}	
		}
		for(rev in (num_values:2)){
			if(abs(adj_pos[rev]-adj_pos[rev-1])<min_spacing){
				adj_pos[rev]=adj_pos[rev]+adj;
				adjusted=T;
				break;
			}	
		}
		if(!adjusted){
			cat("Done adjusting at iteration: ", adj, "\n");
			break;	
		}
	}
	
	#cat("Raw:\n");
	#print(raw);
	#cat("Adj:\n");
	#print(adj_pos);

	extra_sample_space=2;
	predictors_per_plot=5;
	pred_names=colnames(factors);
	raw_adj_range=range(c(raw, adj_pos));

	num_plots=num_pred %/% predictors_per_plot + 1;

	for(plot_ix in 1:num_plots){
		# Plot the samples on the y axis
		plot(0, xlab="", ylab=stat_name, type="n", xaxt="n", main=stat_name,
			ylim=c(raw_adj_range[1], raw_adj_range[2]),
			xlim=c(-1, predictors_per_plot+1+extra_sample_space));

		# Plot 
		for(i in 1:num_values){
			lines(x=c(-1, -.75), y=c(raw[i], raw[i]));
			lines(x=c(-.75, -.25), y=c(raw[i], adj_pos[i]));

		}
		# Label samples
		text(zeros-.20, adj_pos, label=names(raw), pos=4, cex=.5);

		# Label predictors
		predictors_per_plot=min(predictors_per_plot, num_pred);

		if(predictors_per_plot>0){
			for(j in 1:predictors_per_plot){
				pred_ix=((plot_ix-1)*predictors_per_plot)+(j-1)+1;
			
				if(!is.na(pred_arr[pred_ix]) && !length(grep(":", pred_arr[pred_ix]))){

					fact_val=factors[, pred_arr[pred_ix]];
					uniq_val=unique(fact_val);

					if(length(uniq_val)>=bin_cont && !is.factor(fact_val)){
						factor=signif(fact_val, 4);
						coloring=as.numeric(as.factor(bin_continuous_values(fact_val, num_bins=bin_cont)));
					}else{
						factor=as.factor(fact_val);
						coloring=as.numeric(factor);
					}

					if(is.factor(factor)){
						factor=substr(factor, 1, 10);
					}

					text(zeros+j+extra_sample_space, adj_pos, label=factor, 
						col=coloring,
						cex=.5
						);
				}
				# label predictor/factor name
				lab_cex=min(1, 8/nchar(pred_arr[pred_ix]));
				text(j+extra_sample_space, raw_adj_range[2], 
					pred_arr[pred_ix], cex=lab_cex*.75, col="black", family="", font=2, pos=3);
			}
		}
	}

}

###############################################################################

plot_overlapping_histograms=function(raw, factors, model_string, title, bin_cont=5){

	orig.par=par(no.readonly=T);
	palette(c(
		"red",
		"green",
		"blue",
		"orange",
		"purple",
		"black",
		"brown",
		"darkgoldenrod3"
	));

	# Extract out predictors
	predictor_string=strsplit(model_string, "~")[[1]][2];
	pred_arr=strsplit(predictor_string, "\\+")[[1]];
	pred_arr=gsub(" ", "", pred_arr);
	num_pred=length(pred_arr);
	num_values=length(raw);
	raw_range=range(raw);

	cat("Num Predictors:", num_pred, "\n");
	cat("Range: ", raw_range[1], " - ", raw_range[2], "\n");

	factor_names=colnames(factors);

	#par(mfrow=c(4,1));
	layout_mat=matrix(c(1,1,2,3,3,4,5,5,6), ncol=1);
	layout(layout_mat);

	for(pix in 1:num_pred){
		cur_pred=pred_arr[pix];
		cat("  Working on: ", cur_pred, "\n");

		if(any(cur_pred==factor_names)){

			# Get levels for each value
			cur_fact_val=factors[,cur_pred];
			is_factor=is.factor(cur_fact_val);

			num_uniq=length(unique(cur_fact_val));
			if(num_uniq>=bin_cont && !is_factor){
				cur_fact_val=as.factor(bin_continuous_values(cur_fact_val, num_bins=bin_cont));
                        }else{
				cur_fact_val=as.factor(cur_fact_val);
			}

			#num_bins=nclass.Sturges(raw)*2;			
			#overall_hist=hist(raw, breaks=num_bins, plot=F);
			#print(overall_hist);
			#cat("  Num bins: ", num_bins, "\n");

			levels=levels(cur_fact_val);
			num_levels=length(levels);
			print(cur_fact_val);
			cat("Num Levels: ", num_levels, "\n");
			
			# Compute the density for each level
			dens_list=list(num_levels);
			level_samp_size=numeric();
			max_dens=0;
			for(lix in 1:num_levels){
				level_val=raw[cur_fact_val==levels[lix]];
				#hist_list[[levels[lix]]]=hist(level_val, breaks=overall_hist$breaks, plot=F);
				num_samp_at_level=length(level_val);
				level_samp_size[lix]=num_samp_at_level;
				if(num_samp_at_level==0){
					dens_list[[lix]]=NA;
				}else{
					if(num_samp_at_level==1){
						level_val=c(level_val, level_val);
					}

					dens_list[[lix]]=density(level_val);
					max_dens=max(max_dens, dens_list[[lix]]$y);	
				}
			}

			# Open a blank plot
			par(mar=c(4,3,2,1));
			plot(0,0, type="n", xlim=raw_range, ylim=c(0,max_dens), 
				main=cur_pred,
				xlab=title, ylab="Density",
				bty="l");

			# Draw the curves for each factor level
			for(lix in 1:num_levels){
				dens=dens_list[[lix]];
				if(!any(is.na(dens))){
					points(dens, type="l", col=lix);
				}
			}

			# Plot the legend for each factor
			legend_labels=character();
			for(lix in 1:num_levels){
				legend_labels[lix]=paste(levels[lix], " [n=", level_samp_size[lix], "] ", sep="");
			}

			par(mar=c(0,0,0,0));
			plot(0,0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", bty="n");
			legend(0,1, legend=legend_labels, fill=1:num_levels, bty="n");

		}else{
			cat("Not ploting interactions...\n");
		}
	}

	par(orig.par);
	
}

##############################################################################

# Matrix to store R^2's
cat("Allocating R^2 Matrix...\n");
rsqrd_mat=matrix(0, nrow=num_div_idx, ncol=2);
colnames(rsqrd_mat)=c("R^2", "Adj. R^2");
rownames(rsqrd_mat)=div_names;

# Build model matrix so we know how many coefficients to expect.
if(ModelFormula==""){
	model_string= paste("trans ~", 
	paste(factor_names, collapse=" + "));
}

cat("\nFitting this regression model: ", model_string, "\n");

trans=rep(0, nrow(factors));
model_matrix=model.matrix(as.formula(model_string), data=factors);
num_coeff=ncol(model_matrix);

# Matrices to store coefficents and p-values
cat("Allocating coeff and pval matrices...\n");
coeff_matrix=matrix(NA, nrow=num_coeff, ncol=num_div_idx);
pval_matrix=matrix(NA, nrow=num_coeff, ncol=num_div_idx);

colnames(coeff_matrix)=div_names;
rownames(coeff_matrix)=colnames(model_matrix);

colnames(pval_matrix)=div_names;
rownames(pval_matrix)=colnames(model_matrix);


# Fit regression model and ANOVA analysis
for(i in 1:num_div_idx){

	cat("Working on: ", div_names[i], "\n", sep="");

	raw=div_mat[, div_names[i]];

	plot_diversity_with_factors(raw, factors, model_string, div_names[i], 6);
	plot_overlapping_histograms(raw, factors, model_string, title=div_names[i], 6);

	trans=transformed[, div_names[i]];
	print(model.matrix(as.formula(model_string), data=factors));
	fit=lm(as.formula(model_string), data=factors);
	print(fit);

	summ=summary(fit);
	print(summ);

	rsqrd_mat[i, 1]=summ$r.squared;
	rsqrd_mat[i, 2]=summ$adj.r.squared;

	# Output ANOVA
	if(ModelFormula!=""){
		aov_model=paste("trans ~", ModelFormula);
	}else{
		aov_model=model_string;
	}

	aov_sumtext=capture.output(summary(aov(as.formula(aov_model), data=factors)));

	plot_text(c(
		paste("ANOVA for: ", div_names[i]),
		"",
		aov_model,
		"", aov_sumtext));

	# Output regression coefficients and significance
	sumtext=(capture.output(summary(fit)));
	plot_text(c(
		paste("Multiple Regression for: ", div_names[i]), 
		"",
		model_string,
		"",
		paste("Diversity Mean: ", mean(raw)),
		paste("Diversity Stderr: ", sd(raw)/sqrt(length(raw))),
		"", sumtext));

	# Generate marginal model plots
	mmps(fit);



	sum_fit=summary(fit);

	computable_coef_names=names(sum_fit$coefficients[,"Estimate"]);

	coeff_matrix[computable_coef_names,i]=sum_fit$coefficients[computable_coef_names,"Estimate"];
	pval_matrix[computable_coef_names,i]=sum_fit$coefficients[computable_coef_names,"Pr(>|t|)"];
}

# R-squred summary
print(rsqrd_mat);
plot_correl_heatmap(t(rsqrd_mat), title="R^2 Values");

# Remove intercept from matrix
intercept_ix=which(rownames(coeff_matrix)=="(Intercept)");
coeff_matrix=coeff_matrix[-intercept_ix, , drop=F];
pval_matrix=pval_matrix[-intercept_ix, , drop=F];

plot_correl_heatmap(coeff_matrix, title="Coefficients");
plot_correl_heatmap(pval_matrix, title="P-values");

significant_coeff=(pval_matrix<0.05)*coeff_matrix;
plot_correl_heatmap(significant_coeff, title="Significant Coefficients", noPrintZeros=T);


dev.off();

##############################################################################

fh=file(paste(OutputRoot, ".div_as_resp.regr_stats.tsv", sep=""), "w");

# Output Header
cat(file=fh, paste(c("Coefficients", "Estimates:", div_names, "p-values:", div_names), collapse="\t"));
cat(file=fh, "\n");

# Output values
coeff_names=rownames(pval_matrix);
for(i in 1:length(coeff_names)){
	cat(file=fh, 
		coeff_names[i], 
		"",
		paste(coeff_matrix[i,], collapse="\t"),
		"",
		paste(pval_matrix[i,], collapse="\t"),
		sep="\t");
	cat(file=fh, "\n");
}

close(fh);

##############################################################################

# Write coefficient p-values to file
write.table(t(pval_matrix), file=paste(OutputRoot, ".div_as_resp.pvals.tsv", sep=""),
        sep="\t", quote=F, col.names=NA, row.names=T);

write.table(t(coeff_matrix), file=paste(OutputRoot, ".div_as_resp.coefs.tsv", sep=""),
        sep="\t", quote=F, col.names=NA, row.names=T);
	

##############################################################################

# Write diversity matrix to file
write.table(div_mat, file=paste(OutputRoot, ".diversity_indices.tsv", sep=""),
	sep="\t", quote=F, col.names=NA, row.names=T);


##############################################################################

cat("Done.\n");
print(warnings());
q(status=0);
