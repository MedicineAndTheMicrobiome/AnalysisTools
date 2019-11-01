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

plot_diff_from_group=function(diff_rec, ind_cols, repl_ids){
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
		names.arg=repl_ids[ind_nm],
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

kept_sample_depth_vals =list();
replicate_differences_tail=list();
replicate_differences_shan=list();
for(repnm in uniq_repl_names){
	replicate_differences_tail[[repnm]]=numeric();
	replicate_differences_shan[[repnm]]=numeric();
	kept_sample_depth_vals[[repnm]]=numeric();
}

iter=0;
for(cur_sample_name in uniq_samp_names){

	cat("Analyzing Sample Group: ", cur_sample_name, "\n");
	
	samp_ix=(sample_names==cur_sample_name);
	samp_ids=sample_ids[samp_ix];
	repl_ids=replicate_group[samp_ix];
	names(repl_ids)=samp_ids;

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
	plot_diff_from_group(shan_diff, colors, repl_ids);

	plot_combined_rarefaction(depth_range, raref_95ci_tail_list, read_depths, colors);
	title(main="Tail Bootstrapped Rarefaction Curve");
	plot_diff_from_group(tail_diff, colors, repl_ids);

	mtext(cur_sample_name, outer=T, font=2, cex=1.5);

	# Accumulate

	for(s in samp_ids){
		rep_id=repl_ids[s];

		sdff=shan_diff[[s]][["diff"]];
		replicate_differences_shan[[rep_id]]=c(replicate_differences_shan[[rep_id]], sdff);

		sdff=tail_diff[[s]][["diff"]];
		replicate_differences_tail[[rep_id]]=c(replicate_differences_tail[[rep_id]], sdff);

		kept_sample_depth_vals[[rep_id]] =c(kept_sample_depth_vals[[rep_id]], read_depths[s]);
	}

	cat("\n");

	iter=iter+1;
	if(iter==20){
		#break;
	}
}

plot_repl_diff=function(rep_diff_list, divname){
	
	rep_nms=names(rep_diff_list);
	num_rps=length(rep_nms);
	par(mfrow=c(num_rps,1));

	min_diff=NA;
	max_diff=NA;
	for(rpnm in rep_nms){
		min_diff=min(rep_diff_list[[rpnm]], min_diff, na.rm=T);
		max_diff=max(rep_diff_list[[rpnm]], max_diff, na.rm=T);
	}

	cat("Diff Range: ", min_diff, " - ", max_diff, "\n");

	# Plot Histograms
	for(rpnm in rep_nms){
		hist(rep_diff_list[[rpnm]], xlim=c(min_diff, max_diff), breaks=20,
			xlab="Difference of diversity from Other Group Members",
			main=rpnm);
		abline(v=0, col="blue");
	}
	mtext(divname, outer=T, font=2, cex=2);

	# Plot depth vs diversity difference
	for(rpnm in rep_nms){
		plot(kept_sample_depth_vals[[rep_id]], rep_diff_list[[rpnm]], type="n",
			main=rpnm,
			ylim=c(min_diff, max_diff),
			xlab="Reads/Sample", ylab="Difference from Other Group Members"	
			);
		fit=lm(rep_diff_list[[rpnm]]~kept_sample_depth_vals[[rep_id]]);
		sumfit=summary(fit);
		slope=sumfit$coefficients[2,"Estimate"];
		pval=sumfit$coefficients[2,"Pr(>|t|)"];

		title(main=paste("Slope: ", signif(slope, 4), "  p-value: ", signif(pval, 3), sep=""), 
			line=.5, font.main=1);
		
		abline(fit, col="red");
		abline(h=0, col="blue");
		points(kept_sample_depth_vals[[rep_id]], rep_diff_list[[rpnm]]);
	}
	mtext(divname, outer=T, font=2, cex=2);

}

par(mar=c(5,4,4,1));
plot_repl_diff(replicate_differences_shan, "Shannon");
plot_repl_diff(replicate_differences_tail, "Tail");

##############################################################################

cat("Done.\n");
print(warnings());
q(status=0);
