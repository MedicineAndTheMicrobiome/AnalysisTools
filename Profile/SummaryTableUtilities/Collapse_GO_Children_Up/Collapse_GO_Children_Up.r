#!/usr/bin/env Rscript

###############################################################################

library('getopt');

source("~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r");
source("~/git/AnalysisTools/Metadata/OutputFileLibrary/OutputFileLibrary.r");

MIN_PROP_THRES=.95;

params=c(
	"summary_table", "s", 1, "character",
	"obo_table", "b", 1, "character",
	"output_file_root", "o", 1, "character",
	"min_true_threshold", "t", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-s <input summary_table.tsv file>\n",
	"	-b <input OBO table>\n",
	"	-o <output file root>\n",
	"	[-t <min proportion of True Positive/Negative, default=", MIN_PROP_THRES, ">]\n",
	"\n",	
	"This script will read in the summary table and OBO table and identify\n",
	"the leaf nodes that should be collapsed upward to its parent.\n",
	"\n",
	"For each leaf node the number of TN and TP will be calculated.\n",
	"If the proportion: (TP+PN)/NumSamp is >", MIN_PROP_THRES, 
	"%, then the leaf will not be collapsed\n",
	"\n",
	"The TP (true positives) are the number of non-zero data points\n",
	"The TN (true negatives) will be based on the mean abundance of the NZ\n",
	"  and the sequencing depth of that sample.\n",
	"\n",
	"\n");

if(
	!length(opt$summary_table) ||
	!length(opt$obo_table) ||
	!length(opt$output_file_root)
){
	cat(usage);
	q(status=-1);
}

options(width=200, useFancyQuotes=F);

###############################################################################
# Parameters

SummaryTable=opt$summary_table;
OBOTable=opt$ob_table;
OutputFileRoot=opt$output_file_root

cat("\n")
cat("Summary Table: ", SummaryTable, "\n");
cat("OBO Table: ", OBOTable, "\n");
cat("Output File Root: ", OutputFileRoot, "\n");
cat("\n");

pdf(paste(OutputFileRoot, ".clps_zero_chldn.pdf", sep=""), height=11, width=8.5);

###############################################################################
###############################################################################

load_obo_map=function(obo_map_fn){

        cat("Loading OBO Map File...\n");

        start_time=Sys.time();

        mat=read.table(InputOBOMap, quote="", comment.char="", sep="\t",
                row.names=NULL, header=T);

        end_time=Sys.time();
        exec_time=end_time-start_time;

        cat("Time: ", exec_time, " secs\n", sep="");
        cat("Num Rows: ", nrow(mat), " Loaded.\n");
        cat("OBO Map Column Names:\n");
        print(colnames(mat));
        cat("\n");

        return(mat);
}

#------------------------------------------------------------------------------

load_obo_to_tree=function(obo_tab){

        cat("Loading OBO Table into Tree...\n");
        start_time=Sys.time();

        num_ids=nrow(obo_tab);
        parent_to_child_tree=vector("list", length=num_ids);
        names(parent_to_child_tree)=obo_tab[,"id"];

        parent_presplit=strsplit(obo_tab[,"is_a"], ";");

        for(i in 1:num_ids){
                child_id=obo_tab[i,"id"];
                parent_ids=parent_presplit[[i]];
                for(pid in parent_ids){
                        parent_to_child_tree[[pid]]=c(
                                parent_to_child_tree[[pid]],
                                child_id);
                }
        }

        end_time=Sys.time();
        cat("ok.\n");
        exec_time=end_time-start_time;
        cat("Time: ", exec_time, " secs\n", sep="")
        return(parent_to_child_tree);

}

#------------------------------------------------------------------------------

calc_zero_stats=function(cts_mat, nrm_mat){

	# Compute num zeros in each category
	num_zeros=apply(cts_mat, 2, function(x){sum(x==0)});

	# Using the non-zero mean, will make it easier to call a true negative
	nz_mean=apply(nrm_mat, 2, 
		function(x){
			nz_ix=!(x==0);
			nz_mean=mean(x[nz_ix]);
		});

	# Using median of all, will result in many 0's, so the mean is necessary
	# nz_mean=apply(nrm_mat, 2, median);
		
	# Stats table 
	zstats=cbind(nz_mean, num_zeros);
	rownames(zstats)=colnames(nrm_mat);
	colnames(zstats)=c("NonZeroMean", "NumZeros");

	return(zstats);
}

estimate_TruePosNegs=function(zstats, nrm_mat, smp_dpt, false_neg_cutoff=.95, min_true_posneg_cutoff=.95){ 
	num_cat=nrow(zstats);
	cat_names=rownames(zstats);
	num_samples=length(smp_dpt);

	cat("Estimating True Pos+Negs:\n");
	cat("Num Cats: ", num_cat, "\n");
	cat("Num Samples: ", num_samples, "\n");

	tp_stat_names=c("NumTrPos", "NumTrNeg", "PropTrue");
	num_tp_stats=length(tp_stat_names);
	tp_tab=matrix(0, nrow=num_cat, ncol=num_tp_stats);
	colnames(tp_tab)=tp_stat_names;
	rownames(tp_tab)=cat_names;

	for(cat_ix in 1:num_cat){

		cur_cat=cat_names[cat_ix];
		cat("Analyzing: ", cur_cat, "\n");

		# Find zeros across samples for this category
		zero_cts_ix=(nrm_mat[,cur_cat]==0);
		# Get samp depth of zero counts
		smp_depth_of_zeros=smp_dpt[zero_cts_ix];

		# Calculate probabilty of being nonzero, given depth and mean non-zero 
		nz_mean=zstats[cur_cat, "NonZeroMean"];
		calc_prob_nz=function(depth){
			# Probability of being nonzero given depth and mean abundance
			probnz=1-pbinom(0, depth, nz_mean);
			return(probnz);
		}
		prob_non_zero=sapply(smp_depth_of_zeros, calc_prob_nz);

		# If prob of nonzero > cutoff, call it TRUE NEG 
		true_neg_ix=(prob_non_zero>false_neg_cutoff);
		num_true_neg=sum(true_neg_ix);

		# Calc number of nonzero categories
		num_true_pos=num_samples-zstats[cur_cat, "NumZeros"];
		
		# Proportion of samples that we are "sure" that are TRUE POS+NEG
		prop_true_np=(num_true_neg+num_true_pos)/num_samples;
		#cat("Num True Neg: ", num_true_neg, "\n");
		#cat("Num True Pos: ", num_true_pos, "\n");
		#cat("Prop True: ", prop_true_np, "\n");

		tp_tab[cur_cat,]=c(num_true_pos, num_true_neg, prop_true_np);

	}

	return(tp_tab);

}

plot_zero_stats=function(zstats, tp_mat=NULL, cutoff=NULL){

	par(mfrow=c(1,1));

	num_cat=nrow(zstats);
	if(is.null(tp_mat)){
		cols=rep("black", num_cat);
	}else{
		cols=rep("grey", num_cat);
		cols[tp_mat[,"PropTrue"]>=cutoff]="blue";
	}

	plot(log10(zstats[,"NonZeroMean"]), log10(zstats[,"NumZeros"]),
		col=cols,
		main="Num Zeros vs. Mean Abund",
		xlab="log10(Mean Abd of NonZeros)", ylab="log10(Number of Zeros)");

	if(!is.null(cutoff)){
		title(main=paste("Min 'True' Cutoff = ", cutoff), cex.main=.9, col.main="blue", line=.5);
	}
}

find_leaves=function(st, cutoff, child_to_parent){

}

###############################################################################
# Load files

counts_mat=load_summary_file(SummaryTable);

num_sumtab_cat=ncol(counts_mat);
num_sumtab_samp=nrow(counts_mat);

cat("Num categories: ", num_sumtab_cat, "\n", sep="");
cat("Num samples: ", num_sumtab_samp, "\n", sep="");
cat("\n");

norm_mat=normalize(counts_mat);

mean_cat_abd=apply(norm_mat, 2, mean);
cat_order_ix=order(mean_cat_abd, decreasing=T);
mean_cat_abd=mean_cat_abd[cat_order_ix];

samp_depths=apply(counts_mat, 1, sum);
samp_order_ix=order(samp_depths, decreasing=T);
samp_depths=samp_depths[samp_order_ix];

# Order by depth and relative abundance
counts_mat=counts_mat[samp_order_ix, cat_order_ix];
norm_mat=norm_mat[samp_order_ix, cat_order_ix];

#------------------------------------------------------------------------------

zero_stats_mat=calc_zero_stats(counts_mat, norm_mat);

plot_zero_stats(zero_stats_mat);

true_pos_mat=estimate_TruePosNegs(zero_stats_mat, norm_mat, samp_depths);

par(mfrow=c(3,1));
hist(true_pos_mat[,"NumTrPos"]/num_sumtab_samp, breaks=20,
	main="True Positives (Observed)", xlab="Prop True Pos");
hist(true_pos_mat[,"NumTrNeg"]/num_sumtab_samp, breaks=20,
	main="True Negatives (Estimated)", xlab="Prop True Neg");
hist(true_pos_mat[,"PropTrue"], breaks=20,
	main="Proportion True", xlab="Prop True (Pos+Neg)");

plot_zero_stats(zero_stats_mat, true_pos_mat, cutoff=.5);
plot_zero_stats(zero_stats_mat, true_pos_mat, cutoff=.75);
plot_zero_stats(zero_stats_mat, true_pos_mat, cutoff=.90);
plot_zero_stats(zero_stats_mat, true_pos_mat, cutoff=.95);

quit();



###############################################################################

###############################################################################

par(mfrow=c(3,1));
hist(num_zeros, main="0's in Categories", xlab="Number of Zeros");
hist(num_zeros/num_sumtab_samp, main="0's in Categories", xlab="Proportion of Zeros");
hist(log10(num_zeros), main="0's in Categories", xlab="Log10(Num Zeros)");

###############################################################################

quit();


###############################################################################
# Remove zero count categories

orig_counts=apply(counts_mat, 1, sum);
out_counts=apply(out_count_mat, 1, sum);

abs_diff=abs(orig_counts-out_counts);
max_abs_diff=max(abs_diff);
cat("Max differences between original and remapped: ", max_abs_diff, "\n");
non_fract_diff=(abs_diff>1);

if(any(non_fract_diff)){
	cat("WARNING:  non-fractional differences between original and remapped totals.\n");
	#print(abs_diff[non_fract_diff]);
}

###############################################################################

# Export summary table with IDs
outfn=paste(OutputFileRoot, ".ids.summary_table.tsv", sep="");
write_summary_file(out_count_mat, outfn)
#print(out_count_mat);

#-----------------------------------------------------------------------------

# Export summary table with names
outfn=paste(OutputFileRoot, ".names.summary_table.tsv", sep="");
colnames(out_count_mat)=name_map[colnames(out_count_mat)];
write_summary_file(out_count_mat, outfn)
#print(out_count_mat);

#-----------------------------------------------------------------------------

# Export summary stats
outfn=paste(OutputFileRoot, ".stats.tsv", sep="");
write.table(stat_mat, outfn, row.names=F, quote=F, sep="\t");

# Export what is in remaining
outfn=paste(OutputFileRoot, ".remaining_ids.tsv", sep="");
write.table(remaining_ids, outfn, row.names=F, col.names=F, quote=F);

###############################################################################

cat("\nDone.\n")
warns=warnings();
if(!is.null(warns)){
	print(warnings());
}
q(status=0)
