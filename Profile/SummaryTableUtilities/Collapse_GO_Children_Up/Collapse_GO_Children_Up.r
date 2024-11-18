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
OBOTable=opt$obo_table;
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

        mat=read.table(obo_map_fn, quote="", comment.char="", sep="\t",
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

load_obo_to_trees=function(obo_tab){

        cat("Loading OBO Table into (Parent-to-child) and (Child-to-Parent) Tree...\n");
        start_time=Sys.time();

        num_ids=nrow(obo_tab);

	# Initialize a empty list with obo_tab IDs and ""
        empty_tree=vector("list", length=num_ids);
	empty_tree[1:num_ids]="";
	names(empty_tree)=obo_tab[,"id"];

	parent_to_child_tree=empty_tree;
	child_to_parent_tree=empty_tree;

	# Split all the is_a (parent) strings ahead of time for efficiency
        parent_presplit=strsplit(obo_tab[,"is_a"], ";");
	
        for(i in 1:num_ids){
                child_id=obo_tab[i,"id"];
                parent_ids=parent_presplit[[i]];

		# Populate parent_to_child Tree
                for(pid in parent_ids){

			if(all(parent_to_child_tree[[pid]]=="")){
				parent_to_child_tree[[pid]]=child_id;
			}else{
				parent_to_child_tree[[pid]]=c(
					parent_to_child_tree[[pid]],
					child_id);
			}

                }

		# Populate child_to_parent Tree
		child_to_parent_tree[[child_id]]=parent_ids;

        }

	# End and compile timing
        end_time=Sys.time();
        cat("ok.\n");
        exec_time=end_time-start_time;
        cat("Time: ", exec_time, " secs\n", sep="")

	# Return obo_trees
	obo_trs=list();
	obo_trs[["p_to_c"]]=parent_to_child_tree;
	obo_trs[["c_to_p"]]=child_to_parent_tree;;

	#print(result);
        return(obo_trs);

}

#------------------------------------------------------------------------------

find_leaf_children=function(obo_trs){
	# Make a list of all the leaf nodes
	leaf_nodes_ix=unlist(lapply(obo_tree, function(x){x=="";}));
	leaf_names=names(obo_tree[leaf_nodes_ix]);
	return(leaf_names);	
}

#------------------------------------------------------------------------------

collapse_leaves=function(obo_trs, tgts){
	
	# For each target leaf, find it's parents, and remove itself. 

	p_to_c_tree=obo_trs[["p_to_c"]];
	c_to_p_tree=obo_trs[["c_to_p"]];

	for(tgt in tgts){

		cat("Removing: ", tgt, "\n", sep="");

		# Remove child from parent (p_to_c)
		parents_of_target=c_to_p_tree[[tgt]];
		for(prnt in parents_of_target){
			p2c_node=p_to_c_tree[[prnt]];
			p2c_node=setdiff(p2c_node, tgt);
			p_to_c_tree[[prnt]]=p2c_node;
			if(length(p_to_c_tree[[prnt]])==0){
				p_to_c_tree[[prnt]]="";
			}
		}
		p_to_c_tree[[tgt]]=NULL;

		# Remove child
		c_to_p_tree[[tgt]]=NULL;
	}

	p_to_c_tree=obo_trs[["p_to_c"]]=p_to_c_tree;
	c_to_p_tree=obo_trs[["c_to_p"]]=c_to_p_tree;
	return(obo_trs);
}

###############################################################################
###############################################################################

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

#------------------------------------------------------------------------------

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

#------------------------------------------------------------------------------

plot_zero_stats=function(zstats, tp_mat=NULL, cutoff=NULL){

	par(mfrow=c(1,1));

	num_cat=nrow(zstats);
	if(is.null(tp_mat)){
		cols=rep("black", num_cat);
	}else{
		cols=rep("grey", num_cat);
		excd_cut_ix=tp_mat[,"PropTrue"]>=cutoff;
		cols[excd_cut_ix]="blue";
		num_excd_cut=sum(excd_cut_ix);
	}

	plot(log10(zstats[,"NonZeroMean"]), log10(zstats[,"NumZeros"]),
		col=cols,
		main="Num Zeros vs. Mean Abund",
		xlab="log10(Mean Abd of NonZeros)", ylab="log10(Number of Zeros)");

	if(!is.null(cutoff)){
		title(main=paste("Min 'True' Cutoff = ", cutoff, ", Num = ", num_excd_cut), 
			cex.main=.9, col.main="blue", line=.5);
	}
}

#------------------------------------------------------------------------------

extract_lowPropTrue=function(tp_mat, cutoff){
	keep_ix=(tp_mat[,"PropTrue"]>=cutoff);
	out_tp=tp_mat[keep_ix,"PropTrue",drop=F];
	return(out_tp);
}

###############################################################################
###############################################################################

collapse_counts=function(cts_mat, child_to_parent_map, targets){

}



###############################################################################
# Load files

if(0){

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

tp_above_co=extract_lowPropTrue(true_pos_mat, cutoff);

print(tp_above_co);

quit();
}

###############################################################################

obo_map=load_obo_map(OBOTable);
trees=load_obo_to_trees(obo_map);

if(0){
	cat("\nOBO Map:\n");
	print(obo_map);

	cat("\nOBO Parent-to-Child Tree:\n");
	print(trees[["p_to_c"]]);

	cat("\nOBO Child-to_Parent Tree:\n");
	print(trees[["c_to_p"]]);

	cat("\n");
	#trees=collapse_leaves(trees, c("GO:N_2","GO:N_3","GO:N_4"));
	#trees=collapse_leaves(trees, c("GO:N_6"));
	#trees=collapse_leaves(trees, c("GO:N_14"));
	#print(trees);
}

quit();

###############################################################################
#
# 1.) Pull cat that can be collapsed
# 2.) Pull leaves from tree
# 3.) Intersect cat and leaves
# 4.) Merge leave with parent(s)
# 5.) Repeat until cat/leave intersect is empty

orig_counts_mat=counts_mat;
orig_norm_mat=norm_mat;
orig_zero_stats_mat=zero_stats_mat;
orig_true_pos_mat=true_pos_mat;

stop_collapsing=F;
while(!stop_collapsing){

	leaf_children_arr=find_leaf_children(obo_tree);

	low_prop_cat_mat=extract_lowPropTrue(true_pos_mat, TP_Cutoff);
	low_prop_cat=rownames(low_prop_cat_mat);

	collapse_targets=intersect(low_prop_cat, leaf_children_arr);
	num_collapse_targets=length(collapse_targets);

	if(num_collapse_targets>0){

		# Collapse OBO and Counts
		collapse_list=get_target_parents(obo_parent_lookup, collapse_targets);
		obo_tree=collapse_leaves(obo_parent_lookup, collapse_list);
		counts_mat=collapse_counts(counts_mat, obo_parent_lookup, collapse_list);
		norm_mat=normalize(counts_mat);

		# Refresh true_pos_mat
		zero_stats_mat=calc_zero_stats(counts_mat, norm_mat);
		true_pos_mat=estimate_TruePosNegs(zero_stats_mat, norm_mat, samp_depths);

	}else{
		stop_collapsing=T;
	}

}

###############################################################################

cat("\nDone.\n")
warns=warnings();
if(!is.null(warns)){
	print(warnings());
}
q(status=0)
