#!/usr/bin/env Rscript

###############################################################################

library('getopt');

source("~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r");
source("~/git/AnalysisTools/Metadata/OutputFileLibrary/OutputFileLibrary.r");

MIN_TRUE_NEG_PROB=.95;
MIN_TPN_PROP=.95;

params=c(
	"summary_table", "s", 1, "character",
	"obo_table", "b", 1, "character",
	"output_dir", "o", 1, "character",
	"false_neg_cutoff", "f", 2, "numeric",
	"min_tpn_cutoff", "t", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-s <input summary_table.tsv file>\n",
	"	-b <input OBO table>\n",
	"	-o <output directory>\n",
	"	[-f <True Negative Cutoff (probability)>, default=", MIN_TRUE_NEG_PROB, ">]\n",
	"	[-t <min proportion of True Positive/Negative (TPNs), default=", MIN_TPN_PROP, ">]\n",
	"\n",	
	"This script will read in the summary table and OBO table.\n",
	"\n",
	"Recursively, it will look at each parent ID, and collapse\n",
	"  the counts from the children.  If the parent's TPNs exceed\n",
	"  the threshold, then it's profile is exported.\n",
	"\n",
	"For each parent ID, the number of TN and TP will be calculated.\n",
	"If the proportion: (TP+PN)/NumSamp is >", MIN_TPN_PROP, "\n",
	"then the parent's summary table will be exported.\n",
	"\n",
	"The TP (true positives) are the number of non-zero data points\n",
	"The TN (true negatives) will be based on the mean abundance of the NZ\n",
	"  and the sequencing depth of that sample.\n",
	"\n",
	"----------------------------------------------------------------------\n",
	"\n",
	"The True Negative Cutoff is the probability a sample's category is\n",
	"  non-zero, given the read depth of the sample, and mean abundance of the\n",
	"  non-zero samples of that category.\n",
	"\n",
	"  In other words, if the probably TN Cutoff is .95, then 95% of the\n",
	"  time, the category should have a nonzero count, given the sample's depth\n",
	"  This means it's a True Negative, because it's zeroness is not due to\n",
	"  sequencing depth alone, but rather it's really 0 for this sample.\n",
	"\n",
	"The True Positive/Negative (TPN) cutoff is the minimum proportion of samples\n",
	"  that must exceeded, for the category to be exported/kept.\n",
	"\n",
	"\n");

if(
	!length(opt$summary_table) ||
	!length(opt$obo_table) ||
	!length(opt$output_dir)
){
	cat(usage);
	q(status=-1);
}

options(width=200, useFancyQuotes=F);

###############################################################################
# Parameters

SummaryTable=opt$summary_table;
OBOTable=opt$obo_table;
OutputDir=opt$output_dir;

if(length(opt$false_neg_cutoff)){
	MinTrueNegCutoff=opt$false_neg_cutoff;
}else{
	MinTrueNegCutoff=MIN_TRUE_NEG_PROB;
}

if(length(opt$min_tpn_cutoff)){
	MinTrueCutoff=opt$min_tpn_cutoff;
}else{
	MinTrueCutoff=MIN_TPN_PROP;
}

cat("\n")
cat("Summary Table: ", SummaryTable, "\n");
cat("OBO Table: ", OBOTable, "\n");
cat("Output Directory: ", OutputDir, "\n");
cat("Minimum True Cutoff: ", MinTrueCutoff, "\n");
cat("Minimum True Neg Cutoff: ", MinTrueNegCutoff, "\n");
cat("\n");

if(!dir.exists(OutputDir)){
	dir.create(OutputDir);
}else{

}

#pdf(paste(OutputFileRoot, ".clps_zero_chldn.pdf", sep=""), height=11, width=8.5);

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
	obo_trs[["c_wts"]]=sapply(child_to_parent_tree, function(x){1/length(x)});

	#print(result);
        return(obo_trs);

}

load_id_to_name_map=function(obo_tab){
	name_map=list();
	name_map=as.list(obo_tab[,"name"]);
	names(name_map)=obo_tab[,"id"];
	return(name_map);
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

estimate_TruePosNegs=function(zstats, nrm_mat, smp_dpt, true_neg_cutoff=.95){ 
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
		#cat("Analyzing: ", cur_cat, "\n");

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
		true_neg_ix=(prob_non_zero>true_neg_cutoff);
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
# Load files

if(1){

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

#print(counts_mat);

}

#------------------------------------------------------------------------------

if(1){

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

}

###############################################################################

obo_map=load_obo_map(OBOTable);

id_to_name_map=load_id_to_name_map(obo_map);

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

###############################################################################

find_root=function(c_to_p_tree){
	# This function will search through the child-to-parent tree
	# to find all the roots.  There could be mulitple roots.
	# The GO ontology has the roots of molecular function, biological process
	# and cellular component
	res=sapply(c_to_p_tree, function(x){length(x)==0});
	root_ids=names(which(res));
	return(root_ids);	
}

roots=find_root(trees[["c_to_p"]]);
cat("Found the Root(s):\n");
print(roots);

###############################################################################

get_vals=function(targ_id, mat){
	# Return counts if target is in matrix, else return 0's
	cat_names=colnames(mat);
	if(targ_id %in% cat_names){
		return(mat[,targ_id, drop=F]);
	}else{
		null_mat=matrix(0, nrow=nrow(mat), ncol=1);
		colnames(null_mat)=targ_id;
		rownames(null_mat)=rownames(mat);
		return(null_mat);
	}
}

remove_zeros_cols=function(zero_cat_name, mat){
	# Return a cleaned up matrix.  
	#  1. Removes all zero categories
	#  2. If all categories are zero, then use zero_cat_name as a placeholder
	
	num_col=ncol(mat);

	# Count categories
	nonzero_ix=apply(mat, 2, function(x){sum(x)>0;});

	if(sum(nonzero_ix)==0){
		# Return all zeros placeholder
		return(mat[,zero_cat_name,drop=F]);
	}else{
		# Return cleaned up matrix
		return(mat[,nonzero_ix,drop=F]);
	}
	
		
}

accumulate_st_by_id=function(par_st_list, parent_id, p_to_c_tree, child_wts, st_mat){
	# parent_st_list: list (hash) of 2D summary tables accumulate from child up to the
	#	the specified parent_id
	#
	# parent_id: id to start searching from.  If "direct" called, the root_id of the
	#	tree is specified. 
	#
	# p_to_c_tree: the tree, mapping each parent to all it's children
	#
	# child_wts: sometimes a child can have multiple parents, so the weights indicate
	#	the proportion of its counts that should go to each parent.
	#
	# st_mat: the original, accummulated, summary table for each id

	cat("Visiting: ", parent_id, "\n");
	children_ids=p_to_c_tree[[parent_id]];

	avail_st_cat=colnames(st_mat);

	# Base case, if leaf node, return children counts
	num_children=length(children_ids);
	if(num_children==1 && children_ids==""){
		# Leaf node
		cat("\tLeaf node\n");
		kid_arr=get_vals(parent_id, st_mat);
		par_st_list[[parent_id]]=kid_arr;
		return(par_st_list);	
	}

	# Recursive case, if parent node, accumulate children counts and totals

	# Allocate summary table for this parent_id
	children_matrix=matrix(0, nrow=nrow(st_mat), ncol=num_children);
	colnames(children_matrix)=children_ids;
	rownames(children_matrix)=rownames(st_mat);

	# For each child of the parent
	for(kid in children_ids){

		# Accumulate kid's st, update parent_st_list since it acts like
		#   a cache
		par_st_list=accumulate_st_by_id(
			par_st_list, kid, p_to_c_tree, child_wts, st_mat);

		# Pull out this kid's st
		kid_matrix=par_st_list[[kid]];

		# Collapse the kids st into a single array and weight it
		children_matrix[,kid]=apply(kid_matrix, 1, sum)*child_wts[kid];
	}

	# Add self to children matrix table
	self_counts=get_vals(parent_id, st_mat);
	self_and_kids=cbind(self_counts, children_matrix);

	# Remove zero count children.  If no categories left, then have all-zero parents ID.
	self_and_kids=remove_zeros_cols(parent_id, self_and_kids);

	cat("Saving: ", parent_id, "\n");
	par_st_list[[parent_id]]=self_and_kids;

	return(par_st_list);	
}

#------------------------------------------------------------------------------

calculate_tpn_prop=function(props, depths, false_neg_cutoff=0.05){

	#cat("Proportions:\n");
	#print(props);
	#cat("Depths:\n");
	#print(depths);
	zero_ix=(props==0);
	nzero_ix=!zero_ix;

	zero_depths=depths[zero_ix];
	nzero_depths=depths[nzero_ix];

	num_samples=length(props);

	nz_mean=mean(props[nzero_ix]);
	cat("Non-zero Mean: ", nz_mean, "\n");
	
	num_zeros=sum(zero_ix);
	cat("Num zeros: ", num_zeros, "\n");
	if(num_zeros==num_samples || num_zeros==0){
                num_true_neg=0;
	}else{
		calc_prob_nz=function(depth){
			# Probability of being nonzero given depth and mean abundance
			probnz=1-pbinom(0, depth, nz_mean);
			return(probnz);
		}
		prob_non_zero_of_zeros=sapply(zero_depths, calc_prob_nz);
		# If prob of nonzero > cutoff, call it TRUE NEG, because the
		#    the sequencing depth should've been deep enough.
		true_neg_ix=(prob_non_zero_of_zeros>false_neg_cutoff);
		num_true_neg=sum(true_neg_ix);
	}
	

	# Calc number of nonzero categories
	num_true_pos=sum(nzero_ix);
		
	# Proportion of samples that we are "sure" that are TRUE POS+NEG
	prop_true_np=(num_true_neg+num_true_pos)/num_samples;

	cat("Num True Neg: ", num_true_neg, "\n");
	cat("Num True Pos: ", num_true_pos, "\n");
	cat("Prop TNP: ", prop_true_np, "\n");

	return(prop_true_np);

}

calculate_tpn_rates=function(par_st_list, smp_dpt, tpn_cutoff){

	parent_ids=names(par_st_list);
	tnp_out=numeric(length(parent_ids));
	names(tnp_out)=parent_ids;

	for(pid in parent_ids){

		cat("Calculating for: ", pid, "\n");
		proportion_arr=apply(par_st_list[[pid]], 1, sum);	

		if(!all(names(proportion_arr)==names(smp_dpt))){
			cat("Error:  Order of Proportion Array doesn't Match Sample Depths\n");	
		}

		tnp_out[pid]=calculate_tpn_prop(proportion_arr, smp_dpt, false_neg_cutoff=tpn_cutoff);

		cat("\n\n");

	}

	return(tnp_out);

}

#------------------------------------------------------------------------------

export_sumtabs_exc_tpn=function(tpn_rt_lst, cnt_mat_list, smpdpt, tpn_thres, outdir, id2nam_map){

	rename_matrix=function(inmat, map){
		catnames=colnames(inmat);
		outnames=sapply(catnames, function(x){ id2nam_map[[x]]});
		outnames=sapply(outnames, function(x){ make.names(x)});
		outnames=sapply(outnames, function(x){ gsub("\\.","_",x)});
		colnames(inmat)=outnames;
		return(inmat);
	}

	append_remaining=function(inmat, depth){
		sums=apply(inmat, 1, sum);
		Remaining=depth-sums;
		if(any(Remaining!=0)){
			outmat=cbind(inmat,Remaining);
		}else{
			outmat=inmat;
		}
		return(outmat);
	}


	cat("Exporting tables above True PN cutoff: ", tpn_thres, "\n");

	parent_ids=names(tpn_rt_lst);
	for(par_id in parent_ids){
		cat("Parent ID: ", par_id, "\n");
		rate=tpn_rt_lst[par_id];
		cat(" Rate: ", rate, "\n");

		if(rate>tpn_thres){
			parent_name=gsub("\\.", "_", make.names(id2nam_map[[par_id]]));
			cat(" Name: ", parent_name, "\n");
			
			outfn=paste(outdir, "/", par_id, ".", parent_name, 
				".tpn", sprintf("%g", round(100*rate)),
				".summary_table.tsv", sep="");

			outfn=gsub("[:]","_",outfn);
			cat(" Summary Tab Filename: \"", outfn, "\"\n", sep="");

			out_mat=cnt_mat_list[[par_id]];

			renamed_mat=rename_matrix(out_mat, id2nam_map);

			renamed_mat=append_remaining(renamed_mat, smpdpt);

			write_summary_file(renamed_mat, outfn);			
				
			cat("\n\n");
		}
	}

}

###############################################################################

parent_st_norm_list=list();
parent_st_cnts_list=list();

for(root_id in roots){

	cat("Working on: ", root_id, "\n");


	cat("Accumulating Normalized Matrices:\n");
	parent_st_norm_list=accumulate_st_by_id(
			parent_st_norm_list, root_id, trees[["p_to_c"]], trees[["c_wts"]], norm_mat);

	cat("Accumulating Count Matrices:\n");
	parent_st_cnts_list=accumulate_st_by_id(
			parent_st_cnts_list, root_id, trees[["p_to_c"]], trees[["c_wts"]], counts_mat);
	

	cat("Calculating TPN Rates, Pr[True Negative > ", MinTrueNegCutoff, "]\n", sep="");
	tnp_rate_list=calculate_tpn_rates(parent_st_norm_list, samp_depths, MinTrueNegCutoff);

	full_name=id_to_name_map[[root_id]];

	export_sumtabs_exc_tpn(
		tnp_rate_list, parent_st_cnts_list, samp_depths,
		MinTrueCutoff, OutputDir, id_to_name_map);

	
}

quit();

###############################################################################

cat("\nDone.\n")
warns=warnings();
if(!is.null(warns)){
	print(warnings());
}
q(status=0)
