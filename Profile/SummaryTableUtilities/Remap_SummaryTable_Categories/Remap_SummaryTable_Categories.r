#!/usr/bin/env Rscript

###############################################################################

library('getopt');

source("~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r");
source("~/git/AnalysisTools/Metadata/OutputFileLibrary/OutputFileLibrary.r");

params=c(
	"summary_table", "s", 1, "character",
	"mapping_file", "m", 1, "character",
	"output_file_root", "o", 1, "character",
	"name_map_file", "n", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-s <input summary_table.tsv file>\n",
	"	-m <input mapping file>\n",
	"	-o <output file root>\n",
	"\n",
	"	[-n <input name map, e.g. <ID>\\t<name>\\n]",
	"\n",	
	"This script will read in a summary table and then for each\n",
	"of the categories, combine/remap them to a new category\n",
	"based on the mapping file.\n",
	"\n",
	"The input mapping file should have two columns and column names:\n",
	"	OutputIDs\\tInputIDs\\t...\\n\n",	
	"	<out_category>\\t<in_category>\\t...\\n\n",
	"\n",
	"\n",
	"The optional name map should have the columns:\n",
	"	OutputIDs\\tName\\n\n",
	"\n",
	"If the input category has ;'s in it, then the category will be split\n",
	"into a list and each item will be mapped to it's parent.  If the parents\n",
	"are different, the counts will be split.  If the mapping from child\n",
	"to parent is also one-to-many, then this is split too.\n",
	"\n");

if(
	!length(opt$summary_table) ||
	!length(opt$mapping_file) ||
	!length(opt$output_file_root)
){
	cat(usage);
	q(status=-1);
}


###############################################################################
# Parameters

SummaryTable=opt$summary_table;
MappingFile=opt$mapping_file;
OutputFileRoot=opt$output_file_root

NameMapFile="";
if(length(opt$name_map_file)){
	NameMapFile=opt$name_map_file;	
}

cat("\n")
cat("Summary Table: ", SummaryTable, "\n");
cat("Mapping File: ", MappingFile, "\n");
cat("Output File Root: ", OutputFileRoot, "\n");
cat("Name Mapping File: ", NameMapFile, "\n");
cat("\n");

###############################################################################
###############################################################################

load_map_from_matrix=function(mat){

	num_ids=nrow(mat);
	uniq_ids=sort(unique(mat[,2]));
	num_uniq=length(uniq_ids);

	# Allocating list
	mapping_list=vector("list", length=num_uniq);
	names(mapping_list)=uniq_ids;
	
	for(i in 1:num_ids){
		src=mat[i,2];
		dst=mat[i,1];
		#cat(src, "->", dst, "\n");

		mapping_list[[src]]=c(mapping_list[[src]], dst);
	}

	# clean up mapping by removing duplicated child to redundant parent mapping
	for(nm in names(mapping_list)){
		mapping_list[[nm]]=unique(mapping_list[[nm]]);
		if(is.null(mapping_list[[nm]])){
			cat("Error: ", src, " maps to NULL.\n");
			quit();
		}
	}

	results=list();
	results[["mapping"]]=mapping_list;
	results[["out_ids"]]=unique(mat[,1]);

	return(results);

}

load_name_map=function(fn){
	mat=as.matrix(read.table(fn, sep="\t", quote="", header=T));
	mapping_arr=gsub("[^a-zA-Z0-9_]", "_", mat[,2]);

	mapping_arr=c(mapping_arr, "Remaining");
	names(mapping_arr)=c(mat[,1], "Remaining");

	num_mappings=length(unique(mapping_arr));
	if(num_mappings<nrow(mat)){
		cat("Error: After making names variable name safe, they were no longer unique.\n");
		quit(status=-1);
	}

	return(mapping_arr);
}

###############################################################################
# Load files

counts_mat=load_summary_file(SummaryTable);
parent_child_map=as.matrix(read.table(MappingFile, sep="\t", quote="",  header=T));
mappings=load_map_from_matrix(parent_child_map);

child_to_parent_map=mappings[["mapping"]];
parent_ids=mappings[["out_ids"]];
num_parent_ids=length(parent_ids);

cat("Out IDs: [", num_parent_ids, "]\n", sep="");
print(parent_ids);
cat("\n");

if(NameMapFile!=""){
	name_map=load_name_map(NameMapFile);
}else{
	cat("Name Mape File not specified.\n");
}

###############################################################################
# Combines descendant categories together by parent

num_sumtab_cat=ncol(counts_mat);
num_sumtab_samples=nrow(counts_mat);
sumtab_cat=colnames(counts_mat);

sumtab_cat_splits=strsplit(sumtab_cat, ";");

out_count_mat=matrix(0.0, nrow=num_sumtab_samples, ncol=num_parent_ids+1);
rownames(out_count_mat)=rownames(counts_mat);
colnames(out_count_mat)=c(parent_ids, "Remaining");

remaining_ids=c();

for(cat_ix in 1:num_sumtab_cat){

	cur_cat_arr=sumtab_cat_splits[[cat_ix]];
	#cat("Categories splits: \n");
	#print(cur_cat_arr);
	
	# Collect parents across the category components
	cat_parents=c();
	for(cat_comp in cur_cat_arr){
		cat_parents=c(cat_parents, child_to_parent_map[[cat_comp]]);
	}

	# Get num uniq parents
	cat_parents=unique(cat_parents);
	num_parents=length(cat_parents);
	#cat("Parents:\n");
	#print(cat_parents);

	if(num_parents==0){
		# If parents are not found, then put it in Remaining
		out_count_mat[,"Remaining"]=
			out_count_mat[,"Remaining"]+counts_mat[,cat_ix];
		remaining_ids=c(remaining_ids, sumtab_cat[cat_ix]);
		
	}else{
		# Spread out counts across parents
		#print(counts_mat[,cat_ix]);
		adj_counts=counts_mat[,cat_ix]/num_parents;
		#print(adj_counts);
		for(par in cat_parents){
			out_count_mat[,par]=out_count_mat[,par] + adj_counts;
		}
	}

}

###############################################################################
# Remove zero count categories

cat_counts=apply(out_count_mat, 2, sum);
zero_count_cat_ix=(cat_counts==0);
nonzero_count_cat_ix=!zero_count_cat_ix;

out_cat=colnames(out_count_mat);
zero_count_parent_ids=out_cat[zero_count_cat_ix];
num_zc_par_ids=length(zero_count_parent_ids);

nonzero_count_parent_ids=out_cat[nonzero_count_cat_ix];
num_nzc_par_ids=length(nonzero_count_parent_ids);

cat("Zero count Parent IDs: [", num_zc_par_ids, "]\n", sep="");
print(zero_count_parent_ids);
cat("\n");
cat("Non-zero count Parent IDs: [", num_nzc_par_ids, "]\n", sep="");
print(nonzero_count_parent_ids);

out_count_mat=out_count_mat[,nonzero_count_parent_ids,drop=F];
#print(out_count_mat);

###############################################################################
# Generate stats across categories

norm_mat=normalize(out_count_mat);
median_norm=apply(norm_mat, 2, function(x){median(x,na.rm=T)});

sorted_median_norm=sort(median_norm, decreasing=T);
cat("\n");
cat("--------------------------------------------------------------\n");
cat("Top Median Abundances:\n");
stat_mat=matrix(sprintf("%0.4f", sorted_median_norm), ncol=1, 
	dimnames=list(names(sorted_median_norm), "Med_Abund"));
stat_cnames=colnames(stat_mat);
stat_mat=cbind(rownames(stat_mat), stat_mat);
colnames(stat_mat)=c("Category", stat_cnames);

if(NameMapFile!=""){
	Description=name_map[stat_mat[,"Category"]];
	stat_mat=cbind(stat_mat, Description);
}

print(stat_mat);
cat("--------------------------------------------------------------\n");
cat("\n\n");

###############################################################################
# Compare counts to confirm everything still adds up 

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
