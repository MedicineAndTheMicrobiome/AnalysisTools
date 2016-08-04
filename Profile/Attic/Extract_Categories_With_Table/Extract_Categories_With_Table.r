#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"table_file", "t", 1, "character",
	"category_column", "c", 1, "numeric",
	"group_column", "r", 1, "numeric",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	-t <table of categories to extract>\n",
	"	-c <category column number (starting from 1)>\n",
	"	-r <group column number (starting from 1)>\n",
	"	[-o <output reduced summary table file name root>]\n",
	"\n",	
	"Reads in table consisting of groups and category names then\n",
	"uses that information to extract and rebuild a summary table\n",
	"based on the entities in the specified column.\n",
	"\n",
	"If a category is in more than one group, then it's counts\n",
	"will be uniformly distributed across the groups that it belongs\n",
	"to.  If a category is in a group more than once, it's only\n",
	"counted as being in the group once.\n",
	"\n",
	"The output summary table will only contain the groups in the\n",
	"specified group column and a 'Remaining' category.\n",
	"If the script is functioning properly, the totals and remaining\n",
	"will remain constant independent of which groups you specify.\n",
	"\n",
	"\n", sep="");

if(!length(opt$input_file) || !length(opt$table_file) || 
	!length(opt$category_column) || !length(opt$group_column)){
	cat(usage);
	q(status=-1);
}else{
	InputFileName=opt$input_file;
	TableFileName=opt$table_file;
	CategoryColumn=opt$category_column;
	GroupColumn=opt$group_column;
}

if(!length(opt$output_file)){
	OutputFileName = gsub("\\.summary_table\\.xls$", "", opt$input_file);
	OutputFileName = gsub("\\.summary_table\\.tsv$", "", OutputFileName);
}else{
	OutputFileName = opt$output_file;
}

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("Group Table File Name: ", TableFileName, "\n");
cat("\n");

###############################################################################
###############################################################################

# Load data
inmat=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="*", row.names=1))
#cat("Original Matrix:\n")
#print(inmat);

# Grab columns we need into a vector, ignore totals, we won't trust it.
counts_mat=inmat[,2:(ncol(inmat)), drop=F];
#print(counts_mat);

num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);
category_names=colnames(counts_mat);
sample_names=rownames(counts_mat);

cat("\n");
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");
cat("\n");

###############################################################################

load_group_table=function(fn, grp_col, cat_col){

	inmat=as.matrix(read.delim(fn, sep="\t", header=TRUE, check.names=F));
	num_rows=nrow(inmat);
	groups=inmat[,grp_col];
	categories=inmat[,cat_col];

	group_map=cbind(groups, categories);

	# Infer groups when there is empty space
	if(grp_col!=cat_col){
		cur_grp="";
		for(i in 1:num_rows){
			if(group_map[i,1]!=""){
				cur_grp=group_map[i,1];
			}else{
				group_map[i,1]=cur_grp;
			}
		}
	}

	# Identify and remove invalid categories
	invalid_categories=which(is.na(categories) | categories=="");
	if(length(invalid_categories)>0){
		cat("Invalid categories found at: \n");
		invalid_table=(cbind(invalid_categories, as.data.frame(group_map[invalid_categories,])));
		colnames(invalid_table)=c("Row_Index", "Group_Name", "Category_Name");
		print(invalid_table);
		cat("\n\n");
		group_map=group_map[-invalid_categories,];
	}

	colnames(group_map)=colnames(inmat)[c(grp_col, cat_col)];

	return(group_map);
}

###############################################################################

grp_map=load_group_table(TableFileName, GroupColumn, CategoryColumn);

unique_groups=unique(grp_map[,1]);
num_unique_groups=length(unique_groups);

cat("Groups Identified (", num_unique_groups, ") :\n", sep="");
print(unique_groups);


unique_categories=unique(grp_map[,2]);
num_unique_categories=length(unique_categories);
multi_group_category=rep(0, num_unique_categories);
names(multi_group_category)=unique_categories;
for(i in 1:num_unique_groups){
	group_name=unique_groups[i];
	group_categories=unique(grp_map[which(grp_map[,1]==group_name), 2]);
	multi_group_category[group_categories]=multi_group_category[group_categories]+1;
}
print(multi_group_category[multi_group_category>1]);


collapsed_table=matrix(0, nrow=num_samples, ncol=num_unique_groups);

for(i in 1:num_unique_groups){

	group_name=unique_groups[i];
	grp_ix=which(group_name==grp_map[,1]);
	grp_cat=grp_map[grp_ix, 2];

	num_subcat=length(grp_cat);
	cat("Extracting group: ", group_name, "\n");
	cat("\tNum categories: ", num_subcat, "\n");
	
	found_categories=intersect(grp_cat, category_names);
	num_found_cat=length(found_categories);
	cat("\tNum found in summary table: ", num_found_cat, "\n");
	
	grp_counts=counts_mat[,found_categories, drop=F];
	#print(grp_counts);

	num_cat_in_grp=ncol(grp_counts);
	cat_names_in_grp=colnames(grp_counts);
	for(cat_name in cat_names_in_grp){
		grp_counts[,cat_name]=grp_counts[,cat_name]/multi_group_category[cat_name];	
	}
	#print(grp_counts);

	collapsed_table[,i]=apply(grp_counts, 1, sum);
	#print(grp_counts);
}


# Compute remaining category for leftover counts
sample_count_sums=apply(counts_mat, 1, sum);
collapsed_sums=apply(collapsed_table, 1, sum);
remaining=sample_count_sums-collapsed_sums;
collapsed_table=cbind(collapsed_table, remaining);
colnames(collapsed_table)=c(unique_groups, "Remaining");

# Remove categories with 0 count across all samples
category_sums=apply(collapsed_table, 2, sum);
empty_categories=which(category_sums==0);

if(length(empty_categories)>0){
	collapsed_table=collapsed_table[,-empty_categories];
}

###############################################################################

write_summary_table=function(counts_mat, stfname){

        cat("Writing summary file table: ", stfname, "\n", sep="");

        st_fh=file(stfname, "w");

        num_samples=nrow(counts_mat);
        sample_names=rownames(counts_mat);
        category_names=colnames(counts_mat);

        # Output Header
        cat(file=st_fh, paste("sample_id", "total", paste(category_names, collapse="\t"), sep="\t"));
        cat(file=st_fh, "\n");

        # Output rows
        for(samp_idx in 1:num_samples){
                total=sum(counts_mat[samp_idx,]);
                outline=paste(sample_names[samp_idx],total,paste(counts_mat[samp_idx,], collapse="\t"), sep="\t");
                cat(file=st_fh, outline, "\n", sep="");
        }

        close(st_fh);
        cat("Done writing.\n");
}

group_name=colnames(grp_map)[1];
group_name=gsub(" ", "_", group_name);
write_summary_table(collapsed_table, paste(OutputFileName, ".", group_name, ".summary_table.tsv", sep=""));


###############################################################################

cat("Done.\n")
if(!is.null(warnings())){
	print(warnings());
}

q(status=0)
