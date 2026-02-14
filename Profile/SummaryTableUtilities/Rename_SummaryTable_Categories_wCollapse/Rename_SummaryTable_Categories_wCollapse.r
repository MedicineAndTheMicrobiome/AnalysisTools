#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"mapping_table", "m", 1, "character",
	"output_file", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table>\n",
	"	-m <sample mapping table>\n",
	"	-o <output summary table>\n",
	"\n",	
	"This script will read in both the summary table and mapping table\n",
	"and merge the counts between categories that map to the same group/id.\n",
	"The mapping file is assumed to not have headers.\n",
	"The first column is the original ID the second column is the new/destination ID.\n",
	"\n",
	"For example the mapping table:\n",
	"	sample1.1 \\t sample1\n",
	"	sample1.2 \\t sample1\n",
	"	cat \\t animal\n",
	"	dog \\t animal\n",
	"	camel \\t animal\n",
	"	junk \\t na\n",
	"\n",
	"If two categories have the same destination ID, then they will be collapsed together.\n",
	"If there is no mapping for a category, it will be left the same.\n",
	"\n",
	"\n");

if(!length(opt$input_file) || !length(opt$mapping_table)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
MappingTable=opt$mapping_table;
OutputFileName=opt$output_file;

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("Mapping Table Name: ", MappingTable, "\n");

###############################################################################
###############################################################################

# Load data
inmat=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, row.names=1, quote=NULL))
#cat("Original Matrix:\n")

# Grab columns we need into a vector, ignore totals, we won't trust it.
counts_mat=inmat[,2:(ncol(inmat)), drop=F];

num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);
cat("\n");
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");
cat("\n");

sample_ids=rownames(counts_mat);
category_ids=colnames(counts_mat);

###############################################################################

cat("Example Sample IDs:\n");
print(sample_ids[1:10]);
cat("\nExample Category IDs:\n");
print(category_ids[1:10]);
cat("\n");

###############################################################################
###############################################################################

load_mapping_file=function(fname){
	mapping_table=as.matrix(read.table(fname, sep="\t", check.names=FALSE, header=F, row.names=NULL, comment="#"))
	colnames(mapping_table)=c();

	num_rows=nrow(mapping_table);
	
	# Run this explicitly so items later in the list can be overridden
	map=list();
	for(i in 1:num_rows){
		old=mapping_table[i,1];
		new=mapping_table[i,2];
		map[[old]]=new;
	}
	
	return(map);
}

###############################################################################

id_map=load_mapping_file(MappingTable);

cat("Example Mapping Table contents:\n");
print(id_map[1:10]);

###############################################################################

new_categories=character(num_categories);
for(i in 1:num_categories){
	new_categories[i]=id_map[[category_ids[i]]];
}

num_remapped_categories=length(new_categories);
unique_remapped_categories=unique(new_categories);
num_unique_remapped_categories=length(unique_remapped_categories);

cat("Num Remapped:", num_remapped_categories, "\n");
cat("Num Unique Remapped: ", num_unique_remapped_categories, "\n");

out_mat=matrix(0, nrow=num_samples, ncol=num_unique_remapped_categories);
rownames(out_mat)=sample_ids;
colnames(out_mat)=unique_remapped_categories;

collapsed=list();

for(i in 1:num_categories){

	old_cat=category_ids[i];
	new_cat=id_map[[old_cat]];

	if(!is.null(collapsed[[new_cat]])){
		cat("Collapsing: ", old_cat, "/", new_cat, "\n");
	}

	out_mat[sample_ids, new_cat]=out_mat[sample_ids, new_cat] + counts_mat[sample_ids, old_cat];
	
	collapsed[[new_cat]]=new_cat;

}

###############################################################################

write_summary_file=function(out_mat, fname, sample_id_cname){
        fc=file(fname, "w");
        cat(file=fc, paste(sample_id_cname, "total", paste(colnames(out_mat), collapse="\t"), sep="\t"));
        cat(file=fc, "\n");
        sample_names=rownames(out_mat);
        num_samples=nrow(out_mat);
        for(samp_idx in 1:num_samples){
                total=sum(out_mat[samp_idx,]);
                outline=paste(sample_names[samp_idx], total,
                        paste(out_mat[samp_idx,], collapse="\t"), sep="\t");
                cat(file=fc, outline);
                cat(file=fc, "\n");
        }
        close(fc);
}

# Output new summary table
sample_id_cname=strsplit(readLines(InputFileName, n=1), "\t")[[1]][1];
write_summary_file(out_mat, OutputFileName, sample_id_cname);

###############################################################################

cat("\nDone.\n");
if(!is.null(warnings())){
	print(warnings());
}

q(status=0);
