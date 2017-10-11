#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"mapping_table", "m", 1, "character",
	"output_fname_root", "o", 1, "character",
	"force_ignore", "f", 2, "logical",
	"column_num", "c", 2, "numeric",
	"output_dir", "d", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table>\n",
	"	-m <sample mapping table>\n",
	"	[-f (force ignore mismatching samples)]\n",
	"	[-o <output filename root>]\n",
	"	[-c <column number starting from 1, excluding sample IDs, default=1>]\n",
	"	[-d <output directory>]\n",
	"\n",	
	"This script will read in both the summary table and mapping table\n",
	"and merge the counts between samples that map to the same group/id.\n",
	"This script should be used to merge replicates and remove unwanted samples.\n",
	"\n",
	"For example the mapping table:\n",
	"	Sample \\t NewID\n",
	"	sample1.1 \\t sample1\n",
	"	sample1.2 \\t sample1\n",
	"	cat \\t animal\n",
	"	dog \\t animal\n",
	"	camel \\t animal\n",
	"	junk \\t na\n",
	"\n",
	"sample1.1 and sample1.2 will be merged into sample1 and the samples\n",
	"cat, dog and camel will be merged into a new sample called animal.\n",
	"If you don't want a sample included in the new output summary table,\n",
	"you put \"na\" in the group id column.  The above mapping table will\n",
	"produce exactly 2 samples, named sample1 and animal in the new summary table.\n",
	"If you put the same group id as the sample id, then the sample id will\n",
	"be preserved in the new summary table.\n",
	"\n",
	"The output file will have the format:\n",
	"	 <output filename root>.<NewID>.summary_table.tsv\n",
	"\n",
	"The force ignore option will allow the program to continue without\n",
	"aborting if the samples int he mapping table and the summary table\n",
	"do not match up.  It's safer to allow the program to abort so you\n",
	"can try to fix the error.  If you chose to use the force ignore,\n",
	"you should make sure your output looks correct.\n",
	"\n",
	"Based the samples included in the new summary table, a new mapping\n",
	"table will be generated.  If the collapsed samples had different\n",
	"metadata entries, they will be combined into a single field separated\n",
	"with semicolons.\n",
	"\n",
	"For example:\n",
	"	Sample	NewID	Metadata1	Metadata2\n",
	"	abc.1	abc	1		red\n",
	"	abc.2	abc	2		red\n",
	"	abc.3	abc	3		blue\n",
	"	def.1	def	1		green\n",
	"	def.2	def	2		green\n",
	"\n",
	"After collapsing by NewID Will look like:\n",
	"	NewID	Metadata1	Metadata2\n",	
	"	abc	1;2;3		red;blue\n",
	"	def	1;2		green\n",
	"\n",
	"The -d output directory option will override directory specified in the input\n",
	"or output file name root, if specified.\n",
	"\n");

if(!length(opt$input_file) || !length(opt$mapping_table)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
MappingTable=opt$mapping_table;

if(length(opt$output_fname_root)){
	OutputFileNameRoot=opt$output_fname_root;
}else{
	OutputFileNameRoot=gsub("\\.summary_table.tsv$", "", opt$input_file);
	OutputFileNameRoot=gsub("\\.summary_table.xls$", "", OutputFileNameRoot);
}

ColumnNum=1;
if(length(opt$column_num)){
	ColumnNum=opt$column_num;	
}

ForceIgnore=F;
if(length(opt$force_ignore)){
	ForceIgnore=T;
}

OutputDirectory="./";
# If output directory is specified, then override/truncate path of the OutputFileNameRoot

ofnr_components=strsplit(OutputFileNameRoot, "/")[[1]];
num_ofnr_comp=length(ofnr_components);
if(length(opt$output_dir)){
	# If output dir specified:
	OutputDirectory=opt$output_dir;
	OutputFileNameRoot=ofnr_components[num_ofnr_comp];
}else{
	# If output dir not specified
	OutputDirectory=paste(head(ofnr_components, num_ofnr_comp-1), collapse="/");
	OutputFileNameRoot=ofnr_components[num_ofnr_comp];
}

if(OutputDirectory==""){
	OutputDirectory="./";
}

if(!file.exists(OutputDirectory)){
	cat("Creating output directory: ", OutputDirectory, "\n");
	dir.create(OutputDirectory);
}

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Override Output Directory: ", OutputDirectory, "\n");
cat("Output File Name Root: ", OutputFileNameRoot, "\n");       
cat("Mapping Table Name: ", MappingTable, "\n");
cat("Column Number: ", ColumnNum, "\n");
cat("Force Ignore : ", ForceIgnore, "\n");

###############################################################################
###############################################################################

# Load data
inmat=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, row.names=1, quote=NULL))
#cat("Original Matrix:\n")
#print(inmat);

# Grab columns we need into a vector, ignore totals, we won't trust it.
counts_mat=inmat[,2:(ncol(inmat)), drop=F];

num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);
cat("\n");
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");
cat("\n");

###############################################################################

load_mapping_file=function(fname){
	mapping_table=as.matrix(read.table(fname, sep="\t", check.names=FALSE, header=T, row.names=1))
	return(mapping_table);
}

###############################################################################

full_mapping_table=load_mapping_file(MappingTable);
mapping_table=full_mapping_table[, ColumnNum, drop=F];

inmap_samples=rownames(inmat);
mapping_table_samples=rownames(mapping_table);
shared_samples=intersect(inmap_samples, mapping_table_samples);

target_column_name=colnames(mapping_table);
# Get target column's name, and create a directory to put the results in it
cat("Targeted Column's Name is: ", target_column_name, "\n", sep="");
OutputDirectory=paste(OutputDirectory, "/", target_column_name, sep="");
dir.create(OutputDirectory);

if(length(setdiff(shared_samples, mapping_table_samples))){
	cat("WARNING: Missing samples in Mapping Table.\n");
}

if(length(setdiff(shared_samples, inmap_samples))){
	cat("WARNING: Missing samples in Summary Table.\n");
}

if(ForceIgnore){
	inmat=inmat[shared_samples,];
	mapping_table=mapping_table[shared_samples, , drop=F];
}

print(mapping_table);

# Remove samples we want to exclude
unique_groups=sort(unique(mapping_table[,1]));
unusable_ix= is.na(unique_groups) | is.null(unique_groups) | (unique_groups=="") | 
	(unique_groups=="na") | (unique_groups=="NA");
unique_groups=unique_groups[!unusable_ix];

num_uniq_grps=length(unique_groups);
cat("\n");
cat("Num unique groups: ", num_uniq_grps, "\n");
cat("Unique Groups:\n");
print(unique_groups);
cat("\n");

###############################################################################
# Extract samples by group and join them if necessary

new_summary_table=matrix(NA, nrow=num_uniq_grps, ncol=num_categories);

cat("Working on extracting and merging:\n");
mapping_sample_names=rownames(mapping_table);
available_sample_names=rownames(counts_mat);

for(i in 1:num_uniq_grps){
	cur_grp=unique_groups[i];
	
	grp_idx=(mapping_table[,1]==cur_grp);
	grp_sample_names=mapping_sample_names[grp_idx];

	cat("Group: ", cur_grp);
	grp_sample_names=grp_sample_names[!is.na(grp_sample_names)];

	#print(grp_sample_names);
	num_to_collapse=length(grp_sample_names);
	cat(" (", num_to_collapse, ")\n", sep="");
	
	#print(available_sample_names);
	target_sample_names=intersect(available_sample_names, grp_sample_names);
	if(length(target_sample_names) != length(grp_sample_names)){
		cat("\n");
		cat("***************************************************************\n");
		cat("* WARNING: Targeted Sample Names don't match Available!!!     *\n");
		cat("* Targeted:                                                   *\n");
		print(grp_sample_names);
		#cat("* Available:                                                  *\n");
		#print(available_sample_names);
		cat("* Providing:                                                  *\n");
		print(target_sample_names);
		cat("***************************************************************\n");
		cat("\n");
	}	

	sub_matrix=counts_mat[target_sample_names, , drop=F];
	new_summary_table[i,]=apply(sub_matrix, 2, sum);
}

colnames(new_summary_table)=colnames(counts_mat);
rownames(new_summary_table)=unique_groups;

# Remove samples without counts
new_sample_counts=apply(new_summary_table, 1, sum);
non_zero=new_sample_counts>0;
if(any(!non_zero)){
	cat("********************************************************************\n");
	cat("* WARNING: Samples have been removed due to zero counts!!!         *\n");
	print(unique_groups[!non_zero])
	cat("********************************************************************\n");
}
new_summary_table=new_summary_table[non_zero,];

###############################################################################

write_summary_file=function(out_mat, fname){
        fc=file(fname, "w");
        cat(file=fc, paste("sample_id\ttotal", paste(colnames(out_mat), collapse="\t"), sep="\t"));
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
group_name=colnames(mapping_table)[1];
output_fname=paste(OutputDirectory, "/", OutputFileNameRoot, ".", group_name, ".summary_table.tsv", sep="");
write_summary_file(new_summary_table, output_fname);

###############################################################################
###############################################################################

# Generate new mapping table with rows collapsed and keyed by new sample ID
#print(full_mapping_table);
collapse_entries=function(in_mat){
	num_col=ncol(in_mat);
	if(num_col==0){
		return(in_mat);
	}
	collapsed=character(num_col);
	for(i in 1:num_col){
		col_val=gsub("^ *", "", in_mat[,i, drop=F]);
		collapsed[i]=paste(unique(col_val), collapse=";");
	}
	return(collapsed);
}

collpsed_map_table=matrix("", nrow=num_uniq_grps, ncol=ncol(full_mapping_table)-1);
rownames(collpsed_map_table)=unique_groups;
colnames(collpsed_map_table)=colnames(full_mapping_table)[-ColumnNum];
key_cat=colnames(full_mapping_table)[ColumnNum];
for(i in 1:num_uniq_grps){
	cur_grp=unique_groups[i];
        cat("Group: ", cur_grp, "\n");
	grp_idx=(full_mapping_table[,ColumnNum]==cur_grp);
	collpsed_map_table[i,]=collapse_entries(full_mapping_table[grp_idx,-ColumnNum,drop=F]);
}

fc=file(paste(OutputDirectory,"/",OutputFileNameRoot, ".", group_name, ".meta.tsv", sep=""));
cat(file=fc, paste(c(key_cat, colnames(collpsed_map_table)), collapse="\t"), "\n", sep="");
write.table(collpsed_map_table, file=paste(OutputDirectory,"/",OutputFileNameRoot, ".", group_name, ".meta.tsv", sep=""),
		quote=F, sep="\t", row.names=T, col.names=F, append=T);


###############################################################################

cat("\nDone.\n");
if(!is.null(warnings())){
	print(warnings());
}

q(status=0);
