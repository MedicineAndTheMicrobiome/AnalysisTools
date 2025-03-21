#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"mapping_table", "m", 1, "character",
	"output_fname_root", "o", 1, "character",
	"force_ignore", "f", 2, "logical",
	"column_name", "c", 2, "character",
	"output_dir", "d", 2, "character",
	"max_depth", "M", 2, "numeric"
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
	"	[-c <column name of new ID in sample mapping table>]\n",
	"	[-d <output directory>]\n",
	"	[-M <max depth used in merging samples, default=Inf>]\n",
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
	"\n",
	"The -M max depth limits the contribution of a sample's reads.\n",
	"  This allows individual samples to be represented more uniformly, independent of\n",
	"  sequencing depth.\n",
	"  For example: \n",
	"    if max_depth = 0, then the sample with the lowest depth will be used as the max_depth.\n",
	"    if max_depth = Inf, then there is no limit and all samples depths will be summed up.\n",
	"    if max_depth = 3000, if there are samples with <3000 reads, then any sample with >=3000\n",
	"        reads will be 'subsampled' down to 3000.  Otherwise, the depth of the sample with\n",
	"        the least depth will be used.\n",
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

ColumnName="";
if(length(opt$column_name)){
	ColumnName=opt$column_name;	
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

if(length(opt$max_depth)){
	MaxDepth=opt$max_depth;
}else{
	MaxDepth=Inf;
}

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Override Output Directory: ", OutputDirectory, "\n");
cat("Output File Name Root: ", OutputFileNameRoot, "\n");       
cat("Mapping Table Name: ", MappingTable, "\n");
cat("Column Name: ", ColumnName, "\n");
cat("Force Ignore : ", ForceIgnore, "\n");
cat("Max Depth : ", MaxDepth, "\n");

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

if(ColumnName==""){
	ColumnName=colnames(full_mapping_table[,1,drop=F]);
}
mapping_table=full_mapping_table[, ColumnName, drop=F];
nona_ix=!is.na(mapping_table[,1]);
mapping_table=mapping_table[nona_ix,,drop=F];
full_mapping_table=full_mapping_table[nona_ix,,drop=F];

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

combine_samples=function(count_mat, totals, max_depth){

	cat("Totals:\n");
	print(totals);
	cat("Max Depth: ", max_depth, "\n");

	num_samples=nrow(count_mat);

	if(max_depth==Inf || num_samples==1){
		combined_counts=apply(count_mat, 2, sum);
	}else{
		norm=matrix(NA, nrow=nrow(count_mat), ncol=ncol(count_mat));

		# Calculate how many reads to retain per sample
		kept_totals=numeric();

		# If all samples have more depth then the max_depth,
		# then use the min across all samples.
		min_cutoff=max(min(totals), max_depth);

		for(i in 1:num_samples){
			
			# Normalize all samples
			norm[i,]=count_mat[i,]/totals[i];

			# Calculate max depth from each sample to use
			kept_totals[i]=min(totals[i], min_cutoff);
		}

		cat("\nOriginal / Kept Totals:\n");
		print(totals);
		print(kept_totals);

		#cat("Normalized: \n");
		#print(norm);
		
		total_depth=sum(kept_totals);
		cat("Kept Depth: ", total_depth, "\n");
		kept_prop=kept_totals/total_depth;
		cat("Kept Proportions: \n");
		print(kept_prop);

		sum_weighted_normalized=rep(0, ncol(count_mat));

		for(i in 1:num_samples){
			sum_weighted_normalized=sum_weighted_normalized+
				kept_prop[i]*norm[i,];
		}

		#print(sum(sum_weighted_normalized));
		combined_counts=round(sum_weighted_normalized*total_depth);

		cat("\nCombined Depth (with rounding):\n");
		print(sum(combined_counts));

	}

	return(combined_counts);
}

###############################################################################
# Extract samples by group and join them if necessary

new_summary_table=matrix(0, nrow=num_uniq_grps, ncol=num_categories);

cat("Working on extracting and merging:\n");
mapping_sample_names=rownames(mapping_table);
available_sample_names=rownames(counts_mat);

for(i in 1:num_uniq_grps){

	cat("\n\n-----------------------------------------------------------------------------\n");

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

	if(length(target_sample_names)>0){

		sub_matrix=counts_mat[target_sample_names, , drop=F];
		sums=apply(sub_matrix, 1, sum);

		if(MaxDepth==0){
			used_max_depth=min(sums);
			cat("Using depth of: ", used_max_depth, "\n");
		}else{
			used_max_depth=MaxDepth;
		}

		new_summary_table[i,]=combine_samples(sub_matrix, sums, used_max_depth);

	}
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

col_num=which(colnames(full_mapping_table)==ColumnName);
colnames(collpsed_map_table)=colnames(full_mapping_table)[-col_num];

for(i in 1:num_uniq_grps){
	cur_grp=unique_groups[i];
        cat("Group: ", cur_grp, "\n");
	grp_idx=(full_mapping_table[,ColumnName]==cur_grp);
	collpsed_map_table[i,]=collapse_entries(full_mapping_table[grp_idx,-col_num,drop=F]);
}

out_meta_filename=paste(OutputDirectory,"/",OutputFileNameRoot, ".", group_name, ".meta.tsv", sep="");

fc=file(out_meta_filename, "w");
cat(file=fc, paste(c(ColumnName, colnames(collpsed_map_table)), collapse="\t"), "\n", sep="");

if(dim(collpsed_map_table)[2]==0){
	write.table(rownames(collpsed_map_table), 
		file=out_meta_filename,
		quote=F,
		row.names=F, col.names=F, append=T)
}else{
	write.table(collpsed_map_table, 
		file=out_meta_filename,
		quote=F, sep="\t", 
		row.names=T, col.names=F, append=T);
}


###############################################################################

cat("\nDone.\n");
if(!is.null(warnings())){
	print(warnings());
}

q(status=0);
