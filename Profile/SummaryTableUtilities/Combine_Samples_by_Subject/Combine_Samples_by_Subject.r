#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"mapping_table", "m", 1, "character",
	"output_fname", "o", 1, "character",
	"collapse_list", "l", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table>\n",
	"	-m <sample mapping table>\n",
	"	-o <output filename>\n",
	"	[-l <collapse list, comma separated>]\n",
	"\n",	
	"This script will read in both the summary table and mapping table\n",
	"The mapping table has distinct categories for each sample in columns\n",
	"For example, a pairs mapping.\n",
	"\n",
	"The sample mapping table should have the format:\n",
	"subject_id	smp_1	smp_2	...	smp_n\n",	
	"\n",
	"The subject ID can be any unique ID, but the smp should\n",
	"reference a sample ID in the input summary table\n",
	"\n",
	"The collapse list should be a list of columns names, eg. smp_1,smp_2\n",
	"\n",
	"\n");

if(!length(opt$input_file) || !length(opt$mapping_table) || !length(opt$output_fname)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
MappingTable=opt$mapping_table;
OutputFileName=opt$output_fname;

if(length(opt$collapse_list)){
	CollapseList=strsplit(opt$collapse_list,",")[[1]];
}else{
	CollapseList=c();
}

if(!length(grep(".summary_table.tsv$", OutputFileName))){
	OutputFileName=paste(OutputFileName, ".summary_table.tsv", sep="");
}

###############################################################################

cat("\n");
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("Mapping Table Name: ", MappingTable, "\n");
cat("Collapse Column List: \n");
print(CollapseList);
cat("\n");

###############################################################################

load_summary_table=function(summary_table_fn){
        # Load data
        cat("Loading Summary Table Counts (", summary_table_fn, ") ...\n", sep="");
        inmat=as.matrix(read.table(summary_table_fn, sep="\t", header=TRUE, 
		check.names=FALSE, row.names=1, quote=""))

        #cat("\nOriginal Matrix:\n")
        #print(inmat);

        # Grab columns we need into a vector, ignore totals, we won't trust it.
        counts_mat=inmat[,2:(ncol(inmat)), drop=F];
        #cat("\nCounts Matrix:\n");
        #print(counts_mat);

        num_samples=nrow(counts_mat);
        num_categories=ncol(counts_mat);
        sample_names=rownames(counts_mat);

        return(counts_mat);
}

write_summary_table=function(out_mat, fname){
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

load_mapping_file=function(fname){
	cat("Loading Mappings (", fname, ") ...\n", sep="");
	mapping_table=as.matrix(read.table(fname, sep="\t", check.names=FALSE, header=T, row.names=1))
	return(mapping_table);
}

###############################################################################

# Load Mapping
counts_mat=load_summary_table(InputFileName);

num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);

cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");
cat("\n");

###############################################################################

full_mapping_table=load_mapping_file(MappingTable);

num_subjects=nrow(full_mapping_table);
num_samp_types=ncol(full_mapping_table);
table_sample_types=colnames(full_mapping_table);

cat("Num Subjects: ", num_subjects, "\n");
cat("Num Sample Types: ", num_samp_types, "\n");
cat("\n");

###############################################################################

missing_types=setdiff(CollapseList, table_sample_types);
if(length(missing_types)){
	cat("Error, missing requested sample types.\n");
	print(missing_types);
}else{
	cat("Ok, all requested samples found.\n");
}

###############################################################################

cat("\n");
cat("Starting to collapse:\n");

map_sbjs=rownames(full_mapping_table);

combine_all=(length(CollapseList)==0);

combined_list=list();

for(sbj in map_sbjs){

	cat("Collapsing: ", sbj);

	if(combine_all){
		samp_ids_arr=full_mapping_table[sbj, ];
	}else{	
		samp_ids_arr=full_mapping_table[sbj, CollapseList];
	}

	# Extract out sample IDs that are not NA
	samp_ids_arr=samp_ids_arr[!is.na(samp_ids_arr)];
	num_samples_to_combine=length(samp_ids_arr);
	
	cat(" [", num_samples_to_combine, "]\n", sep="");
		
	if(num_samples_to_combine>0){

		# Construct a new sample ID based on subject id and samples that were included
		types=names(samp_ids_arr);
		comb_id=paste(sbj, "_", paste(types, collapse="."), sep="");

		# Combine the counts
		comb=apply(counts_mat[samp_ids_arr,,drop=F], 2, sum);

		# Save the counts into a list
		combined_list[[comb_id]]=comb;

	}
	
}

comb_ids=names(combined_list);
cat("\n");
cat("Combined Sample IDs:\n");
print(comb_ids);

num_combined_samples=length(combined_list);

out_matrix=matrix(NA, nrow=num_combined_samples, ncol=num_categories);
categories=colnames(counts_mat);
colnames(out_matrix)=categories;
rownames(out_matrix)=comb_ids;

for(cmb_id in comb_ids){
	out_matrix[cmb_id,categories]=combined_list[[cmb_id]][categories];
}

###############################################################################

write_summary_table(out_matrix, OutputFileName);

###############################################################################

cat("\nDone.\n");
if(!is.null(warnings())){
	print(warnings());
}

q(status=0);
