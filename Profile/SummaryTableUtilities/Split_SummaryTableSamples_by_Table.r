#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 2, "character",
	"group_file", "p", 1, "character",
	"group_column", "r", 1, "numeric",
	"sample_column", "s", 1, "numeric" 
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
        "\nUsage:\n\n", script_name,
        "\n",
	"	-i <Input summary_table.tsv file>\n",
	"	[-o <Output filename root>]\n",
	"	-p <input grouP table tsv file>\n",
	"	-s <column of Sample identifiers in group table>\n",
	"	-r <Column of gRoup identifiers in group table>\n",
	"\n",
	"This script will read in a summary_table and a group table\n",
	"and then based on the groups specified, will generate new\n",
	"summary_table's by splitting the input summary_table.\n",
	"\n",
	"The sample names will be preserved.  The output files will\n",
	"contain the name of the grouping (i.e. the column name of the\n",
	"groups) and name of the group.\n",
	"\n",
	"For example a table file of:\n",
	"	sample_id \\t animals \n",
	"	sample1 \\t cat \n",
	"	sample2 \\t dog \n",
	"	sample3 \\t elephant \n",
	"	sample4 \\t cat \n",
	"	sample5 \\t na \n",
	"\n",
	"Will generate 3 output summary tables named:\n",
	"	<input summary_table root directory>/animals/cat/<filename root>.cat.summary_table.tsv\n",
	"	<input summary_table root directory>/animals/dog/<filename root>.dog.summary_table.tsv\n",
	"	<input summary_table root directory>/animals/elephant/<filename root>.elephant.summary_table.tsv\n",
	"\n");

if(!length(opt$input_file)
	|| !length(opt$group_file)
	|| !length(opt$group_column)
	|| !length(opt$sample_column)
){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileNameRoot=opt$output_file;
GroupFileName=opt$group_file;
GroupColumn=opt$group_column;
SampleColumn=opt$sample_column;

if(length(OutputFileNameRoot)==0){
	OutputFileNameRoot=opt$input_file;
	OutputFileNameRoot=gsub("\\.summary_table.tsv$", "", OutputFileNameRoot);
	OutputFileNameRoot=gsub("\\.summary_table.xls$", "", OutputFileNameRoot);
}

ofnr_components=strsplit(OutputFileNameRoot, "/")[[1]];
num_ofnr_comp=length(ofnr_components);

OutputRootDir=paste(head(ofnr_components, num_ofnr_comp-1), collapse="/");
OutputRootFname=ofnr_components[num_ofnr_comp];

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name Root: ", OutputFileNameRoot, "\n");       
cat("   Directory: ", OutputRootDir, "\n");
cat("   Filename: ", OutputRootFname, "\n");
cat("Group Table Name: ", GroupFileName, "\n");
cat("Group Column Number: ", GroupColumn, "\n");
cat("Sample Column Number: ", SampleColumn, "\n");

###############################################################################
###############################################################################

# Load data
inmat=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, row.names=1))
#cat("\nOriginal Matrix:\n")
#print(inmat);

# Grab columns we need into a vector, ignore totals, we won't trust it.
counts_mat=inmat[,2:(ncol(inmat))];
#print(counts_mat);

num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);
cat("\n");
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");
cat("\n");

###############################################################################

load_group_file=function(fname, sample_col, group_col){
	mapping_table=as.matrix(read.table(fname, sep="\t", check.names=FALSE, header=T));

	group_map=mapping_table[,group_col, drop=F];	
	rownames(group_map)=mapping_table[,sample_col];

	# Find samples without any group information
	uc_group_map=toupper(group_map);
	null_ix= uc_group_map=="NULL" | 
		 uc_group_map=="NA" |
		 uc_group_map=="." |
		 uc_group_map=="";
	cat("Samples without Group: \n");
	print(mapping_table[null_ix, sample_col]);
	cat("\n");

	# Remove samples without any group information
	null_samples=-which(null_ix);
	if(length(null_samples)>0){
		group_map=group_map[null_samples,,drop=F];
	}

	group_map=group_map[unique(rownames(group_map)),, drop=F];

	return(group_map);
}

###############################################################################

group_map=load_group_file(GroupFileName, SampleColumn, GroupColumn);
cat("Group Map (head only):\n");
print(head(group_map));
grouping_name=colnames(group_map);
cat("\n");

cat("Grouping Name: ", grouping_name, "\n");
dir_wgrpname=paste(OutputRootDir, "/", grouping_name, sep="");
cat("Making directory to store split summary tables: ", dir_wgrpname, "\n");
dir.create(dir_wgrpname);


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

###############################################################################

unique_groups=sort(unique(group_map[,1]));
cat("Unique Groups:\n");
num_groups=length(unique_groups);
print(unique_groups);
cat("\n");

sample_names=rownames(group_map);

for(i in 1:num_groups){
	cur_grp=unique_groups[i];

	# Extract sample members for group
	cat("Extracting samples for group: ", cur_grp, "\n");
	grp_ix=which(cur_grp==group_map[,1]);
	cur_grp_samples=sample_names[grp_ix];
	cat("\nGroup members (head only):\n");
	print(cur_grp_samples);
	cat("\n");

	# Extract counts for group
	missing_in_counts=setdiff(cur_grp_samples,rownames(counts_mat));
	if(length(missing_in_counts)>0){
		cat("ERROR: Missing sample(s) in summary table but requested in group: ", cur_grp, ".\n");
		print(missing_in_counts);
		cat("\n");
		quit(status=-1);
	}
	group_counts=counts_mat[cur_grp_samples,, drop=F];
	#print(group_counts);

	# Make directory to store summary table
	cur_grp=gsub(" ","_", cur_grp);
	dir_wgrp_wtype=paste(dir_wgrpname, "/", cur_grp, sep="");
	dir.create(dir_wgrp_wtype);
	
	# Write out file
	full_output_fname=paste(dir_wgrp_wtype, "/", OutputRootFname, ".", cur_grp, ".summary_table.tsv", sep="");
	cat("Writing Summary Table: ", full_output_fname, "\n");
	write_summary_file(group_counts, full_output_fname);

}

###############################################################################

cat("\nDone.\n");
if(!is.null(warnings())){
	print(warnings());
}

q(status=0);
