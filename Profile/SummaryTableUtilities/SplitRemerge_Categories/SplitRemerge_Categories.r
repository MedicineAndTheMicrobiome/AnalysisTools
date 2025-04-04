#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary_table.tsv file name>\n",
	"	-o <output file name summary table>\n",
	"\n",	
	"This script will read in the summary table and then for each\n",
	"of the categories that are a list of ID's separated by a semi-colon ';'\n",
	"split the counts into new categories.  When all the categories have\n",
	"beeen split, redundant categories are then collapsed into single\n",
	"\n",
	"Note: The output will be a summary table, thus:\n",
	"\n",
	"		 -o output.summary_table.tsv\n",
	"\n",
	"	will generate a file by the same name.\n",
	"\n",
	"\n");

if(!length(opt$input_file) || !length(opt$output_file)){
	cat(usage);
	q(status=-1);
}


###############################################################################

InputFileName=opt$input_file;
OutputFileName=opt$output_file;

cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("\n");

###############################################################################
###############################################################################

load_summary_table=function(summary_table_fn){
        # Load data
        cat("Loading Matrix (", summary_table_fn, ") ...\n", sep="");
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

        cat("\n");
        cat("Num Samples: ", num_samples, "\n");
        cat("Num Categories: ", num_categories, "\n");
        cat("\n");
        return(counts_mat);
}

###############################################################################
# Load counts matrix
counts_mat=load_summary_table(InputFileName);
num_categories=ncol(counts_mat);
num_samples=nrow(counts_mat);

# Parse category names and accumulate all the possible split IDs
split_list=list();
category_names=colnames(counts_mat);
categories_found=c();
for(i in 1:num_categories){
	
	cat("Category: [", i, "]: ", category_names[i], "\n", sep="");

	ids=strsplit(category_names[i], ";")[[1]];
	
	# clean up leading/trailing spaces
	for(id_ix in 1:length(ids)){
		ids[id_ix]=gsub("^\\s+", "", ids[id_ix]);
		ids[id_ix]=gsub("\\s+$", "", ids[id_ix]);
	}
	
	split_list[[i]]=ids;

	# Keep track of all see IDs
	categories_found=c(categories_found, ids);
	if(!(i%%50)){
		categories_found=unique(categories_found);
	}

}

unique_categories=sort(unique(categories_found));
num_unique_cat=length(unique_categories);

# Create output matrix
out_matrix=matrix(0, nrow=num_samples, ncol=num_unique_cat);
colnames(out_matrix)=unique_categories;
rownames(out_matrix)=rownames(counts_mat);

# Split out the categories and accumulate them
for(i in 1:num_categories){

	counts=counts_mat[,i];
	split_cat=split_list[[i]];
	num_splits=length(split_cat);

	# Divide the counts evenly by the number of ids
	split_counts=counts/num_splits;

	for(sp_ix in 1:num_splits){
		split_id=split_cat[sp_ix];
		acc=out_matrix[,split_id];
		acc=acc+split_counts;
		out_matrix[,split_id]=acc;
	}
}

#print(split_list);
#print(unique_categories);
#print(out_matrix);

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
write_summary_file(out_matrix, OutputFileName);

###############################################################################

cat("\nDone.\n")
warns=warnings();
if(!is.null(warns)){
	print(warnings());
}
q(status=0)
