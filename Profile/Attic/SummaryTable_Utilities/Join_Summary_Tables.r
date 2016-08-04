#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file_list", "i", 2, "character",
	"input_file_list_file", "l", 2, "character",
	"output_file", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	[-i <input comma-separated list of summary_table1.tsv,summary_table2.tsv,summary_tableN.tsv,...>]\n",
	"	 or\n",
	"	[-l <file name of a list of summary_table.tsv files>]\n",
	"	-o <output file name root>\n",
	"\n",	
	"This script will read in each summary table, specified in the -i option and combine the counts\n",
	"into a single summary table.  \n",
	"\n",
	"The sample names and counts are preserved.  If a category is found in an input summary\n",
	"table but not in another, 0's will be added for padding.\n",
	"\n",
	"If you use the -i option, then the argument is for example: summary_table1.tsv,summary_table2.tsv,summary_tableN.tsv\n",
	"If you use the -l option, then the argument is a file name, where the file contains a summary table name per line.\n",
	"\n",
	"\n");

if((!length(opt$input_file_list)&&(!length(opt$input_file_listfile))) || !length(opt$output_file)){
	cat(usage);
	q(status=-1);
}

InputFileNameList=opt$input_file_list;
InputFileNameListFile=opt$input_file_list_file;
OutputFileName=opt$output_file;

output_fn=gsub("\\.summary_table\\.tsv$", "", OutputFileName);
output_fn=gsub("\\.summary_table\\.xls$", "", output_fn);

cat("Output Filename Root: ", output_fn, "\n", sep="");

###############################################################################

filenames_list=character();
if(length(InputFileNameListFile)){
	filenames_list=scan(InputFileNameListFile, what=character());
}else{
	filenames_list=strsplit(InputFileNameList,",")[[1]];
}

cat("Input Filenames:\n");
print(filenames_list);
cat("\n");

###############################################################################

load_summary_table=function(summary_table_fn){
        # Load data
        cat("Loading Matrix (", summary_table_fn, ") ...\n", sep="");
        inmat=as.matrix(read.table(summary_table_fn, sep="\t", header=TRUE, check.names=FALSE, row.names=1, comment.char="", quote=""))

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

merge_two_matrices=function(mat1, mat2){

	if(is.null(dim(mat1))){
		return(mat2);
	}
	
	nsamp1=nrow(mat1);
	ncat1=ncol(mat1);

	nsamp2=nrow(mat2);
	ncat2=ncol(mat2);

	catnames1=colnames(mat1);
	catnames2=colnames(mat2);

	sampnames1=rownames(mat1);
	sampnames2=rownames(mat2);

	# Figure out outmat col/row names/dimensions
	all_cat=unique(c(catnames1, catnames2));
	num_new_rows=nsamp1+nsamp2;
	num_new_cols=length(all_cat);

	# Allocate 0 matrix
	cat("Combined Matrix Dim: ", num_new_rows, " x ", num_new_cols, "\n", sep="");
	outmat=matrix(0, nrow=num_new_rows, ncol=num_new_cols);
	colnames(outmat)=sort(all_cat);
	rownames(outmat)=sort(c(sampnames1, sampnames2));

	# Copy data into outmat
	for(samp in sampnames1){
		outmat[samp, catnames1]=mat1[samp, catnames1];
	}
	for(samp in sampnames2){
		outmat[samp, catnames2]=mat2[samp, catnames2];
	}
	
	return(outmat);
}

###############################################################################

out_mat=numeric();
total_counts=0;
for(fn in filenames_list){
	tmp_mat=load_summary_table(fn);
	total_counts=total_counts+sum(tmp_mat);
	out_mat=merge_two_matrices(out_mat, tmp_mat);
}

###############################################################################
# Validate counts

cat("\n");
cat("Total Counts across individual summary tables: ", total_counts, "\n", sep="");
output_counts=sum(out_mat);
cat("Total Counts in merged summary table: ", output_counts, "\n", sep="");

if(output_counts==total_counts){
	cat("Great.  Total counts are equal to sum of individual counts...\n");
}else{
	cat("ERROR!!!  Total vs sum of individual counts not equal!!!\n");
	quit(-1);
}

###############################################################################
# Output
cat("\nWriting New Matrix...\n");
fc=file(paste(output_fn, ".summary_table.tsv", sep=""), "w");

write(paste("sample_id", "total", paste(colnames(out_mat), collapse="\t"), sep="\t"), file=fc);
sample_names=rownames(out_mat);
num_samples=nrow(out_mat);
for(samp_idx in 1:num_samples){
        total=sum(out_mat[samp_idx,]);
        outline=paste(sample_names[samp_idx],total,paste(out_mat[samp_idx,], collapse="\t"), sep="\t");
        write(outline, file=fc);
}
close(fc);

###############################################################################

cat("\nDone.\n")
warns=warnings();
if(!is.null(warns)){
        print(warnings());
}
q(status=0)

