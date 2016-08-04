#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file_list", "i", 1, "character",
	"output_file", "o", 1, "character",
	"sum_flag", "s", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input comma separated list of summary_table.xls,summary_table.xls,summary_table.xls,...>\n",
	"	-o <output file name root>\n",
	"	[-s (sum flag, sum samples together, without normalizing)]\n",
	"\n",	
	"This script will read in each summary table, normalize them, then by sample ID merge the samples together.\n",
	"If the -s flag is used, the individual samples will not be normalized, essentially a weighted average of samples.\n",
	"\n",
	"\n");

if(!length(opt$input_file) || !length(opt$output_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileName=opt$output_file;
SumSamples=TRUE;

if(length(opt$sum_flag)==0){
	SumSamples=FALSE;
}

OutputFileNameSummaryTable=paste(OutputFileName, ".summary_table.xls", sep="");

###############################################################################

filenames_list=strsplit(InputFileName,",")[[1]];

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("\n");

###############################################################################
###############################################################################

# Load data
cat("Loading Matrix...\n");
in_matrices=list();
num_inputs=length(filenames_list);

# To store samples/categories across all matrices
all_samples=character();
all_categories=character();

for(i in 1:num_inputs){

	# Load matrix
	cat("Loading: ", filenames_list[i], "\n", sep="");
	inmat=as.matrix(read.table(filenames_list[i], sep="\t", header=TRUE, check.names=FALSE, comment.char="*", row.names=1,quote=""))

	# Grab columns we need into a vector, ignore totals, we won't trust it.
	counts_mat=inmat[,2:(ncol(inmat)), drop=F];
	if(is.null(nrow(counts_mat))){
		counts_mat=matrix(counts_mat, nrow=1, ncol=length(counts_mat));
		rownames(counts_mat)=rownames(inmat)[1];
	}else{
		rownames(counts_mat)=rownames(inmat);
	}
	colnames(counts_mat)=colnames(inmat)[2:(ncol(inmat))];
	

	# Get the dimensions of the counts
	num_samples=nrow(counts_mat);
	num_categories=ncol(counts_mat);
	cat("Num Samples: ", num_samples, "\n");
	cat("Num Categories: ", num_categories, "\n");
	cat("\n");

	# Keep track of all the samples/categories across all matrices we want to merge
	all_samples=unique(c(all_samples, rownames(counts_mat)));
	all_categories=unique(c(all_categories, colnames(counts_mat)));

	# Compute Totals
	total=apply(counts_mat, 1, sum);
	#cat("Total:\n");
	#print(total);

	# Normalize (if it's requested)	
	if(!SumSamples){
		for(s in 1:num_samples){
			counts_mat[s,]=counts_mat[s,]/total[s];	
		}
	}

	# Save normalized counts
	in_matrices[[i]]=counts_mat;

}

###############################################################################

# Allocate matrix to store all counts
final_sample_count=length(all_samples);
final_category_count=length(all_categories);
accumulated_counts=matrix(0.0, ncol=final_category_count, nrow=final_sample_count);
all_categories_sorted=sort(all_categories);
all_samples_sorted=sort(all_samples);
colnames(accumulated_counts)=all_categories_sorted;
rownames(accumulated_counts)=all_samples_sorted;

# Cycle through all the matrices and sum them up
for(i in 1:num_inputs){

	current_matrix=in_matrices[[i]];
	#print(current_matrix);

	cur_categories=colnames(current_matrix);
	cur_samples=rownames(current_matrix);

	num_col=ncol(current_matrix);
	num_row=nrow(current_matrix);

	for(samp in 1:num_row){
		row=which(cur_samples[samp]==all_samples_sorted);

		for(cat in 1:num_col){
			col=which(cur_categories[cat]==all_categories_sorted);

			accumulated_counts[row,col]=
				accumulated_counts[row,col]+current_matrix[samp, cat];
		}
	}
}

#print(accumulated_counts);

if(!SumSamples){
	# Normalize counts so they add up to one again
	outmat=matrix(0, nrow=final_sample_count, ncol=final_category_count);
	sum=apply(accumulated_counts, 1, sum);
	for(i in 1:final_sample_count){
		outmat[i,]=accumulated_counts[i,]/sum[i];
	}
	rownames(outmat)=rownames(accumulated_counts);
	colnames(outmat)=colnames(accumulated_counts);
}else{
	# Just output the nominal counts
	outmat=accumulated_counts;
}

#print(outmat);

###############################################################################
# Output
cat("Writing New Matrix...\n");
fc=file(OutputFileName, "w");

write(paste("sample_id", "total", paste(colnames(outmat), collapse="\t"), sep="\t"), file=fc);
out_num_samples=nrow(outmat);
sample_names=rownames(outmat);
for(samp_idx in 1:out_num_samples){
	total=sum(outmat[samp_idx,]);
	outline=paste(sample_names[samp_idx],total,paste(outmat[samp_idx,], collapse="\t"), sep="\t");
	write(outline, file=fc);
}
close(fc);	

###############################################################################

writeLines("Done.\n")
warns=warnings();
if(length(warns)>0){
	print(warnings());
}

q(status=0)
