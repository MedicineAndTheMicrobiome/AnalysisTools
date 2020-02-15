#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"minimum_total", "c", 1, "character",
	"output_file", "o", 2, "character",
	"generate_plot", "p", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	-c <cutoff>\n",
	"	[-o <output summary table file name>\n",
	"	[-p (generate plot)]\n",
	"\n",	
	"This script will read in the summary table, and recompute the total for each sample,\n",
	"then only output the samples with total reads greater than each of the specified cutoffs.\n",
	"\n",
	"The the cutoffs may be a single value or a comma separated list of values, such as:\n",
	"	-c 750,1000,2000,3000\n",
	"\n",
	"\n");

if(!length(opt$input_file) || !length(opt$minimum_total)){
	cat(usage);
	q(status=-1);
}

minimums_arr=sort(unique(as.numeric(strsplit(opt$minimum_total, ",")[[1]])));

InputFileName=opt$input_file;
MinimumTotalCutoff=opt$minimum_total;
GeneratePlot=opt$generate_plot;
OutputFileRoot=opt$output_file;

if(length(OutputFileRoot)==0){
	OutputFileRoot=InputFileName;
}

OutputFileRoot=gsub("\\.summary_table\\.tsv$", "", OutputFileRoot);
OutputFileRoot=gsub("\\.summary_table\\.xls$", "", OutputFileRoot);

OutputPDFFileName=paste(OutputFileRoot, ".minFilt.hist.pdf", sep="");

if(length(GeneratePlot)==0){
	GeneratePlot=FALSE;
}else{
	GeneratePlot=TRUE;
}


###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name Root: ", OutputFileRoot, "\n");       
cat("Minimum Total Cutoffs: ", "\n");
print(minimums_arr);
cat("Generate Plot: ", GeneratePlot, "\n");
cat("\n");

###############################################################################
###############################################################################

# Load data
cat("Loading Matrix...\n");
inmat<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="*", row.names=1))
#cat("Original Matrix:\n")
#print(inmat);

# Grab columns we need into a vector, ignore totals, we won't trust it.
counts_mat=inmat[,2:(ncol(inmat))];
#print(counts_mat);

num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");

###############################################################################

# Compute Totals
totals=apply(counts_mat, 1, sum);
max_tot=max(totals);
increments=min(diff(minimums_arr));

cat("Maximum Totals: ", max_tot, "\n", sep="");
cat("Increments: ", increments, "\n", sep="");

# Plot Histogram
if(GeneratePlot){
	pdf(OutputPDFFileName, height=8.5, width=11);

	bins=seq(0, max_tot+increments, increments);
	hist(totals, xlab="Sample Totals", main="Sample Total Distribution (Frequencies)", breaks=bins);
	abline(v=minimums_arr, col="blue");

	bins=sort(unique(c(0, minimums_arr, max_tot)));
	hist(totals, xlab="Sample Totals", main="Sample Total Distribution (Unequal Bin Sizes)", breaks=bins);
	abline(v=minimums_arr, col="blue");
}


###############################################################################
# Output

write_summary_table=function(outname, counts_matrix){
	cat("Writing New Matrix...\n");
	fc=file(outname, "w");

	write(paste("Sample_ID", "Total", paste(colnames(outmat), collapse="\t"), sep="\t"), file=fc);
	out_num_samples=nrow(outmat);
	sample_names=rownames(outmat);
	if(out_num_samples>0){
		for(samp_idx in 1:out_num_samples){
			total=sum(outmat[samp_idx,]);
			outline=paste(sample_names[samp_idx],total,paste(outmat[samp_idx,], collapse="\t"), sep="\t");
			write(outline, file=fc);
		}
	}
	close(fc);	
}

write_list=function(fname, char_arr){
	cat("Writing list...\n");
	fc=file(fname, "w");

	list_len=length(char_arr);
	if(list_len>=1){
		for(i in 1:list_len){
			write(char_arr[i], file=fc);
		}
	}else{
		write("<EOF>", file=fc);
	}
	close(fc);
}

cat("\n\n");

max_cutoff=max(minimums_arr);
for(cutoff in minimums_arr){

	cat("Working on cutoff: ", cutoff, "\n");
	keep_idx=totals>=cutoff;
	num_samples_to_keep=sum(keep_idx);
	cat("Number of Samples to Remove: ", num_samples-num_samples_to_keep, "\n", sep="");
	cat("Number of Samples to Keep  : ", num_samples_to_keep, "\n", sep="");

	#Subset out rows to keep
	samp_ids=rownames(counts_mat);
	outmat=counts_mat[keep_idx,];
	excluded=samp_ids[!keep_idx];

	zpad=paste("%0", floor(log10(max_cutoff))+1, "g", sep="");
	sumtab_name=paste(OutputFileRoot, ".min_", sprintf(zpad, cutoff), ".summary_table.tsv", sep=""); 
	write_summary_table(sumtab_name, outmat);

	list_out_name=paste(OutputFileRoot, ".lt_", sprintf(zpad, cutoff), ".exclusion.tsv", sep="");
	write_list(list_out_name, excluded);

	cat("\n");
}



###############################################################################

writeLines("Done.\n")
warns=warnings();
if(length(warns)>0){
	print(warnings());
}

q(status=0)
