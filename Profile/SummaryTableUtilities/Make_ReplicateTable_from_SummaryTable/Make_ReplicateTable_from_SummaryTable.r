#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_summary_table", "s", 1, "character",
	"output_filename_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-s <input .summary_table.tsv>\n",
	"	-o <output filename root>\n",
	"\n",	
	"This script will generate a mapping table based on the\n",
	"input summary table's sample IDs.  The first column\n",
	"is the original sample ID with replicate extension, the\n",
	"the second column will contain the name of the collapsed or\n",
	"underlying sample ID.\n",
	"\n",
	"For example:\n",
	"0208.004012.20180910.BALR	0208.004012.20180910.BALR\n",
	"0208.004012.20180910.BALR.2	0208.004012.20180910.BALR\n",
	"\n",
	"Mapping File format:\n",
	"	<current sample id>\\t<new sample id>\\n\n",
	"\n");

if(!length(opt$input_summary_table) || !length(opt$output_filename_root)){
	cat(usage);
	q(status=-1);
}

InputSummaryTable=opt$input_summary_table;
OuputFilenameRoot=opt$output_filename_root;

MappingFilename=paste(OuputFilenameRoot, ".mapping.tsv", sep="");
UniqueIDsFilename=paste(OuputFilenameRoot, ".unique_ids.lst", sep="");

###############################################################################

cat("\n")
cat("Input Summary Table Name: ", InputSummaryTable, "\n");
cat("Output Mapping File Name Root: ", OuputFilenameRoot, "\n");       
cat("\n");
cat("Mapping File Name: ", MappingFilename, "\n");
cat("Unique ID list: ", UniqueIDsFilename, "\n");
cat("\n");

###############################################################################
###############################################################################

# Load data
cat("Loading Matrix...\n");
inmat=as.matrix(read.table(InputSummaryTable, sep="\t", header=TRUE, 
	check.names=FALSE, comment.char="*", row.names=1))

#cat("Input Excerpt:\n")
#print(inmat[1:10,1:10]);

###############################################################################

sample_ids=rownames(inmat);

cat("Sample ID Excerpt: \n");
print(sample_ids[1:10]);

num_samp_ids=length(sample_ids);
cat("Number of Sample IDs: ", num_samp_ids, "\n");

collapsed_id=character(num_samp_ids);


cat("\n");
for(i in 1:num_samp_ids){

	cur_samp_id=sample_ids[i];

	toks=strsplit(cur_samp_id, "\\.")[[1]];
	num_toks=length(toks);
	if(!is.na(as.numeric(toks[num_toks]))){
		toks=toks[-num_toks];
	}

	collapsed_id[i]=paste(toks, collapse=".");
}

###############################################################################
# Output

cat("Writing Mapping File...\n");
fc=file(MappingFilename, "w");

cat(file=fc, "Original\tCollapsed\n");
for(samp_idx in 1:num_samp_ids){
	cat(file=fc, sample_ids[samp_idx], "\t", collapsed_id[samp_idx], "\n", sep="");
}
close(fc);	

###############################################################################

cat("Writing unique list of sample IDs...\n");
uniq_ids=sort(unique(collapsed_id));
fc=file(UniqueIDsFilename, "w");
cat(file=fc, "SampleID\n");
for(samp_idx in 1:length(uniq_ids)){
	cat(file=fc, uniq_ids[samp_idx], "\n", sep="");
}
close(fc);

###############################################################################

writeLines("Done.\n")
warns=warnings();
if(!is.null(warns)){
	print(warns);
}

q(status=0)
