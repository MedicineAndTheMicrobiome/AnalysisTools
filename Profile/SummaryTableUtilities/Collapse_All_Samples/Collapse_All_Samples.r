#!/usr/bin/env Rscript

###############################################################################

library('getopt');

source('~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r');
source('~/git/AnalysisTools/Metadata/OutputFileLibrary/OutputFileLibrary.r');


params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 1, "character",
	"sample_id_spec", "s", 2, "character",
	"sample_id_samproot", "r", 2, "character",
	"sample_id_fileroot", "f", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];



usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table>\n",
	"	-o <output summary table>\n",
	"\n",
	"	New Sample ID options:\n",
	"	[-s <new sample ID to put into summary table>]\n",
	"	[-r (Use shared prefix across existing samples)]\n",
	"	[-f (use root of input file name)]\n",
	"\n",
	"This script will read in the summary table and collapse\n",
	"all the samples in the file into a new single sample (row)\n",
	"\n",
	"-s will use whatever you specify as the new sample ID\n",
	"\n",
	"-r will try to find a sample id by looking at your current sample ids\n",
	"  and find a prefix shared by all samples.\n",
	"  For example, if your samples are: samp.1, samp.2, samp.3, then the\n",
	"  the sample ID used will be samp.\n",
	"\n",
	"-f will use the root name of the input summary table.\n",
	"  For example, if your input file is ~/data/my_data.summary_table.tsv\n",
	"  then the sample ID used will be my_data\n",
	"\n");

if(!length(opt$input_file) || !length(opt$output_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileName=opt$output_file;

SampleID_Spec="";
SampleID_SampRoot="";
SampleID_FileRoot="";

if(length(opt$sample_id_spec)){
	SampleID_Spec=opt$sample_id_spec;
}

if(length(opt$sample_id_samproot)){
	SampleID_SampRoot=opt$sample_id_samproot;
}

if(length(opt$sample_id_fileroot)){
	SampleID_FileRoot=opt$sample_id_fileroot;
}

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("\n");
cat("SampleID as Specified: ", SampleID_Spec, "\n");
cat("SampleID as Common SampleID Root: ", SampleID_SampRoot, "\n");
cat("SampleID as File Name Root: ", SampleID_FileRoot, "\n");

###############################################################################
###############################################################################

find_common_root=function(ids_arr){
	max_len=max(nchar(ids_arr));

	for(i in 1:max_len){
		current_char=substr(ids_arr, i, i);
		if(!all(current_char[1]==current_char)){
			break;
		}
	}
	
	common_root=substr(ids_arr[1], 1, i-1);
	common_root=gsub("_*$", "", common_root);
	common_root=gsub("\\.*$", "", common_root);
	common_root=gsub(" *$", "", common_root);

	return(common_root);
}

#------------------------------------------------------------------------------

find_filename_root=function(filepath){
	comp=strsplit(filepath, "/")[[1]];
	num_comp=length(comp);
	filename=comp[num_comp];
	fileroot=gsub(".summary_table.tsv", "", filename);
	return(fileroot);
}

#------------------------------------------------------------------------------

# Load data
counts_mat=load_summary_file(InputFileName);

num_samples=nrow(counts_mat);
num_categories=ncol(counts_mat);

cat("\n");
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");
cat("\n");

samp_ids=rownames(counts_mat);

sample_id=find_common_root(samp_ids);

if(SampleID_Spec!= ""){
	sample_id=SampleID_Spec;
}else if(SampleID_FileRoot != ""){
	sample_id=find_filename_root(InputFileName);
}

summed=apply(counts_mat, 2, sum);
outmat=matrix(summed, ncol=num_categories, nrow=1)

colnames(outmat)=colnames(counts_mat);
rownames(outmat)=sample_id;

cat("Sample ID used: ", sample_id, "\n");
write_summary_file(outmat, OutputFileName);

###############################################################################

cat("\nDone.\n");
if(!is.null(warnings())){
	print(warnings());
}

q(status=0);
