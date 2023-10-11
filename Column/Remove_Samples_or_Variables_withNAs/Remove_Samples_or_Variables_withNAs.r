#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);
options(width=120);

params=c(
	"input", "i", 1, "character",
	"output", "o", 1, "character",
	"req_samples_fn", "s", 2, "character",
	"req_variables_fn", "v", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input tab-separated column data file>\n",
	"\n",
	"   One of: \n",
	"	[-s <required sample IDs list, filename>]\n",
	"	[-v <required variable list, filename]\n",
	"\n",
	"	-o <output filename>\n",
	"\n",
	"This script will either determine which samples are available\n",
	"based on required variables (-v), or determine which variables\n",
	"are available based on required samples (-s).\n",
	"\n",
	"Only one of -s or -v should be chosen.\n",
	"\n");

if(
	!length(opt$input) || 
	!length(opt$output)

){
	cat(usage);
	q(status=-1);
}

InputMetadata=opt$input;
OutputMetadata=opt$output;

RequiredVariablesFn=NULL;
RequiredSamplesFn=NULL;


if(length(opt$req_variables_fn)){
	RequiredVariablesFn=opt$req_variables_fn;
}
if(length(opt$req_samples_fn)){
	RequiredSamplesFn=opt$req_samples_fn;
}


cat("\n");
cat("Input Filename:", InputMetadata, "\n");
cat("Output Filename:", OutputMetadata, "\n");
cat("\n");
cat("Required Variables Filaname:", RequiredVariablesFn, "\n");
cat("Required Samples Filename:", RequiredSamplesFn, "\n");
cat("\n");

##############################################################################

load_factors=function(fname){

	cat("Loading: ", fname, "\n", sep="");

	factors=data.frame(read.delim(fname,  header=TRUE, check.names=FALSE,
		row.names=NULL,
		as.is=T, comment.char="", quote="", sep="\t"));

	rownames(factors)=factors[,1];
	sample_id_name=colnames(factors)[1]
	factors=factors[,-1];

	dimen=dim(factors);
	cat("Rows Loaded: ", dimen[1], "\n");
	cat("Cols Loaded: ", dimen[2], "\n");

	res=list();	
	res[["sample_id"]]=sample_id_name;
	res[["data"]]=factors;

	return(res);
}

write_factors=function(fname, table, sample_id_name){
	
	cat("Writing Metadata: ", fname, "\n", sep="");	

	dimen=dim(table);

	cnames=c(sample_id_name, colnames(table));
	table=cbind(rownames(table), table);
	colnames(table)=cnames;
	
	cat("Rows Exporting: ", dimen[1], "\n");
	cat("Cols Exporting: ", dimen[2], "\n");

	write.table(table, fname, quote=F, row.names=F, sep="\t");
}

load_list=function(fname){

	loaded=data.frame(read.delim(fname, header=F, check.names=FALSE,
		row.names=NULL, as.is=T, comment.char="#", quote="", sep="\t"));

	arr=loaded[,1];

	return(arr);
}

write_list=function(fname, arr){

	cat("Writing List: ", fname, "\n", sep="");	
	write.table(arr, fname, quote=F, row.names=F, col.names=F);

}


##############################################################################

# Load factors
res=load_factors(InputMetadata);
sample_id_cname=res[["sample_id"]];
in_data=res[["data"]];

#print(in_data);

num_samples=nrow(in_data);
num_variables=ncol(in_data);

cat("\n");
cat("Num Samples: ", num_samples, "\n");
cat("Num Variables: ", num_variables, "\n");

sample_ids=rownames(in_data);
variable_names=colnames(in_data);


cat("-----------------------------------------------------------------------------------------------------\n");
cat("Sample IDs:\n");
print(sample_ids);
cat("\n");
cat("Variables:\n");
print(variable_names);
cat("\n");
cat("-----------------------------------------------------------------------------------------------------\n");
cat("\n");
cat("=====================================================================================================\n");
# Load requirements
required_variable_list=NULL;
if(!is.null(RequiredVariablesFn)){
	required_variable_list=load_list(RequiredVariablesFn);
}
if(!is.null(required_variable_list)){
	cat("Required Variables:\n");
	print(required_variable_list);
}

required_sample_list=NULL;
if(!is.null(RequiredSamplesFn)){
	required_sample_list=load_list(RequiredSamplesFn);
}
if(!is.null(required_sample_list)){
	cat("Required Samples:\n");
	print(required_sample_list);
}

cat("=====================================================================================================\n");

##############################################################################

filtered=in_data;

num_nas=sum(is.na(filtered));
num_cells=prod(dim(filtered));
original_num_cells=num_cells;
original_num_nas=num_nas;

cat("\n");
cat("Original Stats:\n");
cat("Num Samples: ", dim(filtered)[1], "\n");
cat("Num Variables: ", dim(filtered)[2], "\n");
cat("Num Cells: ", num_cells, "\n");
cat("Num NAs: ", num_nas, "\n");
cat("Proportion NAs: ", round(num_nas/num_cells, 2), "\n");
cat("\n");


if(!is.null(required_variable_list)){

	original_sample_ids=rownames(filtered);
	filtered=filtered[,required_variable_list, drop=F];

	samples_woNAs=apply(filtered, 1, function(x){ !any(is.na(x));});
	filtered=filtered[samples_woNAs,,drop=F];

	removed_sample_ids=setdiff(original_sample_ids, original_sample_ids[samples_woNAs]);

	out_fname=paste(gsub("tsv", "", OutputMetadata), "removed_sample_ids.lst", sep="");
	write_list(out_fname, removed_sample_ids);

}

if(!is.null(required_sample_list)){

	original_variables=colnames(filtered);
	filtered=filtered[required_sample_list,, drop=F];

	variables_woNAs=apply(filtered, 2, function(x){ !any(is.na(x));});
	filtered=filtered[,variables_woNAs,drop=F];

	removed_variables=setdiff(original_variables, original_variables[variables_woNAs]);

	out_fname=paste(gsub("tsv", "", OutputMetadata), "removed_variables.lst", sep="");
	write_list(out_fname, removed_variables);
}

#print(filtered);

num_nas=sum(is.na(filtered));
num_cells=prod(dim(filtered));

cells_with_info_lost=original_num_cells-num_cells-original_num_nas;

cat("\n");
cat("Filtered Stats:\n");
cat("Num Samples: ", dim(filtered)[1], "\n");
cat("Num Variables: ", dim(filtered)[2], "\n");
cat("Num Cells: ", num_cells, "\n");
cat("Num NAs: ", num_nas, "\n");
cat("Proportion NAs: ", round(num_nas/num_cells, 2), "\n");
cat("Num Cells with information lost: ", cells_with_info_lost, "\n");
cat("\n");

write_factors(OutputMetadata, filtered, sample_id_cname);

##############################################################################

cat("\nDone.\n");

print(warnings());
q(status=0);
