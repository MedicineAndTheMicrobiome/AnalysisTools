#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);
options(width=120);

params=c(
	"input_list_cli", "i", 2, "character",
	"input_list_file", "l", 2, "character",
	"output", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	Use one of the input options:\n",
	"	[-i <input, comma separated list of tsv files>]\n",
	"	[-l <input, text file with list of target files>]\n",
	"\n",
	"	-o <output tab_separated column data file>\n",
	"\n",
	"  Eg. -i /path/filename.tsv,/anotherpath/anotherfile.tsv\n",
	"\n",
	"This script will read in the metadata files specified in the comma-separated\n",
	"list of the -i option, and merge the files based on the shared column names.\n",
	"\n",
	"Some diagnostic information will be generated so you can see which\n",
	"columns could be merged.  \n",
	"\n",
	"The total number of rows output will be:\n",
	"	(total num rows input)-(number of files)+1\n",
	"\n");

if(!length(opt$output)){
	cat(usage);
	q(status=-1);
}

InputFNameList=opt$input_list_cli;
InputFNameFile=opt$input_list_file;

Formulas=opt$formulas;
OutputFName=opt$output;

use_cli=NA;
if(length(InputFNameList)){
	cat("List specified on commmand line.\n");
	use_cli=T;
}
if(length(InputFNameFile)){
	cat("List specified in file.\n");
	use_cli=F;
}

if(is.na(use_cli)){
	cat("Please specify one of the input options.\n");
}

cat("Output Filename: ", OutputFName, "\n");
cat("\n");


target_list=character();
if(use_cli){
	target_list=strsplit(InputFNameList, ",")[[1]];	
}else{
	target_list=read.table(InputFNameFile, as.is=T, sep="\n", header=F, check.names=F,
			comment.char="", quote="")[,1];
}

cat("Target list of files to merge:\n");
print(target_list);


##############################################################################

load_factors=function(fname){
	factors=data.frame(read.table(fname,  header=TRUE, check.names=FALSE, as.is=T, 
		comment.char="", quote="", sep="\t"));

	#print(factors);

	dimen=dim(factors);
	cat("Rows Loaded: ", dimen[1], "\n");
	cat("Cols Loaded: ", dimen[2], "\n");

	return(factors);
}

write_factors=function(fname, table){

	dimen=dim(table);
	cat("Rows Exporting: ", dimen[1], "\n");
	cat("Cols Exporting: ", dimen[2], "\n");
	
	write.table(table, fname, quote=F, row.names=F, sep="\t");

}

##############################################################################

num_target_files=length(target_list);

cat("Files to Merge: \n");
names(target_list)=LETTERS[1:length(target_list)];
print(target_list);

loaded_factors=list();
loaded_columns=list();
loaded_numrows=list();
unique_columns=character();
total_rows=0;
for(i in 1:num_target_files){
	
	cur_fname=target_list[i];
	cat("\nLoading: ", cur_fname, "\n", sep="");
	loaded_factors[[i]]=load_factors(cur_fname);	
	loaded_columns[[i]]=colnames(loaded_factors[[i]]);
	loaded_numrows[[i]]=nrow(loaded_factors[[i]]);
	unique_columns=unique(union(unique_columns, loaded_columns[[i]]));
	total_rows=total_rows+nrow(loaded_numrows[[i]]);

}

cat("\n");
cat("**********************************************************************\n");
cat("Unique column names across files specified:\n\n");
print(unique_columns);

lowercase_uniq_cols=tolower(unique_columns);
formatting_uniq_cols=gsub("[\\.\\_]", "", lowercase_uniq_cols);
#print(formatting_uniq_cols);
#print(lowercase_uniq_cols);

cat("\n");
cat("**********************************************************************\n");

merge_warnings=F;

if(length(unique(lowercase_uniq_cols))!=length(unique(unique_columns))){
	cat("WARNING: Case insensitive conflicts.\n");
	cat("  (This means multiple columns have the same name if we ignore case.)\n");
	counts=table(lowercase_uniq_cols);
	dups=counts>1;
	print(counts[dups]);
	merge_warnings=T;
}else{
	cat("Good. No case insensitive conflicts.\n");
}

if(length(unique(formatting_uniq_cols))!=length(unique(unique_columns))){
	cat("WARNING: Formatting conflicts.\n");
	cat("  (This means multiple columns have the same name if we ignore _'s and .'s)\n");
	counts=table(formatting_uniq_cols);
	dups=counts>1;
	print(counts[dups]);
	merge_warnings=T;
}else{
	cat("Good. No formatting conflicts.\n");
}

cat("**********************************************************************\n");
cat("\n");

unique_columns=unique_columns[order(lowercase_uniq_cols)];
#print(unique_columns);
num_unique_cols=length(unique_columns);

shared_matrix=matrix(0, nrow=num_unique_cols, ncol=num_target_files);
rownames(shared_matrix)=unique_columns;
colnames(shared_matrix)=LETTERS[1:num_target_files];

for(i in 1:num_target_files){
	file_mat=loaded_factors[[i]];
	num_nas=apply(file_mat, 2, function(x){100*sum(!is.na(x))/length(x);});
	
	for(vname in names(num_nas)){
		shared_matrix[vname, i]=num_nas[vname];
	}

}

cat("Non-NA Percentages, sorted alphabetically:\n");
cat("  (Check for variables with similar names with mutually exclusive counts.)\n");
cat("\n");
print(shared_matrix, digits=1);
cat("\n");

cat("Non-NA Percentages, sorted by total non-NAs:\n");
cat("  (Look for variables with missing data due to inconsistent naming.)\n");
cat("\n");
shared_totals=apply(shared_matrix, 1, sum);
ordered_by_total_ix=order(shared_totals, decreasing=T, method="shell");
shared_matrix=shared_matrix[ordered_by_total_ix,,drop=F];
print(shared_matrix, digits=1);
cat("\n");

##############################################################################

output_matrix=matrix(NA, nrow=0, ncol=length(unique_columns));
colnames(output_matrix)=rownames(shared_matrix);

for(i in 1:num_target_files){
	tmp_output_matrix=matrix(NA, nrow=loaded_numrows[[i]], ncol=length(unique_columns));
	colnames(tmp_output_matrix)=rownames(shared_matrix);

	varnames=colnames(loaded_factors[[i]]);

	for(j in 1:length(varnames)){	
		tmp_output_matrix[, varnames[j]]=loaded_factors[[i]][,varnames[j]];
	}

	#print(tmp_output_matrix);

	cat("Accumulating: ", nrow(tmp_output_matrix), " Rows.\n");
	output_matrix=rbind(output_matrix, tmp_output_matrix);
}

#print(output_matrix);

##############################################################################

write_factors(OutputFName, output_matrix);

##############################################################################

if(merge_warnings){
	cat("\n");
	cat("***********************************************\n");
	cat("Check program output.  Merge warnings detected.\n");
	cat("***********************************************\n");
}

cat("\nDone.\n");

print(warnings());
q(status=0);
