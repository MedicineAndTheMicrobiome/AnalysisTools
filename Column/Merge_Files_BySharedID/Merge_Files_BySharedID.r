#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);
options(width=120);

params=c(
	"input_list_file", "l", 1, "character",
	"output", "o", 1, "character",
	"out_primary_key", "C", 1, "character",
	"shared_id", "c", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-l <input, text file with list of target files and target primary key>\n",
	"	-o <output tab_separated column data file>\n",
	"	-C <output shared id / primary key>\n",
	"	[-c <default shared id column name, default=output shared id>]\n",
	"\n",
	"This script will read in the metadata files specified in the comma-separated\n",
	"list of the -i option or file with the -l option, and merge the files based\n",
	"on the shared column names.\n",
	"\n",
	"Only the shared IDs across all the files will be returned.\n",
	"\n",
	"You can use the following command to make the target list file:\n",
	"The list would be stored in the file: targets\n",
	"\n",
	"cat <<EOF > targets\n",
	"<file>\\t<primary_key>\\n\n",
	"<file>\\t<primary_key>\\n\n",
	"<file>\\t<primary_key>\\n\n",
	"...\n",
	"<file>\\t<primary_key>\\n\n",
	"EOF\n",
	"\n\n",
	"If the <primary_key> is not specified, the output shared id column name will be used.\n",
	"\n");

if(!length(opt$output) || !length(opt$out_primary_key) || !length(opt$input_list_file)){
	cat(usage);
	q(status=-1);
}

InputFileListFName=opt$input_list_file;
OutputFName=opt$output;
OutputPrimaryKey=opt$out_primary_key;

if(length(opt$shared_id)){
	DefaultSharedID=opt$shared_id;
}else{
	DefaultSharedID=OutputPrimaryKey;
}

cat("\n");
cat("Input File List: ", InputFileListFName, "\n");
cat("Default Shared ID: ", DefaultSharedID, "\n");
cat("Output Filename: ", OutputFName, "\n");
cat("Output Primary Key: ", OutputPrimaryKey, "\n");
cat("\n");


target_list_rec=read.table(InputFileListFName, as.is=T, sep="\n", header=F, check.names=F,
			comment.char="", quote="")[,1];

cat("Target list contents:\n");
print(target_list_rec);
cat("\n");

target_rec_split=strsplit(target_list_rec, "\t");
#print(target_rec_split);

num_targets=length(target_rec_split);
target_rec_mat=matrix(NA, nrow=num_targets, ncol=2);
colnames(target_rec_mat)=c("FileName", "PrimaryKey");
for(i in 1:num_targets){
	target_rec_mat[i,"FileName"]=target_rec_split[[i]][1];	

	if(!is.na(target_rec_split[[i]][2])){
		target_rec_mat[i,"PrimaryKey"]=target_rec_split[[i]][2];	
	}else{
		target_rec_mat[i,"PrimaryKey"]=DefaultSharedID;
	}
}

print(target_rec_mat);
cat("\n");

##############################################################################

load_factors=function(fname, prim_key){
	factors=data.frame(read.table(fname,  header=TRUE, 
		row.names=prim_key,
		check.names=FALSE, as.is=T, 
		comment.char="#", quote="", 
		sep="\t"));

	#print(factors);

	dimen=dim(factors);
	cat("Primary Key: ", prim_key, "\n", sep="");
	cat("Rows Loaded: ", dimen[1], "\n");
	cat("Cols Loaded: ", dimen[2], "\n");

	return(factors);
}

write_factors=function(fname, table, prim_key){

	dimen=dim(table);
	cat("Primary Key: ", prim_key, "\n", sep="");
	cat("Rows Exporting: ", dimen[1], "\n");
	cat("Cols Exporting: ", dimen[2], "\n");
	
	out_mat=cbind(rownames(table), table);
	colnames(out_mat)=c(prim_key, colnames(table));
	write.table(out_mat, fname, quote=F, row.names=F, sep="\t");

}

##############################################################################

cat("Number of Target Files: ", num_targets, "\n\n", sep="");

files_arr=target_rec_mat[,"FileName"];
factors_list=list();
shared_ids_list=list();

for(i in 1:num_targets){
	fname=files_arr[i];
	shared_id=target_rec_mat[i,"PrimaryKey"];
	cat("Loading: ", fname, "\n", sep="");
	factors_list[[fname]]=load_factors(fname, shared_id);
	shared_ids_list[[fname]]=rownames(factors_list[[fname]]);
}

cat("\n\n");

shared_ids=shared_ids_list[[1]];
for(fname in files_arr){
	cat("Intersecting with IDs from: ", fname, "\n");
	shared_ids=intersect(shared_ids, shared_ids_list[[fname]]);
	cat("Num Shared: ", length(shared_ids), "\n\n");	
}

##############################################################################

shared_ids=sort(shared_ids);
combined_factors=matrix(NA, nrow=length(shared_ids), ncol=0);
rownames(combined_factors)=shared_ids;

for(fname in files_arr){
	combined_factors=cbind(combined_factors, factors_list[[fname]][shared_ids,]);
}

write_factors(OutputFName, combined_factors, OutputPrimaryKey);

##############################################################################

print(warnings());
q(status=0);
