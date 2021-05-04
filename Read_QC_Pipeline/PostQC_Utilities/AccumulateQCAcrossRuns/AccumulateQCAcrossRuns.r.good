#!/usr/bin/env Rscript

###############################################################################

library('openxlsx');
library('getopt');

options(useFancyQuotes=F);
options(width=200);

params=c(
	"input_list", "i", 1, "character",
	"output", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input tab-separated>\n",
	"	-o <output filename root>\n",
	"\n",
	"This script will read in a tab separated list of QC output files:\n",
	"	<file path> \\t <run name> \\n\n",
	"\n",
	"The QC output files should the format:\n",
	"	<sample id> \\t <read depth> \\n\n",
	"\n",
	"The output will be tsv files with the specified root:\n",
	"  Each row will be a unique sample id across the runs with replicate extension removed.\n",
	"  Each column will contain the read depth for each run\n",
	"  The Final column will be the total.\n",
	"\n");

if(!length(opt$input_list) || !length(opt$output)){
	cat(usage);
	q(status=-1);
}

InputFileList=opt$input_list;
OutputFNameRoot=opt$output;

cat("Input Filename: ", InputFileList, "\n");
cat("Output Filename: ", OutputFNameRoot, "\n");
cat("\n");

##############################################################################

load_filelist=function(fname){
	file_table=read.table(fname,  header=F, check.names=FALSE, 
		as.is=T, comment.char="", quote="", sep="\t");

	file_table=apply(file_table, 1:2, as.character);

	dimen=dim(file_table);
	cat("Rows Loaded: ", dimen[1], "\n");
	cat("Cols Loaded: ", dimen[2], "\n");

	if(dimen[2]==1){
		num_rows=dimen[1];
		extracted_name=character(num_rows);
		for(i in 1:num_rows){
			# Remove leading path and file extension
			extracted_name[i]=tail(strsplit(file_table[i,1], "/")[[1]], 1);
			extracted_name[i]=gsub("\\.tsv$", "", extracted_name[i]);
		}
		file_table=cbind(file_table, extracted_name);
	}	
	
	colnames(file_table)=c("Path", "DisplayName");

	return(file_table);
}

#------------------------------------------------------------------------------

load_counts=function(fname){

	count_table=read.table(fname,  header=F, check.names=FALSE, 
		as.is=T, comment.char="", quote="", sep="\t");

	count_table=count_table[,1:2];
	colnames(count_table)=c("SampleID", "Counts");

	if(!is.numeric(count_table[1,2])){
		count_table=count_table[-1,];
	}

	dimen=dim(count_table);
	cat("Rows Loaded: ", dimen[1], "\n");
	cat("Cols Loaded: ", dimen[2], "\n");

	# Convert _'s to .'s
	for(i in 1:dimen[1]){
		count_table[i,1]=gsub("_",".", count_table[i,1]);
	}

	# Remove tail replicate ID
	CollapsedID=character(dimen[1]);
	for(i in 1:dimen[1]){
		CollapsedID[i]=gsub("\\.\\d+$", "", count_table[i,1]);
	}
	count_table=cbind(count_table, CollapsedID);

	return(count_table);
}

#------------------------------------------------------------------------------

sample_ids_to_metadata=function(x, min_cat=1, max_cat=50){
	num_ids=length(x);
	sp=strsplit(x, "\\.");

	split_lengths=c(num_ids);
	for(i in 1:num_ids){
		split_lengths[i]=length(sp[[i]]);
	}

	uniq_split_lengths=unique(split_lengths);
	max_split_lengths=max(uniq_split_lengths);
	min_split_lengths=min(uniq_split_lengths);

	cat("Unique split lengths:\n");
	print(uniq_split_lengths);
	cat("Max split lengths:\n");
	print(max_split_lengths);
	cat("Min split lengths:\n");
	print(min_split_lengths);

	meta_matrix=matrix("NA", nrow=num_ids, ncol=max_split_lengths);
	for(smpix in 1:num_ids){
		cur_len=length(sp[[smpix]]);
		meta_matrix[smpix,]=c(sp[[smpix]], rep("NA", max_split_lengths-cur_len));
	}

	is_unique=logical(max_split_lengths);
	unique_list=list();
	for(cix in 1:max_split_lengths){
		unique_list[[cix]]=unique(meta_matrix[, cix]);
		num_categories=length(unique_list[[cix]]);
		is_unique[cix]=num_categories>=min_cat && num_categories<=max_cat;
	}

	results=list();
	results[["matrix"]]=meta_matrix;
	results[["keep"]]=is_unique;
	results[["unique_list"]]=unique_list;
	results[["num_columns"]]=max_split_lengths;
	
	#print(results);
	return(results);

}

##############################################################################

# Load list of target files
target_list=load_filelist(InputFileList);

cat("Target List:\n");
print(target_list);
cat("\n");

# Load counts into associateive array/list
num_targets=nrow(target_list);
target_values=list();
for(tix in 1:num_targets){
	cat("Loading: ", target_list[tix,1], " (", target_list[tix,2], ")\n");
	target_values[[tix]]=load_counts(target_list[tix,1]);
}
cat("\n");

# Accumulate sample IDs to find unique list
unique_sample_ids=character();
for(tix in 1:num_targets){
	unique_sample_ids=c(unique_sample_ids, as.character(target_values[[tix]][,"CollapsedID"]));
}
unique_sample_ids=sort(unique(unique_sample_ids));
num_unique_sample_ids=length(unique_sample_ids);

cat("Unique Sample IDs:\n");
print(unique_sample_ids);
cat("\n");

# Accumulation Matrices
breakout_matrix=matrix("0", nrow=num_unique_sample_ids, ncol=num_targets);
rownames(breakout_matrix)=unique_sample_ids;
colnames(breakout_matrix)=target_list[,2];

accumout_matrix=matrix(0, nrow=num_unique_sample_ids, ncol=num_targets);
rownames(accumout_matrix)=unique_sample_ids;
colnames(accumout_matrix)=target_list[,2];

#print(breakout_matrix);
#print(accumout_matrix);

attempts=matrix(0, nrow=num_unique_sample_ids, ncol=3);
rownames(attempts)=unique_sample_ids;
colnames(attempts)=c("Attempts", "TotalReads", "AvgReadsPerAttempt");

cat("Accumulating runs into matrix...\n");
for(tix in 1:num_targets){
	count_tab=target_values[[tix]];
	num_samples=nrow(count_tab);
	for(six in 1:num_samples){

		samp_id=as.character(count_tab[six,"CollapsedID"]);
		numeric_counts=as.numeric(count_tab[six,"Counts"]);

		attempts[samp_id, "Attempts"]=attempts[samp_id, "Attempts"]+1;
		attempts[samp_id, "TotalReads"]=attempts[samp_id, "TotalReads"]+numeric_counts;


		if(breakout_matrix[samp_id, tix]=="0"){
			breakout_matrix[samp_id, tix]=count_tab[six,"Counts"];
			accumout_matrix[samp_id, tix]=numeric_counts;
		}else{
			breakout_matrix[samp_id, tix]=
				paste(breakout_matrix[samp_id, tix], "+", count_tab[six,"Counts"]);
			accumout_matrix[samp_id, tix]=
				accumout_matrix[samp_id, tix] + numeric_counts;
		}
	}
}

attempts[,"AvgReadsPerAttempt"]=round(attempts[,"TotalReads"]/attempts[,"Attempts"], 0);
totals=apply(accumout_matrix, 1, sum);
attempts=cbind(rownames(attempts), attempts);
colnames(attempts)=c("SampleID", "Attempts", "TotalReads", "AvgReadsPerAttempt");

###############################################################################

cat("Broken Out:\n");
print(breakout_matrix);
cat("\n");
cat("Accumulated:\n");
print(accumout_matrix);
cat("\n");
cat("Attempts:\n");
print(attempts);

totals=apply(accumout_matrix, 1, sum);

output_text_matrix=matrix("", nrow=num_unique_sample_ids, ncol=num_targets+3);
colnames(output_text_matrix)=c("SampleID", c(target_list[,2]), "=", "Total");
rownames(output_text_matrix)=unique_sample_ids;
output_text_matrix[,"SampleID"]=unique_sample_ids;
output_text_matrix[,target_list[,2]]=breakout_matrix[,target_list[,2]];
output_text_matrix[,"Total"]=paste(totals);
print(output_text_matrix, quote=F);

detailed_text_matrix=output_text_matrix;
write.table(detailed_text_matrix, file=paste(OutputFNameRoot, ".accumulated.detail.tsv", sep=""),
	 quote=F, sep="\t", row.names=F, col.names=T);

output_text_matrix[,target_list[,2]]=accumout_matrix[,target_list[,2]];

cat("\n");
cat("Output as Text:\n");

output_text_matrix[,target_list[,2]]=accumout_matrix[,target_list[,2]];
cat("\n");
print(output_text_matrix, quote=F);

write.table(output_text_matrix, file=paste(OutputFNameRoot, ".accumulated.simplify.tsv", sep=""),
	 quote=F, sep="\t", row.names=F, col.names=T);

fname=paste(OutputFNameRoot, ".accumulated.attempts.tsv", sep="");
write.table(attempts, file=fname, quote=F, sep="\t", row.names=F, col.names=T);

##############################################################################
# Generate 2-way contigency tables

meta_table=sample_ids_to_metadata(unique_sample_ids);

make_contingency_table=function(data, c1, c2){

	arr1=data[["matrix"]][,c1];	
	arr2=data[["matrix"]][,c2];

	uniq1=sort(unique(arr1));
	uniq2=sort(unique(arr2));

	num_cat1=length(uniq1);
	num_cat2=length(uniq2);

	num_samp=nrow(data[["matrix"]]);

	cont_mat=matrix(0, nrow=num_cat1, ncol=num_cat2);
	rownames(cont_mat)=uniq1;
	colnames(cont_mat)=uniq2;

	for(i in 1:num_samp){
		cont_mat[arr1[i], arr2[i]]=cont_mat[arr1[i], arr2[i]]+1;
	}
		
	return(cont_mat);
	
}


contigency_table_list=list();
contigency_table_name=character();
tix=1;

cat("Metadata Columns: ", meta_table[["num_columns"]], "\n");

for(i in 1:meta_table[["num_columns"]]){

	if(!meta_table[["keep"]][i]) next;

	for(j in 1:meta_table[["num_columns"]]){

		if(!meta_table[["keep"]][j]) next;

		if(j>i){
			contigency_table_list[[tix]]=make_contingency_table(meta_table, i, j);
			contigency_table_name[tix]=paste("Table ", i, " x ", j, sep="");
			
			tix=tix+1;
		}
	}
}

num_cont_tabs=tix-1;

##############################################################################
# Prepare Sheets and output Excel spreadsheet

print(contigency_table_list);

workbook=createWorkbook();

#print(detailed_text_matrix);
addWorksheet(wb=workbook, sheetName="Detailed");
writeData(workbook, sheet=1, detailed_text_matrix);

addWorksheet(wb=workbook, sheetName="Simplified");
writeData(workbook, sheet=2, output_text_matrix);

addWorksheet(wb=workbook, sheetName="Attempts");
writeData(workbook, sheet=3, attempts);

if(num_cont_tabs>0){
	for(i in 1:num_cont_tabs){
		cat("Adding Sheet: ", contigency_table_name[i], "\n");
		addWorksheet(wb=workbook, sheetName=contigency_table_name[i]);
		writeData(workbook, sheet=3+i, attempts);
	}
}

saveWorkbook(workbook, paste(OutputFNameRoot, ".xlsx", sep=""), overwrite=T);

##############################################################################

cat("\nDone.\n");

print(warnings());
q(status=0);
