#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);
options(width=120);

params=c(
	"input_list", "i", 1, "character",
	"output", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input tab-separated>\n",
	"	-o <output tab_separated>\n",
	"\n",
	"This script will read in a tab separated list of QC output files:\n",
	"	<file path> \\t <run name> \\n\n",
	"\n",
	"The QC output files should the format:\n",
	"	<sample id> \\t <read depth> \\n\n",
	"\n",
	"The output will be a matrix:\n",
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

##############################################################################

# Load list of target files
target_list=load_filelist(InputFileList);

cat("Target List:\n");
print(target_list);
cat("\n");

num_targets=nrow(target_list);
target_values=list();
for(tix in 1:num_targets){
	cat("Loading: ", target_list[tix,1], " (", target_list[tix,2], ")\n");
	target_values[[tix]]=load_counts(target_list[tix,1]);
}
cat("\n");

unique_sample_ids=character();
for(tix in 1:num_targets){
	unique_sample_ids=as.character(target_values[[tix]][,"CollapsedID"]);
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

for(tix in 1:num_targets){
	count_tab=target_values[[tix]];
	num_samples=nrow(count_tab);
	for(six in 1:num_samples){
		samp_id=count_tab[six,"CollapsedID"];
		if(breakout_matrix[samp_id, tix]=="0"){
			breakout_matrix[samp_id, tix]=count_tab[six,"Counts"];
			accumout_matrix[samp_id, tix]=as.numeric(count_tab[six,"Counts"]);
		}else{
			breakout_matrix[samp_id, tix]=
				paste(breakout_matrix[samp_id, tix], "+", count_tab[six,"Counts"]);
			accumout_matrix[samp_id, tix]=
				accumout_matrix[samp_id, tix] + as.numeric(count_tab[six,"Counts"]);
		}
	}
}

cat("Broken Out:\n");
print(breakout_matrix);
cat("\n");
cat("Accumulated:\n");
print(accumout_matrix);
cat("\n");

totals=apply(accumout_matrix, 1, sum);
#print(totals);

output_text_matrix=matrix("", nrow=num_unique_sample_ids, ncol=num_targets+3);
colnames(output_text_matrix)=c("SampleID", c(target_list[,2]), "=", "Total");
rownames(output_text_matrix)=unique_sample_ids;
output_text_matrix[,"SampleID"]=unique_sample_ids;
output_text_matrix[,target_list[,2]]=breakout_matrix[,target_list[,2]];
output_text_matrix[,"Total"]=paste(totals);
print(output_text_matrix, quote=F);

write.table(output_text_matrix, file=paste(OutputFNameRoot, ".accumulated.detail.tsv", sep=""),
	 quote=F, sep="\t", row.names=F, col.names=T);

output_text_matrix[,target_list[,2]]=accumout_matrix[,target_list[,2]];
cat("\n");
print(output_text_matrix, quote=F);

write.table(output_text_matrix, file=paste(OutputFNameRoot, ".accumulated.simplify.tsv", sep=""),
	 quote=F, sep="\t", row.names=F, col.names=T);

##############################################################################

cat("\nDone.\n");

print(warnings());
q(status=0);
