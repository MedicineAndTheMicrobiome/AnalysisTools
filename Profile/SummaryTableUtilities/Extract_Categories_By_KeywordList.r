#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"keyword_file", "k", 1, "character",
	"case_sensitive", "c", 2, "logical",
	"stand_alone_word", "w", 2, "logical",
	"exact_match", "E", 2, "logical",
	"output_file", "o", 2, "character",
	"exclude_remaining", "R", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary_table.tsv file name>\n",
	"	-k <file name for list of keywords look for>\n",
	"	[-c <case sensitive, default = not>]\n",
	"	[-w <stand alone word, default = not>]\n",
	"	[-E (exact match flag, default = not set)]\n",
	"	[-o <output file name root>]\n",
	"	[-R (Do not include 'Remaining' Category)]\n",
	"\n",	
	"This script will search through the summary table's categories\n",
	"and extract the columns that match the keywords in your keyword\n",
	"file. \n", 
	"\n",
	"If you use the -c option, then the search will not be case sensitive.\n",
	"If you use the -w option, then the word you look for needs to exist by\n",
	"  itself.  It can not be part of another word.",
	"\n",
	"The default settings will get you the most hits.\n",
	"\n",
	"\n");

if(!length(opt$input_file) || !length(opt$keyword_file)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$output_file)){
	outputroot=gsub("\\.summary_table\\.xls", "", opt$input_file);
	outputroot=gsub("\\.summary_table\\.tsv", "", outputroot);
	OutputFileName = paste(outputroot, ".kextrct", sep="");
}else{
	OutputFileName=opt$output_file;
	OutputFileName=gsub("\\.summary_table\\.xls$", "", OutputFileName);
	OutputFileName=gsub("\\.summary_table\\.tsv$", "", OutputFileName);
}


###############################################################################

InputFileName=opt$input_file;
RemoveList=opt$remove_list;
KeywordList=opt$keyword_file;

StandAlone=F;
if(length(opt$stand_alone_word)){
	StandAlone=opt$stand_alone_word;
}

CaseSensitive=F;
if(length(opt$case_sensitive)){
	CaseSensitive=opt$case_sensitive;
}

IncludeRemaining=T;
if(length(opt$exclude_remaining)){
	IncludeRemaining=F
}

ExactMatch=F;
if(length(opt$exact_match)){
	ExactMatch=T;
	cat("Exact Match flag set, ignoring word and case settings.\n");
	StandAlone=F;
	CaseSensitive=T;
}

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Keyword List: ", KeywordList, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("\n");
cat("Case Sensitive: ", CaseSensitive, "\n");
cat("Stand Alone Words: ", StandAlone, "\n");
cat("Include 'Remaining' Category: ", IncludeRemaining, "\n");
cat("\n");

###############################################################################
###############################################################################

load_summary_table=function(summary_table_fn){
        # Load data
        cat("Loading Matrix (", summary_table_fn, ") ...\n", sep="");
        inmat=as.matrix(read.table(summary_table_fn, sep="\t", header=TRUE, check.names=FALSE, row.names=1, quote=""))

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

load_ids=function(list_fn){
        cat("Loading List (", list_fn, ") ...\n", sep="");
        list=scan(file=list_fn, what="complex", sep="\t");
        return(list);
}

###############################################################################
# Load counts matrix
counts_mat=load_summary_table(InputFileName);
num_categories=ncol(counts_mat);
num_samples=nrow(counts_mat);

# Load remove list IDs
keyword_list=load_ids(KeywordList);
num_keywords=length(keyword_list);
cat("Num keywords: ", num_keywords, "\n");
cat("Keywords: \n");
print(keyword_list);
cat("\n");

# Get the category names
category_names=colnames(counts_mat);

if(!CaseSensitive){
	keyword_list=tolower(keyword_list);
	target_category_names=tolower(category_names);
}else{
	target_category_names=category_names;
}

keep_idx=numeric();
for(kwidx in 1:num_keywords){
	cur_keyword=keyword_list[kwidx];
	
	cat("Looking for: ", cur_keyword, "\n");
	
	cur_keyword=gsub("^\\s+", "", cur_keyword);
	cur_keyword=gsub("\\s+$", "", cur_keyword);

	if(!ExactMatch){
		if(StandAlone){
			begin_kw=paste("(^", cur_keyword, "\\W+)", sep="");
			end_kw=paste("(\\W+", cur_keyword, "$)", sep="");
			middle_kw=paste("(\\W+", cur_keyword, "\\W+)", sep="");
			cur_keyword=paste(begin_kw, end_kw, middle_kw, sep="|");
		}

		cat("Pattern:", cur_keyword, "\n");
		
		hits=grep(cur_keyword, target_category_names, perl=T);
		num_hits=length(hits);
		cat("Num hits: ", num_hits, "\n");
	}else{
		find_ix=(cur_keyword==category_names);
		if(any(find_ix)){
			cat("Found.\n");
			hits=which(find_ix);
		}else{
			cat("Not Found.\n");
			hits=c();
		}

	}

	print(category_names[hits]);
	keep_idx=c(keep_idx, hits);

	cat("\n");
}

keep_idx=sort(unique(keep_idx));

# Extract columns
outmat=counts_mat[,keep_idx, drop=F];

if(IncludeRemaining){
	# Compute Remaining
	totals=apply(counts_mat, 1, sum);
	extr_totals=apply(outmat, 1, sum);
	remaining_totals=totals-extr_totals;

	# Append remaining to outmat
	outmat=cbind(outmat, remaining_totals);
	colnames(outmat)=c(category_names[keep_idx], "Remaining");
}

###############################################################################

# Output
cat("\nWriting New Matrix...\n");
fc=file(paste(OutputFileName, ".summary_table.tsv", sep=""), "w");

write(paste("sample_id", "total", paste(colnames(outmat), collapse="\t"), sep="\t"), file=fc);
sample_names=rownames(counts_mat);
for(samp_idx in 1:num_samples){
	total=sum(outmat[samp_idx,]);
	outline=paste(sample_names[samp_idx],total,paste(outmat[samp_idx,], collapse="\t"), sep="\t");
	write(outline, file=fc);
}
close(fc);	

###############################################################################

cat("\nWriting extracted categories...\n");
fc=file(paste(OutputFileName, ".identifed.txt", sep=""), "w");
num_to_keep=length(keep_idx);
for(i in keep_idx){
	cat(file=fc, category_names[i], "\n", sep="");
}
close(fc);

###############################################################################

cat("\nDone.\n")
warns=warnings();
if(!is.null(warns)){
	print(warnings());
}
q(status=0)
