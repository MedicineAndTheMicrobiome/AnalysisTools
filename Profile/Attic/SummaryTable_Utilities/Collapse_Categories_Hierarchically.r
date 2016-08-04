#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"level", "l", 1, "character",
	"output_file_root", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	-l <level to collapse to>\n",
	"	[-o <output summary table file name root, default, based on input file name>]\n",
	"\n",	
	"This script will read in a summary table, and then collapse the rows after removing\n",
	"levels after the one you have specified.\n",
	"\n",
	"Level 1: Kingdom\n",
	"Level 2: Phylum\n",
	"Level 3: Class\n",
	"Level 4: Order\n",
	"Level 5: Family\n",
	"Level 6: Genus\n",
	"\n",
	"For example if the category is: \n",
	"	Bacteria_Actinobacteria_Actinobacteria_Actinomycetales_Kineosporiaceae_Kineococcus\n",
	"\n",
	"Choosing a level of 1, would collapse everything down to Bacteria.\n",
	"Choosing a level of 2, would collapse all the Actinobacteria together.\n",
	"\n",
	"*_incertae_sedis will be kept together.\n",
	"\n",
	"\n");

if(!length(opt$input_file) || !length(opt$level)){
	cat(usage);
	q(status=-1);
}

###############################################################################
# Get paramaters

InputFileName=opt$input_file;
Level=opt$level;

if(!length(opt$output_file)){
	OutputFileName=gsub("\\.summary_table\\.xls$", "", opt$input_file);
}else{
	OutputFileName=opt$output_file_root;
}


###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("Level: ", Level, "\n");

###############################################################################
###############################################################################

# Load data
inmat=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="*", row.names=1))
inmat=inmat[,2:ncol(inmat)];
#cat("Original Matrix:\n")
#print(inmat);

# Summarze what we've loaded
num_samples=nrow(inmat);
num_categories=ncol(inmat);
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");

###############################################################################

categories=colnames(inmat);

# Modify incertae sedis's so they stick together when we split by underscore
categories=gsub("genera_incertae_sedis", "-incertae-sedis", categories);
categories=gsub("_incertae_sedis", "-incertae-sedis", categories);
#print(categories);

# Split by underscore
split_list=strsplit(categories, "_");

# Eliminate levels after the one specified and rejoin them back together
categories_clean=character(num_categories);
for(i in 1:num_categories){
	num_levels=length(split_list[[i]]);
	max_level=min(Level, num_levels);
	
	categories_clean[i]=paste(split_list[[i]][1:max_level], collapse="_");
}

# Create list of unique categories
unique_categories=sort(unique(categories_clean));
num_uniq_cat=length(unique_categories);

# Accumulate the counts together
counts=matrix(0, nrow=num_samples, ncol=num_uniq_cat);
for(i in 1:num_uniq_cat){
	matching=(unique_categories[i]==categories_clean);
	#cat("Collapsing ", unique_categories[i], "\n");
	extracted=inmat[,matching];
	
	if(sum(matching)>1){
		counts[,i]=apply(extracted, 1, sum);
	}else{
		counts[,i]=extracted;
	}
}

# Label the matrix
outmat=counts;
rownames(outmat)=rownames(inmat);
colnames(outmat)=unique_categories;

# Come up with new filename
if(Level==1){
	lstring="Kingdom";
}else if(Level==2){
	lstring="Phylum";
}else if(Level==3){
	lstring="Class";
}else if(Level==4){
	lstring="Order";
}else if(Level==5){
	lstring="Family";
}else if(Level==6){
	lstring="Genus";
}else{
	cat("Unknown level.\n");
	lstring="Unknown";
}
fc=file(paste(OutputFileName, ".", lstring, ".summary_table.xls", sep=""), "w");

# Write out summary file table 
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

cat("Done.\n")
if(!is.null(warnings)){
	print(warnings());
}

q(status=0)
