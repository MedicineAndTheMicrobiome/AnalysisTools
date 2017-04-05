#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 2, "character"
);

MAX_CAT=60;

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input: Mothur output *.cons.taxonomy file>\n",
	"	[-o <output filename root>]\n",
	"\n",	
	"This script will read in the input file and compute the degree of\n",
	"diversity within a particular taxonomic classification.\n",
	"\n",
	"In particular, the usage may be get an estimate of the amount\n",
	"of species level diversity and distribution within a genus.\n",
	"\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

###############################################################################

InputFileName=opt$input_file;
OutputFileName=opt$output_file;

if(!length(OutputFileName)){
	OutputFileName=gsub("\\.cons\\.taxonomy", "", InputFileName);
}

OutputFileName=paste(OutputFileName, ".otu_degree", sep="");

cat("\n");
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name Root: ", OutputFileName, "\n");       
cat("\n");

pdf(paste(OutputFileName, ".pdf", sep=""), height=11, width=8.5);

###############################################################################

pull_otus=function(taxa_arr, name_matrix){

	num_levels=length(taxa_arr);
	num_taxa=nrow(name_matrix);

	return_ix=rep(TRUE, num_taxa);

	for(i in 1:num_levels){
		return_ix=return_ix & taxa_arr[i]==name_matrix[,i];	
	}

	return(name_matrix[return_ix,, drop=F]);
}

plot_otus=function(taxa_name, level_ix, associated_otus, otu_sizes){
	cat("Plotting: ", taxa_name, "\n");
	num_col=ncol(associated_otus);

	if(num_col==level_ix){
		shortened_names=rownames(otu_sizes);
	}else{
		shortened_names=apply(associated_otus[, (level_ix+1):num_col, drop=F], 1, function(x){paste(x, collapse=";")} );
	}

	undetermined_ix=(shortened_names=="");
	if(any(undetermined_ix)){
		shortened_names[undetermined_ix]="(Undetermined)";
	}

	shortened_names=gsub("_unclassified", "(U)", shortened_names);
	
	num_categories=nrow(associated_otus);
	if(num_categories>MAX_CAT){
		num_categories=MAX_CAT;
		associated_otus=associated_otus[1:MAX_CAT, , drop=F];
		otu_sizes=otu_sizes[1:MAX_CAT, , drop=F]
		shortened_names=shortened_names[1:MAX_CAT];
		is_only_top=paste(" (Top ", num_categories,")", sep="");
	}else{
		is_only_top="";
	}

	

	max_otu_size=max(otu_sizes);
	print(shortened_names);
	mids=barplot(height=otu_sizes[,"Size"], names.arg=shortened_names, main=paste(taxa_name, is_only_top, sep=" "), xaxt="n",
		ylab="Counts", ylim=c(0, max(10, max_otu_size)));
	
	bar_width=mids[2]-mids[1];
	plot_range=par()$usr;
	plot_height=plot_range[4];
	label_size=min(c(1,.7*bar_width/par()$cxy[1]));
	cat("Label Size: ", label_size, "\n");
	text(mids-par()$cxy[1]/2, rep(-par()$cxy[2]/2, num_col), shortened_names, srt=-45, xpd=T, pos=4, cex=label_size);

}

###############################################################################

# Load data
inmat=read.delim(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="");

cat("Input Column Names:\n");
print(colnames(inmat));
rownames(inmat)=inmat[,"OTU"];
otu_names=inmat[,"OTU"];
num_otus=nrow(inmat);

# Remove (###) confidence calls from name
inmat[,"Taxonomy"]=gsub("\\(\\d+\\)", "", inmat[,"Taxonomy"]);

# Break down taxa
taxa_mat=matrix("", nrow=num_otus, ncol=6, dimnames=list(inmat[,"OTU"], 
	c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")));

for(row_ix in 1:num_otus){
	splits=strsplit(inmat[row_ix,"Taxonomy"], ";")[[1]];
	num_splits=length(splits);
	taxa_mat[row_ix, 1:num_splits]=splits;
}

#print(taxa_mat);

#otu_ix=pull_otus(c("Bacteria", "Firmicutes", "Clostridia"), taxa_mat);

# Remove taxa that are "unclassified"
taxa_mat_no_unclassified=apply(taxa_mat, c(1,2), function(x){return(gsub("unknown", "", gsub(".*_unclassified$", "", x)))});
# print(taxa_mat_no_unclassified);

# Concat and remove duplicates
unique_concat=sort(unique(apply(taxa_mat_no_unclassified, 1, function(x){return(paste(x, collapse=";"))})));
if(any(unique_concat==";;;;;")){
	unique_concat=unique_concat[-which(unique_concat==";;;;;")];
}
num_unique_taxa=length(unique_concat);

LEVELS=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus");
# Break down taxa again
unique_taxa_mat=matrix("", nrow=num_unique_taxa, ncol=6, dimnames=list(1:num_unique_taxa, 
	LEVELS));
for(row_ix in 1:num_unique_taxa){
	splits=strsplit(unique_concat[row_ix], ";")[[1]];
	num_splits=length(splits);
	unique_taxa_mat[row_ix, 1:num_splits]=splits;
}
#print(unique_taxa_mat);

plots_per_level=c(1,1,3,4,5,6);
space_for_names=c(20, 20, 20, 15, 10, 5);


for(level_ix in 1:6){

	par(mfrow=c(plots_per_level[level_ix], 1));
	par(mar=c(space_for_names[level_ix],4,4,space_for_names[level_ix]/sqrt(2)));
	par(oma=c(0,0,2,0));

	cat("Accumulating OTUs at Level: ", LEVELS[level_ix], "\n");

	cur_taxa=apply(unique_taxa_mat[,1:level_ix, drop=F], 1, function(x){paste(x, collapse=";")});
	null_taxa=grep(";$", cur_taxa);
	if(length(null_taxa)){
		cur_taxa=cur_taxa[-null_taxa];
	}
	unique_cur_taxa=unique(cur_taxa);

	for(taxa in unique_cur_taxa){
		cat("\n  Working on: ", taxa, "\n");
		splits=strsplit(taxa, ";")[[1]];
		#print(splits);
		otus=pull_otus(splits, taxa_mat);
		otu_name=rownames(otus);
		plot_otus(taxa, level_ix, otus, inmat[otu_name,"Size", drop=F]);
		mtext(paste("[",LEVELS[level_ix],"]", sep=""), outer=T, col="blue");
	}

	cat("\n");

}


###############################################################################

dev.off();

cat("Done.\n")
if(!is.null(warnings())){
	print(warnings());
}

q(status=0)
