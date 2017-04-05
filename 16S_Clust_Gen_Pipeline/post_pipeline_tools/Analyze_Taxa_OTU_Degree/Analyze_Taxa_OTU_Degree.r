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
	otu_names=rownames(otu_sizes);

	# Based on taxonomic level, determine what part of name can be truncated
	if(num_col==level_ix){
		shortened_names=otu_names;
	}else{
		shortened_names=apply(associated_otus[, (level_ix+1):num_col, drop=F], 1, function(x){paste(x, collapse=";")} );
		shortened_names=paste(shortened_names, " (", otu_names, ")", sep=""); 
	}

	# Abbreviate/renamed undetermined/unclassified levels
	undetermined_ix=(shortened_names=="");
	if(any(undetermined_ix)){
		shortened_names[undetermined_ix]="(Undetermined)";
	}
	shortened_names=gsub("_unclassified", "(U)", shortened_names);
	
	# If there are too many categories to plot, reduce to MAX_CAT
	num_categories=nrow(associated_otus);
	num_otus=num_categories;
	if(num_categories>MAX_CAT){
		num_categories=MAX_CAT;
		associated_otus=associated_otus[1:MAX_CAT, , drop=F];
		otu_sizes=otu_sizes[1:MAX_CAT, , drop=F]
		shortened_names=shortened_names[1:MAX_CAT];
		is_only_top=paste(" (Top ", num_categories,")", sep="");
	}else{
		is_only_top="";
	}

	# Identify max OTU size, so we can resize the OTU size
	max_y_axis=max(10, max(otu_sizes));
	total_otu_seq=sum(otu_sizes);
	otu_prop=otu_sizes/total_otu_seq;

	# Generate barplot
	mids=barplot(height=otu_sizes[,"Size"], names.arg=shortened_names, main=paste(taxa_name, is_only_top, sep=" "), 
		xaxt="n", yaxt="n",
		ylab="Counts", ylim=c(0, max_y_axis*1.1));

	# y axis (count and proportion)
	y_ticks=round(seq(0,max_y_axis, length.out=5));
	axis(side=2, at=y_ticks, labels=y_ticks, las=2);
	prop_labels=round(y_ticks/total_otu_seq, 2);
	prop_lteq1=prop_labels<=1;
	axis(side=4, at=y_ticks[prop_lteq1], labels=prop_labels[prop_lteq1], las=2, cex.axis=.75);
	mtext("Proportion", side=4, line=3, cex=.75);
	
	# Compute and put down labels
	bar_width=mids[2]-mids[1];
	plot_range=par()$usr;
	plot_height=plot_range[4];
	label_size=min(c(1,.7*bar_width/par()$cxy[1]));
	text(mids-par()$cxy[1]/2, rep(-par()$cxy[2]/2, num_col), shortened_names, srt=-45, xpd=T, pos=4, cex=label_size);

	div=diversity(otu_prop[,"Size"]);	

	# Print other stats under title
	offset=1;
	title(paste("Total OTUs: ", num_otus, sep=""), line=offset, cex.main=.75, font.main=1);
	title(paste("Total Sequences: ", total_otu_seq, sep=""), line=offset-.75, cex.main=.75, font.main=1);
	title(paste("Shannon: ", round(div[1],4), sep=""), line=offset-1.5, cex.main=.75, font.main=1);
	title(paste("Evenness: ", round(div[2],4), sep=""), line=offset-2.25, cex.main=.75, font.main=1);

}

diversity=function(prob){
	# Compute shannon diversity index and evenness
	prob=as.matrix(prob);
	num_cat=length(prob);
	equi_prob=rep(1/num_cat, num_cat);
	shan=-sum(prob*log(prob));
	shan_max=-sum(equi_prob*log(equi_prob));
	return(c(shan, shan/shan_max));
}

plot_text=function(strings, max_lines=50){

        plot_page=function(strings){
                orig_par=par(no.readonly=T);

                par(mfrow=c(1,1));
                par(family="Courier");
                par(oma=rep(.5,4));
                par(mar=rep(0,4));

                num_lines=length(strings);

                top=max(as.integer(num_lines), 40);

                plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
                        xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                        );
                for(i in 1:num_lines){
                        #cat(strings[i], "\n", sep="");
                        text(0, top-i, strings[i], pos=4, cex=.7);
                }

                par(orig_par);
        }

        num_lines=length(strings);
        num_pages=ceiling(num_lines / max_lines);
        #cat("Num Pages: ", num_pages, "\n");
        for(page_ix in 1:num_pages){
                start=(page_ix-1)*max_lines+1;
                end=start+max_lines-1;
                end=min(end, num_lines);
                ##print(c(start,end));
                plot_page(strings[start:end]);
        }
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

plots_per_level=c(1,1,3,3,3,4);
space_for_names=c(15, 25, 20, 15, 10, 5);
options(width=240);

for(level_ix in 1:6){

	cat("Accumulating OTUs at Level: ", LEVELS[level_ix], "\n");

	#-----------------------------------------------------------------------------
	# Get current taxa name by joining components to current level
	cur_taxa=apply(unique_taxa_mat[,1:level_ix, drop=F], 1, function(x){paste(x, collapse=";")});
	null_taxa=grep(";$", cur_taxa);
	if(length(null_taxa)){
		cur_taxa=cur_taxa[-null_taxa];
	}
	unique_cur_taxa=unique(cur_taxa);

	#-----------------------------------------------------------------------------
	# Summary Statistics
	par(mfrow=c(1,1));
	par(oma=c(1,0,1,0));
	lines=character();
	stat_mat=matrix(0, nrow=0, ncol=6, dimnames=list(c(), c("NumSeq", "NumOTUs", "Seq/OTU", "GreatestProp", "Shannon", "Evenness")));
	taxa_names=character();
	for(taxa in unique_cur_taxa){
		splits=strsplit(taxa, ";")[[1]];

		otus=pull_otus(splits, taxa_mat);
		otu_names=rownames(otus);
		otu_sizes=as.vector(inmat[otu_names, "Size", drop=F]);
		
		# Compute stats on this set of OTUs
		total_seq=sum(otu_sizes);
		otu_prop=otu_sizes/total_seq;
		total_otus=nrow(otu_sizes);
		max_abund=max(otu_prop);
		seq_per_otu=total_seq/total_otus;
		div=diversity(otu_prop);

		# Add to line buffer
		stat_mat=rbind(stat_mat, 
			c(total_seq, total_otus, round(seq_per_otu,2), round(max_abund,2), round(div[1],3), round(div[2],3))
		);

		# Only keep up to the 2 lowest taxonomic names
		taxa_names=c(taxa_names, paste(tail(splits,2), collapse=";"));
	}
	rownames(stat_mat)=taxa_names;
	plot_text(capture.output(stat_mat));
	mtext(paste("[",LEVELS[level_ix],"]", sep=""), outer=T, col="blue");

	#-----------------------------------------------------------------------------
	# Plot Barcharts
	par(mfrow=c(plots_per_level[level_ix], 1));
	par(mar=c(space_for_names[level_ix],4,5,space_for_names[level_ix]/sqrt(2)+1));
	par(oma=c(0,0,2,0));
	for(taxa in unique_cur_taxa){
		splits=strsplit(taxa, ";")[[1]];
		#print(splits);
		otus=pull_otus(splits, taxa_mat);
		otu_name=rownames(otus);
		otu_info=inmat[otu_name,"Size", drop=F];
		tot_seq=sum(otu_info[,"Size"]);
		if(nrow(otus)>1){
			plot_otus(taxa, level_ix, otus, otu_info);
		}
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
