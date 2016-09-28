#!/usr/bin/env Rscript

###############################################################################

library(MASS)
library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"offset_file", "t", 1, "character",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	-t <offset file>\n",
	"	[-o <output file root name>]\n",
	"\n",
	"	This script will read in the summary table\n",
	"	and a file describing the time from the first\n",
	"	sample.\n",
	"\n",
	"	The format of the offset file is:\n",
	"\n",
	"	<sample id> \\t <sample grouping id> \\t <time stamp> \\n",
	"\n",
	"\n");

if(!length(opt$input_file) || !length(opt$offset_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OffsetFileName=opt$offset_file;

if(length(opt$output_file)>0){
	OutputFileRoot=opt$output_file;
}else{
	OutputFileRoot=InputFileName;
	OutputFileRoot=gsub("\\.summary_table\\.tsv$", "", OutputFileRoot);
	OutputFileRoot=gsub("\\.summary_table\\.xls$", "", OutputFileRoot);
	cat("No output file root specified.  Using input file name as root.\n");
}

###############################################################################

OutputPDF = paste(OutputFileRoot, ".mds_ts.pdf", sep="");
cat("Output PDF file name: ", OutputPDF, "\n", sep="");
pdf(OutputPDF,width=8.5,height=8.5)

###############################################################################

load_offset=function(fname){
        cat("Loading Offsets: ", fname, "\n");
        offsets_mat=read.delim(fname,  header=FALSE, row.names=1, sep="\t", comment.char="#", quote="");
	colnames(offsets_mat)=c("Group ID", "Offsets");

	# reset offsets
	groups=unique(offsets_mat[,"Group ID"]);
	
	cat("Groups:\n");
	print(groups);
	cat("\n");

	# Reset offsets so they are relative to the first/smallest sample
	for(gid in groups){
		offsets=offsets_mat[gid==offsets_mat[,"Group ID"], "Offsets"];
		min_off=min(offsets);
		offsets_mat[gid==offsets_mat[,"Group ID"], "Offsets"]=offsets-min_off;
	}

	return(offsets_mat);
}

load_summary_file=function(fname){
        cat("Loading Summary Table: ", fname, "\n");
        inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", row.names=1))
        counts_mat=inmat[,2:(ncol(inmat))];
        return(counts_mat);
}

normalize=function(counts){
        totals=apply(counts, 1, sum);
        num_samples=nrow(counts);
        normalized=matrix(0, nrow=nrow(counts), ncol=ncol(counts));

        for(i in 1:num_samples){
                normalized[i,]=counts[i,]/totals[i];
        }

        colnames(normalized)=colnames(counts);
        rownames(normalized)=rownames(counts);
        return(normalized);
}

plot_connected_figure=function(coordinates, offsets_mat, title=""){
	sorted_sids=sort(rownames(offsets_mat));

	print(offsets_mat[sorted_sids,]);
	print(coordinates[sorted_sids,]);

	groups=sort(unique(offsets_mat[,"Group ID"]));
	colors=nrow(offsets_mat);
	num_groups=length(groups);
	for(i in 1:num_groups){
		colors[groups[i]==offsets_mat[sorted_sids,"Group ID"]]=i;
	}
	print(colors);

	palette(rainbow(num_groups));

	plot(coordinates[sorted_sids,], main=title, xlab="Dim 1", ylab="Dim 2", type="n");
	text(coordinates[sorted_sids,], labels=offsets_mat[sorted_sids,"Group ID"], col=colors);
}

###############################################################################
###############################################################################

offset_mat=load_offset(OffsetFileName);
#print(offset_mat);

###############################################################################

counts_mat=load_summary_file(InputFileName);

#print(counts_mat);

normalized_mat=normalize(counts_mat);

#print(normalized_mat);

dist_mat=dist(normalized_mat);
#print(dist_mat);

# Remove 0 distances with very small number
for(i in 1:length(dist_mat)){
	if(dist_mat[i]==0){
		dist_mat[i]=1e-323;
	}
}

mds_coord=cmdscale(dist_mat, k=2);
isomds=isoMDS(dist_mat);
mds2_coord=isomds$points;

###############################################################################


plot_connected_figure(mds_coord, offset_mat, title="Metric MDS");
plot_connected_figure(mds2_coord, offset_mat, title="IsoMetric MDS");

##############################################################################

cat("Done.\n")
dev.off();

q(status=0)
