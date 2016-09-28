#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
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

plot_connected_figure=function(coordinates, offsets_mat, groups_per_plot=3, title=""){
	sorted_sids=sort(rownames(offsets_mat));
	coordinates=coordinates[sorted_sids,];
	offsets_mat=offsets_mat[sorted_sids,];

	print(offsets_mat);
	print(coordinates);

	# Get Unique Groups
	groups=sort(unique(offsets_mat[,"Group ID"]));
	num_groups=length(groups);

	colors=rainbow(num_groups);
	color_mat_dim=ceiling(sqrt(num_groups));
	colors=as.vector(t(matrix(colors, nrow=color_mat_dim, ncol=color_mat_dim)));
	palette(colors);

	# Get limits of points
	extra_margin=.2;
	x_range=range(coordinates[,1]);
	y_range=range(coordinates[,2]);
	x_ext=abs(x_range[2]-x_range[1]);
	y_ext=abs(y_range[2]-y_range[1]);
	xlim=c(x_range[1]-x_ext*extra_margin, x_range[2]+x_ext*extra_margin);
	ylim=c(y_range[1]-y_ext*extra_margin, y_range[2]+y_ext*extra_margin);
	cat("\nPoint ranges:\n");
	cat("X:\n");
	print(x_range);
	cat("Y:\n");
	print(y_range);
	cat("Plot ranges:\n");
	cat("X:\n");
	print(xlim);
	cat("Y:\n");
	print(ylim);
	cat("\n");

	# Plot all samples
	plot(0, main=title, xlab="Dim 1", ylab="Dim 2", type="n", xlim=xlim, ylim=ylim);
	for(i in 1:num_groups){
		grp_subset=which(offsets_mat[,"Group ID"]==groups[i]);
		num_members=length(grp_subset);
		print(grp_subset);

		offsets_subset=offsets_mat[grp_subset,];
		coord_subset=coordinates[grp_subset,];

		sort_ix=order(offsets_subset[,"Offsets"], decreasing=F);

		offsets_subset=offsets_subset[sort_ix,];
		coord_subset=coord_subset[sort_ix,];

		print(offsets_subset);
		print(coord_subset);
			
		points(coord_subset, type="b", col=i, pch=20, cex=.5);
		points(coord_subset[c(1, 1, num_members),], type="p", col=i, pch=c(17, 1, 15), cex=c(1, 2, 1.25));
	}

	# Plot subset of samples
	for(i in 1:num_groups){
		if(((i-1) %% groups_per_plot)==0){
			plot(0, main=title, xlab="Dim 1", ylab="Dim 2", type="n", xlim=xlim, ylim=ylim);
		}

		cat("Plotting: ", groups[i], "\n");
		grp_subset=which(offsets_mat[,"Group ID"]==groups[i]);
		num_members=length(grp_subset);
		print(grp_subset);

		offsets_subset=offsets_mat[grp_subset,];
		coord_subset=coordinates[grp_subset,];

		sort_ix=order(offsets_subset[,"Offsets"], decreasing=F);

		offsets_subset=offsets_subset[sort_ix,];
		coord_subset=coord_subset[sort_ix,];

		print(offsets_subset);
		print(coord_subset);
			
		# Label start, stop, and group id
		points(coord_subset, type="b", col=i, pch=20, cex=.5);
		points(coord_subset[c(1, 1, num_members),], type="p", col=i, pch=c(17, 1, 15), cex=c(1, 2, 1.25));
		text(coord_subset[1,1], coord_subset[1,2], labels=groups[i], col="black", pos=1, cex=.75, font=2);

		# Label offsets
		offset_ix=2:num_members;
		text(coord_subset[offset_ix,1], coord_subset[offset_ix,2], labels=offsets_subset[offset_ix,"Offsets"], col=i, adj=c(.5,-.75), cex=.5, font=3);
	}
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

dist_mat=vegdist(normalized_mat, method="horn");
#dist_mat=dist(normalized_mat);
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
warn=warnings();
if(length(warn)){
	print(warn);
}
q(status=0)
