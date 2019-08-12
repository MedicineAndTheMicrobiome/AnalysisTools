#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"offset_file", "t", 1, "character",
	"output_file", "o", 2, "character",
	"distance_type", "d", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

DEF_DIST="euclidean";

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	-t <offset file>\n",
	"	-d <distance type, def=", DEF_DIST, ">\n",
	"	[-o <output file root name>]\n",
	"\n",
	"	This script will read in the summary table\n",
	"	and a file describing the time from the first\n",
	"	sample.\n",
	"\n",
	"	The format of the offset file is:\n",
	"\n",
	"	<sample id> \\t <sample grouping/individual id> \\t <time stamp> \\t <cohort (treat/group) id>\\n",
	"\n",
	"	Distance Types:\n",	
	"		euclidean, manhattan, wrd, ...\n",
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

DistanceType=DEF_DIST;
if(length(opt$distance_type)){
	DistanceType=opt$distance_type;
}

if(DistanceType=="wrd"){
	source("../../SummaryTableUtilities/WeightedRankDifference.r");
}

###############################################################################

OutputFileRoot=paste(OutputFileRoot, ".", substr(DistanceType, 1,3), sep="");

OutputPDF = paste(OutputFileRoot, ".mds_ts.pdf", sep="");
cat("Output PDF file name: ", OutputPDF, "\n", sep="");
pdf(OutputPDF,width=8.5,height=8.5)

###############################################################################

load_offset=function(fname){

        cat("Loading Offsets: ", fname, "\n");
        offsets_mat=read.delim(fname,  header=TRUE, row.names=1, sep="\t", comment.char="#", quote="");

        num_col=ncol(offsets_mat);
        cat("Num Columns Found: ", num_col, "\n");

        extra_colnames=colnames(offsets_mat);
        print(extra_colnames);
        colnames(offsets_mat)=c("Indiv ID", "Offsets", "Group ID", extra_colnames[4:num_col])[1:num_col];

        # reset offsets
        if(is.numeric(offsets_mat[,"Indiv ID"])){
                offsets_mat[,"Indiv ID"]=paste("#", offsets_mat[,"Indiv ID"], sep="");
        }
        groups=unique(offsets_mat[,"Indiv ID"]);

        cat("Groups:\n");
        print(groups);
        cat("\n");

        # Reset offsets so they are relative to the first/smallest sample
        for(gid in groups){
                group_ix=(gid==offsets_mat[,"Indiv ID"]);
                offsets=offsets_mat[group_ix, "Offsets"];
                min_off=min(offsets);
                offsets_mat[group_ix, "Offsets"]=offsets-min_off;
        }

        offsets_data=list();
        offsets_data[["matrix"]]=offsets_mat;
        offsets_data[["IndivID"]]=extra_colnames[1];
        offsets_data[["Offsets"]]=extra_colnames[2];
        offsets_data[["GroupID"]]=extra_colnames[3];

        return(offsets_data);

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

plot_connected_figure=function(coordinates, offsets_mat, groups_per_plot=3, col_assign, ind_colors, title=""){
	sorted_sids=sort(rownames(offsets_mat));
	coordinates=coordinates[sorted_sids,];
	offsets_mat=offsets_mat[sorted_sids,];

	#print(offsets_mat);
	#print(coordinates);

	# Get Unique Groups
	groups=sort(unique(offsets_mat[,"Indiv ID"]));
	num_groups=length(groups);

	palette(ind_colors);

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
		grp_subset=which(offsets_mat[,"Indiv ID"]==groups[i]);
		num_members=length(grp_subset);
		print(grp_subset);

		offsets_subset=offsets_mat[grp_subset,, drop=F];
		coord_subset=coordinates[grp_subset,, drop=F];

		sort_ix=order(offsets_subset[,"Offsets"], decreasing=F);

		offsets_subset=offsets_subset[sort_ix,, drop=F];
		coord_subset=coord_subset[sort_ix,, drop=F];

		#print(offsets_subset);
		#print(coord_subset);

		#--------------------------------------------------------------------------------
			
		# Draw colored lines
		points(coord_subset, type="l", col=col_assign[groups[i]], pch=20, lwd=2.5);
		# Draw reinforcement black lines
		points(coord_subset, type="b", col="black", pch=20, cex=.1);
		# Draw start/stop glyphs
		points(coord_subset[c(1, 1, num_members),], type="p", col=col_assign[groups[i]], 
			pch=c(17, 1, 15), cex=c(1, 2, 1.25));
	}

	# Plot subset of samples
	for(i in 1:num_groups){
		if(((i-1) %% groups_per_plot)==0){
			plot(0, main=title, xlab="Dim 1", ylab="Dim 2", type="n", xlim=xlim, ylim=ylim);
		}

		#cat("Plotting: ", groups[i], "\n");
		grp_subset=which(offsets_mat[,"Indiv ID"]==groups[i]);
		num_members=length(grp_subset);
		#print(grp_subset);

		offsets_subset=offsets_mat[grp_subset,, drop=F];
		coord_subset=coordinates[grp_subset,, drop=F];

		sort_ix=order(offsets_subset[,"Offsets"], decreasing=F);

		offsets_subset=offsets_subset[sort_ix,, drop=F];
		coord_subset=coord_subset[sort_ix,, drop=F];
			
		#--------------------------------------------------------------------------------
		# Draw colored lines
		points(coord_subset, type="l", col=col_assign[groups[i]], pch=20, cex=.5, lwd=2.5);
		# Draw reinforcement black lines
		points(coord_subset, type="l", col="black", lwd=.1);
		# Draw start/stop glyphs
		points(coord_subset[c(1, 1, num_members),], type="p", col=col_assign[groups[i]], 
			pch=c(17, 1, 15), cex=c(1, 2, 1.25));
		# Label individual id
		text(coord_subset[1,1], coord_subset[1,2], labels=groups[i], col="black", pos=1, cex=.75, font=2);

		# Label offsets
		if(num_members>1){
			offset_ix=2:num_members;
			text(coord_subset[offset_ix,1], coord_subset[offset_ix,2], 
				labels=offsets_subset[offset_ix,"Offsets"], col="black", 
				adj=c(.5,-.75), cex=.5, font=3);
		}
	}
}

###############################################################################

plot_sample_distances=function(distmat, offsets_mat, col_assign, ind_colors, title="", dist_type=""){
	sorted_sids=sort(rownames(offsets_mat));
	offsets_mat=offsets_mat[sorted_sids,];

	# Get Unique Groups
	indiv_ids=sort(unique(offsets_mat[,"Indiv ID"]));
	num_indiv=length(indiv_ids);

	palette(ind_colors);

	def_par=par(no.readonly=T);
	par(mfrow=c(4,1));

	# Get range of offsets
	offset_ranges=range(offsets_mat[,"Offsets"]);
	cat("Offset Range:\n");
	print(offset_ranges);

	#print(offsets_mat);
	distmat2d=as.matrix(distmat);
	dist_ranges=range(distmat2d);
	cat("Distance Ranges:\n");
	print(dist_ranges);

	# Plot subset of samples
	for(i in 1:num_indiv){

		cat("Plotting: ", indiv_ids[i], "\n");
		ind_subset=which(offsets_mat[,"Indiv ID"]==indiv_ids[i]);
		num_samples=length(ind_subset);

		offset_info=offsets_mat[ind_subset,];
		sort_ix=order(offset_info[,"Offsets"]);
		offset_info=offset_info[sort_ix,];
		print(offset_info);

		subset_samples=rownames(offset_info);
		subset_dist=distmat2d[subset_samples[1], subset_samples];
		print(subset_dist);

		# Plot colored lines
		plot(offset_info[,"Offsets"], subset_dist, main=indiv_ids[i],
			 xlab="Time", ylab=paste("Distance (", dist_type, ")", sep=""), type="l", col=col_assign[indiv_ids[i]], lwd=2.5,
			 xlim=offset_ranges, ylim=dist_ranges);
		# Plot ends
		points(offset_info[c(1,1, num_samples),"Offsets"], subset_dist[c(1,1, num_samples)], 
			col=col_assign[indiv_ids[i]],
			type="p", pch=c(17, 1, 15), cex=c(1, 2, 1.25));
		# Plot reinforcement thin black lines
		points(offset_info[,"Offsets"], subset_dist, type="b", pch=16, cex=.1, lwd=.1);
	}
	par(def_par);
}

###############################################################################

plot_sample_dist_by_group=function(dist_mat, offsets_mat, col_assign, ind_colors, dist_type=""){
	
	dist_mat=as.matrix(dist_mat);

        sorted_sids=sort(rownames(offsets_mat));
        offsets_mat=offsets_mat[sorted_sids,, drop=F];

        # Get Num Cohorts
        cohorts=sort(unique(offsets_mat[,"Group ID"]));
        num_cohorts=length(cohorts);
        cat("Number of Cohorts: ", num_cohorts, "\n");
        print(cohorts);
        cat("\n");

        # Get range of offsets
        offset_ranges=range(offsets_mat[,"Offsets"]);
        cat("Offset Range:\n");
        print(offset_ranges);

        # Get range of diversity
        dist_ranges=range(dist_mat);
        cat("Distance Range:\n");
        print(dist_ranges);

        # Set up plots per page
        def_par=par(no.readonly=T);
        par(mfrow=c(num_cohorts,1));

        # Set palette for individuals
        palette(ind_colors);

	x_plot_range=c(offset_ranges[1], offset_ranges[2]+(diff(offset_ranges)/10));

        for(g in 1:num_cohorts){
	
		cat("--------------------------------------------------------------------\n");
                cat("Plotting: ", as.character(cohorts[g]), "\n");
                plot(0, 0, main=cohorts[g],
                         xlab="Time", ylab=paste("Distance (", dist_type, ")", sep=""), type="n",
                         xlim=x_plot_range, ylim=dist_ranges);

                coh_offset_mat=offsets_mat[ offsets_mat[,"Group ID"]==cohorts[g], ];
                print(coh_offset_mat);

                # Get Unique Inidividuals
                indivs=sort(unique(coh_offset_mat[,"Indiv ID"]));
                num_indivs=length(indivs);
                cat("Number of Individuals: ", num_indivs, "\n");
                print(indivs);
                cat("\n");

                # Plot individual samples
                for(i in 1:num_indivs){

                        # Grab from individual cohort 
                        cat("Plotting: ", as.character(indivs[i]), "\n");
                        ind_subset=which(coh_offset_mat[,"Indiv ID"]==indivs[i]);
                        num_timepts=length(ind_subset);

                        # Subset offsets, and sort by offset
                        offset_info=coh_offset_mat[ind_subset,,drop=F];
                        sort_ix=order(offset_info[,"Offsets"]);
                        offset_info=offset_info[sort_ix,];

                        # Subset distances
                        subset_samples=rownames(offset_info);
			subset_dist=dist_mat[subset_samples[1], subset_samples];

                        # Plot distances
                        points(offset_info[c(1,1, num_timepts),"Offsets"], subset_dist[c(1,1, num_timepts)],
                                type="p", pch=c(17, 1, 15), cex=c(1, 2, 1.25), col=col_assign[indivs[i]]);
                        points(offset_info[,"Offsets"], subset_dist, type="l", lwd=2.5, col=col_assign[indivs[i]]);
                        points(offset_info[,"Offsets"], subset_dist, type="l", lwd=.1, col="black");
			text(offset_info[num_timepts, "Offsets"], subset_dist[num_timepts], adj=c(-.5,-1),
				labels=indivs[i], col="black", cex=.5);

                }
        }
        par(def_par);
}

###############################################################################

# Get color assignments
get_colors=function(num_col, alpha=1){
        colors=hsv(seq(0,1,length.out=num_col+1), c(1,.5), c(1,.75,.5), alpha=alpha);
        color_mat_dim=ceiling(sqrt(num_col));
        color_pad=rep("grey", color_mat_dim^2);
        color_pad[1:num_col]=colors[1:num_col];
        color_mat=matrix(color_pad, nrow=color_mat_dim, ncol=color_mat_dim);
        colors=as.vector(t(color_mat));
        colors=colors[colors!="grey"];
}

###############################################################################

offset_data=load_offset(OffsetFileName);
offset_mat=offset_data[["matrix"]];

###############################################################################

counts_mat=load_summary_file(InputFileName);
#print(counts_mat);

###############################################################################

offset_mat_samples=rownames(offset_mat);
counts_mat_samples=rownames(counts_mat);
shared=intersect(offset_mat_samples, counts_mat_samples);

cat("Shared:\n");
print(shared);

cat("\n\n");
cat("Samples not represented in offsets file:\n");
print(setdiff(counts_mat_samples, shared));
cat("Samples not represented in summary table file:\n");
print(setdiff(offset_mat_samples, shared));
cat("\n\n");

offset_mat=offset_mat[shared,];
counts_mat=counts_mat[shared,];

###############################################################################

# Get Cohort info
cohort_names=sort(unique(offset_mat[,"Group ID"]));
num_cohorts=length(cohort_names);
cat("Cohorts:\n");
print(cohort_names);
cat("Num Cohorts: ", num_cohorts, "\n");
cat("\n");

# Get Individuals info
indiv_names=sort(unique(offset_mat[,"Indiv ID"]));
num_indiv=length(indiv_names);
cat("Individuals:\n");
print(indiv_names);
cat("Num Individuals: ", num_indiv, "\n");
cat("\n");

###############################################################################

# Assign colors
ind_colors=get_colors(num_indiv);
col_assign=1:num_indiv;
names(col_assign)=indiv_names;

###############################################################################

normalized_mat=normalize(counts_mat);
#print(normalized_mat);

if(DistanceType=="wrd"){
	dist_mat=weight_rank_dist_opt(normalized_mat, 2);
}else{
	dist_mat=vegdist(normalized_mat, method=DistanceType);
}
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


plot_connected_figure(mds_coord, offset_mat, groups_per_plot=5, col_assign, ind_colors, title=paste("Metric MDS (", DistanceType,")", sep=""));
plot_connected_figure(mds2_coord, offset_mat, groups_per_plot=5, col_assign, ind_colors, title=paste("IsoMetric MDS (", DistanceType, ")", sep=""));

plot_sample_distances(dist_mat, offset_mat, col_assign, ind_colors, dist_type=DistanceType);

plot_sample_dist_by_group(dist_mat, offset_mat, col_assign, ind_colors, dist_type=DistanceType);

##############################################################################

cat("Done.\n")
dev.off();
warn=warnings();
if(length(warn)){
	print(warn);
}
q(status=0)
