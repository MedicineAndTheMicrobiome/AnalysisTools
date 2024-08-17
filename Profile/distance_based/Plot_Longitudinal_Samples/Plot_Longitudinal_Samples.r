#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');

source('~/git/AnalysisTools/Longitudinal/Longitudinal.r');

params=c(
        "summary_file", "s", 1, "character",

        "factor_file", "f", 1, "character",
        "offset_col", "t", 1, "character",
        "subject_id_col", "i", 1, "character",

        "output_root", "o", 1, "character",

        "model_file", "m", 2, "character",
        "group_col", "g", 2, "character",

        "dont_reset_offsets", "n", 2, "logical",
        "begin_offset", "b", 2, "numeric",
        "end_offset", "e", 2, "numeric",

	"distance_type", "d", 2, "character",

	"tag_name", "T", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

DEF_DIST="euclidean";

usage = paste(
	"\nUsage:\n", script_name, "\n",
        "       -s <summary file table>\n",
        "\n",
        "       -f <factor file name>\n",
        "       -t <offset column name>\n",
        "       -i <subject identifier column name>\n",
        "       -o <output root>\n",
        "\n",
        "       [-g <group/cohort variable column name>]\n",
        "\n",
        "       [-n (do not reset earliest offsets to 0 to line up time points, default=reset offsets)]\n",
        "       [-b <begin offset, default=-Inf>]\n",
        "       [-e <end offset, default=Inf>]\n",
	"	[-d <distance type, def=", DEF_DIST, ">]\n",
	"	[-m <model file>]\n",
	"\n",
	"	[-T <tag name>]\n",
	"\n");

if(
        !length(opt$summary_file) ||
        !length(opt$factor_file) ||
        !length(opt$offset_col) ||
        !length(opt$subject_id_col) ||
        !length(opt$output_root)
){
	cat(usage);
	q(status=-1);
}

# Required
SummaryFile=opt$summary_file;
FactorFile=opt$factor_file;
OutputFileRoot=opt$output_root;
SubjectIDCol=opt$subject_id_col;
TimeOffsetCol=opt$offset_col;

# Optional, i.e. with defaults
DistanceType=DEF_DIST;
ResetOffsets=T;
BeginOffset=-Inf;
EndOffset=Inf;
GroupCol="";

if(length(opt$distance_type)){
	DistanceType=opt$distance_type;
}

if(length(opt$dont_reset_offsets)){
        ResetOffsets=F;
}

if(length(opt$begin_offset)){
        BeginOffset=opt$begin_offset;
}

if(length(opt$end_offset)){
        EndOffset=opt$end_offset;
}

if(length(opt$group_col)){
        GroupCol=opt$group_col;
}

if(DistanceType=="wrd"){
	source("../../SummaryTableUtilities/WeightedRankDifference.r");
}

if(length(opt$tag_name)){
        TagName=opt$tag_name;
        cat("Setting TagName Hook: ", TagName, "\n");
        setHook("plot.new",
                function(){
                        #cat("Hook called.\n");
                        if(par()$page==T){
                                oma_orig=par()$oma;
                                exp_oma=oma_orig;
                                exp_oma[1]=max(exp_oma[1], 1.5);
                                par(oma=exp_oma);
                                mtext(paste("[", TagName, "]", sep=""), side=1, line=exp_oma[1]-1,
                                        outer=T, col="steelblue4", font=2, cex=.8, adj=.97);
                                par(oma=oma_orig);
                        }
                },
                "append");

}else{
        TagName="";
}

ModelFile="";
if(length(opt$model_file)){
	ModelFile=opt$model_file;
}

###############################################################################

OutputFileRoot=paste(OutputFileRoot, ".", substr(DistanceType, 1,3), sep="");

OutputPDF = paste(OutputFileRoot, ".mds_ts.pdf", sep="");
cat("Output PDF file name: ", OutputPDF, "\n", sep="");
pdf(OutputPDF,width=8.5,height=8.5)

###############################################################################

load_summary_file=function(fname){
        cat("Loading Summary Table: ", fname, "\n");
        inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", row.names=1))
        counts_mat=inmat[,2:(ncol(inmat))];
        return(counts_mat);
}

load_factors=function(fname){
        factors=data.frame(read.table(fname,  sep="\t", header=TRUE, row.names=1,
                check.names=FALSE, comment.char=""));
        factor_names=colnames(factors);

        ignore_idx=grep("^IGNORE\\.", factor_names);

        if(length(ignore_idx)!=0){
                return(factors[-ignore_idx]);
        }else{
                return(factors);
        }
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

sigchar=function(x){
	if(x<=.0001){
		return("***");
	}else if(x<=.001){
		return("**");
	}else if(x<=.01){
		return("*");
	}else if(x<=.05){
		return(":");
	}else if(x<=.1){
		return(".");
	}else{
		return("");
	}
}

plot_connected_figure=function(coordinates, offsets_rec, subject_grouping_rec, subjects_per_plot=3, 
	col_assign, ind_colors, grp_colors, title=""){

	subject_ids=offsets_rec[["SubjectIDs"]];
	num_subjects=offsets_rec[["NumSubjects"]];
	groups=subject_grouping_rec[["Groups"]];
	num_groups=subject_grouping_rec[["NumGroups"]];

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
	title(main="All Subjects, Colored by Subject", line=.3, cex.main=.75);
	for(i in 1:num_subjects){

		sbj_offset=offsets_rec[["OffsetsBySubject"]][[subject_ids[i]]];
		samp_ids=rownames(sbj_offset);
		num_samples=length(samp_ids);

		coord_subset=coordinates[samp_ids,, drop=F];

		#--------------------------------------------------------------------------------
			
		# Draw colored lines
		points(coord_subset, type="l", col=col_assign[subject_ids[i]], pch=20, lwd=2.5);
		# Draw reinforcement black lines
		points(coord_subset, type="b", col="black", pch=20, cex=.1);
		# Draw start/stop glyphs
		points(coord_subset[c(1, 1, num_samples),], type="p", col=col_assign[subject_ids[i]], 
			pch=c(17, 1, 15), cex=c(1, 2, 1.25));
	}

	# compute smoothed mds points
	smoothed_sbj_mds=list();
	for(i in 1:num_subjects){

		sbj_offset=offsets_rec[["OffsetsBySubject"]][[subject_ids[i]]];
		samp_ids=rownames(sbj_offset);
		num_samples=length(samp_ids);

		coord_subset=coordinates[samp_ids,, drop=F];
		num_coords=nrow(coord_subset);
		smoothed_x=ksmooth(1:num_coords, coord_subset[,1], kernel="normal", bandwidth=1.2);
		smoothed_y=ksmooth(1:num_coords, coord_subset[,2], kernel="normal", bandwidth=1.2);

		# Append first/last point, the smoothed coorindates never reach begin/end points
		lastx=coord_subset[num_coords,1];
		lasty=coord_subset[num_coords,2];
		firstx=coord_subset[1,1];
		firsty=coord_subset[1,2];

		sm_crds=cbind(
			c(firstx, smoothed_x$y, lastx),
			c(firsty, smoothed_y$y, lasty)
		);
		smoothed_sbj_mds[[subject_ids[i]]]=sm_crds;
	}

	# Plot individually colored smoothed individuals
	plot(0, main=title, xlab="Dim 1", ylab="Dim 2", type="n", xlim=xlim, ylim=ylim);
	title(main="All Subjects, Colored by Subject: Lightly Smoothed", line=.3, cex.main=.75);
	for(i in 1:num_subjects){

		cur_sbj=subject_ids[i];
		cur_sm_coords=smoothed_sbj_mds[[cur_sbj]];
		cur_num_coords=nrow(cur_sm_coords);

		# Draw colored lines
		points(cur_sm_coords[,1], cur_sm_coords[,2],
			type="l", col=col_assign[subject_ids[i]], pch=20, lwd=2.5);
		# Draw reinforcement black lines
		points(cur_sm_coords[,1], cur_sm_coords[,2],
			type="l", col="black", pch=20, cex=.1);

		# Draw start/stop glyphs
		points(cur_sm_coords[c(1, 1, cur_num_coords),], type="p", col=col_assign[cur_sbj], 
			pch=c(17, 1, 15), cex=c(1, 2, 1.25));
	}

	# Plot with group colors
	plot(0, main=title, xlab="Dim 1", ylab="Dim 2", type="n", xlim=xlim, ylim=ylim);
	title(main="All Subjects, Colored by Group: Lightly Smoothed", line=.3, cex.main=.75);
	for(i in 1:num_subjects){

		cur_sbj=subject_ids[i];
		cur_sm_coords=smoothed_sbj_mds[[cur_sbj]];
		cur_num_coords=nrow(cur_sm_coords);

		cur_grp_col=grp_colors[subject_grouping_rec[["SbjToGrp"]][[cur_sbj]]];

		# Draw colored lines
		points(cur_sm_coords[,1], cur_sm_coords[,2],
			type="l", col=cur_grp_col, pch=20, lwd=2.5);
		# Draw reinforcement black lines
		points(cur_sm_coords[,1], cur_sm_coords[,2],
			type="l", col="black", pch=20, cex=.1);

		# Draw start/stop glyphs
		points(cur_sm_coords[c(1, 1, cur_num_coords),], type="p", col=cur_grp_col, 
			pch=c(17, 1, 15), cex=c(1, 2, 1.25));
	}
	legend(xlim[1], ylim[2], fill=grp_colors, legend=names(grp_colors), bty="n");

	# Plot subset of samples
	for(i in 1:num_subjects){

		cur_sbj=subject_ids[i]

		cat("Subject: ", subject_ids[i], "\n");	
		sbj_offset=offsets_rec[["OffsetsBySubject"]][[subject_ids[i]]];
		samp_ids=rownames(sbj_offset);
                num_samples=length(samp_ids);
		coord_subset=coordinates[samp_ids,, drop=F];

		plot(0, main=title, xlab="Dim 1", ylab="Dim 2", type="n", xlim=xlim, ylim=ylim);
		title(main=paste("By Subject: ", cur_sbj, ", Num Time Pts Labeled: ", num_samples, sep=""), 
			line=.3, cex.main=.75);

		# Draw colored lines
		points(coord_subset, type="l", col=col_assign[subject_ids[i]], pch=20, cex=.5, lwd=2.5);
		# Draw reinforcement black lines
		points(coord_subset, type="l", col="black", lwd=.1);
		# Draw start/stop glyphs
		points(coord_subset[c(1, 1, num_samples),], type="p", col=col_assign[subject_ids[i]], 
			pch=c(17, 1, 15), cex=c(1, 2, 1.25));
		# Label individual id
		text(coord_subset[1,1], coord_subset[1,2], labels=subject_ids[i], 
			col="black", pos=1, cex=.75, font=2);

		# Label offsets
		if(num_samples>1){
			offset_ix=2:num_samples;
			text(coord_subset[offset_ix,1], coord_subset[offset_ix,2], 
				labels=sbj_offset[offset_ix,"Offsets"], col="black", 
				adj=c(.5,-.75), cex=.5, font=3);
		}
	}
}

###############################################################################

plot_sample_distances=function(distmat, offsets_rec, subject_grouping_rec, col_assign, 
	ind_colors, title="", dist_type=""){

	# Get Unique Groups
	subject_ids=offsets_rec[["SubjectIDs"]];
	num_subjects=offsets_rec[["NumSubjects"]];

	palette(ind_colors);

	def_par=par(no.readonly=T);
	par(mfrow=c(4,1));

	# Get range of offsets
	offset_ranges=range(offsets_rec[["Offsets"]]);

	distmat2d=as.matrix(distmat);
	dist_ranges=range(distmat2d);
	cat("Distance Ranges:\n");
	print(dist_ranges);

	# Plot subset of samples
	for(i in 1:num_subjects){

		cat("Plotting: ", subject_ids[i], "\n");
		cur_sbj=subject_ids[i];
		offset_info=offsets_rec[["OffsetsBySubject"]][[cur_sbj]];
		sample_ids=rownames(offset_info);
		num_samples=length(sample_ids);

		subset_dist=distmat2d[sample_ids[1], sample_ids];
		print(subset_dist);
		dist_ranges=range(subset_dist);

		y_pad=diff(dist_ranges)*.15;

		# Plot colored lines
		plot(offset_info[,"Offsets"], subset_dist, main=subject_ids[i],
			xlab="Time", ylab=paste("Distance (", dist_type, ")", sep=""), 
			type="l", col=col_assign[subject_ids[i]], lwd=2.5,
			xaxt="n",
			xlim=offset_ranges, ylim=c(dist_ranges[1]-y_pad, dist_ranges[2]+y_pad));

		# accentuation where points are
		if(num_samples>2){
			middle_points_ix=2:(num_samples-1);
			points(offset_info[middle_points_ix,"Offsets"], subset_dist[middle_points_ix],
				type="p", pch=20, cex=2);
		}

		# Plot ends
		points(offset_info[c(1,1, num_samples),"Offsets"], subset_dist[c(1,1, num_samples)], 
			col=col_assign[subject_ids[i]],
			type="p", pch=c(17, 1, 15), cex=c(2, 4, 2.5));
		# Plot reinforcement thin black lines
		points(offset_info[,"Offsets"], subset_dist, type="l", pch=16, cex=.1, lwd=.1);

		# Generate axis
		tick_pos=sort(unique(c(offset_ranges, offset_info[,"Offsets"])));
		axis(1, at=tick_pos);
	}
	par(def_par);
}

###############################################################################

plot_sample_dist_by_group=function(dist_mat, offsets_rec, subject_grouping_rec, col_assign, 
	ind_colors, dist_type=""){

	dist_mat=as.matrix(dist_mat);
	subjects_ids=offsets_rec[["SubjectIDs"]];
	num_subjects=offsets_rec[["NumSubjects"]];
	group_ids=subject_grouping_rec[["Groups"]];
	num_groups=subject_grouping_rec[["NumGroups"]];
	offset_ranges=range(offsets_rec[["Offsets"]]);

        # Get range of diversity
        dist_ranges=range(dist_mat);
        cat("Distance Range:\n");
        print(dist_ranges);

        # Set up plots per page
        def_par=par(no.readonly=T);
        par(mfrow=c(num_groups,1));

        # Set palette for individuals
        palette(ind_colors);

	x_plot_range=c(offset_ranges[1], offset_ranges[2]+(diff(offset_ranges)/10));

        for(g in 1:num_groups){
	
		cat("--------------------------------------------------------------------\n");
                cat("Plotting: ", as.character(group_ids[g]), "\n");
                plot(0, 0, main=group_ids[g],
                         xlab="Time", ylab=paste("Distance (", dist_type, ")", sep=""), type="n",
                         xlim=x_plot_range, ylim=dist_ranges);

		cur_grp=group_ids[g];
		subj_in_grp=subject_grouping_rec[["GrpToSbj"]][[cur_grp]];
		num_sbj_in_grp=length(subj_in_grp);

                # Plot individual samples
                for(i in 1:num_sbj_in_grp){

                        # Grab from individual cohort 
                        cat("Plotting: ", as.character(subj_in_grp[i]), "\n");
			offset_info=offsets_rec[["OffsetsBySubject"]][[subj_in_grp[i]]];
                        num_timepts=nrow(offset_info);

                        # Subset distances
                        sample_ids=rownames(offset_info);
			subset_dist=dist_mat[sample_ids[1], sample_ids];

                        # Plot distances
                        points(offset_info[c(1,1, num_timepts),"Offsets"], subset_dist[c(1,1, num_timepts)],
                                type="p", pch=c(17, 1, 15), cex=c(1, 2, 1.25), col=col_assign[subj_in_grp[i]]);
                        points(offset_info[,"Offsets"], subset_dist, type="l", lwd=2.5,
				 col=col_assign[subj_in_grp[i]]);
                        points(offset_info[,"Offsets"], subset_dist, type="l", lwd=.1, col="black");
			text(offset_info[num_timepts, "Offsets"], subset_dist[num_timepts], adj=c(-.5,-1),
				labels=subj_in_grp[i], col="black", cex=.5);

                }
        }
        par(def_par);
}

###############################################################################

plot_sample_dist_by_group_loess=function(dist_mat, offsets_rec, subject_grouping_rec, col_assign, 
	ind_colors, dist_type=""){

	dist_mat=as.matrix(dist_mat);
	group_ids=subject_grouping_rec[["Groups"]];
	num_groups=subject_grouping_rec[["NumGroups"]];
	offset_ranges=range(offsets_rec[["Offsets"]]);
	num_offsets=offsets_rec[["NumUniqOffsets"]];

        # Get range of diversity
        dist_ranges=range(dist_mat);
        cat("Distance Range:\n");
        print(dist_ranges);

        # Set up plots per page
        def_par=par(no.readonly=T);
        par(mfrow=c(num_groups,1));

        # Set palette for individuals
        palette(ind_colors);

	x_plot_range=c(offset_ranges[1], offset_ranges[2]+(diff(offset_ranges)/10));

	lowess_rec=list();
	subj_points_rec=list();
	max_subj_dist=0;
	max_loess=0;
        for(g in 1:num_groups){
	
		# Get subjects in group
		cur_grp=as.character(group_ids[g]);
		subj_in_grp=as.character(subject_grouping_rec[["GrpToSbj"]][[cur_grp]]);
		num_sbj_in_grp=length(subj_in_grp);

		# Prep for loess
		group_xs=c();
		group_ys=c();
		for(i in 1:num_sbj_in_grp){

			offset_info=offsets_rec[["OffsetsBySubject"]][[subj_in_grp[i]]];
                        sample_ids=rownames(offset_info);
			subset_dist=dist_mat[sample_ids[1], sample_ids];
			time_vals=offset_info[,"Offsets"];
		
			group_xs=c(group_xs, time_vals);
			group_ys=c(group_ys, subset_dist);
			subj_points_rec[[subj_in_grp[i]]]=cbind(time_vals, subset_dist);

			max_subj_dist=max(max_subj_dist, subset_dist);
		}

		# Compute loess for group
		lowess_res=lowess(group_xs, group_ys, f=2.75/10);
			
		max_lowess=max(max_loess, lowess_res$y, na.rm=T);
		lowess_rec[[cur_grp]]=cbind(lowess_res$x, lowess_res$y);

	}

	#print(loess_rec);
	#print(subj_points_rec);

	# Create plot
	plot(0,0, xlim=offset_ranges, ylim=c(0, max(max_lowess, max_subj_dist, na.rm=T)), 
		main="Lowess by Group", xlab="Time", ylab="Distance",);

	# Group colors
	groups=names(lowess_rec);
	grp_colors=terrain.colors(num_groups+1)[1:num_groups];
	names(grp_colors)=groups;

	# Plot scatter for all 
	subjects=names(subj_points_rec);
	for(sbj_ix in subjects){
		points(subj_points_rec[[sbj_ix]][,1], subj_points_rec[[sbj_ix]][,2],
			col=grp_colors[subject_grouping_rec[["SbjToGrp"]][sbj_ix]]);	
			#col=col_assign[sbj_ix]);	
	} 

	# Plot loess lines
	for(grp_ix in groups){
		cur_lowess=lowess_rec[[grp_ix]];
		points(cur_lowess[,1], cur_lowess[,2], type="l", lwd=2, col=grp_colors[grp_ix]);
	}	

	# Plot legend
	plot(0,0, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", main="", type="n", bty="n");
	legend(0,1, fill=grp_colors, legend=groups);

        par(def_par);

#quit();
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

plot_text=function(strings){
        par(family="Courier");
        par(oma=rep(.1,4));
        par(mar=rep(0,4));

        num_lines=length(strings);

        top=max(as.integer(num_lines), 52);

        plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                );

        text_size=max(.01, min(.8, .8 - .003*(num_lines-52)));
        #print(text_size);

        for(i in 1:num_lines){
                #cat(strings[i], "\n", sep="");
                strings[i]=gsub("\t", "", strings[i]);
                text(0, top-i, strings[i], pos=4, cex=text_size);
        }
}

###############################################################################

export_distances_from_start=function(fname, offset_rec, dist_mat){

	cat("Exporting distances from start into: ", fname, "\n", sep="");
	
	subject_ids=offset_rec[["SubjectIDs"]];
	num_subjects=length(subject_ids);
	dist_mat=as.matrix(dist_mat);

	fh=file(fname, "w");

	for(cur_id in subject_ids){
		
		cur_offsets=offset_rec[["OffsetsBySubject"]][[cur_id]];
		sample_ids=rownames(cur_offsets);

		# Grab distances from first sample
		cur_dist=dist_mat[sample_ids[1], sample_ids];
		cur_times=cur_offsets[,"Offsets"];

		cat(file=fh, "Subject:,", cur_id, "\n", sep="");
		cat(file=fh, "Offset:,", paste(cur_times, collapse=","), "\n", sep="");
		cat(file=fh, "Sample ID:,", paste(sample_ids, collapse=","), "\n", sep="");
		cat(file=fh, "Distance:,", paste(round(cur_dist,4), collapse=","), "\n", sep="");
		cat(file=fh, "\n");

	}
	
	close(fh);

}

###############################################################################

# Load Factor File
cat("Loading Factor File:\n");
factor_info=load_factors(FactorFile);
factor_samp_ids=rownames(factor_info);

###############################################################################

counts_mat=load_summary_file(SummaryFile);

###############################################################################

factor_mat_samples=rownames(factor_info);
counts_mat_samples=rownames(counts_mat);
shared_samp_ids=intersect(factor_mat_samples, counts_mat_samples);

cat("Shared:\n");
print(shared_samp_ids);

cat("\n\n");
cat("Samples not represented in offsets file:\n");
print(setdiff(counts_mat_samples, shared_samp_ids));
cat("Samples not represented in summary table file:\n");
print(setdiff(factor_mat_samples, shared_samp_ids));
cat("\n\n");

counts_mat=counts_mat[shared_samp_ids,,drop=F];
factor_info=factor_info[shared_samp_ids,,drop=F];

# Set up grouping information
if(GroupCol==""){
        group_names=rep("All", nrow(factor_info));
        GroupCol="Group";
}else{
        group_names=factor_info[,GroupCol];
}
subject_grouping_rec=create_GrpToSbj_map(factor_info[,SubjectIDCol], group_names);
unique_group_names=sort(unique(group_names));

# Extract offset info
offset_rec=extract_offset(factor_info, SubjectIDCol, TimeOffsetCol);

print(offset_rec);
print(subject_grouping_rec);

###############################################################################

groups=subject_grouping_rec[["Groups"]];
num_groups=subject_grouping_rec[["NumGroups"]];
subject_ids=offset_rec[["SubjectIDs"]];
num_subjects=offset_rec[["NumSubjects"]];

###############################################################################

# Assign colors
ind_colors=get_colors(num_subjects);
col_assign=1:num_subjects;
names(col_assign)=subject_ids;

group_colors=terrain.colors(num_groups+1)[1:num_groups];
names(group_colors)=groups;

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

long_stats=calculate_stats_on_series_distance(offset_rec, dist_mat);

plot_connected_figure(mds_coord, offset_rec, subject_grouping_rec, 
	subjects_per_plot=as.integer(num_subjects/5), col_assign, 
	ind_colors, group_colors, title=paste("Metric MDS (", DistanceType,")", sep=""));

plot_connected_figure(mds2_coord, offset_rec, subject_grouping_rec, 
	subjects_per_plot=as.integer(num_subjects/5), col_assign, 
	ind_colors, group_colors, title=paste("IsoMetric MDS (", DistanceType, ")", sep=""));

plot_sample_distances(dist_mat, offset_rec, subject_grouping_rec, col_assign, ind_colors, 
	dist_type=DistanceType);

plot_sample_dist_by_group(dist_mat, offset_rec, subject_grouping_rec, col_assign, ind_colors, 
	dist_type=DistanceType);

plot_sample_dist_by_group_loess(dist_mat, offset_rec, subject_grouping_rec, col_assign, ind_colors, 
	dist_type=DistanceType);

group_stat_comparisons=plot_pairwise_grp_comparisons(long_stats, subject_grouping_rec, plots_pp=1);

##############################################################################

output_stat_table_alternate_ordering(group_stat_comparisons, OutputFileRoot);

##############################################################################

export_distances_from_start(paste(OutputFileRoot,".dist_from_start.tsv", sep=""),
	offset_rec, dist_mat);

##############################################################################
##############################################################################

load_list=function(filename){
        val=scan(filename, what=character(), comment.char="#");
        return(val);
}

if(ModelFile!=""){
        model_var_list=load_list(ModelFile);
        model_var_list=c(model_var_list, GroupCol);
        cat("Model variables in: ", ModelFile, "\n");
}else{
        cat("Model File was not specified. Skipping analyses with factors.\n");
        quit();
}


cat("Collapsing Factors...\n");
colpsd_factors=collapse_factors(factor_info, SubjectIDCol, model_var_list);

cat("Regressing Longitudinal Stats...\n");
regres=regress_longitudinal_stats(long_stats, model_var_list, colpsd_factors);

options(width=200);
print(regres);

cat("Generating Heatmaps by Stat...\n");
title_page(paste(
        "Response:\nLongitudinal Stats\n\nPredictors:\n", paste(model_var_list, collapse="\n"), sep=""));

stat_names=names(regres);
for(stat_ix in stat_names){
        plot_heat_maps(
                regres[[stat_ix]][["coef"]],
                regres[[stat_ix]][["pval"]],
                stat_ix);
}

cat("Summarizing Regression Results into Table...\n");
regr_stat_summary=summarize_regression_results(regres, stat_names);
num_sigf_reg_assoc=nrow(regr_stat_summary);

cat("Writing Stats by Alternative Ordering...\n");
output_long_regression_stats_w_alt_ordering(regr_stat_summary);

##############################################################################

par(mfrow=c(1,1));
plot_text(longit_stat_description_distance);

##############################################################################

cat("Done.\n")
dev.off();
warn=warnings();
if(length(warn)){
	print(warn);
}
q(status=0)
