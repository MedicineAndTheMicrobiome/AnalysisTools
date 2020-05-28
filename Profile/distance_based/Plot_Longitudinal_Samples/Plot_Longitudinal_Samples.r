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

	"distance_type", "d", 2, "character"
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

plot_connected_figure=function(coordinates, offsets_rec, subject_grouping_rec, subjects_per_plot=3, 
	col_assign, ind_colors, title=""){

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
	title(main="All Subjects", line=.3, cex.main=.75);
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

	# Plot subset of samples
	for(i in 1:num_subjects){
		if(((i-1) %% subjects_per_plot)==0){
			cat("Opening new plot...\n");
			plot(0, main=title, xlab="Dim 1", ylab="Dim 2", type="n", xlim=xlim, ylim=ylim);
			title(main="(Subset of Subjects)", line=.3, cex.main=.75);
		}

		cat("Subject: ", subject_ids[i], "\n");	
		sbj_offset=offsets_rec[["OffsetsBySubject"]][[subject_ids[i]]];
		samp_ids=rownames(sbj_offset);
                num_samples=length(samp_ids);
		coord_subset=coordinates[samp_ids,, drop=F];

		#--------------------------------------------------------------------------------
		# Draw colored lines
		points(coord_subset, type="l", col=col_assign[subject_ids[i]], pch=20, cex=.5, lwd=2.5);
		# Draw reinforcement black lines
		points(coord_subset, type="l", col="black", lwd=.1);
		# Draw start/stop glyphs
		points(coord_subset[c(1, 1, num_samples),], type="p", col=col_assign[subject_ids[i]], 
			pch=c(17, 1, 15), cex=c(1, 2, 1.25));
		# Label individual id
		text(coord_subset[1,1], coord_subset[1,2], labels=subject_ids[i], col="black", pos=1, cex=.75, font=2);

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

	loess_rec=list();
	subj_points_rec=list();
	max_subj_dist=0;
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
		loess_res=loess(group_ys~group_xs);
			
		loess_fit_x=seq(offset_ranges[1], offset_ranges[2], length.out=num_offsets*2);
		loess_fit_y=tryCatch({
				predict(loess_res, loess_fit_x);
			}, error=function(e){
				return(NULL);
			});

		if(is.null(loess_fit_y)){
			loess_fit_x=time_vals;
			loess_fit_y=subset_dist;
		}

		loess_rec[[cur_grp]]=cbind(loess_fit_x, loess_fit_y);

	}

	#print(loess_rec);
	#print(subj_points_rec);

	# Create plot
	plot(0,0, xlim=offset_ranges, ylim=c(0, max_subj_dist), xlab="Time", ylab="Distance");

	# Plot scatter for all 
	subjects=names(subj_points_rec);
	for(sbj_ix in subjects){
		points(subj_points_rec[[sbj_ix]][,1], subj_points_rec[[sbj_ix]][,2],
			col=col_assign[sbj_ix]);	
	} 

	# Plot loess lines
	groups=names(loess_rec);
	grp_colors=terrain.colors(num_groups+2)[2:(num_groups-1)];
	names(grp_colors)=groups;

	for(grp_ix in groups){
		cur_loess=loess_rec[[grp_ix]];
		points(cur_loess[,1], cur_loess[,2], type="l", lwd=2, col=grp_colors[grp_ix]);
	}	

	# Plot legend
	plot(0,0, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", main="", type="n", bty="n");
	legend(0,1, fill=grp_colors, legend=groups);

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

calculate_stats_on_series=function(offset_rec, dist_mat){

	avg_dist=function(dist_arr, time_arr){

		# Average distance sample spent away from home (0)
		num_pts=length(dist_arr);

		acc_dist=0;
		for(i in 1:(num_pts-1)){
			avg_dist=(dist_arr[i+1]+dist_arr[i])/2;
			dtime=(time_arr[i+1]-time_arr[i]);
			acc_dist=acc_dist+avg_dist*dtime;	
		}
		
		overall_avg_dist=acc_dist/time_arr[num_pts];
		return(overall_avg_dist);
	}

	avg_speed=function(dist_arr, time_arr){
		# total distance traveled divided by time
		num_pts=length(dist_arr);
		acc_dist=0;
		for(i in 1:(num_pts-1)){
			ddist=abs(dist_arr[i+1]-dist_arr[i]);
			acc_dist=acc_dist+ddist;
		}
		average_speed=acc_dist/time_arr[num_pts];
		return(average_speed);
	}

	tot_dist_travelled=function(dist_arr, time_arr){
		# total distance traveled 
		num_pts=length(dist_arr);
		acc_dist=0;
		for(i in 1:(num_pts-1)){
			ddist=abs(dist_arr[i+1]-dist_arr[i]);
			acc_dist=acc_dist+ddist;
		}
		return(acc_dist);
	}

	mean_reversion=function(dist_arr, time_arr){
		fit=lm(dist_arr~time_arr);
		res=list();
		res[["first_dist"]]=fit$coefficients[["(Intercept)"]];
		res[["slope"]]=fit$coefficients[["time_arr"]];
		res[["last_dist"]]=res[["first_dist"]]+res[["slope"]]*tail(time_arr,1);
		res[["sd_res"]]=sd(fit$residuals);
		return(res);
	}

	closest_travel=function(dist_arr, time_arr){
		
		dist_arr=dist_arr[-1];
		time_arr=time_arr[-1];

		min_dist=min(dist_arr);
		ix=min(which(min_dist==dist_arr));
		
		res=list();
		res[["dist"]]=min_dist;
		res[["time"]]=time_arr[ix];
		return(res);
	}

	furthest_travel=function(dist_arr, time_arr){
		
		dist_arr=dist_arr[-1];
		time_arr=time_arr[-1];

		max_dist=max(dist_arr);
		ix=min(which(max_dist==dist_arr));
		
		res=list();
		res[["dist"]]=max_dist;
		res[["time"]]=time_arr[ix];
		return(res);
	}

	closest_return=function(dist_arr, time_arr){

		while(length(dist_arr)>1 && dist_arr[1]<=dist_arr[2]){
			dist_arr=dist_arr[-1];
			time_arr=time_arr[-1];
		}
		dist_arr=dist_arr[-1];
		time_arr=time_arr[-1];

                res=list();
		if(length(dist_arr)){
			min_dist=min(dist_arr);
			ix=min(which(min_dist==dist_arr));
			res[["dist"]]=min_dist;
			res[["time"]]=time_arr[ix];
		}else{
			res[["dist"]]=NA;
			res[["time"]]=NA;
		}

                return(res);
	}

	first_return=function(dist_arr, time_arr){

		while(length(dist_arr)>1 && dist_arr[1]<=dist_arr[2]){
			dist_arr=dist_arr[-1];
			time_arr=time_arr[-1];
		}
		dist_arr=dist_arr[-1];
		time_arr=time_arr[-1];

                res=list();
		if(length(dist_arr)){
			res[["dist"]]=dist_arr[1];
			res[["time"]]=time_arr[1];
		}else{
			res[["dist"]]=NA;
			res[["time"]]=NA;
		}

                return(res);
	}


	cat("Calculating average distance over time...\n");

	uniq_indiv_ids=offset_rec[["SubjectIDs"]];
	num_ind=offset_rec[["NumSubjects"]];
	
	stat_names=c(
		"last_time", "num_time_pts",
		"average_dist", 
		"average_speed",
		"total_dist_travelled",
		"mean_reversion_first_dist", "mean_reversion_last_dist", 
		"mean_reversion_stdev_residuals", "mean_reversion_slope",
		"closest_travel_dist", "closest_travel_time",
		"furthest_travel_dist", "furthest_travel_time",
		"closest_return_dist", "closest_return_time",
		"first_return_dist", "first_return_time");

	out_mat=matrix(NA, nrow=num_ind, ncol=length(stat_names));
	rownames(out_mat)=uniq_indiv_ids;
	colnames(out_mat)=stat_names;

	dist_mat=as.matrix(dist_mat);

	for(cur_id in uniq_indiv_ids){
		
		cur_offsets=offset_rec[["OffsetsBySubject"]][[cur_id]];
		num_timepts=nrow(cur_offsets);

		out_mat[cur_id, "last_time"]=cur_offsets[num_timepts, "Offsets"];
		out_mat[cur_id, "num_time_pts"]=num_timepts;

		samp_ids=rownames(cur_offsets);

		if(num_timepts>1){
			cur_dist=dist_mat[samp_ids[1], samp_ids];
			cur_times=cur_offsets[,"Offsets"];

			out_mat[cur_id, "average_dist"]=avg_dist(cur_dist, cur_times);
			out_mat[cur_id, "average_speed"]=avg_speed(cur_dist, cur_times);
			out_mat[cur_id, "total_dist_travelled"]=tot_dist_travelled(cur_dist, cur_times);

			res=mean_reversion(cur_dist, cur_times);
			out_mat[cur_id, "mean_reversion_first_dist"]=res[["first_dist"]];
			out_mat[cur_id, "mean_reversion_last_dist"]=res[["last_dist"]];
			out_mat[cur_id, "mean_reversion_stdev_residuals"]=res[["sd_res"]];
			out_mat[cur_id, "mean_reversion_slope"]=res[["slope"]];

			res=closest_travel(cur_dist, cur_times);
			out_mat[cur_id, "closest_travel_dist"]=res[["dist"]];
			out_mat[cur_id, "closest_travel_time"]=res[["time"]];

			res=furthest_travel(cur_dist, cur_times);
			out_mat[cur_id, "furthest_travel_dist"]=res[["dist"]];
			out_mat[cur_id, "furthest_travel_time"]=res[["time"]];

			res=closest_return(cur_dist, cur_times);
			out_mat[cur_id, "closest_return_dist"]=res[["dist"]];
			out_mat[cur_id, "closest_return_time"]=res[["time"]];

			res=first_return(cur_dist, cur_times);
			out_mat[cur_id, "first_return_dist"]=res[["dist"]];
			out_mat[cur_id, "first_return_time"]=res[["time"]];
		}	
	}

	return(out_mat);	
}

###############################################################################

plot_barplot_wsignf_annot=function(title, stat, grps, alpha=0.05, samp_gly=T){
	# Generate a barplot based on stats and groupings
	# Annotat barplot with signficance

	cat("Making Barplot with Significance annotated...\n");
        cat("  Alpha", alpha, "\n");
        group_names=names(grps);
        num_grps=length(group_names);

	# Convert matrix into array, if necessary
	if(!is.null(dim(stat))){
		stat_name=colnames(stat);
		stat=stat[,1];
	}else{
		stat_name="value";
	}

	# Remove NAs
	na_ix=is.na(stat);
	subj=names(stat);
	stat=stat[!na_ix];
	na_subj=names(stat);
	for(grnm in group_names){
		grps[[grnm]]=intersect(grps[[grnm]], na_subj);
		print(stat[grps[[grnm]]]);
	}
	print(grps);

        # Precompute pairwise wilcoxon pvalues
	cat("\n  Precomputing group pairwise p-values...\n");
        pval_mat=matrix(1, nrow=num_grps, ncol=num_grps);
	rownames(pval_mat)=group_names;
	colnames(pval_mat)=group_names;
        signf=numeric();
        for(grp_ix_A in 1:num_grps){
                for(grp_ix_B in 1:num_grps){
                        if(grp_ix_A<grp_ix_B){

				grpAnm=group_names[grp_ix_A];
				grpBnm=group_names[grp_ix_B];

                                res=wilcox.test(stat[grps[[grpAnm]]], stat[grps[[grpBnm]]]);

				if(is.na(res$p.value)){
					pval_mat[grpAnm, grpBnm]=1;
				}else{
					pval_mat[grpAnm, grpBnm]=res$p.value;
				}
                                if(pval_mat[grpAnm, grpBnm]<=alpha){
                                        signf=rbind(signf, c(grpAnm, grpBnm, res$p.value));
                                }
                        }
                }
        }

	cat("p-value matrix:\n");
	print(pval_mat);

	# Count how many rows have significant pairings
        num_signf=nrow(signf);
        cat("  Num Significant: ", num_signf, "\n");
        signf_by_row=apply(pval_mat, 1, function(x){sum(x<alpha)});
        cat("  Num Significant by Row:\n");
        print(signf_by_row);

        num_signf_rows=sum(signf_by_row>0);
        cat("  Num Rows to plot:", num_signf_rows, "\n");

        #signf_mat=apply(pval_mat, 1:2,
        #       function(x){
        #               if(x<.001){return("***")}
        #               if(x<.01){return("**")}
        #               if(x<.05){return("*")}
        #               else{return("")}
        #       }
        #);

        #print(signf_mat, quote=F);

        # Compute 95% CI around mean
	cat("\n  Precomputing group means and 95% CI...\n");
        num_bs=320;

        grp_means=numeric(num_grps);
	names(grp_means)=group_names;

        ci95=matrix(NA, nrow=num_grps, ncol=2);
	rownames(ci95)=group_names;
	colnames(ci95)=c("LB", "UB");
        samp_size=numeric(num_grps);
        for(grp_ix in 1:num_grps){

		grpnm=group_names[grp_ix];
                grp_means[grpnm]=mean(stat[grps[[grpnm]]]);
                num_samp=length(grps[[grpnm]]);

                if(num_samp>=40){
                        meds=numeric(num_bs);
                        for(i in 1:num_bs){
                                meds[i]=mean(sample(stat[grps[[grpnm]]], replace=T));

                        }
                        ci95[grp_ix,]=quantile(meds, c(.025, .975));
                }else{
                        ci95[grp_ix,]=rep(mean(stat[grps[[grpnm]]]),2);
                }

                samp_size[grp_ix]=num_samp;
        }

        cat("Group Means:\n");
        print(grp_means);
	print(length(grp_means));
        cat("Group Median 95% CI:\n");
        print(ci95);

        # Estimate spacing for annotations
        annot_line_prop=1/5; # proportion of pl
        min_95ci=min(c(ci95[,1], stat), na.rm=T);
        max_95ci=max(c(ci95[,2], stat), na.rm=T);
	minmax_span=max_95ci-min_95ci;
        plotdatamax=max_95ci+minmax_span*0.3;
	plotdatamin=min_95ci-minmax_span*0.3;;
        space_for_annotations=minmax_span*annot_line_prop*(num_signf_rows+2);
        horiz_spacing=annot_line_prop*plotdatamax;

        # Start plot
        par(mar=c(8,5,4,3));
	cat("  Plot Limits: (", plotdatamin, ", ", plotdatamax, ")\n"); 
	plot(0, type="n", 
		ylim=c(plotdatamin, plotdatamax+space_for_annotations),
		xlim=c(0, num_grps+1),
                yaxt="n", xaxt="n", xlab="", ylab="", bty="n");
	for(grp_ix in 1:num_grps){
		points(c(grp_ix-.25, grp_ix+.25), rep(grp_means[grp_ix],2), type="l", lwd=3);
	}
	mids=1:num_grps;

	yticks=seq(min_95ci, max_95ci, length.out=5);
	cat("Y ticks:\n");
	print(yticks);
	signf_digits=max(ceiling(abs(log10(abs(yticks)))));
	yticks=signif(yticks, signf_digits);
	
	axis(side=2, at=yticks, labels=sprintf("%3.2f", yticks), cex.axis=.5, las=2);
        title(ylab=paste("Mean ", stat_name, "\nwith Bootstrapped 95% CI", sep=""));
        title(main=title, cex.main=1.5);
        title(main="with Wilcoxon rank sum test (difference between group means) p-values",
                line=.25, cex.main=.7, font.main=3);

        bar_width=mean(diff(mids));
	if(is.na(bar_width)){
		bar_width=.5;
	}
        qbw=bar_width/6;

	# Label x-axis
        text(mids-par()$cxy[1]/2, rep(6*-par()$cxy[2]/2, num_grps),
                group_names, srt=-45, xpd=T, pos=4,
                cex=min(c(1,.7*bar_width/par()$cxy[1])));

        # Scatter
        if(samp_gly){
		cat("Plotting point scatter...\n");

                for(grp_ix in 1:num_grps){
			grpnm=group_names[grp_ix];
                        pts=stat[grps[[grpnm]]];
                        numpts=length(pts);
                        points(
                                #rep(mids[grp_ix], numpts),
                                mids[grp_ix]+rnorm(numpts, 0, bar_width/10),
                                pts, col="darkblue", cex=.5, type="p");
                }
        }

        # label CI's
        for(grp_ix in 1:num_grps){
                if(samp_size[grp_ix]>=40){
			cat("Plotting 95% CI Lines...\n");

                        points(
                                c(mids[grp_ix]-qbw, mids[grp_ix]+qbw),
                                rep(ci95[grp_ix, 2],2), type="l", col="blue");
                        points(
                                c(mids[grp_ix]-qbw, mids[grp_ix]+qbw),
                                rep(ci95[grp_ix, 1],2), type="l", col="blue");
                        points(
                                rep(mids[grp_ix],2),
                                c(ci95[grp_ix, 1], ci95[grp_ix,2]), type="l", col="blue");
                }
        }

        # label sample size
        for(grp_ix in 1:num_grps){
                text(mids[grp_ix], 3*-par()$cxy[2]/2, paste("mean =", round(grp_means[grp_ix], 2)), 
			cex=.95, xpd=T, font=3, adj=c(.5,-1));

                text(mids[grp_ix], 4*-par()$cxy[2]/2, paste("n =",samp_size[grp_ix]), 
			cex=.85, xpd=T, font=3, adj=c(.5,-1));
        }

        connect_significant=function(A, B, ypos, pval){
                abline(h=ypos);
        }

        sigchar=function(x){
                if(x<=.0001){
                        return("***");
                }else if(x<=.001){
                        return("**");
                }else if(x<=.01){
                        return("*");
                }else{
                        return("");
                }
        }

        row_ix=1;
        for(i in 1:(num_grps-1)){

                pvalrow=pval_mat[i,];
                #print(pvalrow);

                signf_pairs=(pvalrow<alpha);
                if(any(signf_pairs)){
                        signf_grps=which(signf_pairs);
                        cat("Pairs: ", i, " to:\n");
                        print(signf_grps);

                        y_offset=plotdatamax+horiz_spacing*row_ix;

                        # Draw line between left/reference to each paired signf grp
                        points(c(
                                mids[i], mids[max(signf_grps)]),
                                rep(y_offset,2),
                                type="l", lend="square"
                        );

                        # Mark left/ref group
                        points(
                                rep(mids[i],2),
                                c(y_offset,y_offset-horiz_spacing/4),
                                type="l", lwd=3, lend="butt");

                        # Mark each signf paired reference group
                        for(pair_ix in signf_grps){
                                points(
                                        rep(mids[pair_ix],2),
                                        c(y_offset,y_offset-horiz_spacing/4),
                                        type="l", lwd=1, lend="butt");


                                # label pvalue
                                paird_pval=sprintf("%5.4f", pvalrow[pair_ix]);
                                text(mids[pair_ix], y_offset, paird_pval,
                                        adj=c(.5, -1), cex=.7);
                                text(mids[pair_ix], y_offset, sigchar(pvalrow[pair_ix]),
                                        adj=c(.5, -1.25), cex=1);
                        }

                        row_ix=row_ix+1;

                }

        }
}

###############################################################################

plot_stats_mat=function(sm, subject_grouping_rec){

	cat("Plotting Stats Matrix...\n");
	num_groups=subject_grouping_rec[["NumGroups"]];
	num_stats=ncol(sm);
	num_indiv=nrow(sm);
	stat_name=colnames(sm);

	cat("Num Groups: ", num_groups, "\n");
	cat("Num Indiv: ", num_indiv, "\n");
	cat("Num Stats: ", num_stats, "\n");

	par(mfrow=c(2,1));

	for(stat_ix in 1:num_stats){
		cat("---------------------------------------------------------\n");
		cat("Plotting: ", stat_name[stat_ix], "\n");
		plot_barplot_wsignf_annot(
			title=stat_name[stat_ix], 
			stat=sm[,stat_ix, drop=F],
			grps=subject_grouping_rec[["GrpToSbj"]]);
	}

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

stats_mat=calculate_stats_on_series(offset_rec, dist_mat);

plot_connected_figure(mds_coord, offset_rec, subject_grouping_rec, 
	subjects_per_plot=as.integer(num_subjects/5), col_assign, 
	ind_colors, title=paste("Metric MDS (", DistanceType,")", sep=""));
plot_connected_figure(mds2_coord, offset_rec, subject_grouping_rec, 
	subjects_per_plot=as.integer(num_subjects/5), col_assign, 
	ind_colors, title=paste("IsoMetric MDS (", DistanceType, ")", sep=""));

plot_sample_distances(dist_mat, offset_rec, subject_grouping_rec, col_assign, ind_colors, 
	dist_type=DistanceType);

plot_sample_dist_by_group(dist_mat, offset_rec, subject_grouping_rec, col_assign, ind_colors, 
	dist_type=DistanceType);

plot_sample_dist_by_group_loess(dist_mat, offset_rec, subject_grouping_rec, col_assign, ind_colors, 
	dist_type=DistanceType);


# Extract individual membership

print(stats_mat);
num_stats=ncol(stats_mat);
cols_per_page=4;
num_pages=ceiling(num_stats/cols_per_page);

cat("Printing Stats to PDF...\n");
cat("Num Stats: ", num_stats, "\n");
cat("Num Pages: ", num_pages, "\n");
cat("Num Stats/Page: ", cols_per_page, "\n");
for(page_ix in 1:num_pages){
	start_col=((page_ix-1)*cols_per_page)+1;
	end_col=min(start_col+cols_per_page-1, num_stats);
	cat("Printing columns: ", start_col, " to ", end_col, "\n", sep="");	
	plot_text(capture.output(print(stats_mat[, start_col:end_col, drop=F])));
}

plot_stats_mat(stats_mat, subject_grouping_rec);

stat_description=c(
	"DESCRIPTION OF STATISTICS:",
	"",
	"",
	"last_time: Last recorded time", 
	"num_time_pts: Number of time points",
	"",
	"Changes, relative to 1st sample:",
	"  average_dist: Average distance samples spent away from 1st sample", 
	"  average_speed: (Total changes in distance)/(Last recorded time)",
	"  total_dist_travelled: Sum of distance travelled",
	"",
	"mean_reversion variables:  Fit linear model across all data points",
	"  mean_reversion_first_dist: expected distance of first sample (y-intercept)", 
	"  mean_reversion_last_dist: expected distance of last sample", 
	"  mean_reversion_stdev_residuals: standard deviation of residuals", 
	"  mean_reversion_slope: slope of linear model",
	"",
	"closest_travel_dist: Distance sample came closest to 1st sample", 
	"closest_travel_time: Time when sample came closest to 1st sample",
	"",
	"furthest_travel_dist: Distance of sample furthest from 1st sample", 
	"furthest_travel_time: Time when sample went furthest from 1st sample",
	"",
	"closest_return_dist: Closest distance sample came to 1st sample after rebounding", 
	"closest_return_time: Time when sample came closest to 1st ample after rebounding",
	"",
	"first_return_dist: Distance when sample first rebounds",
	"first_return_time: Time when sample first rebounds");

par(mfrow=c(1,1));
plot_text(stat_description);

##############################################################################

export_distances_from_start(paste(OutputFileRoot,".dist_from_start.tsv", sep=""),
	offset_rec, dist_mat);

##############################################################################

cat("Done.\n")
dev.off();
warn=warnings();
if(length(warn)){
	print(warn);
}
q(status=0)
