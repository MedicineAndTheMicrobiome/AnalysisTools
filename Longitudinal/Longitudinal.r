
extract_offset=function(factor_mat, sbj_cname, timeoffset_cname, start=-Inf, end=Inf){
	
	var_names=colnames(factor_mat);
	if(!any(sbj_cname==var_names)){
		cat("Error: Could not find ", sbj_cname, " in factors.\n");
		quit(status=-1);
	}
	if(!any(timeoffset_cname==var_names)){
		cat("Error: Could not find ", timeoffset_cname, " in factors.\n");
		quit(status=-1);
	}

	offsets_mat=as.data.frame(factor_mat[,c(sbj_cname, timeoffset_cname)]);
	colnames(offsets_mat)=c("SubjectID", "Offsets");

        # Change number IDs to strings
        if(is.numeric(offsets_mat[,"SubjectID"])){
                numdigits=log10(max(offsets_mat[,"SubjectID"]))+1;
                prtf_str=paste("%0",numdigits,"d", sep="");
                offsets_mat[,"SubjectID"]=paste("#",
                         sprintf(prtf_str, offsets_mat[,"SubjectID"]), sep="");
        }else{
		offsets_mat[,"SubjectID"]=as.character(offsets_mat[,"SubjectID"]);
	}

	# Filter offsets by start and end inclusive
	keep_ix=(offsets_mat[, "Offsets"]>=start & offsets_mat[, "Offsets"]<=end)	
	offsets_mat=offsets_mat[keep_ix,];

	cat("Number of Rows in Offset Matrix: ", nrow(offsets_mat), "\n");

        offsets_data=list();
        offsets_data[["matrix"]]=offsets_mat;
        offsets_data[["SampleIDs"]]=sort(rownames(offsets_mat));
        offsets_data[["SubjectIDs"]]=sort(unique(offsets_mat[, "SubjectID"]));
	offsets_data[["NumSubjects"]]=length(offsets_data[["SubjectIDs"]]);
        offsets_data[["Offsets"]]=sort(unique(offsets_mat[, "Offsets"]));
	
	# Offsets by subject
	offsets_by_sbj=list()
	for(sbj in offsets_data[["SubjectIDs"]]){
		sbj_ix=offsets_mat[,"SubjectID"]==sbj;
		offsets_by_sbj[[sbj]]=offsets_mat[sbj_ix,,drop=F];
	}
	offsets_data[["OffsetsBySubject"]]=offsets_by_sbj;

	# Store range information
	if(start==-Inf){
		stag="start";
	}else{
		stag=as.character(start);
	}

	if(end==Inf){
		etag="end";
	}else{
		etag=as.character(end);
	}

	range_tag=paste(stag, "_to_", etag, sep="");
	range_tag=gsub("-", "n", range_tag);

	offsets_data[["RangeTag"]]=range_tag;
	offsets_data[["Start"]]=start;
	offsets_data[["End"]]=end;
	offsets_data[["Earliest"]]=min(offsets_data[["Offsets"]]);
	offsets_data[["Latest"]]=max(offsets_data[["Offsets"]]);
	offsets_data[["Range"]]=offsets_data[["Latest"]]-offsets_data[["Earliest"]];
	offsets_data[["MinOffsetSep"]]=min(diff(offsets_data[["Offsets"]]));
	offsets_data[["NumUniqOffsets"]]=length(offsets_data[["Offsets"]]);
	offsets_data[["OffsetsAsChar"]]=as.character(offsets_data[["Offsets"]]);

	widths=max(ceiling(log10(offsets_data[["Latest"]])+1));
	format=paste("%0",widths+1,".1f", sep="");
	offsets_data[["OffsetsAsCharZeroPad"]]=sprintf(format,offsets_data[["Offsets"]]);

	print(offsets_data);
        return(offsets_data);

}

get_colors=function(num_col, alpha=1){
        colors=hsv(seq(0,1,length.out=num_col+1), c(1,.5), c(1,.75,.5), alpha=alpha);
        color_mat_dim=ceiling(sqrt(num_col));
        color_pad=rep("grey", color_mat_dim^2);
        color_pad[1:num_col]=colors[1:num_col];
        color_mat=matrix(color_pad, nrow=color_mat_dim, ncol=color_mat_dim);
        colors=as.vector(t(color_mat));
        colors=colors[colors!="grey"];
}

create_GrpToSbj_map=function(subjects_arr, groups_arr){
        if(length(subjects_arr) != length(groups_arr)){
                cat("Error: Subj/Grp array lengths do no match.\n");
                cat("Subjects:\n");
                print(subjects_arr);
                cat("\nGroups:\n");
                print(groups_arr);
                quit(status=-1);
        }

	# Create subject lookup by group
        uniq_grps=sort(unique(groups_arr));
        grp_to_sbj_map=list();
        for(grp_ix in uniq_grps){
                grp_to_sbj_map[[grp_ix]]=unique(sort(as.character(subjects_arr[groups_arr==grp_ix])));
        }

	# Create map from subject to group
	sbj_to_grp_map=groups_arr;
	names(sbj_to_grp_map)=subjects_arr;

	# Assign colors to groups and subjects
	num_grps=length(uniq_grps);
	uniq_sbjs=sort(unique(subjects_arr));
	num_sbjs=length(uniq_sbjs);

	# Opaque
	group_colors=get_colors(num_grps,1);
	names(group_colors)=uniq_grps;
	subject_colors=get_colors(num_sbjs,1);
	names(subject_colors)=uniq_sbjs;

	# Transparent
	group_colors_transparent=get_colors(num_grps,.2);
	names(group_colors_transparent)=uniq_grps;
	subject_colors_transparent=get_colors(num_sbjs,.2);
	names(subject_colors_transparent)=uniq_sbjs;

	rec=list();
	rec[["Groups"]]=uniq_grps;
	rec[["GroupColors"]]=group_colors;
	rec[["GroupColorsTransparent"]]=group_colors_transparent;
	rec[["Subjects"]]=uniq_sbjs;
	rec[["SubjectColors"]]=subject_colors;
	rec[["SubjectColorsTransparent"]]=subject_colors_transparent;
	rec[["NumSubjects"]]=num_sbjs;
	rec[["NumGroups"]]=num_grps;
	rec[["GrpToSbj"]]=grp_to_sbj_map;
	rec[["SbjToGrp"]]=sbj_to_grp_map;	

        return(rec);
}


group_offsets=function(offsets_data){
	
	cat("Reorganizing Raw Offset Data...\n");
	mat=offsets_data[["matrix"]];

	indiv_offsets=list();

	indiv_rows=mat[,"Indiv ID"];

	indivs=sort(unique(indiv_rows));
	
	# Groups offsets by individual
	for(indiv_ix in indivs){

		# Extract rows for individual
		rows_ix=(indiv_rows==indiv_ix)
		cur_samp=mat[rows_ix,,drop=F];
	
		# Reorder by offset
		reorder_ix=order(cur_samp[,"Offsets"]);
		cur_samp=cur_samp[reorder_ix,,drop=F];

		# Add column for sample ID from rowname
		samp_ids=rownames(cur_samp);
		cur_samp=cbind(cur_samp, samp_ids);
		colnames(cur_samp)[4]="Samp ID";

		# Store record in list
		indiv_offsets[[indiv_ix]]=cur_samp;
	}

	# Group individuals by Group ID
	grp_info=list();
	grps=as.character(sort(unique(mat[,"Group ID"])));
	for(gr in grps){
		grp_info[[gr]]=as.character(sort(unique(mat[(mat[,"Group ID"]==gr),"Indiv ID"])));
	}
	

	res=list();
	res[["OffsetsByIndiv"]]=indiv_offsets;
	res[["IndivByGrp"]]=grp_info;
	res[["Individuals"]]=indivs;
	res[["Offsets"]]=sort(unique(mat[,"Offsets"]));
	res[["Groups"]]=grps;

	cat("ok.\n");
	return(res);
}

###############################################################################

calculate_stats_on_series_distance=function(offset_rec, dist_mat){

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
        num_subjects=offset_rec[["NumSubjects"]];

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

        dist_mat=as.matrix(dist_mat);

        tmp_mat=matrix(NA, nrow=num_subjects, ncol=length(stat_names));
        rownames(tmp_mat)=uniq_indiv_ids;
        colnames(tmp_mat)=stat_names;

        for(cur_id in uniq_indiv_ids){

                cur_offsets=offset_rec[["OffsetsBySubject"]][[cur_id]];
                num_timepts=nrow(cur_offsets);

                tmp_mat[cur_id, "last_time"]=cur_offsets[num_timepts, "Offsets"];
                tmp_mat[cur_id, "num_time_pts"]=num_timepts;

                samp_ids=rownames(cur_offsets);

                if(num_timepts>1){
                        cur_dist=dist_mat[samp_ids[1], samp_ids];
                        cur_times=cur_offsets[,"Offsets"];

                        tmp_mat[cur_id, "average_dist"]=avg_dist(cur_dist, cur_times);
                        tmp_mat[cur_id, "average_speed"]=avg_speed(cur_dist, cur_times);
                        tmp_mat[cur_id, "total_dist_travelled"]=tot_dist_travelled(cur_dist, cur_times);

                        res=mean_reversion(cur_dist, cur_times);
                        tmp_mat[cur_id, "mean_reversion_first_dist"]=res[["first_dist"]];
                        tmp_mat[cur_id, "mean_reversion_last_dist"]=res[["last_dist"]];
                        tmp_mat[cur_id, "mean_reversion_stdev_residuals"]=res[["sd_res"]];
                        tmp_mat[cur_id, "mean_reversion_slope"]=res[["slope"]];

                        res=closest_travel(cur_dist, cur_times);
                        tmp_mat[cur_id, "closest_travel_dist"]=res[["dist"]];
                        tmp_mat[cur_id, "closest_travel_time"]=res[["time"]];

                        res=furthest_travel(cur_dist, cur_times);
                        tmp_mat[cur_id, "furthest_travel_dist"]=res[["dist"]];
                        tmp_mat[cur_id, "furthest_travel_time"]=res[["time"]];

                        res=closest_return(cur_dist, cur_times);
                        tmp_mat[cur_id, "closest_return_dist"]=res[["dist"]];
                        tmp_mat[cur_id, "closest_return_time"]=res[["time"]];

                        res=first_return(cur_dist, cur_times);
                        tmp_mat[cur_id, "first_return_dist"]=res[["dist"]];
                        tmp_mat[cur_id, "first_return_time"]=res[["time"]];
                }
        }

	results=list();
	matrix_template=matrix(NA, nrow=num_subjects, ncol=1);
	colnames(matrix_template)="distance";
	rownames(matrix_template)=uniq_indiv_ids;

	for(stat_ix in stat_names){
		results[[stat_ix]]=matrix_template;
		for(sbx in uniq_indiv_ids){
			results[[stat_ix]][sbx, "distance"]=tmp_mat[sbx, stat_ix];
		}	

	}

        return(results);
}

###############################################################################

calc_longitudinal_stats=function(offset_rec, alr_cat_val){

        cat("Calculating Longitudinal Statistics...\n");

        l_min=function(x, y){
                res=min(y);
                return(res);
        }
        l_max=function(x, y){
                res=max(y);
                return(res);
        }
        l_median=function(x, y){
                res=median(y);
                return(res);
        }
        l_mean=function(x, y){
                res=mean(y);
                return(res);
        }
        l_stdev=function(x, y){
                res=sd(y);
                return(res);
        }
        l_range=function(x, y){
                r=range(y);
                span=abs(r[1]-r[2]);
                return(span);
        }
        l_N=function(x, y){
                return(length(x));
        }
        l_last_time=function(x, y){
                return(max(x));
        }

        l_volatility=function(x, y){
                if(length(x)>1){
                        fit=lm(y~x);
                        sumfit=summary(fit);
                        vol=sd(sumfit$residuals);
                        return(vol);
                }else{
                        return(NA);
                }
        }
        l_slope=function(x, y){
                if(length(x)>1){
                        fit=lm(y~x);
                        slope=fit$coefficients["x"];
                        return(slope);
                }else{
                        return(NA);
                }
        }

        l_time_wght_avg=function(x, y){
                npts=length(x);
                if(npts>1){

                        tot_avg=0;

                        for(i in 1:(npts-1)){
                                avg_val=(y[i]+y[i+1])/2;
                                duration=(x[i+1]-x[i]);
                                tot_avg=tot_avg+(avg_val*duration);
                        }

                        norm=tot_avg/(x[npts]-x[1]);
                        return(norm);
                }else{
                        return(NA);
                }
        }

        l_time_at_max=function(x, y){
                max_val=max(y);
                ix=min(which(y==max_val));
                return(x[ix]);
        }

        l_time_at_min=function(x, y){
                min_val=min(y);
                ix=min(which(y==min_val));
                return(x[ix]);
        }

        l_time_closest_to_t0=function(x, y){
                starty=y[1];
                y=y[-1];
                dist=abs(y-starty);
                min_dist=min(dist);
                min_ix=min(which(min_dist==dist));
                return(x[min_ix+1]);
        }

        l_time_furthest_fr_t0=function(x, y){
                starty=y[1];
                y=y[-1];
                dist=abs(y-starty);
                max_dist=max(dist);
                max_ix=min(which(max_dist==dist));
                return(x[max_ix+1]);
        }

        l_start_end_diff=function(x,y){
                start=y[1];
                end=tail(y,1);
                return(end-start);
        }

        l_convexcave=function(x, y){
                num_pts=length(x);
                # y=mx+b
                # b=y-mx

                m=(y[num_pts]-y[1])/(x[num_pts]-x[1]);
                b=y[1]-m*x[1];

                lvl_y=numeric(num_pts);
                for(i in 1:num_pts){
                        lvl_y[i]=y[i]-(m*x[i]+b);
                }

                cum_sum=0;
                for(i in 1:(num_pts-1)){
                        dx=x[i+1]-x[i];
                        avgy=(lvl_y[i+1]+lvl_y[i])/2;
                        avg=dx*avgy/2;
                        cum_sum=cum_sum+avg;
                }
                vexcav=cum_sum/(x[num_pts]-x[1]);

                return(vexcav);
        }


        # statistic:
        #    ALR:
        #       individual:

        stat_name=c(
                "min", "max", "range",
                "volatility", "slope", "time_wght_avg",
                "time_at_max", "time_at_min",
                "time_closest_to_t0", "time_furthest_fr_t0",
                "start_end_diff",
                "convexcave"
        );

        results=list();

        alrcat=colnames(alr_cat_val);
        individuals=as.character(offset_rec[["SubjectIDs"]]);

        cat("\n");
        cat("Individuals:\n");
        print(individuals);

        cat("\n");
        cat("Categories:\n");
        print(alrcat);

        num_cat=ncol(alr_cat_val);
        num_ind=length(individuals);


        for(stat_ix in stat_name){

                results[[stat_ix]]=list();

                tmp_mat=matrix(NA, nrow=num_ind, ncol=num_cat);
                rownames(tmp_mat)=individuals;
                colnames(tmp_mat)=alrcat;

                for(cat_ix in alrcat){

                        for(ind_ix in individuals){

                                indv_offsets=offset_rec[["OffsetsBySubject"]][[ind_ix]];
                                samp_ids=rownames(indv_offsets);

                                time=indv_offsets[,"Offsets"];
                                val=alr_cat_val[samp_ids, cat_ix];

                                funct_name=paste("l_", stat_ix, sep="");

                                call_res=do.call(funct_name, list(x=time, y=val));

                                tmp_mat[ind_ix, cat_ix]=call_res;

                        }

                }

                results[[stat_ix]]=tmp_mat;

        }

        return(results);
}

###############################################################################

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, deci_pts=4,
        label_zeros=T, counts=F, value.cex=1,
        plot_col_dendr=F,
        plot_row_dendr=F
){

        cat("Working on: ", title, "\n");

        num_row=nrow(mat);
        num_col=ncol(mat);

        if(num_row==0 || num_col==0){
                cat("Nothing to plot.\n");
                return();
        }

        any_nas=any(is.na(mat));

        if(num_row==1 || any_nas){
                plot_row_dendr=F;
        }
        if(num_col==1 || any_nas){
                plot_col_dendr=F;
        }

        row_names=rownames(mat);
        col_names=colnames(mat);

        orig.par=par(no.readonly=T);

        cat("Num Rows: ", num_row, "\n");
        cat("Num Cols: ", num_col, "\n");

        # Flips the rows, so becuase origin is bottom left
        mat=mat[rev(1:num_row),, drop=F];

        # Generate a column scheme
        num_colors=50;
        color_arr=rainbow(num_colors, start=0, end=4/6);
        if(high_is_hot){
                color_arr=rev(color_arr);
        }

        # Provide a means to map values to an (color) index
        remap=function(in_val, in_range, out_range){
                in_prop=(in_val-in_range[1])/(in_range[2]-in_range[1])
                out_val=in_prop*(out_range[2]-out_range[1])+out_range[1];
                return(out_val);
        }

        # If range is not specified, find it based on the data
        if(is.na(plot_min)){
                plot_min=min(mat, na.rm=T);
        }
        if(is.na(plot_max)){
                plot_max=max(mat, na.rm=T);
        }

        if(plot_min>=-1 && plot_max<=1){
                fractions_only=T;
        }else{
                fractions_only=F;
        }
        cat("Plot min/max: ", plot_min, "/", plot_max, "\n");

        # Get Label lengths
        row_max_nchar=max(nchar(row_names));
        col_max_nchar=max(nchar(col_names));
        cat("Max Row Names Length: ", row_max_nchar, "\n");
        cat("Max Col Names Length: ", col_max_nchar, "\n");

        ##################################################################################################

        get_dendrogram=function(in_mat, type){

                if(type=="row"){
                        dendist=dist(in_mat);
                }else{
                        dendist=dist(t(in_mat));
                }

                get_clstrd_leaf_names=function(den){
                # Get a list of the leaf names, from left to right
                        den_info=attributes(den);
                        if(!is.null(den_info$leaf) && den_info$leaf==T){
                                return(den_info$label);
                        }else{
                                lf_names=character();
                                for(i in 1:2){
                                        lf_names=c(lf_names, get_clstrd_leaf_names(den[[i]]));
                                }
                                return(lf_names);
                        }
                }

                hcl=hclust(dendist, method="ward.D2");
                dend=list();
                dend[["tree"]]=as.dendrogram(hcl);
                dend[["names"]]=get_clstrd_leaf_names(dend[["tree"]]);
                return(dend);
        }

        ##################################################################################################
        # Comput Layouts
        col_dend_height=ceiling(num_row*.1);
        row_dend_width=ceiling(num_col*.2);

        heatmap_height=num_row;
        heatmap_width=num_col;

        if(plot_col_dendr && plot_row_dendr){
                layoutmat=matrix(
                        c(
                        rep(c(rep(4, row_dend_width), rep(3, heatmap_width)), col_dend_height),
                        rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height)
                        ), byrow=T, ncol=row_dend_width+heatmap_width);

                col_dendr=get_dendrogram(mat, type="col");
                row_dendr=get_dendrogram(mat, type="row");

                mat=mat[row_dendr[["names"]], col_dendr[["names"]], drop=F];

        }else if(plot_col_dendr){
                layoutmat=matrix(
                        c(
                        rep(rep(2, heatmap_width), col_dend_height),
                        rep(rep(1, heatmap_width), heatmap_height)
                        ), byrow=T, ncol=heatmap_width);

                col_dendr=get_dendrogram(mat, type="col");
                mat=mat[, col_dendr[["names"]], drop=F];

        }else if(plot_row_dendr){
                layoutmat=matrix(
                        rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height),
                        byrow=T, ncol=row_dend_width+heatmap_width);

                row_dendr=get_dendrogram(mat, type="row");
                mat=mat[row_dendr[["names"]],, drop=F];
        }else{
                layoutmat=matrix(
                        rep(1, heatmap_height*heatmap_width),
                        byrow=T, ncol=heatmap_width);
        }

        #print(layoutmat);
        layout(layoutmat);

        ##################################################################################################

        par(oma=c(col_max_nchar*.60, 0, 3, row_max_nchar*.60));
        par(mar=c(0,0,0,0));
        plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), xaxt="n", yaxt="n", bty="n", xlab="", ylab="");
        mtext(title, side=3, line=0, outer=T, font=2);

        # x-axis
        axis(side=1, at=seq(.5, num_col-.5, 1), labels=colnames(mat), las=2, line=-1.75, cex.axis=value.cex);
        axis(side=4, at=seq(.5, num_row-.5, 1), labels=rownames(mat), las=2, line=-1.75, cex.axis=value.cex);

        if(log_col){
                plot_min=log10(plot_min+.0125);
                plot_max=log10(plot_max+.0125);
        }

        for(x in 1:num_col){
                for(y in 1:num_row){

                        if(log_col){
                                col_val=log10(mat[y,x]+.0125);
                        }else{
                                col_val=mat[y,x];
                        }

                        remap_val=remap(col_val, c(plot_min, plot_max), c(1, num_colors));
                        col_ix=ceiling(remap_val);

                        rect(x-1, y-1, (x-1)+1, (y-1)+1, border=NA, col=color_arr[col_ix]);

                        if(is.na(mat[y,x]) || mat[y,x]!=0 || label_zeros){
                                if(counts){
                                        text_lab=sprintf("%i", mat[y,x]);
                                }else{
                                        text_lab=sprintf(paste("%0.", deci_pts, "f", sep=""), mat[y,x]);
                                        if(fractions_only){
                                                if(!is.na(mat[y,x])){
                                                        if(mat[y,x]==-1 || mat[y,x]==1){
                                                                text_lab=as.integer(mat[y,x]);
                                                        }else{
                                                                text_lab=gsub("0\\.","\\.", text_lab);
                                                        }
                                                }
                                        }
                                }
                                text(x-.5, y-.5, text_lab, srt=atan(num_col/num_row)/pi*180, cex=value.cex, font=2);
                        }
                }
        }

        ##################################################################################################

        par(mar=c(0, 0, 0, 0));

        if(plot_row_dendr && plot_col_dendr){
                rdh=attributes(row_dendr[["tree"]])$height;
                cdh=attributes(col_dendr[["tree"]])$height;
                plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", bty="n", xlim=c(rdh, 0));
                plot(col_dendr[["tree"]], leaflab="none",xaxt="n", yaxt="n", bty="n", ylim=c(0, cdh));
                plot(0,0, type="n", bty="n", xaxt="n", yaxt="n");
                #text(0,0, "Placeholder");
        }else if(plot_row_dendr){
                rdh=attributes(row_dendr[["tree"]])$height;
                plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", bty="n", xlim=c(rdh, 0));
                #text(0,0, "Row Dendrogram");
        }else if(plot_col_dendr){
                cdh=attributes(col_dendr[["tree"]])$height;
                plot(col_dendr[["tree"]], leaflab="none", xaxt="n", yaxt="n", bty="n", ylim=c(0, cdh));
                #text(0,0, "Col Dendrogram");
        }

        par(orig.par);
}

###############################################################################


plot_barplot_wsignf_annot=function(title, stat, grps, alpha=0.1, samp_gly=T){
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

				grpAval=stat[grps[[grpAnm]]];
				grpBval=stat[grps[[grpBnm]]];

				grpAmean=mean(grpAval);
				grpBmean=mean(grpBval);

				if(!(length(grpAval)==0 || length(grpBval)==0)){
					res=wilcox.test(grpAval, grpBval);
					if(is.nan(res$p.value)){
						pval=1;
					}else{
						pval=res$p.value;
					}
					pval_mat[grpAnm, grpBnm]=pval;
					if(pval<=alpha){
						signf=rbind(signf, 
							c(grpAnm, sprintf("%8.4f",grpAmean), 
							grpBnm, sprintf("%8.4f",grpBmean), 
							sprintf("%8.4f", (grpBmean-grpAmean)),
							sprintf("%5.3f",pval)));
					}
				}
                        }
                }
        }

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
        annot_line_prop=2/5; # proportion of pl
        min_95ci=min(c(ci95[,1], stat), na.rm=T);
        max_95ci=max(c(ci95[,2], stat), na.rm=T);
        minmax_span=max_95ci-min_95ci;
	if(minmax_span==0){minmax_span=min_95ci/5;}

        space_for_annotations=minmax_span*0.4;
        plotdatamax=max_95ci+space_for_annotations;
        plotdatamin=min_95ci-minmax_span*0.4;

	padding_for_mean_n_label=3;
	# Amount of y-space to use per significant different annotation
        horiz_spacing=space_for_annotations/(max(1,num_signf_rows)+padding_for_mean_n_label);

	cat("space_for_annotations:", space_for_annotations, "\n");
	cat("horizon_spacing:", horiz_spacing, "\n");

        # Start plot
        par(mar=c(8,5,4,3));
        cat("  Plot Limits: (", plotdatamin, ", ", plotdatamax, ")\n");
        plot(0, type="n",
                ylim=c(plotdatamin, plotdatamax),
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
        qbw=bar_width/6;

        # Label group names x-axis
        text(mids-par()$cxy[1]/2, rep(plotdatamin, num_grps),
                group_names, srt=-45, xpd=T, pos=4,
		font=2,
                cex=min(c(1.5,.7*bar_width/par()$cxy[1])));

        # Scatter
        if(samp_gly){
                for(grp_ix in 1:num_grps){
                        grpnm=group_names[grp_ix];
                        pts=stat[grps[[grpnm]]];
                        numpts=length(pts);
                        points(
                                #rep(mids[grp_ix], numpts),
                                mids[grp_ix]+rnorm(numpts, 0, bar_width/10),
                                pts, col="darkblue", cex=1, type="p");
                }
        }

        # label CI's
        for(grp_ix in 1:num_grps){
                if(samp_size[grp_ix]>=40){
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
                #text(mids[grp_ix], 3*-par()$cxy[2]/2, paste("mean =", round(grp_means[grp_ix], 2)),
                text(mids[grp_ix], max_95ci, paste("mean =", round(grp_means[grp_ix], 2)),
                        cex=.95, xpd=T, font=3, adj=c(.5,-2));

                #text(mids[grp_ix], 4*-par()$cxy[2]/2, paste("n =",samp_size[grp_ix]),
                text(mids[grp_ix], max_95ci, paste("n =",samp_size[grp_ix]),
                        cex=.95, xpd=T, font=2, adj=c(.5,-3));
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

                        y_offset=max_95ci+(horiz_spacing*(row_ix+padding_for_mean_n_label));

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
	
	return(signf);

}

###############################################################################

plot_pairwise_grp_comparisons=function(longit_stats, grp_to_sbj_info_rec, plots_pp=4){

        stat_colnames=c(
                "Statistic",
                "Category",
                "Gr1",
                "mean(Gr1)",
                "Gr2",
                "mean(Gr2)",
                "diff(Gr2-Gr2)",
                "p-value"
        );

        stat_table=matrix(NA, nrow=0, ncol=length(stat_colnames));
        colnames(stat_table)=stat_colnames;

        unique_group_names=as.character(grp_to_sbj_info_rec[["Groups"]]);
        num_grps=grp_to_sbj_info_rec[["NumGroups"]];
        group_to_subject_map=grp_to_sbj_info_rec[["GrpToSbj"]];

        stat_names=names(longit_stats);
        for(stat_ix in stat_names){

		cat("Working on stat: ", stat_ix, "\n");
                num_target_cat=ncol(longit_stats[[stat_ix]]);
                cat_names=colnames(longit_stats[[stat_ix]]);

                grp_mat=matrix(NA, nrow=num_grps, ncol=num_target_cat);
                rownames(grp_mat)=unique_group_names;
                colnames(grp_mat)=cat_names;

                # Compute means by group
                for(grp_ix in unique_group_names){
                        grp_members=group_to_subject_map[[grp_ix]];
                        grp_mat[grp_ix,]=apply(
                                longit_stats[[stat_ix]][grp_members,,drop=F], 2,
                                function(x){mean(x, na.rm=T)});
                }

                # Plot heatmap
                paint_matrix(grp_mat, paste("mean(", stat_ix, ") for Grouping: ", GroupCol, sep=""));

                par(mfrow=c(plots_pp,1));
                par(oma=c(0,0,2,0));
                label_oma=F
                plot_ix=0;
                for(cat_ix in 1:num_target_cat){
		
			cat("Working on category: ", cat_names[cat_ix], "\n");

                        signf=plot_barplot_wsignf_annot(
                                title=cat_names[cat_ix],
                                longit_stats[[stat_ix]][,cat_ix,drop=F],
                                group_to_subject_map,
                                alpha=0.1
                        );

                        if(length(signf)){
                                num_sig_rows=nrow(signf);
                                newrows=cbind(
                                        rep(stat_ix, num_sig_rows),
                                        rep(cat_names[cat_ix], num_sig_rows),
                                        signf);
                                stat_table=rbind(stat_table, newrows);
                        }

                        plot_ix=plot_ix+1;
                        if(plot_ix==plots_pp){
                                mtext(stat_ix, outer=T, cex=1.5, col="blue", font=2);
                                plot_ix=0;
                        }
                }
                if(plot_ix!=plots_pp){
                        mtext(stat_ix, outer=T, cex=1.5, col="blue", font=2);
                }



        }

        return(stat_table);
}

###############################################################################

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

output_stat_table_alternate_ordering=function(stat_table, output_root){

	if(nrow(stat_table)>0){

		pvals=as.numeric(stat_table[,"p-value"]);
		signf=sapply(pvals, sigchar);
		stat_table=cbind(stat_table, signf);
		row_idx_str=paste(1:nrow(stat_table), ".", sep="");

		options(width=1000);

		#------------------------------------------------
		cat("Ordering by P-value...\n");
		# First order by pvalue, then do stable sort
		stat_order_ix=order(pvals);
		stat_table=stat_table[stat_order_ix,];
		print(stat_table);

		#------------------------------------------------
		cat("Ordering by Stat Name...\n");

		stat_order_ix=order(stat_table[,"Statistic"], method="shell");
		out_stat_table=stat_table[stat_order_ix,];
		rownames(out_stat_table)=row_idx_str;
		out=capture.output(print(out_stat_table, quote=F));

		cat("Writing by statistic:\n");
		plot_text(c(
			"By Statistic:",
			"",
			out));

		#------------------------------------------------
		cat("Ordering by Category...\n");

		stat_order_ix=order(stat_table[,"Category"], method="shell");
		out_stat_table=stat_table[stat_order_ix,];
		rownames(out_stat_table)=row_idx_str;
		out=capture.output(print(out_stat_table, quote=F));

		cat("Writing by category:\n");
		plot_text(c(
			"By Category:",
			"",
			out));

		#------------------------------------------------
		cat("Outputing by P-Value...\n");

		out_stat_table=stat_table;
		rownames(out_stat_table)=row_idx_str;
		out=capture.output(print(out_stat_table, quote=F));
		cat("Writing by p-value:\n");
		plot_text(c(
			"By P-value:",
			"",
			out));

		#------------------------------------------------

		cat("Writing stat table...\n");
		write.table(out_stat_table, 
			file=paste(output_root,".lngt_stats.comp.tsv", sep=""), quote=F, sep="\t");

	}else{
		plot_text(c(
			"No Significant Differences Identified Between Groups."
		));
		fh=file(paste(output_root,".lngt_stats.comp.tsv", sep=""), "w");
		cat(file=fh, "Sorry.  No significant differences identified between groups.\n");
		close(fh);
	}

}

###############################################################################

longit_stat_descriptions=c(
        "Explanation of Longitudinal Statistics:",
        "",
	"",
        "These stats are calculated independent of the time component:",
        "  min: minimum value",
        "  max: maximum value",
        "  median: median value",
        "  mean: average value",
        "  stdev: standard deviation of values",
        "  range: (max-min) of values",
        "  N: number of samples",
        "  last_time: last recorded sample",
        "",
        "These stats are based on fitting a linear model (value = intercept + time + error):",
        "  volatility:  standard deviation of residuals",
        "  slope: slope of line",
        "",
        "These stats take time into account:",
        "  time_weighted_average: average between two time points weighted by time between points",
        "  time_at_max: time at value's maximum",
        "  time_at_min: time at value's minimum",
        "  time_closest_to_t0: time when value is closest to first time value",
        "  time_furthest_fr_t0: time when value is furthest from first time value",
        "  start_end_diff: difference between first and last time points (end-start)",
        "  convexcave: if curve is convex or conave, stat will be <0 or >0, respectively",
        "       A concave curve looks like a hill, a convex curve looks like a pit"
);

longit_stat_description_distance=c(
        "Explanation of Longitudinal Statistics:",
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
        "first_return_time: Time when sample first rebounds"
);


###############################################################################
