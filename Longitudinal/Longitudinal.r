
load_offset=function(fname){

        cat("Loading Offsets: ", fname, "\n");
        offsets_mat=read.delim(fname,  header=TRUE, row.names=1, sep="\t", comment.char="#", quote="");

        num_col=ncol(offsets_mat);
        cat("Num Columns Found: ", num_col, "\n");

	if(num_col==2){
		offsets_mat=cbind(offsets_mat, "NoGroup");
		num_col=3;
	}

        extra_colnames=colnames(offsets_mat);
        print(extra_colnames);
        colnames(offsets_mat)=c("Indiv ID", "Offsets", "Group ID", extra_colnames[4:num_col])[1:num_col];

        # Change number IDs to strings
        if(is.numeric(offsets_mat[,"Indiv ID"])){
                numdigits=log10(max(offsets_mat[,"Indiv ID"]))+1;
                prtf_str=paste("%0",numdigits,"d", sep="");
                offsets_mat[,"Indiv ID"]=paste("#",
                         sprintf(prtf_str, offsets_mat[,"Indiv ID"]), sep="");
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

calculate_stats_on_series=function(offset_mat, dist_mat){

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

        uniq_indiv_ids=sort(unique(offset_mat[,"Indiv ID"]));
        num_ind=length(uniq_indiv_ids);

        cat("IDs:\n");
        print(uniq_indiv_ids);
        cat("Num Individuals: ", num_ind, "\n");

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

                row_ix=(offset_mat[,"Indiv ID"]==cur_id);
                cur_offsets=offset_mat[row_ix,,drop=F];

                # Order offsets
                ord=order(cur_offsets[,"Offsets"]);
                cur_offsets=cur_offsets[ord,,drop=F];

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
                                pval_mat[grpAnm, grpBnm]=res$p.value;
                                if(res$p.value<=alpha){
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
        qbw=bar_width/6;

        # Label x-axis
        text(mids-par()$cxy[1]/2, rep(6*-par()$cxy[2]/2, num_grps),
                group_names, srt=-45, xpd=T, pos=4,
                cex=min(c(1,.7*bar_width/par()$cxy[1])));

        # Scatter
        if(samp_gly){
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

