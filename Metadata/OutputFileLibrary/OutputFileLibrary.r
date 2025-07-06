
source("~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r");

##############################################################################

plot_text=function(strings, max_lines_pp=NULL, oma_tag=NULL, echo=F){

        orig.par=par(no.readonly=T);

        par(mfrow=c(1,1));
        par(mar=rep(0, 4));
        par(family="Courier");

        # Allocate space for putting tag on the bottom outer margin
        if(is.null(oma_tag)){
                par(oma=rep(.5,4));
        }else{
                par(oma=c(1, .5, .5, .5));
        }

	# Estimate the size of the output space
        paper_dim=par("pin");
        cat("Estimated Paper Dim: x= ", paper_dim[1], "in, y= ", paper_dim[2], "in\n", sep="");

        def_lines_per_page=ceiling(paper_dim[2]*6);
        if(is.null(max_lines_pp)){
                #default
                lines_per_page=def_lines_per_page;
                text_size=1;
                cat("  12 Pt Lines per page: ", lines_per_page, "\n", sep="");
        }else{
                lines_per_page=ceiling(max_lines_pp);
                text_size=def_lines_per_page/lines_per_page;
                text_size=max(0.01, text_size);
                text_size=min(1, text_size);
                cat("  User-defined Lines per page: ", lines_per_page, "\n", sep="");
        }
        cat("  cex: ", text_size, "\n");

        num_lines=length(strings);
        cat("Num Lines to Plot: ", num_lines, "\n");

        num_pages=max(1, ceiling(num_lines/lines_per_page));
        cat("Num Pages: ", num_pages, " for ", num_lines, " lines.\n");

        lines_pp=min(num_lines, lines_per_page);
        for(p in 1:num_pages){

                plot(0,0, xlim=c(0,1), ylim=c(0,lines_per_page),
                        type="n",  xaxt="n", yaxt="n",
                        xlab="", ylab="", bty="n"
                        );

                if(!is.null(oma_tag)){
                        mtext(paste("[", oma_tag, "]", sep=""), side=1, col="grey25");
                }

                start=(p-1)*lines_pp+1;
                end=start+lines_pp-1;
                end=min(end, num_lines);
                line=1;
                for(i in start:end){
                        strings[i]=gsub("\t", "", strings[i]);
                        text(0, lines_per_page-line, strings[i], pos=4, cex=text_size);
                        line=line+1;

                        if(echo){
                                cat(strings[i],"\n");
                        }

                }

        }

        par(orig.par);
}


plot_title_page=function(title, subtitle=""){

        orig.par=par(no.readonly=T);
        par(family="serif");
        par(mfrow=c(1,1));

        plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n",  xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                );

        # Title
        title_cex=3;
        title_line=1;
        text(0.5, title_line, title, cex=title_cex, font=2, adj=c(.5,1));

        # Subtitle
        num_subt_lines=length(subtitle);
        cxy=par()$cxy;
        for(i in 1:num_subt_lines){
                text(.5, title_line -title_cex*cxy[2] -i*cxy[2], subtitle[i], adj=.5);
        }

        par(orig.par);
}


paint_matrix=function(mat, title="", subtitle="", plot_min=NA, plot_max=NA, log_col=F, 
	high_is_hot=T, deci_pts=4, 
	label_zeros=T, counts=F, value.cex=2,
	plot_col_dendr=F,
	plot_row_dendr=F,
	suppress_grid_lines=F,
	show_leading_zero=T
){

        num_row=nrow(mat);
        num_col=ncol(mat);

	row_names=rownames(mat);
	col_names=colnames(mat);

	orig.par=par(no.readonly=T);

	cat("Painting Matrix: ", title, "\n");
        cat("Num Rows: ", num_row, "\n");
        cat("Num Cols: ", num_col, "\n");


	if(num_row==0 || num_col==0){
		plot(0, type="n", xlim=c(-1,1), ylim=c(-1,1), xaxt="n", yaxt="n", bty="n", xlab="", ylab="",
			main=title);
		text(0,0, "No data to plot...");
		return();
	}

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
        cat("Plot min/max: ", plot_min, "/", plot_max, "\n");

	# Get Label lengths
	row_max_nchar=max(nchar(row_names), 1, na.rm=T);
	col_max_nchar=max(nchar(col_names), 1, na.rm=T);
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

	if(num_row==1){
		plot_row_dendr=F;
	}
	if(num_col==1){
		plot_col_dendr=F;
	}

	# Don't plot dendrogram if there are any NAs in the matrix
	if(any(is.na(mat))){
		plot_col_dendr=F;
		plot_row_dendr=F;
	}

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
		mat=mat[row_dendr[["names"]],,drop=F];
	}else{
		
		cat("Adjusting Layout:\n");
		cat("  In Heatmap Height: ", heatmap_height, "\n");
		cat("  In Heatmap Width:  ", heatmap_width, "\n");

		hw_ratio=heatmap_height/heatmap_width;

		if(hw_ratio>1){
			# Height is greater than width
			smaller_height=100;
			smaller_width=floor(smaller_height/hw_ratio);
		}else{
			# Width is greater than height
			smaller_width=100;
			smaller_height=floor(smaller_width*hw_ratio);	
		}

		cat("  Layout Heatmap Height: ", smaller_height, "\n");
		cat("  Layout Heatmap Width:  ", smaller_width, "\n");

		layoutmat=matrix(
			rep(1, smaller_height*smaller_width), 
			byrow=T, ncol=smaller_width);
	
		# Original
		#layoutmat=matrix(
		#	rep(1, heatmap_height*heatmap_width), 
		#	byrow=T, ncol=heatmap_width);
	}

	#print(layoutmat);
	layout(layoutmat);

	##################################################################################################
	
	par(oma=c(col_max_nchar*.60, 0, 3, row_max_nchar*.60));
	par(mar=c(0,0,0,0));
        plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), 
		xaxt="n", yaxt="n", bty="n", xlab="", ylab="");
	mtext(title, side=3, line=1, outer=T, cex=1.5, font=2);
	mtext(subtitle, side=3, line=-.5, outer=T, cex=.8, font=3);

	axis_row_cex=min(1, 80/num_row);
	axis_col_cex=min(1, 80/num_col);

        # x-axis
        axis(side=1, at=seq(.5, num_col-.5, 1), labels=colnames(mat), las=2, 
		cex.axis=axis_col_cex, line=-1.75);
        axis(side=4, at=seq(.5, num_row-.5, 1), labels=rownames(mat), las=2, 
		cex.axis=axis_row_cex, line=-1.75);

        if(log_col){
                plot_min=log10(plot_min+.0125);
                plot_max=log10(plot_max+.0125);
        }

	if(value.cex==-1){
		value.cex=min(1, .95/strheight("X"));
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

					if(!show_leading_zero){
						if(!is.na(mat[y,x]) && mat[y,x]>-1 && mat[y,x]<1){
							if(mat[y,x]>0){
								text_lab=gsub("^0\\.", ".", text_lab);
							}else{
								text_lab=gsub("^-0\\.", "-.", text_lab);
							}
						}
					}
                                }
                                text(x-.5, y-.5, text_lab, srt=atan(num_col/num_row)/pi*180, 
					cex=value.cex, font=2);
                        }
                }
        }

	if(!suppress_grid_lines){

		calc_guidelines=function(num_cells, rev=F){
			# This function will calculate major guidelines separation that are
			# more pleasing than going with something fixed

			if(num_cells<8){
				# If fewer than 8 cells, don't bother generating them
				return(NA);
			}else{
				grps=c(4, 5, 6, 7);
				rem=num_cells%%grps

				if(any(rem==0)){
					# Take largest of 0 remainders
					use_grp=grps[max(which(rem==0))];
				}else{
					# Take largest remainder
					grp_ix=max(which(rem==max(rem)));
					use_grp=grps[grp_ix];
				}

				if(rev){
					guide_pos=seq(0,num_cells,use_grp);
				}else{
					guide_pos=seq(num_cells,0,-use_grp);
				}

				guide_pos=setdiff(guide_pos, c(0, num_cells));
				return(guide_pos);
			}
		}

		plt_usr=par()$usr;

		# Guide lines
		# <need to implement here>

		# Grid lines
		for(i in 0:num_row){
			points(x=c(0, num_col), y=c(i,i), type="l", lwd=.5, lty="dashed", col="grey75");
		}
		for(i in 0:num_col){
			points(x=c(i, i), y=c(0, num_row), type="l", lwd=.5, lty="dashed", col="grey75");
		}
	}

	##################################################################################################

	par(mar=c(0, 0, 0, 0));

	if(plot_row_dendr && plot_col_dendr){
		rdh=attributes(row_dendr[["tree"]])$height;
		cdh=attributes(col_dendr[["tree"]])$height;
		plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", 
			bty="n", xlim=c(rdh, 0));
		plot(col_dendr[["tree"]], leaflab="none",xaxt="n", yaxt="n", bty="n", ylim=c(0, cdh));
		plot(0,0, type="n", bty="n", xaxt="n", yaxt="n");
		#text(0,0, "Placeholder");
	}else if(plot_row_dendr){
		rdh=attributes(row_dendr[["tree"]])$height;
		plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", 
			bty="n", xlim=c(rdh, 0));
		#text(0,0, "Row Dendrogram");
	}else if(plot_col_dendr){
		cdh=attributes(col_dendr[["tree"]])$height;
		plot(col_dendr[["tree"]], leaflab="none", xaxt="n", yaxt="n", bty="n", ylim=c(0, cdh));
		#text(0,0, "Col Dendrogram");
	}

	par(orig.par);

}

add_sign_col=function(coeff){

	if(nrow(coeff)==0){
		return(coeff);
	}else{

		cnames=colnames(coeff);
		pval_ix=which(cnames=="Pr(>|t|)");
		pval=coeff[,pval_ix];

		sig_char=function(val){
			if(!is.null(val) && !is.nan(val) && !is.na(val)){
				if(val <= .0001){ return("***");}
				if(val <= .001 ){ return("** ");}
				if(val <= .01  ){ return("*  ");}
				if(val <= .05  ){ return(":  ");}
				if(val <= .1   ){ return(".  ");}
				return(" ");
			}else{
				return(" ");
			}
		}

		sig_arr=sapply(pval, sig_char);
		fdr=round(p.adjust(pval, method="fdr"), 4);
		fdr_sig_char=sapply(fdr, sig_char);
		#out_mat=cbind(coeff, fdr);

		fmt=function(x){
			return(formatC(x, format="f", digits=4, width=7));
		}

		out_mat=cbind(
			fmt(coeff), sig_arr, fmt(fdr), fdr_sig_char);
		colnames(out_mat)=c(cnames, "Signf", "FDR", "Signf");
		
		return(out_mat);
	}	
}

signf_char=function(val){
	if(!is.null(val) && !is.nan(val) && !is.na(val)){
		if(val <= .0001){ return("***");}
		if(val <= .001 ){ return("** ");}
		if(val <= .01  ){ return("*  ");}
		if(val <= .05  ){ return(":  ");}
		if(val <= .1   ){ return(".  ");}
	}
	return(" ");
}

###############################################################################

mask_matrix=function(val_mat, mask_mat, mask_thres=1, mask_val=0){
	# If the mask_mat is greater than the mask_thres, then replace it with mask_val
        masked_matrix=val_mat;
        masked_matrix[mask_mat>mask_thres]=mask_val;
        return(masked_matrix);
}

###############################################################################

plot_histograms=function(data_mat, title="", max_col_plots=3, max_row_plots=4){
	
        orig.par=par(no.readonly=T);

	num_var=ncol(data_mat);
	var_names=colnames(data_mat);

	par(oma=c(1,0,2,0));

	# Calculate num plots per row/col
	sqrt_num_var=ceiling(sqrt(num_var));
	num_col_plots=min(max_col_plots, sqrt_num_var);
	num_row_plots=min(max_row_plots, sqrt_num_var);
	ppp=num_col_plots*num_row_plots
	par(mfrow=c(num_row_plots, num_col_plots));
	
	# Generate plots
	for(i in 1:num_var){

		data=data_mat[,i]
		if(!is.numeric(data)){
			par(mar=c(6,4,5,1));
			cat(var_names[i], " is not numeric.\n");
			tab=table(data);
			num_categories=length(tab);
			categ_names=names(tab);
			mids=barplot(tab, names.arg="");

			bar_width=mids[2]-mids[1];
			plot_range=par()$usr;
        		plot_height=plot_range[4];
        		label_size=min(c(1,.7*bar_width/par()$cxy[1]));
        		text(
				mids-par()$cxy[1]/2, 
				rep(-par()$cxy[2]/2, num_categories), 
				categ_names, srt=-45, xpd=T, pos=4, cex=label_size);	
		}else{
			par(mar=c(3,4,5,1));
			hist(data, main="", xlab="", ylab="Counts");
		
			# Calculate quick descriptive stats
			mean=mean(data);
			qtl=quantile(data, c(0.025,0.5,0.975));
			rng=range(data);

			# Label quick stats
			title(main=sprintf("mean = %.3g, median = %.3g", 
				mean, qtl[2]), line=-.5, cex.main=.9, font.main=3);
			title(main=sprintf("95%% CI (LB = %.3g, UB = %.3g)", 
				qtl[1], qtl[3]), line=-1.4, cex.main=.9, font.main=3);
			title(main=sprintf("min = %.3g, max = %.3g", rng[1], rng[2]), 
				line=-2.3, cex.main=.9, font.main=3);
		}

		title(main=var_names[i], cex.main=1.5, font.main=2);

		# If first plot on page (new page), then label the outer margins
		if(((i-1)%%ppp)==0){
			mtext(title, side=3, outer=T, cex=1.2, font=2);
		}
	}

        par(orig.par);
}

###############################################################################

write_summary_file=function(out_mat, fname){

	cat("Writing: ", fname, "\n");
	cat(" Num Samples: ", nrow(out_mat), "\n");
	cat(" Num Categories: ", ncol(out_mat), "\n");

	totals=apply(out_mat, 1, sum);
	sumtab=cbind(rownames(out_mat), totals, out_mat);
	headers=c("sample_id", "total", colnames(out_mat));
	colnames(sumtab)=headers;
	
	write.table(sumtab, file=fname, quote=F, sep="\t", row.names=F, col.names=T);

}

###############################################################################

plot_comparisons_with_theoretical_alr=function(cnt_mat, alr_mat){
	
	# The goal is to understand whether a bimodal alr distribution is
	# true, or just a result of insufficient sequencing depth.

	# Assume that the expected abundance of each taxa is the mean 
	# across the subjects, instead of just using the mean abundance of the
	# non zero subjects.  I tried using the median, but then it looked
	# like the abundances were over estimated.

	# Calculate the alr using the mean, and 95% PI of sequencing depth
	# as a comparison to the observed alr distribution

	# Compute the multinomial distraction across all taxa, so the .5 adjustment
	# is replicated.
	
	orig.par=par(no.readonly=T);

	cat("Comparing observed alr with theoretical...\n");	
	
	sel_cat=colnames(alr_mat);
	cnt_cat=colnames(cnt_mat);

	num_sel_cat=length(sel_cat);
	num_samples=nrow(cnt_mat);

	cat("Num Samples: ", num_samples, "\n");
	cat("Num Selected Categories: ", num_sel_cat, "\n");

	# Confirm that the count matrix and alr matrix have matching categories
	missing_cat=setdiff(sel_cat, cnt_cat);
	if(length(missing_cat)){
		cat("Error: Selected categories not found in count matrix.\n");

		cat("Examples of Selected:\n");
		print(head(sel_cat));
		cat("\n");
		cat("Examples of Count Matrix:\n");
		print(head(cnt_cat));
		cat("\n");

		quit(status=-1);
	}

	# Remove any fractional adjustments to avoid 0's.
	unadj_cnt_mat=apply(cnt_mat, c(1,2), floor);
	
	# Calculate the sequencing depths
	samp_totals=apply(unadj_cnt_mat, 1, sum);
	totals_ci=quantile(samp_totals, c(.025, .5, .975));
	names(totals_ci)=c("lb", "med", "ub");

	cat("\n");
	cat("Sequencing Depths Prediction Intervals:\n");
	print(totals_ci);

	# Normalize with unadjusted to get unbiased prob
	normalized=unadj_cnt_mat/samp_totals;
	
	cat_norm_means=apply(normalized, 2, mean);

	cat("\n");
	cat("Selected Category Means:\n");
	print(cat_norm_means[sel_cat]);

	#----------------------------------------------------------------------

	generate_counts=function(sample_depth, cat_means){
		# Given means and actually sample depths, generate a list of
		# count matrices, assuming different depth scenarios

		#cat("Sample Depths:\n");
		#print(sample_depth);
		#cat("Category Means:\n");
		#print(cat_means);

		num_samp=length(sample_depth);
		num_cat=length(cat_means);

		mat_list=list();
		depth_type=c("obs_depth", "lb", "med", "ub");

		totals_ci=quantile(sample_depth, c(0.025, .5, 0.975));
		names(totals_ci)=c("lb", "med", "ub");

		for(dt in depth_type){
			if(dt=="obs_depth"){
				count_mat=matrix(0, nrow=num_samp, ncol=num_cat);
				for(i in 1:num_samp){
					count_mat[i,]=
						rmultinom(n=1, size=sample_depth[i], prob=cat_means);
				}
			}else{
				count_mat=t(rmultinom(n=num_samp, 
					size=totals_ci[dt], prob=cat_means));
			}

			colnames(count_mat)=names(cat_means);
			rownames(count_mat)=names(sample_depth);

			mat_list[[dt]]=count_mat;
		}

		return(mat_list);
	}

	#----------------------------------------------------------------------

	generate_alr=function(counts_list, sel_cat){

		lr_list=list();

		mat_names=names(counts_list);	
		for(mn in mat_names){
		
			cat("Working on: ", mn, "\n");

			cur_mat=counts_list[[mn]];
			adj_mat=cur_mat+.5;	

			norm_adj_mat=normalize(adj_mat);
			selected_mat=norm_adj_mat[,sel_cat,drop=F];
			sum_selected=apply(selected_mat, 1, sum);
			remainders=1-sum_selected;

			logratios=log(selected_mat/remainders);
			lr_list[[mn]]=logratios;
		}

		return(lr_list);
	}

	#----------------------------------------------------------------------

	count_mats_list=generate_counts(samp_totals, cat_norm_means);
	alr_list=generate_alr(count_mats_list, sel_cat);

	# Add observed alr values to list for plotting
	alr_list=c(list("observed"=alr_mat), alr_list);

	# Specify full/description for each alr_mat instance
	alr_type_title=list(
		"observed"="Observed/Actual ALR",
		"obs_depth"="Obs. Depths (Varied by Sample)",
		"lb"=paste("95% PI Low Bnd Dep (", totals_ci["lb"], ")", sep=""),
		"ub"=paste("95% PI Upp Bnd Dep (", totals_ci["ub"], ")", sep=""),
		"med"=paste("Median Depth (", totals_ci["med"], ")", sep="")
		);

	cat_range_list=list();

	par(mar=c(3, 3, 4, .5));
	par(mfrow=c(5, length(alr_list)));

	for(catname in sel_cat){

		values=c();
		for(alr_type in names(alr_list)){
			values=c(values, alr_list[[alr_type]][,catname]);
		}
			
		rng=range(values);

		for(alr_type in names(alr_list)){

			cur_val=alr_list[[alr_type]][,catname];
			mean_alr=sprintf("%5.3g", mean(cur_val));
			median_alr=sprintf("%5.3g", median(cur_val));
			mean_abund=sprintf("%4.3g", cat_norm_means[catname]);

			hist(cur_val,
				main="", xlab="", ylab="",
				breaks=seq(rng[1],rng[2], length.out=20));
			title(main=catname, cex.main=1, font.main=2, line=2)
			title(main=alr_type_title[[alr_type]], cex.main=.7, font.main=1, line=1)
			title(main=paste("ALR: mean = ", mean_alr , ", median = ", median_alr, sep=""),
				cex.main=.5, font.main=3, line=.5);

			if(alr_type=="observed"){
				title(main=paste("Abund: mean = ", mean_abund, sep=""),
					cex.main=.5, font.main=3, line=0);
			}
		}
	}

        par(orig.par);
}


###############################################################################

plot_rank_abundance=function(abundances, num_disp_max=50, category_colors=NULL, 
	ymax=NULL, title="Rank Abundance", subtitle="", exclude_Remaining=F){

        # This draws a single Rank Abundance Plot
	# The names on the abundances is how the category colors should be indexed

	# Abundance is a matrix, make it a named vector
        abund_dim=dim(abundances);
	if(!is.null(abund_dim)){
		if(abund_dim[1]==1 && abund_dim[2]>1){
			cat("Taking row as vector.\n");
			abundances=abundances[1,];
		}else if(abund_dim[2]==1 && abund_dim[1]>1){
			cat("Taking column as vector.\n");
			abundances=abundances[,1];
		}
	}

	# Remove "Remaining", if requested
	remaining_removed=F;
	if(exclude_Remaining){
		cat_names=names(abundances);
		rem_ix=("Remaining"==cat_names);
		if(any(rem_ix)){
			cat_names=setdiff(cat_names, "Remaining");
			abundances=abundances[cat_names];
			remaining_removed=T;
		}
	}

        # Sort abundances in decreasing order
        sort_ix=order(abundances, decreasing=T);

        # Grab abundances/names of top categories
        num_top_categories=min(num_disp_max, length(abundances));
        top=abundances[sort_ix[1:num_top_categories]];
        cat_names=names(top);
	if(is.null(cat_names)){
		cat_names=sprintf("cat_%02i",1:num_top_categories);
	}

        # Remove labels if abundance is 0
        cat_names[top==0]="";

        # Grab colors by names in top
	if(!is.null(category_colors)){
		if(length(category_colors)==1){
			reordered_colors=rep(category_colors, num_top_categories);
		}else{
			reordered_colors=category_colors[cat_names];
		}
	}else{
		reordered_colors=rep("white", num_top_categories);	
	}

	# If Remaining not removed, color it
	rem_ix=(cat_names=="Remaining");
	if(any(rem_ix)){
		reordered_colors[rem_ix]="grey";
	}

        # Calc percent represented
        prop_draw=sum(top);
        rep_title=paste("Percentage Represented: ",
                        round(prop_draw*100.0, 2),
                        "%");

	# If ymax not defined, maximize it based on data 
	if(is.null(ymax)){
		ymax=max(top);
	}

	# Generate plot and plot labels
        mids=barplot(top, col=reordered_colors, names.arg="", ylim=c(0, ymax*1.05));
        title(main=title, font.main=2, cex.main=1.75, line=-1);
        title(main=subtitle, font.main=3, cex.main=1, line=-2)
        title(main=rep_title, font.main=1, cex.main=1.1, line=-3)
	if(remaining_removed){
		title(main="(\"Remaining\" Category Found and Removed)", 
			font.main=2, col.main="red", cex.main=.95, line=-4);
	}

        # Compute and label categories beneath barplot
        bar_width=mids[2]-mids[1];
        plot_range=par()$usr;
        plot_height=plot_range[4];
        label_size=min(c(1,.7*bar_width/par()$cxy[1]));
        text(mids-par()$cxy[1]/2, rep(-par()$cxy[2]/2, num_top_categories), cat_names, 
		srt=-45, xpd=T, pos=4, cex=label_size);

}

###############################################################################

modify_filenames=function(infn, modlist, usage=F){
# This function implements an abridged version of the C-shell file name modifier.
# To include the usage information in the usage of a script using it, make the usage=T.
# Otherwise example call:  modify_filenames("dir1/dir2/filename.ext1.ext2", ":t:r:r");

        if(usage==T){

                usage_str=c(
                        "   Only the following C-shell filename modifiers are implemented:\n",
                        "     If dir1/dir2/file.ext1.ext2 is the specified file name then:\n",
                        "       :t      tail    file.ext1.ext2\n",
                        "       :r      root    dir1/dir2/file.ext1\n",
                        "\n",
                        "     To get 'file', you would specify :t:r:r\n"
                );

                return(paste(usage_str, collapse="", sep=""));

        }else{

		if(is.null(modlist) || modlist=="" || modlist=="NA"){
			return(infn);
		}

                mod_arr=strsplit(modlist, ":")[[1]];
                mod_arr=mod_arr[mod_arr!=""];
                num_mods=length(mod_arr);

                for(i in 1:num_mods){

                        if(mod_arr[i]=="r"){
                                toks=head(strsplit(infn, "\\.")[[1]], -1);
                                infn=paste(toks, collapse=".");
                        }else if(mod_arr[i]=="t"){
                                toks=tail(strsplit(infn, "\\/")[[1]], 1);
                                infn=toks;
                        }else{
                                cat("Unsupported modifier: ", mod_arr[i], "\n");
                                quit(status=-1);
                        }

                }
                return(infn);
        }
}

###############################################################################

adjust_positions=function(orig_pos, min_spacing, max_pos=Inf, min_pos=-Inf, max_tries=10000, verbose=F){
# This function will adjust the values in the orig_pos so that the spacing between values is
# greater than min_spacing

	num_pos=length(orig_pos);
	names(orig_pos)=1:num_pos;

	if(verbose){
		cat("Requested Min Spacing: ", min_spacing, "\n");
		cat("Limits: min=", min_pos, ", max=", max_pos, "\n");
		cat("Original:\n");
		print(orig_pos);
		cat("\n");
	}
	
	order_ix=order(orig_pos);
	sorted_pos=orig_pos[order_ix];

	if(verbose){
		cat("Sorted:\n");
		print(sorted_pos);
		cat("\n");
	}
	
	for(tries in 1:max_tries){

		adjmade=F;

		for(ix in 1:(num_pos-1)){
			if((sorted_pos[ix+1]-sorted_pos[ix])<min_spacing){

				sorted_pos[ix]=sorted_pos[ix]-min_spacing/4;
				if(sorted_pos[ix]<min_pos){
					sorted_pos[ix]=min_pos;
				}

				sorted_pos[ix+1]=sorted_pos[ix+1]+min_spacing/4;
				if(sorted_pos[ix+1]>max_pos){
					sorted_pos[ix+1]=max_pos;
				}

				adjmade=T;
			}
		}

		if(adjmade==F){
			if(verbose){
				cat("Adjustment Trials Made: ", tries, "\n");
			}
			break;
		}
	}

	unsort_ix=order(as.numeric(names(sorted_pos)));
	unsorted=sorted_pos[unsort_ix];
	
	if(verbose){
		cat("Adjusted Sorted:\n");
		print(sorted_pos);
		cat("\n");
		cat("Adjusted Unsorted:\n");
		print(unsorted);
		cat("\n");
	}
	return(unsorted);
	

}

connect_orig_to_adj_wline=function(orig_x, orig_y, adj_x, adj_y, lwd=1, col="grey25", lty="solid"){

	for(i in 1:length(orig_x)){
		points(c(orig_x[i], adj_x[i]), c(orig_y[i], adj_y[i]), lwd=lwd, col=col, lty=lty,
			type="l");
	}

}
