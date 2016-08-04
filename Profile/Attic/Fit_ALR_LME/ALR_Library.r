#!/usr/bin/env Rscript

cat("Loading ALR Libraries...\n");

load_factors=function(fname){
	factors=data.frame(read.table(fname, sep="\t", header=TRUE, row.names=1, check.names=FALSE));
	factor_names=colnames(factors);
	cat("Num Samples in Factor File: ", nrow(factors), "\n", sep="");
	cat("Num Factors in Factor File: ", ncol(factors), "\n", sep="");
	cat("\n");
	return(factors);
}

load_summary_file=function(fname){
	inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", quote="", row.names=1))
	counts_mat=inmat[,2:(ncol(inmat))];
	cat("   Num Samples in Summary Table: ", nrow(counts_mat), "\n", sep="");
	cat("Num Categories in Summary Table: ", ncol(counts_mat), "\n", sep="");
	cat("\n");
	return(counts_mat);
}

load_reference_levels_file=function(fname){
        inmat=as.matrix(read.table(fname, sep="\t", header=F, check.names=FALSE, comment.char="#", row.names=1))
        colnames(inmat)=c("ReferenceLevel");
        print(inmat);
        cat("\n");
        if(ncol(inmat)!=1){
                cat("Error reading in reference level file: ", fname, "\n");
                quit(status=-1);
        }
        return(inmat);
}

relevel_factors=function(factors, ref_lev_mat){
        num_factors_to_relevel=nrow(ref_lev_mat);
        relevel_names=rownames(ref_lev_mat);
        for(i in 1:num_factors_to_relevel){
                tmp=factors[,relevel_names[i]];
                #print(tmp);
                tmp=relevel(tmp, ref_lev_mat[i, 1]);
                #print(tmp);
                factors[,relevel_names[i]]=tmp;
        }
        return(factors);
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

extract_top_categories=function(ordered_normalized, top){

	num_samples=nrow(ordered_normalized);
	num_categories=ncol(ordered_normalized);

	cat("Samples: ", num_samples, "\n");
	cat("Categories: ", num_categories, "\n");
	
	num_saved=min(c(num_categories, top+1));

	cat("Top Requested to Extract: ", top, "\n");
	cat("Columns to Extract: ", num_saved, "\n");

	top_cat=matrix(0, nrow=num_samples, ncol=num_saved);
	top=num_saved-1;

	# Extract top categories requested
	top_cat[,1:top]=ordered_normalized[,1:top];

	# Included remaineder as sum of remaining categories
	top_cat[,(top+1)]=apply(
		ordered_normalized[,(top+1):num_categories, drop=F],
		1, sum);

	rownames(top_cat)=rownames(ordered_normalized);
	colnames(top_cat)=c(colnames(ordered_normalized)[1:top], "Remaining");

	return(top_cat);
			
}

collapse_to_top_matrix=function(matrix, num_top, remove_remaining_cat_from_top){
	
	# Order categories from largest abundance to smallest
	norm_mat=normalize(matrix);
	categ_prop=apply(norm_mat, 2, mean);
	order_ix=order(categ_prop, decreasing=T);
	ord_mat=matrix[,order_ix];	

	category_names=colnames(ord_mat);	
	num_categories=length(category_names);

	# Find Remaining column by name
	uc_cat_names=toupper(category_names);
	remaining_ix=which(uc_cat_names=="REMAINDER" | uc_cat_names=="REMAINING");
	if(length(remaining_ix)!=1){
		cat("*******************************************************\n");
		cat("* Remainder/Remaining column not found.               *\n");
		cat("*******************************************************\n");
		cat("\n");
		remaining_col_found=F;
	}else{
		remaining_col_found=T;
	}

	if(remaining_col_found){
		# Move remaining column to end
		remaining_col=ord_mat[, remaining_ix, drop=F];
		matrix_wo_rem=ord_mat[, -remaining_ix, drop=F];
		recomb_mat=cbind(matrix_wo_rem, remaining_col);
		category_names=colnames(recomb_mat);
	}else{
		recomb_mat=ord_mat;
	}
	
	# Keep at least one category for the remainder/denominator
	if(num_top+1>num_categories){
		num_top=num_categories-1;
	}

	# Extract top
	top_col=recomb_mat[, 1:num_top, drop=F];
	rem_col=recomb_mat[, (num_top+1):num_categories, drop=F];
	
	# Merge remainder
	sum_rem_col=apply(rem_col, 1, sum);

	# Recombine
	out_mat=cbind(top_col, sum_rem_col)
	top_col_names=colnames(top_col);
	colnames(out_mat)=c(top_col_names, "Remaining");

	out=list();
	out[["remaining"]]="Remaining";
	out[["matrix"]]=out_mat;

	return(out);

}

additive_log_ratio=function(ordered_matrix, denominator_category){

	num_cat=ncol(ordered_matrix);
	num_samp=nrow(ordered_matrix);

	denominator=ordered_matrix[,denominator_category];
	alr_mat=matrix(0, nrow=num_samp, ncol=(num_cat-1));
	
	for(i in 1:num_samp){
		alr_mat[i,]=log(ordered_matrix[i,1:(num_cat-1)]/denominator[i]);
		#print(alr_mat[i,]);
	}

	rownames(alr_mat)=rownames(ordered_matrix)
	colnames(alr_mat)=head(colnames(ordered_matrix), num_cat-1);

	alr_struct=list();
	alr_struct[["transformed"]]=alr_mat;
	alr_struct[["denominator"]]=denominator;

	return(alr_struct);
}

plot_text=function(strings){
	par(family="Courier");
	par(oma=rep(.5,4));
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

reconcile_samples=function(factor_table, summary_table){

	cat("Reconciling sample IDs between Factor Table and Summary Table...\n");

	factor_sample_ids=rownames(factor_table);
	sumtab_sample_ids=rownames(summary_table);

	num_factor_sample_ids=length(factor_sample_ids);
	num_sumtab_sample_ids=length(sumtab_sample_ids);

	#print(factor_sample_id);
	#print(counts_sample_id);

	shared_sample_ids=intersect(factor_sample_ids, sumtab_sample_ids);
	num_shared_sample_ids=length(shared_sample_ids);

	cat("Num summary table sample IDs: ", num_sumtab_sample_ids, "\n");
	cat("       Num factor sample IDs: ", num_factor_sample_ids, "\n");
	cat("       Num shared sample IDs: ", num_shared_sample_ids, "\n");
	cat("\n");

	cat("Samples missing from Summary Table:\n");
	print(setdiff(factor_sample_ids, sumtab_sample_ids));
	cat("\n");

	cat("Samples missing from Factor Table:\n");
	print(setdiff(sumtab_sample_ids, factor_sample_ids));
	cat("\n");

	cat("Total samples shared: ", num_shared_sample_ids, "\n");
	cat("\n");

	shared_sample_ids=sort(unique(shared_sample_ids));

	reconciled=list();
	reconciled[["counts"]]=summary_table[shared_sample_ids,,drop=F];
	reconciled[["factors"]]=factor_table[shared_sample_ids,,drop=F];

	return(reconciled);
}


# Assign 0's to values smaller than smallest abundance across entire dataset
reset_zeros=function(normalized_st){
	min_assay=min(normalized_st[normalized_st!=0]);
	zero_replacment=min_assay/10;

	cat("Lowest non-zero value: ", min_assay, "\n", sep="");
	cat("Replacing 0's with: ", zero_replacment, "\n", sep="");
	cat("\n");

	nozero_mat=normalized_st;
	nozero_mat[normalized_st==0]=zero_replacment;
	return(nozero_mat);
}

###############################################################################

plot_trellis=function(data, group_data, series_data, trt_levels=NULL, fits=NULL){

        orig_par=par(c("mfrow","mar","oma"));

	# Get unique groups (patient ID)  and series (time points)
        groups=unique(group_data[,1]);
        series=unique(series_data[,1]);

	# Get names of columns, e.g. patient id, time, ALR, respectively
        group_name=colnames(group_data);
        series_name=colnames(series_data);
        data_name=colnames(data);

        #cat("Group Name: ", group_name, "\n");
        #print(groups);
        #cat("\n");
        #cat("Series Name: ", series_name, "\n");
        #print(series);

	# Convert data frames (with column names) to vectors
        series_vect=series_data[,1];
        group_vect=group_data[,1];
        data_vect=data[,1];
        trt_levels_vect=trt_levels[,1];
        sample_names=names(data_vect);

	# Get the x (series/time range) and y (alr) limits for the plot
        data_range=range(data_vect);
        series_range=range(series);

	# Compute frequency of data ticks for the Y axis
        data_ticks=seq(data_range[1], data_range[2], length.out=7);
        data_ticks=data_ticks[2:(length(data_ticks)-1)];
        data_ticks_lab=sprintf("%2.3f", data_ticks);

	# Compute frequency of times series ticks for the X axis
        series_ticks=seq(series_range[1], series_range[2], length.out=5);
        series_ticks_lab=sprintf("%2.0f", series_ticks);

	# Compute the dimensions for the number of trellis cells/subplots
	# based on the number of groups/patient IDs
        num_groups=length(groups);
        plot_cells=sqrt(num_groups);
        trel_rows=floor(plot_cells);
        trel_cols=ceiling(plot_cells);;

	# Set up plot margins
        par(mfrow=c(trel_rows,trel_cols));
        par(mar=c(0,0,0,0));
        par(oma=c(.5,4.5,7,.5));


	# sort groups by treatment level
	groups_order=groups;
	if(!is.null(trt_levels)){
		dedup_ix=!duplicated(group_vect);
		dedup_group=group_vect[dedup_ix];
		dedup_trt_levels=trt_levels_vect[dedup_ix];
		dedep_trt_order=order(dedup_trt_levels);
		groups_order=groups_order[dedep_trt_order];
	}

        plot_ix=0;
	# Generate sub plot on a per group/patient id basis
        for(group in groups_order){
	
		# Get index of values based on group
                grp_ix=which(group==group_data);
                #cat("(", group, ") Group Index: ", grp_ix, "\n");

		# Extract x/y based on group index
                xpts = series_vect[grp_ix];
                ypts = data_vect[grp_ix];

		# Order the points by increasing time so lines can be connected
                ix=order(xpts, decreasing=F);
                xpts=xpts[ix];
                ypts=ypts[ix];

		# Plot the points
                plot(xpts, ypts, type="p",
                        xlab=series_name, ylab=data_name,
                        ylim=data_range, xlim=series_range,
                        xaxt="n", yaxt="n", cex=.5);

		# Plot 0 reference line
                abline(h=0, col="blue", lty=6);

		# Label the name of the group/patient ID in each subplot
                title(main=group, line=-.75, cex.main=.75);

		# Label the Y-axis if its the first cell in the row
                if(!(plot_ix %% trel_cols)){
                        axis(side=2, at=data_ticks, labels=data_ticks_lab, cex.axis=.5, las=2);
                }

		# Label the X-axis if it is the first row
                if(plot_ix < trel_cols){
                        axis(side=3, at=series_ticks, labels=series_ticks_lab, cex.axis=.5, las=2);
                }

                # Fit nlm curve if it is available
                if(!is.null(fits)){
                        points(fits[[group]][,1], fits[[group]][,2], col="blue", type="l");
                }

		# Keep track of the plot index
                plot_ix=plot_ix+1;
        }

	# Label axis overall 
        mtext(data_name, side=3, line=5, outer=T, font=2, cex=.8);
        mtext(series_name, side=3, line=2, outer=T);
        mtext(data_name, side=2, line=3, outer=T);

	# Restore plot parameters
        par(orig_par);
}

###############################################################################

plot_text=function(strings, title="", orientation="P"){

	orig_par=par(c("family", "oma", "mar", "mfrow"));	

        par(family="Courier");
        par(oma=rep(1,4));
        par(mar=rep(0,4));
	par(mfrow=c(1,1));

	# Determine number of lines per page
        num_lines=length(strings);
	font_size=11;
	if(orientation=="L"){
		lines_per_page=40; # Lines per page in lanscape 11pt font
	}else{
		lines_per_page=53; # Lines per page in portrait 11pt font
	}
	if(num_lines > lines_per_page){
		font_size=font_size*lines_per_page/num_lines;
		lines_per_page=num_lines;
	}

	# Set up page dimensions, i.e. lines per page
        top=lines_per_page;
        plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  
		xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(.5,.5,.5,.5), mar=c(0,0,0,0)
                );

	# Determine how to scale font size
	font_to_cex_ratio=10;
	text_size=font_size/font_to_cex_ratio;

	# Write text to page
        for(i in 1:num_lines){
                #cat(strings[i], "\n", sep="");
                strings[i]=gsub("\t", "   ", strings[i]);
                text(0, top-i, strings[i], pos=4, cex=text_size);
        }

        # Write the title 
        mtext(title, outer=TRUE, side=3, font=2, cex=.85, family="");

	# Restore original plotting parameters
	par(orig_par);
}

###############################################################################

plot_heatmap=function(mat, pval_mat=NULL, title="", noPrintZeros=F, guideLines=F, max_name_len=30){

        cat("Plotting: ", title, "\n");

	# Make sure matrix is valid
        if(is.null(dim(mat))){
                cat(title, " Matrix is NULL.  No heatmap generated.\n");
                return;
        }
	mat=t(mat);
	if(!is.null(pval_mat)){
		pval_mat=t(pval_mat);
	}

	# Set plotting parameters
	orig_par=par(c("family", "oma", "mar"));
        par(family="mono");
        par(oma=c(10, 10, 1, .5));
        par(mar=c(6.1, 5.1, .5, .5));

        # Generate colors from red to blue
        colors=(rainbow(2^16, start=0, end=0.65));

	# Plot heatmap.
        # Remember that rows and columns are reversed in the image
        image(1:nrow(mat),1:ncol(mat), mat,
                xaxt="n", yaxt="n",
                xlab="", ylab="",
                col=colors
        );


	#----------------------------------------------------------------------
	# Get names from matrix

	cnames=colnames(mat);
	rnames=rownames(mat);

	# Truncate from end
	cnames=substring(cnames, nchar(cnames)-max_name_len);
	rnames=substring(rnames, nchar(rnames)-max_name_len);

	# Pad a space at end for aesthetics
        cnames=paste(cnames, " ", sep="");
        rnames=paste(rnames, " ", sep="");

	#print(cnames);
	#print(rnames);

        # Find length of longest string in each column or row label
        cname_max_len=max(nchar(cnames));
        rname_max_len=max(nchar(rnames));

        # Get the number of rows and columns
        ncols=ncol(mat);
        nrows=nrow(mat);

	# Estimate scale base on name lengths
        cscale=min(c(30/cname_max_len, 55/ncols));
        rscale=min(c(30/rname_max_len, 55/nrows));

	# Set minimum scaling
        cscale=min(1, cscale);
        rscale=min(1, rscale);

        max_width=max(nchar(sprintf("%.2f",mat)));
        cell_cex=sqrt(min(c(cscale, rscale))^2);
        #cell_cex=(1/max_width)*sqrt(min(c(cscale, rscale))^2);

	# Label the cells with the underlying values
        for(i in 1:nrow(mat)){
                for(j in 1:ncol(mat)){

                        if(!is.na(mat[i,j]) && (noPrintZeros && mat[i,j]==0)){
                                # Skip
                        }else{

				str=sprintf("%.2f",mat[i,j]);
				str=gsub("^0\\.",".", str);

				if(is.null(pval_mat)){
					text(i,j,labels=str, cex=cell_cex, srt=45);
				}else{

					if(pval_mat[i, j] <= .01){
						sig="*";
					}else if(pval_mat[i, j] <= .05){
						sig="\"";
					}else if(pval_mat[i, j] <= .1){
						sig="'";
					}else{
						sig="";
					}

					text(i,j,labels=paste(str, sig, sep=""), cex=cell_cex, srt=45);
	
					
				}
                        }
                }
        }

	#----------------------------------------------------------------------

        # Plot guidelines
        if(guideLines){

                splits=c(2,3,4,5);
                h_remainder=ncols %% splits;
                best_h_split=splits[max(which(h_remainder==min(h_remainder)))];
                if(ncols>best_h_split){
                        h_line_pos=seq(best_h_split, ncols, best_h_split)+.5;
                        abline(h=h_line_pos, col="black", lty="dashed");
                        abline(h=h_line_pos, col="white", lty="dotted");
                }

                v_remainder=nrows %% splits;
                best_v_split=splits[max(which(v_remainder==min(v_remainder)))];
                if(is.na(best_v_split)){
                        best_v_split=max(splits);
                }
                if(nrows>best_v_split){
                        v_line_pos=seq(best_v_split, nrows, best_v_split)+.5;
                        abline(v=v_line_pos, col="black", lty="dashed");
                        abline(v=v_line_pos, col="white", lty="dotted");
                }

        }

        # Plot the labels
        mtext(cnames, at=1:ncols, side=2, las=2, cex=cscale);
        mtext(rnames, at=1:nrows, side=1, las=2, cex=rscale);

        # Plot the title
        mtext(title, line=0, at=nrows*.5, side=3, font=2, family="");

	# Restore graphics parameters
	par(orig_par);

}

cat("done.\n");
