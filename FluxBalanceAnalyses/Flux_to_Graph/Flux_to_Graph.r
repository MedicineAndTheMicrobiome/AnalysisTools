#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

options(useFancyQuotes=F);
options(digits=5)

params=c(
	"fba_infile", "f", 1, "character",
	"outputroot", "o", 1, "character",
	"zero_cutoff", "z", 2, "numeric"
);

ZERO_TOLERANCE=1e-6;

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

script_path=paste(head(strsplit(script_name, "/")[[1]], -1), collapse="/");

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <Combined FBA File>\n",
	"	-o <output filename root.\n",
	"\n",
	"	[-z <cutoff for 0, default=", ZERO_TOLERANCE, ">\n",
	"\n",
	"This script will read in the FBA results.\n",
	"Each row is a different taxa of the community\n",
	"The first column contains the sample ID, and each column\n",
	"Contains the flux values for an exchange reaction.\n",
	"\n");

if(
	!length(opt$fba_infile) || 
	!length(opt$outputroot)
){
	cat(usage);
	q(status=-1);
}

FBA_InFile=opt$fba_infile;
OutputFnameRoot=opt$outputroot;
ZeroCutoff=ZERO_TOLERANCE;

if(length(opt$zero_cutoff)){
	ZeroCutoff=opt$zero_cutoff;
}


param_text=capture.output({
	cat("\n");
	cat("Combined Flux File: ", FBA_InFile, "\n");
	cat("Output File Name Root: ", OutputFnameRoot, "\n");
	cat("Zero-Cutoff: ", ZeroCutoff, "\n");
	cat("\n");
});

print(param_text, quote=F);

###############################################################################

plot_text=function(strings){

	orig.par=par(no.readonly=T);

	par(mfrow=c(1,1));
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

	par(orig.par);
}

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, deci_pts=4,
        label_zeros=T, counts=F, value.cex=1,
        plot_col_dendr=F,
        plot_row_dendr=F
){

        num_row=nrow(mat);
        num_col=ncol(mat);

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

                mat=mat[row_dendr[["names"]], col_dendr[["names"]]];

        }else if(plot_col_dendr){
                layoutmat=matrix(
                        c(
                        rep(rep(2, heatmap_width), col_dend_height),
                        rep(rep(1, heatmap_width), heatmap_height)
                        ), byrow=T, ncol=heatmap_width);

                col_dendr=get_dendrogram(mat, type="col");
                mat=mat[, col_dendr[["names"]]];

        }else if(plot_row_dendr){
                layoutmat=matrix(
                        rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height),
                        byrow=T, ncol=row_dend_width+heatmap_width);

                row_dendr=get_dendrogram(mat, type="row");
                mat=mat[row_dendr[["names"]],];
        }else{

		if(heatmap_height*heatmap_width < 10000){
			layoutmat=matrix(
				rep(1, heatmap_height*heatmap_width),
				byrow=T, ncol=heatmap_width);
		}else{
			larger_dim=max(heatmap_height, heatmap_width);

			hm_height_norm=heatmap_height/larger_dim;	
			hm_width_norm=heatmap_width/larger_dim;

			lo_height=100*hm_height_norm;
			lo_width=100*hm_width_norm;

			layoutmat=matrix(
				rep(1, lo_height*lo_width),
				byrow=T, ncol=lo_width);

		}
			
        }

        #print(layoutmat);
        layout(layoutmat);

        ##################################################################################################

        par(oma=c(col_max_nchar*.60, 0, 3, row_max_nchar*.60));
        par(mar=c(0,0,0,0));
        plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), xaxt="n", yaxt="n", bty="n", xlab="", ylab="");
        mtext(title, side=3, line=0, outer=T, font=2);

        # x-axis
        axis(side=1, at=seq(.5, num_col-.5, 1), labels=colnames(mat), las=2, line=-1.75);
        axis(side=4, at=seq(.5, num_row-.5, 1), labels=rownames(mat), las=2, line=-1.75);

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

	calc_guidelines=function(num_cells, rev=F){
		if(num_cells<8){
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

	abline(h=1:(num_row-1), lwd=.125, lty="dotted", col="grey95");
	abline(h=calc_guidelines(num_row), lty="dotted", lwd=2, col="white");
	abline(h=calc_guidelines(num_row), lty="dotted", lwd=1, col="black");

	abline(v=1:(num_col-1), lwd=.125, lty="dotted", col="grey95");
	abline(v=calc_guidelines(num_col, rev=T), lty="dotted", lwd=2, col="white");
	abline(v=calc_guidelines(num_col, rev=T), lty="dotted", lwd=1, col="black");

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

plot_title_page=function(title, subtitle="", title_cex=3){

        orig.par=par(no.readonly=T);
        par(family="serif");
        par(mfrow=c(1,1));

        plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n",  xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                );

        # Title
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

##############################################################################
# Main Program Starts Here!
##############################################################################

pdf(paste(OutputFnameRoot, ".pca.pdf", sep=""), height=11, width=8.5);

##############################################################################

load_community_fba=function(fname, host="Vm"){

	raw_mat=as.data.frame(read.delim(fname,  header=TRUE, check.names=FALSE, sep="\t"));

	headers=colnames(raw_mat);

	cat("\n");
	cat("First headers:\n");
	print(head(headers));
	cat("...\n");
	cat("\n");

	# Load Sample IDs, 1st Column
	unique_sample_ids=sort(unique(raw_mat[,"SampleID"]));
	num_unique_sample_ids=length(unique_sample_ids);
	cat("Num Sample IDs:", num_unique_sample_ids, "\n");
	cat("Sample IDs (examples):\n");
	print(head(unique_sample_ids));
	cat("...\n");
	cat("\n");

	# Load Taxa IDs, 2nd Column
	unique_taxa_ids=sort(unique(raw_mat[,"TaxaID"]));
	num_unique_taxa_ids=length(unique_taxa_ids);
	cat("Num Taxa IDs:", num_unique_taxa_ids, "\n");
	cat("Taxa IDs (Beginning and End):\n");
	print(head(unique_taxa_ids));
	cat("...\n");
	print(tail(unique_taxa_ids));
	cat("\n");

	# Remove host/environment
	cat("Removing host/environment (", host, ") rows...\n");
	host_ix=raw_mat[,"TaxaID"]==host;	
	no_host_mat=raw_mat[!host_ix,];

	cat("------------------------------------------------------------------\n");
	cat("Num Samples: ", num_unique_sample_ids, "\n");
	cat("Num Taxa: ", num_unique_taxa_ids, "\n");
	cat("Num Reactions: ", ncol(raw_mat)-2, "\n");
	cat("\n");
	cat("Returned dimensions:", paste(dim(no_host_mat), collapse=" x "), "\n");
	cat("------------------------------------------------------------------\n");

	return(no_host_mat);
}


filter_low_flux_categories=function(raw_fb_mat, zero_cutoff=1e-6){

	# The raw_fb_mat still has the first 2 columns as sample_id and taxa 

	cat("Screening flux categories at cutoff: ", zero_cutoff, "\n");

	# Extract flux values
	ncol_raw_mat=ncol(raw_fb_mat);
	flux_val_mat=raw_fb_mat[,3:ncol_raw_mat];

	nrow_flux_val=nrow(flux_val_mat);
	cat("Number of Raw Flux Categories:", nrow_flux_val, "\n");

	# Compute magnitude (absolute value) of abundances and keep reactions > cutoff
	mag_flux_mat=abs(flux_val_mat);
	max_flux_arr=apply(mag_flux_mat, 2, max);
	keep_flux_ix=(max_flux_arr>ZeroCutoff);
	num_kept=sum(keep_flux_ix);

	cat("Num Screened Categories: ", num_kept, "\n");

	# Keep essentially non-zero fluxes
	kept_flux_mat=flux_val_mat[,keep_flux_ix];

	# Set fluxes of kept below cutoff to 0.
	pos_mat=kept_flux_mat;
	neg_mat=kept_flux_mat;
	pos_mat[pos_mat<zero_cutoff]=0;
	neg_mat[neg_mat>-zero_cutoff]=0;
	comb_mat=pos_mat+neg_mat;

	out_mat=cbind(raw_fb_mat[,c(1,2)], comb_mat);

	#par(mfrow=c(1,3));
	#brs=log10(range(mag_flux_mat));
	#hist(log10(as.matrix(mag_flux_mat)), main="Original Unfiltered", breaks=brs);
	#hist(log10(as.matrix(abs(kept_flux_mat))), main="Kept", breaks=brs);
	#hist(log10(as.matrix(abs(out_mat))), main="Kept and Filtered", breaks=brs);
	

#	return(graph_list);
	return(out_mat);
}

threshold_flux_by_magnitude=function(raw_fb_mat, prop_to_keep=.95){
	
	cat("Sorting fluxes in decreasing order.\n");

	ncol_raw_mat=ncol(raw_fb_mat);
	flux_val_mat=raw_fb_mat[,3:ncol_raw_mat];
	num_flux_categories=ncol(flux_val_mat);

	# Calculate magnitues of flux
	mag_flux_mat=abs(flux_val_mat);
	sum_mag_flux_perrxn=apply(mag_flux_mat, 2, sum);
	sum_mag_flux_acrrxn=sum(sum_mag_flux_perrxn);

	# Order the reactions by decreasing flux magnitudes
	mag_flux_dec_ord_ix=order(sum_mag_flux_perrxn, decreasing=T);
	mag_flux_dec_ord=sum_mag_flux_perrxn[mag_flux_dec_ord_ix];
	
	# Normalize by total flux across all reactions
	normalized_flux=mag_flux_dec_ord/sum_mag_flux_acrrxn
	# Calculate cumsum
	cumsum_mag_flux=cumsum(normalized_flux);

	# Create table
	summary_mat_stat_names=c("Index", "FluxMag", "Prop", "Cumul");
	summary_mat=matrix(0, nrow=num_flux_categories, ncol=length(summary_mat_stat_names));
	colnames(summary_mat)=summary_mat_stat_names;

	summary_mat[,"Index"]=1:num_flux_categories;
	summary_mat[,"FluxMag"]=mag_flux_dec_ord;
	summary_mat[,"Prop"]=normalized_flux;
	summary_mat[,"Cumul"]=cumsum_mag_flux;

	rownames(summary_mat)=names(mag_flux_dec_ord);
	print(summary_mat);

	# Extract flux from 1 to less than propt_to_keep
	num_flux_to_keep=sum(cumsum_mag_flux<prop_to_keep)+1;
	flux_names_to_keep=rownames(summary_mat[1:num_flux_to_keep,,drop=F]);

	cat("Flux Names to Keep:\n");
	print(flux_names_to_keep);
	keep_matrix=cbind(raw_fb_mat[,1:2], flux_val_mat[,flux_names_to_keep,drop=F]);

	if(0){
		for(i in 1:num_flux_to_keep){
			hist(flux_val_mat[flux_val_mat[,i]!=0,i],breaks=20);
		}
	}

	return(keep_matrix);
}

#------------------------------------------------------------------------------

nz_vec_cor=function(x, y){
	
	nz_x_ix = (x!=0);
	nz_y_ix = (y!=0);
	shared_ix = (nz_x_ix | nz_y_ix);

	#num_shared=sum(shared_ix);
	#cat(num_shared, "\n");

	cor_val=cor(x[shared_ix], y[shared_ix], method="spearman");

	return(cor_val);

}

nz_mat_cor=function(mat){
	
	num_var=ncol(mat);

	out_mat=matrix(NA, nrow=num_var, ncol=num_var);
	diag(out_mat)=1;
	
	for(i in 1:num_var){
		for(j in 1:(i-1)){
			nzcor=nz_vec_cor(mat[,i], mat[,j]);
			out_mat[i,j]=nzcor;
			out_mat[j,i]=nzcor;
		}
	}

	#print(out_mat);

	return(out_mat);

}

calculate_flux_correlation=function(raw_fb_mat){

	ncol_raw_mat=ncol(raw_fb_mat);
	flux_val_mat=raw_fb_mat[,3:ncol_raw_mat];
	num_flux_categories=ncol(flux_val_mat);

	cor_mat=nz_mat_cor(flux_val_mat);
	
}

project_PSD=function(cor_mat, tol=1e-10){

	make_sym=function(A){
		return((A+t(A))/2);
	}

	# Force symmetry
	cor_mat=make_sym(cor_mat);

	# Eigen Val
	eigen_rec=eigen(cor_mat, symmetric=T);
	evals=eigen_rec$values;
	evecs=eigen_rec$vectors;
	
	# Truncate negative and small values to zero
	evals[evals<tol]=0;

	# Reassemble as PSD matrix
	cov_mat=evecs %*% diag(evals) %*% t(evecs);
	
	# Force symmetry
	cov_mat=make_sym(cov_mat);

	# Convert from cov to cor by normalizing
	diag_scal_mat=diag(1/sqrt(diag(cov_mat)));
	psd_cor_mat=diag_scal_mat %*% cov_mat %*% diag_scal_mat;

	return(psd_cor_mat);

}

#------------------------------------------------------------------------------

run_pca=function(cormat, val, prop_to_keep=0.95){

	# Run PCA, find top PCs
	eigen_rec=eigen(cormat, symmetric=T);

	# Normalize PC variance and find top to keep
	pca_propvar=eigen_rec$values/sum(eigen_rec$values);
	pca_propcumsum=cumsum(pca_propvar);
	keep_ix=pca_propcumsum>prop_to_keep;
	num_pc_at_cutoff=sum(keep_ix);

	# The score are the new sample coordinates: V(n x p) %*% E(p x p) = S (n x p)
	scores=(scale(val[3:ncol(val)], center=T, scale=T) %*% eigen_rec$vectors);

	cat("\n");
	cat("Cumulative Variance:\n");
	print(pca_propcumsum);
	cat("\n");

	num_kept=sum(pca_propcumsum<=prop_to_keep)+1;
	acquired_variance=pca_propcumsum[num_kept];
	
	cat("Num PCs to targeted variance (", prop_to_keep, "):", num_kept, "\n");
	cat("Acquired variance = ", acquired_variance, "\n");

	par(mar=c(4,4,4,5));

	barplot(pca_propcumsum[1:num_kept], xlim=c(0,num_kept+1), ylim=c(0,1), 
		names.arg=1:num_kept,
		main="Cumulative Variance Captured by PCs");


	results=list();
	results[["NumPCsKept"]]=num_kept;
	results[["PropKeptVar"]]=acquired_variance;
	results[["PCPropVar"]]=pca_propvar[1:num_kept];
	results[["Scores"]]=scores[,1:num_kept];

	return(results);

}

#------------------------------------------------------------------------------

select_pca_reps=function(pca_res, raw_fb_mat, max_top=15, report_fn=NULL){
	
	num_pcs=pca_res[["NumPCsKept"]];

	vals=raw_fb_mat[3:ncol(raw_fb_mat)];


	num_fluxes=ncol(vals);
	max_top=min(max_top, num_fluxes);

	cat("Num Fluxes: ", num_fluxes, "\n");
	cat("Num Kept PCs: ", num_pcs, "\n");
	select_list=list();
	selected=c();

	for(i in 1:num_pcs){

		pc_score=pca_res[["Scores"]][,i];
		cor_vect=apply(vals, 2, function(x){ cor(pc_score, x, method="spearman")});

		abs_sort_ix=order(abs(cor_vect), decreasing=T);
		sorted_cor=cor_vect[abs_sort_ix];		

		#print(head(sorted_cor));
		select_list[[i]]=sorted_cor[1:max_top];
		selected=c(selected, names(sorted_cor[1]));

	}
	selected=unique(selected);	


	prop_var=pca_res[["PCPropVar"]];

	if(!is.null(report_fn)){

		cat("Generating report...\n");
	
		report_matrix=matrix("", nrow=2+max_top, ncol=3*num_pcs)
		for(i in 1:num_pcs){

			report_matrix[,((i-1)*3)+1]=
				c(paste("PC", i, sep=""), "", names(select_list[[i]]));
			report_matrix[,((i-1)*3)+2]=
				c(sprintf("%6.4f", prop_var[i]), "", sprintf("%5.3f",select_list[[i]]));

		
			write.table(report_matrix, file=report_fn, quote=F, sep="\t", row.names=F, col.names=F);
		}
		print(report_matrix);
	}

	cat("Selected:\n");
	print(selected);
	return(selected);

}

###############################################################################

raw_FBA_matrix=load_community_fba(FBA_InFile, host="Vm");

nz_FBA_matrix=filter_low_flux_categories(raw_FBA_matrix, zero_cutoff=ZeroCutoff);

thresholded_flux_mat=threshold_flux_by_magnitude(nz_FBA_matrix, prop_to_keep=0.95);

flux_cor_mat=calculate_flux_correlation(thresholded_flux_mat);

#print(flux_cor_mat);
#pca_rec=run_pca(flux_cor_mat, thresholded_flux_mat, prop_to_keep=0.95);

# Run PCA
flux_cor_mat=project_PSD(flux_cor_mat);
pca_rec=run_pca(flux_cor_mat, thresholded_flux_mat, prop_to_keep=0.95);

select_pca_reps(pca_rec, thresholded_flux_mat, report_fn="selected.tsv");











##############################################################################

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
