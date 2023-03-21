#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
library('stats4'); # For mle

options(useFancyQuotes=F);
options(digits=5)

INDIV_PCA_CUTOFF=0.025;
CUMUL_PCA_CUTOFF=0.80;
BETA_CUTOFF=2.0;

params=c(
	"factors", "f", 1, "character",
	"outputroot", "o", 1, "character",
	"betas_fn", "b", 1, "character",

	"beta_cutoff", "x", 2, "numeric",
	"indiv_pca_cutoff", "i", 2, "numeric",
	"cumul_pca_cutoff", "c", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

#script_path=paste(head(strsplit(script_name, "/")[[1]], -1), collapse="/");
#source(paste(script_path, "/Impute_Matrix.r", sep=""));

NORM_PVAL_CUTOFF=.2;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <factors/metadata file name>\n",
	"	-b <beta values file name>\n",
	"	-o <output root file name>\n",
	"	[-x <Beta export cutoff for variable export>]\n",
	"	[-i <Individual PC Cutoff, default=", INDIV_PCA_CUTOFF, ">]\n",
	"	[-c <Cumulative PC Cutoff, default=", CUMUL_PCA_CUTOFF, ">]\n",
	"\n",
	"This script will read in the transformed variables that went into\n",
	"the beta estimation, and a list of the beta values and use a variety\n",
	"of methods to try to suggest a beta value cutoff.\n",
	"\n",
	"Methods will focus on rarefaction where variables are added to the accepted set\n",
	"one at a time in the order of decreasing beta.\n",
	"\n",
	"Beta values will be identified for:\n",
	"	1.) When the number of PCs exceeding the cutoff are maximized.\n",
	"	2.) When the number of PCs necessary to achieve the cumulative cutoff is minimized.\n",
	"	3.) When the entropy of the PCs is maximized.\n",
	"\n");

if(
	!length(opt$factors) || 
	!length(opt$betas_fn) || 
	!length(opt$outputroot)
){
	cat(usage);
	q(status=-1);
}

FactorsFname=opt$factors;
BetasFname=opt$betas_fn;
OutputFnameRoot=opt$outputroot;

Indiv_PC_Cutoff=INDIV_PCA_CUTOFF;
Cumul_PC_Cutoff=CUMUL_PCA_CUTOFF;
MinBeta_Cutoff=BETA_CUTOFF;

if(length(opt$indiv_pca_cutoff)){
	Indiv_PC_Cutoff=opt$indiv_pca_cutoff;
}

if(length(opt$cumul_pca_cutoff)){
	Cumul_PC_Cutoff=opt$cumul_pca_cutoff;
}

if(length(opt$beta_cutoff)){
	MinBeta_Cutoff=opt$beta_cutoff;
}

param_text=capture.output({
	cat("\n");
	cat("Transformed Factor File Name: ", FactorsFname, "\n");
	cat("Fitted betas list filename: ", BetasFname, "\n");
	cat("\n");
	cat("Individual PC Cutoff: ", Indiv_PC_Cutoff, "\n");
	cat("Cumulative PC Cutoff: ", Cumul_PC_Cutoff, "\n"); 
	cat("Min Beta Cutoff: ", MinBeta_Cutoff, "\n");
	cat("\n");
});

print(param_text, quote=F);

###############################################################################

load_factors=function(fname){
	factors=as.data.frame(read.delim(fname,  header=TRUE, row.names=1, check.names=FALSE, sep="\t"));
	return(factors);
}

load_list=function(fname){
	cat("Loading: ", fname, "\n");
	lst=read.delim(fname, header=F, check.names=F, comment.char="#", as.is=T);
	return(lst[,1]);	
}

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

		if(heatmap_height < 200 &&  heatmap_width < 200){
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

##############################################################################
# Main Program Starts Here!
##############################################################################

pdf(paste(OutputFnameRoot, ".beta_analyses.pdf", sep=""), height=11, width=8.5);
plot_text(param_text);

# Load factors
cat("Loading Factors...\n");
loaded_factors=load_factors(FactorsFname);
loaded_factor_names=colnames(loaded_factors);
loaded_sample_names=rownames(loaded_factors);

num_samples=nrow(loaded_factors);
cat("Number of Samples:\n");

cat("Loaded factors:\n");
print(loaded_factor_names);
cat("\n");
cat("Loaded sample ids:\n");
print(loaded_sample_names);
cat("\n");

# Load betas
cat("Loading Betas...\n");
loaded_betas=load_factors(BetasFname);

cat("Loaded Betas:\n");
print(loaded_betas);
cat("\n");

##############################################################################

beta_order=order(loaded_betas, decreasing=T);

sorted_beta=loaded_betas[beta_order,,drop=F];
sorted_factors=loaded_factors[,beta_order];

print(sorted_beta);
log_sorted_beta=log10(sorted_beta[,1]);
names(log_sorted_beta)=rownames(sorted_beta);

##############################################################################

rarefy_for_pca=function(sorted_factors, min_var, max_var){

	# Precalculate the correlation matrix
	cat("Calculating Correlation Matrix...\n");
	cor_mat=cor(sorted_factors[,1:max_var]);

	pc_matrix=matrix(0, nrow=max_var, ncol=max_var);

	cat("Stepping through rarefaction...\n");
	for(vrang in min_var:max_var){
		cat(".");

		sub_cor_mat=cor_mat[1:vrang, 1:vrang];

		eigen=eigen(sub_cor_mat);
		pca_propvar=eigen$values/sum(eigen$values);
		pc_matrix[vrang, 1:vrang]=pca_propvar;
	}
	cat("\n");

	return(pc_matrix);

}

#------------------------------------------------------------------------------

# Calculate exceeding cumulative
calculate_num_PCs_exceeding_cum_cutoff=function(pc_rare_res, cumul_cutoff){
	num_pcs_exc=apply(pc_rare_results, 1, function(x){
		sum(cumsum(x)<cumul_cutoff)+1;
		});	
	return(num_pcs_exc);
}

# Calculate exceeding individual
calculate_num_PCs_exceeding_ind_cutoff=function(pc_rare_res, indiv_cutoff){
	num_pcs_exc=apply(pc_rare_res, 1, function(x){
		sum(x>indiv_cutoff);
		});

	return(num_pcs_exc);
}

calculate_coverage_for_num_PCs_exceeding_cutoff=function(pc_rare_res, numPCs_exc_cutoff){
	num_var=length(numPCs_exc_cutoff);
	coverages=numeric(num_var);
	for(i in 1:num_var){
		coverages[i]=sum(pc_rare_res[i, 1:numPCs_exc_cutoff[i]]);
	}	
	return(coverages);
}

find_min_var_max_pcs=function(num_pcs_exc_ind_cutoff, min_num_variables){
	avail_pcs_range=num_pcs_exc_ind_cutoff;
	avail_pcs_range[1:(min_num_variables-1)]=0;
	max_pcs=max(avail_pcs_range);
	min_var_wmax_pcs=min(which(max_pcs==avail_pcs_range));
	return(min_var_wmax_pcs);
}


#------------------------------------------------------------------------------

max_num_var=max(which(sorted_beta>MinBeta_Cutoff));
min_num_var=min(max_num_var, num_samples);

cat("Range for Number of Variables to consider: ", min_num_var, " - ", max_num_var, "\n");

#max_num_var=600;

pc_rare_results=rarefy_for_pca(loaded_factors, 1, max_num_var);

plot_rare_curves=function(pc_rare_results, start, end, pc_indices, 
	min_nvar_threshold,
	recommended_num_var,
	betas_sorted,
	n_pcs_exc_ind, n_pcs_exc_cum, cov_pcs_exc_ind){

	end=min(end,nrow(pc_rare_results));

	var_ix=start:end;

	par(mfrow=c(3,1));
	par(mar=c(4,4.2,6,4.2));

	#-----------------------------------------------------------------------------
	plot(0, type="n", xlim=c(start, end), ylim=c(0,1.1),
		main=paste("Rarefaction: Range [", start, " - ", end, "]", sep=""),
		ylab="Proportion of Variance", xlab="Num Variables Included"	
		);

	plot_par_usr_x=par()$usr[c(1,2)];
	axis_ticks=c(1, setdiff(axisTicks(plot_par_usr_x, log=F), 0));	
	beta_labels=sprintf("%3.2f", betas_sorted[axis_ticks, "beta"]);
	beta_at_rec=sprintf("%3.2f", betas_sorted[recommended_num_var, "beta"]);

	axis(side=3, at=axis_ticks, labels=beta_labels, 
		padj=.65, cex.axis=.8,
		col="purple", col.tick="purple", col.axis="purple");

	abline(v=recommended_num_var, col="orange", lwd=3);
	abline(v=min_nvar_threshold, col="black", lwd=1.5);
	mtext(paste(" rec=", recommended_num_var, sep=""), 
		side=3, line=-1.5, outer=F, at=recommended_num_var, adj=0, col="orange", cex=.8);
	mtext(paste(" beta=", beta_at_rec, sep=""), 
		side=3, line=-2.5, outer=F, at=recommended_num_var, adj=0, col="purple", cex=.8);
	mtext("min ", side=3, line=-1.5, outer=F, at=min_nvar_threshold, adj=1, col="black", cex=.8);

	# Plot PC coverages for specified PCs
	for(pcix in pc_indices){
		points(var_ix, pc_rare_results[var_ix, pcix], col=pcix, cex=.5);
	}

	#-----------------------------------------------------------------------------
	# Plot num PCs for individual coverage
	max_exceeding=max(n_pcs_exc_ind);
	plot(var_ix, n_pcs_exc_ind[var_ix], xlim=c(start, end), ylim=c(0, max_exceeding+2),
		ylab="Number of PCs", xlab="Num Variables Included",
		main="Number of PCs exceeding Individual Cutoff",
		cex=.5);

	axis(side=3, at=axis_ticks, labels=beta_labels, 
		padj=.65, cex.axis=.8,
		col="purple", col.tick="purple", col.axis="purple");

	abline(v=recommended_num_var, col="orange", lwd=3);
	abline(v=min_nvar_threshold, col="black", lwd=1.5);
	mtext(paste(" rec=", recommended_num_var, sep=""), 
		side=3, line=-1.5, outer=F, at=recommended_num_var, adj=0, col="orange", cex=.8);
	mtext(paste(" beta=", beta_at_rec, sep=""), 
		side=3, line=-2.5, outer=F, at=recommended_num_var, adj=0, col="purple", cex=.8);
	mtext("min ", side=3, line=-1.5, outer=F, at=min_nvar_threshold, adj=1, col="black", cex=.8);

	# Overlay with coverage
	points(cov_pcs_exc_ind[var_ix]*max_exceeding, col="blue", type="l", lwd=2);
	markers=c(0,.25,.50,.75, 1);
	mtext("Cumulative Variance of PCs Exceeding Cutoff", side=4, col="blue", cex= .7, line=3);
	axis(side=4, at=markers*max_exceeding, labels=markers, col="blue", col.axis="blue");


	#-----------------------------------------------------------------------------
	# Plot num PCs for cumulative coverage
	max_exceeding=max(n_pcs_exc_cum);
	plot(var_ix, n_pcs_exc_cum[var_ix], xlim=c(start, end), ylim=c(0, max_exceeding+2),
		ylab="Number of PCs", xlab="Num Variables Included",
		main="Number of PCs to Exceed Cumulative Cutoff",
		cex=.5);

	axis(side=3, at=axis_ticks, labels=beta_labels, 
		padj=.65, cex.axis=.8,
		col="purple", col.tick="purple", col.axis="purple");

	abline(v=recommended_num_var, col="orange", lwd=3);
	abline(v=min_nvar_threshold, col="black", lwd=1.5);
	mtext(paste(" rec=", recommended_num_var, sep=""), 
		side=3, line=-1.5, outer=F, at=recommended_num_var, adj=0, col="orange", cex=.8);
	mtext(paste(" beta=", beta_at_rec, sep=""), 
		side=3, line=-2.5, outer=F, at=recommended_num_var, adj=0, col="purple", cex=.8);
	mtext("min ", side=3, line=-1.5, outer=F, at=min_nvar_threshold, adj=1, col="black", cex=.8);

}

#------------------------------------------------------------------------------

min_num_var=num_samples;

num_pcs_exc_indv_cutoff=calculate_num_PCs_exceeding_ind_cutoff(pc_rare_results, Indiv_PC_Cutoff);
#print(num_pcs_exc_indv_cutoff);

num_pcs_exc_cuml_cutoff=calculate_num_PCs_exceeding_cum_cutoff(pc_rare_results, Cumul_PC_Cutoff);
#print(num_pcs_exc_cuml_cutoff);

cov_of_pcs_exc_indv_cutoff=calculate_coverage_for_num_PCs_exceeding_cutoff(pc_rare_results, num_pcs_exc_indv_cutoff);
#print(cov_of_pcs_exc_indv_cutoff);

recommended_num_var=find_min_var_max_pcs(num_pcs_exc_indv_cutoff, min_num_var);
cat("Recommended Number of Variables: ", recommended_num_var, "\n");

beta_at_rec_num=sorted_beta[recommended_num_var, "betas"];

#------------------------------------------------------------------------------

hist(log_sorted_beta, main="Distribution of Log10(Beta)", xlab="Log10(Beta)", breaks=30, col="purple");
abline(v=log10(2), col="blue", lty="dashed");
abline(v=log_sorted_beta[recommended_num_var], col="orange", lty="dashed");


plot_rare_curves(pc_rare_results, start=1, end=2*min_num_var, pc_indices=1:10, 
	min_num_var,
	recommended_num_var,
	sorted_beta,
	num_pcs_exc_indv_cutoff, num_pcs_exc_cuml_cutoff, cov_of_pcs_exc_indv_cutoff);
	
plot_rare_curves(pc_rare_results, start=1, end=2*recommended_num_var, pc_indices=1:10, 
	min_num_var,
	recommended_num_var,
	sorted_beta,
	num_pcs_exc_indv_cutoff, num_pcs_exc_cuml_cutoff, cov_of_pcs_exc_indv_cutoff);

plot_rare_curves(pc_rare_results, start=1, end=Inf, pc_indices=1:10,
	min_num_var,
	recommended_num_var,
	sorted_beta,
	num_pcs_exc_indv_cutoff, num_pcs_exc_cuml_cutoff, cov_of_pcs_exc_indv_cutoff);

dev.off();

#------------------------------------------------------------------------------
# Export variables

sort_var_names=rownames(sorted_beta);

fname=paste(OutputFnameRoot, ".prioritized.tsv", sep="");
fh=file(fname, "w");
cat(file=fh, "SampleID\t");
close(fh);
write.table(loaded_factors[, sort_var_names[1:recommended_num_var]], 
	file=fname, quote=F, row.names=T, sep="\t", append=T);


##############################################################################

cat("Done.\n");

print(warnings());
q(status=0);
