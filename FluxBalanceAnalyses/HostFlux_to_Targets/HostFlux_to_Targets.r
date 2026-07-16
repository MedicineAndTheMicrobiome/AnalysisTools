#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

source("~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r");
source("~/git/AnalysisTools/Metadata/OutputFileLibrary/OutputFileLibrary.r");

options(useFancyQuotes=F);
options(digits=5)

params=c(
	"fba_infile", "f", 1, "character",
	"outputroot", "o", 1, "character",
	"zero_cutoff", "z", 2, "numeric",
	"prop_flux_mag", "t", 2, "numeric",
	"pca_cutoff", "p", 2, "numeric"
);

ZERO_TOLERANCE=1e-6;
PROP_FLUX_MAG=.99;
PROP_PCA_COV=.95;


opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

script_path=paste(head(strsplit(script_name, "/")[[1]], -1), collapse="/");

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <Combined FBA File (single line per sample)>\n",
	"	-o <output filename root>\n",
	"\n",
	"	[-z <cutoff for 0, default=", ZERO_TOLERANCE, ">]\n",
	"	[-t <cutoff for top |flux| to keep, default=", PROP_FLUX_MAG, ">]\n",
	"	[-p <cutoff for cumulative PCA coverage, default=", PROP_PCA_COV, ">]\n", 
	"\n",
	"This script will read in the FBA results, where each line\n",
	"represents a single sample, and the columns contain\n",
	"flux values.  Each sample/row represents what the host in \n",
	"in a multi-organism community FBA model is experiencing.\n",
	"\n",
	"This script will do the following:\n",
	"	1.) Provide summary stats across the fluxes\n",
	"	2.) Set fluxes below zero_cutoff to 0\n",
	"	3.) Remove flux columns that are 0's across all samples\n",
	"	4.) Compute correlation\n",
	"	5.) Run PCA\n",
	"	6.) Identify top PCA proxies.\n",
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
PropFluxMag=PROP_FLUX_MAG;
PropPCACov=PROP_PCA_COV;

if(length(opt$zero_cutoff)){
	ZeroCutoff=opt$zero_cutoff;
}

if(length(opt$prop_flux_mag)){
	PropFluxMag=opt$prop_flux_mag;
}

if(length(opt$pca_cutoff)){
	PropPCACov=opt$pca_cutoff;
}

##############################################################################

pdf(paste(OutputFnameRoot, ".fba.pca.pdf", sep=""), height=11, width=8.5);
options(width=180);

param_text=capture.output({
	cat(script_name);
	cat("\n\n");
	cat("Flux File by Sample: ", FBA_InFile, "\n");
	cat("Output File Name Root: ", OutputFnameRoot, "\n");
	cat("Zero-Cutoff: ", ZeroCutoff, "\n");
	cat("Flux Coverage: ", PropFluxMag, "\n");
	cat("PCA Coverage: ", PropPCACov, "\n");
	cat("\n");
});
print(param_text, quote=F);
plot_text(param_text);

##############################################################################

load_fba_file=function(fname){

	raw_mat=as.data.frame(read.table(fname,  header=TRUE, check.names=FALSE, sep="\t"));

	rownames(raw_mat)=raw_mat[,1];
	raw_mat=raw_mat[,2:ncol(raw_mat)];

	headers=colnames(raw_mat);

	cat("\n");
	cat("First headers:\n");
	print(head(headers));
	cat("...\n");
	cat("\n");

	num_samples=nrow(raw_mat);
	num_fluxes=ncol(raw_mat);

	cat("------------------------------------------------------------------\n");
	cat("Input Matrix Dimensions:\n");
	cat("Num Samples: ", num_samples, "\n");
	cat("Num Fluxes: ", num_fluxes, "\n");
	cat("------------------------------------------------------------------\n");

	return(raw_mat);
}

#------------------------------------------------------------------------------

filter_low_flux_categories=function(raw_fb_mat, zero_cutoff=1e-6){

	cat("Screening flux categories at cutoff: ", zero_cutoff, "\n");
	zeroed_mat=apply(raw_fb_mat, 1:2, function(x){ ifelse(abs(x)<zero_cutoff, 0, x)});
	zero_cols=apply(zeroed_mat, 2, function(x){ sum(abs(x))==0});

	num_zero_cols=sum(zero_cols);
	cat("Num All Zero columns removed: ", num_zero_cols, "\n");
	clean_mat=zeroed_mat[,!zero_cols];
	return(clean_mat);
}

#------------------------------------------------------------------------------

select_top_fluxes=function(fb_mat, target_flux_coverage=.95){

	sum_flux_mag=apply(fb_mat, 2, function(x){ sum(abs(x))});
	sorted_sfm=sort(sum_flux_mag, decreasing=T);
	sorted_flux_names=names(sorted_sfm);
	print(sorted_sfm);

	total_flux=sum(sorted_sfm);
	normalize_flux=sorted_sfm/total_flux;
	cmsm_nf=cumsum(normalize_flux);	

	num_top_flux_cov=min(which(cmsm_nf>target_flux_coverage));
	cat("Num Fluxes to acquire targeted proportion of flux coverage: ", num_top_flux_cov, "\n");

	par(mfrow=c(1,1));
	par(mar=c(4,4,4,4));
	barplot(cmsm_nf[1:num_top_flux_cov], ylim=c(0,1), las=2, 
		names.arg=1:num_top_flux_cov,
		xlab="Num Variables", ylab="Proportion");

	title(main=paste("Proportion of sum(|flux|) covered at ", target_flux_coverage));
	title(main=paste("Num fluxes kept: ", num_top_flux_cov), line=.5, font.main=3); 

	abline(h=target_flux_coverage, col="red");

	kept_fb_mat=fb_mat[,sorted_flux_names[1:num_top_flux_cov], drop=F];
	return(kept_fb_mat);
}

#------------------------------------------------------------------------------

transform_fluxes=function(fb_mat){

	par(mfrow=c(5,2));

	fb_names=colnames(fb_mat);
	trans_mat=fb_mat;

	plot_title_page(
		title="Inverse Hyperbolic\nSine Transform",
		subtitle="\n\n\n\nasinh(x) = log(x + sqrt(x^2 + 1)"
		);

	for(nm in fb_names){

		vals=fb_mat[,nm];

		hist(vals, breaks=20, main=paste("Original: ", nm));
		abline(v=0, col="red", lwd=2);

		# Apply Inverse Hyberbolic Sine
		tran=log(vals + sqrt(vals^2 + 1));
		trans_mat[,nm]=tran;

		hist(tran, breaks=20, main=paste("Transformed: ", nm));
		abline(v=0, col="red", lwd=2);

	}
	
	return(trans_mat);

}

#------------------------------------------------------------------------------

perform_pca=function(tar_mat, pc_cumthres=0.95){

	cor_mat=cor(tar_mat);
	num_input_var=ncol(tar_mat);
	#print(cor_mat);

	eigen_rec=eigen(cor_mat);
	pca_propvar=eigen_rec$values/sum(eigen_rec$values);
	pca_propcumsum=cumsum(pca_propvar);
	num_pcs_at_cutoff=sum(pca_propcumsum<pc_cumthres)+1;
	scores=(scale(tar_mat, center=T, scale=T) %*% eigen_rec$vectors);

	cat("Cumulative PCA Variance: \n");
	print(pca_propcumsum);
	
	cat("Num PCs at Cutoff: ", num_pcs_at_cutoff, "\n");

	plot_title_page(
		title="PCA",
		);


	par(mfrow=c(2,1));
	barplot(pca_propcumsum[1:num_pcs_at_cutoff], xlab="Num PCs", ylab="Prop Variance",
		names.arg=1:num_pcs_at_cutoff, cex.names=.9, las=2);
	abline(h=pc_cumthres, col="red", lty="dashed");
	title(main=paste(
		num_pcs_at_cutoff, " of ", num_input_var," PCs required to cover ", pc_cumthres, " of variance", sep=""
		));


	barplot(pca_propvar[1:num_pcs_at_cutoff], xlab="Num PCs", ylab="Variance Contrib",
		names.arg=1:num_pcs_at_cutoff, cex.names=.9, las=2);

	results=list();
	results[["prop_var"]]=pca_propvar;
	results[["prop_cumcum"]]=pca_propcumsum;
	results[["num_pcs_at_cutoff"]]=num_pcs_at_cutoff;
	results[["scores"]]=scores;

	return(results);
}

#------------------------------------------------------------------------------

generate_ordination_plots=function(pca_results, pc_proxies=NULL){

	scores=pca_results[["scores"]];

	par(mfrow=c(4,3));
	par(mar=c(4,4,3,4));
	num_samples=nrow(scores);
	
	all_extremes=c();

	num_pcs_at_cutoff=pca_results[["num_pcs_at_cutoff"]];

	if(is.null(pc_proxies)){
		pc_proxies=rep("", num_pcs_at_cutoff);
	}else{
		pc_proxies=paste("[", pc_proxies, "]", sep="");
	}

	plot_title_page(
		title="PCA Score Ordination Plots",
		);

	for(i in 1:(num_pcs_at_cutoff%/%2)){
		ix=((i-1)*2)+1;

		plot(scores[,ix], scores[,ix+1], 
			xlab=paste("PC: ", ix, " ", pc_proxies[ix]),
			ylab=paste("PC: ", ix+1, " ", pc_proxies[ix+1]),
			col="grey"
			);

		topbot=c();
		
		pcX=scores[,ix];
		names(pcX)=1:num_samples;
		pcXsort=sort(pcX);
		topbot=c(topbot, names(pcXsort[c(1,2,num_samples-1, num_samples)]));

		pcY=scores[,ix+1];
		names(pcY)=1:num_samples;
		pcYsort=sort(pcY);
		topbot=c(topbot, names(pcYsort[c(1,2,num_samples-1, num_samples)]));

		topbot_ix=unique(as.numeric(topbot));
		topbot_labels=topbot_ix
		text(scores[topbot_ix,ix], scores[topbot_ix,ix+1],
			col="red", cex=.5, labels=topbot_labels,
			);

		all_extremes=c(all_extremes, topbot);
	
	}

	# Accumulate table of extreme/outlier samples across all PCs

	samp_ids=rownames(scores);
	tab=table(all_extremes);
	tab_num_ids=as.numeric(names(tab));

	#print(tab);
	#print(tab_ids);

	# Sort by numeric id 
	num_sort_ix=order(tab_num_ids);
	tab=tab[num_sort_ix];
	tab_num_ids=tab_num_ids[num_sort_ix];

	extr_samp_ids=samp_ids[tab_num_ids];

	infomat=cbind(extr_samp_ids, tab_num_ids, tab);
	colnames(infomat)=c("SampleID", "NumericID", "Frequency");

	plot_title_page(
		title="Accumulated Extreme\nPC Samples",
		);

	# Output by Numeric ID
	rownames(infomat)=NULL;
	plot_text(c(
		"By Index ID:",
		"",
		capture.output(print(as.data.frame(infomat), quote=F))));
	print(infomat);

	# Output by Frequency
	by_freq_ix=order(infomat[,"Frequency"], decreasing=T, method="shell");
	plot_text(c(
		"By Decreasing Frequencies:",
		"",
		capture.output(print(as.data.frame(infomat[by_freq_ix,]), quote=F))));
	print(infomat[by_freq_ix,]);

}

#------------------------------------------------------------------------------

select_pca_proxies=function(pca_res, tar_mat, max_top_proxy_to_report=10){

	cat("Selecting Proxies for each PC:\n");
	
	# Num PCs to find proxies for
	num_pcs=pca_res[["num_pcs_at_cutoff"]];
	num_flux_var=ncol(tar_mat);
	max_tptr=min(max_top_proxy_to_report, num_flux_var);

	cat("Num PCs Kept: ", num_pcs, "\n");
	cat("Num Flux Variables: ", num_flux_var, "\n");
	cat("Max Top Proxies to Report: ", max_tptr, "\n");

	select_list=list();
	selected=c();
	
	scores=pca_res[["scores"]];

	for(i in 1:num_pcs){

		# For each PC, compute correlation with each flux.
		pc_score=scores[,i];
		cor_vect=apply(tar_mat, 2, function(x){ cor(pc_score, x)});

		# Order the correlatons by greatest magnitude
		abs_sort_ix=order(abs(cor_vect), decreasing=T);
		sorted_cor=cor_vect[abs_sort_ix];		

		#print(head(sorted_cor));
		select_list[[i]]=sorted_cor[1:max_tptr];
		selected=c(selected, names(sorted_cor[1]));

	}

	results=list();
	results[["selected"]]=selected;
	results[["unique_proxies"]]=unique(selected);
	results[["groups"]]=select_list;

	return(results);
}

#------------------------------------------------------------------------------

barplot_top_proxies=function(sel_prox){

	par(mfrow=c(4,3));
	par(mar=c(7,4,3,4));

	sel_prox_list=sel_prox[["groups"]];
	num_pc=length(sel_prox_list);

	plot_title_page(
		title="Top PCs and Closest Proxies",
		);

	cat("Num PCs to Plot: ", num_pc, "\n");

	for(i in 1:num_pc){

		top_proxy=names(sel_prox_list[[i]])[1];


		barplot(sel_prox_list[[i]],
			ylim=c(-1,1), las=2,
			main=paste("PC: ", i)
			);

		title(main=top_proxy, line=.25, font.main=2, cex.main=.7);

	}
	
}

#------------------------------------------------------------------------------

export_proxies=function(targets_arr, values_mat, outfn){

	cat("Num Proxies: ", length(targets_arr), "\n");
	uniq_targ=unique(targets_arr);
	num_uniq_targ=length(uniq_targ);
	cat("Num Unique Proxies: ", num_uniq_targ, "\n");

	cat("Selected Proxies: \n");
	print(uniq_targ);

	selected_mat=values_mat[,uniq_targ,drop=F];

	outmat=cbind(rownames(selected_mat), selected_mat);
	cn=c("SampleID", colnames(selected_mat));
	colnames(outmat)=cn;
	write.table(outmat, file=outfn, quote=F, sep="\t", row.names=F);

}

#------------------------------------------------------------------------------

plot_selected=function(targets_arr, values_mat){

	cat("Num Proxies: ", length(targets_arr), "\n");
	uniq_targ=unique(targets_arr);
	num_uniq_targ=length(uniq_targ);
	cat("Num Unique Proxies: ", num_uniq_targ, "\n");

	cat("Selected Proxies: \n");
	print(uniq_targ);

	par(mfrow=c(4,3));
	for(nm in uniq_targ){
		hist(values_mat[,nm], main=nm, xlab="Values");
	}


}


#------------------------------------------------------------------------------

###############################################################################

# Load file
raw_inmat=load_fba_file(FBA_InFile);
num_samples=nrow(raw_inmat);
num_orig_fluxes=ncol(raw_inmat);

#------------------------------------------------------------------------------

# Remove 0 and very low abundance categories
nolowflux_mat=filter_low_flux_categories(raw_inmat, ZeroCutoff);
num_low_flux_filtered=ncol(nolowflux_mat);

#------------------------------------------------------------------------------

# Keep the top fluxes 
top_fluxes_mat=select_top_fluxes(nolowflux_mat, PropFluxMag);
num_top_fluxes_kept=ncol(top_fluxes_mat);

#------------------------------------------------------------------------------

# Transform the fluxes
trans_top_flux_mat=transform_fluxes(top_fluxes_mat);

#------------------------------------------------------------------------------

# Perform PCA
pca_res=perform_pca(trans_top_flux_mat, PropPCACov);
num_PCs_kept=pca_res[["num_pcs_at_cutoff"]];

#------------------------------------------------------------------------------

# Select the top proxies
selected_proxies=select_pca_proxies(pca_res, top_fluxes_mat);
num_proxies=length(selected_proxies[["unique_proxies"]]);

#------------------------------------------------------------------------------

# Plot ordination plots with proxy names labeled
generate_ordination_plots(pca_res, pc_proxies=selected_proxies[["selected"]]);

#------------------------------------------------------------------------------

# Plot selected proxies and other top close fluxes
barplot_top_proxies(selected_proxies);

#------------------------------------------------------------------------------

proxylist=cbind(selected_proxies[["unique_proxies"]]);
colnames(proxylist)="Unique Selected Proxies";
rownames(proxylist)=1:num_proxies;

plot_text(c(
	paste("FBA In File: ", FBA_InFile),
	"",
	"Original Input Dimensions:",
	paste("Num Samples: ", num_samples),
	paste("Num Fluxes: ", num_orig_fluxes),
	"",
	paste("Min flux trim (|f|<", ZeroCutoff, "):"),
	paste("Num Fluxes: ", num_low_flux_filtered),
	"",
	paste("Min prop of flux to keep (sum(|f|)>", PropFluxMag, "):"),
	paste("Num Fluxes: ", num_top_fluxes_kept),
	"",
	paste("Min PC coverage (sum(var)>", PropPCACov, ")"),
	paste("Num PCs selected: ", num_PCs_kept),
	"",
	paste("Num Unique Proxies Selected: ", num_proxies),
	"",
	capture.output(print(proxylist, quote=F))
));
#---------------------------------------------

plot_selected(
	targets_arr=selected_proxies[["unique_proxies"]], 
	values_mat=trans_top_flux_mat);

export_proxies(
	targets_arr=selected_proxies[["unique_proxies"]], 
	values_mat=trans_top_flux_mat,
	outfn=paste(OutputFnameRoot, ".trans.proxies.tsv", sep=""));


##############################################################################

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
