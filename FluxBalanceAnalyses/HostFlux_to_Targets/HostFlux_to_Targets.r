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
	"zero_cutoff", "z", 2, "numeric"
);

ZERO_TOLERANCE=1e-6;

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

script_path=paste(head(strsplit(script_name, "/")[[1]], -1), collapse="/");

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <Combined FBA File (single line per sample)>\n",
	"	-o <output filename root>\n",
	"\n",
	"	[-z <cutoff for 0, default=", ZERO_TOLERANCE, ">]\n",
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

if(length(opt$zero_cutoff)){
	ZeroCutoff=opt$zero_cutoff;
}

param_text=capture.output({
	cat("\n");
	cat("Flux File by Sample: ", FBA_InFile, "\n");
	cat("Output File Name Root: ", OutputFnameRoot, "\n");
	cat("Zero-Cutoff: ", ZeroCutoff, "\n");
	cat("\n");
});

print(param_text, quote=F);

options(width=180);


##############################################################################
# Main Program Starts Here!
##############################################################################

pdf(paste(OutputFnameRoot, ".fba.pca.pdf", sep=""), height=11, width=8.5);

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


filter_low_flux_categories=function(raw_fb_mat, zero_cutoff=1e-6){

	cat("Screening flux categories at cutoff: ", zero_cutoff, "\n");
	zeroed_mat=apply(raw_fb_mat, 1:2, function(x){ ifelse(x<zero_cutoff, 0, x)});
	zero_cols=apply(zeroed_mat, 2, function(x){ sum(abs(x))==0});

	num_zero_cols=sum(zero_cols);
	cat("Num All Zero columns removed: ", num_zero_cols, "\n");
	clean_mat=zeroed_mat[,!zero_cols];
	return(clean_mat);
}

perform_pca=function(tar_mat, pc_cumthres=0.95){

	cor_mat=cor(tar_mat);
	#print(cor_mat);

	eigen_rec=eigen(cor_mat);
	pca_propvar=eigen_rec$values/sum(eigen_rec$values);
	pca_propcumsum=cumsum(pca_propvar);
	num_pcs_at_cutoff=sum(pca_propcumsum<pc_cumthres)+1;
	scores=(scale(tar_mat, center=T, scale=T) %*% eigen_rec$vectors);

	cat("Cumulative PCA Variance: \n");
	print(pca_propcumsum);
	
	cat("Num PCs at Cutoff: ", num_pcs_at_cutoff, "\n");

	par(mfrow=c(2,1));
	barplot(pca_propcumsum[1:(num_pcs_at_cutoff+1)], xlab="Num PCs", ylab="Prop Variance",
		names.arg=1:(num_pcs_at_cutoff+1), cex.names=.9, las=2);
	barplot(pca_propvar[1:(num_pcs_at_cutoff+1)], xlab="Num PCs", ylab="Variance Contrib",
		names.arg=1:(num_pcs_at_cutoff+1), cex.names=.9, las=2);

	results=list();
	results[["prop_var"]]=pca_propvar;
	results[["prop_cumcum"]]=pca_propcumsum;
	results[["num_pcs_at_cutoff"]]=num_pcs_at_cutoff;
	results[["scores"]]=scores;

	return(results);
}

generate_ordination_plots=function(pca_results, pc_proxies){

	scores=pca_results[["scores"]];

	par(mfrow=c(4,3));
	par(mar=c(4,4,3,4));
	num_samples=nrow(scores);
	
	all_extremes=c();

	num_pcs_at_cutoff=pca_results[["num_pcs_at_cutoff"]];

	for(i in 1:(num_pcs_at_cutoff%/%2)){
		ix=((i-1)*2)+1;

		plot(scores[,ix], scores[,ix+1], 
			xlab=paste("PC: ", ix),
			ylab=paste("PC: ", ix+1),
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

###############################################################################

raw_inmat=load_fba_file(FBA_InFile);
nolowflux_mat=filter_low_flux_categories(raw_inmat, ZeroCutoff);

pca_res=perform_pca(nolowflux_mat);

generate_ordination_plots(pca_res, NULL);


#print(nolowflux_mat);
quit();

	#print(out_mat);

	return(out_mat);

}


#------------------------------------------------------------------------------

select_pca_reps=function(pca_res, raw_fb_mat, max_top=15, report_fn=NULL){
	
	num_pcs=pca_res[["NumPCsKept"]];

	vals=raw_fb_mat[,3:ncol(raw_fb_mat)];


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

#------------------------------------------------------------------------------






##############################################################################

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
