#!/usr/bin/env Rscript

library('getopt');

source("~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r");
source("~/git/AnalysisTools/Metadata/OutputFileLibrary/OutputFileLibrary.r");

###############################################################################

params=c(
	"summary_table", "s", 1, "character",
	"reference_si_fn", "r", 1, "character",
	"output_root", "o", 1, "character",
	"proport_fn", "p", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-s <input summary_table.tsv>\n",
	"	-r <reference spike-in table filename>\n",
	"	-o <output file root>\n",
	"	[-p <proportions tsv filename>]\n",
	"\n",
	"Example Format of Reference Spike In File:\n",
	"\n",
	"First column (ReferenceName):\n",
	"   Bacteria;Deinococcota;Deinococci;Deinococcales;Trueperaceae;Truepera\n",
	"   Bacteria;Bacteroidota;Bacteroidia;Flavobacteriales;Flavobacteriaceae;Imtechella\n",
	"   Bacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillaceae_uncl\n",
	"\n",
	"Second column (CopyCount):\n",	
	"	2e5\n",
	"		(2.0 x 10^5 copies per 20uL)\n",
	"	3e4\n",
	"		(3.0 x 10^4 copies per 20uL)\n",
	"	7e3\n",
	"		(7.0 x 10^3 copies per 20uL)\n",
	"\n",
	"Third column (MinLog10Abund):\n",
	"	-1.0\n",
	"		(10^-1 = .10 = 10% abundance)\n",
	"	-1.25\n",
	"		(10^-1.25 = .0562 = 5.62% abundance)\n",
	"	-2.75\n",
	"		(10^-2.75 = .00178 = 0.178% abundance)\n",
	"\n",
	"Note that the concentration that was spiked in is not important, rather it's the\n",
	"number of 16S (or reference DNA) copies that was mixed into the experimental sample.\n",
	"In the case of a the Zymo protocol, 20uL of reagent was spiked in, so the copies\n",
	"is based on cell density (cells/20 uL) x 16S copies per cell.\n",
	"\n",
	"The volume (or concentration) of the recipient experimental sample is also not\n",
	"important, as long as the sampling unit, e.g. swap, scoop, drop, of the experimental\n",
	"sample is consistent.\n",
	"\n",
	"The third column contains the minimum relative abundance of the reference needed\n",
	"to be identified in each sample in order for the sample's copy counts to be\n",
	"considered reliable.\n",
	"\n",
	"Outputs:\n",
	"	1.) Summary Tables: *.[N].summary_table.tsv\n",
	"	For each of the reference categories and combined, a summary table of the\n",
	"	copy counts per category will be estimated.\n",
	"\n",
	"	2.) Copy Count Statistics: *.copy_count_stats.tsv\n",
	"	This file reports for each sample and each reference (and combined)\n",
	"	statistics including, sample read depth, reference abundance, estimated\n",
	"	total copy counts for each sample, and whether the reference abundance\n",
	"	is above the required minimum (reliable)\n",
	"\n",
	"	3.) Spike-in Analysis PDF: *.abs_spikein.pdf\n",
	"	This PDF contains histograms and figures which try to illustrate the\n",
	"	the reliability and relationships between the reference, sample depth,\n",
	"	etc., for each reference category.\n",	
	"\n",
	"\n");

if(
	!length(opt$summary_table)||
	!length(opt$reference_si_fn)||
	!length(opt$output_root)
	){
	cat(usage);
	q(status=-1);
}

options(width=200);

###############################################################################

InputSummaryTable=opt$summary_table;
ReferenceFilename=opt$reference_si_fn;
OutputRoot=opt$output_root;

ProportionsFile=NULL;
if(length(opt$proport_fn)){
	ProportionsFile=opt$proport_fn;
}

cat("\n");
cat("Input Summary Table: ", InputSummaryTable, "\n");
cat("Reference Category: ", ReferenceFilename, "\n");
cat("Output Root: ", OutputRoot, "\n");
cat("Proportions File: ", ProportionsFile, "\n");
cat("\n");

pdf(paste(OutputRoot, ".abs_spikein.pdf", sep=""), height=8.5, width=11);

###############################################################################

plot_text(c(
	script_name,
	"",
	paste("Input Summary Table: ", InputSummaryTable),
	paste("Reference Category: ", ReferenceFilename),
	paste("Output Root: ", OutputRoot),
	paste("Proportions File: ", ProportionsFile)
	));

###############################################################################
# Load summary table

orig_counts=load_summary_file(InputSummaryTable);
sample_depths=apply(orig_counts, 1, sum);
log_sample_depths=log10(sample_depths);

cat("Normalizing...\n");
normalized=normalize(orig_counts);

categories=colnames(normalized);
#print(categories);

###############################################################################
# Load reference copy count table

reference_info=read.table(ReferenceFilename, header=T, sep="");

cat("Loaded Reference Counts:\n");
print(reference_info);
cat("\n");

reference_names=reference_info[,"ReferenceName"];
reference_copycnts=reference_info[,"CopyCount"];
names(reference_copycnts)=reference_names;

reference_log10minabund=reference_info[,"MinLog10Abund"];
names(reference_log10minabund)=reference_names;

combined_ref_copycount=sum(reference_info[,"CopyCount"]);
combined_min_log10abund=log10(sum(10^reference_log10minabund));

cat("Combined Ref CopyCount: ", combined_ref_copycount, "\n");
cat("Combined Min Log10(Abund): ", combined_min_log10abund, "\n");

plot_text(c(
	paste("Reference File Name: ", ReferenceFilename, sep=""),
	capture.output(print(reference_info, quote=F)),
	"",
	"",
	paste("Combined Copy Count: ", combined_ref_copycount, sep=""),
	paste("Combined Min Log10(Abund): ", combined_min_log10abund, sep="")
	));

###############################################################################
# Remove reference spike ins from normalized table
print(reference_names);
print(categories);

categories_wo_spikeins=setdiff(categories, reference_names);
norm_wo_spikeins=normalized[,categories_wo_spikeins,drop=F];

###############################################################################

num_references=nrow(reference_info);

nz_ix=(normalized!=0);
min_nz=min(normalized[nz_ix]);
log10_min_nz=log10(min_nz);
cat("Min Nonzero Abundance: ", min_nz, "\n");
cat("Log10(Min Nonzero): ", log10_min_nz, "\n");


mode_of_continuous_distribution=function(x){
	x=x[is.finite(x)];
	d=density(x);
	return(d$x[which.max(d$y)]);
}
#------------------------------------------------------------------------------

par(mar=c(5,5,5,1));
par(oma=c(0,0,2,0));

sample_info_list=list();

ref_itn=1;
for(cur_ref_nm in c(reference_names, "_combined_")){
	
	par(mfrow=c(2,2));
	cat("--------------------------------------------------\n");
	cat("Working on: ", cur_ref_nm, "\n", sep="");
	cat("Looking for reference [", ref_itn, "] ...\n");

	if(cur_ref_nm!="_combined_"){
		# Look for and confirm individual references could be found
		reference_col_ix=(categories==cur_ref_nm);

		ref_index=NULL;
		if(!any(reference_col_ix)){
			cat("Error: Could not find reference category.\n");
			quit(status=-1);
		}else{
			ref_index=which(reference_col_ix);
			cat("Found.\n");
		}

		cur_ref_list=cur_ref_nm;	
		copy_num=reference_copycnts[cur_ref_nm];
		l10_minabd=reference_log10minabund[cur_ref_nm];
	}else{
		# Combined reference
		cur_ref_list=reference_names;
		copy_num=combined_ref_copycount;
		l10_minabd=combined_min_log10abund;
	}

	reference_abund=normalized[,cur_ref_list,drop=F];
	#print(reference_abund);

	coll_ref_abd=apply(reference_abund, 1, sum);
	l10_coll_ref_abd=log10(coll_ref_abd+min_nz);
	l10_cn=log10(copy_num);

	min_ref_abd=10^l10_minabd;
	ref_abv_cutoff_ix=(coll_ref_abd >= min_ref_abd);

	#----------------------------------------------------------------------
	# Generate histograms of reference abundances
	hist(coll_ref_abd, breaks=seq(0,1, length.out=20), xlab="Ref Abundance",
		main="Distr. of Reference Abundance");
	title(main=sprintf("--- Specified Min Abund: %5.3f", min_ref_abd), 
		col.main="red", line=.5, cex.main=.85);
	abline(v=10^l10_minabd, col="red", lty="dashed");

	# Generate histograms of log10(reference abundances)
	max_log10=max(c(log10(coll_ref_abd+min_nz),0));
	hist(log10(coll_ref_abd+min_nz), main="Distr of Log10(Reference Abundance)", 
		breaks=seq(log10_min_nz, max_log10, length.out=20), xlab="log10(Ref Abundance)");
	title(main=sprintf("--- Specified Min Log10(Abund): %5.3f", l10_minabd), 
		col.main="red", line=.5, cex.main=.85);
	abline(v=l10_minabd, col="red", lty="dashed");

	#----------------------------------------------------------------------
	# Calculate absolute counts
	abs_wo_refer=round(copy_num*norm_wo_spikeins/as.vector(coll_ref_abd),2);
	samp_copy_counts=apply(abs_wo_refer, 1, sum);
	log10_samp_copy_counts=log10(samp_copy_counts+1);
	max_l10cc=max(log10_samp_copy_counts);
	
	#----------------------------------------------------------------------
	# Plot sample depths vs. proportion of reference
	plot(log_sample_depths, l10_coll_ref_abd,
		xlab="Log10(Sample Depths)", ylab="Log10(Ref Abundance)",
		main="Ref Abundance vs. Sample Depth"
		);
	text(par()$usr[1], l10_minabd+par()$cxy[2]/2, adj=c(0,.4),
		"  Reliably Assayed Ref", col="blue", font=4, cex=.7);
	text(par()$usr[1], l10_minabd-par()$cxy[2]/2, adj=c(0,.4),
		"  Dubiously Assayed Ref", col="red", font=4, cex=.7);
	abline(h=seq(0,-10,-1), col="black", lty="dotted");
	abline(h=l10_minabd, col="red", lty="dashed");
	correl=cor(log_sample_depths, l10_coll_ref_abd, use="na.or.complete");
	cat("Depth vs Ref Abun Cor: ", correl, "\n");
	title(main=sprintf("Correlation: %5.3f", correl), line=1.25, cex.main=.85);
	title(main=sprintf("--- Specified Min Log10(Abund): %5.3f", l10_minabd), 
		col.main="red", line=.5, cex.main=.85);

	# Plot sample depths vs. predicted copy counts
	min_log_smp_dep=min(log_sample_depths);
	plot(log_sample_depths, log10_samp_copy_counts,
		xlab="Log10(Sample Depths)", ylab="Log10(Abs Sample Copy Counts)",
		ylim=c(range(c(log10_samp_copy_counts, l10_cn), na.rm=T)),
		main="Estimated Sample Copy Counts vs. Sample Depth"
		);
	
	text(par()$usr[1]+par()$cxy[1], log10(copy_num)+par()$cxy[2]*.5, "Samp CC > Ref CC ", 
		col="red", adj=c(0,.4), font=4, cex=.8);
	text(par()$usr[1]+par()$cxy[1], log10(copy_num)-par()$cxy[2]*.5, "Samp CC < Ref CC", 
		col="green", adj=c(0,.4), font=4, cex=.8);

	abline(h=log10(copy_num), col="blue", lty="dashed");
	correl=cor(log_sample_depths, log10_samp_copy_counts, use="na.or.complete");
	cat("Depth vs. Sample Copy Counts: ", correl, "\n");
	title(main=sprintf("Correlation: %5.3f", correl), line=1.25, cex.main=.85);
	title(main=sprintf("--- Reference Log10(Copy Number): %5.3f", log10(copy_num)), 
		col.main="blue", line=.5, cex.main=.85);

	# Label reference in outer margin
	mtext(cur_ref_nm, side=3, line=0, outer=T);

	#----------------------------------------------------------------------
	# Plot histogram of Estimated Copy Counts


	plot_hist_and_annot=function(data, title, subtitle, barcol, hxlim, hylim, hbreaks, ref_cn){
		median=median(data, na.rm=T);
		mode=mode_of_continuous_distribution(data);

		widen=function(span, prop=1){
			width=span[2]-span[1];
			margin=width*prop/2;
			return(c(span[1]-margin, span[2]+margin));	
		}

		h=hist(data,
			main="", col=barcol,
			xlim=widen(hxlim, 1.1), ylim=hylim*1.1, breaks=hbreaks,
			xlab="Log10(Copies)");

		#print(data);
		#print(h$counts);
		#print(sum(h$counts));
		num_nan=sum(is.nan(data));

		text(hxlim[1]/2, hylim[2]/2, 
			paste("[Num NaN]\n", format(round(num_nan,0), big.mark=","), sep=""), 
			font=2, col="gray20");

		title(main=title, line=3, cex.main=1.5, font.main=2);

		subtit_smpcnt_str=paste(subtitle, ": Num Samples = ", length(data), sep="");
		title(main=subtit_smpcnt_str, line=1.7, cex.main=1);
		title(main=paste("--- Ref Log10(Copy Number)") , line=.7, cex.main=.85);

		abline(v=median, col="orange", lwd=3);
		abline(v=mode, col="green");
		points(c(ref_cn, ref_cn), hylim, col="black", lty="dashed", type="l");

		if(median>mode){
			med_pos=4;
			mod_pos=2;
		}else{
			med_pos=2;
			mod_pos=4;
		}		

		h=strheight("M");
		text(mode, hylim[2]+h, 
			font=2,
			col="darkgreen",
			paste("[Mode]\n", 
				format(round(10^mode,0), big.mark=","), sep=""), pos=mod_pos);
		text(median, hylim[2]+h, 
			font=2,
			col="darkorange",
			paste("[Median]\n", 
				format(round(10^median,0), big.mark=","), sep=""), pos=med_pos);

	}

	# Compute combined ranges
	comb_hist=hist(log10_samp_copy_counts, plot=F, breaks=20);
	comb_xlim=range(comb_hist$breaks);
	comb_ylim=c(0, max(comb_hist$counts));
	comb_breaks=comb_hist$breaks;

	l10cn=log10(copy_num);

	plot_hist_and_annot(log10_samp_copy_counts, 
		"All Estim. Sample Copy Counts",
		"(All)", 
		"purple",
		comb_xlim, comb_ylim, comb_breaks, l10cn);
	plot_hist_and_annot(log10_samp_copy_counts[ref_abv_cutoff_ix], 
		"'Reliable' Est. Sample Copy Counts",
		"(ref abd >= cutoff)",
		"blue",
		comb_xlim, comb_ylim, comb_breaks, l10cn);
	plot_hist_and_annot(log10_samp_copy_counts[!ref_abv_cutoff_ix],
		"'Dubious' Est. Sample Copy Counts",
		"(ref abd < cutoff)",
		"red",
		comb_xlim, comb_ylim, comb_breaks, l10cn);


	# Label reference in outer margin
	mtext(cur_ref_nm, side=3, line=0, outer=T);
	
	#----------------------------------------------------------------------
	# Write reference-specific summary table 
	out_sumtab_fn=paste(OutputRoot, ".", ref_itn, ".summary_table.tsv", sep="");
	write_summary_file(abs_wo_refer, out_sumtab_fn);

	# Accumulate stats for downstream across references analyses 
	sample_info_list[[cur_ref_nm]]=list();
	sample_info_list[[cur_ref_nm]][["RefAbnd"]]=coll_ref_abd;
	sample_info_list[[cur_ref_nm]][["TotalCopyCounts"]]=samp_copy_counts;
	sample_info_list[[cur_ref_nm]][["RefAbdAboveCutoff"]]=ref_abv_cutoff_ix;

	ref_itn=ref_itn+1;
}

##############################################################################

output_sample_copy_count_summary=function(samp_info_rec, outfn, smp_dep){

	cat("Writing: ", outfn, "\n");

	header=c("SampleID", "SampleDepth");
	out_mat=cbind(names(smp_dep), smp_dep);
	rownames(out_mat)=names(smp_dep);

	# Keep track of the reference name, above the ref specific stats
	superheader=c("", "");

	ix=1;
	for(ref_nm in names(samp_info_rec)){

		cur_info=samp_info_rec[[ref_nm]];

		refmat=cbind(
			"",
			cur_info[["RefAbnd"]],
			cur_info[["TotalCopyCounts"]], 
			cur_info[["RefAbdAboveCutoff"]]
			);

		out_mat=cbind(out_mat, refmat);

		header=c(header, "", "Ref_Abund", "EstTotCpCnt", "AbvRefCutoff");

		index_str=sprintf("[%i] ", ix);
		superheader=c(superheader, index_str, ref_nm, "", "");

		ix=ix+1;
	}

	
	colnames(out_mat)=header;

	# Write Reference Name SuperHeader
	fh=file(outfn, "w");
	cat(file=fh, paste(superheader, collapse="\t"), "\n", sep="");
	close(fh);

	# Write stats
	write.table(out_mat, outfn, append=T, quote=F, sep="\t", row.names=F, col.names=T);
}

output_summary_fn=paste(OutputRoot, ".copy_count_stats.tsv", sep="");
output_sample_copy_count_summary(sample_info_list, output_summary_fn, sample_depths);

##############################################################################

# Generate estimated total plots across samples

generate_sample_prof_line_plot=function(samp_info_rec){

	num_rec=length(samp_info_rec);

	sample_ids=names(samp_info_rec[["_combined_"]][["TotalCopyCounts"]]);
	num_samples=length(sample_ids);
	num_references=num_rec-1;
	ref_names=setdiff(names(samp_info_rec), "_combined_");

	# sort sample ids by combined
	combined_cc=samp_info_rec[["_combined_"]][["TotalCopyCounts"]];
	sorted_combined_samp_ids=names(combined_cc[order(combined_cc, decreasing=T)]);
	combined_abvchar=samp_info_rec[["_combined_"]][["RefAbdAboveCutoff"]];
	combined_abv_char=rep(1, num_samples);
	combined_abv_char[!combined_abvchar]=4;

	# Sum up the NaN so we can annotate the figure
	nan_count_arr=numeric(length(ref_names));
	names(nan_count_arr)=ref_names;
	
	# Find ranges of copy counts across datasets
	max_est=0;
	for(ref_nm in names(samp_info_rec)){
		tcc=samp_info_rec[[ref_nm]][["TotalCopyCounts"]];
		max_est=max(max_est, tcc, na.rm=T);
		nan_count_arr[ref_nm]=sum(is.nan(tcc));
	}
	log10_max_est=log10(max_est);

	cat("Max copy count estimated: ", max_est, "\n");

	par(mfrow=c(2,1));
	par(mar=c(5, 6, 3, 1));
	ref_col=rainbow(num_references, start=0, end=4/6);
	palette(ref_col);

	#------------------------------------------------------------------------
	# Plot counts
	plot(0, type="n", ylim=c(0, max_est*1.1), xlim=c(0,num_samples+1), las=2,
		xlab="Sample", ylab="Estimated Sample Copy Count", 
		main="Variability in Estimated Copy Counts");
	abline(v=num_samples, col="grey25", lty="dashed");

	ref_ix=1;	
	for(ref_nm in ref_names){

		tcc=samp_info_rec[[ref_nm]][["TotalCopyCounts"]];
		rabvc=samp_info_rec[[ref_nm]][["RefAbdAboveCutoff"]];

		sorted_tcc=tcc[sorted_combined_samp_ids];
		sorted_rabvc=rabvc[sorted_combined_samp_ids];
		abv_char=rep(1, num_samples);
		abv_char[!sorted_rabvc]=4;

		points(1:num_samples, sorted_tcc, col=ref_ix, pch=abv_char, cex=.5);
		ref_ix=ref_ix+1;
	}

	points(1:num_samples, combined_cc[sorted_combined_samp_ids], col="black", 
		pch=combined_abv_char, cex=.75);

	legend(num_samples/8, max_est,
		legend=c(ref_names, "COMBINED","", 
			"  > Min Ref Abnd (Reliable)", 
			"  < Min Ref Abnd (Dubious, but computable)",
			"  Min Ref Abnd = 0 (NaN, not shown)"), 
		fill=c(ref_col, "black", NA, NA, NA, NA),
		border=rep(F, 6),
		pch=c(NA, NA, NA, NA, NA, 1, 4, NA), 
		cex=.66);

	#------------------------------------------------------------------------
	# Plot log counts
	plot(0, type="n", ylim=c(0, log10_max_est*1.1), xlim=c(0,num_samples+1), las=2,
		xlab="Sample", ylab="Estimated Log10(Sample Copy Count)", 
		main="Variability in Log10(Estimated Copy Counts)");
	abline(v=num_samples, col="grey25", lty="dashed");

	ref_ix=1;	
	unrel_count_arr=numeric(length(ref_names)+1);
	names(unrel_count_arr)=c(ref_names, "_combined_");

	for(ref_nm in c(ref_names, "_combined_")){

		tcc=samp_info_rec[[ref_nm]][["TotalCopyCounts"]];
		rabvc=samp_info_rec[[ref_nm]][["RefAbdAboveCutoff"]];

		sorted_tcc=tcc[sorted_combined_samp_ids];
		sorted_rabvc=rabvc[sorted_combined_samp_ids];
		abv_char=rep(1, num_samples);
		abv_char[!sorted_rabvc]=4;

		unrel_count_arr[ref_nm]=sum(!sorted_rabvc);

		if(ref_nm!="_combined_"){
			points(1:num_samples, log10(sorted_tcc+1), col=ref_ix, pch=abv_char, cex=.5);
		}
		ref_ix=ref_ix+1;
	}
	points(1:num_samples, log10(combined_cc[sorted_combined_samp_ids]+1), col="black",
		pch=combined_abv_char, cex=.75);

	cnt_stat=cbind(nan_count_arr, unrel_count_arr);

	computable_but_unreliable=apply(cnt_stat, 1, diff);
	ref_stat_str=paste(
		num_samples, " - (", cnt_stat[,2], " = ", cnt_stat[,1], " + ", 
		computable_but_unreliable, ") = ", num_samples-unrel_count_arr, 
		sep="");

	legend(num_samples*2/3, log10_max_est, 
		cex=.66,
		title="Total - (Dubious = NaN + Computable) = Reliable",
		fill=c(ref_col, "black"), 
		legend=c(ref_stat_str),
		bty="n"
		);

	text(num_samples/4, 0, "<-- Samples with more CC", col="red", font=4, cex=1);
	text(num_samples*3/4, 0, "Samples with less CC -->", col="blue", font=4, cex=1);
	
}

generate_sample_prof_line_plot(sample_info_list);

##############################################################################

generate_spaghetti_plots=F;

if(generate_spaghetti_plots){

	cat("Generate spaghetti plots...\n");

	num_samples=nrow(normalized);

	num_color_vals=ceiling(num_samples^(1/3));
	colors=rgb(runif(num_samples), runif(num_samples), runif(num_samples), maxColorValue=1);

	par(mar=c(15,4,4,1));
	par(mfrow=c(1,1));

	plot(0, type="n", xlim=c(1, num_references), ylim=c(0, 1), xaxt="n", xlab="", las=2,
		ylab="Abundance",
		main="Spaghetti Plot: Abundance");
	axis(side=1, at=1:num_references, labels=reference_names, cex.axis=1, las=2);
	for(i in 1:num_samples){
		points(1:num_references, normalized[i, reference_names], col=colors[i], type="l"); 
	}

	plot(0, type="n", xlim=c(1, num_references), ylim=c(log10_min_nz, 0), xaxt="n", xlab="", las=2,
		ylab="log10(Abundance)",
		main="Spaghetti Plot: Log10(Abundance)");
	axis(side=1, at=1:num_references, labels=reference_names, cex.axis=1, las=2);
	for(i in 1:num_samples){
		points(1:num_references, log10(normalized[i, reference_names]+min_nz), 
			col=colors[i], type="l"); 
	}

	plot(0, type="n", xlim=c(1, num_references), ylim=c(1, num_samples), xaxt="n", xlab="", las=2,
		ylab="Rank(Abundance)",
		main="Spaghetti Plot: Rank(Abundance)");
	axis(side=1, at=1:num_references, labels=reference_names, cex.axis=1, las=2);
	rank_norm=apply(normalized[,reference_names], 2, rank);
	colnames(rank_norm)=reference_names;
	for(i in 1:num_samples){
		points(1:num_references, rank_norm[i, reference_names], col=colors[i], type="l"); 
	}

	###############################################################################

	cat("Calculate Correlations among references...\n");
	spearmans_cor=cor(normalized[, reference_names]);

	par(mar=c(7,3,5,7));
	paint_matrix(spearmans_cor, "Rank Correlation between References", plot_min=0, plot_max=1,
		deci_pts=3, show_leading_zero=F);

	###############################################################################

	max_rank_diff=num_samples/10;

	plot(0, type="n", xlim=c(1, num_references), ylim=c(1, num_samples), xaxt="n", xlab="", las=2,
		ylab="Rank(Abundance)",
		main=paste("Spaghetti Plot: Rank(Abundance) Increasing > ", max_rank_diff, sep=""));
	axis(side=1, at=1:num_references, labels=reference_names, cex.axis=1, las=2);
	colnames(rank_norm)=reference_names;
	for(i in 1:num_samples){
		for(r in 1:(num_references-1)){
			if(diff(rank_norm[i, c(r,r+1)])>max_rank_diff){
				points(c(r,r+1), rank_norm[i, c(r,r+1)], col=colors[i], type="l"); 
			}
		}
	}

	plot(0, type="n", xlim=c(1, num_references), ylim=c(1, num_samples), xaxt="n", xlab="", las=2,
		ylab="Rank(Abundance)",
		main=paste("Spaghetti Plot: Rank(Abundance) Decreasing > ", max_rank_diff, sep=""));
	axis(side=1, at=1:num_references, labels=reference_names, cex.axis=1, las=2);
	colnames(rank_norm)=reference_names;
	for(i in 1:num_samples){
		for(r in 1:(num_references-1)){
			if(diff(rank_norm[i, c(r,r+1)])< -max_rank_diff){
				points(c(r,r+1), rank_norm[i, c(r,r+1)], col=colors[i], type="l"); 
			}
		}
	}

}

###############################################################################

if(!is.null(ProportionsFile)){

	cat("---------------------------------------------------------------------------\n");
	cat("Generating adjust copy counts with specified proportions file.\n");

	# Load the proportions file
	prop_tab=read.table(ProportionsFile, header=T, row.names=1, sep="\t");
	prop_arr=prop_tab[,1];
	names(prop_arr)=rownames(prop_tab);
	prop_samp_ids=rownames(prop_tab);
	num_prop_samp_ids=length(prop_samp_ids);

	cat("Num Samples with Proportion Information: ", num_prop_samp_ids, "\n");


	output_prop_adj_copy_counts=function(samp_info_rec, outfn, smp_dep, props){

		cat("Writing: ", outfn, "\n");

		# Copy props into array, and make placeholder for missing props info.
		smp_dep_ids=names(smp_dep);
		num_smp_dep_ids=length(smp_dep_ids);
		prop_to_keep_arr=rep(0, num_smp_dep_ids);
		names(prop_to_keep_arr)=smp_dep_ids;

		props_ids=names(props);
		shared_props_ids=intersect(smp_dep_ids, props_ids);
		prop_to_keep_arr[shared_props_ids]=props[shared_props_ids];

		header=c("SampleID", "SampleDepth", "PropToKeep");
		out_mat=cbind(smp_dep_ids, smp_dep, sprintf("%5.4f",prop_to_keep_arr));
		rownames(out_mat)=smp_dep_ids;
		colnames(out_mat)=header;
		print(out_mat);

		# Keep track of the reference name, above the ref specific stats
		superheader=c("", "", "", "");

		ix=1;
		for(ref_nm in names(samp_info_rec)){

			cur_info=samp_info_rec[[ref_nm]];
			cc_samp_ids=names(cur_info[["TotalCopyCounts"]]);

			shared_samp_ids=intersect(props_ids, cc_samp_ids);

			prop_adj_cc=rep(NA, length(cc_samp_ids));
			names(prop_adj_cc)=cc_samp_ids;
			prop_adj_cc[shared_samp_ids]=
				cur_info[["TotalCopyCounts"]][shared_samp_ids]*
				props[shared_samp_ids];
			
			refmat=cbind(
				"",
				sprintf("%6.4f", cur_info[["RefAbnd"]]),
				sprintf("%4.0f", cur_info[["TotalCopyCounts"]]), 
				cur_info[["RefAbdAboveCutoff"]],
				sprintf("%4.0f",prop_adj_cc)
				);

			out_mat=cbind(out_mat, refmat);

			header=c(header, "", "Ref_Abund", "EstTotCpCnt", "AbvRefCutoff", "PropAdj");

			index_str=sprintf("[%i] ", ix);
			superheader=c(superheader, index_str, ref_nm, "", "", "");

			ix=ix+1;
		}

		
		colnames(out_mat)=header;

		# Write Reference Name SuperHeader
		fh=file(outfn, "w");
		cat(file=fh, paste(superheader, collapse="\t"), "\n", sep="");
		close(fh);

		# Write stats
		write.table(out_mat, outfn, append=T, quote=F, sep="\t", row.names=F, col.names=T);
	}

	scc_prpadj_fname=paste(OutputRoot, ".sample_copy_counts.prop_adj.tsv", sep="");
	output_prop_adj_copy_counts(sample_info_list, scc_prpadj_fname, sample_depths, prop_arr);

	cat("---------------------------------------------------------------------------\n");

}else{
	cat("-----------------------------------------------------\n");
	cat("Proportions File not specified.\n");
	cat("No problem, if that was on purpose.\n");
	cat("-----------------------------------------------------\n");
}
	

###############################################################################

cat("Done.\n")
print(warnings());

q(status=0)
