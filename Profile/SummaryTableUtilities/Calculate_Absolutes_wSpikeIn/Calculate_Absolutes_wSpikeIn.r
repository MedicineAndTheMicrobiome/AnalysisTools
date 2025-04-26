#!/usr/bin/env Rscript

library('getopt');

source("~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r");
source("~/git/AnalysisTools/Metadata/OutputFileLibrary/OutputFileLibrary.r");

###############################################################################

params=c(
	"summary_table", "s", 1, "character",
	"reference_si_fn", "r", 1, "character",
	"output_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-s <input summary_table.tsv>\n",
	"	-r <reference spike-in table filename>\n",
	"	-o <output file root>\n",
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
	"Note that the concentration that was spiked in is not important, rather it's the\n",
	"number of 16S (or reference DNA) copies that was mixed into the experimental sample.\n",
	"In the case of a the Zymo protocol, 20uL of reagent was spiked in, so the copies\n",
	"is based on cell density (cells/20 uL) x 16S copies per cell.\n",
	"\n",
	"The volume (or concentration) of the recipient experimental sample is also not\n",
	"important, as long as the sampling unit, e.g. swap, scoop, drop, of the experimental\n",
	"sample is consistent.\n",
	"\n",
	"For each reference taxa, the taxa counts will be estimated and averaged together.\n",
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

cat("\n");
cat("Input Summary Table: ", InputSummaryTable, "\n");
cat("Reference Category: ", ReferenceFilename, "\n");
cat("Output Root: ", OutputRoot, "\n");
cat("\n");

pdf(paste(OutputRoot, ".abs_spikein.pdf", sep=""), height=8.5, width=11);

###############################################################################

plot_text(c(
	script_name,
	"",
	paste("Input Summary Table: ", InputSummaryTable),
	paste("Reference Category: ", ReferenceFilename),
	paste("Output Root: ", OutputRoot)
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

	# Generate histograms of reference abundances
	hist(coll_ref_abd, breaks=seq(0,1, length.out=20), xlab="Ref Abundance",
		main="Distr. of Reference Abundance");
	title(main=sprintf("--- Specified Min Abund: %5.3f", min_ref_abd), 
		col.main="red", line=.5, cex.main=.85);
	abline(v=10^l10_minabd, col="red", lty="dashed");

	# Generate histograms of log10(reference abundances)
	hist(log10(coll_ref_abd+min_nz), main="Distr of Log10(Reference Abundance)", 
		breaks=seq(log10_min_nz, 0, length.out=20), xlab="log10(Ref Abundance)");
	title(main=sprintf("--- Specified Min Log10(Abund): %5.3f", l10_minabd), 
		col.main="red", line=.5, cex.main=.85);
	abline(v=l10_minabd, col="red", lty="dashed");

	# Calculate absolute counts
	abs_wo_refer=round(copy_num*norm_wo_spikeins/as.vector(coll_ref_abd),2);
	samp_copy_counts=apply(abs_wo_refer, 1, sum);
	log10_samp_copy_counts=log10(samp_copy_counts+1);
	max_l10cc=max(log10_samp_copy_counts);
	
	# Plot sample depths vs. proportion of reference
	plot(log_sample_depths, l10_coll_ref_abd,
		xlab="Log10(Sample Depths)", ylab="Log10(Ref Abundance)",
		main="Ref Abundance vs. Sample Depth"
		);
	abline(h=seq(0,-10,-1), col="black", lty="dotted");
	abline(h=l10_minabd, col="red", lty="dashed");
	correl=cor(log_sample_depths, l10_coll_ref_abd, use="na.or.complete");
	cat("Depth vs Ref Abun Cor: ", correl, "\n");
	title(main=sprintf("Correlation: %5.3f", correl), line=1.25, cex.main=.85);
	title(main=sprintf("--- Specified Min Log10(Abund): %5.3f", l10_minabd), 
		col.main="red", line=.5, cex.main=.85);

	mtext(cur_ref_nm, side=3, line=0, outer=T);


	# Plot sample depths vs. predicted copy counts
	plot(log_sample_depths, log10_samp_copy_counts,
		xlab="Log10(Sample Depths)", ylab="Log10(Abs Sample Copy Counts)",
		ylim=c(range(c(log10_samp_copy_counts, l10_cn), na.rm=T)),
		main="Sample Copy Counts vs. Sample Depth"
		);
	abline(h=log10(copy_num), col="blue", lty="dashed");
	correl=cor(log_sample_depths, log10_samp_copy_counts, use="na.or.complete");
	cat("Depth vs. Sample Copy Counts: ", correl, "\n");
	title(main=sprintf("Correlation: %5.3f", correl), line=1.25, cex.main=.85);
	title(main=sprintf("--- Reference Log10(Copy Number): %5.3f", log10(copy_num)), 
		col.main="blue", line=.5, cex.main=.85);

	#plot_area=par()$usr;
	#legend(plot_area[1], plot_area[4], bty="n",
	#	legend="Ref Copy Num", col="blue", lty="dashed");

	# Plot histogram of Estimated Copy Counts
	comb_hist=hist(log10_samp_copy_counts, main="Distr. of All Estimated Sample Copy Counts",
		xlab="Log10(Copies)", breaks=20);

	comb_plot_lim=par()$usr;
	comb_xlim=comb_plot_lim[c(1,2)];
	comb_ylim=comb_plot_lim[c(3,4)];
	comb_breaks=comb_hist$breaks;

	hist(log10_samp_copy_counts[ref_abv_cutoff_ix], 
		main="", col="blue",
		xlim=comb_xlim, ylim=comb_ylim, breaks=comb_breaks,
		xlab="Log10(Copies)");
	title(main="Distr. of 'Reliable' Est. Sample Copy Counts", line=1.5);
	title(main="(ref abd >= cutoff)", line=.75, cex.main=.85)
	

	hist(log10_samp_copy_counts[!ref_abv_cutoff_ix], 
		main="", col="red",
		xlim=comb_xlim, ylim=comb_ylim, breaks=comb_breaks,
		xlab="Log10(Copies)");
	title(main="Distr. of 'Dubious' Est. Sample Copy Counts", line=1.5);
	title(main="(ref abd < cutoff)", line=.75, cex.main=.85)
	
	# Write reference-specific summary table 
	out_sumtab_fn=paste(OutputRoot, ".", ref_itn, ".summary_table.tsv", sep="");
	write_summary_file(abs_wo_refer, out_sumtab_fn);

	sample_info_list[[cur_ref_nm]]=list();
	sample_info_list[[cur_ref_nm]][["RefAbnd"]]=coll_ref_abd;
	sample_info_list[[cur_ref_nm]][["TotalCopyCounts"]]=samp_copy_counts;
	sample_info_list[[cur_ref_nm]][["RefAbdAboveCutoff"]]=ref_abv_cutoff_ix;

	ref_itn=ref_itn+1;
}

cat("--------------------------------------------------\n");

##############################################################################

output_sample_copy_count_summary=function(samp_info_rec, outfn, smp_dep){

	cat("Writing: ", outfn, "\n");

	header=c("SampleID", "SampleDepth");
	out_mat=cbind(names(smp_dep), smp_dep);
	rownames(out_mat)=names(smp_dep);

	# Keep track of the reference name, above the ref specific stats
	superheader=c("", "");

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

		superheader=c(superheader, "", ref_nm, "", "");
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
	points(1:num_references, log10(normalized[i, reference_names]+min_nz), col=colors[i], type="l"); 
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

###############################################################################

###############################################################################

cat("Done.\n")
print(warnings());

q(status=0)
