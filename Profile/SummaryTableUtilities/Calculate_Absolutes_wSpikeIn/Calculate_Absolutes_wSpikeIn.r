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

plot_text(c(
	paste("Reference File Name: ", ReferenceFilename, sep=""),
	capture.output(print(reference_info, quote=F))
	));

###############################################################################
# Remove reference spike ins from normalized table
print(reference_names);
print(categories);

categories_wo_spikeins=setdiff(categories, reference_names);
norm_wo_spikeins=normalized[,categories_wo_spikeins,drop=F];

###############################################################################

num_references=nrow(reference_info);
par(mfrow=c(3,2));

nz_ix=(normalized!=0);
min_nz=min(normalized[nz_ix]);
log10_min_nz=log10(min_nz);
cat("Min Nonzero Abundance: ", min_nz, "\n");
cat("Log10(Min Nonzero): ", log10_min_nz, "\n");

par(mar=c(5,5,5,1));

ref_ix=1;
for(ref in reference_names){

	cat("Working on: ", ref, "\n", sep="");
	cat("Looking for reference [", ref_ix, "] ...\n");

	reference_ix=(categories==ref);

	ref_index=NULL;
	if(!any(reference_ix)){
		cat("Error: Could not find reference category.\n");
		quit(status=-1);
	}else{
		ref_index=which(reference_ix);
		cat("Found.\n");
	}

	reference_abund=normalized[,ref,drop=F];
	mean_refer_abund=mean(reference_abund);

	cat("Reference Abundances:\n");
	print(reference_abund);
	cat("\n");
	hist(reference_abund, main=ref, breaks=seq(0,1, length.out=20), xlab="Abundance");
	hist(log10(reference_abund+min_nz), main=ref, 
		breaks=seq(log10_min_nz, 0, length.out=20), xlab="log10(Abundance)");

	cat("Mean Reference Abundance:\n");
	print(mean_refer_abund);
	cat("\n");

	cat("Spike-in Copy Count: ", reference_copycnts[ref], "\n");
	reference_counts=reference_abund*reference_copycnts[ref];
	cat("Reference Counts:\n");
	print(reference_counts);
	cat("\n");

	#abs_counts=(normalized/as.vector(reference_abund))*SpikeinCount;
	#print(abs_counts[,ReferenceCategory,drop=F]);

	abs_wo_refer=reference_copycnts[ref]*norm_wo_spikeins/as.vector(reference_abund);
	#print(abs_wo_refer);

	# Write reference-specific summary table 
	out_sumtab_fn=paste(OutputRoot, ".", ref_ix, ".summary_table.tsv", sep="");
	write_summary_file(abs_wo_refer, out_sumtab_fn);

	ref_ix=ref_ix+1;
}

###############################################################################

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
