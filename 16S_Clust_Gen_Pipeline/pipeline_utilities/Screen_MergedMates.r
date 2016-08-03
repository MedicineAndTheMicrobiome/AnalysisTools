#!/usr/bin/env Rscript

###############################################################################

library('getopt');
PROP_OVERLAP_MISM=0.2;
MAX_NS=4;
MIN_OVERLAP=25;

params=c(
	"input_file", "i", 1, "character",
	"output_fname_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input report>\n",
	"	[-o <output root]\n",
	"\n",	
	"This script will read in a report file from Mothur's make.contig()\n",
	"command and generate some statistics and graphs for diagnosis.\n",
	"\n",
	"Generates a read id keeplist that filters by:\n",
	"	Proportion Overlap Mismatch > ", PROP_OVERLAP_MISM, " bp\n",
	"	Maximum N's allowed = ", MAX_NS, "\n",
	"	Minimum Required Overlap = ", MIN_OVERLAP, "\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileNameRoot=opt$output_fname_root;

if(length(opt$output_fname_root)){
	OutputFileNameRoot=opt$output_fname_root;
}else{
	OutputFileNameRoot=gsub("\\.report$","", InputFileName);
}

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name Root: ", OutputFileNameRoot, "\n");       

###############################################################################
###############################################################################
# Load data
cat("Loading Report: ", InputFileName, "\n");
mat=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, row.names=1));

num_contigs=nrow(mat);

# subsamp
subsample=F;
if(subsample){
    ss=200;
    mat=mat[sample(num_contigs, ss),];
    num_contigs=ss;
}

cat("Num Contigs Described: ", num_contigs, "\n");
prop_overlap_mismatch=mat[,"MisMatches"]/mat[,"Overlap_Length"];
mat=cbind(mat, prop_overlap_mismatch);

pdf(paste(OutputFileNameRoot, ".pdf", sep=""), height=11, width=8.5);
par(oma=c(0,0,2,0));

###############################################################################
# Contig Lengths

cat("Plotting Raw Contig Length Distribution...\n");
par(mfrow=c(2,1));
h=hist(mat[,"Length"], breaks=200, plot=F);
lab_freq=seq(1,length(h$counts), length.out=15);

m=barplot(h$counts, space=0, ylab="Frequency", main="Contig Lengths"); 
axis(side=1, at=m[lab_freq], labels=h$mids[lab_freq]);

nz=h$counts;
nz[nz==0]=1;
m=barplot(log10(nz), space=0, ylab="Log10(Frequency)");
axis(side=1, at=m[lab_freq], labels=h$mids[lab_freq]);
mtext("Raw", outer=TRUE);

###############################################################################
# Contig statistics

print(head(mat));

cat("Plotting Report Statistics...\n");

plot_all_stats=function(mat, mar_text){
	par(mfrow=c(4,2));
	plot(mat[,"Length"], mat[,"Num_Ns"], xlab="Contig Length", ylab="Num Ns");
	hist(mat[,"Num_Ns"], main="Num Ns Frequencies", xlab="Mismatch (bp)");
	plot(mat[,"Length"], mat[,"prop_overlap_mismatch"], 
		xlab="Contig Length", ylab="Proportion of Overlap that Mismatch");
	hist(mat[,"prop_overlap_mismatch"], main="Prop Mismatch in Overlap Frequencies", xlab="Proportions");
	plot(mat[,"Length"], mat[,"Overlap_Length"], xlab="Contig Length", ylab="Overlap Length" );
	hist(mat[,"Overlap_Length"], main="Overlap Length Frequencies", xlab="Overlap (bp)");
	plot(mat[,"MisMatches"], mat[,"Num_Ns"], xlab="Num Mismatches", ylab="Num N's", main="Mismatch to N Conversion");

	if(nrow(mat)>2){
		m=lm(mat[,"Num_Ns"]~mat[,"MisMatches"]);
		abline(m, col="red");
		mtext(sprintf("Mismatch Conversion Failure Rate: %3.2f (Smaller is Better)", m$coefficients[2]), cex=.5);
		mtext(mar_text, outer=TRUE);
	}
}

plot_all_stats(mat, mar_text="Raw");

###############################################################################

# Filter by proportion mismatch
cat("Filtering poor quality contigs...\n");
keep_list=rep(TRUE, num_contigs);
keep_list[mat[,"prop_overlap_mismatch"]>PROP_OVERLAP_MISM]=FALSE;
keep_list[mat[,"Num_Ns"]>MAX_NS]=FALSE;
keep_list[mat[,"Overlap_Length"]<MIN_OVERLAP]=FALSE;
mat_filt=mat[keep_list,, drop=F];

###############################################################################

# Replot graphs
cat("Plotting Filtered Contig Length Distribution...\n");
par(mfrow=c(2,1));

if(nrow(mat_filt)==0){
	plot(0,0, type="n", xlab="", ylab="", yaxt="n", xaxt="n");
	text(0,0, "No Data Remaining after filtering.");
}else{
	h=hist(mat_filt[,"Length"], breaks=200, plot=F);

	lab_freq=seq(1,length(h$counts), length.out=15);

	m=barplot(h$counts, space=0, ylab="Frequency", main="Contig Lengths"); 
	axis(side=1, at=m[lab_freq], labels=h$mids[lab_freq]);

	nz=h$counts;
	nz[nz==0]=1;
	m=barplot(log10(nz), space=0, ylab="Log10(Frequency)");
	axis(side=1, at=m[lab_freq], labels=h$mids[lab_freq]);
	mtext("Filtered", outer=TRUE);

	###############################################################################

	cat("Plotting Filtered Report Statistics...\n");
	plot_all_stats(mat_filt, mar_text="Filtered");
}

###############################################################################

# Output read IDs to keep
fh=file(paste(OutputFileNameRoot, ".keep.list", sep=""),"w");
samp_names=rownames(mat_filt);
num_kept=length(samp_names);
cat("Num Kept: ", num_kept, "\n");

if(num_kept>0){
	for(i in 1:num_kept){
		cat(file=fh, samp_names[i], "\n", sep="");
	}
	close(fh);
}

# Output overall summary
fh=file(paste(OutputFileNameRoot, ".summary.report", sep=""), "w");

num_ctgs_raw=nrow(mat);
med_len_raw=median(mat[,"Length"]);
med_overlap_raw=median(mat[,"Overlap_Length"]);

num_ctgs_filt=nrow(mat_filt);
med_len_filt=median(mat_filt[,"Length"]);
med_overlap_filt=median(mat_filt[,"Overlap_Length"]);

cat(file=fh, paste(c("# Output", "Raw_Num_Ctgs", "Raw_Med_Ctg_Len", "Raw_Med_Ovlp_Len", "Filt_Num_Ctgs", "Filt_Med_Ctg_Len", "Filt_Med_Ovlp_Len"), collapse="\t"), "\n", sep="");

cat(file=fh, paste(c(OutputFileNameRoot, num_ctgs_raw, med_len_raw, med_overlap_raw, num_ctgs_filt, med_len_filt, med_overlap_filt), collapse="\t"), "\n", sep="");

close(fh);

###############################################################################

q(status=0);
