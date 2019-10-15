#!/usr/bin/env Rscript

###############################################################################

library('getopt');
library('seqinr');

params=c(
	"input_file_list", "l", 1, "character",
	"output_filename_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-l <input fasta file list>\n",
	"	-o <output filename root>\n",
	"\n",	
	"This script will read in a list of fasta files\n",
	"and generate statistical plots on run quality\n",
	"based on sequencing depth.  For example, the list of the\n",
	"fasta files can be from the same sequencing run.\n",
	"\n");

if(
	!length(opt$input_file_list) ||
	!length(opt$output_filename_root)
){
	cat(usage);
	q(status=-1);
}

###############################################################################

InputFileList=opt$input_file_list;
OutputFilenameRoot=opt$output_filename_root;

cat("\n")
cat("Input File Name: ", InputFileList, "\n");
cat("Output Filename Root: ", OutputFilenameRoot, "\n");
cat("\n");

###############################################################################

read_file_list=function(filelist_fname){

	filelist=as.matrix(read.table(filelist_fname, sep="\t", header=F));

	filelist_dim=dim(filelist);
	cat("File List Dimensions:\n");
	print(filelist_dim);
	cat("\n");

	# If sample names have not been specified,
	# then use the sample paths without the fasta extension
	if(filelist_dim[2]==1){
		cat("Sample names not specified.\n");
		clean_filelist=gsub("\\.fasta", "", filelist);

		num_files=length(clean_filelist);
		for(fix in 1:num_files){
			clean_filelist[fix]=tail(strsplit(clean_filelist[fix], "/")[[1]],1);
		}

		filelist=cbind(clean_filelist, filelist);
	}

	out_list=matrix(filelist[,2], ncol=1);
	rownames(out_list)=filelist[,1];

	return(out_list);
}
###############################################################################

read_fasta_file=function(fasta_fname){

	cat("File Name: ", fasta_fname, "\n");

	seq_data=tryCatch({
		read.fasta(fasta_fname);
	}, error=function(e){
		cat("FASTA file was empty.\n");
		NULL;
	});

	return(seq_data);
}

###############################################################################

plot_reads_per_sample=function(num_reads, sample_names, label.cex=1, sorted=F){

	if(sorted){
		samp_ix=order(num_reads, decreasing=T);
		num_reads=num_reads[samp_ix];
		sample_names=sample_names[samp_ix];
	}

	par(lwd=.05);
	barplot(num_reads, names.arg=sample_names, las=2, col="white",
		ylab="Number of Reads",
		main="Reads per Sample",
		cex.names=label.cex
	);
}

plot_cdf_and_thresholds=function(num_readspersamp){

	ord_ix=order(num_readspersamp, decreasing=T);
	num_readspersamp=num_readspersamp[ord_ix];
	num_samps=length(num_readspersamp);

	key_counts=c(500, 1000, 3000, 10000, 50000, 100000);
	key_colors=c("red", "orange", "green", "blue", "purple", "pink");
	num_key_counts=length(key_counts);
	log_kc=log10(key_counts);
	
	plot_ymax=max(log_kc, log10(num_readspersamp));

	par(mfrow=c(2,1));

	# Plot CDF
	par(mar=c(5,5,10,7));
	plot(100*((1:num_samps)/num_samps), log10(num_readspersamp),
		xlab="Percent of Samples",
		ylab="Log10(Reads/Sample)",
		main=paste("Total Samples: ", num_samps, sep=""),
		ylim=c(0, plot_ymax),
		xlim=c(0, 100),
		type="l",
		);

	abline(h=log_kc, col="blue", lwd=.5, lty="dotted");
	axis(side=4, at=log_kc, labels=key_counts, cex.axis=.7, las=2);

	# Plot a few key bar plots
	perc=numeric(num_key_counts);
	num_abv_keycts=numeric(num_key_counts);
	for(i in 1:num_key_counts){
		num_abv_keycts[i]=sum(num_readspersamp>key_counts[i])
		perc[i]=round(100*num_abv_keycts[i]/num_samps, 2);
		abline(v=perc[i], col=key_colors[i]);
		axis(side=3, at=num_abv_keycts[i]/num_samps*100, labels=num_abv_keycts[i], cex.axis=.7, las=2);
	}

	par(mar=c(5,5,5,7));
	bmids=barplot(perc, names.arg=paste(">",key_counts, sep=""), 
		ylim=c(0,105), col=key_colors,
		xlab="Reads/Sample Thresholds", ylab="Percent Meeting/Exceeding Threshold");
	abline(h=100, col="grey", lty="dotted");
	text(bmids, perc, 
		paste(num_abv_keycts, " rds\n",round(perc, 1), "%", sep=""), 
		pos=3, cex=.75, font=2);

	names(perc)=key_counts;

	return(perc);

}

###############################################################################

pdf(paste(OutputFilenameRoot, ".fa.depth.pdf", sep=""), height=11, width=8.5);

# Load sample list
file_list=read_file_list(InputFileList);
sample_names=rownames(file_list);
filepaths=file_list[,1];
num_samples=length(sample_names);

num_reads_per_sample=numeric(num_samples);
num_reads_per_sample=numeric(num_samples);

# Load FASTA file
for(i in 1:num_samples){
	cat("Working on loading: ", sample_names[i], "\n");
	cat("   Path:", filepaths[i], "\n");

	data=read_fasta_file(filepaths[i]);
	if(is.null(data)){
		num_reads=0;
	}else{
		num_reads=length(data);
	}
	cat("Num_Reads: ", num_reads, "\n");
	num_reads_per_sample[i]=num_reads
}

###############################################################################

# Plot depth stats
par(oma=c(1,1,3,1));
key_cts=plot_cdf_and_thresholds(num_reads_per_sample);
mtext(OutputFilenameRoot, side=3, outer=T, cex=2);


# Output tab file
out_tab=t(c(OutputFilenameRoot, key_cts));
write.table(out_tab, paste(OutputFilenameRoot, ".fa.depth.key_counts.tsv", sep=""), quote=F, 
	sep="\t", row.names=F);

###############################################################################

cat("Done.\n")
print(warnings());

q(status=0);
