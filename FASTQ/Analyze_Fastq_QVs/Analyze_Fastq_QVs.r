#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file_list", "l", 1, "character",
	"offset", "f", 1, "numeric",
	"output_filename_root", "o", 1, "character",
	"head", "h", 2, "numeric",
	"max_qv", "m", 2, "numeric",
	"validation", "v", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

DEFAULT_MAX_QV=41;

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-l <input fastq file list>\n",
	"	-f <QV offset, i.e. 33 or 64>\n",
	"	-o <output filename root>\n",
	"\n",
	"	[-m <maximum 'reasonable' QV, default=", DEFAULT_MAX_QV, "]\n",
	"	[-h <number of top reads to select (to reduce PDF files size)>]\n",
	"	[-v (validation flag)]\n",
	"\n",	
	"This script will read in a list of fastq files\n",
	"and generate statistics and plots in order to\n",
	"compare the quality of reads between samples.\n",
	"\n",
	"You can specify a list with one of two formats:\n",
	"  1.) A list of fastq files (1 column)\n",
	"  2.) A mapping of sample ID and fastq files (2 columns)\n",
	"\n",
	"Note that Phred QV's should be between 0 and 40ish.\n",
	"If the wrong offset encoding is specified, you will\n",
	"get QV's out of that range, and the problem should\n",
	"be looked into.\n",
	"\n",
	"\n");

if(
	!length(opt$input_file_list) ||
	!length(opt$offset) ||
	!length(opt$output_filename_root)
){
	cat(usage);
	q(status=-1);
}

###############################################################################

InputFileList=opt$input_file_list;
Offset=opt$offset;
OutputFilenameRoot=opt$output_filename_root;
NoValidation=!length(opt$validation);

if(length(opt$head)>0){
	Head=opt$head;
	head_ext=paste(".h", Head, sep="");
}else{
	Head=-1;
	head_ext="";
}
OutputFilenameRoot=paste(OutputFilenameRoot, head_ext, sep="");

if(length(opt$max_qv)>0){
	MaxQV=opt$max_qv;
}else{
	MaxQV=DEFAULT_MAX_QV;
}

cat("\n")
cat("Input File Name: ", InputFileList, "\n");
cat("Offset: ", Offset, "\n");
cat("Output Filename Root: ", OutputFilenameRoot, "\n");
cat("Head Values: ", Head, "\n");
cat("No Validation: ", NoValidation, "\n");
cat("Max QV allowed: ", MaxQV, "\n");
cat("\n");

###############################################################################

read_file_list=function(filelist_fname){

	filelist=as.matrix(read.table(filelist_fname, sep="\t", header=F));

	filelist_dim=dim(filelist);
	cat("File List Dimensions:\n");
	print(filelist_dim);
	cat("\n");

	# If sample names have not been specified,
	# then use the sample paths without the fastq extension
	if(filelist_dim[2]==1){
		cat("Sample names not specified.\n");
		clean_filelist=gsub("\\.fastq", "", filelist);

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

read_fastq_file=function(fastq_fname, keep_ids=F, keep_seq=F, keep_qual=F, quick=T, head=-1){

	cat("File Name: ", fastq_fname, "\n");
	
	fh=file(fastq_fname);
	open(fh);

	data=list();
	if(keep_ids){
		data[["ids"]]=character();
	}
	if(keep_seq){
		data[["seq"]]=character();
	}
	if(keep_qual){
		data[["qual"]]=character();
	}

	num_records=0;
	min_len=.Machine$integer.max;
	max_len=0;

	keep_reading=T;
	while(keep_reading){

		lines=readLines(fh, n=4, warn=F);
		num_lines_read=length(lines);

		if(num_lines_read==0){
			keep_reading=F;
		}else if(num_lines_read!=4){
			cat("Error:  Expecting 4 lines per fastq record.\n");
			cat("Around Line: ", num_records*4, "\n");
			quit(status=-1);
		}else{
			#print(lines);
			len=nchar(lines[2]);
			min_len=min(min_len, len);
			max_len=max(max_len, len);

			if(!quick){
				if(nchar(lines[2])!=nchar(lines[4])){
					cat("Error:  Sequence and quality lengths not equal.\n");
					cat("Around Line: ", num_records*4, "\n");
					cat("   Seq length: ", nchar(lines[2]), "\n");
					cat("   Qual length: ", nchar(lines[4]), "\n");
					quit(status=-1);
				}
				if(lines[3] != "+"){
					cat("Error:  Expecting '+' between sequence and quality, but not found.\n");
					cat("Around Line: ", num_records*4, "\n");
					quit(status=-1);
				}
			}

			idx=num_records+1;
			if(keep_ids){
				data[["ids"]][idx]=lines[1];
			}
			if(keep_seq){
				data[["seq"]][idx]=lines[2]
			}
			if(keep_qual){
				data[["qual"]][idx]=lines[4]
			}

			num_records=num_records+1;
			if(num_records==head){
				keep_reading=F;
			}
		}
	}

	cat("Num Records Read: ", num_records, "\n");
	cat("Min Length: ", min_len, " / Max Length: ", max_len, "\n");

	close(fh);

	data[["num_recs"]]=num_records;
	data[["min_len"]]=min_len;
	data[["max_len"]]=max_len;

	return(data);
}

###############################################################################

qvstr_to_int_str=function(qvstr, offset){
	ascii=as.integer(charToRaw(qvstr));
	return(ascii-offset);
}

qvstr_to_int_list=function(qvstr_list, offset){
	res=list();
	num_recs=length(qvstr_list);
	for(i in 1:num_recs){
		res[[i]]=qvstr_to_int_str(qvstr_list[[i]], offset);
		is_reasonable=is_reasonable_qv(res[[i]], max_qv=MaxQV);
		#is_reasonable=T; 
		if(!is_reasonable){
			qv_range=range(res[[i]]);
			cat("Error:  QV's are out of range (<0 or >", MaxQV, "). Are you sure offsets are correct?\n", sep="");
			cat("Record: ", i, "\n");
			cat("QV Range: ", qv_range[1], " - ", qv_range[2], "\n", sep="");
			quit(status=-1);
		}
	}
	return(res);
}

is_reasonable_qv=function(qvint, max_qv=41){
	if(all(qvint>=0 & qvint<=max_qv)){
		return(T);
	}else{
		return(F);
	}
}

###############################################################################

plot_qv_by_pos=function(pos_med_ci, max_len){

	max_disp_qv=65;

	plot(pos_med_ci[2,],
		xlim=c(0, max_len),
		ylim=c(0, max_disp_qv), 
		type="p", cex=.15, pch=16, col="blue",
		xlab="Position (bp)",
		ylab="Quality Value",
		main="Quality Value across Positions",
		bty="l"
	);

	# Plot bounds
	points(pos_med_ci[1,], pch=2, cex=.15, col="red");
	points(pos_med_ci[3,], pch=6, cex=.15, col="green");

	# Compute smoothing
	med_lowess=lowess(pos_med_ci[2,], f=1/6);
	lb_lowess=lowess(pos_med_ci[1,], f=1/6);
	ub_lowess=lowess(pos_med_ci[3,], f=1/6);

	# Smooth lines
	points(med_lowess, col="blue", type="l", lwd=3, lty="solid");
	points(lb_lowess, col="red", type="l", lwd=1, lty="dashed");
	points(ub_lowess, col="green", type="l", lwd=1, lty="dashed");

	legend(
		0, 65,
		fill=c("blue", "red", "green"),
		legend=c("Median", "95% Lower Bound", "95% Upper Bound"),
		cex=.7
	)

	# Compute median QVs cross positions
	med_lb=median(pos_med_ci[1,]);
	med_med=median(pos_med_ci[2,]);
	med_ub=median(pos_med_ci[3,]);

	# If QVs appear to be low, print banner
	if(med_lb<=20 || med_med<=30 || med_ub<=35){
		cat("Median LB: ", med_lb, "\n");
		cat("Median Med: ", med_med, "\n");
		cat("Median UB: ", med_ub, "\n");
		transparent_red=rgb(1,0,0, alpha=.25);
		text(max_len/2, max_disp_qv/2, 
			cex=5.5,
			col=transparent_red, font=2,
			labels="WARNING: Low QVs!"
		);
	}
	
}

plot_counts_by_pos=function(cts, max_len){

	max_cts=max(cts);
	plot(cts,
		xlim=c(0, max_len),
		ylim=c(0, max_cts),
		xlab="Position (bp)",
		ylab="Counts",
		pch=15,
		cex=.25,
		type="b",
		main="Length Falloff of Reads",
		bty="l"
	);

	if(sd(cts)/mean(cts)<.05){
		transparent_red=rgb(1,0,0, alpha=.25);
		text(max_len/2, max_cts/2, 
			cex=3.5,
			col=transparent_red, font=2,
			labels="WARNING: Homogenous Lengths!"
		);
	}

}

plot_frequency_of_qv=function(all_qv){
	qv_max=50;
	h=hist(all_qv, breaks=c(0:qv_max),
		main="QV Frequencies",
		xlab="Quality Value",
	);

	lowqv_dens=sum(h$density[0:15])
	cat("Proportion of Low QVs (<=15): ", lowqv_dens, "\n");
	if(lowqv_dens>(1/5)){
		transparent_red=rgb(1,0,0, alpha=.25);
		text(qv_max/2, max(h$counts)/2, 
			cex=3,
			col=transparent_red, font=2,
			labels="WARNING: High Proportion of Low QV!"
		);
	}
	

}

plot_len_vs_qv=function(lengths, med_qv, max_len){
	plot(lengths, med_qv,
		xlim=c(0, max_len),
		ylim=c(0, 60),
		xlab="Length",
		ylab="Median QV",
		main="Median QV vs. Length"
	);

}


#------------------------------------------------------------------------------

analyze_qv_list=function(qvint_list, max_len, name){

	num_recs=length(qvint_list);

	combined=matrix(NA, nrow=num_recs, ncol=max_len);	
	seq_len=numeric(num_recs);

	for(i in 1:num_recs){
		seq_len[i]=length(qvint_list[[i]]);	
		combined[i,1:seq_len[i]]=qvint_list[[i]][1:seq_len[i]];
	}

	# Generate per read statistics
	med_ci=apply(combined, 2, function(x){quantile(x, probs=c(0.025, .5, 0.975), na.rm=T)});	
	cts=apply(combined, 2, function(x){sum(!is.na(x))});
	med_per_seq=apply(combined, 1, function(x){median(x,na.rm=T)});
	all_qv=as.vector(combined);
	all_qv=all_qv[!is.na(all_qv)];

	# Plot per read statistics
	par(oma=c(0,0,3,0));
	par(mfrow=c(4,1));
	plot_qv_by_pos(med_ci, max_len);
	plot_counts_by_pos(cts, max_len);
	plot_frequency_of_qv(all_qv);
	plot_len_vs_qv(seq_len, med_per_seq, max_len);

	title_size=2-nchar(name)/40;
	mtext(name, side=3, outer=T, cex=title_size);
		
	# Generate per sample statistics
	sample_qv_ci=quantile(all_qv, probs=c(0.025, .5, 0.975));
	sample_len_ci=quantile(seq_len, probs=c(0.025, .5, 0.975));
	
	samp_summ=list();
	samp_summ[["qv_ci"]]=sample_qv_ci;
	samp_summ[["len_ci"]]=sample_len_ci;
	return(samp_summ);
}

#------------------------------------------------------------------------------

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

plot_qv=function(qv_ci, sample_names, label.cex=1, sorted=F){

	max_qv=max(qv_ci, na.rm=T);

	num_samples=nrow(qv_ci);
	mids=barplot(qv_ci[,2], plot=F);
	spacing=mids[2]-mids[1];
	
	spacing=ifelse(is.na(spacing), 2, spacing);

	plot(0,0, type="n", 
		bty="l",
		xlim=c(0, (num_samples+.25)*spacing), ylim=c(0, max_qv*1.2),
		xaxt="n", yaxt="n", xlab="", ylab="");

	if(sorted){
		samp_ix=order(qv_ci[,2], decreasing=T);
		qv_ci=qv_ci[samp_ix,];
		sample_names=sample_names[samp_ix];
	}


	par(lwd=.05);
	barplot(qv_ci[,2], names.arg=sample_names, las=2, col="white",
		ylab="Median QV",
		main="Median Quality Values",
		add=T,
		cex.names=label.cex
	);
	points(mids, qv_ci[,1], pch="+", cex=2, col="red");
	points(mids, qv_ci[,3], pch="+", cex=2, col="green");

}

plot_len=function(len_ci, sample_names, label.cex=1, sorted=F){
	max_len=max(len_ci, na.rm=T);

	num_samples=nrow(len_ci);
	mids=barplot(len_ci[,2], plot=F);
	spacing=mids[2]-mids[1];
	spacing=ifelse(is.na(spacing), 2, spacing);

	plot(0,0, type="n", 
		bty="l",
		xlim=c(0, (num_samples+.25)*spacing), ylim=c(0, max_len*1.2),
		xaxt="n", yaxt="n", xlab="", ylab="");

	if(sorted){
		samp_ix=order(len_ci[,2]);
		len_ci=len_ci[samp_ix,];
		sample_names=sample_names[samp_ix];
	}

	par(lwd=.05);
	barplot(len_ci[,2], names.arg=sample_names, las=2, col="white",
		ylab="Length (bps)",
		main="Median Lengths",
		add=T,
		cex.names=label.cex
	);
	points(mids, len_ci[,1], pch="+", cex=2, col="red");
	points(mids, len_ci[,3], pch="+", cex=2, col="green");
}

plot_qnorm=function(values, sample_names, data_name){
	pts=qqnorm(values, plot.it=F);
	plot(pts, 
		xlab="Normal Standard Deviations",
		ylab=data_name,
		main=paste("Sample ", data_name, " Normal Q-Q Plot", sep=""),
		xlim=range(pts$x, na.rm=T)*c(1.2, 1.2), ylim=range(pts$y, na.rm=T)*c(1/1.2, 1.2), cex=.5, col="red");
	qqline(values, col="blue");
	text(pts$x, pts$y, sample_names, col="black", cex=.5, pos=1);
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

pdf(paste(OutputFilenameRoot, ".indiv.qv.pdf", sep=""), height=11, width=8.5);

file_list=read_file_list(InputFileList);
sample_names=rownames(file_list);
filepaths=file_list[,1];

num_samples=length(sample_names);

num_reads_per_sample=numeric(num_samples);
qv_ci=matrix(NA, nrow=num_samples, ncol=3);
len_ci=matrix(NA, nrow=num_samples, ncol=3);
max_samp_name_len=0;

for(i in 1:num_samples){
	cat("Working on loading: ", sample_names[i], "\n");
	max_samp_name_len=max(max_samp_name_len, nchar(sample_names[i]));

	try_res=try({data=read_fastq_file(filepaths[i], keep_qual=T, head=Head, quick=NoValidation)});
	if(class(try_res)=="try-error"){
		cat("Could not open:", filepaths[i], "\n");
		cat("Skipping...\n\n");
		next;
	}

	if(data[["num_recs"]]>0){
		qv_int=qvstr_to_int_list(data[["qual"]], Offset);
		summary=analyze_qv_list(qv_int, data[["max_len"]], sample_names[i]);
		qv_ci[i,]=summary[["qv_ci"]];
		len_ci[i,]=summary[["len_ci"]];
	}else{
		qv_ci[i,]=c(NA,NA,NA);
		len_ci[i,]=c(NA,NA,NA);
	}

	num_reads_per_sample[i]=data[["num_recs"]];
	
	cat("\n");
}

dev.off();

###############################################################################

pdf(paste(OutputFilenameRoot, ".overall.qv.pdf", sep=""), height=11, width=8.5);

#print(num_reads_per_sample);
#print(qv_ci)
#print(len_ci);

par(mfrow=c(3,1));
par(oma=c(0,0,3,0));
def_mar=par()$mar;
fat_bottom=def_mar;

fat_bottom[1]=fat_bottom[1]*max(1, max_samp_name_len/7);
if(fat_bottom[1]>10){
	label.cex=(cex=9/fat_bottom[1]);
	fat_bottom[1]=10;
}else{
	label.cex=1;
}

par(mar=fat_bottom);

# Plot by name order
plot_reads_per_sample(num_reads_per_sample, sample_names, label.cex=label.cex);
plot_qv(qv_ci, sample_names, label.cex=label.cex);
plot_len(len_ci, sample_names, label.cex=label.cex);
mtext(OutputFilenameRoot, side=3, outer=T, cex=2);

#-------------------------------------------------------------------------------

# Plot sorted by metric
plot_reads_per_sample(num_reads_per_sample, sample_names, label.cex=label.cex, sorted=T);
plot_qv(qv_ci, sample_names, label.cex=label.cex, sorted=T);
plot_len(len_ci, sample_names, label.cex=label.cex, sorted=T);
mtext(OutputFilenameRoot, side=3, outer=T, cex=2);

#-------------------------------------------------------------------------------

# Generate QQ plot for median read lengths
par(mfrow=c(3,1));
plot_qnorm(num_reads_per_sample, sample_names, "Reads/Sample");
plot_qnorm(qv_ci[,2], sample_names, "Read QVs");
plot_qnorm(len_ci[,2], sample_names, "Read Lengths");

mtext(OutputFilenameRoot, side=3, outer=T, cex=2);

#-------------------------------------------------------------------------------

# Generate histograms across samples
#print(num_reads_per_sample);
h=hist(num_reads_per_sample, main="Distribution of Reads/Sample", xlab="Reads/Sample", ylab="Counts");
#print(h);
cxy=par()$cxy;
if(Head!=-1){
	abline(v=Head, col="blue");
	max_bin_height=max(h$counts);
	text(Head-cxy[1]/4, max_bin_height/2, paste("Req Max Subs Size: ", Head, sep=""), cex=.7, srt=90, pos=3, col="blue");
}

#print(qv_ci);
hist(qv_ci[,2], main="Distribution of Median Sample QVs", xlab="Median Sample QV", ylab="Counts");

#print(len_ci);
hist(len_ci[,2], main="Distribution of Median Sample Lengths", xlab="Median Sample Lengths", ylab="Counts");

mtext(OutputFilenameRoot, side=3, outer=T, cex=2);

#-------------------------------------------------------------------------------

key_cts=plot_cdf_and_thresholds(num_reads_per_sample);
mtext(OutputFilenameRoot, side=3, outer=T, cex=2);

out_tab=t(c(OutputFilenameRoot, key_cts));

write.table(out_tab, paste(OutputFilenameRoot, ".key_counts.tsv", sep=""), quote=F, 
	sep="\t", row.names=F);

###############################################################################
# Output Spreadsheets

combined=cbind(
	qv_ci[,c(2,1,3), drop=F], 
	len_ci[,c(2,1,3), drop=F], 
	matrix(num_reads_per_sample, ncol=1));

colnames(combined)=c(
	"QV Median", "QV 95% LB", "QV 95% UB",
	"Length Median", "Length 95% LB", "Length 95% UB",
	"Num Reads Analyzed"
);

rownames(combined)=sample_names;

write.table(combined, paste(OutputFilenameRoot, ".qv.tsv", sep=""), quote=F, sep="\t", 
	row.names=T, col.names=NA);

###############################################################################

cat("Done.\n")
print(warnings());

q(status=0)
