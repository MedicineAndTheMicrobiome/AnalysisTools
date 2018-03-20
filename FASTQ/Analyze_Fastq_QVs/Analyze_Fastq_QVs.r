#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file_list", "l", 1, "character",
	"offset", "f", 1, "numeric",
	"output_filename_root", "o", 1, "character",
	"head", "h", 2, "numeric",
	"validation", "v", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-l <input fastq file list>\n",
	"	-f <QV offset, i.e. 33 or 64>\n",
	"	-o <output filename root>\n",
	"	[-h <number of top reads to select>]\n",
	"	[-v (validation flag)]\n",
	"\n",	
	"This script will read in a list of fastq files\n",
	"and generate statistics and plots in order to\n",
	"compare the quality of reads between samples.\n",
	"\n",
	"Note that QV's should be between 0 and 40ish.\n",
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
}else{
	Head=-1;
}

cat("\n")
cat("Input File Name: ", InputFileList, "\n");
cat("Offset: ", Offset, "\n");
cat("Output Filename Root: ", OutputFilenameRoot, "\n");
cat("Head Values: ", Head, "\n");
cat("No Validation: ", NoValidation, "\n");
cat("\n");

###############################################################################

read_file_list=function(filelist_fname){
	filelist=as.matrix(read.table(filelist_fname, sep="\t", header=T, row.names=1));
	return(filelist);
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
		is_reasonable=is_reasonable_qv(res[[i]]);
		if(!is_reasonable){
			cat("Error:  QV's are out of range. Are you sure offsets are correct?\n");
			cat("Record: ", i, "\n");
			quit(status=-1);
		}
	}
	return(res);
}

is_reasonable_qv=function(qvint){
	if(all(qvint>=0 & qvint<=41)){
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

	plot(cts,
		xlim=c(0, max_len),
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
		text(max_len/2, mean(cts), 
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
	mtext(name, side=3, outer=T, cex=2);
		
	# Generate per sample statistics
	sample_qv_ci=quantile(all_qv, probs=c(0.025, .5, 0.975));
	sample_len_ci=quantile(seq_len, probs=c(0.025, .5, 0.975));
	
	samp_summ=list();
	samp_summ[["qv_ci"]]=sample_qv_ci;
	samp_summ[["len_ci"]]=sample_len_ci;
	return(samp_summ);
}

#------------------------------------------------------------------------------

plot_reads_per_sample=function(num_reads, sample_names, label.cex=1){
	barplot(num_reads, names.arg=sample_names, las=2, col="white",
		ylab="Number of Reads",
		main="Reads per Sample",
		cex.names=label.cex
	);
}

plot_qv=function(qv_ci, sample_names, label.cex=1){

	max_qv=max(qv_ci);

	num_samples=nrow(qv_ci);
	mids=barplot(qv_ci[,2], plot=F);
	spacing=mids[2]-mids[1];
	
	spacing=ifelse(is.na(spacing), 2, spacing);

	plot(0,0, type="n", 
		bty="l",
		xlim=c(0, (num_samples+.25)*spacing), ylim=c(0, max_qv*1.2),
		xaxt="n", yaxt="n", xlab="", ylab="");

	barplot(qv_ci[,2], names.arg=sample_names, las=2, col="white",
		ylab="Median QV",
		main="Median Quality Values",
		add=T,
		cex.names=label.cex
	);
	points(mids, qv_ci[,1], pch="+", cex=2, col="red");
	points(mids, qv_ci[,3], pch="+", cex=2, col="green");

}

plot_len=function(len_ci, sample_names, label.cex=1){
	max_len=max(len_ci);

	num_samples=nrow(len_ci);
	mids=barplot(len_ci[,2], plot=F);
	spacing=mids[2]-mids[1];
	spacing=ifelse(is.na(spacing), 2, spacing);

	plot(0,0, type="n", 
		bty="l",
		xlim=c(0, (num_samples+.25)*spacing), ylim=c(0, max_len*1.2),
		xaxt="n", yaxt="n", xlab="", ylab="");

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
		xlim=range(pts$x)*c(1.2, 1.2), ylim=range(pts$y)*c(1/1.2, 1.2), cex=.5, col="red");
	qqline(values, col="blue");
	text(pts$x, pts$y, sample_names, col="black", cex=.5, pos=1);
}

###############################################################################

pdf(paste(OutputFilenameRoot, ".qv.pdf", sep=""), height=11, width=8.5);

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

	data=read_fastq_file(filepaths[i], keep_qual=T, head=Head, quick=NoValidation);
	qv_int=qvstr_to_int_list(data[["qual"]], Offset);

	summary=analyze_qv_list(qv_int, data[["max_len"]], sample_names[i]);

	num_reads_per_sample[i]=data[["num_recs"]];
	qv_ci[i,]=summary[["qv_ci"]];
	len_ci[i,]=summary[["len_ci"]];
	
	cat("\n");
}

#print(num_reads_per_sample);
#print(qv_ci)
#print(len_ci);

par(mfrow=c(3,1));
par(oma=c(0,0,3,0));
def_mar=par()$mar;
fat_bottom=def_mar;

fat_bottom[1]=fat_bottom[1]*max(1, max_samp_name_len/7);
if(fat_bottom[1]>10){
	label.cex=(cex=10/fat_bottom[1]);
	fat_bottom[1]=10;
}else{
	label.cex=1;
}

par(mar=fat_bottom);

plot_reads_per_sample(num_reads_per_sample, sample_names, label.cex=label.cex);
plot_qv(qv_ci, sample_names, label.cex=label.cex);
plot_len(len_ci, sample_names, label.cex=label.cex);

# Generate QQ plot for median read lengths
par(mfrow=c(3,1));
plot_qnorm(num_reads_per_sample, sample_names, "Reads/Sample");
plot_qnorm(qv_ci[,2], sample_names, "Read QVs");
plot_qnorm(len_ci[,2], sample_names, "Read Lengths");

mtext("Across All Samples", side=3, outer=T, cex=2);

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
