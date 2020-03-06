#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"directory", "d", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

LOG_FILE_NAME="fastq_qc_log";
LOG_FILE_EXT="tsv";
CTL_REGEX="^00[0-9][0-9]\\.";
POSTQC_DIR="PostQC";

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-d <FASTQ QC Pipeline Output Directory>\n",
	"\n",
	"This script should be run after the QC pipeline has\n",
	"completed.\n",
	"\n",
	"The files: ", LOG_FILE_NAME, ".*.", LOG_FILE_EXT, "\n",
	"will be concatenated together.\n",
	"\n",
	"Some plots looking at how sequence statistics have changed\n",
	"after each QC step will be generated, for both forward\n",
	"and reverse reads, separately.\n",
	"\n",
	"Paired merged statistics will also be generated and output\n",
	"of sample IDs that underperformed according to various cutoffs\n",
	"will also be generated.\n",
	"\n",
	"Output will go into the specific -d directory, under a new\n",
	"directory named, ", POSTQC_DIR, "\n",
	"\n",
	"The control names, specified with the regex: ", CTL_REGEX, "\n",
	"will be split out before the analyses begin.\n",
	"\n",
	"\n");

if(
	!length(opt$directory)
){
	cat(usage);
	q(status=-1);
}

###############################################################################

plot_text=function(strings){
        par(family="Courier");
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
}

###############################################################################

QC_LogDir=opt$directory;
post_qc_res_dir=paste(QC_LogDir, "/PostQC", sep="");

cat("\n")
cat("QC Result Directory: ", QC_LogDir, "\n");
cat("Post QC Result Directory: ", post_qc_res_dir, "\n");

dir.create(post_qc_res_dir);

###############################################################################

# Find list of files in directory
pattern=paste(LOG_FILE_NAME, "\\.\\d+\\.", LOG_FILE_EXT, sep="");
cat("Log Pattern: ", pattern, "\n");
qc_log_files_arr=sort(list.files(QC_LogDir, pattern=pattern));
qc_log_fullpath_arr=paste(QC_LogDir, "/", qc_log_files_arr, sep="");
print(qc_log_fullpath_arr);

cat("\n");
cat("Found these log files:\n");
print(qc_log_files_arr);

num_log_files=length(qc_log_files_arr);

cat("\nNumber of log files: ", num_log_files, "\n", sep="");

# Parse files names
name_parse=strsplit(qc_log_files_arr, "\\.");
indices_found=numeric();
for(i in 1:num_log_files){
	indices_found[i]=as.numeric(name_parse[[i]][2]);
}
indices_found=sort(indices_found);

cat("Indices Found:\n");
print(indices_found);

# Confirm all files in range exist
if(all(indices_found==(1:num_log_files))){
	cat("Indices found across expected range. (1 - ", num_log_files, ")\n", sep="");
}else{
	cat("Error:  Unexpected differences between found and expected log file indices.\n");
}

###############################################################################

full_file_fname=paste(post_qc_res_dir, "/",LOG_FILE_NAME, "._ALL_.", LOG_FILE_EXT, sep="");
cat("Concatenating all files into: ", full_file_fname, "\n", sep="");

full_tab=data.frame();

for(i in 1:num_log_files){
	cur_file=qc_log_files_arr[i];
	cat("Loading: ", cur_file, "\n"); 
	tab=read.table(qc_log_fullpath_arr[i], header=T, sep="\t", comment.char="");
	
	cnam=colnames(tab);
	cnam[1:3]=c("SampleID", "Direction", "QCStep");
	colnames(tab)=cnam;

	full_tab=rbind(full_tab, tab);
}

uniq_samp_ids=sort(unique(as.character(full_tab[,"SampleID"])));
num_samp_ids=length(uniq_samp_ids);

cat("Number of Samples Read:\n", num_samp_ids, "\n");

write.table(full_tab, file=full_file_fname, sep="\t", quote=F, row.names=F);

###############################################################################

cat("Splitting out controls based on: ", CTL_REGEX, "\n");
ctl_ix=grep(CTL_REGEX, full_tab[,"SampleID"]);

ctl_tab=full_tab[ctl_ix,,drop=F];
exp_tab=full_tab[-ctl_ix,,drop=F];

ctl_ids=as.character(ctl_tab[,"SampleID"]);
exp_ids=as.character(exp_tab[,"SampleID"]);

uniq_ctl_ids=sort(unique(ctl_ids));
uniq_exp_ids=sort(unique(exp_ids));

num_uniq_ctl_ids=length(uniq_ctl_ids);
num_uniq_exp_ids=length(uniq_exp_ids);

cat("Control IDs:\n");
print(uniq_ctl_ids);
cat("\n");
cat("Experimental IDs:\n");
print(uniq_exp_ids);
cat("\n");

# Extract out first sample's steps as a basis
first_exp=uniq_exp_ids[1];
first_entry=exp_tab[first_exp==exp_ids,];

forw_steps=as.character(first_entry[first_entry[,"Direction"]=="F","QCStep"]);
reve_steps=as.character(first_entry[first_entry[,"Direction"]=="R","QCStep"]);

rev_reads_found=T;
if(length(reve_steps)==0){
	rev_reads_found=F;	
}

if(length(forw_steps)==length(reve_steps) && !all(forw_steps==reve_steps)){
	cat("Error: Forward and Reverse steps don't match.\n");
	cat("Forward:\n");
	print(forw_steps);
	cat("Reverse:\n");
	print(reve_steps);
	cat("\n");
	quit(-1);
}

qc_steps=forw_steps;

cat("Step Order Identified:\n");
print(qc_steps);

plot_step_histograms=function(table, directions, steps, target_stats){

	num_steps=length(steps);

	maxs=apply(table[,target_stats,drop=F], 2, function(x){max(x, na.rm=T)} );
	cat("Maxes:\n");
	print(maxs);

	ac=function(x){
		return(prettyNum(x, big.mark=",", scientific=F));
	}

	for(stat_ix in target_stats){

		if(length(intersect(stat_ix, c("NumRecords", "NumBases")))){
			disp_stat=paste("log10(", stat_ix, ")",sep="");
			log_trans=T;
		}else{
			disp_stat=paste(stat_ix, sep="");
			log_trans=F;
		}

		cat("Stat:", stat_ix, "\n", sep="");

		if(num_steps>1){
			par(mfcol=c(num_steps,length(directions)));
		}else{
			par(mfcol=c(length(directions), 1));
		}

		par(mar=c(3,4,4,1));
		par(oma=c(0,0,3,0));
	
		for(cur_dir in directions){

			prev_med=NA;
			first_med=NA;

			for(step_ix in steps){
				cat("\tStep:", step_ix, "(", cur_dir, ")\n", sep="");

				srows=(table[,"QCStep"]==step_ix & table[,"Direction"]==cur_dir);
				values=table[srows, stat_ix];
				med_val=median(values, na.rm=T);

				if(is.na(first_med)){
					first_med=med_val;	
				}
				if(is.na(prev_med)){
					prev_med=med_val;
				}

				perc_tot=paste(round(100*med_val/first_med, 2), "%", sep="");
				perc_prev=paste(round(100*med_val/prev_med, 2), "%", sep="");

				if(log_trans){
					hist(log10(values+1), 
						breaks=seq(0, log10(maxs[stat_ix])*1.025, length.out=40), 
						xlab="", main=""); 
					abline(v=median(log10(values+1)), col="blue");
				}else{
					hist(values, 
						breaks=seq(0, maxs[stat_ix]*1.025, length.out=40), 
						xlab="", main="");
					abline(v=median(values, na.rm=T), col="blue");
				}
				title(main=paste(cur_dir, ": ", step_ix, ", median = ", ac(med_val), sep=""),
					cex.main=1);
				title(main=paste("Of All: ", perc_tot, "   Of Previous Step: ", perc_prev, sep=""), 
					font=3, line=1, cex.main=.8);

				prev_med=med_val;
			}
		}
		mtext(disp_stat, side=3, outer=T, cex=2, font=2);

	}	
}


pdf(paste(post_qc_res_dir, "/", LOG_FILE_NAME, ".summary.pdf", sep=""), height=11, width=8.5);

partial_tab=exp_tab[,
	c("SampleID", "Direction", "QCStep", "NumRecords", "NumBases", "LB95Len", "LB95QV")];

target_stats=c("NumRecords", "NumBases", "LB95Len", "LB95QV");

plot_text(c(
	"QC Run Post-Analysis:",
	"",
	paste("Run on: ", date(), sep=""),
	paste("Run by: ", system("whoami", intern=T), sep=""),
	"",
	paste("Command: "),
	paste("\t", commandArgs()),
	"",
	paste("Working Directory: ", getwd()),
	paste("Target Directory: ", QC_LogDir, sep=""),
	paste("Output Directory: ", post_qc_res_dir, sep=""),
	"",
	paste("Num Log Files Detected: ", num_log_files, sep=""),
	paste("Num Sample IDs Detected: ", num_samp_ids, sep=""),
	paste("Control Regular Expression: ", CTL_REGEX, sep=""),
	paste("Num Controls Detected: ", num_uniq_ctl_ids, sep=""),
	paste("Num Experimentals Detected: ", num_uniq_exp_ids, sep=""),
	"",
	paste("Num QC Steps Detected: ", length(qc_steps)),
	paste(capture.output(print(qc_steps, quote=F)))
));

###############################################################################
	
par(mfrow=c(1,1));
plot(0, xlim=c(-1,1), ylim=c(-1,1), type="n", bty="n", xaxt="n", yaxt="n", main="", xlab="", ylab="");
text(0,0, "Forward vs. Reverse Reads", font=2, cex=3);

plot_step_histograms(partial_tab, c("F", "R"), qc_steps, target_stats);

par(mfrow=c(1,1));
plot(0, xlim=c(-1,1), ylim=c(-1,1), type="n", bty="n", xaxt="n", yaxt="n", main="", xlab="", ylab="");
text(0,0, "Pairability Results", font=2, cex=3);
plot_step_histograms(partial_tab, c("for_frag", "rev_frag", "paired"), "merged", target_stats);

###############################################################################

paired_merged_ix=((partial_tab[,"Direction"]=="paired") & (partial_tab[,"QCStep"]=="merged"));
paired_merged_num_records=partial_tab[paired_merged_ix, "NumRecords"];
paired_merged_samp_id=partial_tab[paired_merged_ix, "SampleID"];
names(paired_merged_num_records)=paired_merged_samp_id;

plot_bar_cutoffs=function(values, cutoffs){
	
	num_cutoffs=length(cutoffs);

	num_exceeding=numeric();
	num_under=numeric();
	for(i in 1:num_cutoffs){
		num_exceeding[i]=sum(values>=cutoffs[i])
		num_under[i]=sum(values<cutoffs[i]);
	}

	max_counts=length(values);;	

	par(mfrow=c(2,1));
	par(oma=c(0,0,0,0));
	mids=barplot(num_exceeding, names.arg=paste(">=", cutoffs), 
		ylim=c(0, max_counts*1.1),
		main="Samples (NumRecords/2) Exceeding Depth Cutoffs",
		xlab="Paired Cutoffs", ylab="Num Samples", cex.names=.8);
	perc_excd=round(100*num_exceeding/max_counts, 2);
	text(mids, num_exceeding, paste(num_exceeding, "\n(", perc_excd, "%)", sep=""), font=3, cex=.7, pos=3);
	

	mids=barplot(num_under, names.arg=paste("<", cutoffs), 
		ylim=c(0, max_counts*1.1),
		main="Samples (NumRecords/2) Below Depth Cutoffs",
		xlab="Paired Cutoffs", ylab="Num Samples", cex.names=.8);
	perc_under=round(100*num_under/max_counts, 2);
	text(mids, num_under, paste(num_under, "\n(", perc_under, "%)", sep=""), font=3, cex=.7, pos=3);

}


num_pairable_sequences=paired_merged_num_records/2;

targeted_cutoffs=c(750, 1000, 2000, 3000, 10000, 20000, 50000, 10^ceiling(log10(max(num_pairable_sequences))));
plot_bar_cutoffs(num_pairable_sequences, targeted_cutoffs);

###############################################################################

export_redos=function(fname_root, counts_wnames, cutoffs){

	num_cutoff=length(cutoffs);
	cutoffs=sort(cutoffs);
	maxcutoff=max(cutoffs[is.finite(cutoffs)]);
	
	numdigits=floor(log10(maxcutoff))+1;
	spf_str=paste("%0", numdigits, "g", sep="");
	
	for(cutix in cutoffs){

		if(is.infinite(cutix)){
			fname=paste(fname_root, ".paird_cnts.ALL", ".tsv", sep="");
		}else{
			fname=paste(fname_root, ".paird_cnts.lt_", sprintf(spf_str, cutix), ".tsv", sep="");
		}

		below=as.matrix(counts_wnames[counts_wnames<cutix]);
		cat("Exporting to: ", fname, "\n");
		cat("   Num Samples: ", length(below), "\n");
		write.table(below, fname, quote=F, sep="\t", row.names=T, col.names=F);
	}

}

cat("\n");

export_redos(paste(post_qc_res_dir, "/", LOG_FILE_NAME, sep=""), 
	sort(num_pairable_sequences), c(750,1000,2000,3000,Inf));

dev.off();

###############################################################################

cat("Tar'ing Output...\n");
system(paste("tar -C ", QC_LogDir, " -cvzf ", QC_LogDir, "/QC.tgz PostQC", sep=""));

###############################################################################

cat("Done.\n")
print(warnings());

q(status=0)
