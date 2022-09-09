#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
library(stringr);

options(useFancyQuotes=F);
options(digits=5);
options(width=200);

params=c(
	"targets", "t", 1, "character",
	"count_tab", "c", 1, "character",
	"outputroot", "o", 1, "character",
	"abund_cutoff", "A", 2, "numeric",
	"prevl_cutoff", "P", 2, "numeric",
	"fasta_file", "f", 2, "character",
	"mothur_bin", "m", 2, "character",
	"fa_extract_bin", "e", 2, "character",
	"ref_taxa_align", "n", 2, "character",
	"ref_taxa_tax", "x", 2, "character"
);

ABUND_CUTOFF=0.001;
PREVL_CUTOFF=0.15;

MOTHUR_DEF_BIN="/usr/bin/mothur_1.44.1/mothur/mothur";
FASTA_EXTR_BIN="~/git/AnalysisTools/FASTA/Extract_Record_From_FASTA_By_List.pl";
REF_ALGN="/home/kelvinli/git/AnalysisTools/16S_Clust_Gen_Pipeline/silva_reference/silva.nr_v138.align";
REF_TAX="/home/kelvinli/git/AnalysisTools/16S_Clust_Gen_Pipeline/silva_reference/silva.nr_v138.tax";

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-t <Target Sequence List Filename>\n",
	"	-c <Count Table Filename>\n",
	"	-o <Output Filename Root>\n",	
	"\n",
	"	[-A <abundance cutoff, default=", ABUND_CUTOFF, ">]\n",
	"	[-P <prevalence cutoff, default=", PREVL_CUTOFF, ">]\n",
	"\n",
	"	Classification paramaters\n",
	"	[-f <fasta of representative sequences>]\n",
	"	[-m <mothur bin, \n\t\tdefault=", MOTHUR_DEF_BIN, ">]\n",
	"	[-e <extraction script, \n\t\tdefault=", FASTA_EXTR_BIN, ">]\n",
	"	[-n <reference taxa .align file, \n\t\tdefault=", REF_ALGN, ">]\n",
	"	[-x <reference taxa .tax file, \n\t\tdefault=", REF_TAX, ">]\n",
	"\n",
	"This script will evaluate the prevalence and abundance of each\n",
	"unique sequence.\n",
	"\n",
	"Within the Mothur pipeline the files should be applied\n",
	"to the following parameters:\n",
	" Target Sequence List: *precluster.denovo.vsearch.accnos\n",
	" Count Table:          *precluster.count_table\n",
	" FASTA File:		*precluster.fasta\n",
	"\n");

if(
	!length(opt$targets) || 
	!length(opt$count_tab) ||
	!length(opt$outputroot)
){
	cat(usage);
	q(status=-1);
}

TargetsFile=opt$targets
CountsFile=opt$count_tab;
OutputFnameRoot=opt$outputroot;

AbundanceCutoff=ABUND_CUTOFF;
if(length(opt$abund_cutoff)){
	AbundanceCutoff=opt$abund_cutoff;
}

PrevalenceCutoff=PREVL_CUTOFF;
if(length(opt$prevl_cutoff)){
	PrevalenceCutoff=opt$prevl_cutoff;
}

MothurBin=MOTHUR_DEF_BIN;
FastaExtrBin=FASTA_EXTR_BIN;
RefAlign=REF_ALGN;
RefTax=REF_TAX;
FastaFile="";


if(length(opt$fasta_file)){
	FastaFile=opt$fasta_file;
}
if(length(opt$mothur_bin)){
	MothurBin=opt$mothur_bin;
}
if(length(opt$fa_extract_bin)){
	FastaExtrBin=opt$fa_extract_bin;
}
if(length(opt$ref_taxa_align)){
	RefAlign=opt$ref_taxa_align;
}
if(length(opt$ref_taxa_tax)){
	RefTax=opt$ref_taxa_tax;
}

cat("\n");
cat("Targets Filename: ", TargetsFile, "\n", sep="");
cat("Counts Filename: ", CountsFile, "\n", sep="");
cat("Output Filename Root: ", OutputFnameRoot, "\n", sep="");
cat("\n");
cat("Abundance Cutoff: ", AbundanceCutoff, "\n", sep="");
cat("Prevalence Cutoff: ", PrevalenceCutoff, "\n", sep="");
cat("\n");
cat("FASTA Filename: ", FastaFile, "\n", sep="");
cat("Mothur Bin Path: ", MothurBin, "\n", sep="");
cat("FASTA Extr Bin Path: ", FastaExtrBin, "\n", sep="");
cat("\n");

##############################################################################

load_count_table_prev_abund_only=function(fname){

	num_lines=as.numeric(system(paste("wc -l ", fname, " | cut -f 1 -d ' '"), intern=T, wait=T));
	cat("Num lines in : ", fname, ": ", num_lines, "\n");
	num_repres=num_lines-3;

	# Read in first few lines, then table
	fh=file(fname, "r");
	comments=readLines(fh,n=1);			# Skip comments
	sampleoffsets_example=readLines(fh,n=1);	# Skip examples

	# Read in order/offset of each sample id being represented
	sampleoffset_line=readLines(fh,n=1);
	sampleoffsets=strsplit(sampleoffset_line, "\\t")[[1]];
	sample_ids=sampleoffsets[3:length(sampleoffsets)];	# First 2 fields are headers

	# Accumulate to return
	rep_to_sample_map=vector(mode="list", length=num_repres);
	
	cat("Loading Rep to Sample Mapping...\n");
	line=readLines(fh,n=1);
	rep_ix=1;
	seq_rep_ids=character(num_repres);
	cents=floor(num_lines/20);
	while(length(line)>0){
	
		# Split the lines into representative/num_reps/represented,counts,etc.
		lineinfo=strsplit(line, "\\t")[[1]];

		representative=lineinfo[1];
		num_represented=lineinfo[2];
		source_ids=lineinfo[3:length(lineinfo)];

		rep_to_sample_map[[rep_ix]]=source_ids;
		seq_rep_ids[rep_ix]=representative;
		line=readLines(fh,n=1);

		# Report progress
		rep_ix=rep_ix+1
		if(!(rep_ix%% cents)){
			cat("Loaded: ", rep_ix, " lines (", round(rep_ix/num_repres*100,2), "%).\n", sep="");
		}
	}
	close(fh);
	cat("done.\n");

	names(rep_to_sample_map)=seq_rep_ids;

	cat("Coverting Rep to Sample Map to Prevalence/Abundance Table...\n");
	# Convert to count and sample
	num_reps=length(seq_rep_ids);
	
	cn=c("NumSamples", "NumReads", "Prevalence", "Abundance");
	stat_mat=matrix(NA, nrow=num_reps, ncol=length(cn));
	colnames(stat_mat)=cn;
	rownames(stat_mat)=seq_rep_ids;

	cat("Num Sequence Representatives: ", num_reps, "\n");
	status_ix=0;
	for(sr_ix in seq_rep_ids){

		#cat("\n\nExamining: ", sr_ix, "\n");

		samp_info=rep_to_sample_map[[sr_ix]];
		#print(samp_info);
		split_list=strsplit(samp_info, ",");
		#print(split_list);

		num_samples=length(samp_info);
		#cat("  Num Samples: ", num_samples, "\n");

		rep_count=0;
		for(i in 1:num_samples){
			samp_id=split_list[[i]][1];
			counts=as.numeric(split_list[[i]][2]);

			rep_count=rep_count+counts;
		}

		stat_mat[sr_ix, "NumReads"]=rep_count;
		stat_mat[sr_ix, "NumSamples"]=num_samples;

		# Report progress
		status_ix=status_ix+1;
		if(!(status_ix%%10000)){
			cat("Processing: ", status_ix, " records (", round(100*status_ix/num_reps,2), "%)\n", sep="");
		}

	}
	cat("done.\n");

	sample_ids=unique(sample_ids);
	num_samples=length(sample_ids);
	cat("Num Samples Identified: ", num_samples, "\n");

	num_reads=sum(stat_mat[,"NumReads"]);
	cat("Num Reads Identified: ", num_reads, " represented.\n");

	stat_mat[,"Prevalence"]=stat_mat[,"NumSamples"]/num_samples;
	stat_mat[,"Abundance"]=stat_mat[,"NumReads"]/num_reads;

	return(stat_mat);
}

##############################################################################

plot_text=function(strings, max_lines=75){

	plot_page=function(strings){
		orig_par=par(no.readonly=T);

		par(mfrow=c(1,1));
		par(family="Courier");
		par(oma=rep(.5,4));
		par(mar=rep(0,4));

		num_lines=length(strings);

		top=max_lines;

		plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
			xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
			);
		for(i in 1:num_lines){
			#cat(strings[i], "\n", sep="");
			text(0, top-i, strings[i], pos=4, cex=.8);
		}

		par(orig_par);
	}

	num_lines=length(strings);
	num_pages=ceiling(num_lines / max_lines);
	#cat("Num Pages: ", num_pages, "\n");
	for(page_ix in 1:num_pages){
		start=(page_ix-1)*max_lines+1;
		end=start+max_lines-1;
		end=min(end, num_lines);
		##print(c(start,end));
		plot_page(strings[start:end]);
	}
}

##############################################################################

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, counts=F){

	orig_par=par(no.readonly=T);
	par(mfrow=c(1,1));

        num_row=nrow(mat);
        num_col=ncol(mat);

        cat("Num Rows: ", num_row, "\n");
        cat("Num Cols: ", num_col, "\n");

        mat=mat[rev(1:num_row),, drop=F];

        num_colors=50;
        color_arr=rainbow(num_colors, start=0, end=4/6);
        if(high_is_hot){
                color_arr=rev(color_arr);
        }

        remap=function(in_val, in_range, out_range){
                in_prop=(in_val-in_range[1])/(in_range[2]-in_range[1])
                out_val=in_prop*(out_range[2]-out_range[1])+out_range[1];
                return(out_val);
        }

        if(is.na(plot_min)){
                plot_min=min(mat);
        }
        if(is.na(plot_max)){
                plot_max=max(mat);
        }
        cat("Plot min/max: ", plot_min, "/", plot_max, "\n");
	par(mar=c(10,10,1,1));
        plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), 
		xaxt="n", yaxt="n", bty="n", xlab="", ylab="", main=title);

        # x-axis
        axis(side=1, at=seq(.5, num_col-.5, 1), labels=colnames(mat), las=2, cex.axis=.7);
        axis(side=2, at=seq(.5, num_row-.5, 1), labels=rownames(mat), las=2, cex.axis=.7);

        if(log_col){
                plot_min=log10(plot_min+.0125);
                plot_max=log10(plot_max+.0125);
        }

	text_size=min(1, 17.5/num_row);
	cat("Text Size: ", text_size, "\n");

        for(x in 1:num_col){
                for(y in 1:num_row){

                        if(log_col){
                                col_val=log10(mat[y,x]+.0125);
                        }else{
                                col_val=mat[y,x];
                        }

                        remap_val=remap(col_val, c(plot_min, plot_max), c(1, num_colors));
                        col_ix=ceiling(remap_val);

                        rect(x-1, y-1, (x-1)+1, (y-1)+1, border=NA, col=color_arr[col_ix]);

                        if(counts){
                                text_lab=sprintf("%i", mat[y,x]);
                        }else{
                                text_lab=mat[y,x];
                        }
                        text(x-.5, y-.5, text_lab, srt=45, cex=text_size, font=2);
                }
        }

	par(orig_par);

}

##############################################################################

load_list=function(filename){
        val=scan(filename, what=character(), comment.char="#");
        return(val);
}

##############################################################################
# Main Program Starts Here!
##############################################################################

pdf(paste(OutputFnameRoot, ".seq_check.pdf", sep=""), height=11, width=10);

param_msg=capture.output({
	cat("Targets Filename:\n");
	cat("  ", TargetsFile, "\n");
	cat("\n");
	cat("Counts Filename:\n");
	cat("  ", CountsFile, "\n");
	cat("\n");
	cat("Output Filename Root: ", OutputFnameRoot, "\n");
	cat("\n");
	cat("Abundance Cutoff: ", AbundanceCutoff, "\n");
	cat("Prevalence Cutoff: ", PrevalenceCutoff, "\n");
});

plot_text(param_msg);

# Load target list
cat("Loading Targets File...\n");
targets=load_list(TargetsFile);
num_targets=length(targets);
cat("Number of targets loaded: ", num_targets, "\n");

# Load Counts
cat("Loading Counts File...\n");
stat_mat=load_count_table_prev_abund_only(CountsFile);

##############################################################################
# Sort by prev/abund

prev_sort_ix=order(stat_mat[,"Prevalence"], decreasing=T);
temp_stat_mat_by_prev=stat_mat[prev_sort_ix,];
abund_sort_ix=order(temp_stat_mat_by_prev[,"Abundance"], decreasing=T, method="shell");
stat_matrix_by_abundance=temp_stat_mat_by_prev[abund_sort_ix,];

abun_sort_ix=order(stat_mat[,"Abundance"], decreasing=T);
temp_stat_mat_by_abun=stat_mat[abun_sort_ix,];
prev_sort_ix=order(temp_stat_mat_by_abun[,"Prevalence"], decreasing=T, method="shell");
stat_matrix_by_prevalence=temp_stat_mat_by_abun[prev_sort_ix,];

max_head=min(20, nrow(stat_mat));

table_heads_msg=capture.output({
	cat("Top ", max_head, " Representative Sequences: \n\n", sep="");
	cat("Sorted by Abundance:\n");
	print(stat_matrix_by_abundance[1:max_head,]);
	cat("\n\n");
	cat("Sorted by Prevalence:\n");
	print(stat_matrix_by_prevalence[1:max_head,]);
});

print(table_heads_msg);
plot_text(table_heads_msg);

abund_range=range(stat_mat[,"Abundance"]);
log_min_abund=log10(abund_range[1]);
log_max_abund=log10(abund_range[2]);

##############################################################################

num_rep_sequences=nrow(stat_mat);
cat("Num Total Representatives: ", num_rep_sequences, "\n");

CUTOFF=0.95;
abund_cutoff_ix=sum(cumsum(stat_matrix_by_abundance[,"Abundance"]));
cat("Num Representives to go to cutoff:", abund_cutoff_ix, "\n");

NUM_PLOT=1500;
top_targ_num=min(NUM_PLOT, num_rep_sequences);
top_reps=unique(c(
	rownames(stat_matrix_by_abundance[1:top_targ_num,]),
	rownames(stat_matrix_by_prevalence[1:top_targ_num,])
));

num_top_reps=length(top_reps);
cat("Number of Top Reps to Plot: ", num_top_reps, "\n");
plot(
	log10(stat_matrix_by_abundance[top_reps, "Abundance"]), 
	stat_matrix_by_prevalence[top_reps, "Prevalence"],
	xlab="Log10(Abundance)",
	ylab="Prevalence",
	xlim=c(log_min_abund, log_max_abund),
	ylim=c(0,1),
	main=paste("Top ", NUM_PLOT, " Representative Sequnces", sep="")
);

abline(v=log10(AbundanceCutoff), col="green", lty="dashed");
abline(h=PrevalenceCutoff, col="blue", lty="dashed");


# Count singletons
num_singletons=sum(stat_mat[,"NumReads"]==1);
num_sequences_represented=sum(stat_mat[,"NumReads"]);
cat("Number of Reads Represented: ", num_sequences_represented, "\n");

##############################################################################

targets_by_abund=intersect(rownames(stat_matrix_by_abundance), targets);
targ_stat_mat=stat_matrix_by_abundance[targets_by_abund,];

max_head=min(20, nrow(targ_stat_mat));
#print(targ_stat_mat[1:max_head,]);

cat("Num Targets: ", length(targets_by_abund), "\n");
head(targets_by_abund)

# Mark targets on scatterplot, don't plot singletons, there are too many
non_sing_ix=(targ_stat_mat[,"NumReads"]>1);
num_singletons_target=sum(!non_sing_ix);
min_prev=min(targ_stat_mat[,"Prevalence"]);

top_targ=intersect(top_reps, targets);

points(log10(targ_stat_mat[top_targ,"Abundance"]), 
	targ_stat_mat[top_targ,"Prevalence"], pch=4, cex=1.2, col="red");

if(num_singletons_target){
	points(log_min_abund, min_prev, pch=4, col="red", cex=1.5);
	text(log_min_abund, min_prev, 
		paste(num_singletons_target, " singletons", sep=""),
		srt=90, adj=c(-.5,.5),
		cex=.7, font=2
		); 
}

num_targets_represented=sum(targ_stat_mat[,"NumReads"]);

##############################################################################
# Extract targets that exceeded cutoff

extract_targets_ix=(targ_stat_mat[,"Abundance"]>=AbundanceCutoff) & 
	(targ_stat_mat[,"Prevalence"]>=PrevalenceCutoff);
extr_targ_stat_mat=targ_stat_mat[extract_targets_ix,,drop=F];
extr_targ_seq_ids=rownames(extr_targ_stat_mat);

msg=c(
	"The following sequence IDs were targeted, but they exceed both",
	"the abundance and prevalance cutoffs:",
	"",
	paste("Abundance: ", AbundanceCutoff),
	paste("Prevelance: ", PrevalenceCutoff),
	"",
	"Sample IDs:",
	capture.output(print(extr_targ_stat_mat, quote="F"))
);

print(msg);
plot_text(msg);

# Export extracted targets
output_extr_lst_fn=paste(OutputFnameRoot, ".extr_targets.lst", sep="");
fh=file(output_extr_lst_fn, "w");
cat(file=fh, paste(extr_targ_seq_ids, collapse="\n"), "\n", sep="");
close(fh);

# Export all seq above cutoffs
extract_all_above_ix=
	(stat_mat[,"Abundance"]>=AbundanceCutoff) & 
	(stat_mat[,"Prevalence"]>=PrevalenceCutoff);
all_above_stat_mat=stat_mat[extract_all_above_ix,];
all_above_seq_ids=rownames(all_above_stat_mat);
cat("All Above Sequence IDs:\n");
print(all_above_seq_ids);
output_above_lst_fn=paste(OutputFnameRoot, ".above_cutoff.lst", sep="");
fh=file(output_above_lst_fn, "w");
cat(file=fh, paste(all_above_seq_ids, collapse="\n"), "\n", sep="");
close(fh);

##############################################################################

msg=c(
	"Note: The 'representative' sequences are a set of unique sequences",
	"that have been selected to represent the underlying sequences in each sample.",
	"A representative sequence that only represents a single underlying sequence",
	"is called a 'singleton'.  Singletons could be very low abundance taxa",
	"or artifacts due to sequencing error. The 'represented number of sequences'",
	"is the number of redundant sequences the representive sequences represents.",
	"This value is used to estimate the (expanded) abundance of that sequence across",
	"all the samples.",
	"",
	"OVERALL:",
	paste("Number Representative Sequences (unique): ", num_rep_sequences, " read from count table", sep=""),
	paste("Number Sequences Represented (expanded): ", num_sequences_represented, sep=""),
	"",
	paste("Number of Singletons: ", num_singletons, "/", num_rep_sequences,
		" (", round(num_singletons/num_rep_sequences*100,2), "% of representatives are singletons)", sep=""),
	paste("Cumulative Singleton Contribution: ", num_singletons, "/", num_sequences_represented, 
		" (", round(num_singletons/num_sequences_represented*100,2), 
		"% of abundance are singletons)", sep=""),
	"",
	"",
	"TARGETED:",
	paste("Number Representatives Targets Loaded: ", num_targets, sep=""),
	paste("Representatives Targeted: ", num_targets, "/", num_rep_sequences,
		" (", round(num_targets/num_rep_sequences*100, 2), "% of all representatives)", sep=""),
	paste("Number of Represented (abundance) Targeted: ", 
		num_targets_represented, "/", num_sequences_represented,
		" (", round(num_targets_represented/num_sequences_represented*100, 2), 
		"% of all abundance)", sep=""),
	"",
	paste("Number of Singletons Targeted: ", num_singletons_target),
	paste("Percentage of All Singletons Targeted: ",
		num_singletons_target, "/", num_singletons,
		"(", round(num_singletons_target/num_singletons*100, 2), "% of Singleton representatives)", sep=""), 
	paste("Cumulative Targeted Singleton Contribution: ",
		num_singletons_target, "/", num_sequences_represented, 
		"(", round(num_singletons_target/num_sequences_represented*100,2), "% of abundance)", sep=""), 

	""
);

print(msg);
plot_text(msg);

##############################################################################

par(oma=c(0,0,2,0));
par(mfrow=c(2,1));

# Abundance
abund_breaks=seq(log_min_abund, 0, length.out=20);
hdall=hist(log10(stat_mat[,"Abundance"]), breaks=abund_breaks, 
	xlab="Log10(Abundance)", main="All Reads");
hdtar=hist(log10(targ_stat_mat[,"Abundance"]), breaks=abund_breaks, 
	xlab="Log10(Abundance)", main="Targeted");
mtext("Comparison of Log10(Abundance) Distributions", outer=T, cex=1.5, font=2);

# Log scale abundance
log_hdall=hdall;
log_hdtar=hdtar;
log_hdall$counts=log10(log_hdall$counts+1);
log_hdtar$counts=log10(log_hdtar$counts+1);
plot(log_hdall, ylab="Log10(Frequency)", main="All Reads", xlab="Log10(Abundance)");
plot(log_hdtar, ylab="Log10(Frequency)", main="Targeted", xlab="Log10(Abundance)");
mtext("Log Scale Comparison of Log10(Abundance) Distributions", outer=T, cex=1.5, font=2);

# Prevalence
prevl_breaks=seq(0,1,length.out=20);
hist(stat_mat[,"Prevalence"], breaks=prevl_breaks, main="All Reads",
	xlab="Prevalence");
hist(targ_stat_mat[,"Prevalence"], breaks=prevl_breaks, main="Targeted",
	xlab="Prevalence");
mtext("Comparison of Prevalence Distributions", outer=T, cex=1.5, font=2);

# Log scale Prevalence
max_prev=max(stat_mat[,"Prevalence"]);
prevl_breaks=seq(0,max_prev,length.out=20);
hdall=hist(stat_mat[,"Prevalence"], breaks=prevl_breaks, plot=F);
hdtar=hist(targ_stat_mat[,"Prevalence"], breaks=prevl_breaks, plot=F);
log_hdall=hdall;
log_hdtar=hdtar;
log_hdall$counts=log10(log_hdall$counts+1);
log_hdtar$counts=log10(log_hdtar$counts+1);
plot(log_hdall, ylab="Log10(Frequency)", main="All Reads", xlab="Log10(Abundance)");
plot(log_hdtar, ylab="Log10(Frequency)", main="Targeted", xlab="Log10(Abundance)");
mtext("Log Scale Comparison of Prevalence Distributions", outer=T, cex=1.5, font=2);

##############################################################################

if(FastaFile != ""){

	# seq ids are from: all_above_stat_mat
	cat("All Seq IDs above cutoffs matrix:\n");
	print(all_above_stat_mat);
	cat("Targeted Seqs above cutoff:\n");
	print(extr_targ_seq_ids);
	all_above_seq_ids=rownames(all_above_stat_mat);

	# Extract Fasta file
	extracted_fasta_filename=paste(OutputFnameRoot, ".above_cutoff.fasta", sep="");
	cmd=paste(
		FastaExtrBin, 
			" -f ", FastaFile, 
			" -l ", output_above_lst_fn, 
			" > ", extracted_fasta_filename,
		sep="");

	cat(cmd, "\n");
	system(cmd, wait=T);

	# Run mothur's RDP classifier on it
	cmd=paste(
		MothurBin,
			"\"#classify.seqs(",
				"fasta=", extracted_fasta_filename, ", ",
				"template=", RefAlign, ", ",
				"taxonomy=", RefTax, ", ",
				"cutoff=50", ", ",
				"processors=16",
			")\"", sep=" ");

	cat(cmd, "\n");
	system(cmd, wait=T);

	# Get extension used in mothur's output
	split_res=strsplit(RefTax, "\\.")[[1]];
	num_tok=length(split_res);
	vers_ext=split_res[num_tok-1];
	class_result_file=paste(OutputFnameRoot, ".above_cutoff.", vers_ext, ".wang.taxonomy", sep="");
	class_data=read.table(class_result_file, stringsAsFactors=F);

	# Merge classification results with stat matrix
	colnames(class_data)=c("TargetSeqID", "FullTaxonomy");
	print(class_data);
	num_class=nrow(class_data);
	genus=character(num_class);
	tsid=character(num_class);

	for(i in 1:num_class){
		tsid[i]=class_data[i, "TargetSeqID"];
		genus[i]=strsplit(class_data[i,"FullTaxonomy"], ";")[[1]][6];
	}
	names(genus)=tsid;


	# Format output
	Genus=rep("", nrow(all_above_stat_mat));
	all_above_stat_mat_char=cbind(all_above_stat_mat, Genus);
	rownames(all_above_stat_mat_char)=rownames(all_above_stat_mat);

	all_above_ids=rownames(all_above_stat_mat);
	all_above_stat_mat_char[all_above_ids, "Genus"]=genus[all_above_ids];

	all_above_stat_mat_char[,"Prevalence"]=
		sprintf("%3.3f", as.numeric(all_above_stat_mat[,"Prevalence"]));
	all_above_stat_mat_char[,"Abundance"]=
		sprintf("%6.6f", as.numeric(all_above_stat_mat[,"Abundance"]));

	# Mark the rows that were targeted
	Targets=rep("", nrow(all_above_stat_mat));
	all_above_stat_mat_char=cbind(all_above_stat_mat_char, Targets);
	all_above_stat_mat_char[extr_targ_seq_ids, "Targets"]="*";

	# Output stats
	msg=c(
		"Classifications of Targeted Sequences:",
		"",
		capture.output(print(all_above_stat_mat_char, quote=F))
	);
	
	print(msg);
	plot_text(msg);
	
	#######################################################################
	# Regenerate scatter plot and label points with classifications. 

	par(mfrow=c(1,1));
	plot(0,0, type="n", xlab="Log10(Abundance)", ylab="Prevalence",
		xlim=c(log10(AbundanceCutoff),0), 
		ylim=c(PrevalenceCutoff, max(stat_mat[,"Prevalence"])*1.1),
		main="Sequences Above Cutoffs"
	);

	# Draw lines where the cutoffs are
	abline(v=log10(AbundanceCutoff), col="green", lty="dashed");
	abline(h=PrevalenceCutoff, col="blue", lty="dashed");
	
	# Plot all above cutoff
	points(
		log10(all_above_stat_mat[,"Abundance"]),
		all_above_stat_mat[,"Prevalence"]);

	# Mark targets above cutoff
	points(
		log10(all_above_stat_mat[extr_targ_seq_ids,"Abundance"]),
		all_above_stat_mat[extr_targ_seq_ids,"Prevalence"],
		pch=4, col="red", cex=1.2
		);
	
	# Label all the classifications
	for(i in 1:num_class){
		x=log10(as.numeric(all_above_stat_mat[i,"Abundance"]));
		y=as.numeric(all_above_stat_mat[i,"Prevalence"]);
		label=all_above_stat_mat_char[i,"Genus"];	
		if(any(all_above_seq_ids[i]==extr_targ_seq_ids)){
			labcol="black";
			labcex=1;
		}else{
			labcol="grey40";
			labcex=.8;
		}
		text(x, y, label, pos=3, font=2, cex=labcex, col=labcol);
	}
	

}

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
