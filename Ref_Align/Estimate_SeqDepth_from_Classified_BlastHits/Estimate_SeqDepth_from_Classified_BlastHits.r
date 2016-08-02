#!/usr/bin/env Rscript

###############################################################################
options(width=200);

library('getopt');

params=c(
	"blast_in", "b", 1, "character",
	"total_raw", "r", 1, "numeric",
	"total_blasted", "q", 1, "numeric",
	"output_root", "o", 1, "character",
	"acc_map", "m", "2", "character",
	"excl_only", "e", "2", "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"\n",
	"	-b <blast output, >=2 subjects, sorted by query/read ID>\n",
	"	-r <total raw fragments before QC>\n",
	"	-q <total good fragments blasted, divide by two if F/R included, i.e. interleaved paired fasta was used>\n",
	"	-o <output filename root>\n",
	"	[-m <accession to name map>]\n",
	"	[-e (flag, only report exclusive hits)]\n",
	"\n",
	"This script will read in a blast file that contains hits\n",
	"to more than one organism, sorted by read ID.\n",
	"\n",
	"It will then compute the number of reads necessary to uniquely identify\n",
	"a specific number of reads from that organism (to that subject ID) with a 95% confidence.\n",
	"\n",
	"Note that if the reads are paired, i.e., the fasta file contains\n",
	"forward and reverse reads from the same fragment, then for this script\n",
	"to work properly, the sequence ids must be the same.  If the forward\n",
	"and reverse reads have different sequence IDs, then they will be double\n",
	"counted.\n",
	"\n",
	"At low sequencing depths the relationship between sequencing and \n",
	"minimum acquired fragments from an organism will be non-linear. However,\n",
	"at greater sequencing depth the relationship becomes linear because even\n",
	"though the variance of the binomial distribution increases with trials,\n",
	"as a proportion of trials, the variance is very low.\n",
	"\n",
	"This script will generate the following output in the PDF file:\n",
	"	1.) Summary of set overlap counts and proportions relative to\n",
	"		total raw counts and total reads that were blasted.\n",
	"	2.) A curve for each of the set overlap possibilities across\n",
	"		various percent identity cutoffs.\n",
	"	3.) Rarefaction lines/curves (for each set overlap possibility)\n",
	"		predicting the number of reads necessary to detect\n",
	"		various levels of overlap counts.\n",
	"\n",
	"The -m option, specifies a map file so that the accessions used in\n",
	"blast input can be displayed with a more easy to remember\n",
	"sequence name.\n",
	"\n",
	"The -e option, specifies that only the sequences that are exclusively\n",
	"mapping will be displayed.\n",
	"\n", sep="");

if(
	!length(opt$blast_in) || 
	!length(opt$total_raw) || 
	!length(opt$total_blasted) ||
	!length(opt$output_root) 
){
	cat(usage);
	q(status=-1);
}

BlastInput=opt$blast_in;
TotalRaw=opt$total_raw;
TotalBlasted=opt$total_blasted;
OutputRoot=opt$output_root;
AccessionMap=opt$acc_map;
ExclusiveOnly=length(opt$excl_only)>0;

cat("\n");
cat("Blast input: ", BlastInput, "\n");
cat("Total number of raw: ", TotalRaw, "\n");
cat("Total number of blasted: ", TotalBlasted, "\n");
cat("Output Root: ", OutputRoot, "\n");
cat("\n");
cat("Accession Map: ", AccessionMap, "\n");
cat("Exclusive Only: ", ExclusiveOnly, "\n");


################################################################################

cutoffs=c(85, 87.5, 90, 92.5, 95, 97.5, 100);
num_cutoffs=length(cutoffs);
subj_comb=character();	# Subject Combinations

hits=list();

# Example:
# hits[["A"]]      = c(3, 3, 2, 1, 1);
# hits[["A & B"]]  = c(3, 3, 3, 2, 2);

filter=function(lines, cutoff){
	values=as.integer(lines[,3]);
	return(lines[values>=cutoff,, drop=F]);
}

#------------------------------------------------------------------------------

process_buffer=function(lines, hits){

	#print(lines);

	cutoff=character();
	for(i in 1:num_cutoffs){
		lines_abv_cutoff=filter(lines, cutoffs[i]);	
		subjects=sort(unique(lines_abv_cutoff[,2]));

		if(length(subjects)>0){
			subjs_str=paste(subjects, collapse=" & ");

			# If new combination found, initialize new row
			if(!any(names(hits)==subjs_str)){
				hits[[subjs_str]]=rep(0, num_cutoffs);
				subj_comb=c(subj_comb, subjs_str);
			}

			# Add counts to the combination
			tmp=hits[[subjs_str]];
			tmp[i]=tmp[i]+1;
			hits[[subjs_str]]=tmp;

		}
	}

	#print(hits);

	return(hits);

}

#------------------------------------------------------------------------------

read_and_process_blast=function(blast_in, hits){

	fh=file(blast_in, open="r");

	lines=character();
	prev="";
	buf="";

	hits_procs=0;
	qrys_procs=0;

	eof=F;
	while(!eof){

		buf=readLines(con=fh, n=1);

		if(length(buf)>0){

			comp=unlist(strsplit(buf, "\t"));

			# Process records if new query detected
			if(comp[1]!=prev && length(lines)>0){
				hits=process_buffer(lines, hits);
				qrys_procs=qrys_procs+1;
				lines=character();
			}

			# Add new line to buffer
			lines=rbind(lines,comp);
			prev=comp[1];
			hits_procs=hits_procs+1;

		}else{
			eof=T;
		}

	}


	# Process remaining rows
	if(nrow(lines)>0){
		hits=process_buffer(lines, hits);
		qrys_procs=qrys_procs+1;
	}

	close(fh);

	# Report processing statistics
	cat("Hits Processed: ", hits_procs, "\n");
	cat("Queries Processed: ", qrys_procs, "\n");

	return(hits);
}

################################################################################

plot_text=function(strings){

	def=par(no.readonly=T);
        par(family="Courier");
        par(oma=rep(.5,4));
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

	par(def);
}

################################################################################

plot_rarefaction=function(proportions, title, sbj, vbar){

	subject=rownames(proportions);
	cutoffs=as.numeric(colnames(proportions));
	num_cutoffs=length(cutoffs);

	cat("Working on: ", title, ": ", subject, "\n");

	max_x=1.2e7;
	increments=seq(0, max_x, 5e5);
	num_increments=length(increments);
	
	mat=matrix(0, nrow=num_cutoffs, ncol=num_increments);

	for(i in 1:num_cutoffs){
		min_hits_at95=qbinom(0.05, increments, proportions[i]);
		mat[i,]=min_hits_at95;
	}

	max_y=max(mat);
	par(mar=c(10, 4.1,4.1,2.1));
	plot(0, type="n", 
		ylim=c(0,max_y), xlim=c(0, max_x),
		xlab="", ylab="Minimum Expected Hits >95% Confidence",
		main=sbj,
		xaxt="n",
		yaxt="n"
	);

	# Draw axes
	# X
	x_ticks=seq(0, max_x, round(max_x/10));
	axis(side=1, at=x_ticks, label=prettyNum(x_ticks, big.mark=","), las=2);
	# Y
	y_ticks=seq(0, max_y, ceiling(max_y/10));
	if(length(y_ticks)<2){
		y_ticks=c(0,2);	
	}
	axis(side=2, at=y_ticks, label=y_ticks);

	# Mark observed/empirical depths
	abline(v=vbar, col="grey", lwd=2);

	for(i in 1:num_cutoffs){
		points(increments, mat[i,], col=i, type="b", pch=i, lty=i);
	}


	mtext(title, side=3, line=.5, cex=.7);

	legend(0, max(y_ticks), paste(">=", cutoffs, "% Identity"), pch=1:num_cutoffs, col=1:num_cutoffs);


}

################################################################################

plot_identity_curves=function(counts, title){


	num_curves=nrow(counts);
	num_ids=ncol(counts);
	hit_names=rownames(counts);
	idents=as.numeric(colnames(counts));

	def=par(no.readonly=T);

	par(mfrow=c(1,1));

	y_range=range(counts);

	plot(0, type="n", 
		xlim=range(idents), ylim=c(y_range[1], y_range[2]*1.2),
		xlab="Percent Identity", ylab="Number of Hits",
		main=title
	);

	for(i in 1:num_curves){
		points(idents, counts[i,], col=i, pch=i, type="b");
	}

	max_y=max(counts);
	legend(95, max_y*1.2, hit_names, pch=1:num_curves, col=1:num_curves);

	par(def);

}

################################################################################

load_map=function(mapfname){
	inmap=as.matrix(read.table(mapfname, sep="\t", header=FALSE));
	map=list();
	for(i in 1:nrow(inmap)){
		map[[inmap[i,1]]]=inmap[i,2];
	}
	return(map);
}

################################################################################

hits=read_and_process_blast(BlastInput, hits);

# Combine results into matrix
combinations=names(hits);
num_comb=length(combinations);

count_mat=matrix(0, nrow=num_comb, ncol=num_cutoffs, 
	dimnames=list(combinations, cutoffs));

for(i in 1:num_comb){
	count_mat[i,]=hits[[combinations[i]]];
}

if(ExclusiveOnly){
	cat("Removing overlapping hits...\n");
	hit_names=rownames(count_mat);
	non_exc=grep(" & ", hit_names);
	count_mat=count_mat[-non_exc,, drop=F];
}

if(length(AccessionMap)>0){
	map=load_map(AccessionMap);
	hit_names=rownames(count_mat);
	for(i in 1:num_comb){
		new_name=map[[hit_names[i]]];
		if(length(new_name)){
			hit_names[i]=new_name;
		}
	}
	rownames(count_mat)=hit_names;
	combinations=hit_names;
	num_comb=length(combinations);
}


print(count_mat);
# Sort counts by decreasing order
count_mat=count_mat[order(count_mat[,1], decreasing=T), , drop=F];

cat("Counts:\n");
print(count_mat);

blasted_read_prob=count_mat/TotalBlasted;
raw_read_prob=count_mat/TotalRaw;

blasted_read_ppm=round(blasted_read_prob*1e6);
raw_read_ppm=round(raw_read_prob*1e6);

cat("\n");
cat("Blasted PPM (parts per million):\n");
print(blasted_read_ppm);
cat("\n");
cat("Raw PPM (parts per million):\n");
print(raw_read_ppm);

cat("\n");
cat("Blasted Proportions:\n");
print(blasted_read_prob);
cat("\n");
cat("Raw Proportions:\n");
print(raw_read_prob);

################################################################################
# Output spreadsheets

fh=file(paste(OutputRoot, ".csv", sep=""), "w");

cat(file=fh, paste("Blast input:,", BlastInput), "\n");
cat(file=fh, paste("Output Root:,", OutputRoot), "\n");
cat(file=fh, "\n");

cat(file=fh, paste("Total number of raw:,", TotalRaw, "\n"));
cat(file=fh, paste("Total number of blasted:,", TotalBlasted, "\n"));
cat(file=fh, "\n");

cat(file=fh, "Total Hits:\n");
cat(file=fh, "Subject\\%Id,");
write.table(count_mat, quote=F, sep=",", file=fh);
cat(file=fh, "\n");

cat(file=fh, "Blasted PPM:\n");
cat(file=fh, "Subject\\%Id,");
write.table(blasted_read_ppm, quote=F, sep=",", file=fh);
cat(file=fh, "\n");

cat(file=fh, "Raw PPM:\n");
cat(file=fh, "Subject\\%Id,");
write.table(raw_read_ppm, quote=F, sep=",", file=fh);
cat(file=fh, "\n");

cat(file=fh, "Blasted Proportions:\n");
cat(file=fh, "Subject\\%Id,");
write.table(blasted_read_prob, quote=F, sep=",", file=fh);
cat(file=fh, "\n");

cat(file=fh, "Raw Proportions:\n");
cat(file=fh, "Subject\\%Id,");
write.table(raw_read_prob, quote=F, sep=",", file=fh);
cat(file=fh, "\n");

close(fh);


################################################################################

pdf(paste(OutputRoot, ".pdf", sep=""), height=8.5, width=11);

#-------------------------------------------------------------------------------

par(mfrow=c(1,1));

ct_txt=capture.output(print(count_mat));
bl_ppm_txt=capture.output(print(blasted_read_ppm));
bl_txt=capture.output(print(blasted_read_prob));
rw_ppm_txt=capture.output(print(raw_read_ppm));
rw_txt=capture.output(print(raw_read_prob));


plot_text(c(

	paste("Blast input: ", BlastInput),
	paste("Output Root: ", OutputRoot),
	"\n",
	paste("Total number of raw: ", prettyNum(TotalRaw, big.mark=",")),
	paste("Total number of blasted: ", prettyNum(TotalBlasted, big.mark=",")),
	"\n",
	"Blast Counts:",
	ct_txt,
	"\n",
	"Blast PPM:",
	bl_ppm_txt,
	"\n",
	"Raw PPM:",
	rw_ppm_txt,
	"\n"

));

plot_text(c(

	"Raw Proportions:", 
	rw_txt,
	"\n",
	"Blasted Proportions:", 
	bl_txt

));

#-------------------------------------------------------------------------------
# Plot percent identity curves

plot_identity_curves(count_mat, title=OutputRoot);

#-------------------------------------------------------------------------------
# Plot rarefaction for Blasted Reads

cat("\n");

par(mfrow=c(1,3));
for(i in 1:num_comb){
	plot_rarefaction(blasted_read_prob[i,, drop=F], title="Blasted/Good Reads", sbj=combinations[i], vbar=TotalBlasted);
}

#-------------------------------------------------------------------------------
# Plot rarefaction for Raw Reads

par(mfrow=c(1,3));
for(i in 1:num_comb){
	plot_rarefaction(raw_read_prob[i,, drop=F], title="Raw Reads", sbj=combinations[i], vbar=TotalRaw);
}

################################################################################

dev.off();

w=warnings();
if(length(w)){
	print(w);
}

cat("\nDone.\n");

q(status=0)
