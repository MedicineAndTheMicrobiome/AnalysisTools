#!/usr/bin/env Rscript

###############################################################################
#                                                                             #
#       Copyright (c) 2014 J. Craig Venter Institute.                         #
#       All rights reserved.                                                  #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                             #
###############################################################################
###############################################################################

library('getopt');
library(sn);		# Skew normal libraries
library(VGAM);		# Additiona skewed normal libraries

params=c(
	"profA", "a", 1, "character",	# Transriptomics
	"profB", "b", 1, "character",	# Genomics
	"output_root", "o", 1, "character",
	"nameA", "A", 2, "character",
	"nameB", "B", 2, "character",
	"idmap", "m", 2, "character",
	"housekeeping", "h", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-a <input a (eg. transcriptomics) summary_table.tsv>\n",
	"	[-b <input b (eg. metagenomics) summary_table.tsv>]\n",
	"	-o <output filename root>\n",
	"	-h <housekeeping gene list file>\n",
	"	[-A <Prof A's name, default=Transcriptomics>]\n",
	"	[-B <Prof B's name, default=Metagenomics>]\n",
	"	[-m <ID mapping file>]\n",
	"\n",	
	"	If the metagenomics summary table (-b) is not specified, then over expression\n",
	"	will have to be based on relative abundance only.  If the metagenomics\n",
	"	is specified, then the analysis will be based on comparing the log odds ratio.\n",
	"	If only transcriptomics is specified, the analysis will use the log odds.\n",
	"	In both cases, over expression is computed relative to the housekeeping genes.",
	"\n",
	"	The ID mapping file should contain a mapping of the accessions\n",
	"	used in the summary_table to a more meaningful description.\n",
	"\n",
	"	The housekeeping gene list should be described with the same\n",
	"	accessions as in the summary_table, of course.  So if your\n",
	"	summary_table is full of GO ID's, then your housekeeping gene list should\n",
	"	consist of GO ID's and your ID mapping file should map from\n",
	"	GO ID to description.\n",
	"\n");

if(!length(opt$profA) || !length(opt$output_root) || !length(opt$housekeeping)){
	cat(usage);
	q(status=-1);
}

###############################################################################

ProfAFname=opt$profA;
ProfBFname=opt$profB;

GenomicsAvailable=TRUE;
if(length(ProfBFname)==0){
	GenomicsAvailable=FALSE;	
}

OutputRoot=opt$output_root;
nameA=opt$nameA;
nameB=opt$nameB;
MapFname=opt$idmap;
HousekeepingFname=opt$housekeeping;


cat("\n")
cat("Prof A: ", ProfAFname, "\n");

if(GenomicsAvailable){
	cat("Prof B: ", ProfBFname, "\n");
}

cat("Housekeeping Filename: ", HousekeepingFname, "\n");
cat("Output Filename Root: ", OutputRoot, "\n");
cat("\n");

if(length(MapFname)){
	cat("Using ID Mapping from: ", MapFname, "\n", sep="");
}


if(length(nameA)==0){
	nameA="Transcriptomics";	
}

if(length(nameB)==0){
	nameB="Metagenomics";
}

pdf(paste(OutputRoot, ".pdf", sep=""), height=11, width=8.5);
par(oma=c(0,0,1,0));
MAR_DEF=par()$mar;

cat("\n");

###############################################################################
###############################################################################

load_idmap=function(MapFname){
	
	cat("Loading ID Mapping: ", MapFname, "\n");

	mat=as.matrix(read.delim(MapFname, header=F, sep="\t", na.strings=""));
	nrows=nrow(mat);

	cat("Num rows: ", nrows, "\n");

	mat[,1]=gsub("^ +", "", mat[,1]);
	mat[,1]=gsub(" +$", "", mat[,1]);

	mapping=list(nrows);

	for(i in 1:nrows){
		id=mat[i,1];
		name=mat[i,2];
		mapping[[id]]=name;
	}
	#print(head(mapping));
	cat("Done loading ID Map.\n");
	return(mapping);
}

###############################################################################

load_list=function(ListFname){
	cat("Loading list: ", ListFname, "\n");
	mat=as.matrix(read.delim(ListFname, header=F, sep="\t", na.strings=""));
	cat("Num items loaded: ", nrow(mat), "\n");
	return(mat[,1]);
}

###############################################################################

rename=function(names, mapping){
	num_names=length(names);
	mapnames=names(mapping);
	
	new_names=character(num_names);
	for(i in 1:num_names){
		idx=which(mapnames==names[i]);
		#cat("idx: ", idx, "\n");
		if(length(idx)!=0){
			new_names[i]=mapping[[names[i]]];
		}else{
			new_names[i]=names[i];
		}
	}
	#print(new_names);
	return(new_names);

}

###############################################################################

load_summary_file=function(InputFileName){

	cat("Loading: ", InputFileName, "\n");

	# Load data
	inmat=as.matrix(read.delim(InputFileName, sep="\t", header=TRUE, check.names=FALSE,
		na.strings="", comment.char="", quote="", row.names=1))

	#cat("Original Matrix:\n")
	#print(inmat);

	# Grab columns we need into a vector, ignore totals, we won't trust it.
	counts_mat=inmat[,2:(ncol(inmat)), drop=FALSE];
	#print(counts_mat);

	rownames(counts_mat)=rownames(inmat);
	colnames(counts_mat)=colnames(inmat)[2:(ncol(inmat))];
	return(counts_mat);
}

###############################################################################

normalize=function(counts_mat){
	sample_totals=apply(counts_mat, 1, sum);
	normalized=matrix(0, nrow=nrow(counts_mat), ncol=ncol(counts_mat));
	for(i in 1:nrow(counts_mat)){
		normalized[i,]=counts_mat[i,]/sample_totals[i];
	}
	colnames(normalized)=colnames(counts_mat);
	rownames(normalized)=rownames(counts_mat);
	return(normalized);	
}

###############################################################################

compute_logratio=function(trans_norm_mat, genom_norm_mat){
	ncol=ncol(trans_norm_mat);
	nrow=nrow(trans_norm_mat);
	lr_mat=log(trans_norm_mat/genom_norm_mat);
	colnames(lr_mat)=colnames(trans_norm_mat);
	rownames(lr_mat)=rownames(trans_norm_mat);
	return(lr_mat);
}

###############################################################################

estimate_pvalue=function(x, empirical_distribution){
	num_pts=length(empirical_distribution);
	pval=1-sum(empirical_distribution<x)/num_pts;
	return(pval);
}

###############################################################################
###############################################################################

if(length(MapFname)>0){
	idmap=load_idmap(MapFname);
	cat("\n");
}


###############################################################################
# Load summary files

prof_a_count_mat=load_summary_file(ProfAFname);

num_a_samples=nrow(prof_a_count_mat);
num_a_categories=ncol(prof_a_count_mat);

cat("A: Num Samples:", num_a_samples, "\n");
cat("A: Num Categories: ", num_a_categories, "\n");
cat("\n");

if(GenomicsAvailable){
	prof_b_count_mat=load_summary_file(ProfBFname);

	num_b_samples=nrow(prof_b_count_mat);
	num_b_categories=ncol(prof_b_count_mat);
	cat("B: Num Samples:", num_b_samples, "\n");
	cat("B: Num Categories: ", num_b_categories, "\n");
	cat("\n");
}else{
	prof_b_count_mat=NULL;
}

###############################################################################
# Normalize, identify 0's, keep track of minimum non-zero abundance

prof_a_norm_mat=normalize(prof_a_count_mat)[1,];
nonzero_a=prof_a_norm_mat>0;
num_nonzero_a=sum(nonzero_a);
cat("Num non zero genes in A: ", num_nonzero_a, "\n");
min_a_prop=min(prof_a_norm_mat[nonzero_a]);
set_min_a_prop=min_a_prop/10;
genes_in_a=names(prof_a_norm_mat[nonzero_a]);

num_genes_in_a=length(genes_in_a);
prof_a_norm_mat=prof_a_norm_mat[nonzero_a];

cat("Min in A: ", min_a_prop, "\n");
cat("Reset A's 0's to: ", set_min_a_prop, "\n");
cat("Num peptides in A: ", num_genes_in_a, "\n", sep="");
cat("\n");

if(GenomicsAvailable){
	prof_b_norm_mat=normalize(prof_b_count_mat)[1,];
	nonzero_b=prof_b_norm_mat>0;
	num_nonzero_b=sum(nonzero_b);
	cat("Num non zero genes in B: ", num_nonzero_b, "\n");
	min_b_prop=min(prof_b_norm_mat[nonzero_b]);
	set_min_b_prop=min_b_prop/10;
	genes_in_b=names(prof_b_norm_mat[nonzero_b]);
	num_genes_in_b=length(genes_in_b);
	prof_b_norm_mat=prof_b_norm_mat[nonzero_b];

	cat("Min in B: ", min_b_prop, "\n");
	cat("Reset B's 0's to: ", set_min_b_prop, "\n");
	cat("Num peptides in B: ", num_genes_in_b, "\n", sep="");

	genes_in_b=names(prof_b_norm_mat[nonzero_b, drop=F]);
	genes_in_ab=intersect(genes_in_a, genes_in_b);
	num_genes_in_b=length(genes_in_b);
	num_genes_in_ab=length(genes_in_ab);
	cat("Num peptides in A & B: ", num_genes_in_ab, "\n", sep="");
}

cat("\n");

###############################################################################
# Set the abundances for genomics if transcriptomics is available.

num_genes_in_a=length(genes_in_a);

if(GenomicsAvailable){
		
	# Independent of if any of the genes are housekeeping,
	# If the genes are in the transriptomics, they must be
	# in the genomics, so set the min genomics for any
	# transcriptomic gene to the 1/10 min nonzero genomics

	# Trans		Genomics	Status?
	# 0		0		no evidence for gene 
	# 0		1		not expressed, though capability exists
	# 1		0		over expressed, so gene must be there but not enough sequencing depth
	# 1		1		log ratio could be computed.

	cat("Setting the gene abundances for no hit genomics, if there were transcriptomic hits.\n");
	
	# Pad out transcriptomics with 0, if there is metagenomics evidence.
	for(gene in genes_in_b){
		if(is.na(prof_a_norm_mat[gene])){
			prof_a_norm_mat[gene]=0;
		}
	}

	# Update metagenomics, if there is transcriptomic evidence
	for(gene in genes_in_a){
		if(is.na(prof_b_norm_mat[gene]) || prof_b_norm_mat[gene]==0){
			prof_b_norm_mat[gene]=set_min_b_prop;
		}
	}

}

###############################################################################
###############################################################################
# Load housekeeping genes and set minimums for these

housekeeping_list=load_list(HousekeepingFname);

hk_in_a=intersect(housekeeping_list, genes_in_a);
num_hk_in_a=length(hk_in_a);
cat("Housekeeping genes in A: ", num_hk_in_a, "\n");

# For housekeeping genes, it is assumed that they must be present since they
# should be.  So automatically set these genes for transcriptomics and genomics
# to 1/10 the mimimum non zero for each respectively.

# Trans		Genomics	Status?
# 0		0		Set Tr and Gen, to 1/10 min
# 0		1		Underexpressed, set Tr to 1/10 min
# 1		0		Overexpressed, set Gen to 1/10 min
# 1		1		log ratio could be computed.

for(hk_gene in housekeeping_list){
	if(is.na(prof_a_norm_mat[hk_gene]) || prof_a_norm_mat[hk_gene]==0){
		prof_a_norm_mat[hk_gene]=set_min_a_prop;
	}
}

if(GenomicsAvailable){

	hk_in_b=intersect(housekeeping_list, genes_in_b);
	num_hk_in_b=length(hk_in_b);
	cat("Housekeeping genes in B: ", num_hk_in_b, "\n");

	for(hk_gene in housekeeping_list){
		if(is.na(prof_b_norm_mat[hk_gene]) || prof_b_norm_mat[hk_gene]==0){
			prof_b_norm_mat[hk_gene]=set_min_b_prop;
		}
	}
}

###############################################################################

all_gene_names=sort(names(prof_a_norm_mat));

prof_a_norm_mat=prof_a_norm_mat[all_gene_names];

odda=(prof_a_norm_mat/(1-prof_a_norm_mat))

if(GenomicsAvailable){
	profAlen=length(prof_a_norm_mat);
	profBlen=length(prof_b_norm_mat);
	if(profAlen!=profBlen){
		cat("Note that the length of each profile is not equal:\n");
		cat("Prof A Length: ", profAlen, "\n");
		cat("Prof B Length: ", profBlen, "\n");
	}
	
	prof_b_norm_mat=prof_b_norm_mat[all_gene_names];
	oddb=(prof_b_norm_mat/(1-prof_b_norm_mat))

	stat=log(odda/oddb,10);

	stat_label=paste("Log10[odds(", nameA,")/odds(", nameB, ")]", sep="");
}else{
	stat=log(odda,10);
	stat_label=paste("Log10[odds(", nameA, ")]");

	#stat=log(prof_a_norm_mat,10);
	#stat_label=paste("Log10(Pr(", nameA, "))");
}

# Remove negative infinities (ie. 0 expression)
notNegInf=(stat!=-Inf);
stat=stat[notNegInf];

print(stat);
cat("Labeling x-axis as: ", stat_label, "\n");

all_genes=names(stat);
housekeeping_stats=stat[housekeeping_list];
remaining_stats=stat[setdiff(all_genes, housekeeping_list)];

###############################################################################

###############################################################################

###############################################################################
###############################################################################

compute_MLE=function(stats, title, xlabel, color, x_range){
	
	cat("Working on: ", title, "\n");

	x_fit_pts=seq(x_range[1], x_range[2], length.out=100);

	# Define negative log likelihoood object function
	nLL=function(loc, scal, shape){
		if(scal <= 0){
			return(1e305);
		}
		nll=(-sum(dsn(stats, loc, scal, shape, log=T)));
		#cat(sprintf("%f: %3.2f, %3.2f, %3.4f\n", nll, loc, scal, shape));
		return(nll);
	}

	# Use MLE to estimate parameters
	fit_loglik=-1e305;

	# Try multiple starting points to avoid local minima
	cat("Selecting starting points:\n");
	loc_range=seq(-10, 10, length.out=3);
	scale_range=seq(.001, 10, length.out=3);
	shape_range=seq(-50, 50, length.out=3);

	best_coef=list();
	for(loc_start in loc_range){
		for(scal_start in scale_range){
			for(shape_start in shape_range){

				cat(sprintf("Trying: loc=%3.2f, scale=%3.2f, shape=%3.2f ",
					 loc_start, scal_start, shape_start));

				tc_res=tryCatch({
					skew_fit=mle(nLL, start=list(loc=loc_start, scal=scal_start, shape=shape_start));
				}, error=function(e){
					cat("x\n");
				});

				if(!is.null(tc_res)){
					cat(".\n");
					ll=logLik(tc_res);
					if(ll>fit_loglik){
						best_coef=attr(tc_res, "coef");
						fit_loglik=ll;
						cat("\nBest LogLik so far: ", ll, "\n");
						print(best_coef);
						cat("\n");
					}
				}
			}
		}
	}

	cat("Lowest LogLike = ", fit_loglik, "\n");
	fit_coef=c(best_coef["loc"], best_coef["scal"], best_coef["shape"]);
	print(best_coef);

	# Plot gene distribution
	num_breaks=nclass.Sturges(stats)*2;
	h=hist(stats, breaks=num_breaks, plot=F);
	hist(stats, freq=F,  xlim=x_range, ylim=c(0, max(h$density)*1.2),
		breaks=num_breaks,
		xlab=stat_label, col="grey",
		main=paste("Distribution for ", title, sep=""),
		axes=F
	);
	
	stat_range=range(stats);
	axis(side=1, las=1, at=seq(floor(stat_range[1]), ceiling(stat_range[2]), 1));
	y_axis_limit=max(h$density)*1.2;
	axis(side=2, las=2, at=seq(0, y_axis_limit, round(y_axis_limit/5, digits=1)));

	# Get pdf for fit
	y_fit_pts=dsn(x_fit_pts,
		location=fit_coef[1],
		scale=fit_coef[2],
		shape=fit_coef[3]
	);

	# Overlay fitted pdf
	points(x_fit_pts, y_fit_pts, col=color, type="b", cex=.3);

	legend(x=par()$usr[1], y=par()$usr[4], 
		legend=c("Observed", "Fitted"),
		fill=c("black", color), bty="n",
		cex=.7);

	# Report fitted skewed normal parameters
	text(x=par()$usr[1], y=par()$usr[4]/2, labels=
		sprintf("MLE-based Skewed Normal Parameters:\n  Location = %3.2f\n  Scale = %3.2f\n  Shape = %3.2f",
		  best_coef[1], best_coef[2], best_coef[3]),
		pos=4,
		cex=.7
	);

	# Report number of observations used in fit
	text(x=par()$usr[1], y=par()$usr[4]*3/4, labels=
		sprintf("Num Observations: %i", length(stats)),
		pos=4,
		cex=.7
	);

	# Return values in a list
	return_list=list();
	return_list[["fit_coef"]]=fit_coef;
	return_list[["mle_xdens"]]=x_fit_pts;
	return_list[["mle_ydens"]]=y_fit_pts;
	return(return_list);
}

###############################################################################
###############################################################################
###############################################################################

xranges=range(stat);
cat("Ranges for statistic: ", xranges[1], " - ", xranges[2], "\n");

delta=xranges[2]-xranges[1];
xranges=c(xranges[1]-delta*.2, xranges[2]+delta*.2);
par(mfrow=c(3,1));

housekeeping_fit=compute_MLE(housekeeping_stats, "Housekeeping", stat_label, "red", xranges);
remaining_fit=compute_MLE(remaining_stats, "Remaining", stat_label, "blue", xranges);

# Write filenames in output
mtext(paste("Transcriptomics: ", ProfAFname), outer=T, col="blue", cex=.5);
if(GenomicsAvailable){
	mtext(paste("Genomics: ", ProfBFname), outer=T, col="blue", cex=.5, line=-1);
}

###############################################################################
# Plot Overlay of all and HK log ratios distributions

yranges=range(c(housekeeping_fit$mle_ydens, remaining_fit$mle_ydens));

plot(0,0, type="n", 
	xlim=xranges,
	ylim=yranges,
	main="Overlay of Remaining and Housekeeping Distributions",
	xlab=stat_label,
	ylab="Density"
);

points(housekeeping_fit$mle_xdens, housekeeping_fit$mle_ydens, type="b", col="red", cex=.3);
points(remaining_fit$mle_xdens, remaining_fit$mle_ydens, type="b", col="blue", cex=.3);

legend(x=par()$usr[1], y=par()$usr[4], 
	legend=c("Fitted Housekeeping", "Fitted Remaining"),
	fill=c("red", "blue"), bty="n",
	cex=.7);

###############################################################################
# Generate bar plot for the housekeeping genes

sorted_housekeeping_stats=sort(housekeeping_stats, decreasing=T);
housekeeping_stat_range=range(sorted_housekeeping_stats);

plot_range=housekeeping_stat_range;
if(plot_range[2]<0){
	plot_range[2]=0;
}
axis_ticks=seq(floor(plot_range[1]), ceiling(plot_range[2]), .5);

# Generate bar plot
par(mfrow=c(2,1));
par(mar=c(15, 4.1, 4.1, 5.1));
barmids=barplot(sorted_housekeeping_stats, border=NA, axisnames=F,
	main=paste("Sorted Distribution for Housekeeping Genes", sep=""),
	ylim=range(axis_ticks),
	ylab=stat_label,
	xlab="",
	yaxt="n"
);

# Place axis labels on right hand side

axis(side=2, at=axis_ticks, labels=axis_ticks, las=2, cex.axis=1);

if(length(MapFname)>0){
	shkablr_ids=names(sorted_housekeeping_stats);
	shkablr_names=rename(shkablr_ids, idmap);
	names(sorted_housekeeping_stats)=shkablr_names;
}

# Label genes underneath bar plot
num_hk_genes=length(sorted_housekeeping_stats);
gene_desc=names(sorted_housekeeping_stats);
spacing=mean(diff(barmids));
label_size=min(c(1, 24/num_hk_genes));

for(i in 1:num_hk_genes){
	text(barmids[i]-1*spacing/2, -.1, gene_desc[i], srt=-45, xpd=T, cex=label_size, pos=4);
}

###############################################################################
# Plot the distribution of p-values across the positive expression genes

cat("Plotting p-values across all genes.\n");

sorted_remaining_stats=sort(remaining_stats, decreasing=T);

pvalues=1-psn(sorted_remaining_stats, 
	housekeeping_fit$fit_coef[1],
	housekeeping_fit$fit_coef[2],
	housekeeping_fit$fit_coef[3]
);

if(GenomicsAvailable){
	positive_idx=sorted_remaining_stats>0;
}else{
	positive_idx=rep(T, length(sorted_remaining_stats));
}

positive_pvals=pvalues[positive_idx];
all_positive=sorted_remaining_stats[positive_idx];

all_positive_colors=rep("blue", length(all_positive));
all_positive_colors[positive_pvals<.1]="green";
all_positive_colors[positive_pvals<.05]="yellow";
all_positive_colors[positive_pvals<.01]="orange";
all_positive_colors[positive_pvals<.001]="red";

bar_mids=barplot(all_positive, col=all_positive_colors, border=NA, axisnames=F,
	main="Sorted Positively Expressed Genes",
	ylab=stat_label,
	xlab="Genes"
);

max_stat=floor(max(all_positive));
min_stat=ceiling(min(all_positive));
x_marks=numeric();
marker_cnt=1;
for(i in min_stat:max_stat){
	x_marks[marker_cnt]=max(which(all_positive>=i));
	marker_cnt=marker_cnt+1;
}

axis(side=1, at=bar_mids[x_marks], label=min_stat:max_stat, cex.axis=.75, col="grey40", col.axis="grey40");

abline(v=bar_mids[x_marks], col="grey50", lwd=.5);

legend(x=par()$usr[2], y=par()$usr[4], 
	xjust=1, # Right Justified
	title="P-values Color Encoding:",
	legend=c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "Not Significant"),
	fill=c("red", "orange", "yellow", "green", "blue"), 
	bty="n",
	cex=.7);

###############################################################################
# Plot the genes with significant p-values labeled

cat("Plotting genes with significant overexpression p-values.\n");

if(length(MapFname)>0){
	aplr_ids=names(all_positive);
	aplr_names=rename(aplr_ids, idmap);
	names(all_positive)=aplr_names;
}
	
num_significant_genes=sum(pvalues<0.1);
cat("Number of significant genes: ", num_significant_genes, "\n", sep="");

par(mfrow=c(2,1));
par(mar=c(20, 1.1, 4.1, 7.1));
genes_per_plot=50;
num_plots=ceiling(num_significant_genes/genes_per_plot);
max_pos=max(all_positive);
min_pos=min(all_positive);

if(min_pos<0 && max_pos<0){
	max_pos=0;
}

for(i in 1:num_plots){
	start_ix=(i-1)*genes_per_plot + 1;
	end_ix=i*genes_per_plot;

	cat("Plotting Top Significantly Expressed Genes from ", start_ix, " to ", end_ix, ".\n", sep="");
	
	# Generate bar plot
	
	
	barmids=barplot(
		rep(0, end_ix-start_ix+1), 
		border=NA, 
		axisnames=F,
		main=sprintf("Sorted Positively Expressed Genes: %i to %i", start_ix, end_ix),
		ylab=stat_label,
		xlab="",
		ylim=c(min_pos, max_pos),
		yaxt="n",
	);

	abline(h=0:max_pos, col="grey40", lty="dotted");

	barmids=barplot(all_positive[start_ix:end_ix], 
		col=all_positive_colors[start_ix:end_ix],
		axisnames=F,
		add=T,
		yaxt="n"
	);

#	barmids=barplot(all_positive[start_ix:end_ix], 
#		col=all_positive_colors[start_ix:end_ix],
#		border=NA, 
#		axisnames=F,
#		main=sprintf("Sorted Positively Expressed Genes: %i to %i", start_ix, end_ix),
#		ylab=stat_label,
#		xlab="",
#		ylim=c(min_pos, max_pos),
#		yaxt="n"
#	);

	# Place axis labels on right hand side
	axis_ticks=seq(min_pos, max_pos, .5);
	axis(side=4, at=axis_ticks, labels=sprintf("%3.1f", axis_ticks), las=2, cex.axis= .5);

	# Label genes underneath bar plot
	gene_desc=names(all_positive[start_ix:end_ix]);
	spacing=mean(diff(barmids));
	for(i in 1:genes_per_plot){
                text(barmids[i]-1.2*spacing, -.1, gene_desc[i], srt=-45, xpd=T, cex=32/genes_per_plot, pos=4);
        }

	# Place legend
	legend(x=par()$usr[2], y=par()$usr[4], 
		xjust=1, # Right Justified
		title="P-values Color Encoding:",
		legend=c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "Not Significant"),
		fill=c("red", "orange", "yellow", "green", "blue"), 
		bty="n",
		cex=.5);

}


###############################################################################

dev.off();
cat("Done.\n")
if(!is.null(warnings())){
	print(warnings());
}
q(status=0)
