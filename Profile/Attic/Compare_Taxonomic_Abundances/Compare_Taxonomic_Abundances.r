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

COMBO_TOP=100;
IND_TOP=50;
REMAINING="Remaining";

library('getopt');

params=c(
	"profA", "a", 1, "character",
	"profB", "b", 1, "character",
	"output_root", "o", 1, "character",
	"nameA", "A", 2, "character",
	"nameB", "B", 2, "character",
	"idmap", "m", 2, "character",
	"sum_up", "s", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-a <input a (eg. transcriptomics) summary_table.csv>\n",
	"	-b <input b (eg. genomiacs) summary_table.csv>\n",
	"	-o <output filename root>\n",
	"	[-A <Prof A's name>]\n",
	"	[-B <Prof B's name>]\n",
	"	[-m <ID mapping file>]\n",
	"	[-s (sum up across samples first)]\n",
	"\n",	
	"This script will take in two summary tables, A (eg. metatranscriptomic data)\n",
	"and B (eg. metagenomic data), and generate a several plots which compare their\n",
	"relative abundances.\n",
	"\n",
	"Accumulation of Shared Taxa: Plots how many taxa are shared between the two\n",
	"profiles as the number of top taxa in each profile are expanded.  For example,\n",
	"if at the top 40 (x-axis), the number of shared taxa are 7 (y-axis), then it\n",
	"means that when looking at the top 40 taxa of each profile, only 7 taxa are shared.\n",
	"\n",
	"Log ratio of Top Shared Taxa: A bar plot of the shared taxa from the two profiles.\n",
	"The order of the taxa in the plot are sorted by order in which they were accumulated\n",
	"as the top X taxa are expanded.\n",
	"\n",
	"Log ratio of Ranked Abundance: A plot of the ratio of the ranked abundances.\n",
	"This plots the relative evenness of the two profiles.  If two profiles have\n",
	"a similar distribution of abundances then the line will be flat along 0.\n",
	"Note that this plot does not match taxa.  It is purely based on comparing the\n",
	"ranked abundance curves relative shapes.\n",
	"\n",
	"Rank abundance curve:  A plot of the usual rank abundance curve.\n",
	"\n",
	"Log ratio of Prof A/Prof B sorted by A: A plot of the Log(A/B), sorted by the abundance of A.\n",
	"This plot provides the LR's from the perspective of A's top taxa.  Since these\n",
	"are based on the top taxa of A which is in the numerator, there should not be\n",
	"any negative infinities, and the differences will typically be statistically\n",
	"significant.  This avoids comparing low abundance vs low abundance which may\n",
	"produce a wide range of ratios.\n",
	"\n",
	"Log ratio of Prof B/Prof A sorted by B: similar to above graph, but from perspective\n",
	"of B.",
	"\n",
	"\n",
	"The number of top taxa compared: ", COMBO_TOP, "\n",
	"The number of top taxa graphed for individual graphs: ", IND_TOP, "\n",
	"The taxa name that will be removed from the graphs (but included in normalization): ", REMAINING, "\n",
	"\n");

if(!length(opt$profA) || !length(opt$profB) || !length(opt$output_root)){
	cat(usage);
	q(status=-1);
}

###############################################################################

ProfAFname=opt$profA;
ProfBFname=opt$profB;
OutputRoot=opt$output_root;

nameA=opt$nameA;
nameB=opt$nameB;

if(length(nameA)==0){
	nameA=tail(strsplit(ProfAFname, "/")[[1]], 1);
}
if(length(nameB)==0){
	nameB=tail(strsplit(ProfBFname, "/")[[1]], 1);
}


MapFname=opt$idmap;
SumUp=opt$sum_up;
if(length(SumUp)==0){
	SumUp=FALSE;
}else{
	SumUp=TRUE;
}

cat("\n")
cat("Prof A: ", ProfAFname, "\n");
cat("Prof B: ", ProfBFname, "\n");
cat("Output Filename Root: ", OutputRoot, "\n");
cat("Sum up across samples?: ", SumUp, "\n");
cat("\n");

if(length(MapFname)){
	cat("Using ID Mapping from: ", MapFname, "\n", sep="");
}

pdf(paste(OutputRoot, ".pdf", sep=""), height=11, width=8.5);
MAR_DEF=par()$mar;

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

prof_a_count_mat=load_summary_file(ProfAFname);
prof_b_count_mat=load_summary_file(ProfBFname);

if(SumUp){
	sum_a=t(apply(prof_a_count_mat, 2, sum));
	sum_b=t(apply(prof_b_count_mat, 2, sum));
	colnames(sum_a)=colnames(prof_a_count_mat);
	colnames(sum_b)=colnames(prof_b_count_mat);
	prof_a_count_mat=sum_a;
	prof_b_count_mat=sum_b;
	cat("NOTE: Column values have been summed up across samples.\n");
}else{
	if(nrow(prof_a_count_mat)>1 || nrow(prof_b_count_mat)>1){
		cat("WARNING: You have multiple samples in this summary table file.\n");
	}
}

if(length(MapFname)>0){
	cat("Performing remapping of IDs...\n")
	idmap=load_idmap(MapFname);
	colnames(prof_a_count_mat)=rename(colnames(prof_a_count_mat), idmap);
	colnames(prof_b_count_mat)=rename(colnames(prof_b_count_mat), idmap);
	cat("Ok.\n");

}

cat("\n");
num_a_samples=nrow(prof_a_count_mat);
num_a_categories=ncol(prof_a_count_mat);
cat("A: Num Samples:", num_a_samples, "\n");
cat("A: Num Categories: ", num_a_categories, "\n");
cat("\n");
num_b_samples=nrow(prof_b_count_mat);
num_b_categories=ncol(prof_b_count_mat);
cat("B: Num Samples:", num_b_samples, "\n");
cat("B: Num Categories: ", num_b_categories, "\n");
cat("\n");

if(num_a_samples!=num_b_samples){
	cat("Error, number of samples (rows) do not match.\n");
	quit(status=-1);
}else{
	num_samples=num_a_samples;
}

###############################################################################
# Normalize
prof_a_norm_mat=normalize(prof_a_count_mat);
prof_b_norm_mat=normalize(prof_b_count_mat);

###############################################################################
###############################################################################

plot_shared_top=function(profA, profB, num_top, nameA="A", nameB="B"){
	a_sorted=sort(profA, decreasing=T);
	b_sorted=sort(profB, decreasing=T);

	a_gtzero=a_sorted>0;
	b_gtzero=b_sorted>0;

	# Do not compare ranks if they have zero abundance
	a_sorted=a_sorted[a_gtzero];
	b_sorted=b_sorted[b_gtzero];

	top_A=names(a_sorted);
	top_B=names(b_sorted);

	lenA=length(top_A);
	lenB=length(top_B);
	
	minlen=min(c(lenA, lenB));

	if(num_top>minlen){
		num_top=minlen;
	}

	num_shared=numeric();
	for(i in 1:num_top){
		memA=top_A[1:i];	
		memB=top_B[1:i];
		shared=intersect(memA, memB);
		num_shared[i]=length(shared);
	}

	max_shared=num_shared[num_top];
	plot(1:num_top, num_shared,
		main="Accumulation of Shared Among Top Categories",
		xlab="Top X Categories",
		ylab="Number of Shared Categories",
		type="n"
	)
	
	abline(h=seq(0, max_shared, 5), col="grey");
	points(1:num_top, 1:num_top, col="darkgreen", type="l", lty="dashed");
	points(1:num_top, num_shared);

	mtext("Green line represents 100% concordance.", side=3, at=0, adj=0,
		cex=.5, col="darkgreen");

	median_concordance=median(num_shared/(1:num_top));
	
	mtext(sprintf("Median concordance across top categories: %3.2f%%", 100*median_concordance), side=3, at=num_top, adj=1,
		cex=.5);
	

	return(shared);
}

###############################################################################

plot_ratio=function(profA, profB, whichTaxa, nameA="A", nameB="B"){

	keptA=profA[whichTaxa];
	keptB=profB[whichTaxa];
	num_taxa=length(whichTaxa);

	log_ratio=log10(keptA/keptB);

	ratio=paste(nameA, "/", nameB, sep="");

	par(mar=c(13.2, 4.1, 4.1, 5.1));
	barmids=barplot(log_ratio, las=2, 
		main=paste("Log10(", ratio, ") of Top Shared Categories Abundances", sep=""),
		axisnames=F);

	graph_bottom=min(log_ratio);
	spacing=mean(diff(barmids));

        for(i in 1:num_taxa){
                text(barmids[i]-spacing*.25, graph_bottom-.1, whichTaxa[i], srt=-45, xpd=T, cex=35/num_taxa, pos=4);
        }

	for(i in 1:num_taxa){
		bottom=min(log_ratio[i], 0);
		points(
			c(barmids[i], barmids[i]),
			c(bottom, graph_bottom),
			type="l", col="black", lty=3, lwd=.5
		);
        }

	par(mar=MAR_DEF);
}

###############################################################################

plot_evenness=function(profA, profB, num_top, nameA="A", nameB="B"){
	profA_sorted=sort(profA, decreasing=T)[1:num_top];
	profB_sorted=sort(profB, decreasing=T)[1:num_top];
	
	profA_cumsum=cumsum(profA_sorted);
	profB_cumsum=cumsum(profB_sorted);

	#plot(0,0, type="n", ylim=c(0,1), xlim=c(0, num_top*1.1));

	ratio=paste(nameA, "/", nameB, sep="");

	plot(log10(profA_sorted/profB_sorted), type="b",
		main=paste("Log10(", ratio, ") of Ranked Abundances\nComparison of Evenness", sep=""),
		xlab="Category Rank",
		ylab=paste("Log10(", ratio, ")", sep="")
		);
	abline(h=0, col="grey");

}

###############################################################################

plot_topratios=function(profA, profB, num_top, nameA="A", nameB="B"){

	profA_sorted=sort(profA, decreasing=T)[1:num_top];

	Asorted_names=names(profA_sorted);
	profB_names=names(profB);

	profB_byA=numeric();

	for(i in 1:num_top){
		a_name=Asorted_names[i];
		if(any(a_name==profB_names)){
			profB_byA[i]=profB[a_name];
		}else{
			profB_byA[i]=0;
		}
	}
	
	log_ratio=log10(profA_sorted/profB_byA);
	finite=is.finite(log_ratio);
	finite_range=range(log_ratio[finite]);
	finite_max=finite_range[2];

	log_ratio[!finite]=finite_max*1.1;

	colors=rep("grey", num_top);
	colors[!finite]="black";

	ratio=paste(nameA, "/", nameB, sep="");

	par(mar=c(11.2, 4.1, 4.1, 4.1));
	barmids=barplot(log_ratio, las=2, 
		main=paste("Log10(", ratio, ") of ", nameA, "'s Top Categorical Abundances", sep=""),
		ylim=finite_range*1.1,
		ylab=paste("Log10(", ratio, ")", sep=""),
		col=colors,
		axisnames=F);

	graph_bottom=min(log_ratio);
	spacing=mean(diff(barmids));

        for(i in 1:num_top){
		if(!finite[i]){
			text(barmids[i], finite_max*1.1, quote(infinity) , 
				xpd=T, cex=70/num_top, pos=3);
		}
	}

	for(i in 1:num_top){
		bottom=min(log_ratio[i], 0);
		points(
			c(barmids[i], barmids[i]),
			c(bottom, graph_bottom),
			type="l", col="black", lty=3, lwd=.5
		);
        }

        for(i in 1:num_top){
                text(barmids[i]-spacing, graph_bottom-.1, Asorted_names[i], 
			srt=-45, xpd=T, cex=35/num_top, pos=4);
        }

	par(mar=MAR_DEF);

}

###############################################################################

plot_rankabund=function(prof, num_top, name="Sample"){

	prof_sorted=sort(prof, decreasing=T)[1:num_top];
	sorted_names=names(prof_sorted);

	cumsum=cumsum(prof_sorted);
	
	par(mar=c(11.2, 4.1, 4.1, 4.1));
	barmids=barplot(prof_sorted,axisnames=F,
		ylab="Abundance",
		main=paste("Rank Abundance of ", name, sep="")
		);
	spacing=mean(diff(barmids));
	graph_bottom=0;

	max_val=max(prof_sorted);

        for(i in 1:num_top){
                text(barmids[i]-spacing, graph_bottom-(max_val*.05), sorted_names[i], 
			srt=-45, xpd=T, cex=35/num_top, pos=4);
        }

	par(mar=MAR_DEF);

	return(sorted_names);

}

###############################################################################
###############################################################################

# remove low abundance that has been group together as 'Remaining'
prof_a_colnames=colnames(prof_a_norm_mat);
prof_b_colnames=colnames(prof_b_norm_mat);

prof_a_rem_idx=which(prof_a_colnames==REMAINING);
prof_b_rem_idx=which(prof_b_colnames==REMAINING);

if(length(prof_a_rem_idx)!=0){
	cat("Index of 'Remaining' in Profile A: ", prof_a_rem_idx, sep="", "\n");
	prof_a_norm_mat=prof_a_norm_mat[,-prof_a_rem_idx, drop=F];
}

if(length(prof_b_rem_idx)!=0){
	cat("Index of 'Remaining' in Profile B: ", prof_b_rem_idx, sep="", "\n");
	prof_b_norm_mat=prof_b_norm_mat[,-prof_b_rem_idx, drop=F];
}

used_taxa=c();

IND_TOP=min(c(IND_TOP, ncol(prof_a_norm_mat), ncol(prof_b_norm_mat)));
cat("Num top actually used: ", IND_TOP, "\n");

for(i in 1:num_samples){

	par(mfrow=c(3,1));
	shared_taxa_list=plot_shared_top(prof_a_norm_mat[i,], prof_b_norm_mat[i,], COMBO_TOP, nameA, nameB);
	plot_ratio(prof_a_norm_mat[i,], prof_b_norm_mat[i,], shared_taxa_list, nameA, nameB);
	plot_evenness(prof_a_norm_mat[i,], prof_b_norm_mat[i,], IND_TOP, nameA, nameB);

	par(mfrow=c(2,1));
	trans_top_taxa=plot_rankabund(prof_a_norm_mat[i,], IND_TOP, nameA);
	plot_topratios(prof_a_norm_mat[i,], prof_b_norm_mat[i,], IND_TOP, nameA, nameB);
	genom_top_taxa=plot_rankabund(prof_b_norm_mat[i,], IND_TOP, nameB);
	plot_topratios(prof_b_norm_mat[i,], prof_a_norm_mat[i,], IND_TOP, nameB, nameA);

	used_taxa=unique(c(used_taxa, shared_taxa_list, trans_top_taxa, genom_top_taxa));
}

###############################################################################

dev.off();
cat("Done.\n")
if(!is.null(warnings())){
	print(warnings());
}
q(status=0)
