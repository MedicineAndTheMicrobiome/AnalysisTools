#!/usr/bin/env Rscript

###############################################################################
#                                                                             #
#       Copyright (c) 2013 J. Craig Venter Institute.                         #
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

params=c(
	"input_file", "i", 1, "character",
	"output_file_root", "o", 2, "character",
	"color_map", "c", 2, "character",
	"height", "h", 2, "numeric",
	"width", "w", 2, "numeric",
	"ubiquity", "u", 2, "numeric",	
	"abundance", "a", 2, "numeric",
	"keep_list_file", "k", 2, "character",
	"grey_subcutoff_taxa", "r", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

PaperHeight=8.5;
PaperWidth=11;


usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"	-i <input CDF file>\n",
	"	[-c <color scheme, taxa to color map>]\n",
	"\n",
	"	[-u <ubiquity cutoff, eg. .5 is 50%>]\n",
	"	[-a <abundance cutoff, eg. .01 is 1%>]\n",
	"	[-k <keep list file (default all taxa)>\n",
	"	[-r (flag to plot all taxa, but use grey and do not label taxa not exceeding cutoff)]\n",
	"\n",
	"	[-o <output filename root>]\n",
	"	[-h <paper height, default=", PaperHeight, ">]\n",
	"	[-w <paper width>]\n",
	"\n",
	"This script will plot the Ub-Ab CDF curves for the non-Major core taxa.\n",
	"\n",
	"The non-core are:\n",
	"	Low Abundance but High Ubiquity (These are the Minor Core)\n",
	"	High Abundance but Low Ubiquity (These may be sample outliers)\n",
	"\n",
	"\n", sep="");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;

if(length(opt$output_file_root)){
	OutputFileRoot=opt$output_file_root;
}else{
	OutputFileRoot=gsub(".cdf", "", opt$input_file);
}
OutputFileRoot=paste(OutputFileRoot, ".otherCore", sep="");

if(length(opt$height)){
	PaperHeight=opt$height;
}

if(length(opt$width)){
	PaperWidth=opt$width;
}

if(length(opt$ubiquity)){
	Ubiquity=opt$ubiquity;
}else{
	Ubiquity=.5;
}

if(length(opt$abundance)){
	Abundance=opt$abundance;
}else{
	Abundance=.01;
}

if(length(opt$keep_list_file)){
	KeepList=opt$keep_list_file;
}else{
	KeepList="";
}

if(length(opt$grey_subcutoff_taxa)){
	GreySubcutoffTaxa=TRUE;
}else{
	GreySubcutoffTaxa=FALSE;
}

ColorMap=NULL;
if(length(opt$color_map)){
	ColorMap=opt$color_map;
}

cat("Input File List: ", InputFileName, "\n");
cat("Output File Root: ", OutputFileRoot, "\n");
cat("\n");
cat("Ubiquity Cutoff: ", Ubiquity, "\n");
cat("Abundance Cutoff: ", Abundance, "\n");
cat("Grey Subcutoff Taxa? ", GreySubcutoffTaxa, "\n");

cat("\nPaper Dimensions: h: ", PaperHeight, " x w: ", PaperWidth, "\n", sep="");

if(!(KeepList=="")){
	cat("Keep List: ", KeepList, "\n");
}

################################################################################

read_keeplist_from_file=function(keeplist_fn){
	keep_list=as.character(as.matrix(read.table(keeplist_fn, sep="\t", header=F)));
	return(as.vector(keep_list));
}

################################################################################

read_cdf_from_file=function(cdf_fn){
# Read precomputed CDF presence values from file

	cdf_info=list();
	cdfs=as.matrix(read.table(cdf_fn, sep=",", header=TRUE, row.names=1, check.names=FALSE));
	
	cdf_info$fname=cdf_fn;
	cdf_info$cutoffs=as.numeric(colnames(cdfs));
	cdf_info$taxa=gsub("X","",rownames(cdfs));
	cdf_info$cdf_val=cdfs;
	cdf_info$num_taxa=length(cdf_info$taxa);

	return(cdf_info);
}

################################################################################

optimize_lines=function(x, y){
# Given two vectors describing lines, it will remove the lines that are redundant.

        diff_f=abs(c(1,diff(y)))>0;
        diff_r=abs(c(rev(diff(rev(y))),1))>0;
        keep=diff_f | diff_r;

        ol=list();
        ol$x=x[keep];
        ol$y=y[keep];

        return(ol);
}

################################################################################

plot_cdf=function(cdf_info, max_abund=NULL, max_ubiq=NULL, mark_ubiquity=0, mark_abundance=0, type){
# Plot the CDF

	num_taxa=length(cdf_info$taxa);
	cat("Num taxa to plot: ", num_taxa, "\n");
	if(num_taxa==0){
		cat("Have all taxa been filtered?  There's nothing left to plot.\n");	
		return();
	}

	non_zero_min=sort(cdf_info$cutoffs)[2];
	log_min=log(non_zero_min,10);

	cat("Using ", non_zero_min, " as minimum non zero abundance.\n", sep="");

	if(is.null(max_abund)){
		max_abund=0;
	}
	if(is.null(max_ubiq)){
		max_ubiq=1;
	}

	# Compute space to give for taxa labels
	max_taxa_name_len=max(nchar(cdf_info$taxa));
	cat("Longest taxa name length: ", max_taxa_name_len, "\n");
	left_name_space=((max_abund-log_min)/4)*max_taxa_name_len/30;
	bottom_name_space=1*max_taxa_name_len/100;
	cat("Left Name Space: ", left_name_space, "\n");
	cat("Bottom Name Space: ", bottom_name_space, "\n");

	# Set up blank plot
	par(mar=c(5,5,5,1));
	plot(0,0, type="n", xlab=expression(Log[10](Abundance)), ylab="Percent Ubiquity", main=InputFileName,
		xlim=c(log_min-left_name_space,max_abund), ylim=c(-bottom_name_space,max_ubiq),
		yaxt="n", xaxt="n", 
		);

	# y-axis labels
	ubiquity_ticks=(seq(0,100,10))/100;
	ubiquity_labels=sprintf("%i", ubiquity_ticks*100);
	axis(side=2, at=ubiquity_ticks, labels=ubiquity_labels, tick=TRUE, outer=FALSE, las=2);

	# x-axis labels
	abundance_ticks=(log_min:0);
	abundance_labels=sprintf("%i", abundance_ticks);
	#axis(side=3, at=abundance_ticks, labels=abundance_labels, tick=TRUE, outer=FALSE, las=1);
	axis(side=1, at=abundance_ticks, labels=abundance_labels, tick=TRUE, outer=FALSE, las=2);


	if(cdf_info$cutoffs[1]==0){
		cdf_info$cutoffs[1]=non_zero_min;
	}
	log_cutoffs=log(cdf_info$cutoffs,10);

	# Plot each taxa's line
	if(num_taxa>0){

		grey_idx=which(cdf_info$colors=="grey");
		non_grey_idx=which(cdf_info$colors!="grey");	
		draw_order=c(grey_idx, non_grey_idx);	

		for(t in draw_order){

			ol=optimize_lines(log_cutoffs, cdf_info$cdf_val[t,]);

			lines(ol$x, ol$y, 
				col=cdf_info$colors[t], 
				lty=cdf_info$line_type[t], lwd=1.5);

		}

		idx=numeric(num_taxa);
		for(t in 1:num_taxa){
			# Determine where curve touches the x-axis
			ubiquities=cdf_info$cdf_val[t,];
			min_nonzero=min(ubiquities[ubiquities>0]);
			idx[t]=max(which(min_nonzero==ubiquities));
		}

		# Label below x axis
		text(log_cutoffs[idx]+.08, bottom_name_space*-.02, cdf_info$taxa,
			 cex=.3, col=cdf_info$colors, srt=90, pos=2);
		
	
		# Label along y axis on the left
		text(log_min,  cdf_info$cdf_val[,2], cdf_info$taxa, cex=.3, col=cdf_info$colors, pos=2);

	}
	
	# Draw the point where the cutoffs are
	log_mark_abundance=log(mark_abundance,10)

	range_col="midnightblue";
	subtitle="";
	if(type=="HiUb_LoAb"){
		points(log_min,mark_ubiquity, col=range_col, pch=16);
		points(log_min,1, col=range_col, pch=16);
		lines(x=c(log_min, log_min),y=c(mark_ubiquity,1), lwd=3, col=range_col); #Top left bar

		points(log_min, 0, col=range_col, pch=16);
		points(log_mark_abundance, 0, col=range_col, pch=16);
		lines(x=c(log_min, log_mark_abundance),y=c(0,0), lwd=3, col=range_col); #Bottom left bar

		subtitle="The Minor Core: High Ubiquity w/ Low Abundance";

	}else if(type=="LoUb_HiAb"){
		points(log_min,0, col=range_col, pch=16);
		points(log_min,mark_ubiquity, col=range_col, pch=16);
		lines(x=c(log_min, log_min),y=c(0, mark_ubiquity), lwd=3, col=range_col); #Bottom left bar

		points(0, 0, col=range_col, pch=16);
		points(log_mark_abundance, 0, col=range_col, pch=16);
		lines(x=c(log_mark_abundance,0),y=c(0,0), lwd=3, col=range_col); #Bottom right bar

		subtitle="The Outliers: Low Ubiquity w/ High Abundance";
	}

	mtext(subtitle, line=1);

}

################################################################################

load_color_map=function(colormapfn){
# Load a color mapping from a file

	color_map=list();
	map=as.matrix(read.table(colormapfn, sep="\t", header=FALSE, comment.char=""));
	num_entries=nrow(map);
	for(i in 1:num_entries){
		clean=gsub("^ *", "", map[i,1]);
		clean=gsub(" *$", "", clean);
		color_map[[clean]]=map[i,2];
	}
	return(color_map);
}

################################################################################

taxa_to_color=function(taxa_list, color_map){
# Assign a color to each taxa, in the same order

	num_taxa=length(taxa_list);
	colors=rep("", num_taxa);
	if(num_taxa>0){
		for(i in 1:num_taxa){
			color=color_map[[taxa_list[i]]];
			if(!is.null(color)){
				colors[i]=color;
			}	
		}
	}
	return(colors);	
}

################################################################################

filter_cdfs_byKeeplist=function(cdf_info, keeplist_vector){
# Removes taxa not in the keeplist
	num_taxa=cdf_info$num_taxa;
	keep_index=rep(FALSE, num_taxa);
	for(i in 1:num_taxa){
		if(any(cdf_info$taxa[i]==keeplist_vector)){
			keep_index[i]=TRUE;
		}
	}
	cdf_info$cdf_val=cdf_info$cdf_val[keep_index,];
	cdf_info$taxa=cdf_info$taxa[keep_index];
	cdf_info$num_taxa=sum(keep_index);

	return(cdf_info);
}

################################################################################

filter_cdfs=function(cdf_info, ubiquity, abundance, type){
# Removes taxa that do not meet ubiquity/abundance cutoff

	# Get lower bound
	non_zero_lowerbound=min(cdf_info$cutoffs[cdf_info$cutoffs>0]);
	non_zero_lowerbound_idx=min(which(cdf_info$cutoffs>0));
	cat("Non-zero lowerbound: ", non_zero_lowerbound, "\n", sep="");
	cat("Non-zero lowerbound idx: ", non_zero_lowerbound_idx, "\n", sep="");

	# Get ubiquity at abundance cutoff
	ub_at_ab=max(which(cdf_info$cutoffs<=abundance));
	ab_cutoff=cdf_info$cutoffs[ub_at_ab];
	cat("Abundance Cutoff: ", ab_cutoff, "\n", sep="");
	cat("Abundance Cutoff idx: ", ub_at_ab, "\n", sep="");

	keep_idx=rep(TRUE,cdf_info$num_taxa);

	if(type=="HiUb_LoAb"){

		for(i in 1:cdf_info$num_taxa){
			if(cdf_info$cdf_val[i, non_zero_lowerbound_idx]<ubiquity){
				keep_idx[i]=FALSE;
			}
			if(cdf_info$cdf_val[i, ub_at_ab]>0){
				keep_idx[i]=FALSE;
			}
		}

	}else if(type=="LoUb_HiAb"){

		for(i in 1:cdf_info$num_taxa){
			if(cdf_info$cdf_val[i, non_zero_lowerbound_idx]>ubiquity){
				keep_idx[i]=FALSE;
			}
			if(cdf_info$cdf_val[i, ub_at_ab]==0){
				keep_idx[i]=FALSE;
			}
		}

	}
	num_kept=sum(keep_idx);

	# Store results back into cdf_info
	cdf_info$cdf_val=cdf_info$cdf_val[keep_idx,];
	cdf_info$taxa=cdf_info$taxa[keep_idx];
	cdf_info$num_taxa=num_kept;

	if(cdf_info$num_taxa==1){
		# Annoying R feature that automatically changes single row matrices to vectors.
		cdf_info$cdf_val=matrix(cdf_info$cdf_val, nrow=1, ncol=length(cdf_info$cutoffs));
	}
	
	return(cdf_info);
}

################################################################################

# Read in color map
color_map=NULL;
if(!is.null(ColorMap)){
	color_map=load_color_map(ColorMap);
}

# Read in CDF info
cdf_info=read_cdf_from_file(InputFileName);

# Filter by list
if(KeepList!=""){
	keeplist_vector=read_keeplist_from_file(KeepList);
	cat("Keep List Members: \n");
	print(keeplist_vector);
	cdf_info=filter_cdfs_byKeeplist(cdf_info, keeplist_vector);
}


ub_ext=sprintf(".u%0.2f", Ubiquity);
ab_ext=sprintf(".a%0.5f", Abundance);
pdf(paste(OutputFileRoot, ub_ext, ab_ext, ".pdf", sep=""), height=PaperHeight, width=PaperWidth);

for(filter_type in c("HiUb_LoAb", "LoUb_HiAb")){

	cat("\nWorking on Filter Type: ", filter_type, "\n");

	type_ext=paste(".", filter_type, sep="");

	if(!GreySubcutoffTaxa){
		# Filter by cutoffs
		f_cdf_info=filter_cdfs(cdf_info, Ubiquity, Abundance, filter_type);

		# Assign colors based on order in taxa of cdf_info
		colors=taxa_to_color(f_cdf_info$taxa, color_map);

		# Assign grey scale and line types
		line_types=rep(1, f_cdf_info$num_taxa);
		oc=0;
		if(f_cdf_info$num_taxa>0){
			for(i in 1:f_cdf_info$num_taxa){
				if(colors[i]==""){
					grey_idx=(oc%%3)+1;
					line_types[i]=(oc%%5)+1;
					if(grey_idx==1){
						colors[i]="gray0";
					}else if(grey_idx==2){
						colors[i]="gray30";
					}else if(grey_idx==3){
						colors[i]="grey60";
					}
					oc=oc+1;
				}
			}
		}
		f_cdf_info$colors=colors;
		f_cdf_info$line_type=line_types;
	}else{

		# Get colors for all taxa
		colors=taxa_to_color(cdf_info$taxa, color_map);

		for(i in 1:cdf_info$num_taxa){
			if(colors[i]==""){
				colors[i]="black";
			}
		}

		# Get a list of taxa that exceed the cutoff
		exceeding_cdf_info=filter_cdfs(cdf_info, Ubiquity, Abundance, filter_type);
		above_cutoff_taxa=exceeding_cdf_info$taxa;

		# If the taxa exceeds the cutoff, leave it alone, remove the name and color it grey
		f_cdf_info=cdf_info;
		for(i in 1:f_cdf_info$num_taxa){
			if(!any(f_cdf_info$taxa[i]==above_cutoff_taxa)){
				colors[i]="grey";
				f_cdf_info$taxa[i]="";
			}
		}
		
		f_cdf_info$colors=colors;
		f_cdf_info$line_type=rep(1, cdf_info$num_taxa);
	}


	################################################################################


	plot_cdf(f_cdf_info, mark_ubiquity=Ubiquity, mark_abundance=Abundance, type=filter_type);


	################################################################################

	if(Ubiquity>0 || Abundance>0){
		kept_taxa=filter_cdfs(cdf_info, Ubiquity, Abundance, filter_type);
	
		peak_ubiquity_idx=min(which(cdf_info$cutoffs>0));

		fh=file(paste(OutputFileRoot, ub_ext, ab_ext, ".", filter_type, ".csv", sep=""), "w");

		num_exc_cutoff=length(order);

		cat(file=fh, "Filename:,", OutputFileRoot, "\n", sep="");
		cat(file=fh, "Ubiquity:,", Ubiquity, "\n", sep="");
		cat(file=fh, "Abundance:,", Abundance, "\n", sep="");
		cat(file=fh, "Type:,", filter_type, "\n", sep="");
		cat(file=fh, "NumTaxa Exceeding Cutoffs:,", kept_taxa$num_taxa, "\n", sep="");
		cat(file=fh, "\n");

		cat(file=fh, "Index,TaxaName,PeakUbiquity,PeakAbundance\n", sep="");

		if(kept_taxa$num_taxa>0){
			for(i in 1:kept_taxa$num_taxa){
				peak_abundance_idx=min(which(kept_taxa$cdf_val[i,]==0));

				lineout=paste(c(
					i, 
					kept_taxa$taxa[i], 
					kept_taxa$cdf_val[i,peak_ubiquity_idx], 
					kept_taxa$cutoffs[peak_abundance_idx]), 
					collapse=",");

				cat(file=fh, lineout, "\n", sep="");
			}
		}

		close(fh);
	}

}

################################################################################

dev.off();

w=warnings();
if(length(w)){
	print(w);
}
cat("Done.\n")

q(status=0)
