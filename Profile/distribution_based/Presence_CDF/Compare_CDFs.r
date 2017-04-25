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
	"input_file_A", "A", 1, "character",
	"input_file_B", "B", 1, "character",
	"id_map_file", "m", 2, "character",

	"min_ubiq_diff", "d", 2, "numeric",
	"min_avg_abund", "a", 2, "numeric",

	"output_file_root", "o", 2, "character",
	"color_map", "c", 2, "character",

	"height", "h", 2, "numeric",
	"width", "w", 2, "numeric",
	"keep_list_file", "k", 2, "character",
	"top_names", "t", 2, "numeric",

	"shorten", "s", 2, "logical",
	"suppress_coordinates", "q", 2, "logical",
	"label_size", "l", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

PaperHeight=8.5;
PaperWidth=8.5;

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"	-A <input CDF A file>\n",
	"	-B <input CDF B file>\n",
	"	[-m <id map file>]\n",
	"	[-c <color scheme, taxa to color map>]\n",
	"\n",
	"      Filtering options:\n",
	"	[-d <min difference in ubiquity to plot (between 0.0 and 1.0), default=0.0>]\n",
	"	[-a <min abundance at max ubiquity difference to plot (between 0.0 and 1.0), default=0.0>]\n",
	"	[-k <keep list file, default=all)]>\n",
	"\n",
	"      Labelling options:\n",
	"	[-t <number (integer) of top names to plot, default=all>]\n",
	"\n",
	"      Output options:\n",
	"	[-o <output filename root>]\n",
	"	[-h <paper height, default=", PaperHeight, " in>]\n",
	"	[-w <paper width, default=", PaperWidth, " in>]\n",
	"	[-s (shorten category/taxa names flag, default=F)]\n",
	"	[-q (suppress printing coordinates next to category/taxa, default=F)]\n",
	"	[-l <label size, default=1>]\n",
	"\n",
	"This script will generate a Ubiquit-Ubiquity (U-U) Plot, that\n",
	"compares the ubiquities between two cohorts at matching abundances.\n",
	"Each curve represents a different taxa.\n",
	"\n",
	"Because the number of taxa may produce a cluttered plot, you\n",
	"may wish to try one of the filtering options.\n",
	"  1.) Min difference filtering removes taxa that are closest to the diagonal reference line.\n",
	"  2.) Min abundance filtering removes taxa that are below the specified abundance cutoff.\n",
	"      The abundance examined is at the max difference of ubiquity between the two cohorts.\n",
	"  3.) Keep list only plots curves specified in the list file.\n",
	"\n",
	"Labelling options:\n",
	"  Top names only keeps the top N taxa ordered by differences in ubiquity in both magnitudes.\n",
	"\n", sep="");

if(!length(opt$input_file_A) || !length(opt$input_file_B)){
	cat(usage);
	q(status=-1);
}

InputFileNameA=opt$input_file_A;
InputFileNameB=opt$input_file_B;

if(length(opt$output_file_root)){
	OutputFileRoot=opt$output_file_root;
}else{
	# Reduce file name length
	Aroot=gsub("\\.cdf$","",InputFileNameA);
	Aroot=gsub("\\.summary_table\\.xls", "", Aroot);
	Aroot=gsub("\\.summary_table\\.tsv", "", Aroot);

	Broot=gsub("\\.cdf$","",InputFileNameB);
	Broot=gsub("\\.summary_table\\.xls", "", Broot);
	Broot=gsub("\\.summary_table\\.tsv", "", Broot);

	OutputFileRoot=paste(Aroot, "_vs_", Broot, sep="");
}

if(length(opt$height)){
	PaperHeight=opt$height;
}

if(length(opt$width)){
	PaperWidth=opt$width;
}

if(length(opt$min_ubiq_diff)){
	MinUbiqDiff=opt$min_ubiq_diff;
}else{
	MinUbiqDiff=0;
}

if(length(opt$min_avg_abund)){
	MinAvgAbund=opt$min_avg_abund;
}else{
	MinAvgAbund=0;
}

if(length(opt$keep_list_file)){
	KeepList=opt$keep_list_file;
}else{
	KeepList="";
}

if(length(opt$id_map_file)){
	IDMapFile=opt$id_map_file;
}else{
	IDMapFile="";
}

ColorMap=NULL;
if(length(opt$color_map)){
	ColorMap=opt$color_map;
}

NumTopNamesToPlot=-1;
if(length(opt$top_names)){
	NumTopNamesToPlot=opt$top_names;
}

SuppressCoordinates=F;
if(length(opt$suppress_coordinates)){
	SuppressCoordinates=T;
}

ShortenNames=F;
if(length(opt$shorten)){
	ShortenNames=T;
}

LabelSize=1;
if(length(opt$label_size)){
	LabelSize=opt$label_size;
}

################################################################################

cat("Input PDF File A: ", InputFileNameA, "\n");
cat("Input PDF File B: ", InputFileNameB, "\n");
cat("Output File Root: ", OutputFileRoot, "\n");
cat("\n");
cat("Min Ubiquity Diff: ", MinUbiqDiff, "\n");
cat("Min Abundance at Diff: ", MinAvgAbund, "\n");

cat("\nPaper Dimensions: h: ", PaperHeight, " x w: ", PaperWidth, "\n\n", sep="");

if(!(KeepList=="")){
	cat("Keep List: ", KeepList, "\n");
}
if(!(IDMapFile=="")){
	cat("ID Map File: ", IDMapFile, "\n");
}	

################################################################################

read_keeplist_from_file=function(keeplist_fn){

	clean=character();

	tryCatch({
		keep_list=as.character(as.matrix(read.table(keeplist_fn, sep="\t", header=F)));
		clean=gsub("^ *","", keep_list);
		clean=gsub(" *$","", clean);
		#print(clean);
	}, error = function(e){
		cat("Error Reading Keep List: ", keeplist_fn, "\n");
	});

	return(as.vector(clean));
}

################################################################################

read_map_file=function(map_fn){

	map_info=as.matrix(read.table(map_fn, sep="\t", header=F));
	num_maps=nrow(map_info);
	num_cols=ncol(map_info);

	if(num_cols < 2){
		cat("Error:  Not detecting enough columns in mapping file: ", map_fn, "\n", sep="");
	}
	mapping=list();
	for(i in 1:num_maps){
		mapping[[map_info[i,1]]]=map_info[i, 2];
	}
	return(mapping);
}


################################################################################

read_cdf_from_file=function(cdf_fn){
# Read precomputed CDF presence values from file

	cdf_info=list();
	cdfs=as.matrix(read.table(cdf_fn, sep=",", header=TRUE, row.names=1, check.names=FALSE));
	
	cdf_info$fname=cdf_fn;
	cdf_info$cutoffs=as.numeric(colnames(cdfs));
	cdf_info$taxa=rownames(cdfs);
	cdf_info$cdf_val=cdfs;
	cdf_info$num_taxa=length(cdf_info$taxa);
	
	return(cdf_info);
}

################################################################################

plot_compare_cdf=function(cdf_infoA, cdf_infoB, min_ubiq_diff=0, min_avg_abund=0, 
	interpolation_method="average", id_mapping_info=NULL,
	label_size=1, shorten_names=F){
# Plot the CDF

	num_taxaA=length(cdf_infoA$taxa);
	num_taxaB=length(cdf_infoB$taxa);
	cat("Num Taxa in A: ", num_taxaA, "\n");
	cat("Num Taxa in B: ", num_taxaB, "\n");

	# Compute what is shared between A and B
	shared_taxa=sort(intersect(cdf_infoA$taxa, cdf_infoB$taxa));
	num_shared=length(shared_taxa);
	cat("Shared Taxa: ", num_shared, "\n");
	
	# Get indices of shared taxa
	sharedA_idx=numeric();
	sharedB_idx=numeric();
	for(i in 1:num_shared){
		sharedA_idx[i]=which(cdf_infoA$taxa==shared_taxa[i])	
		sharedB_idx[i]=which(cdf_infoB$taxa==shared_taxa[i])	
	}

	#print(sharedA_idx);
	#print(sharedB_idx);

	#----------------------------------------------------------------------
	# Set up blank plot

	A_name_len=nchar(cdf_infoA$fname);
	B_name_len=nchar(cdf_infoB$fname);
	MAX_DISP_LEN=65;
	
	A_disp=substr(cdf_infoA$fname, A_name_len-MAX_DISP_LEN+1, A_name_len);
	B_disp=substr(cdf_infoB$fname, B_name_len-MAX_DISP_LEN+1, B_name_len);

	if(A_name_len-MAX_DISP_LEN+1 > 1){
		A_disp=paste("...", A_disp, sep="");
	}

	if(B_name_len-MAX_DISP_LEN+1 > 1){
		B_disp=paste("...", B_disp, sep="");
	}

	plot(0,0, type="n", 
		main="Comparison of Ubiquities",
		xlab=paste("Ubiquity of (", A_disp, ")", sep=""),
		ylab=paste("Ubiquity of (", B_disp, ")", sep=""),
		xlim=c(0,1), ylim=c(0,1));

	# Draw common line
	abline(a=0, b=1, lwd=10, col="grey");
	num_plotted=0;

	if(num_shared>0){
		#----------------------------------------------------------------------
		# Compute Draw order so less interesting is draw first (underneath)
		ks_pos=numeric();
		pos_idc=numeric();
		ks_neg=numeric();
		neg_idc=numeric();
		all_diff=numeric();
		pos_idx=1;
		neg_idx=1;
	
		for(t in 1:num_shared){
			# Get shared indices
			idx_A=sharedA_idx[t];
			idx_B=sharedB_idx[t];

			# Get difference between all cdf points for taxa t
			tsA=cdf_infoA$cdf_val[idx_A,];
			tsB=cdf_infoB$cdf_val[idx_B,];
			diff=tsA-tsB;

			# Compute abs, min and max differences
			all_diff[t]=max(abs(diff));
			max_diff=max(diff);
			min_diff=abs(min(diff));

			# Keep track of min and max differences
			if(max_diff>min_diff){
				ks_pos[pos_idx]=all_diff[t];
				pos_idc[pos_idx]=t;
				pos_idx=pos_idx+1;
			}else{
				ks_neg[neg_idx]=all_diff[t];
				neg_idc[neg_idx]=t;
				neg_idx=neg_idx+1;
			}
		}

		# Determine which taxa have the greatest shift
		ks_pos_top=order(ks_pos, decreasing=TRUE);
		ks_neg_top=order(ks_neg, decreasing=TRUE);

		# Get indices of taxa that have the greatest shift
		if(NumTopNamesToPlot==-1){
			NumTopNamesToPlot=num_shared;
		}
		label_indices=numeric();
		if(pos_idx>0){
			num_lab=min(c(pos_idx-1,NumTopNamesToPlot));
			label_indices=c(label_indices, pos_idc[ks_pos_top[1:num_lab]]);
		}
		if(neg_idx>0){
			num_lab=min(c(neg_idx-1,NumTopNamesToPlot));
			label_indices=c(label_indices, neg_idc[ks_neg_top[1:num_lab]]);
		}

		# Determine order of how to draw lines. Draw lines with greatest difference on top, ie. later
		draw_order=order(abs(all_diff),decreasing=FALSE);

		#----------------------------------------------------------------------
		# Plot each taxa's line

		for(t in draw_order){

			#cat("Plotting: ", shared_taxa[t], "\n");
			idx_A=sharedA_idx[t];
			idx_B=sharedB_idx[t];

			ptsA=cdf_infoA$cdf_val[idx_A,];
			ptsB=cdf_infoB$cdf_val[idx_B,];

			# Compute a good place to place labels (where there is a max difference)
			differences=abs(ptsA-ptsB);
			max=max(differences);
			max_diff_idx=max(which(differences==max));
			abund_at_diff=cdf_infoA$cutoffs[max_diff_idx]

			if((max > min_ubiq_diff) && (abund_at_diff > min_avg_abund)){

				# Draw a thicker line, if we're going to label it
				if(any(t==label_indices)){
					line_width=3;	
				}else{
					line_width=1;
				}	

				# Draw lines
				lines(ptsA, ptsB, 
					col=cdf_infoA$colors[idx_A], 
					lty=cdf_infoA$line_type[idx_A],
					lwd=line_width);

				# Change ID's to descriptions
				if(!is.null(id_mapping_info)){
					display_name=id_mapping_info[[shared_taxa[t]]];
				}else{
					display_name=shared_taxa[t];
					if(shorten_names){
						display_name=tail(strsplit(display_name, ";")[[1]],1);
					}
				}

				if(any(t==label_indices)){

					cat("Labeling: ", display_name);

					# Decide if label goes above or below the line
					if(ptsA[max_diff_idx]<ptsB[max_diff_idx]){
						rel_pos=-.75;
						cat(" (above)\n");
					}else{
						rel_pos=2
						cat(" (below)\n");
					}

					# Place a "pin" over the line of the same color as the text
					points(ptsA[max_diff_idx], ptsB[max_diff_idx], 
						col=cdf_infoA$colors[idx_A], pch=19, cex=label_size); 

					# Label the line
					coordinate_string="";
					if(!SuppressCoordinates){
						coordinate_string=sprintf(" (%2.2f,%2.2f)", ptsA[max_diff_idx], ptsB[max_diff_idx]);
					}
					text(ptsA[max_diff_idx], ptsB[max_diff_idx], 
						label=paste(sep="", display_name, coordinate_string), srt=45, 
						col=cdf_infoA$colors[idx_A], adj=c(.5, rel_pos), 
						cex=label_size, font=3);
				}

				num_plotted=num_plotted+1;
			}
		}
	}else{
		text(.5, .5, "Nothing to Plot...", cex=3);
	}

	cat("Num Plotted: ", num_plotted, "\n");
	
}

################################################################################

load_color_map=function(colormapfn){
# Load a color mapping from a file

	color_map=list();
	map=as.matrix(read.table(colormapfn, sep="\t", header=FALSE, comment.char=""));
	num_entries=nrow(map);
	for(i in 1:num_entries){
		clean=gsub("^ *","", as.character(map[i,1]));
		clean=gsub(" *$","", clean);
		color_map[[clean]]=map[i,2];
	}
	#print(color_map);
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
	cdf_info$cdf_val=cdf_info$cdf_val[keep_index,,drop=F];
	cdf_info$taxa=cdf_info$taxa[keep_index];
	cdf_info$num_taxa=sum(keep_index);

	return(cdf_info);
}

################################################################################

filter_cdfs=function(cdf_info, ubiquity, abundance){
# Removes taxa that do not meet ubiquity/abundance cutoff

	# Find index into the abundances
	cat("Filtering for Abundance: ", abundance, "\n");
	abund_idx=min(which(abundance<=cdf_info$cutoffs));
	cat("Index of abundance cutoff: ", abund_idx, "\n");

	# Find taxa that exceed ubiquity
	cat("Filtering for Ubiquity: ", ubiquity, "\n");
	taxa_exc_ubiquity=which(cdf_info$cdf_val[,abund_idx]>=ubiquity);
	num_kept=length(taxa_exc_ubiquity);
	cat("Num taxa kept: ", num_kept, "\n");

	# Store results back into cdf_info
	cdf_info$cdf_val=cdf_info$cdf_val[taxa_exc_ubiquity, , drop=F];
	cdf_info$taxa=cdf_info$taxa[taxa_exc_ubiquity];
	cdf_info$num_taxa=num_kept;

	return(cdf_info);
}

################################################################################

interpolate=function(x, y, new_x, method="linear"){

        # Given known x and y, will try to generate a list of new_y's for each
        # value of new_x by linear extrapolation between the two closest points

        num_newx=length(new_x);
        new_y=numeric(num_newx);

        order=order(x);
        x=x[order];
        y=y[order];

        for(i in 1:num_newx){

                lb_ix=max(which(new_x[i]>=x));
                ub_ix=min(which(new_x[i]<=x));

                ylb=y[lb_ix];   xlb=x[lb_ix];
                yub=y[ub_ix];   xub=x[ub_ix];

                if(method=="linear"){
                        delta_y=yub-ylb;
                        delta_x=xub-xlb;
                        if(delta_y==0){
                                new_y[i]=ylb;   # Slope = 0
                        }else if(delta_x==0){
                                new_y[i]=(ylb+yub)/2;   # Two y's for same x, just average
                        }else{
                                slope=delta_y/(xub-xlb);
                                new_y[i]=(new_x[i]-xlb)*slope + ylb; # if both x and y change, compute slope and interpolate
                        }
                }else if(method=="average"){
                        new_y[i]=(yub+ylb)/2;
		}else{
			cat("Unrecognized interpolation method: ", method, "\n");
			quit(status=-1);
		}

                #cat("Target X:", new_x[i], " idx:", lb_ix, "-", ub_ix, " : X:", xlb, "-", xub, " / Y:", ylb, "-", yub, ": Est ", new_y[i], "\n", sep="");
        }
        return(new_y);
}

################################################################################

# Read in color map
color_map=NULL;
if(!is.null(ColorMap)){
	color_map=load_color_map(ColorMap);
}

# Read id mapping file
idmapping=NULL;
if(IDMapFile!=""){
	idmapping=read_map_file(IDMapFile);
}

# Read in CDF info
cdf_infoA=read_cdf_from_file(InputFileNameA);
cdf_infoB=read_cdf_from_file(InputFileNameB);

# Filter by list
if(KeepList!=""){
	keeplist_vector=read_keeplist_from_file(KeepList);
	cat("Keep List Members: \n");
	#print(keeplist_vector);

	cdf_infoA=filter_cdfs_byKeeplist(cdf_infoA, keeplist_vector);
	cdf_infoB=filter_cdfs_byKeeplist(cdf_infoB, keeplist_vector);
}

# Assign colors based on order in taxa of cdf_info
colors=taxa_to_color(cdf_infoA$taxa, color_map);

cat("Assigning grey scales to uncolored taxa.\n");
# Assign grey scale and line types
line_types=rep(1, cdf_infoA$num_taxa);
oc=0;
if(cdf_infoA$num_taxa>0){
	for(i in 1:cdf_infoA$num_taxa){
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


cdf_infoA$colors=colors;
cdf_infoA$line_type=line_types;
cdf_infoB$colors=colors;
cdf_infoB$line_type=line_types;

################################################################################

pdf(paste(OutputFileRoot, ".UU.pdf", sep=""), height=PaperHeight, width=PaperWidth);

cat("Plotting comparisons:\n");

plot_compare_cdf(cdf_infoA, cdf_infoB, 
	min_ubiq_diff=MinUbiqDiff, min_avg_abund=MinAvgAbund, id_mapping_info=idmapping,
	label_size=LabelSize, shorten_names=ShortenNames);

dev.off();

################################################################################

w=warnings();
if(length(w)){
	print(w);
}
cat("Done.\n")

q(status=0)
