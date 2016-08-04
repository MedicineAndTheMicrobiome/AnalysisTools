#!/usr/bin/env Rscript

###############################################################################
#                                                                             # 
#       Copyright (c) 2009 J. Craig Venter Institute.                         #     
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
NoiseCutoff=.001;
DEFAULT_COMBINED="combined_entropies";

params=c(
        "input_file", "e", 2, "character",
	"output_file", "o", 2, "character",
        "topn", "n", 2, "numeric",
        "threshold", "t", 2, "numeric",
	"noise", "s", 2, "numeric",
	"input_list", "l", 2, "character",
	"annotation", "a", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];


usage = paste (
        "\nUsage:\n\n", script_name,
        "\n",
	"	-e <entropy files, separated by commas and no spaces>\n",
	"		or\n",
	"	-l <list of entropy files>\n",
	"\n",
	"	[-o <output file name>]\n",
	"	[-n <label top n points>]\n",
	"	[-t <label points above t, where 0 < t < 1 >]\n",
	"	[-s <noise threshold, 0 < s < 1, default=", NoiseCutoff, " >]\n",
	"	[-a <annotation file (begin\\tend\\tmask\\tcomments\\n)>]\n",
        "\n",
	"This script will read in one to many entropy file(s) and plot the points into a pdf file.\n",
	"\n",
	"By default the output file name will be based on the input file.\n",
	"If more than one input file is used, the name will be: ", DEFAULT_COMBINED, "\n",
	"\n",
	"Specify the -n option if you want to label the highest scoring positions.\n",
	"Specify the -t option if you want to label all points above thes specified threshold.\n",
	"Specify the -s option to not plot points below the noise threshold.\n",
	"\n",
	"Note: Output and labeled indices start from 1, although input entropies should start from position 0.\n",
	"Format of the input is: <position>\\t<entropy>\\n\n",
	"\n",
	"Use the -a option to highlight regions on the plot with an annotation file.\n",
	"The annotation file should be a list of tab separated coordinates, begin and end.\n",
	"If there is more than 2 columns in this annotation file, the 3rd column is considered\n",
	"a mask, where a value greater than 1 is considered on (highlight the regon).\n",
	"If the 3rd column is 0, the region is not highlighted on the plot.\n",
	"The 4th column is optional and is ignored by this script.\n",
        "\n",
	"Example Usage:\n",
	"	", script_name, " -e data1.ent,data2.ent \n",
        "\n", sep="");

if(!length(opt$input_file) && !length(opt$input_list)){
        cat(usage);
        q(status=-1);
}

#----------------------------------------------------------
# Get input files

input_file_names=c();
if(length(opt$input_file)){
	input_file_names=c(input_file_names, strsplit(opt$input_file, ",")[[1]]);
	num_inputs=length(input_file_names);
	if(num_inputs>1){
		OutputFileName=DEFAULT_COMBINED;
	}else{
		OutputFileName=opt$input_file;
	}
}
if(length(opt$input_list)){
	input_file_names=c(input_file_names, as.matrix(read.delim(opt$input_list, sep="\n", header=FALSE)));
	OutputFileName=opt$input_list;
	num_inputs=length(input_file_names);
}

#----------------------------------------------------------

TopN=0;
if(length(opt$topn)){
	TopN=opt$topn;
}

Threshold=1;
if(length(opt$threshold)){
	Threshold=opt$threshold;
}

if(length(opt$noise)){
	NoiseCutoff=opt$noise;
}

if(length(opt$output_file)){
	OutputFileName=opt$output_file;
}

Annotation="";
if(length(opt$annotation)){
	Annotation=opt$annotation;
}

###############################################################################

# Load specified dataset
entropies=list();
entropy_length=-1;

# For each file, read in the position and the entropy value
max_pos=0;
max_entropy=0;
for(i in 1:num_inputs){
	fname=input_file_names[i];
	cat("Loading '", fname, "'\n", sep="");
	entropies[[fname]]=as.matrix(read.table(fname, check.names=FALSE, sep="\t"));
	colnames(entropies[[fname]])=c("Position", "Entropy");
	max_pos=max(entropies[[fname]][,1], max_pos);
	max_entropy=max(max_entropy, entropies[[fname]][,2]);
}
cat("Max position: ", max_pos, "\n");


# Load annotation file
if(Annotation!=""){
	table=read.table(Annotation, sep="\t", comment.char="#");
	ncols=ncol(table);
	if(ncols>2){
		# If more than 2 columns, assume 3rd column is the mask
		annot_coord_wMask=apply(as.matrix(table[,1:3]), c(1,2), as.numeric);
		mask=annot_coord_wMask[,3]>=1;
		annot_coord=annot_coord_wMask[mask,1:2];
	}else if(ncols==2){
		# If only 2 columns, assume all will be annotated
		annot_coord=apply(as.matrix(table[,1:2]), c(1,2), as.numeric);
	}else{
		cat("Error: Could not understand annotation file: ", Annotation, "\n");
		quit(status=-1);
	}

	if (is.vector(annot_coord)) {
		annot_coord=as.matrix(t(annot_coord));
	} 
	num_annotations=nrow(annot_coord);
	
	cat("Number of regions to annotate: ", num_annotations, "\n");
}

# Apply filter to each entropy profile
for(i in 1:num_inputs){
	fname=input_file_names[i];
	cur_entropy=entropies[[fname]];
	gt_cutoff_idx=which(cur_entropy[,2]<NoiseCutoff);
	cur_entropy[gt_cutoff_idx,2]=0;
	entropies[[fname]]=cur_entropy;
}

#------------------------------------------------------------------------------
# Produce the output

# Change this if you want a different output format
entropy_PDF=paste(OutputFileName, ".pdf", sep="")
pdf(entropy_PDF,width=11,height=8.5)


layout_matrix=t(matrix(c(
	1,1,1,1,1,1,2,2,
	1,1,1,1,1,1,2,2,
	1,1,1,1,1,1,2,2
	), nrow=8));
layout(layout_matrix);
colors=rainbow(num_inputs);
cat("Noise cutoff: ", NoiseCutoff, "\n", sep="");

y_limits=c(1, .5, max_entropy*1.1, .1);

for(plot_idx in 1:length(y_limits)){

	plot_ymax=y_limits[plot_idx];

	par(oma=c(0,0,0,0));
	par(mar=c(10,5,10,0));

	# Initiate empty plot
	plot(0,0, xlim=c(0,max_pos), ylim=c(0,plot_ymax), type="n",
			xlab="Position", ylab="Normalized Entropy");
	title(main=OutputFileName, cex.main=2, font.main=2)
	mtext(text=paste("y zoom: ", sprintf("%3.1f%%", 100/plot_ymax), sep=""), side=3, line=0, cex=.75, font=3);

	if(Annotation!=""){
		for(annot_idx in 1:num_annotations){
			rect(
				annot_coord[annot_idx, 1],
				0,
				annot_coord[annot_idx, 2],
				1,
				border=NA,
				col="grey85");
		}
	}

	abline(h=NoiseCutoff,col="gray50", lwd=.5, lty=2)

	# Plot the points
	for(i in 1:num_inputs){
		fname=input_file_names[i];
		cur_entropy=entropies[[fname]];

		# Don't bother plotting points equal to zero
		gtz=which(cur_entropy[,2]>0);
		points(cur_entropy[gtz,1], cur_entropy[gtz,2], pch=1, cex=.6, col=colors[i]);
	}

	# If we're going to annotate the top points
	if(TopN>0 || Threshold<1){


		# Combine all the entropies together
		all_entropy=numeric();
		for(i in 1:num_inputs){
			fname=input_file_names[i];
			cur_entropy=entropies[[fname]];
			all_entropy=rbind(all_entropy, cur_entropy);
		}

		# Label the top N points
		if(TopN>0){
			sort_rec=sort(all_entropy[,2], index.return=TRUE, decreasing=TRUE);
			all_entropy_sorted=all_entropy[sort_rec$ix,]
			topN_entropies=all_entropy_sorted[1:TopN,];

			x=topN_entropies[,1];
			y=topN_entropies[,2];

			text(x, y, labels=x, cex=.7, adj=c(.5,-.3), col="black");
		}

		# Label the points above the threshold
		if(Threshold<1){
			cat("Labeling above ", Threshold, " entropies.\n", sep="");
			top_indices=which(all_entropy[,2]>=Threshold);

			x=all_entropy[top_indices,1];
			y=all_entropy[top_indices,2];

			text(x, y, labels=x, cex=.7, adj=c(.5,-.3), col="black");

		}

		if(plot_idx==1){
			# Output points in the order of their position
			labeled_txt=paste(OutputFileName, ".labeled.txt", sep="")
			cat("Outputing top indices to: ", labeled_txt, "\n");
			fh=file(labeled_txt, "w");
			num_points=length(x);
			for(i in order(x)){
				cat(file=fh, x[i], y[i], sep="\t");
				cat(file=fh, "\n");
			}
			close(fh);
		}
	}

	# Plot the legend
	par(mar=c(10,1,10,0));
	plot(0,0, xlim=c(-10,100), ylim=c(-100, 0), ylab="", xlab="", type="n", xaxt="n", yaxt="n", bty="n");
	clean_names=gsub(".ent$", "", input_file_names);
	legend(-5,0, legend=clean_names, fill=colors);

}

dev.off();

cat("Done.\n");

