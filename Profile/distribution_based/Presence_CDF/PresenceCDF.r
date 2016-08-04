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
	"output_file", "o", 2, "character",
	"remove_weak_classifications", "w", 0, "logical",
	"in_sample_abundance", "a", 2, "numeric",
	"across_sample_prevalence", "p", 2, "numeric",
	"draw_legend", "L", 0, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"	-i <input summary_table.xls file>\n",
	"	[-o <output cdf file name>]\n",
	"	[-w] (remove weak classifications flag)\n",
	"	[-a <in sample abundance cutoff, eg .10 is 10% within sample]\n",
	"	[-p <across sample prevalence cutoff, eg. .10 is 10% across samples]\n",
	"	[-L] (draw legend instead of trying to label each line)\n",
	"\n",
	"Generates CDF-like graph for the members in the summary table.\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputRoot=gsub(".summary_table.xls", "", InputFileName);

remove_weak_classifications=FALSE;
if(length(opt$remove_weak_classifications)){
	remove_weak_classifications=TRUE;
}

cat("Removing weak classifications: ", remove_weak_classifications, "\n", sep="");

if(length(opt$output_file)){
	OutputRoot=opt$output_file;
}

if(remove_weak_classifications){
	OutputCDF_PDF=paste(OutputRoot, ".no_weak_class.presence_cdf.pdf", sep="");
	OutputDepthFileName=paste(OutputRoot, ".no_weak_class.depth_count.tsv", sep="");
	OutputCDF_Values=paste(OutputRoot, ".no_weak_class.cdf_values.tsv", sep="");
}else{
	OutputCDF_PDF=paste(OutputRoot, ".presence_cdf.pdf", sep="");
	OutputDepthFileName=paste(OutputRoot, ".depth_count.tsv", sep="");
	OutputCDF_Values=paste(OutputRoot, ".cdf_values.tsv", sep="");
}

if(length(opt$in_sample_abundance)){
	in_sample_abundance=opt$in_sample_abundance;
}else{
	in_sample_abundance=0.005;
}

if(length(opt$across_sample_prevalence)){
	across_sample_prevalence=opt$across_sample_prevalence;
}else{
	across_sample_prevalence=0.005;
}

DrawLegend=FALSE;
if(length(opt$draw_legend)){
	DrawLegend=TRUE;
}

###############################################################################
# Main program loop

	cat("Working on ", InputFileName, "\n", sep="");

	# Load summary_table.xls
	mat<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, row.names=1))

	# Exclude total counts column
	count_mat=mat[,2:ncol(mat)]

	categories=as.vector(colnames(count_mat));
	num_categories=length(categories);

	sample_names=rownames(count_mat);

	#print("Original category names:\n");
	#print(categories);

	# Compute shorted names
	short_names=character(num_categories);
	for(i in 1:num_categories){
		taxonomy=unlist(strsplit(categories[i], " "));
		short_names[i]=taxonomy[length(taxonomy)];
	}

	#print("Shortened category names:\n");
	#print(short_names);

	NumSamples=nrow(count_mat);
	NumCategories=ncol(count_mat);

	cat("\n");
	cat("Num Samples: ", NumSamples, "\n");
	cat("Num Categories: ", NumCategories, "\n");
	cat("\n");

	# Normalize
	sample_totals=numeric(NumSamples);
	prob_mat=matrix(nrow=NumSamples,ncol=NumCategories);
	warn=FALSE;
	for(i in 1:NumSamples){
		sample_totals[i]=sum(count_mat[i,]);
		if(sample_totals[i]==0){
			cat("Sample: ", sample_names[i], " had a zero total.\n", sep="");
			warn=TRUE;
		}
		prob_mat[i,]=count_mat[i,]/sample_totals[i];
	}
	if(warn){
		cat("Error:  You can't have 0's as your total, or else your across/within values won't make sense.\n");
		q(status=-1);
	}

	# Compute CDF
	if(in_sample_abundance < .005){
		increm=seq(0,1, in_sample_abundance);
	}else{
		increm=seq(0,1,.005);
	}
	num_increm=length(increm);

	all_cdf_presence=matrix(nrow=NumCategories, ncol=num_increm);
	for(cat_id in 1:NumCategories){

		#cat("Calculating for Sample: ", short_names[cat_id], "\n");
		probs=prob_mat[,cat_id];

		cdf_presence=numeric(0);
		for(i in 1:num_increm){
			presence=sum(probs>=increm[i]);
			cdf_presence[i]=presence/NumSamples;
		}

		#print(cdf_presence);
		all_cdf_presence[cat_id,]=cdf_presence;
	}

	# Define our cutoff
	anno_idx=which(increm==in_sample_abundance);	

	# Precompute which categories will be drawn so we can allocate a colors
	num_samples_to_draw=0;
	for(cat_id in 1:NumCategories){
		if(all_cdf_presence[cat_id, anno_idx] > across_sample_prevalence){
			num_samples_to_draw=num_samples_to_draw+1;
		}
	}
	cat("Number of samples (above cutoff) to draw: ", num_samples_to_draw, "\n");
	if(num_samples_to_draw <= 8){
		colors=c("red", "blue", "green", "orange", "purple", "black", "brown", "yellow4");
	}else{
		colors=rainbow(num_samples_to_draw, start=0, end=0.65);
	}
	cat("Num samples to draw:", num_samples_to_draw, "\n");
	cat("Number of colors generated:", length(colors), "\n");

	# Define which zooms we want to draw
	ranges=c(1, .5, in_sample_abundance*8, in_sample_abundance*4, in_sample_abundance*2);
	ranges=ranges[ranges<=1];
	ranges=unique(ranges);
	cat("Drawing graphs for zoom: ", paste(ranges, collapse=", "), "\n");

	# Open up PDF file to write
	cat("Opening up PDF file for writing:", OutputCDF_PDF, "\n");
	pdf(OutputCDF_PDF, width=11,height=8.5)

	cdf_values_fh=file(OutputCDF_Values, "wt");
	cat("InSampleCutoff", in_sample_abundance, "\n", file=cdf_values_fh, sep="\t");
	cat("AcrossSampleCutoff", across_sample_prevalence, "\n", file=cdf_values_fh, sep="\t");
	cat("WithinSampleIncrement", file=cdf_values_fh, increm, "\n", sep="\t");

	# Iterate through each zoom
	for(ridx in 1:length(ranges)){

		x_max=ranges[ridx];
		cat("Working on zoom: ", x_max, "\n");

		# Open up empty plot
		plot(0,0,col="white", ylim=c(0,1), xlim=c(0,x_max), xlab="Percent of Organism in Sample", ylab="Percent of Samples with Organism", 
			main=paste("Core: ", InputFileName, sep=""),
			yaxt="n", xaxt="n"
		);

		axis(1,at=seq(0,1,x_max/10) ,labels=seq(0,100,x_max*10));
		axis(2,at=seq(0,1,.1) ,labels=seq(0,100,10));

		# Zoom (in title)
		range_txt=paste("In sample zoom: 0 - ", x_max, sep="");
		mtext(range_txt, side=3, cex=.8, line=.7);

		# Cutoff (in title)
		cutoff_txt=paste("(In sample cutoff = ", in_sample_abundance*100, "%, Across sample cutoff = ", across_sample_prevalence*100, "%)", sep="");
		mtext(cutoff_txt, side=3, cex=.8, line=0);

		# Plot cutoff point
		points(in_sample_abundance, across_sample_prevalence, col="gray", cex=2.3);
		text(in_sample_abundance, across_sample_prevalence, "cutoff", cex=.5);

		# Draw each of the CDF's lines and labels
		num_samples_drawn=0;
		legend_labels=character();
		legend_colors=character();
		for(cat_id in 1:NumCategories){
			
			draw=TRUE;

			if(remove_weak_classifications){
				if(length(grep(":", short_names[cat_id]))>0){
					cat("Not drawing:  ",  short_names[cat_id], "\n", sep="");
					draw=FALSE;
				}
			}

			if(all_cdf_presence[cat_id, anno_idx] < across_sample_prevalence){
				draw=FALSE;
			}else{
				num_samples_drawn=num_samples_drawn+1;
			}
				

			# Only plot categories that exceed the cutoff
			if(draw){

				# Draw CDF
				lines(increm, all_cdf_presence[cat_id,], col=colors[num_samples_drawn]);
				
				# add line to output values
				if(ridx==1){
					cat(short_names[cat_id], "\t", file=cdf_values_fh);
					cat(all_cdf_presence[cat_id, ], "\n", sep="\t", file=cdf_values_fh);
				}

				if(!DrawLegend){
				# Label each line
					text(increm[anno_idx], all_cdf_presence[cat_id, anno_idx], 
						short_names[cat_id], cex=.7, pos=4, col=colors[num_samples_drawn]);
				}else{
					legend_labels[num_samples_drawn]=short_names[cat_id];
					legend_colors[num_samples_drawn]=colors[num_samples_drawn];
				}

				# Keep track of number of lines drawn so we can cycle through the predefined colors
			}
		}

		if(DrawLegend){
			legend(x_max/2, 1, legend=legend_labels, fill=legend_colors, cex=.6);
		}

	}

################################################################################

close(cdf_values_fh);

################################################################################
# Output counts/depths of categories so we can have a heatmap later.

# match the in sample granularity with the cdf granularity

in_sample_values=ncol(all_cdf_presence); # columns
in_sample_dim=29; # columns
across_sample_dim=29; # rows
depth_matrix=matrix(0,nrow=across_sample_dim, ncol=in_sample_dim);

# Sum up area under all curves
for(cat_idx in 1:NumCategories){

	for(in_sample_idx in 1:in_sample_dim){
		across_sample=all_cdf_presence[cat_idx, as.integer(in_sample_idx/in_sample_dim*in_sample_values)];
		across_sample_max_idx=as.integer(across_sample_dim*across_sample);

		if(across_sample_max_idx>across_sample_dim){
			across_sample_max_idx=across_sample_dim;
		}
		if(across_sample_max_idx>0){
			for(across_sample_idx in  1:across_sample_max_idx){
				depth_matrix[across_sample_idx, in_sample_idx]=depth_matrix[across_sample_idx, in_sample_idx]+1;
			}
		}
	}
}

# Compute labels
in_sample_label=sprintf("%i", seq(0,100,length.out=in_sample_dim));
across_sample_label=sprintf("%i",seq(0,100,length.out=across_sample_dim));

# Reverse rows so that bottom left is 0,0, instead of top left
across_sample_label=rev(across_sample_label);
depth_matrix=depth_matrix[across_sample_dim:1,];

# Write out depth table
# Write out column labels
write.table(matrix(c("",in_sample_label), nrow=1, ncol=(in_sample_dim+1)), OutputDepthFileName,  sep="\t", quote=F, row.names=F, col.names=F, append=F);
# Write out data and row labels
write.table(depth_matrix,OutputDepthFileName, sep="\t", quote=F, row.names=across_sample_label, col.names=F, append=T);

################################################################################

writeLines("Done.\n")

q(status=0)
