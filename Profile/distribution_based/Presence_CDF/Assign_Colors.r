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
	"output_file_root", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"	-i <input list of taxa>\n",
	"	[-o <output filename root>]\n",
	"\n",
	"This script will assign colors to the list of taxa\n",
	"which can be used as a global color assignment.\n",
	"\n",
	"The input file is just a single column of taxa (alphanumeric) followed by carriage return.\n",
	"	For example:\n",
	"		apples\\n\n",
	"		bananas\\n\n",
	"		crocodiles\\n\n",
	"		...\n",
	"		zebras\\n\n",
	"\n",
	"The first output file is just a the input file's first column, plus a color assignment in hexadecimal RGB.\n",
	"	For example:\n",
	"		apples\\t#804040\n",
	"		bananas\\t#FF9900\n",
	"		crocodiles\\t#99FF80\n",
	"		...\n",
	"		zebras\\t#008080\n",
	"\n",
	"The second output is a pdf file, which contains a legend, so you can see what color was assigned.\n",
	"\n", sep="");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;

if(length(opt$output_file_root)){
	OutputFileRoot=opt$output_file_root;
}else{
	OutputFileRoot=InputFileName;
}

cat("Input File List: ", InputFileName, "\n");
cat("Output File Root: ", OutputFileRoot, "\n");

################################################################################

# Read in list of taxa
taxa=as.matrix(read.table(InputFileName, sep="\t"))[,1];
#print(taxa);
num_taxa=length(taxa);
cat("Num taxa: ", num_taxa, "\n");

# Step through possible colors
colors=rep(NaN, 250);
counter=0;
for(hue in seq(0,1-(1/10),1/10)){
	if(hue>.12 && hue < .18){ # exclude yellow
		next;
	}
	for(val in seq(.5,1, length.out=2)){
		for(sat in seq(.5,1, length.out=2)){
			counter=counter+1;
			#cat(hue, sat, val, sep=" ", "\n");
			colors[counter]=hsv(hue,sat,val);
		}
	}
}
orig_colors=colors;
cat("Num colors generated: ", counter, "\n", sep="");

# Shuffle the colors so neighbors are not close to each other anymore
dim=ceiling(sqrt(counter));
colors=as.vector(t(matrix((colors[1:(dim^2)]), nrow=dim, ncol=dim)));
colors=colors[colors!="NaN"];

################################################################################

# Output color legend
pdf(paste(OutputFileRoot, ".colormap.pdf", sep=""), height=((num_taxa/5.5)+.5), width=8.5);
#barplot(1:counter, 1, col=orig_colors[1:counter]);
par(oma=c(0,0,0,0), mar=c(0,0,0,0));
plot(0,0, xlim=c(0, 100), ylim=c(0,100), type="n", yaxt="n", xaxt="n");
clean_taxa=gsub("_"," ", taxa);
legend(0,100, legend=clean_taxa, fill=colors, cex=.8);
dev.off();

# Output color assignment as a map
fh=file(paste(OutputFileRoot, ".colormap", sep=""), "w");
for(i in 1:num_taxa){
	col_idx=i%%counter;
	if(col_idx==0){
		col_idx=counter;
	}
	cat(file=fh, taxa[i], "\t", colors[col_idx],"\n", sep="");
}
close(fh);

################################################################################

cat("Done.\n")
quit(status=0)
