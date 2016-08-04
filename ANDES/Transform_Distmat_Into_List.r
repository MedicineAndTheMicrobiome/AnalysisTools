#!/usr/bin/env Rscript

###############################################################################
#                                                                             # 
#       Copyright (c) 2011 J. Craig Venter Institute.                         #     
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
library(MASS);

params=c(
	"input_file", "i", 1, "character",
	"hide_filename", "f", 2, "logical",
	"standard_deviation", "s", 2, "numeric",
	"output_filename", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage <- paste (
		"\nUsage:\n", script_name, "\n",
		"	-i <input distance matrix>\n",
		"	[-f (hide filename in output option)]\n",
		"	[-s <standard deviation>]\n",
		"	[-o <output file name root>]\n",
		"\n",
		"Reads in a distance matrix, sorts its pairs by increasing distance.\n",
		"By default, the first column is the source file name.\n",
		"If you don't want that output in the first column, then use the -f\n",
		"option to suppress that.\n",
		"\n",
		"The second and third column are the two objects that the 4th column's\n",
		"distance is referring to.\n",
		"If the -s option is specified, the 5th column is the standard deviation\n",
		"of the distances for that distance matrix.  It will be used for all the\n",
		"the distances.\n",
		"\n",
		"Example:\n",
		"\n",
		"    [filename]\\t<ObjectA>\\t<ObjectB>\\t<Distance>\\t[standard deviation]\\n\n",
		"    myDist\tNewYork\tNewJersey\t66.8\t6.18\n",
		"\n",
		"The output file name will be based on the input file, unless it is specified\n",
		"by the -o option.\n",
		"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileRoot=InputFileName;
OutputFileRoot=gsub("\\.r_distmat$","", OutputFileRoot);
OutputFileRoot=gsub("\\.distmat$","", OutputFileRoot);
OutputFileRoot=gsub("\\.dist_mat$","", OutputFileRoot);

path=strsplit(OutputFileRoot, "\\/");
infname=tail(path[[1]],1);

if(length(opt$output_filename)){
	OutputFileRoot=opt$output_filename;
}

HideFilename=FALSE;
if(length(opt$hide_filename)){
	HideFilename=TRUE;
}

StandardDeviation=-1;
if(length(opt$standard_deviation)){
	StandardDeviation=opt$standard_deviation;
}

###############################################################################


load_distmat=function(distmat_fname){

	seps=c(",", "\t", " ");

	for(i in 1:length(seps)){

		dist_mat=as.matrix(read.table(InputFileName, header=TRUE, sep=seps[i], row.names=1, check.names=FALSE));

		nrow=nrow(dist_mat);
		ncol=ncol(dist_mat);

		if(nrow==ncol){
			cat("Using: '", seps[i], "' as separator.\n", sep="");
			break;
		}
	}

	if(nrow!=ncol){
		cat("rows: ", nrow, "  cols: ", ncol, "\n");
		cat("Error:  Your distance matrix is not square.\n");
		quit(status=-1);
	}

	cat("Dist Matrix Dimensions: ", nrow, "(r) x ", ncol, "(c)\n");
	return(dist_mat);
}

dist_mat=load_distmat(distmat_fname);
nrow=nrow(dist_mat);
ncol=ncol(dist_mat);
colnames=colnames(dist_mat);

###############################################################################

src_name=character();
dst_name=character();
distances=numeric();
num_dist=1;

for(i in 1:nrow){
	cat(".");
	for(j in 1:nrow){
		if(i<=j){
			next;
		}
		distances[num_dist]=dist_mat[i,j];
		src_name[num_dist]=colnames[i];
		dst_name[num_dist]=colnames[j];
		num_dist=num_dist+1;
	}
}
cat("\n");

###############################################################################

sort_info=sort(distances, index.return=1, decreasing=FALSE);

fh=file(paste(OutputFileRoot, ".dist_flat", sep=""), "w");

for(i in sort_info$ix){

	# Build output line
	line=character();
	if(!HideFilename){
		line=c(line, infname);
	}
	
	line=c(line, src_name[i], dst_name[i], distances[i]);

	if(StandardDeviation>=0){
		line=c(line, StandardDeviation);
	}

	cat(file=fh, paste(line, collapse="\t"), "\n", sep="");
}

close(fh);

###############################################################################

cat("Printing warnings.\n");
print(warnings());
cat("Done.\n");
