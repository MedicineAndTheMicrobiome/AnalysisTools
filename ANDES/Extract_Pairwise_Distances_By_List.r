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
		"matrix_list", "m", 1, "character",
		"subject_list", "s", 1, "character",
		"extracted", "e", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage <- paste (
		"\nUsage:\n", script_name, "\n",
		"	-m <matrix list filename>\n",
		"	-s <subject list filename>\n",
		"	-e <extracted output filename>\n",
		"\n",
		"Given a list of subjects and a list of matrices, will pull\n",
		"out of the distance matrices all the distances among the subjects.\n",
		"\n");

if(!length(opt$matrix_list) || !length(opt$subject_list) || !length(opt$extracted)){
	cat(usage);
	q(status=-1);
}

MatrixList=opt$matrix_list;
SubjectList=opt$subject_list;
ExtractedOut=opt$extracted;

###############################################################################

load_distmat=function(distmat_fname){

	dist_mat=as.matrix(read.table(distmat_fname, header=TRUE, check.names=FALSE));

	nrow=nrow(dist_mat);
	ncol=ncol(dist_mat);

	cat("Dist Matrix Dimensions: ", nrow, "(r) x ", ncol, "(c)\n");

	if(nrow!=ncol){
		cat("Error:  Your distance matrix is not square.\n");
		quit(status=-1);
	}

	return(dist_mat);
}

###############################################################################

load_list=function(list_fname){
	list=as.matrix(read.table(list_fname, header=FALSE));
	return(as.vector(list[,1]));
}

###############################################################################

# Load subjects to extract
subject_list=load_list(SubjectList);
subject_list_length=length(subject_list);
cat("Subject List Length: ", subject_list_length, "\n");
print(subject_list);
cat("\n");

# Load list of matrices to extract from
matrix_list=load_list(MatrixList);
matrix_list_length=length(matrix_list);
cat("Matrix List Length: ", matrix_list_length, "\n");
print(matrix_list);

cat("\n\n");

###############################################################################

src=character();
dst=character();
dist=numeric();
dist_counter=1;
for(m_idx in 1:matrix_list_length){		
	
	# Load matrix
	cat("Loading: ", matrix_list[m_idx], "\n");
	distmat=load_distmat(matrix_list[m_idx]);
	matrix_subjects=colnames(distmat);
	mat_dim=ncol(distmat);

	# Determine which subjects are common between matrix and list
	common_subjects=intersect(subject_list, matrix_subjects);
	num_common=length(common_subjects);
	cat("Num Common: ", num_common, "\n");

	# Extract out only subjects in list
	common_matrix=distmat[common_subjects, common_subjects];
	distances_to_extract=((num_common*(num_common-1))/2);
	cat("Num Dist to extract: ", distances_to_extract, "\n");

	# Store bottom triangle
	for(i in 1:num_common){
		for(j in i:num_common){
			if(i<=j){
				dist[dist_counter]=common_matrix[i,j];
				src[dist_counter]=common_subjects[j];
				dst[dist_counter]=common_subjects[i];
				dist_counter=dist_counter+1;
			}
		}
	}

	cat("\n");

}

###############################################################################

fh=file(ExtractedOut, "w");

for(i in 1:(dist_counter-1)){
	cat(file=fh, paste(src[i], dst[i], dist[i], sep="\t"), "\n", sep="");
}

close(fh);
###############################################################################

cat("Printing warnings.\n");
print(warnings());
cat("Done.\n");
