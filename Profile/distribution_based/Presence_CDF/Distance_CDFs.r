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
	"input_filename", "i", 1, "character",
	"output_file_root", "o", 2, "character",
	"ubiquity_cutoff", "u", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

UbiquityCutoff=.2;

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"	-i <input filename cdf>\n",
	"	[-o <output filename root>]\n",
	"	[-u <ubiquity cutoff, default=", UbiquityCutoff, ">]\n",
	"\n",
	"Reads in a presence CDF file and computes a pairwise distance\n",
	"between each line.  Generates a distance matrix.\n",
	"\n", sep="");

if(!length(opt$input_filename)){ 
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_filename;
OutputFileRoot=gsub(".cdf", "", InputFileName);

if(length(opt$output_file_root)){
	OutputFileRoot=opt$output_file_root;	
}

if(length(opt$ubiquity_cutoff)){
	UbiquityCutoff=opt$ubiquity_cutoff;
}	


cat("Input File Name: ", InputFileName, "\n");
cat("Output File Root: ", OutputFileRoot, "\n");
cat("Ubiquity Cutoff: ", UbiquityCutoff, "\n");

################################################################################

read_cdf_from_file=function(cdf_fn){
	cdf_info=list();
	cdfs=as.matrix(read.table(cdf_fn, sep=",", header=TRUE, row.names=1, check.names=FALSE));
	
	cdf_info$fname=cdf_fn;
	cdf_info$cutoffs=as.numeric(colnames(cdfs));
	cdf_info$taxa=rownames(cdfs);
	cdf_info$cdf_val=cdfs;
	cdf_info$num_taxa=nrow(cdfs);
	
	return(cdf_info);
}

write_distmat_to_file=function(dist, fh){
	dist=as.matrix(dist);
	nrow=nrow(dist);
	taxa=rownames(dist);
	cat(file=fh, " ", paste(taxa, collapse=" ", sep=" "), "\n", sep="");
	for(i in 1:nrow){
		cat(file=fh, taxa[i], " ", paste(dist[i,], collapse=" "), "\n", sep="");	
	}
}


max_diff=function(a,b){
	diff=max(abs(a-b))
	return(diff);
}


max_diff_dist=function(mat){
	num_samples=nrow(mat);
	distmat=matrix(0,nrow=num_samples, ncol=num_samples);
	for(i in 1:num_samples){
		for(j in 1:i){
			distmat[i,j]=max_diff(mat[i,], mat[j,]);
		}
	}
	colnames(distmat)=rownames(mat);
	rownames(distmat)=rownames(mat);
	return(as.dist(distmat));
}

################################################################################

cdf_info=read_cdf_from_file(InputFileName);

cat("Num Taxa: ", cdf_info$num_taxa, "\n");

taxa_to_keep=which(cdf_info$cdf_val[,2]>=UbiquityCutoff);
num_kept=length(taxa_to_keep);

cat("Num Taxa exceeding ubiquity cutoff: ", num_kept, "\n");

#distances=dist(cdf_info$cdf_val[taxa_to_keep,]);
distances=max_diff_dist(cdf_info$cdf_val[taxa_to_keep,]);

################################################################################
# Output results

fh=file(paste(OutputFileRoot, ".dist_mat", sep=""), "w");
write_distmat_to_file(distances, fh);
close(fh);

################################################################################

cat("Done.\n")

q(status=0)
