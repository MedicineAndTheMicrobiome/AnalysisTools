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
	"input_file_list", "i", 1, "character",
	"output_file_root", "o", 2, "character",
	"average_curves", "a", 2, "logical",
	"multiply_curves", "m", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"	-i <list of input cdf files, comma separated, no spaces>\n",
	"	-o <output filename root>\n",
	"	[-a (flag for averaging curves together)]\n",
	"	[-m (flag for multiply (ANDing) curves together)]\n",
	"\n",
	"This script will take in a list of CDF files (-i) and generate a combined\n",
	"CDF file by either averaging (-a) or multiplying (-m) the probabilities at each\n",
	"abundance cutoff together across all the specified cdf files in order to arrive\n",
	"at new combined ubiquity values across all abundances.\n",
	"\n",
	"This script is used for combining the CDF values across multiple environments\n",
	"in order to get a composite analysis of what core may be.  For example, if you\n",
	"have 5 skin location, you can combine them together to get a single skin CDF file.\n",
	"which you may compare against a oral region composed of 10 oral locations.  This way,\n",
	"you won't have to worry as much about the fact that you have twice as many oral\n",
	"locations that might bias your analysis towards the core of your oral regions.\n",
	"\n",

	"\n", sep="");

if(!length(opt$input_file_list) && !length(opt$output_file_root)){
	cat(usage);
	q(status=-1);
}

InputFileNameList=opt$input_file_list;
OutputFileRoot=opt$output_file_root;

AverageCurves=FALSE;
if(length(opt$average_curves)){
	AverageCurves=TRUE;	
}

MultiplyCurves=FALSE;
if(length(opt$multiply_curves)){
	MultiplyCurves=TRUE;	
}

if(MultiplyCurves & AverageCurves){
	cat("Error, you can't perform both operations.\n");
	quit(status=-1);
}

if(!(MultiplyCurves | AverageCurves)){
	cat("Error, you have to pick an operation.\n");
	quit(status=-1);
}

cat("Input File List: ", InputFileNameList, "\n");
cat("Output File Root: ", OutputFileRoot, "\n");
cat("Multiply? ", MultiplyCurves, "\n");
cat("Average? ", AverageCurves, "\n");


################################################################################
# Split input filename list string into array
input_files_vector=strsplit(InputFileNameList, ",")[[1]];
num_input_cdfs=length(input_files_vector);
cat("Num input CDFs in list: ", num_input_cdfs, "\n");
for(i in 1:num_input_cdfs){
	cat("\t", input_files_vector[i], "\n");
}

################################################################################
read_cdf_from_file=function(cdf_fn){
	cdf_info=list();
	cdfs=as.matrix(read.table(cdf_fn, sep=",", header=TRUE, row.names=1, check.names=FALSE));
	
	cdf_info$fname=cdf_fn;
	cdf_info$cutoffs=as.numeric(colnames(cdfs));
	cdf_info$taxa=rownames(cdfs);
	cdf_info$cdf_val=cdfs;
	
	return(cdf_info);
}

#-------------------------------------------------------------------------------

write_cdf_to_file=function(matrix, fh){
	ncol=ncol(matrix);
	nrow=nrow(matrix);
	increments=colnames(matrix);
	taxa=rownames(matrix);
	cat(file=fh, ",");
	cat(file=fh, paste(sprintf("%s",increments), collapse=","), "\n");
#	cat(file=fh,paste(increments, collapse=","), "\n");
	for(i in 1:nrow){
		cat(file=fh, taxa[i], paste(matrix[i,], collapse=","), sep=",");
		cat(file=fh,"\n");
	}	
}

################################################################################
# Load files into list
cdf_list=list()
for(i in 1:num_input_cdfs){
	cur_filename=input_files_vector[i];
	cdf_list[[i]]=read_cdf_from_file(input_files_vector[i]);
}

#-------------------------------------------------------------------------------
# Confirm that all inputs have the same set of taxa and cutoffs
cutoffs=cdf_list[[1]]$cutoffs;
taxa_list=cdf_list[[1]]$taxa;
for(i in 2:num_input_cdfs){
	if(!all(cdf_list[[i]]$cutoffs==cutoffs)){
		cat("Error! Cutoffs don't match.\n");
		quit(status=-1);
	}
	if(!all(cdf_list[[i]]$taxa==taxa_list)){
		cat("Error! Taxas don't match.\n");
		quit(status=-1);
	}
}
num_taxa=length(taxa_list);
num_cutoffs=length(cutoffs);

#-------------------------------------------------------------------------------
# Perform computations

if(AverageCurves){
	accumulator=matrix(0, nrow=num_taxa, ncol=num_cutoffs);
	for(i in 1:num_input_cdfs){
		cur_cdf=cdf_list[[i]]$cdf_val;
		for(j in 1:num_cutoffs){
			for(i in 1:num_taxa){
				accumulator[i,j]=accumulator[i,j]+cur_cdf[i,j];
			}
		}
	}
	accumulator=accumulator/num_input_cdfs;
	ext="avg";
}

if(MultiplyCurves){
	accumulator=matrix(1, nrow=num_taxa, ncol=num_cutoffs);
	for(i in 1:num_input_cdfs){
		cur_cdf=cdf_list[[i]]$cdf_val;
		for(j in 1:num_cutoffs){
			for(i in 1:num_taxa){
				accumulator[i,j]=accumulator[i,j]*cur_cdf[i,j];
			}
		}
	}
	ext="and";
}

colnames(accumulator)=sprintf("%1.16f",cutoffs);
rownames(accumulator)=taxa_list;

################################################################################
# Output results

fh=file(paste(OutputFileRoot, ".", ext, ".cdf", sep=""), "w");
write_cdf_to_file(accumulator, fh);
close(fh);

################################################################################

writeLines("Done.\n")

q(status=0)
