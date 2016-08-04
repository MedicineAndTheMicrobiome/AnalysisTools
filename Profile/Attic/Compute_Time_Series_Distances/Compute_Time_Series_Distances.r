#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"association_file", "a", 1, "character",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	-a <input association file>\n",
	"	[-o <output file name>]\n",
	"\n",	
	"This script will compute the distance between consecutive\n",
	"samples specified by the association file.\n",
	"\n",
	"An output file will be generated for each consecutive pair.\n",
	"\n",
	"\n");

if(!length(opt$input_file)||!(length(opt$association_file))){
	cat(usage);
	q(status=-1);
}

###############################################################################

# Figure out where we are running, this script so we can find the "WeightedRankDifference.r" code.
path_comp=strsplit(script_name, "/")[[1]];
bin_comp_idx=length(path_comp);
bin_path=paste(path_comp[-bin_comp_idx], collapse="/", sep="");
if(nchar(bin_path)==0){
        bin_path=".";
}
cat("Binary path: '", bin_path, "'\n", sep="");
source(paste(bin_path, "WeightedRankDifference.r", sep="/"));

###############################################################################

InputFileName=opt$input_file;
OutputFilename=opt$output_file;
AssociationFile=opt$association_file;

if(length(opt$output_file)==0){
	OutputFilenameRoot=paste(gsub("\\.summary_table\\.xls$", "", InputFileName));
}

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name Root: ", OutputFilenameRoot, "\n");
cat("Association File Name: ", AssociationFile, "\n");
cat("\n");

###############################################################################
###############################################################################

load_summary_file_table=function(input_filename){
        # Load data
        cat("Loading: ", input_filename, "\n");
        inmat=as.matrix(read.table(input_filename, sep="\t", header=TRUE, check.names=FALSE, comment.char="*", row.names=1))
        #cat("Original Matrix:\n")
        #print(inmat);

        # Grab columns we need into a vector, ignore totals, we won't trust it.
        counts_mat=inmat[,2:(ncol(inmat))];
        #print(counts_mat);

        # Summary what we've loaded
        num_samples=nrow(counts_mat);
        num_categories=ncol(counts_mat);
        cat("Num Samples: ", num_samples, "\n");
        cat("Num Categories: ", num_categories, "\n");

        return(counts_mat);
}

#-----------------------------------------------------------------------------

normalize=function(inmat){
        samp_sums=apply(inmat, 1, sum);
        norm_mat=matrix(0, nrow=nrow(inmat), ncol=ncol(inmat));
        for(i in 1:length(samp_sums)){
                norm_mat[i,]=inmat[i,]/samp_sums[i];
        }

        rownames(norm_mat)=rownames(inmat);
        colnames(norm_mat)=colnames(inmat);

        #cat("Normalized:\n");
        #print(norm_mat);
        return(norm_mat);

}

#-----------------------------------------------------------------------------

load_associations=function(input_filename){
        cat("Loading association file: ", input_filename, "\n", sep="");

        inmat=as.matrix(read.table(input_filename, sep=",", header=TRUE, check.names=FALSE, row.names=1))

        associations=list();
        associations$sets=inmat;
        associations$num_groups=nrow(inmat);
        associations$num_types=ncol(inmat);
        associations$type_names=colnames(inmat);
        associations$group_ids=rownames(inmat);

        return(associations);
}

###############################################################################

# Load associations
associations=load_associations(AssociationFile);
cat("Num Time Points: ", associations$num_types, "\n");
cat("Num Time Sequences: ", associations$num_groups, "\n");
#print(associations);

# Load data
in_matrix=load_summary_file_table(InputFileName);
#cat("Original Matrix:\n")
#print(in_matrix);
category_names=colnames(in_matrix);
sample_names=rownames(in_matrix);
num_samples=nrow(in_matrix);
num_categories=ncol(in_matrix);

# Normalize
normalized=normalize(in_matrix);
#cat("Normalized Matrix:\n");
#print(normalized);

for(timepoints in 1:(associations$num_types-1)){

	cat("Working on: ", associations$type_names[timepoints], " to ", associations$type_names[timepoints+1], "\n", sep="");
	fc=file(paste(OutputFilenameRoot, ".", associations$type_names[timepoints], ".dist", sep=""), "w");
	cat("sample_id,from_", associations$type_names[timepoints], "\n", sep="", file=fc);
	for(sequences in 1:associations$num_groups){
		pt1=associations$sets[sequences, timepoints];
		pt2=associations$sets[sequences, timepoints+1];

		#cat("Computing distance from ", pt1, " to ", pt2, "\n", sep="");
		pt1_distrib=normalized[pt1,];
		pt2_distrib=normalized[pt2,];
		dist=order_dist(pt1_distrib, pt2_distrib, deg=4);
		cat(pt1, ",", dist, "\n", sep="", file=fc);
	}
	close(fc);
}

###############################################################################
#fc=file(OutputFilenameCSV, "w");

cat("Done.\n")
#print(warnings());

q(status=0)
