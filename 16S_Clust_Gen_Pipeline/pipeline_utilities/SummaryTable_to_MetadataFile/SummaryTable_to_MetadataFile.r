#!/usr/bin/env Rscript

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_file_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary_table.tsv>\n",
	"	-o <output filename root>\n",
	"\n",	
	"This script will read in the summary table and generate\n",
	"metadata/factor file based on the sample ID.\n",
	"\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

###############################################################################

InputFileName=opt$input_file;
OutputFileNameRoot=opt$output_file_root;

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name Root: ", OutputFileNameRoot, "\n");

###############################################################################
###############################################################################

project_ids_map=list(
	"0000"="MoBio.Powersoil.DNA.Ext.Neg",
	"0001"="PCR.Negative",
	"0002"="PCR.Positive.Zymo", 
	"0003"="PCR.Positive.Ecoli.DH5a",
	"0004"="Stool.Test.Extraction",
	"0005"="Qiagen.FastStool.DNA.Ext.Neg",
	"0006"="MoBio.Ultraclean.DNA.Ext.Neg",
	"0007"="PCR.Positive.Other",
	"0008"="Vehicle.Control.Saline",
	"0009"="Vehicle.Control.Swab",
	"0010"="MoBio.PowerMicrobiome.DNA.Ext.Neg",
	"0011"="Zymo.IntctCells.Pos.Ctrl.For.Extr",
	"0012"="Ecoli.IntctCells.Pos.Ctrl.For.Extr",
	"0013"="Pos.Ctrl.For.Extr.Other",
	"0015"="Neg.Ctrl.Genomic.DNA.by.Invstgtr"
);

cat("Mapped Project IDs:\n");
print(project_ids_map);

###############################################################################
# Load data

inmat=as.matrix(read.delim(InputFileName, sep="\t", header=TRUE, 
	check.names=FALSE, comment.char="", quote="", row.names=1))

sample_names=rownames(inmat);
num_samples=nrow(inmat);
read_depth=inmat[,1];

cat("Sample Names:\n");
print(sample_names);

samp_name_split=strsplit(sample_names, "\\.");

# Translate project ID to description
proj_ids=character(num_samples);
proj_desc=character(num_samples);
for(i in 1:num_samples){
	proj_ids[i]=samp_name_split[[i]][1];
	desc=project_ids_map[[proj_ids[i]]];

	if(!is.null(desc)){
		proj_desc[i]=project_ids_map[[proj_ids[i]]];
	}else{
		proj_desc[i]=paste("P_", proj_ids[i], sep="");
	}
}
names(proj_desc)=proj_ids;

cat("\nProject IDs:\n");
print(proj_desc);

###############################################################################
# Output metadata

outfh=file(paste(OutputFileNameRoot, ".metadata.tsv", sep=""), "w");

cat(file=outfh, "Sample_ID\tSample_Type\tRead_Depth\n");

for(row_ix in 1:num_samples){
	cat(file=outfh, 
		paste(sample_names[row_ix], proj_desc[row_ix], read_depth[row_ix], sep="\t"), "\n", sep="");
}

close(outfh);

###############################################################################

cat("Done.\n")
print(warnings());

q(status=0)
