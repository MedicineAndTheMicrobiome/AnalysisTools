#!/usr/bin/env Rscript

###############################################################################

library('getopt');
library('seqinr');

params=c(
	"clustal_align", "c", 1, "character",
	"target_seq", "s", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
        "\nUsage:\n\n", script_name,
        "\n",
	"	-c <clustal alignment file>\n",
	"	-s <target sequences to extract>\n",
	"\n",
	"This script will read in an alignment, and then determine\n",
	"the range of positions that are covered by the target sequences.\n",
	"The sub region's alignment across all sequences will then be\n",
	"extracted.\n",
	"\n",
	"For example, if v1v3 regions were aligned to a set of full.\n",
	"length 16s references, then the script will determine the.\n",
	"range of the alignment where v1v3 was relevant, and only.\n",
	"extract out that region.  If v3v5 regions were included,\n",
	"in the alignment, they will not be extracted since they\n",
	"would not overlap the targeted v1v3 sequences.\n",
        "\n", sep="");

if(!length(opt$clustal_align) || !length(opt$target_seq)){
        cat(usage);
        q(status=-1);
}

#----------------------------------------------------------

AlignmentFilename=opt$clustal_align;
TargetSequences=opt$target_seq;

OutputFilenameRoot=paste(gsub(".aln", "", AlignmentFilename));

cat("Alignment File: ", AlignmentFilename, "\n");
cat("Target List: ", TargetSequences, "\n");
cat("Output Filename Root: ", OutputFilenameRoot, "\n");

################################################################################
# Load target list
targets_id_list=scan(TargetSequences, what=character());
num_target_ids=length(targets_id_list);
#print(targets_id_list);
cat("Num target IDs read: ", num_target_ids, "\n");

################################################################################
# Load alignments
aligns=read.alignment(AlignmentFilename, format="clustal");

#aligns$nb   num alignments
#aligns$nam  names
#aligns$seq  gapped sequences
#aligns$com  commends
#print(aligns);

# Get information on the alignments
num_columns=nchar(aligns$seq[[1]]);
num_rows=aligns$nb;
cat("Num sequences aligned: ", aligns$nb, "\n");
cat("Total aligned length: ", num_columns, "\n");

# Determine which positions are residues, ie. not gaps or dots
align_matrix=as.matrix(aligns);
gap_matrix=align_matrix=="-";
dot_matrix=align_matrix==".";
residues=!(gap_matrix | dot_matrix);
#print(residues)

# Get the row indices for the target sequences
target_idx=numeric(num_target_ids);
for(i in 1:num_target_ids){
	target_idx[i]=which(aligns$nam==targets_id_list[i]);
}

# Get the bounds for the residues for the targeted sequences
lb=numeric(num_target_ids);
ub=numeric(num_target_ids);
for(i in 1:num_target_ids){
	residue_range=range(which(residues[target_idx[i],]));
	lb[i]=residue_range[1];
	ub[i]=residue_range[2];
	#cat(aligns$nam[i], ": ", residue_range[1], "-", residue_range[2], "\n", sep="");
}
target_lb=min(lb);
target_ub=max(ub);
range_string=paste(target_lb, "-", target_ub, sep="");
cat("Range for target sequences: ", target_lb, " - ", target_ub, "\n");

################################################################################

cat("Outputing alignments.\n");
# Open output FASTA file
subaln_fh=file(paste(OutputFilenameRoot, ".", range_string, ".fasta", sep=""), "w");
# Open output list of excluded sequences
excluded_fh=file(paste(OutputFilenameRoot, ".", range_string, ".excluded", sep=""), "w");

# Extract out region across all sequences
num_excluded_ids=0;
num_kept_ids=0;
for(i in 1:num_rows){

	# Extract the sequence for that targeted region
	gapped_array=align_matrix[i,  target_lb:target_ub];
	
	# Test for empty string
	gaps=gapped_array=="-";
	dots=gapped_array==".";
	row_residues=which(!(gaps | dots))

	if(length(row_residues)==0){
		# If the alignment for that region is empty, do not output it
		cat(file=excluded_fh, aligns$nam[i], "\n", sep="");
		num_excluded_ids=num_excluded_ids+1;
	}else{
		# If the region is not empty, then output it
		gapped_string=paste(gapped_array, collapse="");
		cat(file=subaln_fh, ">", aligns$nam[i], "\n", sep="");
		cat(file=subaln_fh, gapped_string, "\n", sep="");
		num_kept_ids=num_kept_ids+1;
	}
}

close(subaln_fh);
close(excluded_fh);

cat("\n");
cat("Number of excluded sequences: ", num_excluded_ids, "\n", sep="");
cat("Number of extracted/kept sequences: ", num_kept_ids, "\n", sep="");
cat("\n");

################################################################################

cat("Outputing graphics.\n");
pdf(paste(OutputFilenameRoot, ".", range_string, ".pdf", sep=""), height=8.5, width=11);

out_map=matrix(1, nrow=num_rows, ncol=num_columns);
for(i in 1:num_rows){
	if(any(i==target_idx)){
		color=1;
	}else{
		color=.5;
	}
	for(j in 1:num_columns){
		if(residues[i,j]){
			out_map[i,j]=color;
		}else{
			out_map[i,j]=0;
		}
	}
}
image(t(out_map[num_rows:1,]), col=c("white", "black", "red"),
	xaxt="n", yaxt="n",
	ylab="sequences", xlab="residues",
	main="Targeted Multiple Sequence Alignments"
	);

axis(side=2, at=seq(0,1, length.out=10), labels=rev(round(seq(1,num_rows, length.out=10))), las=2, cex.axis=.5);
axis(side=1, at=seq(0,1, length.out=10), labels=round(seq(1,num_columns, length.out=10)), las=2, cex.axis=.5);

dev.off();

################################################################################

cat("Done.\n");
