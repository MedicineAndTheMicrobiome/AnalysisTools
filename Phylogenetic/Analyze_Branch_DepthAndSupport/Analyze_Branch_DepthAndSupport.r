#!/usr/bin/env Rscript

###############################################################################

library('getopt');
library('ape');

PAPER_WIDTH=8.5;
PAPER_HEIGHT=-1;

params=c(
	"tree_fn", "t", 1, "character",
	"output_fn_root", "o", 2, "character",
	"reference_accession", "r", 2, "character",
	"height", "H", 2, "numeric",
	"width", "W", 2, "numeric",
	"color_fn", "c", 2, "character",
	"description_fn", "d", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
        "\nUsage:\n\n", script_name,
        "\n",
	"	-t <tree file name>\n",
	"	[-o <output filename root>]\n",
	"	[-r <reference accession>]\n",
	"	[-H <paper height, default is dynamic>]\n",
	"	[-W <paper width, default is 8.5>]\n",
	"	[-c <color map>]\n",
	"	[-d <description map>]\n",
	"\n",
	"This script will read in a bootstrapped trees in Newick format, and\n",
	"then report some results.\n",
	"\n",
	"	1.) Generate boxplots for distances from reference.\n",
	"	2.) Generate text output for 95% CI from reference.\n",
	"	3.) Generate consensus tree topologies from .5 to 1.0 support.\n",
        "\n", sep="");

if(!length(opt$tree_fn)){
        cat(usage);
        q(status=-1);
}

#----------------------------------------------------------

if(length(opt$tree_fn)){
	TreesFilename=opt$tree_fn;
}

if(length(opt$output_fn_root)){
	OutputRoot=opt$output_fn_root;
}else{
	OutputRoot=TreesFilename;
}

if(length(opt$reference_accession)){
	ReferenceAccession=opt$reference_accession;
}else{
	ReferenceAccession="";
}

if(length(opt$height)){
	Height=opt$height;
}else{
	Height=-1;
}

if(length(opt$width)){
	Width=opt$width;
}else{
	Width=8.5;
}

if(length(opt$color_fn)){
	ColorMapFname=opt$color_fn;
}else{
	ColorMapFname="";
}

if(length(opt$description_fn)){
	DescriptionMapFname=opt$description_fn;
}else{
	DescriptionMapFname="";
}

###############################################################################
# Functions for loading accession to color or description maps

load_map=function(map_filename){
        if(map_filename!=""){
                mat=as.matrix(read.delim(map_filename, sep="\t", header=FALSE));
                map=list();
                for(i in 1:nrow(mat)){
                        map[[mat[i,1]]]=mat[i,2];
                }
                return(map);
        }else{
                return(NULL);
        }
}

rename=function(accessions, map, keep_accession=TRUE, default=""){

        num_accessions=length(accessions);
        new_names=character(num_accessions);
        for(i in 1:num_accessions){

                cur_acc=accessions[i];
                new_name=map[[cur_acc]];
                if(length(new_name)>0){
                        if(keep_accession){
                                new_names[i]=paste(cur_acc, ": ", new_name, sep="");
                        }else{
                                new_names[i]=new_name;
                        }
                }else{
                        if(default==""){
                                new_names[i]=cur_acc;
                        }else{
                                new_names[i]=default;
                        }
                }

        }
        return(new_names);

}

description_map=load_map(DescriptionMapFname);
#print(description_map);
color_map=load_map(ColorMapFname);
#print(color_map);

###############################################################################
# Load tree from file
cat("Loading trees: ", TreesFilename, "\n");
trees=read.tree(file=TreesFilename);
num_trees=length(trees);

cat("Number of trees loaded: ", num_trees, "\n");
if(num_trees==1){
	cat("Error:  Only one tree was specified.\n");
	quit(status=-1);
}

sorted_sample_names=sort(trees[[1]]$tip.label);
num_samples=length(sorted_sample_names);

cat("Tree leaf names:\n");
print(sorted_sample_names);

cat("Num samples: ", num_samples, "\n");

################################################################################
# Open PDF file

if(Height==-1){
        Height=1/10*num_samples;
}
cat("Height Used: ", Height, "\n");

pdf(paste(OutputRoot, ".pdf", sep=""), height=Height, width=Width);

################################################################################
# Compute cophenetic distance matrix

coph_mats=array(0, dim=c(num_samples, num_samples, num_trees));
for(i in 1:num_trees){
	#plot(trees[[i]]);
	coph_mats[,,i]=cophenetic(trees[[i]])[sorted_sample_names,sorted_sample_names];
}

# Compute bound indices
alpha=.05;
ubix=ceiling(num_trees*(1-alpha/2));
lbix=max(1,floor(num_trees*(alpha/2)))+1;

cat("Alpha: ", alpha, "\n");
cat("Index Bounds: (", lbix, ", ", ubix, ")\n");

max_dist=max(coph_mats);

# Store bootstrapped statistics
median=matrix(0, nrow=num_samples, ncol=num_samples);
lb=matrix(0, nrow=num_samples, ncol=num_samples);
ub=matrix(0, nrow=num_samples, ncol=num_samples);

# Go through each sample pairwise to compute median and bounds
for(i in 1:num_samples){
	cat("Working on: ", sorted_sample_names[i], "\n");
	for(j in i:num_samples){
		sorted=sort(coph_mats[i,j,]);
		median[i,j]=median(sorted);
		lb[i,j]=sorted[lbix];
		ub[i,j]=sorted[ubix];

		median[j,i]=median[i,j];
		lb[j,i]=lb[i,j];
		ub[j,i]=ub[i,j]
		#hist(sorted, xlim=c(0, 10), ylim=c(0, 30), breaks=seq(0, 10, length.out=20));
	}
}

################################################################################
# Get distance from reference

if(ReferenceAccession!=""){
	reference_idx=which(ReferenceAccession==sorted_sample_names);
	cat("Reference Accession: ", ReferenceAccession, "\n");
	cat("Reference Index: ", reference_idx, "\n");
	
	if(length(reference_idx)==0){
		cat("Error: Reference accession specified (", ReferenceAccession, ").  Not found.\n", sep="");
		quit(status=-1);
	}

	#----------------------------------------------------------------------
	# Output distances to reference into text file, and report 95% CI

	fc=file(paste(OutputRoot, ".dist_to_ref.csv", sep=""), "w");
	to_ref_med=median[reference_idx,];
	to_ref_lb=lb[reference_idx,];
	to_ref_ub=ub[reference_idx,];
	med_sort_ix=order(to_ref_med);

	cat(file=fc, paste(
		"Sample Name",
		"Median",
		"LowerBound",
		"UpperBound",
		sep=","),
		"\n", sep="");
	
	for(i in med_sort_ix){
		cat(file=fc,
			sorted_sample_names[i],
			to_ref_med[i],
			to_ref_lb[i],
			to_ref_ub[i],
			sep=","
		)
		cat(file=fc, "\n");
	}
	close(fc);

	#----------------------------------------------------------------------
	# Perpare date for generating boxplots

	# Make matrix of bootstrap values to reference
	dist_to_ref=matrix(0, ncol=num_samples, nrow=num_trees);
	for(i in 1:num_samples){
		dist_to_ref[,i]=coph_mats[i, reference_idx,];	
	}

	ssns=sorted_sample_names[rev(med_sort_ix)];
	# Rename accessions to descriptions
	if(!is.null(description_map)){
		boxplot_names=rename(ssns, description_map, keep_accession=TRUE);
	}

	# Generate colors for each accession
	if(!is.null(color_map)){
		colors=rename(ssns, color_map, keep_accession=FALSE, default="grey");
	}else{
		colors=rep(1, length(boxplot_names));
	}

	# Generate boxplot
	par(mar=c(5,25,1,1));
	par(oma=c(0,0,3,0));
	boxplot(dist_to_ref[,rev(med_sort_ix)], horizontal=T, cex=.5, 
		names=boxplot_names, las=2, cex.axis=.6,
		col=colors,
		xlim=c(num_samples*.035,num_samples*(1-.035))
		);

	# Label
	xaxp=par()$xaxp;
	ticks=seq(xaxp[1], xaxp[2], length.out=xaxp[3]+1);
	abline(v=ticks, col="pink", lwd=1, lty=3);
	mtext("Ordered Cophenetic Distances to Reference");
	
}

################################################################################
# Sum Trees and output

par(mfrow=c(1,1));
par(mar=c(5.1,1,1,1));
par(oma=c(0,0,2,0));

cat("Computing consensus tree...\n");

# Specify which cutoffs to compute
# p=.5, majority rule
proportion_cutoffs=seq(.5,1,.1);

# Generate plot at each cutoff
for(p in proportion_cutoffs){
	cat("Working on: ", p, "\n");
	cons_tree=consensus(trees, p=p);
	num_tips=length(cons_tree$tip.label);
	num_int_nodes=cons_tree$Nnode
	edges_dist=rep(1, num_tips+num_int_nodes-1);
	cons_tree$edge.length=edges_dist;
	if(ReferenceAccession!=""){
		cons_tree=root(cons_tree, ReferenceAccession);
		cat("Tree rooted with: ", ReferenceAccession, "\n", sep="");
	}

	# Descriptions
	tiplabels=cons_tree$tip.label;
	if(!is.null(description_map)){
		cons_tree$tip.label=rename(tiplabels, description_map, keep_accession=TRUE);
	}

	# Colors
	if(!is.null(color_map)){
		colors=rename(tiplabels, color_map, keep_accession=FALSE, default="black");
	}else{
		colors=1;
	}

	# Plot
	plot(cons_tree, cex=.5, tip.color=colors);
	mtext(sprintf("Support Cutoff = %2.2f",p), cex=1.5, font=2);
}

################################################################################

dev.off()

cat("Done.\n");
