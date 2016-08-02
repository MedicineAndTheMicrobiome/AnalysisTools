#!/usr/bin/env Rscript

###############################################################################

library('getopt');
library('ape');

PAPER_WIDTH=11;
PAPER_HEIGHT=-1;

params=c(
	"tree_fn", "t", 1, "character",

	"output_fn_root", "o", 2, "character",
	"color_fn", "c", 2, "character",
	"description_fn", "d", 2, "character",
	"width", "w", 2, "numeric",
	"height", "h", 2, "numeric",
	"root_acc", "r", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
        "\nUsage:\n\n", script_name,
        "\n",
	"	-t <tree file name>\n",
	"	[-o <output filename root>]\n",
	"	[-c <color map>]\n",
	"	[-d <description map>]\n",
	"	[-w <paper width, default=", PAPER_WIDTH, ">]\n",
	"	[-h <paper height, default is dynamic>]\n",
	"	[-r <accession with which to root the tree>]\n",
	"\n",
	"This script will read in a tree in Newick or Nexus format, and\n",
	"then generate a plot using the paper size specified\n",
	"and the colors/descriptions.\n",
	"\n",
	"The tree may be rooted with the -r option.\n",
	"\n",
	"This script plots only one tree.  If you have several trees\n",
	"from bootstrapping, then a consensus tree (eg. from sumtrees.py)\n",
	"needs to be generated first.\n",
        "\n", sep="");

if(!length(opt$tree_fn)){
        cat(usage);
        q(status=-1);
}

#----------------------------------------------------------

ColorFilename="";
DescriptionFilename="";

if(length(opt$tree_fn)){
	TreeFilename=opt$tree_fn;
}

if(length(opt$output_fn_root)){
	OutputRoot=opt$output_fn_root;
}else{
	OutputRoot=TreeFilename;
}

if(length(opt$width)){
	Width=opt$width;
}else{
	Width=PAPER_WIDTH;
}

if(length(opt$height)){
	Height=opt$height;
}else{
	Height=-1;
}

if(length(opt$description_fn)){
	DescriptionMapFname=opt$description_fn;
}else{
	DescriptionMapFname="";
}

if(length(opt$color_fn)){
	ColorMapFname=opt$color_fn;
}else{
	ColorMapFname="";
}

if(length(opt$root_acc)){
	RootAcc=opt$root_acc;
}else{
	RootAcc="";
}

###############################################################################

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

###############################################################################

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

###############################################################################

description_map=load_map(DescriptionMapFname);
#print(description_map);
color_map=load_map(ColorMapFname);
#print(color_map);

# Load tree from file

magic_number=scan(TreeFilename, what="character", nmax=1);
if(magic_number=="#NEXUS"){
	cat("Reading tree as Nexus format.\n");
	tree=read.nexus(file=TreeFilename);
}else{
	cat("Reading tree as Newick format.\n");
	tree=read.tree(file=TreeFilename);
}
#tree$edge
#tree$Nnode
#tree$tip.label
#tree$edge.length
#tree$node.label

# Reroot
if(RootAcc!=""){
	is_rooted=is.rooted(tree);
	cat("Is rooted? ", is_rooted, "\n", sep="");
	if(!is_rooted){
		tree=root(tree, RootAcc);
		cat("Tree rooted with: ", RootAcc, "\n", sep="");
	}else{
		cat("Tree is already rooted.\n");
	}
}

num_labels=length(tree$tip.label);

# Annotate accessions with descriptions
original_accessions=tree$tip.label;
tree$tip.label=rename(original_accessions, description_map, keep_accession=TRUE);

# Get colors in the order of the tip labels 
if(!is.null(color_map)){
	colors=rename(original_accessions, color_map, keep_accession=FALSE, default="black");
}else{
	colors=rep(1, num_labels);
}

# Get the number of labels
cat("\n");
cat("Number of leaves (labels): ", num_labels, "\n");

# Size the output pdf height based on the size of the tree
if(Height==-1){
	Height=1/10*num_labels;
}
cat("Height Used: ", Height, "\n");

# Open PDF
pdf(paste(OutputRoot, ".pdf", sep=""), height=Height, width=Width);
par(mar=c(0,0,0,0));
par(oma=c(0,0,0,0));

# Plot tree
plot(tree, cex=.5, tip.color=colors);
add.scale.bar(0,-2,lcol="blue", cex=.6);

###############################################################################

dev.off();

cat("Done.\n");

