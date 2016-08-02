#!/usr/bin/env Rscript

###############################################################################

library('getopt');
library('ape');

PAPER_WIDTH=8.5;
PAPER_HEIGHT=-1;

params=c(
	"tree_fn", "t", 1, "character",
	"num_clusters", "k", 1, "numeric",
	"color_fn", "c", 2, "character",
	"description_fn", "d", 2, "character",
	"output_fn_root", "o", 2, "character",
	"width", "w", 2, "numeric",
	"height", "h", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
        "\nUsage:\n\n", script_name,
        "\n",
	"	-t <tree file name>\n",
	"	-k <number of clusters to split>\n",
	"	[-o <output filename root>]\n",
	"	[-d <description map>]\n",
	"	[-w <paper width, default=", PAPER_WIDTH, ">]\n",
	"	[-h <paper height, default is dynamic>]\n",
	"\n",
	"This script will read in a tree in Newick format, and\n",
	"generates a cophenetic distance matrix.  The distance\n",
	"matrix is then hierarchically clustered and cut based\n",
	"on the specified -k paramter.  Both the phylogenetic\n",
	"and Ward's based dendrogram is plotted with the nodes\n",
	"colored based on the hierarchical clusters.\n",
	"\n",
	"A file with a list of members for each cluster is generated.\n",
        "\n", sep="");

if(!length(opt$tree_fn) || !length(opt$num_clusters)){
        cat(usage);
        q(status=-1);
}

#----------------------------------------------------------

DescriptionFilename="";

TreeFilename=opt$tree_fn;
NumClusters=opt$num_clusters;

if(length(opt$output_fn_root)){
	OutputRoot=opt$output_fn_root;
}else{
	OutputRoot=TreeFilename;
}
OutputRoot=paste(OutputRoot, sprintf(".%02i", NumClusters), sep="");

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

# Load tree from file
tree=read.tree(file=TreeFilename);
#tree$edge
#tree$Nnode
#tree$tip.label
#tree$edge.length
#tree$node.label

num_labels=length(tree$tip.label);

# Annotate accessions with descriptions
original_accessions=tree$tip.label;
tree$tip.label=rename(original_accessions, description_map, keep_accession=TRUE);

# Get the number of labels
cat("\n");
cat("Number of leaves (labels): ", num_labels, "\n");

# Size the output pdf height based on the size of the tree
if(Height==-1){
	Height=1/10*num_labels;
}
cat("Height Used: ", Height, "\n");


###############################################################################

set_lab=function(n, colormap) {
	if(is.leaf(n)) {
		leaf_attr=attributes(n);
		leaf_name=leaf_attr$label;
		attr(n, "nodePar")$lab.col=colormap[[leaf_name]];
		attr(n, "nodePar")$lab.cex=.5;
		attr(n, "nodePar")$cex=0;
	}
	return(n);	
}

###############################################################################

# Open PDF
pdf(paste(OutputRoot, ".pdf", sep=""), height=Height, width=Width);
par(mar=c(0,0,0,5));

# Get cophenetic distance among all members in the tree
cophenetic.dist=cophenetic(tree);
hcl=hclust(as.dist(cophenetic.dist), method="ward");
den=as.dendrogram(hcl);

# Assign colors to groups
if(NumClusters<=10){
	colorlist=c("black","red","green","orange", "blue", "pink", "brown", "grey", "cyan", "salmon");
}else{
	colorlist=rainbow(NumClusters, end=.85);
}
palette(colorlist);

# Cut tree based on user specified number of clusters
parts=cutree(hcl, k=NumClusters);

# Plot hierarchical clsuter 
den=dendrapply(den, set_lab, parts);
plot(den, horiz=T);

# Output members into file
num_clusters=max(parts);
for(cl in 1:num_clusters){
        outputfilename=paste(OutputRoot, ".", sprintf("%02i",cl), sep="");
        fc=file(outputfilename, "w");
        members=names(parts[parts==cl]);
        write(members,file=fc);
        close(fc);
}

# Plot phylogenetic tree
plot(tree, cex=.5, tip.color=parts[tree$tip.label]);

###############################################################################

dev.off();

cat("Done.\n");

