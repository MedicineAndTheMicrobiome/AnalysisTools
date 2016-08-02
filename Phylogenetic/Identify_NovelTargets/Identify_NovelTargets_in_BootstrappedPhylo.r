#!/usr/bin/env Rscript

###############################################################################

library('getopt');
library('ape');

PAPER_WIDTH=8.5;
PAPER_HEIGHT=-1;

params=c(
	"tree_fn", "t", 1, "character",
	"ref_fn", "r", 1, "character",
	"target_fn", "x", 1, "character",
	"output_fn_root", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
        "\nUsage:\n\n", script_name,
        "\n",
	"	-t <bootstrapped trees file name>\n",
	"	-r <references filename (accessions for characterized taxa, eg. SILVA or GenBank)>\n",
	"	-x <targets filename (accessions for uncharacterized taxa, ie. OTUs)>\n",
	"	[-o <output filename root>]\n",
	"\n",
	"This script will read in a set of bootstrapped trees, then based on which IDs are\n",
	"targets (OTUs) determine if they are different enough to a reference sequence\n",
	"to be considered novel.\n",
	"\n",
	"Novelty is based on determining whether >95% of the bootstrapped cophenetic distances\n",
	"between the closest reference and the targeted OTU are greater than 0.\n",
	"\n",
	"Uncorrected p-values are reported for each effect size (branch length).\n",
	"\n",
        "\n", sep="");

if(!length(opt$tree_fn) || !length(opt$ref_fn) || !length(opt$target_fn)){
        cat(usage);
        q(status=-1);
}

#------------------------------------------------------------------------------

if(length(opt$tree_fn)){
	TreeFilename=opt$tree_fn;
}

if(length(opt$output_fn_root)){
	OutputRoot=opt$output_fn_root;
}else{
	OutputRoot=TreeFilename;
}

if(length(opt$ref_fn)){
	ReferenceFname=opt$ref_fn;
}

if(length(opt$target_fn)){
	TargetFname=opt$target_fn;
}
###############################################################################

# Load tree from file
cat("Loading tree: ", TreeFilename, "\n");
tree=read.tree(file=TreeFilename);

tree_or_list=attributes(tree);
cat("File Type: ", tree_or_list$class, "\n");
if(tree_or_list$class=="phylo"){
	cat("\n\nError:  You only have one tree.  You need bootstraps of trees, not the ML tree.\n\n");
	quit(status=-1);
}	

num_bootstraps=length(tree);
cat("Num trees/bootstraps loaded: ", num_bootstraps, "\n", sep="");

# Load reference accessions
ref_ids=sort(scan(ReferenceFname, what=character(), sep="\n"));
num_references=length(ref_ids);
cat("Num references loaded: ", num_references, "\n", sep="");

# Load target accessions
target_ids=sort(scan(TargetFname, what=character(), sep="\n"));
num_targets=length(target_ids);
cat("Num targets loaded: ", num_targets, "\n", sep="");

################################################################################

# Generate cophenetic distance matrices and subset targets x references
cat("Computing cophenetic distances...\n");
cophenetic_dist=array(0, dim=c(num_targets, num_references, num_bootstraps));
colnames(cophenetic_dist)=ref_ids;
rownames(cophenetic_dist)=target_ids;

for(i in 1:num_bootstraps){
    cophenetic_dist[,,i]=cophenetic.phylo(tree[[i]])[target_ids, ref_ids];
}
cat("done.\n\n");

# Note: the rows and columns are ordered alphebetically
#print(cophenetic_dist)

################################################################################

# Specify range and interval of effects to calculate
effects_values=seq(0, .3, .05);
num_effects=length(effects_values);

# Compute indices for upper and lower confidence intervals
alpha=0.05;
ub_ix=ceiling(num_bootstraps*(1-alpha/2));
lb_ix=ceiling(num_bootstraps*(alpha/2));

# Open output text file
fh=file(paste(OutputRoot, ".novelty.txt", sep=""), "w");

# Output header
effects_str=sprintf("Effects_%0.2f", effects_values);
cat(file=fh, paste("Target", "Closest_Reference", "Median", "95%_LowerBound", "95%_UpperBound", paste(effects_str, collapse="\t"), sep="\t"), "\n");

# For each target, identify closest reference by median cophenetic distance
for(i in 1:num_targets){

	cat("\nComputing stats for: ", target_ids[i], "\n");

	# Compute median distance to every reference
	med_coph_dists=numeric(num_references);
	for(j in 1:num_references){
		med_coph_dists[j]=median(cophenetic_dist[i,j,]);
	}

	# Identify closest reference
	closest_ref_ix=which(med_coph_dists==min(med_coph_dists))[1];	
	closest_ref=ref_ids[closest_ref_ix];
	closest_dist=med_coph_dists[closest_ref_ix];

	# Extract out bootstrap distances for closest reference
	cat("\tClosest reference: ", closest_ref, " (", closest_dist, ") \n", sep="");
	sorted_coph=sort(cophenetic_dist[i, closest_ref_ix, ]);

	# Compute 95% confidence intervals
	ub=sorted_coph[ub_ix];
	lb=sorted_coph[lb_ix];

	# Compute p-value for each effect
	effects_pvalues=numeric(num_effects);
	for(e in 1:num_effects){
		effects_pvalues[e]=sum(sorted_coph<effects_values[e])/num_bootstraps;
		#print(effects_pvalues[e]);
	}

	# Report p-values to each effect size
	cat(file=fh, paste(
		target_ids[i],
		closest_ref,
		closest_dist,
		lb,
		ub,
		paste(effects_pvalues, collapse="\t"),
		sep="\t"
		)
	);
	cat(file=fh,"\n");
		
}

cat("\n\n");

################################################################################

close(fh);
#dev.off()
cat("Done.\n");
