#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"pathway_map", "p", 1, "character",
	"output_file", "o", 1, "character",
	"pathway_names", "n", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input imputed summary_table.tsv file>\n",
	"	-o <output root>\n",
	"	-p <pathway map>\n",
	"	[-n <pathway names>]\n",
	"\n",
	"This script will read in the abundances for each of the ECs\n",
	"from the summary table.  Based on the pathway map, a log(likelihood)\n",
	"for each pathway will be calculated by summing together all the log10(probabilities)\n",
	"(abundances) for each of the steps. \n",
	"\n",
	"\n", sep="");

if(
	!length(opt$input_file) || 
	!length(opt$output_file)){
	cat(usage);
	q(status=-1);
}

InputFilename=opt$input_file;
OutputFilename=opt$output_file;
PathwayMap=opt$pathway_map;

PathwayNames="";
if(length(opt$pathway_names)){
	PathwayNames=opt$pathway_names;
}

###############################################################################

load_factors=function(fname){
        factors=data.frame(read.table(fname,  header=TRUE, check.names=FALSE, row.names=1, 
		comment.char="", quote="", sep="\t", stringsAsFactors=TRUE));
        dimen=dim(factors);
        cat("Rows Loaded: ", dimen[1], "\n");
        cat("Cols Loaded: ", dimen[2], "\n");
        return(factors);
}

load_summary_file=function(fname){
        cat("Loading Summary Table: ", fname, "\n");
        inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", row.names=1))
        counts_mat=inmat[,2:(ncol(inmat))];
        return(counts_mat);
}

plot_text=function(strings){

	orig.par=par(no.readonly=T);
        par(family="Courier");
        par(oma=rep(.1,4));
        par(mar=rep(0,4));

        num_lines=length(strings);

        top=max(as.integer(num_lines), 52);

        plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                );

        text_size=max(.01, min(.8, .8 - .003*(num_lines-52)));
        #print(text_size);

        for(i in 1:num_lines){
                #cat(strings[i], "\n", sep="");
                strings[i]=gsub("\t", "", strings[i]);
                text(0, top-i, strings[i], pos=4, cex=text_size);
        }

	par(orig.par);
}

load_pathway_map=function(file){
	# File Format:
	file_data=read.delim(file, header=F, check.names=F, 
		stringsAsFactors=F, comment.char="#", row.names=NULL);
	ec_pathway_data=file_data[,c(1,2)];
	colnames(ec_pathway_data)=c("EC", "Pathway");
	pathway_sort_ix=order(ec_pathway_data[,"Pathway"]);
	
	ec_pathway_data=ec_pathway_data[pathway_sort_ix,];

	unique_pathways=sort(unique(ec_pathway_data[,"Pathway"]));
	num_pathways=length(unique_pathways);
	
	cat("Unique Pathways:", num_pathways, "\n");

	map_as_list=list();
	for(i in 1:num_pathways){
		cur_path=unique_pathways[i];
		path_ix=which(cur_path==ec_pathway_data[,"Pathway"]);
		map_as_list[[cur_path]]=sort(unique(ec_pathway_data[path_ix, "EC"]));
	}

	#aprint(map_as_list);
	return(map_as_list);
	
}

load_pathway_names=function(file){
	
	file_data=read.delim(file, header=T, check.names=F, 
		stringsAsFactors=F, comment.char="#", row.names=NULL);
	map_mat=file_data[,c(1,2)];
	num_names=nrow(map_mat);
	
	hash=list();
	for(i in 1:num_names){
		hash[[map_mat[i, 1]]]=map_mat[i,2];
	}

	return(hash);
}


###############################################################################

pathway_name_hash=NULL;
if(PathwayNames!=""){
	pathway_name_hash=load_pathway_names(PathwayNames);
}

map_list=load_pathway_map(PathwayMap);
abund_mat=load_summary_file(InputFilename);

num_categories=ncol(abund_mat);
num_samples=nrow(abund_mat);
sample_ids=rownames(abund_mat);

cat("Num Categories: ", num_categories, "\n", sep="");
cat("Num Samples: ", num_samples, "\n", sep="");


pdf(paste(OutputFilename, ".pwy_loglik.hist.pdf", sep=""), height=11, width=8.5);


log_prob_mat=log10(abund_mat);
avail_categories=colnames(abund_mat);

num_pathways=length(map_list);
pathway_names=names(map_list);

plot_wmesg=function(msg, main, tsize){
	plot(0,0, type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", main=main, cex.main=tsize);
	text(0,0, msg);
}


fullfilled_pathway_lengths=rep(NA, num_pathways);

par(mfrow=c(5,2));

title_size=.9;

saved_logliks=list();
for(i in 1:num_pathways){

	cur_pathway_rxns=map_list[[i]];
	cur_pathway_name=pathway_names[i];
	num_rxns_in_pathway=length(cur_pathway_rxns);

	# Look up full pathway name if available
	if(!is.null(pathway_name_hash)){
		full_name=pathway_name_hash[[cur_pathway_name]];
	}else{
		full_name="";
	}

	# Pathway info
	pathway_info_str=paste("[", i, "/", num_pathways, "] ", cur_pathway_name, 
		ifelse(full_name!="", paste("\n", full_name, sep=""), ""),
		"\n(", num_rxns_in_pathway, " Reactions)", sep="");
	cat("\n", pathway_info_str, "\n", sep="");
	print(cur_pathway_rxns);

	# Check for any missing rxns
	missing_rxns=setdiff(cur_pathway_rxns, avail_categories);
	num_missing=length(missing_rxns);

	# If missing any rxns, skip
	if(num_missing){
		cat("Missing ", num_missing, " reactions.\n", sep=""); 
		cat("Skipping...\n");

		if(num_missing==num_rxns_in_pathway){
			num_missing="ALL";
		}

		plot_wmesg(paste(
			"Missing ", num_missing, " Reaction(s) from Pathway.", sep=""), 
			main=pathway_info_str, tsize=title_size);
		next;
	}

	# Extract log probabilities
	extracted_log_prob=log_prob_mat[,cur_pathway_rxns,drop=F];
	#print(extracted_log_prob);

	# Sum up log probabilities to get likelikhoods
	pathway_loglikelihood=apply(extracted_log_prob, 1, sum);

	# Check make sure no rxns are still zeros (failed previous imputation?)
	if(!is.finite(min(pathway_loglikelihood))){
		cat("Too many zeros... (-Inf Log Likelihoods)\n", sep=""); 
		cat("Skipping...\n");
		plot_wmesg("Zeros/-Inf for reaction(s) in Pathway", main=pathway_info_str, tsize=title_size);
		next;
	}

	# Plot hist 
	hist(pathway_loglikelihood, main=pathway_info_str, xlab="Pathway Log Likelihood", 
		cex.main=title_size,
		breaks=50);
	
	saved_logliks[[cur_pathway_name]]=pathway_loglikelihood;

	# Keep track of which pathways were fullfilled
	fullfilled_pathway_lengths[i]=num_rxns_in_pathway;
	#print(pathway_likelihood);
}

###############################################################################

fullfilled_pathway_lengths=fullfilled_pathway_lengths[!is.na(fullfilled_pathway_lengths)];

par(mfrow=c(1,1));
hist(fullfilled_pathway_lengths, 
	breaks=30,
	main="Distribution of Fullfilled Pathway Lengths", xlab="Num Reactions");

###############################################################################
# Export computed log liks

cat("\n");
cat("Saving log likelihoods to matrix...\n");
num_fullfilled_pwy=length(saved_logliks);
fullfilled_pwy_names=names(saved_logliks);
outmat=matrix(numeric(), nrow=num_samples, ncol=num_fullfilled_pwy);
colnames(outmat)=fullfilled_pwy_names;
rownames(outmat)=sample_ids;

for(fpwy in fullfilled_pwy_names){
	outmat[,fpwy]=saved_logliks[[fpwy]];
}

cat("Writing log likelihood matrix to file...\n");
SampleID=sample_ids;
outmat=cbind(SampleID, outmat);
write.table(outmat, file=paste(OutputFilename, ".pwy_loglik.tsv", sep=""),
	sep="\t", row.names=F, col.names=T, quote=F);


###############################################################################

dev.off();

###############################################################################

cat("Done.\n")
warn=warnings();
if(length(warn)){
	print(warn);
}
q(status=0)
