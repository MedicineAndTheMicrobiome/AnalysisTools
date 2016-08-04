#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file_gene_list", "e", 1, "character",
	"input_file_taxa_list", "t", 1, "character",
	"output_filename", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n   ", script_name, "\n",
	"	-e <rarefaction curve file genes list>\n",
	"	-t <rarefaction curve file taxa list>\n",
	"	[-o <output file name>\n",
	"\n",
	"This script will read in two files, each of which is a list of rarefaction curve files\n",
	"with medians, means, and 95% confidence intervals. \n",
	"\n",
	"The rows should go like this:\n",
	"\n",
	"	<Label>\\t<medians starting with 1 sample>\n",
	"	<Label>\\t<means starting with 1 sample>\n",
	"	<Label>\\t<95% CI lowerbound starting with 1 sample>\n",
	"	<Label>\\t<95% CI upperbound starting with 1 sample>\n",
	"\n",
	"\n");

if(!length(opt$input_file_gene_list) || !length(opt$input_file_taxa_list)){
	cat(usage);
	q(status=-1);
}

InputFileNameListGenes=opt$input_file_gene_list;
InputFileNameListTaxa=opt$input_file_taxa_list;
OutputFilename=opt$output_filename;

if(length(OutputFilename)==0){
	gene_file_name=tail(strsplit(InputFileNameListGenes, "/")[[1]],1);
	taxa_file_name=tail(strsplit(InputFileNameListTaxa, "/")[[1]],1);
	OutputFilename=paste(gene_file_name, "_vs_", taxa_file_name, sep="");
}
cat("Output Filename Root: ", OutputFilename, "\n");

###############################################################################
# Load counts from file

load_file_list=function(filename){
	cat("Loading list: ", filename, "...\n");
	dat=as.matrix(read.table(filename));
	cat("Done.\n");
	return(as.vector(dat[,1]));
}

load_rarefaction_curves=function(filename){
	cat("Loading rarefaction file: ", filename, "...\n");
	dat=as.matrix(read.table(filename, sep="\t", check.names=FALSE, row.names=1, fill=TRUE));
	rownames(dat)=c();
	return(dat);
}

load_curves_into_list=function(filename_list){
	rare_list=list();
	list_len=length(filename_list);
	cat("Number of files to load: ", list_len, "\n");
	#print(filename_list);
	for(i in 1:list_len){
		rare_list[[filename_list[i]]]=load_rarefaction_curves(filename_list[i]);
	}
	return(rare_list);
}

###############################################################################

genes_list=load_file_list(InputFileNameListGenes);
taxa_list=load_file_list(InputFileNameListTaxa);

genes_list_length=length(genes_list);
taxa_list_length=length(taxa_list);
if(genes_list_length!=taxa_list_length){
	cat("Error: The number of samples in your genes list does not equal the number in your taxa list.\n");
	quit(status=-1);
}else{
	cat("Ok, number of samples between genes and taxa match.\n");
	num_sample_pairs=genes_list_length;
}

genes_rare_list=load_curves_into_list(genes_list)
taxa_rare_list=load_curves_into_list(taxa_list)

###############################################################################

num_samples_in_genes_rarefaction=numeric();
num_samples_in_taxa_rarefaction=numeric();
max_discovered=0;

for(i in 1:num_sample_pairs){
	upperbound_genes_rarefaction_curve=(genes_rare_list[[i]][4,]); # 4 is 4th row (upperbound) in matrix 
	upperbound_taxa_rarefaction_curve=(taxa_rare_list[[i]][4,]); # 4 is 4th row (upperbound) in matrix 

	num_samples_in_genes_rarefaction[i]=length(upperbound_genes_rarefaction_curve);
	num_samples_in_taxa_rarefaction[i]=length(upperbound_taxa_rarefaction_curve);

	if(max_discovered<max(upperbound_genes_rarefaction_curve)){
		max_discovered=max(upperbound_genes_rarefaction_curve);
	}
	if(max_discovered<max(upperbound_taxa_rarefaction_curve)){
		max_discovered=max(upperbound_taxa_rarefaction_curve);
	}
}

max_sampled=max(c(num_samples_in_genes_rarefaction, num_samples_in_taxa_rarefaction));
cat("\n");
cat("Max Discovered (among all): ", max_discovered, "\n");
cat("Num samples in genes: ", num_samples_in_genes_rarefaction, "\n");
cat("Num samples in taxa: ", num_samples_in_taxa_rarefaction, "\n");
cat("Max Sampled (among all): ", max_sampled, "\n");

cat("\n");

# Clean names
clean_name=character();
for(i in 1:genes_list_length){
	sample_name=genes_list[i];
	clean_name[i]=gsub("_", " ", tail(strsplit(sample_name,"/")[[1]],1));
	cat("Working on ", sample_name, "\n");
	cat("\tCleaned name: ", clean_name[i], "\n");
}

###############################################################################

library(plotrix);
pdf(paste(OutputFilename, ".pdf", sep=""), height=11, width=8.5);

###############################################################################

cat("\n");
cat("Working on combined curves...\n");
# Plot combined rarefaction curves
plot(0,0, type="n", xlim=c(0,max_sampled*1.1), ylim=c(0,log10(max_discovered)), ylab="Log(Discovered)", xlab="Number of Samples", main="Samples Taken vs. Discovered");
for(i in 1:num_sample_pairs){
	median_genes_curve=as.vector(genes_rare_list[[i]][1,]); # median
	lb_genes_curve=as.vector(genes_rare_list[[i]][3,]); # lowerbound
	ub_genes_curve=as.vector(genes_rare_list[[i]][4,]); # upperbound

	median_taxa_curve=as.vector(taxa_rare_list[[i]][1,]); # median
	lb_taxa_curve=as.vector(taxa_rare_list[[i]][3,]); # lowerbound
	ub_taxa_curve=as.vector(taxa_rare_list[[i]][4,]); # upperbound

	median_genes_curve=log10(median_genes_curve);
	lb_genes_curve=log10(lb_genes_curve);
	ub_genes_curve=log10(ub_genes_curve);
	median_taxa_curve=log10(median_taxa_curve);
	lb_taxa_curve=log10(lb_taxa_curve);
	ub_taxa_curve=log10(ub_taxa_curve);

	ebarwidth=.001;
	plotCI(1:num_samples_in_taxa_rarefaction[i], median_taxa_curve, ui=ub_taxa_curve, li=lb_taxa_curve, add=TRUE, cex=.3, col=i, sfrac=ebarwidth);
	plotCI(1:num_samples_in_genes_rarefaction[i], median_genes_curve, ui=ub_genes_curve, li=lb_genes_curve, add=TRUE, cex=.3, col=i, sfrac=ebarwidth);

	taxalastpoint=tail(median_taxa_curve,1);
	geneslastpoint=tail(median_genes_curve,1);
	
	text(num_samples_in_taxa_rarefaction[i], taxalastpoint, 10^taxalastpoint, cex=.6, col=i, pos=4);
	text(num_samples_in_genes_rarefaction[i], geneslastpoint, 10^geneslastpoint, cex=.6, col=i, pos=4);

	points(num_samples_in_taxa_rarefaction[i]+1, taxalastpoint, pch=15, col=i, cex=1);
	points(num_samples_in_genes_rarefaction[i]+1, geneslastpoint, pch=16, col=i, cex=1);

}
legend(max_sampled*.70, log10(max_discovered)*.25, 
	legend=c("Genes", "Taxa", clean_name), 
	col=c("black", "black", 1:num_sample_pairs), 
	pch=c(1, 0, rep(15, num_sample_pairs)),
	cex=1
	);

####################################################################
cat("\n");
cat("Working on individual...\n");
for(i in 1:num_sample_pairs){
	cat("Working on: ", clean_name[i], "\n");
	median_genes_curve=as.vector(genes_rare_list[[i]][1,]); # median
	lb_genes_curve=as.vector(genes_rare_list[[i]][3,]); # lowerbound
	ub_genes_curve=as.vector(genes_rare_list[[i]][4,]); # upperbound

	median_taxa_curve=as.vector(taxa_rare_list[[i]][1,]); # median
	lb_taxa_curve=as.vector(taxa_rare_list[[i]][3,]); # lowerbound
	ub_taxa_curve=as.vector(taxa_rare_list[[i]][4,]); # upperbound

	median_genes_curve=log10(median_genes_curve);
	lb_genes_curve=log10(lb_genes_curve);
	ub_genes_curve=log10(ub_genes_curve);
	median_taxa_curve=log10(median_taxa_curve);
	lb_taxa_curve=log10(lb_taxa_curve);
	ub_taxa_curve=log10(ub_taxa_curve);

	max_local_sampled=max(length(median_genes_curve), length(median_taxa_curve));
	max_local_discovered=max(max(ub_genes_curve), max(ub_taxa_curve));	

	plot(0,0, type="n", xlim=c(0,max_local_sampled*1.1), ylim=c(0,max_local_discovered), 
		ylab="Log(Discovered)", xlab="Num Samples", 
		main=paste(clean_name[i],": Samples Taken vs. Discovered"));

	ebarwidth=.002
	plotCI(1:length(median_taxa_curve), median_taxa_curve, ui=ub_taxa_curve, li=lb_taxa_curve, add=TRUE, cex=.3, col=i, sfrac=ebarwidth);
	plotCI(1:length(median_genes_curve), median_genes_curve, ui=ub_genes_curve, li=lb_genes_curve, add=TRUE, cex=.3, col=i, sfrac=ebarwidth);

	taxalastpoint=tail(median_taxa_curve,1);
	geneslastpoint=tail(median_genes_curve,1);
	
	max_local_genes_sampled=length(median_genes_curve);
	max_local_taxa_sampled=length(median_taxa_curve);

	text(max_local_taxa_sampled+1, taxalastpoint, 10^taxalastpoint, cex=.6, col=i, pos=4);
	text(max_local_genes_sampled+1, geneslastpoint, 10^geneslastpoint, cex=.6, col=i, pos=4);

	points(max_local_taxa_sampled+1, taxalastpoint, pch=15, col=i, cex=1);
	points(max_local_genes_sampled+1, geneslastpoint, pch=16, col=i, cex=1);

	legend(max_local_sampled*.70, log10(max_local_discovered)*.45, 
		legend=c("Genes", "Taxa"), 
		col=c(i,i), 
		pch=c(16, 15),
		cex=1
		);
}

####################################################################
cat("\n");
cat("Working on proportions...\n");
logProp=list();
maxProp=0;
minProp=0;
for(i in 1:num_sample_pairs){
	median_genes_curve=as.vector(genes_rare_list[[i]][1,]); # median
	median_taxa_curve=as.vector(taxa_rare_list[[i]][1,]); # median
	min_length=min(length(median_genes_curve), length(median_taxa_curve));
	logProp[[i]]=log10(median_genes_curve[1:min_length]/median_taxa_curve[1:min_length]);
	maxProp=max(maxProp, max(logProp[[i]]));
	minProp=min(minProp, min(logProp[[i]]));

}

plot(0,0, type="n", xlim=c(0,max_local_sampled), ylim=c(minProp,maxProp), ylab="Log(Genes/Taxa)", xlab="Num Samples", main=paste("Samples Taken vs. Log(Gene/Taxa)"));
abline(h=0, col="grey");
for(i in 1:num_sample_pairs){
	points(1:length(logProp[[i]]), logProp[[i]], cex=.6, col=i);
	lastpoint=tail(logProp[[i]],1);
	
	if(lastpoint>0){
		ratio_string=sprintf("%i : 1", round(10^lastpoint));
	}else{
		ratio_string=sprintf("1 : %i", round(10^(-lastpoint)));
	}

	text(length(logProp[[i]]), lastpoint, ratio_string,  cex=.6, col=i, pos=4);
}

legend(max_local_sampled*.75, maxProp*.25, legend=clean_name, fill=1:num_sample_pairs, cex=.8);











cat("Done.\n");
q(status=0);
