#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file_i", "i", 1, "character",
	"input_file_j", "j", 1, "character",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
		"\nUsage:\n\t", script_name, "\n",
		"\t\t-i <Input summary_table.xls 1 FileName, ie. taxon-based>\n",
		"\t\t-j <Input summary_table.xls 2 FileName, ie. otu-based>\n",
		"\t\t[-o <Output summary_table.xls FileName>]\n",
		"\n",
		"Compares two summary_table.xls files with the same sample.\n",
		"\n\n"
);

if(!length(opt$input_file_i) || !length(opt$input_file_j)){
	cat(usage);
	q(status=-1);
}

###############################################################################
# Main program loop

InputFileName1=opt$input_file_i;
InputFileName2=opt$input_file_j;

OutputFilename="Comparison";
if(length(opt$output_file)){
	OutputFilename=opt$output_file;	
}

###############################################################################
###############################################################################

load_counts=function(summary_file_fn){
	inmat<-as.matrix(read.table(summary_file_fn, sep="\t", header=TRUE, check.names=FALSE, row.names=1));
	count_mat=inmat[,2:ncol(inmat)];

	# Sort rows by names
	names=rownames(inmat);
	names_sorted=sort(names, index.return=TRUE);
	count_mat=count_mat[names_sorted$ix,];

	return(count_mat);
}

num_taxa=function(count_mat){
	num_samples=nrow(count_mat);
	num_taxa=rep(0, num_samples);
	for(i in 1:num_samples){
		num_taxa[i]=sum(count_mat[i,]>0);
	}
	names(num_taxa)=rownames(count_mat);
	return(num_taxa);
}

confirm_matching=function(count_mat1, count_mat2){
	names1=rownames(count_mat1);
	names2=rownames(count_mat2);

	if(length(names1)!=length(names2)){
		return(FALSE);
	}

	for(i in 1:length(names1)){
		if(names1[i]!=names2[i]){
			return(FALSE);
		}
	}
	return(TRUE);
}

###############################################################################

# Load summary tables
in1mat=load_counts(InputFileName1);
in2mat=load_counts(InputFileName2);

# Confirma both samples line up
if(!confirm_matching(in1mat, in2mat)){
	cat("Error:  Samples do not match between summary_files.xls.\n");
	q(status=-1);
}
sample_names=rownames(in1mat);

# Counts num taxa per sample
num_taxa1=num_taxa(in1mat);
num_taxa2=num_taxa(in2mat);

# Compute ratios
ratios=num_taxa2/num_taxa1;
names(ratios)=sample_names;
#print(ratios);

sorted_ratios=sort(ratios, decreasing=TRUE, index.return=TRUE);
print(sorted_ratios);

# Plot num_taxas
pdf(paste(OutputFilename, ".pdf", sep=""), height=11, width=8.5);
barheights=rbind(num_taxa1, num_taxa2);
barheights=barheights[,sorted_ratios$ix];

adj_barheights=barheights;
adj_barheights[2,]=barheights[2,]-barheights[1,];

par(oma=c(10,.5,.5,.5));
xpos=barplot(adj_barheights,names.arg=names(barheights), las=2, 
	ylim=c(0, 1.2*max(barheights)),
	ylab="Taxonomic/OTU counts"
	);
text(xpos,barheights[1,]+10, barheights[1,], pos=1, offset=.2, cex=.7, col="green");
text(xpos,barheights[2,], barheights[2,], pos=1, offset=.2, cex=.7, col="blue");
text(xpos,barheights[2,], labels=sprintf("%3.2f", sorted_ratios$x), pos=3, offset=.5, cex=.9, col="red");
legend(xpos[length(xpos)]*.8, max(barheights)*.8, 
	legend=c("Ratio/Est. Deg", "OTUs", "Taxa"),
	fill=c("red", "blue", "green"));


#print(num_taxa1);
#cat("\n");
#print(num_taxa2);

q(status=1);
