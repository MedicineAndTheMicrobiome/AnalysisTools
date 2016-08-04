#!/usr/bin/env Rscript

###############################################################################

library('getopt');

UbiquityCutoff=.8;
AbundanceCutoff=1e-4;

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 2, "character",
	"abundance", "a", 2, "numeric",
	"ubiquity", "u", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"	-i <input summary_table.xls file>\n",
	"	[-u <ubiquity cutoff, default >= ", UbiquityCutoff, ">\n",
	"	[-a <abundance cutoff, default >= ", AbundanceCutoff, ">\n",
	"	[-o <output cdf matrix name>]\n",
	"\n",
	"Computes the correlation among the top taxa across all donors.\n",
	"\n",
	"\n", sep="");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputRoot=gsub(".summary_table.xls", "", InputFileName);

if(length(opt$output_file)){
	OutputRoot=opt$output_file;
}

if(length(opt$abundance)){
	AbundanceCutoff=opt$abundance;
}

if(length(opt$ubiquity)){
	UbiquityCutoff=opt$ubiquity;
}

cat("Input File: ", InputFileName, "\n");
cat("Output File Root: ", OutputRoot, "\n");
cat("\n");
cat("Ubiquity Cutoff: ", UbiquityCutoff, "\n");
cat("Abundance Cutoff: ", AbundanceCutoff, "\n");

################################################################################

load_summary_table=function(filename){
	st=as.matrix(read.table(filename, header=TRUE, row.names=1));
	return(st[,2:ncol(st)]);
}

#-------------------------------------------------------------------------------

normalize=function(st){
	sums=apply(st, 1, sum);
	n=matrix(0,nrow=nrow(st), ncol=ncol(st));
	rownames(n)=rownames(st);
	colnames(n)=colnames(st);
	for(i in 1:nrow(st)){
		n[i,]=st[i,]/sums[i];
	}
	return(n);
}

#-------------------------------------------------------------------------------

sig=function(r, n){
	t=r*sqrt((n-2)/(1-(r^2)));
	pval=1-pt(t, n-2);
	return(pval*2); 	# Two tailed adjustment
}

#-------------------------------------------------------------------------------

comp_pval=function(cor_mat, num_samples){
	ncol=ncol(cor_mat);
	pval_mat=matrix(0,ncol=ncol,nrow=ncol);
	rownames(pval_mat)=rownames(cor_mat);
	colnames(pval_mat)=colnames(cor_mat);
	for(i in 1:ncol){
		for(j in 1:ncol){
			pval_mat[i,j]=sig(cor_mat[i,j], num_samples);
		}
	}
	return(pval_mat);
}

#-------------------------------------------------------------------------------

write_matrix=function(mat, fh){
	nrow=nrow(mat);
        sample_names=rownames(mat);
        cat(file=fh, ",", paste(sample_names, collapse=",", sep=","), "\n", sep="");
        for(i in 1:nrow){
                cat(file=fh, sample_names[i], ",", paste(mat[i,], collapse=","), "\n", sep="");
        }
}

################################################################################

# Load summary table
st=load_summary_table(InputFileName);
num_taxa=ncol(st);
num_samples=nrow(st);
taxa_names=colnames(st);
#print(st);

# Normalize counts
nst=normalize(st);
#print(nst);

# Identify taxas exceeding cutoffs
taxas_exc_cutoff=numeric();
num_taxas_exc_cutoff=0;
for(t in 1:num_taxa){
	num_ubiq=sum(nst[,t]>=AbundanceCutoff);
	perc_ubiq=num_ubiq/num_samples;
	if(perc_ubiq>=UbiquityCutoff){
		num_taxas_exc_cutoff=num_taxas_exc_cutoff+1;	
		taxas_exc_cutoff[num_taxas_exc_cutoff]=t;
	}
}
#print(taxas_exc_cutoff);

# Warn if there will be no output
if(length(taxas_exc_cutoff)==0){
	cat("No taxas identified that exceed both cutoffs.\n");
	quit(status=1);
}

# Report which taxas we're going to calculate correlations for
cat("Taxas exceeding cutoff:\n");
for(i in taxas_exc_cutoff){
	cat("\t", taxa_names[i], "\n");
}

# Compute Correlation Coefficients
cor_mat=cor(nst[,taxas_exc_cutoff]);

# Compute P-values
pval_mat=comp_pval(cor_mat, num_samples);

################################################################################

# Output Correlations
fh=file(paste(OutputRoot, UbiquityCutoff, AbundanceCutoff, ".cor.csv", sep=""), "w");
write_matrix(cor_mat, fh);
close(fh);

# Output p-values
fh=file(paste(OutputRoot, UbiquityCutoff, AbundanceCutoff, ".pval.csv", sep=""), "w");
write_matrix(pval_mat, fh);
close(fh);

################################################################################

cat("Done.\n")
q(status=0)
