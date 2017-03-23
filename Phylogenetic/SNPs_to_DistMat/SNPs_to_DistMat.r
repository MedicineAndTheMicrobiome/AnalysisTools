#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"snps_matrix_filename", "s", 1, "character",
	"keep_list", "k", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
        "\nUsage:\n\n", script_name,
        "\n",
	"	-s <SNP Matrix Filename>\n",
	"	[-k <keep list of samples>]\n",
	"\n",
	"This script will read in a the SNP Matrix file and\n",
	"generate a distance matrix.\n",
	"\n",
	"The input file can be generated from: https://msub.csbio.unc.edu/\n",
	"\n",
	"Input Format:\n",
	"	<loci_title> , <samp_id1> , <samp_id2> , ...  <samp_idN> \\n\n",
	"	<loci_1>     , <allele1,1>  , <allele2,1>  , ...  <alleleN,1> \\n\n",
	"	<loci_2>     , <allele1,2>  , <allele2,2>  , ...  <alleleN,2> \\n\n",
	"	<loci_3>     , <allele1,3>  , <allele2,3>  , ...  <alleleN,3> \\n\n",
	"	...\n",
	"	<loci_M>     , <allele1,M>  , <allele2,M>  , ...  <alleleN,M> \\n\n",
	"\n",
	"Where an <allelei,j> can be {A, T, G, C}\n",
        "\n", sep="");

if(!length(opt$snps_matrix_filename)){
        cat(usage);
        q(status=-1);
}

if(length(opt$keep_list)){
	KeepList=opt$keep_list;
}else{
	KeepList="";
}

#----------------------------------------------------------

InputFilename=opt$snps_matrix_filename;

OutputFilenameRoot=paste(gsub(".csv", "", InputFilename), sep="");

cat("SNP Input File: ", InputFilename, "\n");
cat("Output File: ", OutputFilenameRoot, "\n");
cat("\n");

################################################################################
# Read SNPs File:

snps_matrix=read.delim(InputFilename, header=T, sep=",", check.names=F, row.names=1);
#print(snps_matrix);
num_loci=nrow(snps_matrix);
num_samp=ncol(snps_matrix);

cat("Num Loci: ", num_loci, "\n");
cat("Num Samples: ", num_samp, "\n");
cat("\n");
prof_matrix=t(snps_matrix);

################################################################################
# Read Keep List:

samp_ids=rownames(prof_matrix);
cat("Samples in SNPs Matrix:\n");
print(sort(samp_ids));
cat("\n");

if(KeepList!=""){
	sample_keep_list=scan(KeepList, what=character());
	cat("Keep List samples:\n");
	print(sample_keep_list);
	cat("\n");
	overlapping_samples=intersect(sample_keep_list, samp_ids);

	cat("Overlapping samples with keep list:\n");
	print(sort(overlapping_samples));
	cat("\n");
	prof_matrix=prof_matrix[overlapping_samples,];

	cat("Samples in Keep List missing from SNP Matrix:\n");
	print(sort(setdiff(sample_keep_list, overlapping_samples)));
	cat("\n");
}

samp_ids=sort(rownames(prof_matrix));
prof_matrix=prof_matrix[samp_ids,];

num_samp=nrow(prof_matrix);

################################################################################

print(prof_matrix);


nuc_dist=function(a, b){
	len=length(a);
	
	tot=0;
	for(i in 1:len){
		if(a[i]!=b[i]){
			tot=tot+1;
		}
	}

	dist=tot/len;
}

distmat=matrix(0, nrow=num_samp, ncol=num_samp, dimnames=list(samp_ids, samp_ids));

for(i in 1:num_samp){
	for(j in 1:i){
		distmat[i,j]=nuc_dist(prof_matrix[i,], prof_matrix[j,]);
		distmat[j,i]=distmat[i,j];
	}
}

dist=as.dist(distmat);

################################################################################

pdf(paste(OutputFilenameRoot, ".snp.dist.pdf", sep=""), height=8.5, width=11);

hcl=hclust(dist, method="ward.D");
plot(hcl, xlab="", main=OutputFilenameRoot);

dev.off();

################################################################################

cat("Done.\n");
