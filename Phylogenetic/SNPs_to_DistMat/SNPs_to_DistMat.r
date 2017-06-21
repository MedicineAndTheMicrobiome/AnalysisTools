#!/usr/bin/env Rscript

###############################################################################

library('getopt');
library('MASS');

DEF_PEN_TYPE="liber";

params=c(
	"snps_matrix_filename", "s", 1, "character",
	"keep_list", "k", 2, "character",
	"nuc_penalty", "p", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
        "\nUsage:\n\n", script_name,
        "\n",
	"	-s <SNP Matrix Filename>\n",
	"	[-k <keep list of samples>]\n",
	"	[-p <penalty type, default=",DEF_PEN_TYPE,">]\n",
	"\n",
	"This script will read in a the SNP Matrix file and\n",
	"generate a distance matrix.\n",
	"\n",
	"The input file can be generated from: https://msub.csbio.unc.edu/\n",
	"\n",
	"Input Format:\n",
	"	<samp_id1> , <samp_id2> , ...  <samp_idN> \\n\n",
	"	<allele1,1>  , <allele2,1>  , ...  <alleleN,1> \\n\n",
	"	<allele1,2>  , <allele2,2>  , ...  <alleleN,2> \\n\n",
	"	<allele1,3>  , <allele2,3>  , ...  <alleleN,3> \\n\n",
	"	...\n",
	"	<allele1,M>  , <allele2,M>  , ...  <alleleN,M> \\n\n",
	"\n",
	"Where an <allelei,j> can be {A, T, G, C, etc.}\n",
	"The alleles can also include the IUPAC ambiguity codes.\n",
	"\n",
	"The choices for penalty are:\n",
	"        liber: This will lean towards smaller distances\n",
	"		because if there is any overlap between\n",
	"		the nucleotides of an ambiguity code\n",
	"		the distance will be 0.\n",
	"       probab: This is based on probability.\n",
	"		For example, if K={G,T} and Y={C,T}, then\n",
	"		the truth table is:\n",
	"		    G  T\n",
	"		  C 0  0\n",
	"		  T 0  1\n",
	"\n",
	"		so the distance = 1-1/4 = 3/4\n",	
	"        simil: This depends on the amount of overlap\n",
	"		between two ambiguity codes.\n",
	"		For example if K={G,T} and Y={C,T}, the\n",
	"		possible overlap={C,G,T}, with a length of 3\n",
	"		and the actual overlap={T}, so the distance is\n",
	"		1-1/3 = 2/3.\n",
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

if(length(opt$nuc_penalty)){
	PenaltyType=opt$nuc_penalty;
}else{
	PenaltyType=DEF_PEN_TYPE;
}

#----------------------------------------------------------

InputFilename=opt$snps_matrix_filename;

OutputFilenameRoot=paste(gsub(".csv", "", InputFilename), ".", PenaltyType, sep="");

cat("SNP Input File: ", InputFilename, "\n");
cat("Output File: ", OutputFilenameRoot, "\n");
cat("\n");

################################################################################
# Read SNPs File:

snps_matrix=read.delim(InputFilename, header=T, sep=",", check.names=F, row.names=NULL);
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

compute_nuc_distmat=function(type=""){
	codes=c("A","C","G","T","U","M","R","W","S","Y","K","V","H","D","B","N");
	num_codes=length(codes);
	
	map=list();
	map[["A"]]=c("A");
	map[["C"]]=c("C");
	map[["G"]]=c("G");
	map[["T"]]=c("T");
	map[["U"]]=c("T");

	map[["M"]]=c("A","C");
	map[["R"]]=c("A","G");
	map[["W"]]=c("A","T");
	map[["S"]]=c("C","G");
	map[["Y"]]=c("C","T");
	map[["K"]]=c("G","T");

	map[["V"]]=c("A","C","G");
	map[["H"]]=c("A","C","T");
	map[["D"]]=c("A","G","T");
	map[["B"]]=c("C","G","T");

	map[["N"]]=c("G","A","T","C");

	nuc_distmat=matrix(0, nrow=num_codes, ncol=num_codes, dimnames=list(codes, codes));

	if(type=="simil"){
		for(nuc1 in codes){
			for(nuc2 in codes){
				possibilities=unique(c(map[[nuc1]], map[[nuc2]]));
				intersecting=intersect(map[[nuc1]], map[[nuc2]]);
				nuc_distmat[nuc1, nuc2]=1-length(intersecting)/length(possibilities);	
			}
		}
	}else if(type=="liber"){
		for(nuc1 in codes){
			for(nuc2 in codes){
				intersecting=intersect(map[[nuc1]], map[[nuc2]]);
				nuc_distmat[nuc1, nuc2]=(length(intersecting)==0);
			}
		}
	}else if(type=="probab"){
		for(nuc1 in codes){
			for(nuc2 in codes){
				intersecting=intersect(map[[nuc1]], map[[nuc2]]);
				len1=length(map[[nuc1]]);
				len2=length(map[[nuc2]]);
				nuc_distmat[nuc1, nuc2]=1-length(intersecting)/(len1*len2);
			}
		}
	}else{
		cat("Error: Unrecognized penalty type: ", type, "\n");
		quit(status=-1);
	}

	return(nuc_distmat);

}

penalty_matrix=compute_nuc_distmat(PenaltyType);
cat("Penalty Type: ", PenaltyType, "\n");
print(penalty_matrix);

nuc_dist=function(a, b){
	len=length(a);
	
	tot=0;
	for(i in 1:len){
		tot=tot+penalty_matrix[a[i], b[i]];
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

distmat[distmat==0]=1e-323;
diag(distmat)=0;
dist=as.dist(distmat);

################################################################################

pdf(paste(OutputFilenameRoot, ".snp.dist.pdf", sep=""), height=8.5, width=11);

hcl=hclust(dist, method="ward.D");
plot(hcl, xlab="", main="Ward's Minimum Variance");

hcl=hclust(dist, method="average");
plot(hcl, xlab="", main="Average (i.e. UPGMA)");

print(dist);
isomds=isoMDS(dist, k=2);
print(isomds);
plot(isomds$points[,1], isomds$points[,2], type="n");
text(isomds$points[,1], isomds$points[,2], labels=row.names(distmat));

dev.off();

################################################################################

# Output distance matrix
asFull=as.matrix(dist);
fh=file(paste(OutputFilenameRoot, ".", PenaltyType, ".distmat", sep=""), "w");

sample_names=row.names(distmat);
num_samples=length(sample_names);
for(i in 1:num_samples){
        cat(file=fh, " ", sample_names[i], sep="");
}
cat(file=fh,"\n");
for(i in 1:num_samples){
        cat(file=fh, sample_names[i], asFull[i,], sep=" ");
        cat(file=fh, "\n");
}
close(fh);

################################################################################

cat("Done.\n");
