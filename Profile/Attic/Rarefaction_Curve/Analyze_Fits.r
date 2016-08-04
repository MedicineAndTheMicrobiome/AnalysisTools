#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"	-i <input fit_and_estimates.csv >\n",
	"\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFileRoot=gsub(".csv", "", InputFileName);

cat("Output file name root:", OutputFileRoot, "\n", sep="");

###############################################################################
# Load counts from file

confid=function(x){
	x=sort(x);
	bidx=.05/2*length(x);
	lb=x[1+bidx];
	ub=x[length(x)-bidx];
	ret=list();
	ret$lb=lb;
	ret$ub=ub;
	ret$med=median(x);
	return(ret);
}

cat("Working on ", InputFileName, "\n", sep="");

# Load summary_table.xls
mat<-as.matrix(read.table(InputFileName, sep=",", header=TRUE, check.names=FALSE))

disttype=mat[,4];
rmsd=as.numeric(mat[,5]);
maxtaxa=as.numeric(mat[,6]);

lnorm_idc=(disttype=="lnorm")
gamma_idc=(disttype=="gamma")
pareto_idc=(disttype=="pareto")
frechet_idc=(disttype=="frechet")

lnorm=median(rmsd[lnorm_idc]);
gamma=median(rmsd[gamma_idc]);
pareto=median(rmsd[pareto_idc]);
frechet=median(rmsd[frechet_idc]);

lnorm_range=confid(maxtaxa[lnorm_idc]);
gamma_range=confid(maxtaxa[gamma_idc]);
pareto_range=confid(maxtaxa[pareto_idc]);
frechet_range=confid(maxtaxa[frechet_idc]);

fh=file(paste(OutputFileRoot, ".ci", sep=""), "wt");

cat("Type\tSmax (95% CI)\tRMSD\n", file=fh);
cat("lnorm\t", lnorm_range$med, " (", lnorm_range$lb, ", ", lnorm_range$ub, ")\t", lnorm, "\n", sep="", file=fh);
cat("gamma\t", gamma_range$med, " (", gamma_range$lb, ", ", gamma_range$ub, ")\t", gamma, "\n", sep="", file=fh);
cat("pareto\t", pareto_range$med, " (", pareto_range$lb, ", ", pareto_range$ub, ")\t", pareto, "\n", sep="", file=fh);
cat("frechet\t", frechet_range$med, " (", frechet_range$lb, ", ", frechet_range$ub, ")\t", frechet, "\n", sep="", file=fh);

###############################################################################

cat("Done.\n");
q(status=0);
