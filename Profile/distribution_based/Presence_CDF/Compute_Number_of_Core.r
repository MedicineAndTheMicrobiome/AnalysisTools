#!/usr/bin/env Rscript

###############################################################################
#                                                                             #
#       Copyright (c) 2013 J. Craig Venter Institute.                         #
#       All rights reserved.                                                  #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                             #
###############################################################################
###############################################################################

ALPHA=.05

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 2, "character",
	"ubiquity_cutoff", "u", 1, "numeric",
	"abundance_cutoff", "a", 1, "numeric",
	"alpha", "c", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"	-i <input summary_table.xls file>\n",
	"	-u <minimum ubiquity cutoff (between 0.0 and 1.0)>\n",
	"	-a <minimum abundance cutoff (between 0.0 and 1.0)>\n",
	"	[-c <confidence interval (1-alpha)*100%, default=", ALPHA, ">]\n",
	"	[-o <output core numbers filename>]\n",
	"\n",
	"For the specified summary table, computes the number of\n",
	"core taxa given the specified ubiquity and abundance cutoff\n",
	"The median and specified confidence intervals are also reported.\n",
	"\n",
	"These were computed by bootstrapping the donors and reads per sample.\n",
	"\n",
	"\n", sep="");

if(!length(opt$input_file) || !length(opt$ubiquity_cutoff) || !length(opt$abundance_cutoff)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputRoot=gsub(".summary_table.xls", "", InputFileName);

if(length(opt$output_file)){
	OutputRoot=opt$output_file;
}

if(length(opt$alpha)){
	Alpha=opt$alpha;
}else{
	Alpha=ALPHA;
}
# 2 times the minimum we need
NumBoots=max(c((2/Alpha)*2, 2000));

UbiquityCutoff=opt$ubiquity_cutoff;
AbundanceCutoff=opt$abundance_cutoff;

cat("Input File: ", InputFileName, "\n");
cat("Output File Root: ", OutputRoot, "\n");
cat("Ubiquity Cutoff: ", UbiquityCutoff, "\n");
cat("Abundance Cutoff: ", AbundanceCutoff, "\n");
cat("\n");
cat("Alpha: ", Alpha, "\n");
cat("Num bootstraps to perform: ", NumBoots, "\n");

################################################################################

load_summary_table=function(filename){
	st=as.matrix(read.table(filename, header=TRUE, sep="\t", row.names=1, check.names=F));
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

compute_num_core=function(nst, ubiq_cutoff, abund_cutoff){

	above_abund=matrix(0, nrow=nrow(nst), ncol=ncol(nst));
	above_abund[]=nst>abund_cutoff;
	ubiq=apply(above_abund, 2, sum)/nrow(nst);
	is_core=ubiq>ubiq_cutoff;
	num_core=sum(is_core);

	#print(nst);
	#print(above_abund);
	#print(ubiq)
	#print(is_core);
	#print(num_core);
	core_rec=list();
	core_rec[["num_core"]]=num_core;
	core_rec[["members"]]=is_core;
	return(core_rec);
}

#-------------------------------------------------------------------------------

ci = function(x, alpha){
        n=length(x);
        med=median(x);
        ba=1-(n-2)/n;
        #cat("Best alpha = ", ba, "\n");
        if(ba <= (alpha+.0000000000000001)){
                sorted=sort(x);
                lb=sorted[floor(n*(alpha/2))+1];
                ub=sorted[ceiling(n*(1-(alpha/2)))];
                return(c(med,lb,ub));
        }else{
                return(c(med,NA,NA))
        }
}

################################################################################

resample_summary_table=function(normalized, counts){
        num_samp=nrow(normalized);
        seq=sample(1:num_samp, replace=T);
        new_counts=numeric();
        for(i in seq){
                new_counts=cbind(new_counts, rmultinom(1, counts[i], prob=normalized[i,]));
        }
        colnames(new_counts)=sprintf("%s.%i", rownames(normalized)[seq], 1:num_samp);
        return(t(new_counts));
}

################################################################################

# Load summary table and get basic info
cat("Loading Summary Table: ", InputFileName, "\n");
st=load_summary_table(InputFileName);
num_taxa=ncol(st);
num_samples=nrow(st);
taxa_names=colnames(st);
sample_names=rownames(st);

cat("\n");
cat("Num Taxa: ", num_taxa, "\n");
cat("Num Samples: ", num_samples, "\n");
cat("\n");

# Normalize counts
nst=normalize(st);
sample_tot_cts=apply(st, 1, sum);

# Bootstrap
num_core=numeric(NumBoots);
core_members=rep(0, num_taxa);
for(i in 1:NumBoots){
	
	#cat("Iteration ", i, " of ", NumBoots, "\n", sep="");

	instance=resample_summary_table(nst, sample_tot_cts);
	n_instance=normalize(instance);

	core_rec=compute_num_core(n_instance, UbiquityCutoff, AbundanceCutoff);

	num_core[i]=core_rec$num_core;
	core_members=core_members+core_rec$members;
}

prop_core_members=core_members/NumBoots;
core_idx=which(prop_core_members>=(1-Alpha));

cat("Significant Core Taxa: \n");
print(taxa_names[core_idx]);
cat("\n");

# Compute observed core
num_instance_core_rec=compute_num_core(nst, UbiquityCutoff, AbundanceCutoff);
num_instance_core=num_instance_core_rec$num_core;

# Compute confidence intervals
#print(num_core);
confidence=ci(num_core, Alpha);
med=confidence[1];
lb=confidence[2];
ub=confidence[3];

cat("Num observed core (in input): ", num_instance_core, "\n");
cat("Bootstrapped results:\n");
cat("  Median Core: ", med, "\n");
cat("  Lowerbound: ", lb, "\n");
cat("  Upperbound: ", ub, "\n");

################################################################################
# Output core statistics into file

fh=file(paste(OutputRoot, ".core_statistics.tsv", sep=""), "w");
cat(file=fh, paste("FilenameRoot", "ObsNumCore", "Median", "LowerBound", "UpperBound", "alpha", "UbCutoff", "AbCutoff", sep="\t"), "\n");
cat(file=fh, paste(OutputRoot, num_instance_core, med, lb, ub, Alpha, UbiquityCutoff, AbundanceCutoff, sep="\t"), "\n");
close(fh);


fh=file(paste(OutputRoot, ".core_members.tsv", sep=""), "w");
cat(file=fh, paste("Index", "MemberName", "Pr", sep="\t"), "\n");
for(i in 1:length(core_idx)){
	ix=core_idx[i];
	cat(file=fh, paste(i, taxa_names[ix], prop_core_members[ix], sep="\t"), "\n");
}
close(fh);

# Observed Core members
num_obs_core=num_instance_core_rec$num_core;
fh=file(paste(OutputRoot, ".observed_core.tsv", sep=""), "w");
obs_core_members=taxa_names[num_instance_core_rec$members];
for(i in 1:num_obs_core){
        cat(file=fh, obs_core_members[i], "\n");
}
close(fh);

################################################################################

cat("\nDone.\n")

q(status=0)
