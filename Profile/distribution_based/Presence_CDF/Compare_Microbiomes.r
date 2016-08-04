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

library('getopt');

params=c(
	"input_A_file", "a", 1, "character",
	"input_B_file", "b", 1, "character",
	"target_taxa_list", "t", 2, "character",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"	-a <input summary_table.xls file A>\n",
	"	-b <input summary_table.xls file B>\n",
	"	[-t <target taxa of interest list>]\n",
	"	[-o <output file name>]\n",
	"\n",
	"Reads in two summary files and uses the Kolmogrov-Smirnov test to determine if there \n",
	"is a difference in the ubiquity of two organisms between the two groups.\n",
	"\n",
	"If you specify the -t option, then only the taxa in the target taxa list will be displayed/output.\n",
	"\n",
	"Produces a (.compare.pdf) PDF file, with one page per taxa:\n",
	"	Contains a histogram of abundances for each group, and then a QQ plot between abundances.\n",
	"	Uses KS test to determine if the abundance of two taxa are difference.\n",
	"	Also contains a abundance vs. ubiquity plot for each group, and QQ blot between ubiquities.\n",
	"\n",
	"Generates a (.sig_dif.txt) text file:\n",
	"	Contains the bonferroni corrected p-value for the KS test between the two groups\n",
	"\n",
	"Generates a (.ubiq_diff.txt) text file:\n",
	"	Contains a list of max differences between the two groups. For down stream comparisons.\n",
	"\n",
	"Generates a (.uncor_names.txt) text file:\n",
	"	Contains a list of names that were significant with a uncorrect p-value of >0.05.\n",
	"\n",
	"\n", sep="");

if(!length(opt$input_A_file) || !length(opt$input_B_file)){
	cat(usage);
	q(status=-1);
}

InputFileAName=opt$input_A_file;
InputFileANameClean=gsub(".summary_table.xls", "", InputFileAName);
InputFileBName=opt$input_B_file;
InputFileBNameClean=gsub(".summary_table.xls", "", InputFileBName);

RootA=gsub(".summary_table.xls", "", InputFileAName);
RootB=gsub(".summary_table.xls", "", InputFileBName);
OutputRoot=paste(RootA, "_vs_", RootB, sep="");

if(length(opt$output_file)){
	OutputRoot=opt$output_file;
}

TargetTaxaList="";
if(length(opt$target_taxa_list)){
	TargetTaxaList=opt$target_taxa_list;
}

cat("Input File A: ", InputFileAName, "\n");
cat("Input File B: ", InputFileBName, "\n");
cat("Target Taxa List: ", TargetTaxaList, "\n");
cat("Output File Root: ", OutputRoot, "\n");


################################################################################

load_summary_table=function(filename){
	st=as.matrix(read.table(filename, header=TRUE, sep="\t", row.names=1, check.names=F));
	
	# Only return values taxa that are greater than zero
	counts=st[,2:ncol(st)];
	sum=apply(counts, 2, sum);
	return(counts[,sum>0]);
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

read_list=function(filename){
	taxa=as.character(as.matrix(read.table(filename, sep="\t"))[,1]);
	num_taxa=length(taxa);
	cat("Num Taxa Loaded: ", num_taxa, "\n");
	return(taxa);
}


################################################################################

st_A=load_summary_table(InputFileAName);
st_B=load_summary_table(InputFileBName);

st_A=normalize(st_A);
st_B=normalize(st_B);

if(TargetTaxaList!=""){
	target_taxa_list=read_list(TargetTaxaList);
	list_length=length(target_taxa_list);
	cat("Target Taxa List Length: ", list_length, "\n");

	# Get number of rows
	num_samples_inA=nrow(st_A);
	num_samples_inB=nrow(st_B);

	# Get available columns
	inA=colnames(st_A);
	inB=colnames(st_B);
	
	# Create zero count columns
	notInA=setdiff(target_taxa_list, inA);
	notInB=setdiff(target_taxa_list, inB);

	# Give taxa names to zero count columns
	zerosA=matrix(0, nrow=num_samples_inA, ncol=length(notInA));
	colnames(zerosA)=notInA;
	zerosB=matrix(0, nrow=num_samples_inB, ncol=length(notInB));
	colnames(zerosB)=notInB;

	# Attach zero count columns to dataset
	st_A=cbind(st_A, zerosA);
	st_B=cbind(st_B, zerosB);
	
	# Subset out taxa we want
	st_A=st_A[,target_taxa_list];
	st_B=st_B[,target_taxa_list];

	#print(st_A);
	#print(st_B);
}

num_taxa_A=ncol(st_A);
num_taxa_B=ncol(st_B);

cat("Num Taxa: \n");
cat("	File A: ", num_taxa_A, "\n");
cat("	File B: ", num_taxa_B, "\n");

num_samples_A=nrow(st_A);
num_samples_B=nrow(st_B);

cat("Num Samples: \n");
cat("	File A: ", num_samples_A, "\n");
cat("	File B: ", num_samples_B, "\n");


################################################################################

taxa_list_A=colnames(st_A);
taxa_list_B=colnames(st_B);

shared_taxa=intersect(taxa_list_A, taxa_list_B);
num_shared=length(shared_taxa);
only_A_taxa=setdiff(taxa_list_A, shared_taxa);
only_B_taxa=setdiff(taxa_list_B, shared_taxa);
cat("Num Shared: ", num_shared, "\n");
all_taxa=union(taxa_list_A, taxa_list_B);
num_all_taxa=length(all_taxa);

################################################################################

fh=file(paste(OutputRoot, ".exc_dif.txt", sep=""), "w");
pdf(paste(OutputRoot, ".compare.pdf", sep=""), height=8.5, width=11);

cat(file=fh, "A: ", InputFileAName, "\n");
cat(file=fh, "B: ", InputFileBName, "\n");
cat(file=fh, "\n");

cat(file=fh, "Exclusive in A:\n");
cat(file=fh, paste(only_A_taxa, collapse="\n"));
cat(file=fh, "\n");
cat(file=fh, "\n");

cat(file=fh, "Exclusive in B:\n");
cat(file=fh, paste(only_B_taxa, collapse="\n"));
cat(file=fh, "\n");
cat(file=fh, "\n");

################################################################################

par(mfrow=c(2,3));
par(oma=c(0,0,4,0));

uncorrected_pvalue=numeric(num_all_taxa);

num_A_samples=nrow(st_A);
num_B_samples=nrow(st_B);

all_ubiquity_diff=numeric(num_all_taxa);
average_abundance=numeric(num_all_taxa);

for(i in 1:num_all_taxa){
	taxa_name=all_taxa[i];
	#cat("Working on: ", taxa_name, "\n");
	#print(st_A[,taxa_name]);
	#print(st_B[,taxa_name]);

	if(!(any(taxa_name==taxa_list_A))){
		A_val=rep(0, num_A_samples);
	}else{
		A_val=st_A[,taxa_name];
	}	

	if(!(any(taxa_name==taxa_list_B))){
		B_val=rep(0, num_B_samples);
	}else{
		B_val=st_B[,taxa_name];
	}	

	# Compute average abundance in A and in B, then average them together
	A_mean_ab=mean(A_val);
	B_mean_ab=mean(B_val);
	average_abundance[i]=(A_mean_ab+B_mean_ab)/2;

	# Compute KS test
	ks_result=ks.test(A_val, B_val);
	#cat("p-value: ", ks_result$p.value, "\n");
	uncorrected_pvalue[i]=ks_result$p.value;
	bon_cor=uncorrected_pvalue[i]*num_all_taxa;
	
	# Plot abundances in A and in B
	comb_hist=hist(c(A_val, B_val), plot=FALSE);
	a_hist=hist(A_val, breaks=comb_hist$breaks, plot=FALSE);
	b_hist=hist(B_val, breaks=comb_hist$breaks, plot=FALSE);
	max_counts=max(a_hist$counts, b_hist$counts);
	hist(A_val, ylim=c(0,max_counts), main="", breaks=comb_hist$breaks, xlab="Abundances in A");
	hist(B_val, ylim=c(0,max_counts), main="", breaks=comb_hist$breaks, xlab="Abundances in B");
	max_abundance=max(c(A_val, B_val));

	# Plot QQ plot of A Abund vs B Abund
	qq=qqplot(A_val, B_val, plot.it=FALSE);
	qqlim=max(c(qq$x,qq$y))
	plot(qq$x, qq$y, xlim=c(0,qqlim*1.1), ylim=c(0,qqlim*1.1), type="n", xlab="A Abundance", ylab="B Abundance");
	mtext(sprintf("p-value: %3.3f", ks_result$p.value), cex=.6);
	mtext(sprintf("Bonferroni Corrected p-value: %3.3f", bon_cor), cex=.5, line=-.75 );
	abline(a=0,b=1, col="grey", lwd=3, lty=1);	# Mark equal line
	points(qq$x, qq$y, type="b", cex=.3);		# Draw qq lines

	# Plot Abundance vs Ubiquity
	num_ticks=100;
	abundance_ticks=seq(0, max_abundance, length.out=num_ticks);
	cdf_A=numeric(0);
	cdf_B=numeric(0);
	for(j in 1:num_ticks){
		cdf_A[j]=sum(A_val>abundance_ticks[j])/num_A_samples;
		cdf_B[j]=sum(B_val>abundance_ticks[j])/num_B_samples;
	}

	max_diff=max(abs(cdf_A-cdf_B));	
	all_ubiquity_diff[i]=max_diff;

	plot(abundance_ticks, cdf_A, type="l", ylab="Ubiquity of A", xlab="Abundance of A", ylim=c(0,1));
	plot(abundance_ticks, cdf_B, type="l", ylab="Ubiquity of B", xlab="Abundance of B", ylim=c(0,1));

	# Plot QQ plot of A ubiquity vs B ubiquity
	plot(cdf_A, cdf_B, xlab="Ubiquity of A", ylab="Ubiquity of B", type="n", ylim=c(0,1), xlim=c(0,1));
	mtext(sprintf("Max Diff: %3.3f", max_diff), cex=.6);
	abline(a=0,b=1, col="grey", lwd=3, lty=1);	# Mark equal line
	points(cdf_A, cdf_B, type="b", cex=.3);

	# Label plots with taxa name
	if(bon_cor<.05){
		color="red";
	}else{
		color="black";
	}
	mtext(taxa_name, cex=.85, font=2, outer=TRUE, col=color);
	
	#cat("\n\n");

}

par(mfrow=c(1,1));
median_diff=median(all_ubiquity_diff);
diff_hist=hist(all_ubiquity_diff, xlab="Ubiquity Differences", main="Distribution of Ubiquity Differences", breaks=20);
abline(v=median_diff, col="blue", lty=2);
text(median_diff, max(diff_hist$counts), labels=sprintf("median = %3.3f", median_diff), col="blue", pos=4);

################################################################################
# Output B&H corrected p-values

ix=order(uncorrected_pvalue, decreasing=FALSE);
sorted_taxa_name=all_taxa[ix];
sorted_uncor_pvalue=uncorrected_pvalue[ix];

bh_qval=numeric(num_all_taxa);
bon_qval=numeric(num_all_taxa);

bh_qval=p.adjust(sorted_uncor_pvalue, method="BH");
bon_qval=p.adjust(sorted_uncor_pvalue, method="bonferroni");

fh=file(paste(OutputRoot, ".p-values.txt", sep=""), "w");
cat(file=fh, "Uncor\tFDR\tBonferroni\tTaxa Name\n");
for(i in 1:num_all_taxa){
	cat(file=fh, sorted_uncor_pvalue[i], "\t", bh_qval[i], "\t", bon_qval[i], "\t", sorted_taxa_name[i], "\n");
}
close(fh);

################################################################################
# Output differences in ubiquity across all taxa
fh=file(paste(OutputRoot, ".ubiq_diff.txt", sep=""), "w");
cat(file=fh, OutputRoot, ",", sep="");
cat(file=fh, paste(sprintf("%3.4f",all_ubiquity_diff), collapse=","),"\n");
close(fh);

################################################################################
# Output Weighted KS-Statistic

weighted_ks_stat=sum(average_abundance*all_ubiquity_diff);

fh=file(paste(OutputRoot, ".weighted_ks.txt", sep=""), "w");
cat(file=fh, "OutputRoot", "\t", "Abundance_Weighted_KS_Statistic", "\t", weighted_ks_stat, "\n", sep="");
close(fh);

################################################################################
# Output list of names that were significant before correcting

fh=file(paste(OutputRoot, ".uncor_pval_lt05", sep=""), "w");

uncor_sig_names=sorted_taxa_name[sorted_uncor_pvalue<=0.05];
if(length(uncor_sig_names)>0){
	for(i in 1:length(uncor_sig_names)){
		cat(file=fh, uncor_sig_names[i], "\n", sep="");
	}
}

close(fh);

fh=file(paste(OutputRoot, ".uncor_pval_lt10", sep=""), "w");

uncor_sig_names=sorted_taxa_name[sorted_uncor_pvalue<=0.10];
if(length(uncor_sig_names)>0){
	for(i in 1:length(uncor_sig_names)){
		cat(file=fh, uncor_sig_names[i], "\n", sep="");
	}
}

close(fh);

################################################################################

cat("Done.\n")
print(warnings());

q(status=0)
