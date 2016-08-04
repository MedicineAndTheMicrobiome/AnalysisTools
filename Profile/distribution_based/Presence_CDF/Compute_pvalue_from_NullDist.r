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
	"null_distribution", "n", 1, "character",
	"test_statistics", "t", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

PaperHeight=11;
PaperWidth=8.5;

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"	-n <null distribution file>\n",
	"	-t <test statistics file>\n",
	"\n",
	"This script will read in a null distribution and a test statistics file\n",
	"and compute the p-value for each test statistic.\n",
	"\n",
	"For the 1-tail, only testing for the upper tail is performed.\n",
	"\n", sep="");

if(!length(opt$null_distribution) || !length(opt$test_statistics)){
	cat(usage);
}


InputNullDistribution=opt$null_distribution;
InputTestStatistics=opt$test_statistics;
TwoTailedFlag=FALSE;

if(length(opt$two_tailed)){
	TwoTailedFlag=TRUE;
}

cat("Input Null Distribution: ", InputNullDistribution, "\n");
cat("Input Test Statistics: ", InputTestStatistics, "\n");
cat("Two Tailed?", TwoTailedFlag, "\n");
cat("\n");

################################################################################

load_values=function(filename){
	values=as.vector(as.matrix(read.table(filename, sep="\t")));
	return(values);
}

################################################################################

null_dist=load_values(InputNullDistribution);
num_null_dist_val=length(null_dist);

test_stats=load_values(InputTestStatistics);
num_test_stat_val=length(test_stats);

cat("Num null distribution values: ", num_null_dist_val, "\n");
cat("Num test statistics: ", num_test_stat_val, "\n");

hist_rec=hist(null_dist, plot=FALSE);
null_median=median(null_dist);
max_counts=max(hist_rec$counts);

null_sorted=sort(null_dist);

outfile=file(paste(InputTestStatistics, ".pval", sep=""), "w");

cat(file=outfile, "#Index\tTestStatistic\t1-tailed\t2-tailed\n");

pdf(paste(InputTestStatistics, ".pval.pdf", sep=""), height=11, width=8.5);
par(mfrow=c(2,3));

for(i in 1:num_test_stat_val){
	cur_test_stat=test_stats[i];
	
	num_greater_than_ts=sum(null_sorted>cur_test_stat);

	pval=num_greater_than_ts/num_null_dist_val;
	
	one_tailed_pval=pval;

	if(pval>.5){
		two_tailed_pval=(1-pval)*2;
	}else{
		two_tailed_pval=pval*2;
	}

	if(pval==0){
		one_tailed_pval=sprintf("<%f", 1/num_null_dist_val);
		two_tailed_pval=sprintf("<%f", 2/num_null_dist_val);
	}	

	cat(file=outfile, i, "\t", cur_test_stat, "\t", one_tailed_pval, "\t", two_tailed_pval, "\n");

	r=range(c(cur_test_stat, hist_rec$breaks));
	plot(hist_rec, main=sprintf("Test Statistic = %3.3f", cur_test_stat), xlim=r, xlab="Null Distribution");
	abline(v=null_median, col="blue", lty=2);
	abline(v=cur_test_stat, col="red");

}

################################################################################

w=warnings();
if(length(w)){
	print(w);
}
cat("Done.\n")

q(status=0)
