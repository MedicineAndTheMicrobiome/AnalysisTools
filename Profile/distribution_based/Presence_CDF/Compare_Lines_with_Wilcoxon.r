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
	"stats", "s", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"	-s <Filename containing statistics to compare medians of.>\n",
	"\n",
	"Reads in a list of groups to compare.\n",
	"\n",
	"Format is: \n",
	" <group name1>,<test statistics (comma separated)>\\n\n",
	" <group name2>,<test statistics (comma separated)>\\n\n",
	" ...\n",
	"\n", sep="");

if(!length(opt$stats)){
	cat(usage);
	q(status=-1);
}

InputStatsFile=opt$stats;
OutputFileName=gsub(".txt", "", opt$stats);

cat("Input Stats File: ", InputStatsFile, "\n");
cat("Output Stats File: ", OutputFileName, "\n");

################################################################################

fh=file(InputStatsFile);
stats=readLines(fh);
close(fh);

num_groups=length(stats);
cat("Num Groups to compare: ", num_groups, "\n", sep="");

stat_list=list();
names=character();
for(i in 1:num_groups){
	array=strsplit(stats[i], ",")[[1]];
	stat_list[[i]]=as.numeric(array[2:length(array)]);
	names[i]=array[1];
}

################################################################################

plot_compare_two=function(a_data, b_data, a_name, b_name){
	comb_hist=hist(c(a_data, b_data), plot=FALSE);
	a_hist=hist(a_data, breaks=comb_hist$breaks, plot=FALSE);
	b_hist=hist(b_data, breaks=comb_hist$breaks, plot=FALSE);
	max_counts=max(a_hist$counts, b_hist$counts);

	par(mar=c(5.1,4.1,2.1,2.1));

	hist(a_data, ylim=c(0,max_counts*1.1), main="", breaks=comb_hist$breaks, xlab=a_name);
	a_median=median(a_data);
	abline(v=a_median, col="blue", lty=2);
	text(a_median, max_counts, sprintf("median=%3.4f", a_median), pos=3, cex=.8);

	hist(b_data, ylim=c(0,max_counts*1.2), main="", breaks=comb_hist$breaks, xlab=b_name);
	b_median=median(b_data);
	abline(v=b_median, col="blue", lty=2);
	text(b_median, max_counts, sprintf("median=%3.4f", b_median), pos=3, cex=.8);
}

################################################################################

pdf(paste(OutputFileName, ".wilcox_comp.pdf", sep=""), height=11, width=8.5);

pval_matrix=matrix(0, nrow=num_groups, ncol=num_groups);
rownames(pval_matrix)=names;
colnames(pval_matrix)=names;
par(mfrow=c(1,3));
for(i in 1:num_groups){
	for(j in 1:num_groups){
		if(i<j){
			wc=wilcox.test(stat_list[[i]], stat_list[[j]]);
			pval_matrix[i,j]=wc$p.val;

			plot_compare_two(stat_list[[i]], stat_list[[j]], a_name=names[i], b_name=names[j]);

			med_i=median(stat_list[[i]]);
			med_j=median(stat_list[[j]]);

			plot(0,0, ylim=c(0,10), xlim=c(0,1), xlab="", ylab="", main="", yaxt="n", xaxt="n", bty="n", type="n");
			text(0,8, label=sprintf("\n\n\n\n\nWilcoxon Rank Sum Test:\n  p-value = %3.4f\n Left Median = %3.4f\n Right Median = %3.4f\n",
				 pval_matrix[i,j], med_i, med_j), 
				pos=4, cex=1);
		}
	}
}

################################################################################

cat("Done.\n")
print(warnings());

q(status=0)
