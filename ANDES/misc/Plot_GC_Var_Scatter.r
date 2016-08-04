#!/usr/bin/env Rscript

###############################################################################
#                                                                             # 
#       Copyright (c) 2009 J. Craig Venter Institute.                         #     
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

progname <- commandArgs(FALSE)[4]
args <- commandArgs(TRUE)

arg_count=length(args);

if(arg_count<2){

	script_name <- unlist(strsplit(progname,"="))[2]
	usage <- paste (
		"\nUsage:\n\t", script_name, "<outputfilename root> <scatter plot data1>... <scatter plot datan>\n",
		"Plots multiple scatter plots.",
		"\n")

	writeLines(usage)
	writeLines("Input FileName not defined.\n")
	quit(status=0)
}

###############################################################################

OutputFilePDF=paste(args[1], ".scat.pdf", sep="")
pdf(OutputFilePDF, width=11, height=8.5);

num_files=arg_count-1;

rain=rainbow(num_files);

infile_cnt=1;

for(filecount in 2:arg_count){
	filename=args[filecount];
	writeLines(paste("Plotting: ", filename, sep=""));
	writeLines(paste("Color: ", toString(rain[infile_cnt]), sep=""));
	data=read.delim(filename);
	if(filecount == 2){
		plot(data, col=rain[infile_cnt], xlim=c(0,1), ylim=c(0,1));
	}else{
		points(data, col=rain[infile_cnt]);
	}
	legend(0.0, 1-(infile_cnt/30), filename, bty="n", text.col=rain[infile_cnt]);
	infile_cnt=infile_cnt+1;
}


#len=length(A[,1])
#num_sets=ncol(A[1,])
#
#for(x in 2:num_sets){
#	p=A[,x];
#	
#	xbuf=rep(0,len);
#	ybuf=rep(0,len);
#	non_zero=1;
#	for(i in 1:len){
#		if(p[i]>0){
#			xbuf[non_zero]=i;
#			ybuf[non_zero]=p[i];
#			non_zero=non_zero+1;
#		}
#	}
#	xx=xbuf[1:non_zero-1];
#	yy=ybuf[1:non_zero-1]
#
#	points(xx,yy, col=rain[x], pch=1, cex=.6);
#	#points(xx,yy, col=rain[x], pch=20);
#}
#
#


#q(status=0)
