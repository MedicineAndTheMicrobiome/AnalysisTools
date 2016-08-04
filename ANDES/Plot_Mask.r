#!/usr/bin/env Rscript

###############################################################################
#                                                                             # 
#       Copyright (c) 2011 J. Craig Venter Institute.                         #     
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
		"input_mask", "m", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage <- paste (
		"\nUsage:\n", script_name, "\n",
		"	-m <input file>\n",
		"\n",
		"Plots a mask so you can see how the weights are distributed.\n",
		"\n");

if(!length(opt$input_mask)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_mask;

MaskPlotPDF=paste(InputFileName, ".pdf", sep="")

###############################################################################

is.polar_pos=function(x){
	return((x=="R" | x=="H" | x=="K"));
}
is.polar_neg=function(x){
	return((x=="D" | x=="E"));
}
is.polar_neut=function(x){
	return((x=="N" | x=="C" | x=="Q" | x=="S" | x=="T" | x=="Y"));
}
is.nonpolar=function(x){
	return((x=="A" | x=="G" | x=="I" | x=="L" | x=="M" | x=="F" | x=="P" | x=="V" | x=="W"));
}

mole_weight=function(x){
	# Molecular weights
	if(x=="I"){ return(131.1736)}
	if(x=="L"){ return(131.1736)}
	if(x=="K"){ return(146.1882)}
	if(x=="M"){ return(149.2124)}
	if(x=="F"){ return(165.1900)}
	if(x=="T"){ return(119.1197)}
	if(x=="W"){ return(204.2262)}
	if(x=="V"){ return(117.1469)}
	if(x=="R"){ return(174.2017)}
	if(x=="H"){ return(155.1552)}
	if(x=="A"){ return(89.0935)}
	if(x=="N"){ return(132.1184)}
	if(x=="D"){ return(133.1032)}
	if(x=="C"){ return(121.1590)}
	if(x=="E"){ return(147.1299)}
	if(x=="Q"){ return(146.1451)}
	if(x=="G"){ return(75.0669)}
	if(x=="P"){ return(115.1310)}
	if(x=="S"){ return(105.0930)}
	if(x=="Y"){ return(181.1894)}
	cat("Error!  Unknown AA. :", x, "\n");
}

######################################################################################################
# Load data
mask=as.matrix(read.table(InputFileName, header=TRUE, check.names=FALSE));

residues=mask[,1];
values=as.numeric(mask[,2]);
mask_len=nrow(mask);

cat("Length of Mask: ", mask_len, "\n");

#-----------------------------------------------------------------------------------------------------

# Assign colors
#print(residues);
colors=character(mask_len);
colors[is.polar_pos(residues)]="red";
colors[is.polar_neg(residues)]="blue";
colors[is.polar_neut(residues)]="purple";
colors[is.nonpolar(residues)]="mediumpurple";

# Assign point sizes
sizes=numeric(mask_len);
for(i in 1:mask_len){
	sizes[i]=mole_weight(residues[i]);
}
max_size=max(sizes);
sizes=sizes/max_size;
#print(sizes);

#-----------------------------------------------------------------------------------------------------

# Generate plot
pdf(MaskPlotPDF, height=8.5, width=11);
plot(values, pch=15, col=colors, cex=sizes, xlab="Position", ylab="Weight");
dev.off()

######################################################################################################

cat("Done.\n");
