#!/usr/bin/env Rscript

###############################################################################

library('getopt');
library(plotrix);

params=c(
	"input_file_list", "i", 1, "character",
	"output_filename", "o", 2, "character",
	"svg", "s", 2, "logical",
	"landscape", "l", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n   ", script_name, "\n",
	"	-i <rarefaction curve file list>\n",
	"	[-o <output file name>]\n",
	"\n",
	"	[-s (output file an SVG file)]\n",
	"	[-l (output in Landscape orientation)]\n",
	"\n",
	"This script will read in one file.  It is a list of rarefaction curve files\n",
	"with medians, means, and 95% confidence intervals. \n",
	"\n",
	"The rows should go like this:\n",
	"\n",
	"	<Label>\\t<medians starting with 1 sample>\n",
	"	<Label>\\t<means starting with 1 sample>\n",
	"	<Label>\\t<95% CI lowerbound starting with 1 sample>\n",
	"	<Label>\\t<95% CI upperbound starting with 1 sample>\n",
	"\n",
	"\n",
	"\n", sep="");

InputFileNameList=opt$input_file_list;
OutputFilename=opt$output_filename;

if(length(opt$input_file_list)==0){
	cat(usage);
	quit(status=-1);
}

if(length(OutputFilename)==0){
	OutputFilename=tail(strsplit(InputFileNameList, "/")[[1]],1);
}

DoSVG=FALSE;
if(length(opt$svg)!=0){
	DoSVG=TRUE;	
	cat("Generating SVG instead of PDF.\n");
}

DoLandscape=FALSE;
if(length(opt$landscape)!=0){
	DoLandscape=TRUE;	
	cat("Generating output in Landscape (11 x 8.5) format instead of Portrait (8.5 x 11).\n");
}

#OutHeight=11;
#OutWidth=8.5;
OutHeight=11;
OutWidth=4.25;
if(DoLandscape){
	OutHeight=8.5;
	OutWidth=11;
}

cat("Output Filename Root: ", OutputFilename, "\n\n");
if(DoSVG){
	svg(paste(OutputFilename, ".1.svg", sep=""), height=OutHeight, width=OutWidth);
}else{
	pdf(paste(OutputFilename, ".pdf", sep=""), height=OutHeight, width=OutWidth);
}


###############################################################################
# Load counts from file

load_file_list=function(filename){
	cat("Loading list: ", filename, "...\n");
	dat=as.matrix(read.table(filename));
	cat("Done.\n");
	return(as.vector(dat[,1]));
}

load_rarefaction_curves=function(filename){
	cat("Loading rarefaction file: ", filename, "...\n");
	dat=as.matrix(read.table(filename, sep="\t", check.names=FALSE, row.names=1, fill=TRUE));
	rownames(dat)=c();
	return(dat);
}

load_curves_into_list=function(filename_list){
	rare_list=list();
	list_len=length(filename_list);
	cat("Number of files to load: ", list_len, "\n");
	#print(filename_list);
	for(i in 1:list_len){
		rare_list[[filename_list[i]]]=load_rarefaction_curves(filename_list[i]);
	}
	return(rare_list);
}

###############################################################################
# Load Genes/Taxa curves into memory

curves_filelist=load_file_list(InputFileNameList);
curves_list=load_curves_into_list(curves_filelist)
num_curves=length(curves_list);
cat("\nNum Curves Loaded: ", num_curves, "\n");

###############################################################################
# Compute and assign color scheme

get_colors=function(num_groups){
	# Allocate colors
        rainbow_natural=rainbow(num_groups)
        rainbow_lowval=rainbow(num_groups,v=.5)
        rainbow_lowsat=rainbow(num_groups,s=.75)
        odd=seq(1,num_groups,2);
        even=seq(2,num_groups,2);
        rainbow_colors=character();
        rainbow_colors[odd]=rainbow_natural[odd];
        rainbow_colors[even]=rainbow_lowsat[even];
        yll=round((.5/6)*(num_groups))+1;
        yul=round((1.5/6)*(num_groups))+1;
        yellows=seq(yll,yul,1);
        rainbow_colors[yellows]=rainbow_lowval[yellows]

        # Shuffle colors
        matrix_dim=ceiling(sqrt(num_groups));
        shuf_vect=rep(0,matrix_dim^2);
        shuf_vect[1:num_groups]=1:num_groups;
        shuf_vect=as.vector(t(matrix(shuf_vect,ncol=matrix_dim)));
        shuf_vect=shuf_vect[shuf_vect>0];
        rainbow_colors=rainbow_colors[shuf_vect];

	return(rainbow_colors);
}

cat("Number of Colors Needed:", num_curves, "\n");
color_scheme=get_colors(num_curves);

###############################################################################
# Compute max samples/discovered

max_discovered=0;
min_discovered=Inf;
max_sampled=0;
num_samples_in_curve=numeric();
for(i in 1:num_curves){
	upperbound_curve=(curves_list[[i]][4,]); # 4 is 4th row (upperbound) in matrix 
	lowerbound_curve=(curves_list[[i]][3,]); # 3 is 3th row (lowerbound) in matrix 
	median_curve=(curves_list[[i]][1,]); # 1 is 1th row (median) in matrix 
	num_samples_in_curve[i]=length(upperbound_curve);

	max_of_curve=max(upperbound_curve);
	if(max_discovered<max_of_curve){
		max_discovered=max_of_curve;
	}
	
	min_of_curve=min(lowerbound_curve);
	if(min_discovered>min_of_curve){
		min_discovered=min_of_curve;
	}
}
max_sampled=max(num_samples_in_curve);
min_sampled=min(num_samples_in_curve);

cat("Max Discovered (among all): ", max_discovered, "\n");
cat("Max Sampled (among all): ", max_sampled, "\n");
cat("Min Sampled (among all): ", min_sampled, "\n");

cat("\n");

###############################################################################
# Get the number discovered at min among all
num_discovered_at_min=numeric();
for(i in 1:num_curves){
	num_discovered_at_min[i]=curves_list[[i]][1,min_sampled]; # Median
}


min_discovered_at_min_sample=Inf;
for(i in 1:num_curves){
	min_discovered_at_min_sample=min(min_discovered_at_min_sample, curves_list[[i]][3,min_sampled]); # lowerbound
}

ordered_by_discovery=sort(num_discovered_at_min, decreasing=TRUE, index.return=TRUE);
obd=ordered_by_discovery$ix;
odc=numeric();
odc[obd]=1:num_curves;

###############################################################################
# Clean names
clean_name=character();
for(i in 1:num_curves){
	clean_name[i]=gsub("_", " ", tail(strsplit(curves_filelist[i],"/")[[1]],1));
	cat("\t", i, ": Cleaning name: ", curves_filelist[i], " -> ", clean_name[i], "\n");
}

###############################################################################

# Plot combined rarefaction curves
plot(0,0, type="n", xlim=c(0,max_sampled*1.1), ylim=c(0,max_discovered), ylab="Discovered", xlab="Number of Samples", main="Samples Taken vs. Units Discovered");
ebarwidth=.001;
for(i in 1:num_curves){

	median_curve=as.vector(curves_list[[i]][1,]); # median
	lb_curve=as.vector(curves_list[[i]][3,]); # lowerbound
	ub_curve=as.vector(curves_list[[i]][4,]); # upperbound

	# Down sample the CI's to plot
	mask=-seq(i,num_samples_in_curve[i],num_curves);
	masked_ub_curve=ub_curve;
	masked_lb_curve=lb_curve;
	masked_ub_curve[mask]=median_curve[mask];
	masked_lb_curve[mask]=median_curve[mask];
	masked_lb_curve[num_samples_in_curve[i]]=lb_curve[num_samples_in_curve[i]];
	masked_ub_curve[num_samples_in_curve[i]]=ub_curve[num_samples_in_curve[i]];

	# Plot individual
	plotCI(1:num_samples_in_curve[i], median_curve, ui=masked_ub_curve, li=masked_lb_curve, add=TRUE, cex=.3, col=color_scheme[odc[i]], sfrac=ebarwidth);
	last_point=tail(median_curve,1);
	text(num_samples_in_curve[i], last_point, last_point, cex=.6, col=color_scheme[odc[i]], pos=4);
	#text(num_samples_in_curve[i], last_point, clean_name[i], cex=.6, col=color_scheme[odc[i]], pos=4);
	points(num_samples_in_curve[i], last_point, pch=16, col=color_scheme[odc[i]], cex=1);

}
legend(0, max_discovered, 
	legend=clean_name[obd], 
	fill=color_scheme,
	cex=.8
	);

###############################################################################

cat("\n");
if(DoSVG){
	svg(paste(OutputFilename, ".2.svg", sep=""), height=OutHeight, width=OutWidth);
}
# Plot combined rarefaction curves with labels
plot(0,0, type="n", xlim=c(0,max_sampled*1.1), ylim=c(min_discovered_at_min_sample*.9,max_discovered), ylab="Discovered", xlab="Number of Samples", main="Samples Taken vs. Units Discovered");
ebarwidth=.001;
for(i in 1:num_curves){

	median_curve=as.vector(curves_list[[i]][1,]); # median
	lb_curve=as.vector(curves_list[[i]][3,]); # lowerbound
	ub_curve=as.vector(curves_list[[i]][4,]); # upperbound

	# Down sample the CI's to plot
	mask=-seq(i,num_samples_in_curve[i],num_curves);
	masked_ub_curve=ub_curve;
	masked_lb_curve=lb_curve;
	masked_ub_curve[mask]=median_curve[mask];
	masked_lb_curve[mask]=median_curve[mask];
	#masked_lb_curve[num_samples_in_curve[i]]=lb_curve[num_samples_in_curve[i]];
	#masked_ub_curve[num_samples_in_curve[i]]=ub_curve[num_samples_in_curve[i]];

	# Plot individual
	plotCI(1:num_samples_in_curve[i], median_curve, ui=masked_ub_curve, li=masked_lb_curve, add=TRUE, cex=.3, col=color_scheme[odc[i]], sfrac=ebarwidth);
	last_point=tail(median_curve,1);
	#text(num_samples_in_curve[i], last_point, last_point, cex=.6, col=color_scheme[odc[i]], pos=1, offset=-.5);
	text(num_samples_in_curve[i], last_point, clean_name[i], cex=.6, col=color_scheme[odc[i]], pos=4, offset=.5);
	points(num_samples_in_curve[i], last_point, pch=16, col=color_scheme[odc[i]], cex=1);

}

####################################################################

max_discovered=log10(max_discovered);
min_discovered=log10(min_discovered);
min_discovered_at_min_sample=log10(min_discovered_at_min_sample);

if(DoSVG){
	svg(paste(OutputFilename, ".3.svg", sep=""), height=OutHeight, width=OutWidth);
}
plot(0,0, type="n", xlim=c(0,max_sampled*1.15), ylim=c(min_discovered,max_discovered), ylab="Log(Discovered)", xlab="Number of Samples", main="Samples Taken vs. Units Discovered");
ebarwidth=.001;
for(i in 1:num_curves){

	median_curve=as.vector(curves_list[[i]][1,]); # median
	lb_curve=as.vector(curves_list[[i]][3,]); # lowerbound
	ub_curve=as.vector(curves_list[[i]][4,]); # upperbound

	# Log base 10 transform
	median_curve=log10(median_curve);
	lb_curve=log10(lb_curve);
	ub_curve=log10(ub_curve);

	# Down sample the CI's to plot
	mask=-seq(i,num_samples_in_curve[i],num_curves);
	masked_ub_curve=ub_curve;
	masked_lb_curve=lb_curve;
	masked_ub_curve[mask]=median_curve[mask];
	masked_lb_curve[mask]=median_curve[mask];
	masked_lb_curve[num_samples_in_curve[i]]=lb_curve[num_samples_in_curve[i]];
	masked_ub_curve[num_samples_in_curve[i]]=ub_curve[num_samples_in_curve[i]];
	
	# Plot individual
	plotCI(1:num_samples_in_curve[i], median_curve, ui=masked_ub_curve, li=masked_lb_curve, add=TRUE, cex=.3, col=color_scheme[odc[i]], sfrac=ebarwidth);
	last_point=tail(median_curve,1);
	text(num_samples_in_curve[i], last_point, round(10^last_point), cex=.6, col=color_scheme[odc[i]], pos=4);
	points(num_samples_in_curve[i], last_point, pch=16, col=color_scheme[odc[i]], cex=1);

}
print(max_sampled);
print(max_discovered);
legend(max_sampled*.45, max_discovered*.5, 
	legend=clean_name[obd], 
	fill=color_scheme,
	cex=.7
	);

####################################################################

if(DoSVG){
	svg(paste(OutputFilename, ".4.svg", sep=""), height=OutHeight, width=OutWidth);
}
plot(0,0, type="n", xlim=c(0,max_sampled*1.35), ylim=c(min_discovered_at_min_sample*.9,max_discovered), ylab="Log(Discovered)", xlab="Number of Samples", main="Samples Taken vs. Units Discovered");
ebarwidth=.001;
for(i in 1:num_curves){

	median_curve=as.vector(curves_list[[i]][1,]); # median
	lb_curve=as.vector(curves_list[[i]][3,]); # lowerbound
	ub_curve=as.vector(curves_list[[i]][4,]); # upperbound

	# Log base 10 transform
	median_curve=log10(median_curve);
	lb_curve=log10(lb_curve);
	ub_curve=log10(ub_curve);

	# Down sample the CI's to plot
	mask=-seq(i,num_samples_in_curve[i],num_curves);
	masked_ub_curve=ub_curve;
	masked_lb_curve=lb_curve;
	masked_ub_curve[mask]=median_curve[mask];
	masked_lb_curve[mask]=median_curve[mask];
	#masked_lb_curve[num_samples_in_curve[i]]=lb_curve[num_samples_in_curve[i]];
	#masked_ub_curve[num_samples_in_curve[i]]=ub_curve[num_samples_in_curve[i]];
	
	# Plot individual
	plotCI(1:num_samples_in_curve[i], median_curve, ui=masked_ub_curve, li=masked_lb_curve, add=TRUE, cex=.3, col=color_scheme[odc[i]], sfrac=ebarwidth);
	last_point=tail(median_curve,1);
	#text(num_samples_in_curve[i]+1, last_point, 10^last_point, cex=.6, col=color_scheme[odc[i]], pos=1, offset=-.5);
	name_split=gsub(" ", "\n", clean_name[i]);
	text(num_samples_in_curve[i]+1, last_point, name_split, cex=.6, col=color_scheme[odc[i]], pos=4, offset=.5);
	points(num_samples_in_curve[i], last_point, pch=16, col=color_scheme[odc[i]], cex=1);

}

####################################################################

cat("Done.\n");
warn_msg=warnings();
if(!is.null(warn_msg)){
	print(warnings());
}
q(status=0);
