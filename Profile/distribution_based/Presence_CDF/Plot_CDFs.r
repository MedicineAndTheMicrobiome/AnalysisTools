#!/usr/bin/env Rscript

###############################################################################
PaperHeight=11;
PaperWidth=8.5;

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_file_root", "o", 2, "character",
	"color_map", "c", 2, "character",
	"height", "h", 2, "numeric",
	"width", "w", 2, "numeric",
	"ubiquity", "u", 2, "numeric",	
	"abundance", "a", 2, "numeric",
	"keep_list_file", "k", 2, "character",
	"ignore_bulls_eye", "b", 2, "logical",
	"grey_subcutoff_taxa", "r", 2, "logical",
	"id_map", "m", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"	-i <input CDF file>\n",
	"	[-u <ubiquity cutoff, eg. .8 is 80%>]\n",
	"	[-a <abundance cutoff, eg. .01 is 1%>]\n",
	"	[-c <color scheme, taxa to color map, recommended.>]\n",
	"	[-m <id map>]\n",
	"\n",
	"	[-k <keep list file (default all taxa)>]\n",
	"	[-r (flag to plot all taxa, but use grey and do not label taxa not exceeding cutoff)]\n",
	"	[-b (flag to ignore the plot of bulls eye in the graph)]\n",
	"\n",
	"	[-o <output filename root>]\n",
	"	[-h <paper height, default=", PaperHeight, " in>]\n",
	"	[-w <paper width, default=", PaperWidth, " in>]\n",
	"\n",
	"This script will generate a Ub-Ab Plot.\n",
	"The output file is a Adobe PDF file.\n",
	"\n",
	"Note:  The input file is a .cdf fle, not a summary table.\n",
	"\n",
	"Using the color scheme (-c) is highly recommended.  This can be generated with the\n",
	"Assign_Colors.r script.  This assigns a color to each taxa that you specify by creating\n",
	"a mapping of taxa name to RGB value.\n",
	"\n",
	"The -u and -a options, specify cutoffs for drawing the bull's eye.\n",
	"\n",
	"By default, only the taxa above both ubiquity and abundance cutoffs will be plotted.\n",
	"If the -r flag is set, then all taxa will be plotted, but only the taxa above the cutoff\n",	
	"	will be colored based on your specified color scheme.\n",
	"\n",
	"The keep list (-k) may be used to specify which taxa you want to draw.  You may want to use\n",
	"this if you have very many spurious taxa that you want to exclude.\n",
	"\n",
	
	"\n", sep="");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;

if(length(opt$output_file_root)){
	OutputFileRoot=opt$output_file_root;
}else{
	OutputFileRoot=gsub(".cdf", "", opt$input_file);
}

if(length(opt$height)){
	PaperHeight=opt$height;
}

if(length(opt$width)){
	PaperWidth=opt$width;
}

if(length(opt$ubiquity)){
	Ubiquity=opt$ubiquity;
}else{
	Ubiquity=0;
}

if(length(opt$abundance)){
	Abundance=opt$abundance;
}else{
	Abundance=0;
}

if(length(opt$keep_list_file)){
	KeepList=opt$keep_list_file;
}else{
	KeepList="";
}

if(length(opt$grey_subcutoff_taxa)){
	GreySubcutoffTaxa=TRUE;
}else{
	GreySubcutoffTaxa=FALSE;
}

if(length(opt$ignore_bulls_eye)){
        IgnoreBullsEye=TRUE;
}else{
        IgnoreBullsEye=FALSE;
}

ColorMap=NULL;
if(length(opt$color_map)){
	ColorMap=opt$color_map;
}

IDMap=NULL;
if(length(opt$id_map)){
	IDMap=opt$id_map;
}

cat("Input File List: ", InputFileName, "\n");
cat("Output File Root: ", OutputFileRoot, "\n");
cat("\n");
cat("Ubiquity Cutoff: ", Ubiquity, "\n");
cat("Abundance Cutoff: ", Abundance, "\n");
cat("Grey Subcutoff Taxa? ", GreySubcutoffTaxa, "\n");
cat("ignore bulls eye? ", IgnoreBullsEye, "\n");

cat("\nPaper Dimensions: h: ", PaperHeight, " x w: ", PaperWidth, "\n", sep="");

if(!(KeepList=="")){
	cat("Keep List: ", KeepList, "\n");
}

################################################################################

read_keeplist_from_file=function(keeplist_fn){
	keep_list=as.character(as.matrix(read.table(keeplist_fn, sep="\t", header=F)));
	return(as.vector(keep_list));
}

################################################################################

read_cdf_from_file=function(cdf_fn){
# Read precomputed CDF presence values from file

	cdf_info=list();
	cdfs=as.matrix(read.table(cdf_fn, sep=",", header=TRUE, row.names=1, check.names=FALSE, quote=""));
	
	cdf_info$fname=cdf_fn;
	cdf_info$cutoffs=as.numeric(colnames(cdfs));
	cdf_info$taxa=gsub("X","",rownames(cdfs));
	cdf_info$cdf_val=cdfs;
	cdf_info$num_taxa=length(cdf_info$taxa);
	cdf_info$filtering=character();

	return(cdf_info);
}

################################################################################

optimize_lines=function(x, y){
# Given two vectors describing lines, it will remove the lines that are redundant.

	diff_f=abs(c(1,diff(y)))>0;
	diff_r=abs(c(rev(diff(rev(y))),1))>0;
	keep=diff_f | diff_r;

	ol=list();
	ol$x=x[keep];
	ol$y=y[keep];

	return(ol);
}

################################################################################

plot_cdf=function(cdf_info, log=FALSE, max_abund=NULL, max_ubiq=NULL, mark_ubiquity=0, mark_abundance=0, top_right_axis=FALSE){
# Plot the CDF

	num_taxa=length(cdf_info$taxa);
	cat("Num taxa to plot: ", num_taxa, "\n");
	if(num_taxa==0){
		cat("Have all taxa been filtered?  There's nothing left to plot.\n");	
		return();
	}

	# Compute drawing priority.  If it's gray, then draw it later.
	# The larger the number, the later it is drawn, so it is on top.
	not_gray=substr(cdf_info$colors,1,4)!="gray";
	is_gray=!not_gray;
	priority=rep(2,num_taxa);
	priority[is_gray]=1;
	draw_order=order(priority);

	if(log==TRUE){
		non_zero_min=sort(cdf_info$cutoffs)[2];
		log_min=log(non_zero_min,10);

		cat("Using ", non_zero_min, " as minimum non zero abundance.\n", sep="");

		if(is.null(max_abund)){
			max_abund=0;
		}
		if(is.null(max_ubiq)){
			max_ubiq=1;
		}

		# Compute space to give for taxa labels
		max_taxa_name_len=max(nchar(cdf_info$taxa));
		cat("Longest taxa name length: ", max_taxa_name_len, "\n");
		left_name_space=((max_abund-log_min)/4)*max_taxa_name_len/50;
		bottom_name_space=1*max_taxa_name_len/190;
		cat("Left Name Space: ", left_name_space, "\n");
		cat("Bottom Name Space: ", bottom_name_space, "\n");

		# Set up blank plot
		if(top_right_axis){
			par(mar=c(1,1,9,5));
		}else{
			par(mar=c(5,5,5,1));
		}
		plot(0,0, type="n", xlab="", ylab="", main=InputFileName,
			xlim=c(log_min-left_name_space,max_abund), ylim=c(-bottom_name_space,max_ubiq),
			yaxt="n", xaxt="n", 
			);

		# y-axis labels
		ubiquity_ticks=(seq(0,100,10))/100;
		ubiquity_labels=sprintf("%i", ubiquity_ticks*100);

		# x-axis labels
		abundance_ticks=(log_min:0);
		abundance_labels=sprintf("%i", abundance_ticks);

		if(top_right_axis){
			#X
			axis(side=3, at=abundance_ticks, labels=abundance_labels, tick=TRUE, outer=FALSE, las=1);
			mtext("Percent Ubiquity", side=4, line=2.75, at=.5);
			#Y
			axis(side=4, at=ubiquity_ticks, labels=ubiquity_labels, tick=TRUE, outer=FALSE, las=2);
			mtext(expression(Log[10](Abundance)), side=3, line=2.25, at=-2);
		}else{
			#X
			axis(side=1, at=abundance_ticks, labels=abundance_labels, tick=TRUE, outer=FALSE, las=1);
			mtext(expression(Log[10](Abundance)), side=1, line=2.5, at=-2);
			#Y
			axis(side=2, at=ubiquity_ticks, labels=ubiquity_labels, tick=TRUE, outer=FALSE, las=2);
			mtext("Percent Ubiquity", side=2, line=2.75, at=.5);
		}


		if(cdf_info$cutoffs[1]==0){
			cdf_info$cutoffs[1]=non_zero_min;
		}
		log_cutoffs=log(cdf_info$cutoffs,10);


		# Plot each taxa's line
		if(num_taxa>0){

			for(t in draw_order){
				if(cdf_info$cdf_val[t,2]>0){
					ol=optimize_lines(log_cutoffs, cdf_info$cdf_val[t,]);


					if(not_gray[t]){
						thick=2;
					}else{
						thick=1;
					}
			
					lines(ol$x, ol$y, col=cdf_info$colors[t], 
						lty=cdf_info$line_type[t], lwd=thick);
				}

			}

			idx=numeric(num_taxa);
			for(t in 1:num_taxa){
				# Determine where curve touches the x-axis
				ubiquities=cdf_info$cdf_val[t,];
				min_nonzero=min(ubiquities[ubiquities>0]);
				idx[t]=max(which(min_nonzero==ubiquities));
			}

			# Label below x axis
			x=log_cutoffs[idx]+.08;
			labels=sprintf("%s  ", cdf_info$taxa);
			colors=cdf_info$colors;

			if(all(!not_gray)){
				selected_labels=!not_gray;
			}else{
				selected_labels=not_gray;
			}
			#print(labels[selected_labels]);
			text(x[selected_labels], bottom_name_space*-.02, 
				labels[selected_labels], cex=.3, col=colors[selected_labels], srt=90, pos=2);
		
			# Label along y axis on the left
			y=cdf_info$cdf_val[,2];
			labels=cdf_info$taxa;
			colors=cdf_info$colors;
			text(log_min, y[selected_labels], 
				labels[selected_labels], cex=.3, col=colors[selected_labels], pos=2);

		}
		
		# Draw the point where the cutoffs are
		if(!IgnoreBullsEye){
			log_mark_abundance=log(mark_abundance,10)
			points(log_mark_abundance, mark_ubiquity, pch=1, cex=1.5, col="red")
			points(log_mark_abundance, mark_ubiquity, pch=16, cex=.75, col="red")
		}

	}else{
		# Graph non-log based plot
	
		if(is.null(max_abund)){
			max_abund=1;
		}
		if(is.null(max_ubiq)){
			max_ubiq=1;
		}

		# Set up blank plot
		plot(0,0, type="n", main=InputFileName,
			xlab="Abundance", ylab="Ubiquity (%)", xlim=c(0,max_abund), ylim=c(0,max_ubiq), yaxt="n");

		# y-axis labels
		ubiquity_ticks=(seq(0,100,10))/100;
		ubiquity_labels=sprintf("%i", ubiquity_ticks*100);
		axis(side=2, at=ubiquity_ticks, labels=ubiquity_labels, tick=TRUE, outer=FALSE, las=2);

		# Plot each taxa's line
		if(num_taxa>0){
			for(t in draw_order){
				if(cdf_info$cdf_val[t,2]>0){
					ol=optimize_lines(cdf_info$cutoffs, cdf_info$cdf_val[t,]);
					lines(ol$x, ol$y,
						col=cdf_info$colors[t], lty=cdf_info$line_type[t]);
				}
			}
		}

		# Draw the point where the cutoffs are
		if(!IgnoreBullsEye){
			points(mark_abundance, mark_ubiquity, pch=1, cex=1.5, col="red")
			points(mark_abundance, mark_ubiquity, pch=16, cex=.75,  col="red")
		}
	}
}

################################################################################

load_color_map=function(colormapfn){
# Load a color mapping from a file

	color_map=list();
	map=as.matrix(read.table(colormapfn, sep="\t", header=FALSE, comment.char="", quote=""));
	num_entries=nrow(map);
	for(i in 1:num_entries){
		clean=gsub("^ *", "", map[i,1]);
		clean=gsub(" *$", "", clean);
		color_map[[clean]]=map[i,2];
	}
	return(color_map);
}

################################################################################

load_id_map=function(idmapfn){
	id_map=list();
	map=as.matrix(read.table(idmapfn, sep="\t", header=FALSE, comment.char="", quote=""));
	num_entries=nrow(map);
	cat("Number of entries loaded from ", idmapfn, " : ", num_entries, "\n", sep="");

	map_hash=list();
	for(i in 1:num_entries){
		map_hash[[map[i,1]]]=map[i,2];
	}
	return(map_hash);
}

################################################################################

taxa_to_color=function(taxa_list, color_map){
# Assign a color to each taxa, in the same order

	num_taxa=length(taxa_list);
	colors=rep("", num_taxa);
	if(num_taxa>0){
		for(i in 1:num_taxa){
			color=color_map[[taxa_list[i]]];
			if(!is.null(color)){
				colors[i]=color;
			}	
		}
	}
	return(colors);	
}

################################################################################

filter_cdfs_byKeeplist=function(cdf_info, keeplist_vector){
# Removes taxa not in the keeplist
	num_taxa=cdf_info$num_taxa;
	keep_index=rep(FALSE, num_taxa);
	for(i in 1:num_taxa){
		if(any(cdf_info$taxa[i]==keeplist_vector)){
			keep_index[i]=TRUE;
		}
	}
	cdf_info$cdf_val=cdf_info$cdf_val[keep_index,];
	cdf_info$taxa=cdf_info$taxa[keep_index];
	cdf_info$num_taxa=sum(keep_index);
	cdf_info$filtering=c(cdf_info$filtering, "KeepList");

	return(cdf_info);
}

################################################################################

filter_cdfs=function(cdf_info, ubiquity, abundance){
# Removes taxa that do not meet ubiquity/abundance cutoff

	# Find index into the abundances
	cat("Filtering for Abundance: ", abundance, "\n");
	abund_idx=min(which(abundance<=cdf_info$cutoffs));
	cat("Index of abundance cutoff: ", abund_idx, "\n");

	# Find taxa that exceed ubiquity
	cat("Filtering for Ubiquity: ", ubiquity, "\n");
	taxa_exc_ubiquity=which(cdf_info$cdf_val[,abund_idx]>=ubiquity);
	num_kept=length(taxa_exc_ubiquity);
	cat("Num taxa kept: ", num_kept, "\n");
	
	# Store results back into cdf_info
	cdf_info$cdf_val=cdf_info$cdf_val[taxa_exc_ubiquity,];
	cdf_info$taxa=cdf_info$taxa[taxa_exc_ubiquity];
	cdf_info$num_taxa=num_kept;
	cdf_info$filtering=c(cdf_info$filtering, "Core");

	if(cdf_info$num_taxa==1){
		# Annoying R feature that automatically changes single row matrices to vectors.
		cdf_info$cdf_val=matrix(cdf_info$cdf_val, nrow=1, ncol=length(cdf_info$cutoffs));
	}
	
	return(cdf_info);
}

################################################################################

remap_cdf_info=function(cdf_info, id_map){

	orig_names=cdf_info$taxa;
	new_names=character(length(orig_names));
	num_names=length(orig_names);
	for(i in 1:num_names){
		splits=strsplit(orig_names[i], ";")[[1]];
		for(j in 1:length(splits)){
			old_id=splits[j];
			new_id=id_map[[old_id]];
			if(is.null(new_id)){
				new_id=old_id;
			}
			splits[j]=new_id;
		}
		new_names[i]=paste(splits, collapse=";");	
	}
	cdf_info$taxa=new_names;
	return(cdf_info);
}

################################################################################

# Read in color map
color_map=NULL;
if(!is.null(ColorMap)){
	color_map=load_color_map(ColorMap);
}

# Read in id map
id_map=NULL;
if(!is.null(IDMap)){
	id_map=load_id_map(IDMap);
}

# Read in CDF info
cdf_info=read_cdf_from_file(InputFileName);

# Filter by list
if(KeepList!=""){
	keeplist_vector=read_keeplist_from_file(KeepList);
	cat("Keep List Members: \n");
	print(keeplist_vector);
	cdf_info=filter_cdfs_byKeeplist(cdf_info, keeplist_vector);
}

if(!is.null(id_map)){
	cdf_info=remap_cdf_info(cdf_info, id_map);
}

if(!GreySubcutoffTaxa){
	# Filter by cutoffs
	cdf_info=filter_cdfs(cdf_info, Ubiquity, Abundance);

	# Assign colors based on order in taxa of cdf_info
	colors=taxa_to_color(cdf_info$taxa, color_map);

	# Assign grey scale and line types
	line_types=rep(1, cdf_info$num_taxa);
	oc=0;
	if(cdf_info$num_taxa>0){
		for(i in 1:cdf_info$num_taxa){
			if(colors[i]==""){
				grey_idx=(oc%%3)+1;
				line_types[i]=(oc%%5)+1;
				if(grey_idx==1){
					colors[i]="gray0";
				}else if(grey_idx==2){
					colors[i]="gray30";
				}else if(grey_idx==3){
					colors[i]="gray60";
				}
				oc=oc+1;
			}
		}
	}
	cdf_info$colors=colors;
	cdf_info$line_type=line_types;

}else{
	# Do this if you want the taxa below cutoff in grey

	# Get colors for all taxa
	colors=taxa_to_color(cdf_info$taxa, color_map);

	for(i in 1:cdf_info$num_taxa){
		if(colors[i]==""){
			colors[i]="gray";
		}
	}

	# Get a list of taxa that exceed the cutoff
	above_cutoff_cdf_info=filter_cdfs(cdf_info, Ubiquity, Abundance);
	above_cutoff_taxa=above_cutoff_cdf_info$taxa;

	# If the taxa exceeds the cutoff, leave it alone, remove the name and color it grey
	for(i in 1:cdf_info$num_taxa){
		if(!any(cdf_info$taxa[i]==above_cutoff_taxa)){
			colors[i]="gray";
			cdf_info$taxa[i]="";
		}
	}
	
	cdf_info$colors=colors;
	cdf_info$line_type=rep(1, cdf_info$num_taxa);
}

ub_ext=sprintf(".u%0.2f", Ubiquity);
ab_ext=sprintf(".a%0.5f", Abundance);

################################################################################

pdf(paste(OutputFileRoot, ub_ext, ab_ext, ".pdf", sep=""), height=PaperHeight, width=PaperWidth);

#plot_cdf(cdf_info, log=FALSE, mark_ubiquity=Ubiquity, mark_abundance=Abundance);
plot_cdf(cdf_info, log=TRUE, mark_ubiquity=Ubiquity, mark_abundance=Abundance, top_right_axis=TRUE);
#plot_cdf(cdf_info, log=TRUE, mark_ubiquity=Ubiquity, mark_abundance=Abundance, top_right_axis=FALSE);
#plot_cdf(cdf_info, log=TRUE, max_ubiq=.6, mark_ubiquity=Ubiquity, mark_abundance=Abundance);
#plot_cdf(cdf_info, log=TRUE, max_ubiq=.25, mark_ubiquity=Ubiquity, mark_abundance=Abundance);

dev.off();

################################################################################

if(Ubiquity>0 || Abundance>0){
	cdf_info=filter_cdfs(cdf_info, Ubiquity, Abundance);

	fh=file(paste(OutputFileRoot, ub_ext, ab_ext, ".csv", sep=""), "w");
	abund_idx=min(which(Abundance<=cdf_info$cutoffs));
	order=order(cdf_info$cdf_val[,abund_idx], decreasing=TRUE);

	num_exc_cutoff=length(order);

	cat(file=fh, "Filename:,", OutputFileRoot, "\n", sep="");
	cat(file=fh, "Ubiquity:,", Ubiquity, "\n", sep="");
	cat(file=fh, "Abundance:,", Abundance, "\n", sep="");
	cat(file=fh, "NumTaxa Exceeding Cutoffs:,", num_exc_cutoff, "\n", sep="");
	cat(file=fh, "\n");

	cat(file=fh, "Index,TaxaName,Ubiquity at ",Abundance," Abundance\n", sep="");
	count=1;
	for(i in order){
		cat(file=fh, count, ",", cdf_info$taxa[i], ",", cdf_info$cdf_val[i,abund_idx], "\n", sep="");
		count=count+1;
	}
	close(fh);
}

################################################################################

w=warnings();
if(length(w)){
	print(w);
}
cat("Done.\n")

q(status=0)
