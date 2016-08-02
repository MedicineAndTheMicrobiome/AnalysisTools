#!/usr/bin/env Rscript

###############################################################################
PaperHeight=9.5;
PaperWidth=20;
options(width=200);

library('getopt');

params=c(
	"input_file_list", "i", 1, "character",
	"is_simple", "I", 2, "logical",
	"subsample", "b", 2, "numeric",
	"percent_cutoff", "p", 2, "numeric",
	"output_file_root", "o", 2, "character",
	"regions_file", "r", 2, "character",
	"height", "h", 2, "numeric",
	"width", "w", 2, "numeric",
	"sets", "s", 2, "character",
	"test", "t", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n",
	"\n",
	"	-i <Input File List>\n",
	"	[-I (simple file list tag)]\n",
	"	[-b <subsample, eg. 10000, default is none performed>]\n",
	"	[-p <minimum percent identity cutoff, eg. 95, default=none>]\n",
	"\n",
	"	[-o <Output Filename Root>]\n",
	"	[-r <Regions File>]\n",
	"	[-h <height, default=", PaperHeight, " in>]\n",
	"	[-w <width, default=", PaperWidth, " in>]\n",
	"	[-s <Sets File>]\n",
	"\n",
	"This script will generate a fragment recruitment plot based on the\n",
	"work of Doug Rusch.\n",
	"\n",
	"The input file list is a list of fragment recruitment positions:\n",
	"\n",
	"	<display name>\\t<path of recruitment file>\\t<R color>\\n\n",
	"\n",
	"The -I flag allows you to input a very simple input file list.\n",
	"	If you just pass in a list of file paths, it will automatically\n",
	"	assume the file name is the sample ID.\n",
	"\n",
	"The format of the recruitment file is:\n",
	"	<begin>\\t<end>\\t<percent identity>\\n\n",
	"\n",
	"Where the begin and end are the coordinates of a read aligned to the referenece\n",
	"genome.  These can be cut out of the -m 8 formatted blast results.  The percent\n",
	"identity can also be used directly from blast, as it is a number between 0 and 100.\n",
	"It's not a proportion or probability, it's a percentage.\n",
	"\n",
	"Generally speaking, using 0 or 1 space or residue coordinates won't\n",
	"make a difference since we are so zoomed out.\n",
	"\n",
	"The sets file will contain a list of sample IDs and whether or not\n",
	"to include the datasets.\n",
	"\n",
	"\n",
	"\n",
	"\n", sep="");

if(!length(opt$input_file_list)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;

if(length(opt$is_simple)){
	IsSimple=T;
}else{
	IsSimple=F;
}

if(length(opt$subsample)){
	Subsample=opt$subsample;
}else{
	Subsample=-1;
}

if(length(opt$percent_cutoff)){
	PercentIDCutoff=opt$percent_cutoff;
}else{
	PercentIDCutoff=0;
}

if(length(opt$output_file_root)){
	OutputFileRoot=opt$output_file_root;
}else{
	OutputFileRoot=gsub(".txt$", "", opt$input_file);
	OutputFileRoot=gsub(".tsv$", "", OutputFileRoot);
}

if(length(opt$height)){
	PaperHeight=opt$height;
}

if(length(opt$width)){
	PaperWidth=opt$width;
}


if(length(opt$regions_file)){
	RegionsFile=opt$regions_file;
}else{
	RegionsFile=NULL;
}

if(length(opt$sets)){
	SetsFile=opt$sets;
}else{
	SetsFile=NULL;
}

cat("Input File List: ", InputFileName, "\n");
cat("Output File Root: ", OutputFileRoot, "\n");

if(!is.null(RegionsFile)){
	cat("Regions File: ", RegionsFile, "\n");
}


################################################################################

guess_sample_ids=function(filelist){

		cat("Guessing sample ids...\n");

		# Extract filename from path
		file_name=character();
		num_samp=nrow(filelist);
		for(i in 1:num_samp){
			file_name[i]=tail(strsplit(filelist[i,], "/")[[1]], 1);
		}

		# If there is only one file, let it be the sample id
		if(length(filelist)==1){
			return(file_name);
		}


		# Split each file name into components by periods
		samp_ids_list=list();
		max_len=100;
		comp=matrix("", nrow=num_samp, ncol=max_len);
		for(i in 1:length(file_name)){
			samp_ids_list[[i]]=strsplit(file_name[i], "\\.")[[1]];			
			num_comp=length(samp_ids_list[[i]]);
			if(max_len<num_comp){
				cat("Error:  Too many periods in file name.\n");
				quit(status=-1);
			}

			# Save the reverse 
			comp[i, 1:num_comp]=rev(samp_ids_list[[i]]);
		}

		# If all the components are the same, remove it.
		for(i in 1:ncol(comp)){
			if(length(unique(comp[,i]))==1){
				comp[,i]="";
			}
		}


		# Reconstruct a sample id based on the kept components
		samp_ids=numeric();
		for(i in 1:num_samp){
			samp_comp=comp[i, comp[i,]!=""];
			samp_ids[i]=paste(rev(samp_comp), collapse=".");
		}

		cat("Guessed sample ids:\n");
		print(samp_ids);
		
		return(samp_ids);

}

read_file_list=function(filelist_name, is_simple){

	if(!is_simple){
		file_info=as.matrix(read.delim(filelist_name, sep="\t", header=F, row.names=1, comment.char="#"));
		if(ncol(file_info)==1){
			file_info=cbind(file_info,rep(NA, nrow(file_info)));
		}
		colnames(file_info)=c("Path", "Color");
		#print(file_info);
	}else{
		# Read in list and guess what the sample ID should be
		cat("Loading read file list as simple...\n");
		file_info=as.matrix(read.delim(filelist_name, sep="\t", header=F, comment.char="#"));
		print(file_info);

		# Extract file name
		samp_ids=guess_sample_ids(file_info);

		file_info=cbind(file_info, c(1:length(samp_ids)));
		colnames(file_info)=c("Path", "Color");
		rownames(file_info)=samp_ids;

	}
	return(file_info);
}

read_hits=function(filename){
	cat("Loading: ", filename, "\n", sep="");
	file_info=as.matrix(read.delim(filename, sep="\t", header=F));
	
	# Removed emptylines
	isna=is.na(file_info[,3]);
	file_info=file_info[!isna,, drop=F];

	if(ncol(file_info)!=3){
		cat("Error:  Input hits (", filename, ") do not have 3 columns.\n", sep="");
		quit(status=-1);
	}

	colnames(file_info)=c("Begin", "End", "PercID");
	return(file_info);
}

read_regions=function(filename){
	cat("Loading: ", filename, "\n", sep="");
	file_info=as.matrix(read.delim(filename, sep="\t", header=F, row.names=1));
	colnames(file_info)=c("Begin", "End");
	return(file_info);
}

read_sets_file=function(filename){

	# Load entire table in
	sets_table=read.table(filename, sep="\t", header=T, check.names=F,
		comment.char="", quote="", row.names=1);

	for(i in 1:ncol(sets_table)){
		sets_table[,i]=ifelse(sets_table[,i]=="", NA, as.character(sets_table[,i]));
	}

	return(sets_table);
}

################################################################################

file_info=read_file_list(InputFileName, IsSimple);

num_datasets=nrow(file_info);
hits=list();
max_coord=0;
for(i in 1:num_datasets){

	# Load hits
	hits_info=read_hits(file_info[i, "Path"]);

	# Try to identify size of genome
	max_coord=max(c(max_coord, 
		max(
			c(hits_info[,1], hits_info[,2])
		)));

	sample_name=rownames(file_info)[i];

	cat("Num Positions Read: ", nrow(hits_info), "\n");

	# Extract only hits above cutoff
	if(PercentIDCutoff>0){
		hits_info=hits_info[hits_info[,3]>=PercentIDCutoff,];	
		cat("Num Positions >=", PercentIDCutoff, ": ", nrow(hits_info), "\n", sep="");
	}

	# Perform subsampling if requested
	if(Subsample!=-1){
		num_pos=nrow(hits_info);
		targ_ss=min(c(Subsample, num_pos));
		cat("Keeping subsample of size: ", targ_ss, " of ", num_pos, "\n", sep="");
		ss_ix=sample(1:nrow(hits_info), targ_ss, replace=F);
		hits_info=hits_info[ss_ix,, drop=F];
	}

	hits[[sample_name]]=hits_info;
}

cat("Largest coordinate position: ", max_coord, "\n");

################################################################################

if(!is.null(SetsFile)){
	sets_table=read_sets_file(SetsFile);
	cat("\nSets table: \n");
	print(sets_table);
	cat("\n");
}else{
	sets_table=matrix(1:num_datasets, nrow=num_datasets, ncol=1);
	rownames(sets_table)=rownames(file_info);
}

################################################################################

if(!is.null(RegionsFile)){
	regions_info=read_regions(RegionsFile);
	print(regions_info);
}else{
	regions_info=NULL;
}

################################################################################

log_trans=function(x, lb){
	x[x==0] = lb/10;
	trans=log(x,10);
	return(trans);
}

pick_colors=function(num_colors){
	if(num_colors==1){
		num_colors=2;
	}
	if(num_colors<=12){
		colors=c("blue", "red", "green",
			"orange", "purple", "brown",
			"#7366BD", "pink", "#0D98BA",
			"#FFAE42", "#C0448F", "#C5E384"
			)
	}else{
		subcolors=ceiling(num_colors/4);
		colors=c(
			rainbow(subcolors, s=1, v=.9),
			rainbow(subcolors, s=1, v=.66),
			rainbow(subcolors, s=.66, v=.9),
			rainbow(subcolors, s=.66, v=.66)
		);
	}	
	return(colors[1:num_colors]);
#	barplot(rep(1,num_colors), col=colors);
}


plot_hits=function(file_info, hits, genome_size, min_perc_id, regions_info=NULL, plot_name=NULL){

	num_datasets=nrow(file_info);

	subsample_names=rownames(file_info);

	if(all(is.na(file_info[,"Color"]))){
		cat("No colors were assigned.  Using rainbow.\n");
		cat("Number of colors to assign: ", num_datasets, "\n");
		file_info[,"Color"]=pick_colors(num_datasets);
	}else{
		num_color_groups=length(unique(file_info[,"Color"]));
		palette(pick_colors(num_color_groups));
	}

	layout_matrix=matrix(c(
		2,2,2,2,2,2,2,2,2,2,4,
		2,2,2,2,2,2,2,2,2,2,4,
		1,1,1,1,1,1,1,1,1,1,3,
		1,1,1,1,1,1,1,1,1,1,3,
		1,1,1,1,1,1,1,1,1,1,3,
		1,1,1,1,1,1,1,1,1,1,3,
		1,1,1,1,1,1,1,1,1,1,3,
		1,1,1,1,1,1,1,1,1,1,3,
		5,5,5,5,5,5,5,5,5,5,5,
		5,5,5,5,5,5,5,5,5,5,5
	), ncol=11, byrow=T);
	layout(layout_matrix);	
	ml=4.1;
	mr=0;
	mt=4.1;
	mb=5.1;
	perc_id_range=c(min_perc_id,100);
	
	#------------------------------------------------------
	# Plot 1 Fragment recruitment
	par(mar=c(4,ml,0,mr));
	plot(0,0, 
		main="",
		ylab="Percent Identity (%)",
		xlab="Genomic Position (bp)",
		bty="o",
		xaxt="n",
		yaxt="n",
		xlim=c(0, genome_size), ylim=perc_id_range);


	# Keep track of regions to grey out
	non_regions=matrix(0, nrow=2, ncol=2);
	non_regions[1,]=c(-genome_size/10, 0);
	non_regions[2,]=c(genome_size, genome_size+genome_size/10);

	cat("\n");
	if(!is.null(regions_info)){
		for(i in 2:(nrow(regions_info))){
			non_regions=rbind(non_regions,
					c(regions_info[i-1,2], regions_info[i, 1]));
		}

		num_ticks=genome_size/10;
		region_names=rownames(regions_info);
		for(i in 1:nrow(regions_info)){

			cat("Region Name: ", region_names[i], "\n", sep="");

			region_width=(regions_info[i, 2] - regions_info[i, 1]);
			cat("  Region Width: ", region_width, "\n");

			ticks_in_region = ceiling((region_width/genome_size)*10);
			cat("  Ticks in Region: ", ticks_in_region, "\n");

			region_ticks=seq(0, region_width, signif(region_width/ticks_in_region, 1));
			global_ticks=regions_info[i, 1] + region_ticks;
			labels=sprintf("%g", region_ticks);
			axis(side=1, at=global_ticks, labels=labels);

			# Label the region name
			text(regions_info[i, 1], 80, pos=4, region_names[i]);
		}
	
	}else{
		ticks=seq(0, genome_size, signif(genome_size/10, 1));
		labels=sprintf("%g", ticks);
		axis(side=1, at=ticks, labels=labels);
	}
	cat("\n");

	# Grey out non-regions
	for(i in 1:nrow(non_regions)){
		cat("Greying from: ", non_regions[i,1], " - ", non_regions[i,2], "\n", sep="");
		rect(non_regions[i,1], 0 , non_regions[i,2], 110, col="grey85", border=NA);
	}

	# Plot y axis
	yticks=seq(0, 100, 5);
	axis(side=2, at=yticks, labels=sprintf("%3.0f", yticks));

	comb_begins=numeric(0);	
	comb_ends=numeric(0);	
	comb_perc_id=numeric(0);	
	comb_color=character(0);

	for(i in 1:num_datasets){
		ssn=subsample_names[i]
		comb_begins=c(comb_begins, hits[[ssn]][,1]);
		comb_ends=c(comb_ends, hits[[ssn]][,2]);
		comb_perc_id=c(comb_perc_id, hits[[ssn]][,3]);
		comb_color=c(comb_color, rep(file_info[i, "Color"], length(hits[[ssn]][,1])));
	}

	num_total_points=length(comb_begins);
	cat("Total hits: ", num_total_points, "\n");

	rand_seq=sample(1:num_total_points, num_total_points);

	num_heartbeats=10000;
	iter=0;
	cat(paste(rep(".", num_total_points/num_heartbeats), collapse=""), "\n");
	for(hix in rand_seq){	
		lines(
			c(comb_begins[hix], comb_ends[hix]), 
			c(comb_perc_id[hix], comb_perc_id[hix]), 
			col=comb_color[hix], lwd=3,
			lend="butt"
		);	
		iter=iter+1;
		if(!(iter%%num_heartbeats)){
			cat(">");
		}
	}
	cat("\n");

	#------------------------------------------------------
	# Plot 2 Genomic Density
DENSITY=F;

	if(DENSITY==T){
		hist_rec=list();
		density=list();
		min_den=1;
		for(i in 1:num_datasets){
			ssn=subsample_names[i]
			hist_rec[[i]]=hist(hits[[ssn]][,1], breaks=100, plot=F);
			density[[i]]=hist_rec[[i]]$counts/sum(hist_rec[[i]]$counts);
			min_den=min(min_den, density[[i]][density[[i]]!=0]);
		}
		cat("Genomic min non zero density: ", min_den, "\n");

		for(i in 1:num_datasets){
			d=density[[i]];
			density[[i]]=log_trans(d, min_den);
		}
		min_den=log(min_den,10);
			
		par(mar=c(0,ml,mt,mr));
		plot(0,0, xlim=c(0,genome_size), ylim=c(min_den,0),  
			main=plot_name,
			bty="n",
			type="n", 
			xaxt="n",
			#yaxt="n",
			ylab="Log(Prop)",
			xlab="Genomic Position (bp)"
			);

		for(i in 1:num_datasets){
			points(hist_rec[[i]]$mids, density[[i]], col=file_info[i, "Color"], type="l");
		}
	}else{
		cdfs=list();
		hist_rec=list();
		divisions=seq(0, genome_size, length.out=200);
		for(i in 1:num_datasets){
			ssn=subsample_names[i]
			hist_rec[[i]]=hist(hits[[ssn]][,1], breaks=divisions, plot=F);
			total=sum(hist_rec[[i]]$counts);
			cdfs[[i]]=cumsum(c(0, hist_rec[[i]]$counts))/total;
		}

		par(mar=c(0, ml, mt, mr));
		plot(0,0, xlim=c(0,genome_size), ylim=c(0,1),  
			main=plot_name,
			bty="n",
			type="n", 
			xaxt="n",
			#yaxt="n",
			ylab="CDF",
			xlab="Genomic Position (bp)"
			);

		for(i in 1:num_datasets){
			points(divisions, cdfs[[i]], type="l", col=file_info[i, "Color"]);
		}

	}

	#------------------------------------------------------
	# Plot 3 Identity Density


	if(DENSITY==T){
		hist_rec=list();
		density=list();
		min_den=1;
		for(i in 1:num_datasets){
			ssn=subsample_names[i];
			hist_rec[[i]]=hist(hits[[ssn]][,3], breaks=25, plot=F);
			density[[i]]=hist_rec[[i]]$counts/sum(hist_rec[[i]]$counts);
			min_den=min(min_den, density[[i]][density[[i]]!=0]);
		}
		cat("Identity min non zero density: ", min_den, "\n");

		for(i in 1:num_datasets){
			d=density[[i]];
			density[[i]]=log_trans(d, min_den);
		}
		min_den=log(min_den,10);

		par(mar=c(4,0,0,1));
		plot(0,0, type="n", 
			xlab="Log(Prop)", ylab="", 
			#xaxt="n",
			yaxt="n",
			bty="n",
			xlim=c(min_den,0),
			ylim=perc_id_range);

		for(i in 1:num_datasets){
			points(density[[i]], hist_rec[[i]]$mids, col=file_info[i, "Color"], type="l");
		}
	}else{
		cdfs=list();
		hist_rec=list();
		divisions=seq(min(perc_id_range), max(perc_id_range), length.out=100);
		for(i in 1:num_datasets){
			ssn=subsample_names[i]
			hist_rec[[i]]=hist(hits[[ssn]][,3], breaks=divisions, plot=F);
			total=sum(hist_rec[[i]]$counts);
			cdfs[[i]]=cumsum(c(0, hist_rec[[i]]$counts))/total;
		}

		par(mar=c(4, 0, 0, 1));
		plot(0,0, 
			xlim=c(0,1), 
			ylim=perc_id_range,  
			bty="n",
			type="n", 
			xaxt="n",
			ylab="CDF",
			xlab=""
			);

		for(i in 1:num_datasets){
			points(cdfs[[i]], divisions, type="l", col=file_info[i, "Color"]);
		}
	}

	#------------------------------------------------------
	# Plot 4 Empty Space
	par(mar=c(0,0,mt,mr));
	plot(0,0,type="n", bty="n", yaxt="n", xaxt="n", xlab="", ylab="");

	if(DENSITY==T){
		text(0,0,label="Probability", srt=-45);
	}else{
		text(0,0,label="Cumulative", srt=-45);
	}

	#------------------------------------------------------
	# Plot 5 Labels
	
	par(mar=c(0,ml,0,mr));
	plot(0,0, xlim=c(0,50), ylim=c(-1,0), type="n",
		xlab="", ylab="",
		bty="n", 
		yaxt="n", 
		xaxt="n"
		);
	for(i in 1:num_datasets){
		text(i-1,0,
			labels=(file_info[i,"Name"]), 
			col=file_info[i, "Color"], 
			pos=4,
			srt=-45
		);
	}
	
}

################################################################################

combine_hits=function(subset_hits, sample_names){
	combined_hits=numeric();
	for(smp_nm in sample_names){
		combined_hits=rbind(combined_hits, subset_hits[[smp_nm]]);
	}
	cat("Num data points: ", nrow(combined_hits), "\n");
	return(combined_hits);
}

get_ks=function(x, y){
	res=ks.test(x,y);
	pval=res$p.value;
	str=sprintf("Kolmogorov-Smirnov Test: p-value = %g", pval);
	return(str);
}

plot_qqplots=function(subset_file_info, subset_hits, plot_set_name){

	cat("\n*** Starting QQ Plots: \"", plot_set_name, "\"\n", sep="");

	# Get sample to group info
	group_names=unique(subset_file_info[,"Group"]);
	sample_names=rownames(subset_file_info);
	
	samp_names_a=sample_names[subset_file_info[,"Group"]==group_names[1]];
	samp_names_b=sample_names[subset_file_info[,"Group"]==group_names[2]];

	# Combine the hits across samples in the same group
	cat("\n");
	cat("Sample A (", group_names[1], ") Members: \n");
	print(samp_names_a);
	sample_a_hits=combine_hits(subset_hits, samp_names_a);	

	cat("\n");
	cat("Sample B (", group_names[2], ") Members:\n");
	print(samp_names_b);
	sample_b_hits=combine_hits(subset_hits, samp_names_b);	

	# Compute the K-S statistics
	ks_res_str_pos=get_ks(sample_a_hits[,"Begin"], sample_b_hits[,"Begin"]);
	ks_res_str_perc=get_ks(sample_a_hits[,"PercID"], sample_b_hits[,"PercID"]);

	par(mfrow=c(1,2));
	par(mar=c(6,6,8,4));
	
	# Plot positional QQ Plot
	qqplot(sample_a_hits[,"Begin"], sample_b_hits[,"Begin"], 
		main="Positional Read Distribution Q-Q Plot",
		xlab=paste(group_names[1], " bps", sep=""),
		ylab=paste(group_names[2], " bps", sep="")
	);
	abline(coef=c(0,1), col="blue");
	mtext(paste(group_names[2], " vs. ", group_names[1]), side=3, line=2.5);
	mtext(ks_res_str_pos, side=3, line=1.5);

	# Plot perc identity QQ Plot
	qqplot(sample_a_hits[,"PercID"], sample_b_hits[,"PercID"],
		main="Percent Identity Read Distribution Q-Q Plot",
		xlab=paste(group_names[1], " % Identity", sep=""),
		ylab=paste(group_names[2], " % Identity", sep="")
	);
	abline(coef=c(0,1), col="blue");
	mtext(paste(group_names[2], " vs. ", group_names[1]), side=3, line=2.5);
	mtext(ks_res_str_perc, side=3, line=1.5);

	cat("\n");	
	cat("*** Ending QQ Plots\n");
	return;
}

################################################################################

plot_stats=function(stat, colors, groups){
	
	# Order samples by statistic
	order_ix=order(stat[,1]);
	stat=stat[order_ix,, drop=F];	

	# Get range of stat and reorder group colors
	stat_range=range(stat[,1]);
	stat_name=colnames(stat);	
	samp_names=rownames(stat);
	groups=groups[samp_names, , drop=F];
	colors=colors[samp_names, , drop=F];
	grp_col_map=list();
	for(i in 1:length(groups)){
		grp_col_map[[groups[i,1,drop=F]]]=colors[i,1,drop=F];
	}

	# Prepare to adjust positions so they don't collide in the plot
	num_samples=nrow(stat);
	stat_range=range(stat);
	min_spacing=(stat_range[2]-stat_range[1])/50;	
	if(min_spacing==0){
		min_spacing=1;
	}
	adj=min_spacing/20;
	cat("Min spacing: ", min_spacing, "\n");
	cat("Spacing adjustment increment: ", adj, "\n");

	# Repeatedly adjust the positions until they don't collide
	adjusted=F;
	adj_pos=stat;
	for(adj_ix in 1:10000){
                adjusted=F
                for(forw in 1:(num_samples-1)){
                        if(abs(adj_pos[forw]-adj_pos[forw+1])<min_spacing){
                                adj_pos[forw]=adj_pos[forw]-adj;
                                adjusted=T;
                                break;
                        }
                }
                for(rev in (num_samples:2)){
                        if(abs(adj_pos[rev]-adj_pos[rev-1])<min_spacing){
                                adj_pos[rev]=adj_pos[rev]+adj;
                                adjusted=T;
                                break;
                        }
                }
                if(!adjusted){
                        cat("Done adjusting at iteration: ", adj, "\n");
                        break;
                }
        }

	# Create an empty plot based on the ranges of the adjusted positions
	pos_range=range(adj_pos);
	plot(0,0, type="n", ylim=c(pos_range[1], pos_range[2]), xlim=c(0,10),
		main=stat_name, xlab="", ylab=stat_name, xaxt="n"
	);

	# Draw hash line which marks the actual value
	segments(.5, stat[,1], 1, stat[,1], lwd=2);
	# Draw connecting line to the label that was adjusted
	segments(1, stat[,1], 2, adj_pos[,1]);
	# Draw label using colors that were specified
	text(2, adj_pos, labels=samp_names, pos=4, col=colors[,1]);

	# Find group medians
	group_asgn=groups[,1];
	names(group_asgn)=rownames(groups);
	unique_groups=sort(unique(group_asgn));
	num_groups=length(unique_groups);
	cat("Num Groups: ", num_groups, "\n");
	group_medians=numeric(num_groups);
	names(group_medians)=unique_groups;
	group_values=list();
	for(grp in unique_groups){
		cat("Group: ", grp, "\n");
		members=names(group_asgn[group_asgn==grp]);
		print(stat[members,1]);
		group_medians[grp]=median(stat[members,1]);
		group_values[[grp]]=stat[members,1];
	}

	# Plot group medians
	for(grp in unique_groups){
		points(0, group_medians[grp], pch=19, cex=3, col=grp_col_map[[grp]]);
		text(7, group_medians[grp], label=sprintf("%3.4f", group_medians[grp]), cex=2, pos=4, col=grp_col_map[[grp]]);
	}

	ret_rs=rep(NA, 8);
	names(ret_rs)=c(
		"GroupA",
		"MedianA",
		"GroupB",
		"MedianB",
		"Difference",
		"nA",
		"nB",
		"p-value"
	);

	# Compute wilcox if there are only two groups
	if(length(unique_groups)==2){
		res=wilcox.test(group_values[[unique_groups[1]]], group_values[[unique_groups[2]]]);
		mtext(sprintf("Median: %3.4f", group_medians[1]), line=2, col=grp_col_map[[unique_groups[1]]], cex=.6);
		mtext(sprintf("Median: %3.4f", group_medians[2]), line=1, col=grp_col_map[[unique_groups[2]]], cex=.6);
		mtext(sprintf("Wilcoxon Rank Sum Test, p-value: %g", res$p.value), line=0, cex=.6);

		ret_rs[]=c(
			unique_groups[1],
			group_medians[1], 
			unique_groups[2],
			group_medians[2], 
			abs(group_medians[1]-group_medians[2]),
			length(group_values[[unique_groups[1]]]),
			length(group_values[[unique_groups[2]]]),
			res$p.value
		);
	}
print(ret_rs);
	return(ret_rs);
}

################################################################################

co_str=sprintf(".%02g", PercentIDCutoff);

pdf(paste(OutputFileRoot, co_str, ".pdf", sep=""), height=PaperHeight, width=PaperWidth);

num_records=nrow(file_info);

min_perc_id=100;
for(i in 1:num_records){
	ranges=range(c(hits[[i]][,1], hits[[i]][,2]));
	cat(file_info[i, 1], ": ", ranges[1], " - ", ranges[2], " \n");

	min_perc_id=min(min_perc_id, hits[[i]][,3]);
}

cat("Minimum percent identify found: ", min_perc_id, "\n");
min_disp_perc_id=(min_perc_id %/% 5)*5;
cat("Minimum displayed percent identity: ", min_disp_perc_id, "\n");


cat("\n");
cat("Samples loaded with hits:\n");
samples_with_hits=names(hits);
print(samples_with_hits);
cat("\n");

cat("Samples with set information:\n");
samples_with_setinfo=rownames(sets_table);
print(samples_with_setinfo);
cat("\n");

cat("Shared samples:\n");
shared_samples=sort(intersect(samples_with_hits, samples_with_setinfo));
print(shared_samples);
cat("\n");

cat("Subsetting shared samples\n");
file_info=file_info[shared_samples,, drop=F];
sets_table=sets_table[shared_samples,, drop=F];

print(file_info);
print(sets_table);

num_sets=ncol(sets_table);
sets_names=colnames(sets_table);
sample_names=rownames(sets_table);


# Compute statistics for position and percent identity.
marg_stats_of_interest=c(1,2,3,4);
cat("Precomputing sample specific statistics:\n");
num_samples=length(sample_names);
marginal_stats=matrix(0, nrow=num_samples, ncol=length(marg_stats_of_interest));
rownames(marginal_stats)=sample_names;
colnames(marginal_stats)=c("Percent Identity", "K-S Statistic", "Coverage (reads/bp)", "Coverage (bases/bp)");

perc_id_range=c(min_disp_perc_id, 100);
par(mfrow=c(2,4));

pgenome=function(q, lower.tail=TRUE, log.p=FALSE){
	return(punif(q, min=0, max=max_coord, lower.tail, log.p));
}

total_reads=0;
total_bases=0;
for(sample_id in sample_names){	
	samp_hits=hits[[sample_id]];

	med=median(samp_hits[,"PercID"]);
	total_reads=total_reads+nrow(samp_hits);
	hist(
		samp_hits[,"PercID"],
		breaks=seq(perc_id_range[1], perc_id_range[2], length.out=10),
		main=paste("Perc Id Hit Distribution: ", sample_id),
		xlab="Percent Identity (%)" 
	);
	
	# Compute KS
	ks_res=ks.test(samp_hits[,"Begin"], pgenome);

	# Estimate Coverage
	num_hits=length(samp_hits[,"PercID"]);
	num_bases=sum(abs(samp_hits[,"Begin"]-samp_hits[,"End"])+1);
	total_bases=total_bases+num_bases;
	read_depth=num_hits/max_coord;
	base_cov_depth=total_bases/max_coord;

	marginal_stats[sample_id, ]=c(med, ks_res$statistic, read_depth, base_cov_depth);

	mtext(paste("Median Perc ID: ", med), cex=.75, line=0);
	mtext(paste("KS Statistic: ", ks_res$statistic), cex=.75, line=-1);
	mtext(paste("KS p-value: ", ks_res$p.value), cex=.75, line=-2);
}

overall_read_coverage=total_reads/max_coord;
overall_base_coverage=total_bases/max_coord;
plot(0,0, type="n", bty="n", ylim=c(0,40), xlim=c(0,40), yaxt="n", xaxt="n", xlab="", ylab="");
text(0, 35, sprintf("Overall Reads/Genomic Base: %g", overall_read_coverage), pos=4);
text(0, 33, sprintf("Overall Bases/Genomic Base: %g", overall_base_coverage), pos=4);
text(0, 31, sprintf("Total Reads: %g", total_reads), pos=4);

# Variables to store statistical tests
cumulative_stat_res=matrix(NA, nrow=1, ncol=8);
cumulative_stat_name=matrix("", nrow=1, ncol=2);
colnames(cumulative_stat_name)=c("Set Name", "Stat Name");

# Generate plots for each set variation
for(i in 1:num_sets){

	cat("\n\n*********************************************************************\n\n");
	cat("Plotting: ", sets_names[i], "\n");

	categories=as.character(sets_table[,i]);

	group_info=cbind(sample_names, categories);

	included=which(!is.na(group_info[,2]));
	group_info=group_info[included, ];

	cat("Included: \n");
	print(included);
	
	# Build subset
	new_name=paste(sample_names[included], " (", categories[included], ")", sep ="");
	subset_file_info=cbind(new_name, file_info[included, "Path"], 
		as.factor(categories[included]), as.character(categories[included]));
	colnames(subset_file_info)=c("Name", "Path", "Color", "Group");

	plot_set_name=paste(OutputFileRoot, " : ", sets_names[i], sep="");

	plot_hits(subset_file_info, hits, max_coord, min_disp_perc_id, regions_info, 
		plot_name=plot_set_name
	);

	# Plot rankings of position and percent identity
	par(mfrow=c(1,length(marg_stats_of_interest)));
	par(mar=c(4,6,8,2));
	for(statix in marg_stats_of_interest){
		used_samples=rownames(subset_file_info);
		stat_res=plot_stats(
			marginal_stats[used_samples,statix, drop=F], 
			subset_file_info[,"Color", drop=F],
			subset_file_info[,"Group", drop=F]
		);

		cumulative_stat_res=rbind(cumulative_stat_res, stat_res);
		cumulative_stat_name=rbind(cumulative_stat_name, 
			c(sets_names[i], colnames(marginal_stats[used_samples,statix, drop=F])));
	}

	# Plot statistical tests
	groups=unique(categories[included]);
	print(groups);
	if(length(groups)==2){
		plot_qqplots(subset_file_info, hits, plot_set_name);
	}

}

cumulative_stat_res=cumulative_stat_res[-1,];
cumulative_stat_name=cumulative_stat_name[-1,];
#print(cumulative_stat_res);
#print(cumulative_stat_name);

dev.off();

################################################################################
# Output statistics in a table

table_fname=paste(OutputFileRoot, ".frg_rec_stats.tsv", sep="");
fh=file(table_fname, "w");

out_fields=c(colnames(cumulative_stat_name), colnames(cumulative_stat_res));
cat(file=fh, paste(out_fields, collapse="\t"), "\n", sep="");

num_stats=nrow(cumulative_stat_name);
for(i in 1:num_stats){
	out_fields=c(cumulative_stat_name[i,], cumulative_stat_res[i,]);
	cat(file=fh, paste(out_fields, collapse="\t"), "\n", sep="");
}

close(fh);

################################################################################

w=warnings();
if(length(w)){
	print(w);
}
cat("Done.\n")

q(status=0)
