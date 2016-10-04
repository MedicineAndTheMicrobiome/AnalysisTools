#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"offset_file", "t", 1, "character",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	-t <offset file>\n",
	"	[-o <output file root name>]\n",
	"\n",
	"	This script will read in the summary table\n",
	"	and a file describing the time from the first\n",
	"	sample.\n",
	"\n",
	"	For each sample a time series bar plot and\n",
	"	and diversity will be generated.\n",
	"\n",
	"	The format of the offset file is:\n",
	"\n",
	"	<sample id> \\t <sample (individual) grouping id> \\t <time stamp> \\t <sample (cohort) group id>\\n",
	"\n",
	"\n");

if(!length(opt$input_file) || !length(opt$offset_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OffsetFileName=opt$offset_file;

if(length(opt$output_file)>0){
	OutputFileRoot=opt$output_file;
}else{
	OutputFileRoot=InputFileName;
	OutputFileRoot=gsub("\\.summary_table\\.tsv$", "", OutputFileRoot);
	OutputFileRoot=gsub("\\.summary_table\\.xls$", "", OutputFileRoot);
	cat("No output file root specified.  Using input file name as root.\n");
}

###############################################################################

OutputPDF = paste(OutputFileRoot, ".div_ts.pdf", sep="");
cat("Output PDF file name: ", OutputPDF, "\n", sep="");
pdf(OutputPDF,width=8.5,height=11)

###############################################################################

load_offset=function(fname){
        cat("Loading Offsets: ", fname, "\n");
        offsets_mat=read.delim(fname,  header=FALSE, row.names=1, sep="\t", comment.char="#", quote="");
	colnames(offsets_mat)=c("Indiv ID", "Offsets", "Group ID");

	# reset offsets
	groups=unique(offsets_mat[,"Indiv ID"]);
	
	cat("Groups:\n");
	print(groups);
	cat("\n");

	# Reset offsets so they are relative to the first/smallest sample
	for(gid in groups){
		offsets=offsets_mat[gid==offsets_mat[,"Indiv ID"], "Offsets"];
		min_off=min(offsets);
		offsets_mat[gid==offsets_mat[,"Indiv ID"], "Offsets"]=offsets-min_off;
	}

	return(offsets_mat);
}

load_summary_file=function(fname){
        cat("Loading Summary Table: ", fname, "\n");
        inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", row.names=1))
        counts_mat=inmat[,2:(ncol(inmat))];
        return(counts_mat);
}

normalize=function(counts){
        totals=apply(counts, 1, sum);
        num_samples=nrow(counts);
        normalized=matrix(0, nrow=nrow(counts), ncol=ncol(counts));

        for(i in 1:num_samples){
                normalized[i,]=counts[i,]/totals[i];
        }

        colnames(normalized)=colnames(counts);
        rownames(normalized)=rownames(counts);
        return(normalized);
}

simplify_matrix_categories=function(normalized_mat, top=4){
	keep_cat=character();
	num_samp=nrow(normalized_mat);
	samp_names=rownames(normalized_mat);
	
	for(i in 1:num_samp){
		#cat(samp_names[i], "\n");
		abund=sort(normalized_mat[i,], decreasing=T);
		top_cat=(names(abund)[1:top]);	
		#print(top_cat);
		#cat("\n");
		keep_cat=c(keep_cat, top_cat);
	}

	uniq_keep_cat=unique(keep_cat);
	cat("Top ", top, " across all samples:\n", sep="");
	print(uniq_keep_cat);

	# Keep categories across top categories
	keep_mat=normalized_mat[,uniq_keep_cat];

	# Sort categories in decreasing order
	avg_abund=apply(keep_mat, 2, mean);
	sort_ix=order(avg_abund, decreasing=T);
	keep_mat=keep_mat[, sort_ix];
	return(keep_mat);
	
}

###############################################################################

plot_dist=function(x, y, width=20, abundances){
	
	rect(
		xleft=x-width/2,
		ybottom=0,
		xright=x+width/2,
		ytop=1,
		lwd=.01,
		col="grey"
	);

	num_abund=length(abundances);
	prev=0;
	for(i in 1:num_abund){
		rect(
			xleft=x-width/2,
			ybottom=prev,
			xright=x+width/2,
			ytop=prev+abundances[i],
			lwd=.01,
			col=i
		);	
		prev=prev+abundances[i];
	}
		
}

plot_sample_distances_by_individual=function(diversity_arr, normalized_mat, offsets_mat, col_assign, category_colors, ind_colors){
	sorted_sids=sort(rownames(offsets_mat));
	offsets_mat=offsets_mat[sorted_sids,];

	# Get Unique Groups
	groups=sort(unique(offsets_mat[,"Indiv ID"]));
	num_groups=length(groups);


	def_par=par(no.readonly=T);
	par(mfrow=c(3,1));

	# Get range of offsets
	offset_ranges=range(offsets_mat[,"Offsets"]);
	cat("Offset Range:\n");
	print(offset_ranges);

	diversity_ranges=range(diversity_arr);
	cat("Diversity Range:\n");
	print(diversity_ranges);

	# Get closest between offsets
	min_period=numeric();
	for(i in 1:num_groups){
		grp_subset=which(offsets_mat[,"Indiv ID"]==groups[i]);
                offsets=offsets_mat[grp_subset, "Offsets"];
		min_period=min(diff(offsets), min_period)
	}
	cat("Minimum period: ", min_period, "\n");

	cat_names=colnames(normalized_mat);

	# Plot individual samples
	for(i in 1:num_groups){

		# Grab members from group
		cat("Plotting: ", groups[i], "\n");
		grp_subset=which(offsets_mat[,"Indiv ID"]==groups[i]);
		num_members=length(grp_subset);

		# Subset offsets, and sort by offset
		offset_info=offsets_mat[grp_subset,];
		sort_ix=order(offset_info[,"Offsets"]);
		offset_info=offset_info[sort_ix,];
		print(offset_info);

		###############################################################

		# Subset diversity
		subset_samples=rownames(offset_info);
		subset_diversity=diversity_arr[subset_samples];
		print(subset_diversity);


		# Plot Diversity
		palette(ind_colors);
		plot(offset_info[,"Offsets"], subset_diversity, main=groups[i],
			 xlab="Time", ylab="Diversity", type="l", col=col_assign[groups[i]], lwd=2.5,
			 xlim=offset_ranges, ylim=diversity_ranges);

		points(offset_info[c(1,1, num_members),"Offsets"], subset_diversity[c(1,1, num_members)], 
			type="p", pch=c(17, 1, 15), cex=c(1, 2, 1.25));
		points(offset_info[,"Offsets"], subset_diversity, type="b", pch=16, cex=.5, lwd=.1);

		###############################################################

		# Plot abundances
		palette(category_colors);
		subset_norm=normalized_mat[subset_samples,];
		plot(offset_info[,"Offsets"], subset_diversity, main=groups[i],
			 xlab="Time", ylab="Taxa", type="n", col=i, lwd=2,
			 xlim=offset_ranges, ylim=c(0,1));

		for(t in 1:num_members){
			plot_dist(offset_info[subset_samples[t],"Offsets"], y=0, 
				abundances=normalized_mat[subset_samples[t],], width=min_period);
		}

		###############################################################

		# Plot color key/legend
		num_in_key=min(35, length(cat_names));
		plot(0, 0, main=groups[i],
			 xlab="", ylab="", type="n", col=i, lwd=2,
			 xlim=c(0,10), ylim=c(0,num_in_key), xaxt="n", yaxt="n");

		for(j in 1:num_in_key){
			rect(0, j-1, .9, j+.9-1, col=j, lwd=.1);
			text(1, j-1+.4, labels=cat_names[j], pos=4, cex=.8);
		}
		
	
	}
	par(def_par);
}

###############################################################################

plot_sample_distances_by_group=function(diversity_arr, normalized_mat, offsets_mat, col_assign, ind_colors){
	sorted_sids=sort(rownames(offsets_mat));
	offsets_mat=offsets_mat[sorted_sids,];


	# Get Num Cohorts
	cohorts=sort(unique(offsets_mat[,"Group ID"]));	
	num_cohorts=length(cohorts);
	cat("Number of Cohorts: ", num_cohorts, "\n");
	print(cohorts);
	cat("\n");

	# Get range of offsets
	offset_ranges=range(offsets_mat[,"Offsets"]);
	cat("Offset Range:\n");
	print(offset_ranges);

	# Get range of diversity
	diversity_ranges=range(diversity_arr);
	cat("Diversity Range:\n");
	print(diversity_ranges);

	# Set up plots per page
	def_par=par(no.readonly=T);
	par(mfrow=c(num_cohorts,1));

	# Set palette for individuals
	palette(ind_colors);

	for(g in 1:num_cohorts){

		cat("Plotting: ", cohorts[g], "\n");
		plot(0, 0, main=cohorts[g],
			 xlab="Time", ylab="Diversity", type="n", 
			 xlim=offset_ranges, ylim=diversity_ranges);

		coh_offset_mat=offsets_mat[ offsets_mat[,"Group ID"]==cohorts[g], ];
		print(coh_offset_mat);

		# Get Unique Inidividuals
		indivs=sort(unique(coh_offset_mat[,"Indiv ID"]));
		num_indivs=length(indivs);
		cat("Number of Individuals: ", num_indivs, "\n");
		print(indivs);
		cat("\n");

		# Plot individual samples
		for(i in 1:num_indivs){

			# Grab from group
			cat("Plotting: ", indivs[i], "\n");
			ind_subset=which(coh_offset_mat[,"Indiv ID"]==indivs[i]);
			num_timepts=length(ind_subset);

			# Subset offsets, and sort by offset
			offset_info=coh_offset_mat[ind_subset,];
			sort_ix=order(offset_info[,"Offsets"]);
			offset_info=offset_info[sort_ix,];
			print(offset_info);

			###############################################################

			# Subset diversity
			subset_samples=rownames(offset_info);
			subset_diversity=diversity_arr[subset_samples];
			print(subset_diversity);

			# Plot Diversity

			points(offset_info[c(1,1, num_timepts),"Offsets"], subset_diversity[c(1,1, num_timepts)], 
				type="p", pch=c(17, 1, 15), cex=c(1, 2, 1.25), col=col_assign[indivs[i]]);
			points(offset_info[,"Offsets"], subset_diversity, type="l", lwd=2.5, col=col_assign[indivs[i]]);
			points(offset_info[,"Offsets"], subset_diversity, type="l", lwd=.1, col="black");

		}
	}
	par(def_par);
}

###############################################################################
###############################################################################

offset_mat=load_offset(OffsetFileName);
print(offset_mat);

###############################################################################

counts_mat=load_summary_file(InputFileName);
#print(counts_mat);

###############################################################################

offset_mat_samples=rownames(offset_mat);
counts_mat_samples=rownames(counts_mat);
shared=intersect(offset_mat_samples, counts_mat_samples);

cat("\n\n");
cat("Samples not represented in summary table file:\n");
print(setdiff(counts_mat_samples, shared));
cat("Samples not represented in offsets file:\n");
print(setdiff(offset_mat_samples, shared));
cat("\n\n");

offset_mat=offset_mat[shared,];
counts_mat=counts_mat[shared,];
indiv_ids=unique(offset_mat[,"Indiv ID"]);
num_indiv=length(indiv_ids);

###############################################################################

normalized_mat=normalize(counts_mat);
#print(normalized_mat);

diversity_arr=diversity(normalized_mat, "shannon");
#print(diversity_arr);

###############################################################################

# simplify matrix
simplified_mat=simplify_matrix_categories(normalized_mat, top=3);
num_simp_cat=ncol(simplified_mat);

# Get color assignments
get_colors=function(num_col){
	colors=hsv(seq(0,1,length.out=num_col), c(1,.5), c(1,.75,.5));
	color_mat_dim=ceiling(sqrt(num_col));
	color_pad=rep("grey", color_mat_dim^2);
	color_pad[1:num_col]=colors;
	color_mat=matrix(color_pad, nrow=color_mat_dim, ncol=color_mat_dim);
	colors=as.vector(t(color_mat));
	colors=colors[colors!="grey"];
}

category_colors=get_colors(num_simp_cat);
ind_colors=get_colors(num_indiv);

col_assign=1:num_indiv;
names(col_assign)=indiv_ids;


###############################################################################

plot_sample_distances_by_individual(diversity_arr, simplified_mat, offset_mat, col_assign, category_colors, ind_colors);
plot_sample_distances_by_group(diversity_arr, simplified_mat, offset_mat, col_assign, ind_colors);

##############################################################################

cat("Done.\n")
dev.off();
warn=warnings();
if(length(warn)){
	print(warn);
}
q(status=0)
