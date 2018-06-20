#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"factor_file", "f", 2, "character",
	"factor_subset", "M", 2, "character",
	"top_categories", "t", 2, "character",
	"output_file", "o", 2, "character",
	"diversity_type", "d", 2, "character",
	"shorten_category_names", "s", 2, "character",
	"crossing_string", "c", 2, "character",
	"label_threshold", "l", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

DEF_DIVERSITY="tail";
TOP_CATEGORIES=4;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	[-f <factor file>]\n",
	"	[-M <list of factors/variables to focus on>]\n",
	"	[-t <top categories to display, default=", TOP_CATEGORIES, ">]\n",
	"	[-o <output file root name>]\n",
	"	[-d <diversity, default=", DEF_DIVERSITY, ".]\n",
	"	[-s <shorten category names, with separator in double quotes (default=\"\")>]\n",
	"	[-l <label abundances greater than specified threshold, default=1.0, recommended 0.01]\n",
	"\n",
	"	[-c <crossing/interactions list, e.g., \"var1,var2,var3\" >]\n",
	"\n",
	"	This script will read in the summary table\n",
	"	and the factor file.\n",
	"\n",
	"	For each factor (column), a set of stacked\n",
	"	barplots will be generated for each factor level\n",
	"	If there is no factor file specified, then\n",
	"	then the average distribution will be plotted\n",
	"\n",
	"	Diversity types include:\n",
	"		shannon, tail, simpson, invsimpson\n",
	"\n",
	"	Each stacked barplot will be labeled with median and 95%CI\n",
	"	of the diversity index\n",
	"\n",
	"\n",
	"	Specifying the -c cross/interactions list, creates an\n",
	"	extra pdf file that slices the data according to the\n",
	"	underlying factor levels. For example, if you had variables\n",
	"	A, B, and C, with levels Aa, Ab, Ac, and Ba, Bb, and Ca, Cb,\n",
	"	respectively, plots for AxB, AxC, and BxC, would be generated.\n",
	"	Two sets of plots (PDFs) will be produced:\n",
	"		1.) AxB, AxC, BxC: Assuming the missing variable\n",
	"			doesn't matter.  i.e. in AxB, C is collapsed\n",
	"		2.) C(AxB), B(AxC), A(BxC): First split by the missing\n",
	"			variable's levels, then produced the cross plots.\n",
	"			i.e. if C has 5 levels, then there will be 5.\n",
	"			AxB plots.\n",
	"\n",
	"\n", sep="");

if(!length(opt$input_file) || !length(opt$factor_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
FactorFileName=opt$factor_file;

DiversityType=DEF_DIVERSITY;
if(length(opt$diversity_type)){
	DiversityType=opt$diversity_type;
}

NumTopCategories=TOP_CATEGORIES;
if(length(opt$top_categories)){
	NumTopCategories=opt$top_categories;
}

FactorSubset=NULL;
if(length(opt$factor_subset)){
	FactorSubset=opt$factor_subset;
}

if(length(opt$output_file)>0){
	OutputFileRoot=opt$output_file;
}else{
	OutputFileRoot=InputFileName;
	OutputFileRoot=gsub("\\.summary_table\\.tsv$", "", OutputFileRoot);
	OutputFileRoot=gsub("\\.summary_table\\.xls$", "", OutputFileRoot);
	cat("No output file root specified.  Using input file name as root.\n");
}

if(length(opt$crossing_string)>0){
	CrossingString=opt$crossing_string;
}else{
	CrossingString=NA;
}

if(length(opt$shorten_category_names)){
        ShortenCategoryNames=opt$shorten_category_names;
}else{
        ShortenCategoryNames="";
}

if(ShortenCategoryNames==TRUE){
        cat("Error:  You need to specify a delimitor to split the category names.\n");
        cat("        i.e., this -x option is not a flag, it requires a parameter.\n");
        quit(status=-1);
}

LabelThreshold=1.0;
if(length(opt$label_threshold)){
	LabelThreshold=opt$label_threshold;
	if(LabelThreshold<1.0){
		OutputFileRoot=paste(OutputFileRoot, ".lbl", sep="");
	}
}

cat("Label Threshold: ", LabelThreshold, "\n");


###############################################################################

OutputFileRoot=paste(OutputFileRoot, ".", substr(DiversityType, 1, 4), ".t", NumTopCategories,  sep="");
OutputPDF = paste(OutputFileRoot, ".stckd_bp.pdf", sep="");
cat("Output PDF file name: ", OutputPDF, "\n", sep="");
pdf(OutputPDF,width=8.5,height=14)

if(!is.na(CrossingString)){
	cat("Crossing String Specified: ", CrossingString, "\n");
	crossing_var=strsplit(CrossingString, ",")[[1]];
	num_crossings=length(crossing_var);
	for(i in 1:num_crossings){
		crossing_var[i]=gsub("^\\s+", "", crossing_var[i]);
		crossing_var[i]=gsub("\\s+$", "", crossing_var[i]);
	}
}else{
	num_crossings=0;
}

###############################################################################

load_factors=function(fname){
        factors=data.frame(read.table(fname,  header=TRUE, check.names=FALSE, row.names=1, comment.char="", quote="", sep="\t"));
        dimen=dim(factors);
        cat("Rows Loaded: ", dimen[1], "\n");
        cat("Cols Loaded: ", dimen[2], "\n");
        return(factors);
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
	# This function will reduce the matrix to only contain
	# the categories that are in the top n across all of the
	# samples.  This way, all of the categories will be represented
	# to some extent.  The returned matrix will have its
	# columns sorted by decreasing abundance.

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

plot_text=function(strings){

	orig.par=par(no.readonly=T);
        par(family="Courier");
        par(oma=rep(.1,4));
        par(mar=rep(0,4));

        num_lines=length(strings);

        top=max(as.integer(num_lines), 52);

        plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                );

        text_size=max(.01, min(.8, .8 - .003*(num_lines-52)));
        #print(text_size);

        for(i in 1:num_lines){
                #cat(strings[i], "\n", sep="");
                strings[i]=gsub("\t", "", strings[i]);
                text(0, top-i, strings[i], pos=4, cex=text_size);
        }

	par(orig.par);
}

###############################################################################

# Get color assignments
get_colors=function(num_col, alpha=1){
	colors=hsv(seq(0,1,length.out=num_col+1), c(1,.5), c(1,.75,.5), alpha=alpha);
	color_mat_dim=ceiling(sqrt(num_col));
	color_pad=rep("grey", color_mat_dim^2);
	color_pad[1:num_col]=colors[1:num_col];
	color_mat=matrix(color_pad, nrow=color_mat_dim, ncol=color_mat_dim);
	colors=as.vector(t(color_mat));
	colors=colors[colors!="grey"];
}

###############################################################################

plot_dist=function(x, y, width=20, abundances, num_ticks=3, label_abund=0){
	# This function will plot a stack box plot
	# The location is center around x, and over y, with a bar height of 1
	# If the abundance is less than label_abund, it will be labeled
	
	if(all(is.na(abundances))){

		points(
			c(x-width/2, x+width/2),
			c(y, y),
			lwd=.1, type="l"
		);

		points(
			c(x-width/2, x-width/2),
			c(y, y+.1),
			lwd=.1, type="l"
		);

		points(
			c(x+width/2, x+width/2),
			c(y, y+.1),
			lwd=.1, type="l"
		);

	}else{
		
		rect(
			xleft=x-width/2,
			ybottom=y,
			xright=x+width/2,
			ytop=y+1,
			lwd=.01,
			col="grey"
		);

		num_abund=length(abundances);
		if(is.null(dim(abundances))){
			cat_name=names(abundances);
		}else{
			cat_name=colnames(abundances);
		}

		tick_pos=seq(1, num_abund, length.out=num_ticks+2);
		tick_pos=ceiling(tick_pos[2:(num_ticks+1)]);

		prev=y;
		for(i in 1:num_abund){
			rect(
				xleft=x-width/2,
				ybottom=prev,
				xright=x+width/2,
				ytop=prev+abundances[i],
				lwd=.01,
				col=i
			);	

			# Plot tick marks on the left
			if(any(i==tick_pos)){
				xleft=x-width/2;
				ymid=prev+abundances[i]/2;

				points(
					c(xleft-width/20, xleft), 
					c(ymid, ymid),
					type="l", lwd=.5, col="black"
				);
			}

			# Label abundance if greater than label_abund
			if(abundances[i]>label_abund){
				text(x, prev+abundances[i]/2, 
					paste(cat_name[i], " [", round(abundances[i], 4), "]" , sep=""), cex=.1);
			}

			prev=prev+abundances[i];
		}

	}
		
}

plot_legend=function(categories, size=.7, num_ticks=3){

	orig.par=par(no.readonly=T);

	par(mar=c(0,0,0,0));
	num_cat=length(categories);
	plot(0,0, type="n", ylim=c(-10,0), xlim=c(0,30), bty="n", xaxt="n", yaxt="n");
	leg_info=legend(0,0, legend=rev(c(categories, "Remaining")), 
		fill=rev(c(1:num_cat, "grey")), cex=size, pt.lwd=.1);

	# Compute tick positions
	tick_pos=seq(1, num_cat, length.out=num_ticks+2);
	tick_pos=ceiling(tick_pos[2:(num_ticks+1)])
	xleft=leg_info$rect$left;
	xright=(xleft+leg_info$text$x[1])/4;

	for(i in 1:num_ticks){
		ypos=leg_info$text$y[tick_pos[num_ticks-i+1]];	
		points(c(xleft, xright), c(ypos, ypos), type="l", lwd=1);
	}

	par(mar=orig.par$mar);
}



###############################################################################

tail_statistic=function(x){
        sorted=sort(x, decreasing=TRUE);
        norm=sorted/sum(x);
        n=length(norm);
        tail=0;
        for(i in 1:n){
                tail=tail + norm[i]*((i-1)^2);
        }
        return(sqrt(tail));
}

###############################################################################

plot_abundance_matrix=function(abd_mat, title="", plot_cols=8, plot_rows=4, 
	samp_size=c(), divname="diversity", median_diversity=c(), mean_diversity=c(),
	label_threshold=0){
	# This function will plot a sample x abundance (summary table)
	# There will be one plot for row (sample) in the matrix

	num_cat=ncol(abd_mat);
	num_samples=nrow(abd_mat);
	sample_names=rownames(abd_mat);
	cat_names=colnames(abd_mat);
	label_samp_size=(length(samp_size)>0);
	label_diversity=(length(median_diversity)>0);

	# Set up layout
	tot_plots_per_page=plot_cols*plot_rows;
	layout_mat=matrix(1:tot_plots_per_page, byrow=T, nrow=plot_rows, ncol=plot_cols);
	layout_mat=rbind(layout_mat, rep(tot_plots_per_page+1, plot_cols));
	layout_mat=rbind(layout_mat, rep(tot_plots_per_page+1, plot_cols));
	#cat("Layout Matrix:\n");
	#print(layout_mat);
	layout(layout_mat);

	orig.par=par(no.readonly=T);
	par(oma=c(.5,.5,3.5,.5));
	par(mar=c(1,1,1,1));

	dist_bar_width=2;
	label_offset1=-(dist_bar_width/2 + .20);
	label_offset2=-(dist_bar_width/2 + .5);

	i=0;
	while(i<num_samples){
		for(y in 1:plot_rows){
			for(x in 1:plot_cols){
				if(i<num_samples){
					sample=sample_names[i+1];
					plot(0,0, xlim=c(-1.5,1.5), ylim=c(0,1), type="n", bty="n", xaxt="n", yaxt="n");

					mtext(sample, line=0, cex=.5, font=2);

					if(label_samp_size){
						n=samp_size[i+1];
						mtext(paste("n=",n,sep=""), line=-.5, cex=.4, font=3);
					}
					if(length(samp_size) && samp_size[i+1]>1){
						text(label_offset1, 0, paste("Median ", divname, " = ",
							signif(median_diversity[i+1], 4),sep=""),
							srt=90, adj=0, cex=.7);
						text(label_offset2, 0, paste("Mean ", divname, " = ",
							signif(mean_diversity[i+1], 4),sep=""),
							srt=90, adj=0, cex=.7);
					}else{
						text(label_offset1, 0, paste(divname, " = ",
							signif(median_diversity[i+1], 4),sep=""),
							srt=90, adj=0, cex=.7);
					}

					abundances=abd_mat[sample,,drop=F];
					plot_dist(0, 0, width=dist_bar_width, abundances, label_abund=label_threshold);
				}else{
					plot(0,0, type="n", bty="n", xaxt="n", yaxt="n");
				}
				i=i+1;
			}
		}
		cat("Plotting legend...\n");
		plot_legend(cat_names);
		mtext(text=title, side=3, outer=T, cex=2, font=2, line=.5);
	}
	par(orig.par);
}

plot_diversity_barplot=function(title, diversity_name, samp_size,
	mean_diversity, diversity_95lb, diversity_95ub){

	grp_names=names(mean_diversity);
	num_grps=length(grp_names);

	# Calculate max diversity for ylim
	div_max=max(diversity_95ub[!is.na(diversity_95ub)]);
	if(is.na(div_max) || !is.finite(div_max)){
		cat("95% CI Upperbound not finite...\n");
		div_max=max(mean_diversity);
	}	


	par(mar=c(10, 5, 4, 1));
	mids=barplot(mean_diversity, main=title, las=2, 
		ylim=c(0, div_max*1.1),
		ylab=paste("Mean ", diversity_name, " w/ 95% CI", sep=""), 
		names.arg=rep("", num_grps));
	bar_width=mids[2]-mids[1];

	# Label sample sizes per group
	ylim=par()$usr[4];
	text(mids, ylim/20, paste("n=", samp_size, sep=""), font=3,
		cex=min(c(1,.5*bar_width/par()$cxy[1])));

	# Label x-axis
	text(mids-par()$cxy[1]/2, rep(-par()$cxy[2]/2, num_grps),
		grp_names, srt=-45, xpd=T, pos=4,
		cex=min(c(1,.7*bar_width/par()$cxy[1])));

	# mark 95lb/95ub
	intvl_col="blue";
	for(i in 1:num_grps){
		if(!is.na(diversity_95lb[i])){
			m=mids[i];
			lb=diversity_95lb[i];
			ub=diversity_95ub[i];
			points(
				c(m-bar_width/6, m+bar_width/6),
				c(lb, lb),
				type="l", col=intvl_col, lwd=.75);
			points(
				c(m-bar_width/6, m+bar_width/6),
				c(ub, ub),
				type="l", col=intvl_col, lwd=.75);
			points(
				c(m, m),
				c(lb, ub),
				type="l", col=intvl_col, lwd=.5);
		}
	}


}


###############################################################################

orig_factors_mat=load_factors(FactorFileName);
#print(factors_mat);

###############################################################################

orig_counts_mat=load_summary_file(InputFileName);
#print(counts_mat);

###############################################################################

orig_factors_samples=rownames(orig_factors_mat);
orig_counts_samples=rownames(orig_counts_mat);
shared=intersect(orig_factors_samples, orig_counts_samples);

cat("\n\n");
cat("Samples not represented in summary table file:\n");
excl_to_st=setdiff(orig_counts_samples, shared);
print(excl_to_st);
cat("Samples not represented in offsets file:\n");
excl_to_fct=setdiff(orig_factors_samples, shared);
print(excl_to_fct);
cat("\n\n");

num_shared=length(shared);
cat("Number of Shared Samples: ", num_shared, "\n");

factors_mat=orig_factors_mat[shared,];
counts_mat=orig_counts_mat[shared,];

###############################################################################

if(ShortenCategoryNames!=""){
        full_names=colnames(counts_mat);
        splits=strsplit(full_names, ShortenCategoryNames);
        short_names=character();
        for(i in 1:length(full_names)){
                short_names[i]=tail(splits[[i]], 1);

                short_names[i]=gsub("_unclassified$", "_uncl", short_names[i]);
                short_names[i]=gsub("_group", "_grp", short_names[i]);
        }
        colnames(counts_mat)=short_names;

        cat("Names have been shortened.\n");
}

normalized_mat=normalize(counts_mat);
#print(normalized_mat);

# simplify matrix
simplified_mat=simplify_matrix_categories(normalized_mat, top=NumTopCategories);
num_simp_cat=ncol(simplified_mat);
cat("Number of Simplified Abundances: ", num_simp_cat, "\n");

if(DiversityType=="tail"){
	diversity_arr=apply(normalized_mat, 1, tail_statistic);
}else{
	diversity_arr=diversity(normalized_mat, DiversityType);
}
cat("Diversity:\n");
print(diversity_arr);

plot_text(c(
	paste("Summary Table File: ", InputFileName),
	paste("Factor File: ", FactorFileName),
	"",
	paste("Diversity Index:", DiversityType),
	"",
	"Summary Table:",
	paste("    Num Samples:", nrow(orig_counts_mat)),
	paste(" Num Categories:", ncol(orig_counts_mat)),
	"",
	"Factor Table:",
	paste("    Num Samples:", nrow(orig_factors_mat)),
	paste("    Num Factors:", ncol(orig_factors_mat)),
	"",
	paste("Shared Samples:", num_shared),
	"",
	paste("Number of Top Categories from each sample to summarize:", NumTopCategories),
	paste("Number of Unique Categories across top categories extracted:", num_simp_cat),
	"",
	"Samples exclusive to Summary Table:",
	capture.output(print(excl_to_st)),
	"",
	"Samples exclusive to Factor Table:",
	capture.output(print(excl_to_fct))
));

###############################################################################

category_colors=get_colors(num_simp_cat);
palette(category_colors);

plot_abundance_matrix(simplified_mat, title="By Sample ID", 
	divname=DiversityType, 
	median_diversity=diversity_arr, mean_diversity=diversity_arr,
	label_threshold=LabelThreshold);

###############################################################################

map_val_to_grp=function(fact_mat){
	# This function will convert a factor matrix, into
	# a grouping matrix to reduce the number of continous values

	num_factors=ncol(factors_mat);
	num_values=nrow(factors_mat);
	fact_names=colnames(factors_mat);
	map_mat=as.data.frame(fact_mat);

	for(fidx in 1:num_factors){
		fact_name=fact_names[fidx];
		cat("\nMapping on: ", fact_name, "\n");

		fact_val=factors_mat[,fidx];
		print(fact_val);
		non_na_ix=!is.na(fact_val);
		fact_val=fact_val[non_na_ix];
		num_fact_val=length(fact_val);

		if(is.factor(fact_val)){
			cat(fact_name, " is a factor.\n", sep="");
			fact_lev=levels(fact_val);
			print(fact_lev);
		}else{
			unique_val=unique(fact_val);
			num_unique=length(unique_val);

			if(num_unique<=2){
				cat(fact_name, ": few enough unique values, NOT grouping\n", sep="");
				map_mat[,fidx]=as.character(map_mat[,fidx]);
			}else{
				cat(fact_name, ": too many unique values, grouping...\n", sep="");
				hist_res=hist(fact_val,breaks=nclass.Sturges(fact_val), plot=F);
				cat("Values\n");
				print(fact_val);
				num_grps=length(hist_res$breaks);
				cat("Num Groups: ", num_grps, "\n");

				# Attach an order to the grp level to keep sorting correct
				order_pad=sprintf(paste("%0",floor(log10(num_grps))+1, "i", sep=""), 1:(num_grps-1));
				print(order_pad);
			
				grp_levels=paste(
					paste("(", order_pad, ") ", sep=""),
					hist_res$breaks[1:(num_grps-1)], 
					"-", hist_res$breaks[2:num_grps], sep="");
				cat("Group:\n");
				print(grp_levels);

				grp_asn=character(num_fact_val);
				lowerbounds=hist_res$breaks[1:(num_grps-1)];
				for(i in 1:num_fact_val){
					grp_asn[i]=grp_levels[max(which(fact_val[i]>=lowerbounds))];
					#cat(fact_val[i], "->", grp_levels[grp_asn[i]],"\n");	
				}
				cat("Assigned Groups:\n");
				print(grp_asn);

				# Convert strings to factors
				grp_as_factor=factor(grp_asn, levels=grp_levels, ordered=F);
				# Initialize an array
				tmp=rep(grp_as_factor[1],num_values);
				# Copy values over
				tmp[non_na_ix]=grp_asn;
				# Replace NAs
				tmp[setdiff(1:num_values, which(non_na_ix))]=NA;
				map_mat[,fidx]=tmp;
			}

			cat("Unique Val:", unique_val, "\n");
		
		}	
	}
	return(map_mat);
}

###############################################################################
# Plot stacked bar plots across each of the factors

if(!is.null(FactorSubset)){
	fact_subset_arr=scan(FactorSubset, what=character(), comment.char="#");
	cat("Focusing on these factors:\n");
	print(fact_subset_arr);
	factors_mat=factors_mat[,fact_subset_arr];
}

grp_mat=map_val_to_grp(factors_mat);
print(grp_mat);

sample_names=rownames(grp_mat);
grp_names=colnames(grp_mat);

for(i in 1:ncol(grp_mat)){
		
	values=(grp_mat[,i]);
	all_levels=levels(values);
	groups=sort(unique(values[!is.na(values)]));
	grp_name=grp_names[i];
	num_grps=length(groups);

	cat("Plotting: ", grp_name, "\n");
	cat("Available Groups: \n");
	print(groups);

	cat("Num Available Groups: ", num_grps, "\n");
	cat("Possible Groups (levels): \n");
	print(all_levels);

	if(num_grps==0){
		cat("No informations... skipping...\n");
		next;
	}

	combined_abd=matrix(0, nrow=num_grps, ncol=num_simp_cat);
	rownames(combined_abd)=groups;
	colnames(combined_abd)=colnames(simplified_mat);
	sample_sizes=numeric(num_grps);
	diversity_median=numeric(num_grps);
	diversity_95lb=numeric(num_grps);
	diversity_95ub=numeric(num_grps);
	diversity_mean=numeric(num_grps);
	grp_i=1;
	for(grp in groups){
		cat("Extracting: ", grp, "\n");
		samp_ix=(which(grp==values));
		sample_sizes[grp_i]=length(samp_ix);
		sample_arr=sample_names[samp_ix];

		diversity_median[grp_i]=median(diversity_arr[sample_arr]);
		diversity_mean[grp_i]=mean(diversity_arr[sample_arr]);
		if(sample_sizes[grp_i]>1){
			se=sd(diversity_arr[sample_arr])/sqrt(length(sample_arr));
			diversity_95lb[grp_i]=diversity_mean[grp_i]-se*1.96;
			diversity_95ub[grp_i]=diversity_mean[grp_i]+se*1.96;
		}else{
			diversity_95lb[grp_i]=NA;
			diversity_95ub[grp_i]=NA;
		}

		combined_abd[grp,]=apply(simplified_mat[sample_arr,,drop=F], 2, mean);
		#print(combined_abd);	
		grp_i=grp_i+1;
	}
	names(sample_sizes)=groups;
	names(diversity_median)=groups;
	names(diversity_mean)=groups;

	#print(combined_abd);
	#print(sample_sizes);
	plot_abundance_matrix(combined_abd, title=grp_name, samp_size=sample_sizes, 
		divname=DiversityType, median_diversity=diversity_median, mean_diversity=diversity_mean,
		label_threshold=LabelThreshold);


	par(mfrow=c(3,1));
	plot_diversity_barplot(title=paste(grp_name, ": Ordered by Category", sep=""), 
		diversity_name=DiversityType, samp_size=sample_sizes,
                mean_diversity=diversity_mean, 
		diversity_95lb, diversity_95ub);

	order_ix=order(diversity_mean);
	plot_diversity_barplot(title=paste(grp_name, ": Ordered by Diversity", sep=""), 
		diversity_name=DiversityType, samp_size=sample_sizes[order_ix],
                mean_diversity=diversity_mean[order_ix], 
		diversity_95lb[order_ix], diversity_95ub[order_ix]);

	order_ix=order(sample_sizes);
	plot_diversity_barplot(title=paste(grp_name, ": Ordered by Sample Size", sep=""), 
		diversity_name=DiversityType, samp_size=sample_sizes[order_ix],
                mean_diversity=diversity_mean[order_ix], 
		diversity_95lb[order_ix], diversity_95ub[order_ix]);


	cat("\n");
	
}
dev.off();

cat("Completed per variable plots.\n");

###############################################################################

if(num_crossings>0){

	plot_2D_stacked=function(var1, var2, title_arr, grpd_factors, simplfd_sumtab, label_threshold=0){
		# This function will plot a matrix of stacked bar plots
		# Var1 levels will be the number of columns
		# Var2 levels will be the nubmer of rows

		fact_nonnas=!is.na(grpd_factors[,var1]) & !is.na(grpd_factors[,var2]);
		nonna_samp_id=rownames(grpd_factors)[fact_nonnas];
		grpd_factors=grpd_factors[nonna_samp_id,];
		simplfd_sumtab=simplfd_sumtab[nonna_samp_id,];

		v1_levels=sort(levels(grpd_factors[,var1]));
		v2_levels=sort(levels(grpd_factors[,var2]));

		num_v1_levels=length(v1_levels);
		num_v2_levels=length(v2_levels);

		samp_ids=rownames(grpd_factors);

		par(mar=c(.25,.25,.25,.25));
		par(oma=c(3,3,5,1));

		if(num_v1_levels==0){
			cat(var1, ": does not contain levels.\n");
			v1_levels=unique(grpd_factors[,var1]);
			num_v1_levels=length(v1_levels);
		}
		if(num_v2_levels==0){
			cat(var2, ": does not contain levels.\n");
			v2_levels=unique(grpd_factors[,var2]);
			num_v2_levels=length(v2_levels);
		}
		
		par(mfrow=c(num_v2_levels, num_v1_levels));

		for(j in 1:num_v2_levels){
			v2_lev=v2_levels[j];
			cat("Rows: ", v2_lev, "\n");
			v2_samp_id=(grpd_factors[,var2]==v2_lev);

			for(i in 1:num_v1_levels){
				v1_lev=v1_levels[i];
				cat("Cols: ", v1_lev, "\n");
			
				# Pull IDs
				v1_samp_id=(grpd_factors[,var1]==v1_lev);

				plot(0,0, type="n", 
					bty="n", yaxt="n", xaxt="n",
					xlim=c(-.5,.5), ylim=c(-.15,1));

				v1_x_v2_samp_id=(v1_samp_id & v2_samp_id);
				num_samp=sum(v1_x_v2_samp_id);

				combined_abd=apply(simplfd_sumtab[samp_ids[v1_x_v2_samp_id],,drop=F], 2, mean);
				plot_dist(0, 0, width=1, combined_abd, label_abund=label_threshold);

				# Label number of samples
				text(0, 0, pos=1, paste("n=", num_samp, sep=""), cex=.6);
				

				# Label levels on side and top
				if(i==1){
					axis(2, at=.5, labels=v2_lev, tick=F, font=2);
				}
				if(j==num_v2_levels){
					axis(1, at=0, labels=v1_lev, tick=F, font=2);
				}
				if(j==1 && i==1){
					mtext(title_arr[1], side=3, line=3, outer=T, font=2, cex=.7);
					mtext("x",          side=3, line=2, outer=T, font=1, cex=.5);
					mtext(title_arr[2], side=3, line=1, outer=T, font=2, cex=.7);
					mtext(title_arr[3], side=3, line=0, outer=T, font=1, cex=.5);
				}

			}
		}
	}

	cat("Starting crossed variable plots:\n");
	print(crossing_var);


	num_uniq=integer(num_crossings);
	cat("Num of Levels:\n");
	for(i in 1:num_crossings){
		#print(grp_mat[,crossing_var[i]]);
		num_uniq[i]=length(unique(grp_mat[,crossing_var[i]]));
		cat(crossing_var[i], ": ", num_uniq[i], "\n");
	}
	names(num_uniq)=crossing_var;
	
	if(num_crossings==3){

		samp_ids=rownames(simplified_mat);
		
		for(excl_var in crossing_var){
			rem_cros_var=setdiff(crossing_var, excl_var);

			cat("PDF Width: ", num_uniq[1], " Height: ", num_uniq[2], "\n");
			pdf(paste(OutputFileRoot, ".", rem_cros_var[1], "_x_", rem_cros_var[2], "_x_", excl_var, ".pdf", sep=""),
				width=max(num_uniq[rem_cros_var[1]], 3),
				height=num_uniq[rem_cros_var[2]]*2
			);

			# Plot combined
			plot_2D_stacked(rem_cros_var[1], rem_cros_var[2], 
				c(crossing_var[1], crossing_var[2]),
				grp_mat, simplified_mat, LabelThreshold);	

			
			# Plot split by excluded variable
			levels=sort(unique(grp_mat[,excl_var]));
			cat("Splitting by: ", levels, "\n");
			for(lev in levels){
				keep_ix=(grp_mat[,excl_var]==lev);
				plot_2D_stacked(rem_cros_var[1], rem_cros_var[2], 
					c(rem_cros_var[1], rem_cros_var[2], lev),
					grp_mat[keep_ix,,drop=F], simplified_mat, LabelThreshold);	
			}

			# Plot legend	
			par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0));
			plot_legend(colnames(simplified_mat), size=.5);

			dev.off();

		}

	}else if(num_crossings==2){
		# do AxB crossings
		cat("PDF Width: ", num_uniq[1], " Height: ", num_uniq[2], "\n");
		pdf(paste(OutputFileRoot, ".", crossing_var[1], "_x_", crossing_var[2], ".pdf", sep=""),
			width=max(num_uniq[1], 3),
			height=num_uniq[2]*2
		);
		
		cat("Generating plot for 2-way crossings...\n");
		plot_2D_stacked(crossing_var[1], crossing_var[2], 
			c(crossing_var[1],  crossing_var[2], ""),
			grp_mat, simplified_mat, LabelThreshold);	

		# Plot legend
		par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0));
		plot_legend(colnames(simplified_mat), size=.5);

		dev.off();
	}
}

###############################################################################

cat("Done.\n")
warn=warnings();
if(length(warn)){
	print(warn);
}
q(status=0)
