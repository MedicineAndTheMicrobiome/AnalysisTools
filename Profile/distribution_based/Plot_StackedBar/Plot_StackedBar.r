#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"factor_file", "f", 2, "character",
	"factor_subset", "M", 2, "character",
	"top_categories", "T", 2, "character",
	"output_file", "o", 2, "character",
	"diversity_type", "d", 2, "character",
	"shorten_category_names", "s", 2, "character",
	"crossing_string", "c", 2, "character",
	"label_threshold", "l", 2, "numeric",
	"tag_name", "t", 2, "character",
	"harmonizing_colormap", "h", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

DEF_DIVERSITY="tail";
TOP_CATEGORIES=4;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	[-f <factor file>]\n",
	"	[-M <list of factors/variables to focus on (filename)>]\n",
	"	[-T <top categories to display, default=", TOP_CATEGORIES, ">]\n",
	"	[-o <output file root name>]\n",
	"	[-d <diversity, default=", DEF_DIVERSITY, ".]\n",
	"	[-s <shorten category names, with separator in double quotes (default=\"\")>]\n",
	"	[-l <label abundances greater than specified threshold, default=1.0, recommended 0.01]\n",
	"\n",
	"	[-c <crossing/interactions list, e.g., \"var1,var2,var3\" >]\n",
	"\n",
	"	[-t <tag name>]\n",
	"\n",
	"	[-h <harmonizing color map file, from previous run>]\n",
	"\n",
	"	This script will read in the summary table\n",
	"	and the factor file.\n",
	"\n",
	"	For each factor (column), a set of stacked\n",
	"	barplots will be generated for each factor level\n",
	"	If there is no factor file specified, then\n",
	"	then the average distribution will be plotted\n",
	"\n",
	"	Note if a factor file is not specified, the code will\n",
	"	try to guess the categories by parsing the sample IDs based\n",
	"	on periods (.) and underscore (_).",
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

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;

if(length(opt$factor_file)){
	FactorFileName=opt$factor_file;
}else{
	FactorFileName=NULL;
}

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

if(length(opt$tag_name)){
        TagName=opt$tag_name;
        cat("Setting TagName Hook: ", TagName, "\n");
        setHook("plot.new",
                function(){
                        cat("Hook called.\n");
                        if(par()$page==T){
                                oma_orig=par()$oma;
                                exp_oma=oma_orig;
                                exp_oma[1]=max(exp_oma[1], 1);
                                par(oma=exp_oma);
                                mtext(paste("[", TagName, "]", sep=""), side=1, line=exp_oma[1]-1,
                                        outer=T, col="steelblue4", font=2, cex=.8, adj=.97);
                                par(oma=oma_orig);
                        }
                }, "append");

}else{
        TagName="";
}

if(length(opt$harmonizing_colormap)){
	HarmonizingColorMap=opt$harmonizing_colormap;
}else{
	HarmonizingColorMap="";
}

###############################################################################

OutputFileRoot=paste(OutputFileRoot, ".", substr(DiversityType, 1, 4), ".t", NumTopCategories,  sep="");
OutputPDF = paste(OutputFileRoot, ".stckd_bp.pdf", sep="");
cat("Output PDF file name: ", OutputPDF, "\n", sep="");
pdf(OutputPDF,width=8.5*1.25,height=14)

output_grouped_abundances_fh=file(paste(OutputFileRoot, ".stckd_bp.abd.tsv", sep=""), "w");

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

	# If there is a "Remaining" category, remove it
	cat_names=colnames(normalized_mat);
	rem_bool=(tolower(cat_names)=="remaining");
	rem_ix=which(rem_bool);	
	if(length(rem_ix)==1){
		cat("Detected \"Remaining\" category...\n");
		cat("Excluding from the top categories.\n");
		normalized_mat=normalized_mat[,-rem_ix, drop=F];
	}
	
	top=min(top, ncol(normalized_mat));

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
	keep_mat=normalized_mat[,uniq_keep_cat, drop=F];

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

plot_dist=function(x, y, width=20, abundances, num_ticks=3, label_abund=0, color_map){
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
				col=color_map[cat_name[i]]
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
					paste(cat_name[i], " [", round(abundances[i], 4), "]" , sep=""), 
					cex=20*label_abund);
			}

			prev=prev+abundances[i];
		}

	}
		
}

plot_legend=function(categories, size=NULL, num_ticks=3, max_labels=40, color_map){

	orig.par=par(no.readonly=T);

	par(mar=c(0,2,0,0));
	num_cat=length(categories);
	orig_num_cat=num_cat;
	orig_categories=categories;

	cat("Legend: Num categories:", num_cat, "\n");
	orig_v_plot_diff=0;
	if(num_cat>=max_labels){
		categories=categories[1:(max_labels-1)];
		num_cat=length(categories);
		orig_v_plot_diff=(num_cat-max_labels);
	}

	if(is.null(size)){
		size=min(1, 45/num_cat);
	}

	plot(0,0, type="n", ylim=c(15,0), xlim=c(0,30), 
		bty="n", 
		xaxt="n", yaxt="n");

	leg_info=legend(0, 0, legend=rev(c(categories, "Remaining")), 
		fill=rev(c(color_map[categories[1:num_cat]], "grey")), cex=size, pt.lwd=.1);

	if(num_cat>=num_ticks){
		# Compute tick positions
		tick_pos=seq(1, orig_num_cat, length.out=num_ticks+2);
		tick_pos=ceiling(tick_pos[2:(num_ticks+1)])
		xleft=leg_info$rect$left;
		xright=(xleft+leg_info$text$x[1])/4;

		# Get categories of original tick pos
		ticked_cat_names=character(num_ticks);

		for(i in 1:num_ticks){
			ticked_cat_names[i]=orig_categories[tick_pos[num_ticks-i+1]];
		}
		cat("Ticked Category Names:\n");
		print(ticked_cat_names);
		
		# Tick the categories
		for(i in 1:num_cat){
			if(any(categories[i]==ticked_cat_names)){
				ypos=leg_info$text$y[i];	
				points(c(xleft, xright), c(ypos, ypos), type="l", lwd=1);
			}
		}
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

plot_abundance_matrix=function(abd_mat, title="", plot_cols=8, plot_rows=5, 
	samp_size=c(), divname="diversity", median_diversity=c(), mean_diversity=c(),
	label_threshold=0, color_map){
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

	layout_mat=cbind(layout_mat, matrix(tot_plots_per_page+1, nrow=plot_rows, ncol=4));

	#layout_mat=rbind(layout_mat, rep(tot_plots_per_page+1, plot_cols));
	#cat("Layout Matrix:\n");
	#print(layout_mat);
	layout(layout_mat);

	orig.par=par(no.readonly=T);
	par(oma=c(.5,.5,3.5,.5));
	par(mar=c(1,1,1,1));

	dist_bar_width=2;
	label_offset1=-(dist_bar_width/2 + .20);
	label_offset2=-(dist_bar_width/2 + .5);

	mean_sample_name_length=mean(nchar(sample_names));
	stagger_labels=(mean_sample_name_length>20)

	i=0;
	while(i<num_samples){
		for(y in 1:plot_rows){
			for(x in 1:plot_cols){
				if(i<num_samples){
					sample=sample_names[i+1];
					plot(0,0, xlim=c(-1.5,1.5), ylim=c(0,1), type="n", bty="n", xaxt="n", yaxt="n");

					smp_nm_len=nchar(sample);
					smp_nm_cex=.5*min(1, 21/smp_nm_len);
					
					if(stagger_labels){
						if(i%%2){
							l=.3;
						}else{
							l=-.3;
						}
					}else{
						l=0;
					}
		
					mtext(sample, line=l, cex=smp_nm_cex, font=2);
					

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
					plot_dist(0, 0, width=dist_bar_width, abundances, 
						num_ticks=3,
						label_abund=label_threshold, color_map=color_map);
				}else{
					plot(0,0, type="n", bty="n", xaxt="n", yaxt="n");
				}
				i=i+1;
			}
		}
		cat("Plotting legend...\n");

		cat("Num Categories in Legend: ", num_cat, "\n");

		plot_legend(cat_names, size=NULL, color_map=color_map);
		mtext(text=title, side=3, outer=T, cex=2, font=2, line=.5);
	}
	par(orig.par);
}

plot_diversity_barplot=function(title, diversity_name, samp_size,
	mean_diversity, diversity_95lb, diversity_95ub){

	grp_names=names(mean_diversity);
	num_grps=length(grp_names);

	# Calculate max diversity for ylim
	div_max=max(diversity_95ub[!is.na(diversity_95ub)], mean_diversity);
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

plot_diversity_barplot_signf=function(title, diversity_name, grpd_div_list, alpha=.05, sample_glyps=F){

	
	group_names=names(grpd_div_list);
	num_grps=length(group_names);

	# Precompute pairwise wilcoxon pvalues
	pval_mat=matrix(1, nrow=num_grps, ncol=num_grps);
	#diag(pval_mat)=1;
	signf=numeric();
	for(grp_ix_A in 1:num_grps){
		for(grp_ix_B in 1:num_grps){
			if(grp_ix_A<grp_ix_B){
				res=wilcox.test(grpd_div_list[[grp_ix_A]], grpd_div_list[[grp_ix_B]]);

				pv=ifelse(is.na(res$p.value), 1, res$p.value);

				pval_mat[grp_ix_A, grp_ix_B]=pv;
				#pval_mat[grp_ix_B, grp_ix_A]=res$p.value;

				if(pv<=alpha){
					signf=rbind(signf, c(grp_ix_A, grp_ix_B, res$p.value));
				}
			}
		}	
	}

	num_signf=nrow(signf);

	signf_by_row=apply(pval_mat, 1, function(x){sum(x<alpha)});
	
	cat("Alpha", alpha, "\n");
	cat("Significance by Row:\n");
	print(signf_by_row);
	
	num_signf_rows=sum(signf_by_row>0);
	cat("Num Rows to plot:", num_signf_rows, "\n");
	
	#print(signf);
	cat("Num Significant: ", num_signf, "\n");

	#signf_mat=apply(pval_mat, 1:2, 
	#	function(x){ 
	#		if(x<.001){return("***")}
	#		if(x<.01){return("**")}
	#		if(x<.05){return("*")}
	#		else{return("")}
	#	}
	#);
	
	#print(signf_mat, quote=F);	

	# Compute 95% CI around mean
	num_bs=320;
	grp_means=numeric(num_grps);
	ci95=matrix(NA, nrow=num_grps, ncol=2);
	samp_size=numeric(num_grps);
	for(grp_ix in 1:num_grps){
		grp_means[grp_ix]=mean(grpd_div_list[[grp_ix]]);

		num_samp=length(grpd_div_list[[grp_ix]]);
		
		if(num_samp>=40){	
			meds=numeric(num_bs);
			for(i in 1:num_bs){
				meds[i]=mean(sample(grpd_div_list[[grp_ix]], replace=T));
				
			}
			ci95[grp_ix,]=quantile(meds, c(.025, .975));
		}else{
			ci95[grp_ix,]=rep(mean(grpd_div_list[[grp_ix]]),2);
		}

		samp_size[grp_ix]=num_samp;
	}

	cat("Group Means:\n");
	print(grp_means);
	cat("Group Median 95% CI:\n");
	print(ci95);

	# Estimate spacing for annotations
	annot_line_prop=1/10; # proportion of pl
	max_95ci=max(ci95[,2], na.rm=T);
	datamax=max_95ci;
	space_for_annotations=datamax*annot_line_prop*(num_signf_rows+1);
	horiz_spacing=annot_line_prop*datamax;

	# Start plot
	par(mar=c(8,5,4,3));
	mids=barplot(grp_means, ylim=c(0, datamax+space_for_annotations));
	title(ylab=paste("Mean ", diversity_name, " with Bootstrapped 95% CI", sep=""));
	title(main=title, cex.main=1.5);
	title(main="with Wilcoxon rank sum test (difference between group means) p-values", 
		line=.25, cex.main=1, font.main=3);

	bar_width=mean(diff(mids));
	qbw=bar_width/4;

	# Label x-axis
	text(mids-par()$cxy[1]/2, rep(-par()$cxy[2]/2, num_grps),
		group_names, srt=-45, xpd=T, pos=4,
		cex=min(c(1,.7*bar_width/par()$cxy[1])));
	
	# Scatter
	if(sample_glyps){
		for(grp_ix in 1:num_grps){
			pts=grpd_div_list[[grp_ix]];
			numpts=length(pts);
			points(
				#rep(mids[grp_ix], numpts),
				mids[grp_ix]+rnorm(numpts, 0, bar_width/8),
				pts, col="grey25", cex=.5, type="p");
		}
	}

	# label CI's
	for(grp_ix in 1:num_grps){
		if(samp_size[grp_ix]>=40){
			points(
				c(mids[grp_ix]-qbw, mids[grp_ix]+qbw),
				rep(ci95[grp_ix, 2],2), type="l", col="blue");
			points(
				c(mids[grp_ix]-qbw, mids[grp_ix]+qbw),
				rep(ci95[grp_ix, 1],2), type="l", col="blue");
			points(
				rep(mids[grp_ix],2),
				c(ci95[grp_ix, 1], ci95[grp_ix,2]), type="l", col="blue");
		}
	}

	# label sample size
	for(grp_ix in 1:num_grps){
		text(mids[grp_ix], 0, paste("n =",samp_size[grp_ix]), cex=.85, font=3, adj=c(.5,-1));
	}

	connect_significant=function(A, B, ypos, pval){
		abline(h=ypos);
	}

	print(pval_mat);

	#abline(h=datamax+(horiz_spacing*1:num_signf_rows), col="grey80");

	sigchar=function(x){
		if(x<=.0001){
			return("***");
		}else if(x<=.001){
			return("**");
		}else if(x<=.01){
			return("*");
		}else{
			return("");
		}
	}
	
	row_ix=1;
	for(i in 1:(num_grps-1)){
	
		pvalrow=pval_mat[i,];
		#print(pvalrow);
		
		signf_pairs=(pvalrow<alpha);
		if(any(signf_pairs)){
			signf_grps=which(signf_pairs);
			cat("Pairs: ", i, " to:\n");
			print(signf_grps);

			y_offset=datamax+horiz_spacing*row_ix;

			# Draw line between left/reference to each paired signf grp
			points(c(
				mids[i], mids[max(signf_grps)]), 
				rep(y_offset,2),
				type="l", lend="square"
			);

			# Mark left/ref group
			points(
				rep(mids[i],2),
				c(y_offset,y_offset-horiz_spacing/4),
				type="l", lwd=3, lend="butt");

			# Mark each signf paired reference group
			for(pair_ix in signf_grps){
				points(
					rep(mids[pair_ix],2),
					c(y_offset,y_offset-horiz_spacing/4),
					type="l", lwd=1, lend="butt");

				
				# label pvalue
				paird_pval=sprintf("%5.4f", pvalrow[pair_ix]);
				text(mids[pair_ix], y_offset, paird_pval, 
					adj=c(.5, -1), cex=.7);
				text(mids[pair_ix], y_offset, sigchar(pvalrow[pair_ix]), 
					adj=c(.5, -1.25), cex=1);
			}

			row_ix=row_ix+1;

		}

	}

}

#------------------------------------------------------------------------------

write_abundances_to_fh=function(abundances_fh, abd_mat, title="",
	samp_size=c(), divname="diversity", median_diversity=c(), mean_diversity=c()){

	cat(file=abundances_fh, title, "\n", sep="");
	cat("Writing: ", title, "\n", sep="");

	Remaining=apply(abd_mat, 1, function(x){1-sum(x)});
	abd_mat=cbind(abd_mat, Remaining);

	num_grps=nrow(abd_mat);
	num_cats=ncol(abd_mat);
	grp_names=rownames(abd_mat);
	cat_names=colnames(abd_mat);

	# Spacer
	cat(file=abundances_fh, "\n");

	# Column Headers
	cat(file=abundances_fh, "\t", paste(c(
		paste(cat_names, collapse="\t"), 
		"n", 
		paste("median(", divname, ")", sep=""),
		paste("mean(", divname, ")", sep="")
	 ), collapse="\t"), "\n", sep="");

	# Sample Group values
	for(i in 1:num_grps){
		cat(file=abundances_fh, 
			paste(c(grp_names[i], sprintf("%.4g", abd_mat[i,])), collapse="\t"), sep="");

		cat(file=abundances_fh, 
			"\t",
			paste(c(
				sprintf("%i", samp_size[i]), 
				sprintf("%5.3f", median_diversity[i]), 
				sprintf("%5.3f", mean_diversity[i])
			),  collapse="\t"),
			"\n", sep="");
	}
	
	cat(file=abundances_fh, "\n");

}

###############################################################################

orig_counts_mat=load_summary_file(InputFileName);
sample_depths=apply(orig_counts_mat, 1, sum);
#print(counts_mat);

###############################################################################

if(!is.null(FactorFileName)){
	orig_factors_mat=load_factors(FactorFileName);
}else{
	cat("WARNING:  No factor file specified.  Just guessing at groupings now.\n");
	
	num_sumtab_samp=nrow(orig_counts_mat);
	samp_ids=rownames(orig_counts_mat);
	max_tokens=0;
	
	split_list=list();
	for(i in 1:num_sumtab_samp){
		split_list[[i]]=strsplit(samp_ids[i], "[\\.\\_]")[[1]];
		max_tokens=max(max_tokens, length(split_list[[i]]));
	}

	guess_mat=as.data.frame(matrix(NA, nrow=num_sumtab_samp, ncol=max_tokens));
	rownames(guess_mat)=samp_ids;

	for(i in 1:num_sumtab_samp){
		sp=split_list[[i]];
		guess_mat[i, 1:length(sp)]=sp;
	}

	colnames(guess_mat)=paste("ID", 1:max_tokens, sep="");
	orig_factors_mat=guess_mat;
	print(orig_factors_mat);

		
}
#print(factors_mat);

# Append depths to factors
tmp_samp_ids=rownames(orig_factors_mat);
log10_sample_depths=log10(sample_depths[tmp_samp_ids]);
orig_factors_mat=cbind(log10_sample_depths, orig_factors_mat);

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

factors_mat=orig_factors_mat[shared,, drop=F];
counts_mat=orig_counts_mat[shared,, drop=F];

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

max_excl_show=50;

if(length(excl_to_st)>max_excl_show){
	excl_to_st_andmore="...(more)";
}else{
	excl_to_st_andmore="";
}

if(length(excl_to_fct)>max_excl_show){
	excl_to_fct_andmore="...(more)";
}else{
	excl_to_fct_andmore="";
}

plot_text(c(
	paste("Summary Table File: ", InputFileName),
	paste("Factor File: ", FactorFileName),
	paste("Output File Root: ", OutputFileRoot),
	paste("Color Map: ", HarmonizingColorMap),
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
	capture.output(print(head(excl_to_st,max_excl_show))),
	excl_to_st_andmore,
	"",
	"Samples exclusive to Factor Table:",
	capture.output(print(head(excl_to_fct,max_excl_show))),
	excl_to_fct_andmore

));	

###############################################################################
# Handle Colors

if(HarmonizingColorMap==""){
	category_colors=get_colors(num_simp_cat);

	names(category_colors)=colnames(simplified_mat);

	write.table(category_colors, file=paste(OutputFileRoot, ".color_map.tsv", sep=""), quote=F, sep="\t",
		row.names=T, col.names=F);
}else{
	colors_table=read.table(HarmonizingColorMap, as.is=T, comment.char="");
	category_colors=colors_table[,2];
	names(category_colors)=colors_table[,1];
}

print(category_colors);
#palette(category_colors);

###############################################################################

# Plot diversity vs. depth

plot_depth_vs_diversiy=function(cts_mat, div_arr){

	#print(cts_mat);
	#print(div_arr);	

	samp_ids=rownames(cts_mat);
	div_arr=div_arr[samp_ids];

	samp_depth=apply(cts_mat, 1, sum);

	par(mfrow=c(3,1));

	sqrt_samp_depth=sqrt(samp_depth);
	lg10_div_arr=log10(div_arr+1);	

	fit=lm(lg10_div_arr~sqrt_samp_depth);

	hist(sqrt_samp_depth, main="Sample Depth", xlab="sqrt(Sample Depth)", breaks=20);
	hist(lg10_div_arr, main="Diversity", xlab="Log10(Diversity)", breaks=20);
	
	par(mar=c(5,5,5,5));

	sqrt_depth_pretty=sort(c(unique(c(pretty(sqrt_samp_depth), sqrt(c(100, 500, 1000))))));
	log10_diversity_pretty=pretty(lg10_div_arr);

	plot(sqrt_samp_depth, lg10_div_arr, main="Diversity vs. Sample Depth",
		xlab="sqrt(Sample Depth)", ylab="Log10(Diversity)", xaxt="n", yaxt="n");

	abline(fit, col="blue");

	axis(side=1, at=sqrt_depth_pretty, labels=round(sqrt_depth_pretty,2));
	axis(side=3, at=sqrt_depth_pretty, labels=round(sqrt_depth_pretty^2,2));

	axis(side=2, at=log10_diversity_pretty, labels=round(log10_diversity_pretty,2));
	axis(side=4, at=log10_diversity_pretty, labels=round((10^log10_diversity_pretty)-1,0));

	lines(lowess(lg10_div_arr~sqrt_samp_depth), col="darkorange", lty=2, lwd=2);
	

}

plot_depth_vs_diversiy(counts_mat, diversity_arr);

###############################################################################

plot_abundance_matrix(simplified_mat, title="By Sample ID", 
	divname=DiversityType, 
	median_diversity=diversity_arr, mean_diversity=diversity_arr,
	label_threshold=LabelThreshold, color_map=category_colors);

write_abundances_to_fh(output_grouped_abundances_fh, simplified_mat, title="By Sample ID",
	samp_size=rep(1, nrow(simplified_mat)), divname=DiversityType, 
	median_diversity=diversity_arr, mean_diversity=diversity_arr);

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
		non_na_ix=!is.na(fact_val);
		fact_val=fact_val[non_na_ix];
		num_fact_val=length(fact_val);

		if(is.factor(fact_val)){
			cat(fact_name, " is a factor.\n", sep="");
			fact_lev=levels(fact_val);
		}else{
			unique_val=unique(fact_val);
			num_unique=length(unique_val);

			if(num_unique<=2 || is.character(unique_val)){
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
					" - ", hist_res$breaks[2:num_grps], sep="");
				cat("Group:\n");
				print(grp_levels);

				grp_asn=character(num_fact_val);
				lowerbounds=hist_res$breaks[1:(num_grps-1)];
				for(i in 1:num_fact_val){
					grp_asn[i]=grp_levels[max(which(fact_val[i]>=lowerbounds))];
					#cat(fact_val[i], "->", grp_levels[grp_asn[i]],"\n");	
				}
				#cat("Assigned Groups:\n");
				#print(grp_asn);

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
	factors_mat=factors_mat[,fact_subset_arr, drop=F];
}

grp_mat=map_val_to_grp(factors_mat);

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

	diversity_by_group=list();

	grp_i=1;
	for(grp in groups){
		cat("Extracting: ", grp, "\n");
		samp_ix=(which(grp==values));
		sample_sizes[grp_i]=length(samp_ix);
		sample_arr=sample_names[samp_ix];

		diversity_by_group[[grp]]=diversity_arr[sample_arr];

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
	cat("Plotting Abundances:\n");
	plot_abundance_matrix(combined_abd, title=grp_name, samp_size=sample_sizes, 
		divname=DiversityType, median_diversity=diversity_median, mean_diversity=diversity_mean,
		label_threshold=LabelThreshold, color_map=category_colors);

	cat("Writing Abundances to File:\n");
	write_abundances_to_fh(output_grouped_abundances_fh, combined_abd, title=grp_name, samp_size=sample_sizes, 
		divname=DiversityType, median_diversity=diversity_median, mean_diversity=diversity_mean);
	
	cat("Plotting Diversity Bar Plots:\n");
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


	if(num_grps<=20){
		par(mfrow=c(1,1));
		plot_diversity_barplot_signf(
			title=paste(grp_name, ": Ordered by Category", sep=""),
			diversity_name=DiversityType,
			diversity_by_group
		);
	}else{
		cat("Not generating for such a large number of groups.\n");
	}
		


	cat("\n");
	
}

close(output_grouped_abundances_fh);
dev.off();

cat("Completed per variable plots.\n");

###############################################################################

if(num_crossings>0){

	plot_2D_stacked=function(var1, var2, title_arr, grpd_factors, simplfd_sumtab, label_threshold=0, color_map){
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
				plot_dist(0, 0, width=1, combined_abd, num_ticks=3,
					label_abund=label_threshold, color_map=color_map);

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

	#----------------------------------------------------------------------

	plot_2D_diversity=function(var1, var2, title_arr, grpd_factors, diversity_arr){

		cat("\n\n");
		cat("Generating Diversity Plot for: ", var1, " x ", var2, "\n");
		num_rows=nrow(grpd_factors);
		cat("Number of samples: ", num_rows, "\n");

		var1_val=grpd_factors[,var1];
		var2_val=grpd_factors[,var2];

		var1_lev=sort(unique(var1_val[!is.na(var1_val)]));
		var2_lev=sort(unique(var2_val[!is.na(var2_val)]));

		var1_num_lev=length(var1_lev);
		var2_num_lev=length(var2_lev);

		cat(var1, ": ", var1_num_lev, " levels.\n", sep="");
		cat(var2, ": ", var2_num_lev, " levels.\n", sep="");

		samp_names=rownames(grpd_factors);
		
		cum_med=c();
		cum_val=c();

		grpd_members=list();
		grp_idx=1;

		for(j in 1:var2_num_lev){
			v2_samp_ix=(var2_val==var2_lev[j]);

			for(i in 1:var1_num_lev){
				v1_samp_ix=(var1_val==var1_lev[i]);
				
				mutual_samp_ix=(v1_samp_ix & v2_samp_ix);
				mutual_samp_names=samp_names[mutual_samp_ix];

				#cat("V1:", var1_val, "\n")
				#print(v1_samp_ix);
				#cat("V2:", var2_val, "\n")
				#print(v2_samp_ix);

				#cat("Mutual:\n");
				#print(mutual_samp_ix);
				#print(mutual_samp_names);

				mutual_samp_names=mutual_samp_names[!is.na(mutual_samp_names)];
				diversity_subset=diversity_arr[mutual_samp_names];
				quant=quantile(diversity_subset, c(.025,.5,.975), na.rm=T);

				#print(diversity_subset);

				memb_info=list();
				memb_info[["cross"]]=paste(var1_lev[i], " x ", var2_val[j], sep="");
				memb_info[["diversity"]]=diversity_subset;
				memb_info[["lb95"]]=quant[1];
				memb_info[["median"]]=quant[2];
				memb_info[["ub95"]]=quant[3]
				memb_info[["n"]]=length(diversity_subset);

				cum_med=c(cum_med, quant[2]);
				cum_val=c(cum_val, diversity_subset);

				grpd_members[[grp_idx]]=memb_info;
				grp_idx=grp_idx+1;

			}
		}

		num_plots=var2_num_lev*var1_num_lev;
		colors=rev(rainbow(num_plots, start=0, end=4/6));
		cum_med_rank=rank(cum_med);

		cat("Group medians:\n");
		cum_med=cum_med[!is.na(cum_med)];
		print(cum_med);
		cat("\n");

		cum_val=cum_val[!is.na(cum_val)];
		div_rang=range(cum_val);
		mid_y=mean(div_rang);
		cat("Plot diversity range:\n");
		print(div_rang);

		grp_idx=1;
		par(mfrow=c(var2_num_lev, var1_num_lev));
		par(mar=c(1,1,1,2));
		yticks=pretty(c(0, div_rang[2]), 8);

		for(j in 1:var2_num_lev){
			for(i in 1:var1_num_lev){
				
				plot(0, type="n", xlim=c(0,1.75), ylim=c(0, div_rang[2]),
					xaxt="n", yaxt="n", bty="n", xlab="", ylab=""
				);

				abline(h=sort(cum_med), col=colors, lwd=.25);
				abline(h=0, col="black", lwd=1.5);

				mid=barplot(grpd_members[[grp_idx]][["median"]],
					add=T, xaxt="n", yaxt="n", col=colors[cum_med_rank[grp_idx]]
				);

				num_cross_level=grpd_members[[grp_idx]][["n"]];

				str=paste("n=", num_cross_level, "  median=", round(grpd_members[[grp_idx]][["median"]],2), sep="");
				mtext(str, side=1, line=-.5, cex=.5, font=3);

				if(num_cross_level>=20){
					str=paste("95% PI = (", 
						round(grpd_members[[grp_idx]][["lb95"]],2),
						", ",
						round(grpd_members[[grp_idx]][["ub95"]],2),
					")", sep="");
					mtext(str, side=1, line=0,  cex=.5, font=3);

					# draw CI bands
					endlen=1/16;
					points(
						c(mid+3/4, mid+3/4), 
						c(grpd_members[[grp_idx]][["lb95"]], grpd_members[[grp_idx]][["ub95"]]),
						type="l", lwd=.5
					);
					points(
						c(mid+3/4-endlen, mid+3/4+endlen), 
						c(grpd_members[[grp_idx]][["lb95"]], grpd_members[[grp_idx]][["lb95"]]),
						type="l", lwd=.5
					);
					points(
						c(mid+3/4-endlen, mid+3/4+endlen), 
						c(grpd_members[[grp_idx]][["ub95"]], grpd_members[[grp_idx]][["ub95"]]),
						type="l", lwd=.5
					);

					text(mid+3/4, grpd_members[[grp_idx]][["ub95"]],
						"95% PI", cex=.35, font=3, adj=c(.5, -.2));
		
				}


				jitter=rnorm(num_cross_level, 0, .06);
				jitter=jitter-mean(jitter);

				sym_size=.7*seq(.90,1.1, length.out=num_cross_level);
				points(jitter + rep(mid, num_cross_level),
					 grpd_members[[grp_idx]][["diversity"]],
					bg="white", col="black", pch=21, cex=sym_size, lwd=.5);

				grp_idx=grp_idx+1;

				# Label levels on side and top
				if(i==1){
					axis(2, at=div_rang[2]/2, labels=var2_lev[j], tick=F, font=2);
				}
				if(i==var1_num_lev){
					axis(4, at=yticks, labels=yticks, tick=T, cex.axis=.75, las=2);
				}
				if(j==var2_num_lev){
					axis(1, at=.75, labels=var1_lev[i], tick=F, font=2);
				}
				if(j==1 && i==1){
					mtext(title_arr[1], side=3, line=3, outer=T, font=2, cex=.7);
					mtext("x",          side=3, line=2, outer=T, font=1, cex=.5);
					mtext(title_arr[2], side=3, line=1, outer=T, font=2, cex=.7);
					mtext(title_arr[3], side=3, line=0, outer=T, font=1, cex=.5);
				}
			}
		}

		cat("\n\n");
	}


	compute_95CI_median=function(x, num_bs=200){
		if(length(x)==0){return(c(NA,NA))};

		bs_med=numeric(num_bs);
		for(i in 1:num_bs){
			resamp=sample(x, replace=T);
			bs_med[i]=median(resamp);
		}	
		bnds=quantile(bs_med, c(.025, .975));
		return(bnds);
	}

	plot_2D_diversity_comparisons=function(var1, var2, title_arr, grpd_factors, diversity_arr){

		cat("\n\n");
		cat("Generating Diversity Comparison Plots for: (x-axis) ", var1, " x (y-axis)", var2, "\n");
		num_rows=nrow(grpd_factors);
		cat("Number of samples: ", num_rows, "\n");

		var1_val=grpd_factors[,var1];
		var2_val=grpd_factors[,var2];

		var1_lev=sort(unique(var1_val[!is.na(var1_val)]));
		var2_lev=sort(unique(var2_val[!is.na(var2_val)]));

		var1_num_lev=length(var1_lev);
		var2_num_lev=length(var2_lev);

		cat(var1, ": ", var1_num_lev, " levels.\n", sep="");
		cat(var2, ": ", var2_num_lev, " levels.\n", sep="");

		samp_names=rownames(grpd_factors);
		
		par(mfrow=c(3,var1_num_lev));
		par(oma=c(.5, .5, 4, 3));

		glob_margins=c(3,.15,.15,.15);
		bottom_margin=4;
		left_margin=1.75;
		right_margin=1;
	
		grpd_info=list();
		max_div=0;
		# Extract data
		for(i in 1:var1_num_lev){
			v1_samp_ix=(var1_val==var1_lev[i]);
			cat("Var1 X-axis:", var1_lev[i], "\n");

			for(j in 1:var2_num_lev){
				v2_samp_ix=(var2_val==var2_lev[j]);
				cat("  Var2 Y-axis:", var2_lev[j], "\n");

				mutual_samp_ix=(v1_samp_ix & v2_samp_ix);
				#print(mutual_samp_ix);

				mutual_samp_names=samp_names[mutual_samp_ix];
				mutual_samp_names=mutual_samp_names[!is.na(mutual_samp_names)];
				#print(mutual_samp_names);

				div=diversity_arr[mutual_samp_names];
				memb_info=list();
				memb_info[["diversity"]]=div;
				max_div=max(max_div, div);
		
				key=paste(var1_lev[i], "_x_", var2_lev[j], sep="");
				grpd_info[[key]]=memb_info;		

			}
		}

		cat("Generating Barplots...\n");
		for(i in 1:var1_num_lev){

			med_div=numeric();
			div=list();
			
			for(j in 1:var2_num_lev){
				key=paste(var1_lev[i], "_x_", var2_lev[j], sep="");
				div[[var2_lev[j]]]=grpd_info[[key]]$diversity;
				med_div[j]=median(grpd_info[[key]]$diversity);
			}

			print(div);

			#bar_mids=barplot(med_div, xaxt="n", yaxt="n", ylim=c(0, max_div*1.1));
			mlw=.75;

			# Set up margins for
			tmp_mar=glob_margins;
			if(i==var1_num_lev){
				tmp_mar[4]=right_margin;
			}else if(i==1){
				tmp_mar[2]=left_margin;
			}
			par(mar=tmp_mar);

			# Generate plot area
			plot(0,0, type="n", bty="n", 
				xaxt="n", yaxt="n", 
				xlab="", ylab="",
				xlim=c(0, var2_num_lev), ylim=c(0, max_div*1.1));

			# Plot median bar, 95%CI of Median, and scatter
			for(j in 1:var2_num_lev){
				points(	c(j-.5-mlw/2, j-.5+mlw/2),
					c(med_div[j], med_div[j]), type="l", lwd=1, lend=2, col="blue"); 

				ci=compute_95CI_median(div[[var2_lev[j]]]);
				points(	c(j-.5-mlw/3, j-.5+mlw/3),
					c(ci[1], ci[1]), type="l", lwd=.75, lend=2, col="red"); 
				points(	c(j-.5-mlw/3, j-.5+mlw/3),
					c(ci[2], ci[2]), type="l", lwd=.75, lend=2, col="red"); 

				jitter=rnorm(length(div[[var2_lev[j]]]), 0, .10);
				jitter-mean(jitter);
				points(j-.5+jitter,
					div[[var2_lev[j]]],
					bg="white", col="grey15", pch=21, cex=.2, lwd=.25);

				abline(h=0);

				
			}
	
			# Label x-axis
			for(j in 1:var2_num_lev){
				text(
					(j-.5)-par()$cxy[1]/2, 
					0-par()$cxy[2]/2,
					labels=var2_lev[j], srt=-45, xpd=T, pos=4, cex=.5);
			}

			if(i==1){
				title(ylab="Diversity", line=.75, font.lab=2, cex.lab=.75);
				title(ylab="(Medians w/ 95%CI)", line=0, font.lab=1, cex.lab=.5);
			}
			if(i==var1_num_lev){
				ats=pretty(c(0, max_div));
				axis(4, at=ats, labels=ats, las=2, cex.axis=.7);	
			}
		}

		cat("Generating Pairwise Differences...\n");
		pos=(1:var2_num_lev)-.5;
		for(i in 1:var1_num_lev){
			for(j in 1:var2_num_lev){
				key=paste(var1_lev[i], "_x_", var2_lev[j], sep="");
				div[[var2_lev[j]]]=grpd_info[[key]]$diversity;
				med_div[j]=median(grpd_info[[key]]$diversity);
			}

			# Set up margins for
			tmp_mar=glob_margins;
			if(i==var1_num_lev){
				tmp_mar[4]=right_margin;
			}else if(i==1){
				tmp_mar[2]=left_margin;
			}
			par(mar=tmp_mar);

			# Generate plot area
			plot(0,0, type="n", bty="n", xaxt="n", yaxt="n",
				xlim=c(0, var2_num_lev), ylim=c(0, var2_num_lev),
				main="", xlab="", ylab="");

			# Draw grid lines
			abline(h=1:(var2_num_lev-1), col="grey80");
			abline(v=1:(var2_num_lev-1), col="grey80");

			for(x in 1:var2_num_lev){
				for(y in 1:var2_num_lev){
					if(x==y){
						next;
					}
					text(x-.5, y-.5, labels=round(med_div[x]-med_div[y], 2), cex=.75, srt=45);
				}
			}

			
			for(j in 1:var2_num_lev){
				text(
					pos[j]-par()$cxy[1]/2, 
					0-par()$cxy[2]/2,
					labels=var2_lev[j], srt=-45, xpd=T, pos=4, cex=.5);
			}

			if(i==1){
				title(ylab="Differences", line=.75, font.lab=2, cex.lab=.75);
				title(ylab="if(Bottom-Right)<0, then R>B", line=0, font.lab=1, cex.lab=.5);
			}
			if(i==var1_num_lev){
				axis(4, at=pos, labels=var2_lev, las=2, cex.axis=.5, lty=c(0,1));
			}
		}

			
		cat("Generating Significances Matrix...\n");
		for(i in 1:var1_num_lev){
			for(j in 1:var2_num_lev){
				key=paste(var1_lev[i], "_x_", var2_lev[j], sep="");
				div[[var2_lev[j]]]=grpd_info[[key]]$diversity;
			}

			# Set up margins for
			tmp_mar=glob_margins;
			tmp_mar[1]=bottom_margin;
			if(i==var1_num_lev){
				tmp_mar[4]=right_margin;
			}else if(i==1){
				tmp_mar[2]=left_margin;
			}
			par(mar=tmp_mar/2);
			cat("Margins for Plot:\n");
			print(tmp_mar);
			# Generate plot area
			plot(0,0, type="n", bty="n", xaxt="n", yaxt="n",
				xlim=c(0, var2_num_lev), ylim=c(0, var2_num_lev),
				main="", xlab="", ylab="");

			# Draw grid lines
			abline(h=1:(var2_num_lev-1), col="grey80");
			abline(v=1:(var2_num_lev-1), col="grey80");

			for(x in 1:var2_num_lev){
				for(y in 1:var2_num_lev){

					# Skip self-self comparisons
					if(x==y){
						next;
					}

					# Compute wilcoxon difference in medians
					if(length(div[[var2_lev[x]]])==0 || length(div[[var2_lev[y]]])==0){
						pvalue=1;
					}else{
						wcres=wilcox.test(div[[var2_lev[x]]], div[[var2_lev[y]]]);
						pvalue=wcres$p.value;
					}

					# Assign text color based on p-value
					color="grey50";
					font=1;
					if(pvalue<.1){
						color="black";
					}
					if(pvalue<.05){
						color="blue";
						font=2;
					}
					if(pvalue<.01){
						color="red";
						font=2;
					}
					
					text(x-.5, y-.5, labels=round(pvalue, 2), cex=.75, srt=45,
						font=font, col=color);
				}
			}
			
			for(j in 1:var2_num_lev){
				text(
					pos[j]-par()$cxy[1]/2, 
					0-par()$cxy[2]/2,
					labels=var2_lev[j], srt=-45, xpd=T, pos=4, cex=.5);
			}

			if(i==var1_num_lev){
				axis(4, at=pos, labels=var2_lev, las=2, cex.axis=.5, lty=c(0,1));
			}


			# Bottom margin variable 1 factor levels
			title(xlab=var1_lev[i], font.lab=2, line=3);

			# Left margin labels
			if(i==1){
				title(ylab="Significances", line=.75, font.lab=2, cex.lab=.75);
				title(ylab="(Wilcoxon R.S.T.)", line=0, font.lab=1, cex.lab=.5);
			}
		}
		

	        mtext(title_arr[1], side=3, line=3, outer=T, font=2, cex=.7);
		mtext("x",          side=3, line=2, outer=T, font=1, cex=.5);
		mtext(title_arr[2], side=3, line=1, outer=T, font=2, cex=.7);
		mtext(title_arr[3], side=3, line=0, outer=T, font=1, cex=.5);

		cat("\n\n");
	}

	plot_2D_diversity_comparisons_signf_annot=function(var1, var2, title_arr, grpd_factors, diversity_arr){
		#print(var1);
		#print(var2);
		#print(title_arr);
		#print(grpd_factors);
		#print(diversity_arr);

		var1_levels=sort(setdiff(unique(grpd_factors[,var1]), NA));
		var2_levels=sort(setdiff(unique(grpd_factors[,var2]), NA));

		num_var1_levels=length(var1_levels);
		num_var2_levels=length(var2_levels);

		cat(var1, " levels (", num_var1_levels, "):\n");
		print(var1_levels);

		cat(var2, " levels (", num_var2_levels, "):\n");
		print(var2_levels);

		stat_matrix=matrix(NA, nrow=num_var1_levels, ncol=num_var2_levels, 
			dimnames=list(var1_levels, var2_levels));

		mean_mat=stat_matrix;	
		lb95_mat=stat_matrix;	
		ub95_mat=stat_matrix;	
		cnt_mat=stat_matrix;

		# Split diversity into 2D list of lists
		grpd_div_12=list();
		grpd_div_21=list();
		samp_ids=rownames(grpd_factors);
		num_bs=320;
	
		for(v1_lvl in var1_levels){
			grpd_div_12[[v1_lvl]]=list();
			l1_ix=grpd_factors[,var1]==v1_lvl;
			for(v2_lvl in var2_levels){
				l2_ix=grpd_factors[,var2]==v2_lvl;
				samp_ix=setdiff(unique(samp_ids[l1_ix&l2_ix]), NA);
				cat(v1_lvl, " x ", v2_lvl, ":\n", sep="");
				print(samp_ix);
				div=diversity_arr[samp_ix];
				grpd_div_12[[v1_lvl]][[v2_lvl]]=div;
				
				grp_size=length(div);
				mean_mat[v1_lvl, v2_lvl]=mean(div);
				cnt_mat[v1_lvl, v2_lvl]=grp_size;
				
				# Bootstrap calculate 95% CI of mean
				if(grp_size>=40){
					bs_means=numeric(num_bs)
					for(bs_ix in 1:num_bs){
						bs_means[bs_ix]=mean(sample(div, replace=T));
					}

					ci95=quantile(bs_means, c(.025, .975));
					lb95_mat[v1_lvl, v2_lvl]=ci95[1];
					ub95_mat[v1_lvl, v2_lvl]=ci95[2];
				}
				

			}

		}

		for(v2_lvl in var2_levels){
			grpd_div_21[[v2_lvl]]=list();
			for(v1_lvl in var1_levels){
				grpd_div_21[[v2_lvl]][[v1_lvl]]=grpd_div_12[[v1_lvl]][[v2_lvl]];
			}
		}

		print(mean_mat);
		print(lb95_mat);
		print(ub95_mat);
		print(cnt_mat);
		#print(grpd_div);

		plot_bar_annot=function(v1, v2, div, means, lb95s, ub95s, cnts, scat=F){

			cat("Plotting barplots:\n");

			v1_levels=rownames(means);
			v2_levels=colnames(means);
			num_v1_lvls=length(v1_levels);
			num_v2_lvls=length(v2_levels);

			tot_bars=num_v1_lvls*num_v2_lvls;
			cat("Number of Bars: ", tot_bars, "\n");

			stat_mat=matrix(NA, nrow=tot_bars, ncol=4);
			colnames(stat_mat)=c("means", "lb95s", "ub95s", "cnts");

			inner_levels=character(tot_bars);
			outer_levels=character(tot_bars);

			print(v1_levels);
			print(v2_levels);
			
			# Serialize the 2D matrix into Array
			i=1;
			for(v1_lvl in v1_levels){
				for(v2_lvl in v2_levels){
					outer_levels[i]=v1_lvl;
					inner_levels[i]=v2_lvl;
					stat_mat[i, "means"]=means[v1_lvl, v2_lvl];		
					stat_mat[i, "lb95s"]=lb95s[v1_lvl, v2_lvl];		
					stat_mat[i, "ub95s"]=ub95s[v1_lvl, v2_lvl];		
					stat_mat[i, "cnts"]=cnts[v1_lvl, v2_lvl];		
					i=i+1;
				}
			}

			pval_list=list();
			for(v1_lvl in v1_levels){

				pval_list[[v1_lvl]]=matrix(1, nrow=num_v2_lvls, ncol=num_v2_lvls,
					dimnames=list(v2_levels, v2_levels));
				
				for(i in 1:num_v2_lvls){
					idiv=div[[v1_lvl]][[v2_levels[i]]];

					for(j in 1:num_v2_lvls){

						if(i<j){
							jdiv=div[[v1_lvl]][[v2_levels[j]]];

							cat(v1_lvl, ":", v2_levels[i], " vs. ", v2_levels[j], "\n");

							if(length(idiv)>0 && length(jdiv>0)){
								wcres=wilcox.test(idiv, jdiv);
								pval_list[[v1_lvl]][i,j]=wcres$p.value;
							}else{
								pval_list[[v1_lvl]][i,j]=1;
							}

						}
				
					}
				}
			}

			maxplot_val=max(stat_mat[, c("ub95s", "means")], na.rm=T);

			if(scat==T){
				serial_div_list=list();
				i=1;
				max_div=0;
				for(v1_lvl in v1_levels){
					for(v2_lvl in v2_levels){
						serial_div_list[[i]]=div[[v1_lvl]][[v2_lvl]];		
						max_div=max(max_div, serial_div_list[[i]]);
						i=i+1;
					}
				}
				maxplot_val=max(maxplot_val, max_div);
			}

			annot_space=maxplot_val/2;
			annot_start=maxplot_val*1.05;
			ymax_wannot=annot_start+annot_space;

			# Start the plot
			par(mar=c(8, 4, 1, 1));
			mids=barplot(stat_mat[, "means"], ylim=c(0, ymax_wannot), xlab="", ylab="Diversity");

			bar_spacing=(mids[2]-mids[1]);

			# Add dashed lines between grouped bars
			between_bars=mids-bar_spacing/2;	
			if((1+num_v2_lvls)<tot_bars){
				abline(v=between_bars[seq(1+num_v2_lvls,tot_bars,num_v2_lvls)], col="grey50", lty=2);
			}

			# Proportion of intervals annotation width
			bs_div6=bar_spacing/6;
		
			# Draw scatter
			if(scat==T){
				for(i in 1:tot_bars){
					num_samp=length(serial_div_list[[i]]);
					scat_size=.7*seq(.9,1.1, length.out=num_samp);
					points(rnorm(num_samp, mids[i], bs_div6/2), 
						serial_div_list[[i]],
						cex=scat_size, bg="white", pch=21, col="darkgreen", lwd=.9);
				}
			}
			
			# Label plots
			for(i in 1:tot_bars){
				# LBs
				points(
					c(mids[i]-bs_div6, mids[i]+bs_div6),
					rep(stat_mat[i, "lb95s"],2),
					col="blue", type="l"
				);

				# UBs
				points(
					c(mids[i]-bs_div6, mids[i]+bs_div6),
					rep(stat_mat[i, "ub95s"],2),
					col="blue", type="l"
				);
		
				# mids
				points(
					rep(mids[i],2),
					c(stat_mat[i, "lb95s"],stat_mat[i, "ub95s"]),
					col="blue", type="l"
				);

				# counts
				text(mids[i], 0, 
					labels=paste(
						"n=", stat_mat[i, "cnts"], 
						"  mean=", round(stat_mat[i, "means"], 2), 
						sep=""),
					adj=c(.5,-.5), font=3, cex=.5
				);
			}
			
			# Calculate begin/end of outer group bar spacing
			outer_mids=bar_spacing*num_v2_lvls/2 + num_v2_lvls*bar_spacing*(0:(num_v1_lvls-1));

			# Label outer levels
			outer_label_cex=3/num_v1_lvls;
			text(
				outer_mids, 
				rep(0-par()$cxy[2], num_v1_lvls), 
				labels=v1_levels, xpd=T, cex=outer_label_cex, font=2);

			# Label inner levels
			text(
				mids-par()$cxy[1]/2, 
				0-par()$cxy[2]*2,
				labels=inner_levels, srt=-45, xpd=T, pos=4, cex=.8);
			
			print(pval_list);


			# Draw annotation for significance between members of same group
			for(v1_ix in 1:num_v1_lvls){

				v1_lvl=v1_levels[v1_ix];
				cat("Within group:", v1_lvl, "\n");

				for(v2_ix in 1:num_v2_lvls){
					from=v2_levels[v2_ix];
					pvalrow=pval_list[[v1_lvl]][v2_ix,];
					signfrow=pvalrow<.1;
					print(signfrow);

					if(any(signfrow)){
						to=which(signfrow);

						fromx=mids[(v1_ix-1)*num_v2_lvls+v2_ix];
						ypos=annot_start+annot_space*v2_ix/num_v2_lvls;

						# Draw lines between from/to
						points(
							c(fromx,
							mids[(v1_ix-1)*num_v2_lvls+max(to)]), 
							rep(ypos, 2),
							type="l");

						# From tick
						points(
							rep(fromx,2),
							c(ypos, ypos*.97), type="l", lwd=1.5, lend="butt");


						# Draw 'to' ticks and pvalues
						for(i in 1:num_v2_lvls){
							if(pvalrow[i]<.1){

								xpos=mids[(v1_ix-1)*num_v2_lvls+i];
							
								points(
									rep(xpos,2),
									c(ypos, ypos*.98),
									type="l", lwd=1, lend="butt");

								sigch="";
								if(pvalrow[i]<=.001){
									sigch="***";
								}else if(pvalrow[i]<=.01){
									sigch="**";
								}else if(pvalrow[i]<=.05){
									sigch="*";
								}
								text(xpos, ypos, sigch, adj=c(.5,.25), cex=1);
							}
						}
					}

				}
			}

		}

		# plot grouped by var 1
		par(oma=c(1,1,1,1));
		par(mfrow=c(1,1));
		plot_bar_annot(var2, var1, grpd_div_21, t(mean_mat), t(lb95_mat), t(ub95_mat), t(cnt_mat), scat=F);
		plot_bar_annot(var1, var2, grpd_div_12, mean_mat, lb95_mat, ub95_mat, cnt_mat, scat=F);

		plot_bar_annot(var2, var1, grpd_div_21, t(mean_mat), t(lb95_mat), t(ub95_mat), t(cnt_mat), scat=T);
		plot_bar_annot(var1, var2, grpd_div_12, mean_mat, lb95_mat, ub95_mat, cnt_mat, scat=T);
		
		plot_text(c(
			"Notes:",
			"",
			"p-values:",
			"          ***: p <= .001",
			"           **: p <= .01",
			"            *: p <= .05",
			" (tick only) : p <= .10",
			"",
			"P-values are calculated pairwise using Wilcoxon Rank Sum Test.",
			"",
			"Bar heights are mean diversity.",
			"Intervals are bootstrapped 95% CI around mean."
		));
	
	}

	#----------------------------------------------------------------------

	cat("Starting crossed variable plots:\n");
	print(crossing_var);


	num_uniq=integer(num_crossings);
	cat("Num of Levels:\n");
	for(i in 1:num_crossings){
		num_uniq[i]=length(unique(grp_mat[,crossing_var[i]]));
		cat(crossing_var[i], ": ", num_uniq[i], "\n");
	}
	names(num_uniq)=crossing_var;
	
	if(num_crossings==3){

		samp_ids=rownames(simplified_mat);
		
		for(excl_var in crossing_var){
			rem_cros_var=setdiff(crossing_var, excl_var);

			cat("PDF Width: ", num_uniq[1], " Height: ", num_uniq[2], "\n");
			pdf(paste(OutputFileRoot, ".", 
				rem_cros_var[1], "_x_", rem_cros_var[2], "_x_", excl_var, ".pdf", sep=""),
				width=num_uniq[rem_cros_var[1]]*3,
				height=num_uniq[rem_cros_var[2]]*3
			);

			# Plot combined
			plot_2D_stacked(rem_cros_var[1], rem_cros_var[2], 
				c(crossing_var[1], crossing_var[2]),
				grp_mat, simplified_mat, LabelThreshold, category_colors);	

			plot_2D_diversity(rem_cros_var[1], rem_cros_var[2], 
				c(crossing_var[1], crossing_var[2]),
				grp_mat, diversity_arr);

			plot_2D_diversity_comparisons(rem_cros_var[1], rem_cros_var[2], 
				c(crossing_var[1], crossing_var[2]),
				grp_mat, diversity_arr);

			plot_2D_diversity_comparisons_signf_annot(rem_cros_var[1], rem_cros_var[2], 
				c(crossing_var[1], crossing_var[2]),
				grp_mat, diversity_arr);

			
			# Plot split by excluded variable
			levels=sort(unique(grp_mat[,excl_var]));
			cat("Splitting by: ", levels, "\n");
			for(lev in levels){
				keep_ix=(grp_mat[,excl_var]==lev);
				plot_2D_stacked(rem_cros_var[1], rem_cros_var[2], 
					c(rem_cros_var[1], rem_cros_var[2], lev),
					grp_mat[keep_ix,,drop=F], simplified_mat, LabelThreshold, category_colors);	
				plot_2D_diversity(rem_cros_var[1], rem_cros_var[2], 
					c(rem_cros_var[1], rem_cros_var[2], lev),
					grp_mat[keep_ix,,drop=F], diversity_arr);
				plot_2D_diversity_comparisons(rem_cros_var[1], rem_cros_var[2], 
					c(rem_cros_var[1], rem_cros_var[2], lev),
					grp_mat[keep_ix,,drop=F], diversity_arr);
				plot_2D_diversity_comparisons_signf_annot(rem_cros_var[1], rem_cros_var[2], 
					c(rem_cros_var[1], rem_cros_var[2], lev),
					grp_mat[keep_ix,,drop=F], diversity_arr);
			}

			# Plot legend	
			par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0));
			plot_legend(colnames(simplified_mat), size=.5, color_map=color_map);

			dev.off();

		}

	}else if(num_crossings==2){
		# do AxB crossings
		cat("PDF Width: ", num_uniq[1], " Height: ", num_uniq[2], "\n");
		pdf(paste(OutputFileRoot, ".", crossing_var[1], "_x_", crossing_var[2], ".pdf", sep=""),
			width=num_uniq[1]*3,
			height=num_uniq[2]*3
		);
		
		cat("Generating plot for 2-way crossings...\n");
		plot_2D_stacked(crossing_var[1], crossing_var[2], 
			c(crossing_var[1],  crossing_var[2], ""),
			grp_mat, simplified_mat, LabelThreshold, category_colors);	

		plot_2D_diversity(crossing_var[1], crossing_var[2], 
			c(crossing_var[1],  crossing_var[2], ""), grp_mat, diversity_arr);

		plot_2D_diversity_comparisons(crossing_var[1], crossing_var[2], 
			c(crossing_var[1],  crossing_var[2], ""), grp_mat, diversity_arr);

		plot_2D_diversity_comparisons_signf_annot(crossing_var[1], crossing_var[2], 
			c(crossing_var[1],  crossing_var[2], ""), grp_mat, diversity_arr);

		# Plot legend
		par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0));
		plot_legend(colnames(simplified_mat), size=.5, color_map=category_colors);

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
