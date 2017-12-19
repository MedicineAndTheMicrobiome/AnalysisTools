#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"factor_file", "f", 2, "character",
	"top_categories", "t", 2, "character",
	"output_file", "o", 2, "character",
	"diversity_type", "d", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

DEF_DIVERSITY="tail";
TOP_CATEGORIES=4;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	[-f <factor file>]\n",
	"	[-t <top categories to display, default=", TOP_CATEGORIES, ">]\n",
	"	[-o <output file root name>]\n",
	"	[-d <diversity, default=", DEF_DIVERSITY, ".]\n",
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

if(length(opt$output_file)>0){
	OutputFileRoot=opt$output_file;
}else{
	OutputFileRoot=InputFileName;
	OutputFileRoot=gsub("\\.summary_table\\.tsv$", "", OutputFileRoot);
	OutputFileRoot=gsub("\\.summary_table\\.xls$", "", OutputFileRoot);
	cat("No output file root specified.  Using input file name as root.\n");
}

###############################################################################

OutputFileRoot=paste(OutputFileRoot, ".", substr(DiversityType, 1, 4), sep="");
OutputPDF = paste(OutputFileRoot, ".div_ts.pdf", sep="");
cat("Output PDF file name: ", OutputPDF, "\n", sep="");
pdf(OutputPDF,width=8.5,height=14)

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

plot_legend=function(categories){

	orig.par=par(no.readonly=T);
	par(mar=c(0,0,0,0));
	num_cat=length(categories);
	plot(0,0, type="n", ylim=c(-10,0), xlim=c(0,30), bty="n", xaxt="n", yaxt="n");
	legend(0,0, legend=rev(categories), fill=rev(1:num_cat), cex=.7);
	par(mar=orig.par$mar);
}



###############################################################################

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T){
	num_row=nrow(mat);
	num_col=ncol(mat);

	cat("Num Rows: ", num_row, "\n");
	cat("Num Cols: ", num_col, "\n");

	mat=mat[rev(1:num_row),];

	num_colors=50;
	color_arr=rainbow(num_colors, start=0, end=4/6);
	if(high_is_hot){
		color_arr=rev(color_arr);
	}
	
	remap=function(in_val, in_range, out_range){
		in_prop=(in_val-in_range[1])/(in_range[2]-in_range[1])	
		out_val=in_prop*(out_range[2]-out_range[1])+out_range[1];
		return(out_val);
	}

	if(is.na(plot_min)){
		plot_min=min(mat);	
	}
	if(is.na(plot_max)){
		plot_max=max(mat);
	}
	cat("Plot min/max: ", plot_min, "/", plot_max, "\n");

	par(oma=c(12,12,0,1));
	par(mar=c(0, 0, 2, 0));
	plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), xaxt="n", yaxt="n", bty="n", xlab="", ylab="", main=title);
	
	# x-axis
	axis(side=1, at=seq(.5, num_col-.5, 1), labels=colnames(mat), las=2);
	axis(side=2, at=seq(.5, num_row-.5, 1), labels=rownames(mat), las=2);

	if(log_col){
		plot_min=log10(plot_min+.0125);
		plot_max=log10(plot_max+.0125);
	}

	for(x in 1:num_col){
		for(y in 1:num_row){
			
			if(log_col){
				col_val=log10(mat[y,x]+.0125);
			}else{
				col_val=mat[y,x];
			}

			remap_val=remap(col_val, c(plot_min, plot_max), c(1, num_colors));
			col_ix=ceiling(remap_val);
	
			rect(x-1, y-1, (x-1)+1, (y-1)+1, border=NA, col=color_arr[col_ix]);
			text(x-.5, y-.5, sprintf("%0.4f", mat[y,x]), srt=45);
		}
	}
	
}

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
###############################################################################

factors_mat=load_factors(FactorFileName);
#print(factors_mat);

###############################################################################

counts_mat=load_summary_file(InputFileName);
#print(counts_mat);

###############################################################################

factors_samples=rownames(factors_mat);
counts_samples=rownames(counts_mat);
shared=intersect(factors_samples, counts_samples);

cat("\n\n");
cat("Samples not represented in summary table file:\n");
print(setdiff(counts_samples, shared));
cat("Samples not represented in offsets file:\n");
print(setdiff(factors_samples, shared));
cat("\n\n");

num_shared=length(shared);
cat("Number of Shared Samples: ", num_shared, "\n");

factors_mat=factors_mat[shared,];
counts_mat=counts_mat[shared,];

###############################################################################

normalized_mat=normalize(counts_mat);
#print(normalized_mat);

if(DiversityType=="tail"){
	diversity_arr=apply(normalized_mat, 1, tail_statistic);
}else{
	diversity_arr=diversity(normalized_mat, DiversityType);
}
#print(diversity_arr);

###############################################################################

# simplify matrix
simplified_mat=simplify_matrix_categories(normalized_mat, top=NumTopCategories);
num_simp_cat=ncol(simplified_mat);
simp_categories=colnames(simplified_mat);

cat("Number of Simplified Abundances: ", num_simp_cat, "\n");

category_colors=get_colors(num_simp_cat);
palette(category_colors);

plot_cols=8;
plot_rows=4;
tot_plots_per_page=plot_cols*plot_rows;
layout_mat=matrix(1:tot_plots_per_page, byrow=T, nrow=plot_rows, ncol=plot_cols);
layout_mat=rbind(layout_mat, rep(tot_plots_per_page+1, plot_cols));
layout_mat=rbind(layout_mat, rep(tot_plots_per_page+1, plot_cols));
print(layout_mat);
layout(layout_mat);

par(oma=c(.5,.5,.5,.5));
par(mar=c(1,1,1,1));

i=0;
while(i<num_shared){
	for(y in 1:plot_rows){
		for(x in 1:plot_cols){
			if(i<num_shared){
				sample=shared[i+1];
				plot(0,0, c(-1,1), ylim=c(0,1), type="n", bty="n", xaxt="n", yaxt="n");
				title(main=sample, cex.main=.7);
				abundances=simplified_mat[sample,,drop=F];
				plot_dist(0, 0, width=1, abundances);
			}else{
				plot(0,0, c(-1,1), ylim=c(0,1), type="n", bty="n", xaxt="n", yaxt="n");
			}
			i=i+1;
		}
	}
	cat("Plotting legend...\n");
	plot_legend(simp_categories);
}

###############################################################################

cat("Done.\n")
dev.off();
warn=warnings();
if(length(warn)){
	print(warn);
}
q(status=0)
