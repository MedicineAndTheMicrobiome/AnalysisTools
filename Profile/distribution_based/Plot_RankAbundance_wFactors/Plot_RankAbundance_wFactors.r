#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"factor_file", "f", 1, "character",
	"top_categories", "t", 2, "numeric",
	"output_file", "o", 2, "character",
	"diversity_type", "d", 2, "character",
	"shorten_category_names", "s", 2, "character",
	"crossing_string", "c", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

DEF_DIVERSITY="tail";
TOP_CATEGORIES=20;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	-f <factor file>\n",
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
	"	For each factor in the factor/metadata file, a set of plots will generated, with\n",
	"	the samples grouped by the levels in the factor.\n",
	"\n",
	"	Diversity types include:\n",
	"		shannon, tail, simpson, invsimpson\n",
	"\n",
	"	Each rank abundance barplot will be labeled with median and 95%CI\n",
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

NumTopCategories=TOP_CATEGORIES;
if(length(opt$top_categories)){
	NumTopCategories=opt$top_categories;
}


###############################################################################

OutputFileRoot=paste(OutputFileRoot, ".", substr(DiversityType, 1, 4), ".t", NumTopCategories,  sep="");
OutputPDF = paste(OutputFileRoot, ".rnk_abnd.pdf", sep="");
cat("Output PDF file name: ", OutputPDF, "\n", sep="");

page_width=max(8.5, (1+3+NumTopCategories/5));
cat("Page with based on number of categories to plot: ", page_width, "\n");

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

plot_ra=function(abundances, num_top_categories=10, category_colors, ymax){
	# This draws a single Rank Abundance Plot

	category_names=rownames(abundances);	
	print(abundances);

	#barplot(abundances);

}

###############################################################################

plot_rank_abundance_matrix=function(abd_mat, title="", plot_cols=8, plot_rows=4, 
	num_top_categories=10,
	samp_size=c(), divname="diversity", median_diversity=c(), mean_diversity=c()
	){
	# This function will plot a sample x abundance (summary table)
	# There will be one plot for row (sample) in the matrix

	num_cat=ncol(abd_mat);
	num_samples=nrow(abd_mat);
	sample_names=rownames(abd_mat);
	cat_names=colnames(abd_mat);
	label_samp_size=(length(samp_size)>0);
	label_diversity=(length(median_diversity)>0);

	colors=rainbow(num_cat, start=0, end=4/6);

	orig.par=par(no.readonly=T);

	par(mfrow=c(plot_rows, plot_cols));
	par(oma=c(.5,.5,3.5,.5));
	par(mar=c(1,1,1,1));

	abd_avg=apply(abd_mat, 2, mean);
	avg_ix=order(abd_avg, decreasing=T);
	

	i=0;
	while(i<num_samples){
		for(y in 1:plot_rows){
			for(x in 1:plot_cols){
				if(i<num_samples){
					sample=sample_names[i+1];


					mtext(sample, line=0, cex=.5, font=2);

					if(label_samp_size){
						n=samp_size[i+1];
						#mtext(paste("n=",n,sep=""), line=-.5, cex=.4, font=3);
					}
					if(length(samp_size) && samp_size[i+1]>1){
						#text(label_offset1, 0, paste("Median ", divname, " = ",
						#	signif(median_diversity[i+1], 4),sep=""),
						#	srt=90, adj=0, cex=.7);
						#text(label_offset2, 0, paste("Mean ", divname, " = ",
						#	signif(mean_diversity[i+1], 4),sep=""),
						#	srt=90, adj=0, cex=.7);
					}else{
						#text(label_offset1, 0, paste(divname, " = ",
						#	signif(median_diversity[i+1], 4),sep=""),
						#	srt=90, adj=0, cex=.7);
					}

					abundances=abd_mat[sample,,drop=F];

					plot_ra(abd_mat[i,], num_top_categories, colors, ymax);
				}else{
					plot(0,0, type="n", bty="n", xaxt="n", yaxt="n");
				}
				i=i+1;
			}
		}
		mtext(text=title, side=3, outer=T, cex=2, font=2, line=.5);
	}
	par(orig.par);
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

###############################################################################

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
	"",
	"Samples exclusive to Summary Table:",
	capture.output(print(excl_to_st)),
	"",
	"Samples exclusive to Factor Table:",
	capture.output(print(excl_to_fct))
));

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
				order_pad=sprintf(paste("%0",floor(log10(num_grps))+1, "i", sep=""), 
					1:(num_grps-1));
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

grp_mat=map_val_to_grp(factors_mat);
print(grp_mat);

sample_names=rownames(grp_mat);
grp_names=colnames(grp_mat);

num_categories=ncol(normalized_mat);

for(i in 1:ncol(grp_mat)){
		
	values=grp_mat[,i];
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


	# Allocate group/level merged stats and categories
	combined_abd=matrix(0, nrow=num_grps, ncol=num_categories);
	rownames(combined_abd)=groups;
	colnames(combined_abd)=colnames(normalized_mat);
	sample_sizes=numeric(num_grps);
	diversity_median=numeric(num_grps);
	diversity_95lb=numeric(num_grps);
	diversity_95ub=numeric(num_grps);
	diversity_mean=numeric(num_grps);

	# For each group/level combine categories
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

		combined_abd[grp,]=apply(normalized_mat[sample_arr,,drop=F], 2, mean);
		#print(combined_abd);	
		grp_i=grp_i+1;
	}
	names(sample_sizes)=groups;
	names(diversity_median)=groups;
	names(diversity_mean)=groups;

	#print(combined_abd);
	#print(sample_sizes);
	plot_rank_abundance_matrix(
		abd_mat=combined_abd, 
		title=grp_name, 
		num_top_categories=NumCategoriesToPlot,
		samp_size=sample_sizes, 
		divname=DiversityType, 
		median_diversity=diversity_median, 
		mean_diversity=diversity_mean
	);


	cat("\n");
	
}

dev.off();


###############################################################################

cat("Done.\n")
warn=warnings();
if(length(warn)){
	print(warn);
}
q(status=0)
