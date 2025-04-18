#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

options(useFancyQuotes=F);

params=c(
	"summary_file", "s", 1, "character",
	"num_variables", "p", 2, "numeric",
	"contains_remaining", "R", 2, "logical",
	"additional_categories", "a", 2, "character",
	"output_file_root", "o", 2, "character",
	"shorten_category_names", "x", 2, "character",
	"alpha", "f", 2, "numeric",
	"prefix", "P", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

NUM_TOP_CAT=35;
ALPHA=.05;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table for taxa/function>\n",
	"	[-p <number of top taxonomic/categorical variables, default=", NUM_TOP_CAT, ">]\n",
	"	[-R (pay attention to 'remaining' category)]\n",
	"	[-a <filename of additional categories to include>] Not implemented yet\n",
	"	[-o <output filename root>]\n",
	"	[-x <shorten categories names with delimitor>]\n",
	"	[-f <Alpha for filtering correlation heatmap, default=", ALPHA, ">]\n",
	"	[-P <prefix for category names>]\n",
	"\n",
	"Compute ALR transformed values and export to tsv file.\n",
	"Some diagnostic plots will be generated.\n",
	"\n",
	"The -P option will allow the user to specify a prefix to the category names\n",
	"in case more than one summary table is eventually going to be used and the\n",
	"exported ALR values will be used together in the same analysis.  To keep the\n",
	"taxa/category names unique.\n",
	"\n", sep="");

if(!length(opt$summary_file)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$output_file_root)){
	OutputRoot=gsub(".summary_table.xls", "", opt$summary_file);
	OutputRoot=gsub(".summary_table.tsv", "", OutputRoot);
}else{
	OutputRoot=opt$output_file_root;
}

if(!length(opt$num_variables)){
	NumVariables=NUM_TOP_CAT;
}else{
	NumVariables=opt$num_variables;
}

if(length(opt$contains_remaining)){
	UseRemaining=T;
}else{
	UseRemaining=F;
}

if(length(opt$shorten_category_names)){
	ShortenCategoryNames=opt$shorten_category_names;
}else{
	ShortenCategoryNames="";
}

if(length(opt$alpha)){
	Alpha=opt$alpha;
}else{
	Alpha=ALPHA;
}

if(length(opt$prefix)){
	Prefix=opt$prefix;
}else{
	Prefix="";
}

SummaryFile=opt$summary_file;

OutputRoot=paste(OutputRoot, ".alr", sep="");

cat("\n");
cat("Summary File: ", SummaryFile, "\n", sep="");
cat("Output File: ", OutputRoot, "\n", sep="");
cat("\n");
cat("Number of MALR Variables: ", NumVariables, "\n", sep="");
cat("Use Remaining? ", UseRemaining, "\n");
cat("Shorten Category Names: '", ShortenCategoryNames, "'\n", sep="");
cat("Category Prefix: '", Prefix, "'\n", sep="");
cat("Alpha: ", Alpha, "\n");
cat("\n");

if(ShortenCategoryNames==TRUE){
	cat("Error:  You need to specify a delimitor to split the category names.\n");
	cat("        i.e., this -x option is not a flag, it requires a parameter.\n");
	quit(status=-1);
}

options(width=100);
cat("Text Line Width: ", options()$width, "\n", sep="");

##############################################################################

load_summary_file=function(fname){
	inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", quote="", row.names=1))
	counts_mat=inmat[,2:(ncol(inmat))];

	# Clean category names a little
	cat_names=colnames(counts_mat);
	cat_names=gsub("-", "_", cat_names);
	colnames(counts_mat)=cat_names;
	
	cat("Num Categories in Summary Table: ", ncol(counts_mat), "\n", sep="");
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

extract_top_categories=function(ordered_normalized, top, additional_cat=c()){

        num_samples=nrow(ordered_normalized);
        num_categories=ncol(ordered_normalized);

        cat("Samples: ", num_samples, "\n");
        cat("Categories: ", num_categories, "\n");

        num_top_to_extract=min(num_categories-1, top);

        cat("Top Requested to Extract: ", top, "\n");
        cat("Columns to Extract: ", num_top_to_extract, "\n");

        # Extract top categories requested
        top_cat=ordered_normalized[,1:num_top_to_extract];

        if(length(additional_cat)){
                cat("Additional Categories to Include:\n");
                print(additional_cat);
        }else{
                cat("No Additional Categories to Extract.\n");
        }

        # Extract additional categories
        # :: Make sure we can find the categories
        available_cat=colnames(ordered_normalized);
        missing_cat=setdiff(additional_cat, available_cat);
        if(length(missing_cat)){
                cat("Error: Could not find categories: \n");
                print(missing_cat);
                quit(status=-1);
        }

        # :: Remove categories we have already extracted in the top N
        already_extracted_cat=colnames(top_cat);
        extra_cat=setdiff(additional_cat, already_extracted_cat);

        num_extra_to_extract=length(extra_cat);
        cat("Num Extra Categories to Extract: ", num_extra_to_extract, "\n");

        # Allocate/Prepare output matrix
        num_out_mat_cols=num_top_to_extract+num_extra_to_extract+1;
        out_mat=matrix(0, nrow=num_samples, ncol=num_out_mat_cols);
        rownames(out_mat)=rownames(ordered_normalized);
        colnames(out_mat)=c(already_extracted_cat, extra_cat, "Remaining");

        # Copy over top and additional categories, and compute remaining
        all_cat_names=c(already_extracted_cat, extra_cat);
        out_mat[,all_cat_names]=ordered_normalized[,all_cat_names];
        out_mat[,"Remaining"]=apply(out_mat, 1, function(x){1-sum(x)});

        return(out_mat);

}

additive_log_rato=function(ordered_matrix){
# Assumes last column will be the denominator

	num_cat=ncol(ordered_matrix);
	num_samp=nrow(ordered_matrix);

	denominator=ordered_matrix[,num_cat];
	alr_mat=matrix(0, nrow=num_samp, ncol=(num_cat-1));

	for(i in 1:num_samp){
		alr_mat[i,]=log(ordered_matrix[i,1:(num_cat-1)]/denominator[i]);
	}

	rownames(alr_mat)=rownames(ordered_matrix)
	colnames(alr_mat)=head(colnames(ordered_matrix), num_cat-1);

	alr_struct=list();
	alr_struct[["transformed"]]=alr_mat;
	alr_struct[["denominator"]]=denominator;

	return(alr_struct);
}

plot_text=function(strings){
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

	cat(strings, sep="\n");
}

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, deci_pts=4, 
	label_zeros=T, counts=F, value.cex=2, 
	plot_col_dendr=F,
	plot_row_dendr=F
){

        num_row=nrow(mat);
        num_col=ncol(mat);

	row_names=rownames(mat);
	col_names=colnames(mat);

	orig.par=par(no.readonly=T);

	cat("Painting Matrix: ", title, "\n");
        cat("Num Rows: ", num_row, "\n");
        cat("Num Cols: ", num_col, "\n");


	if(num_row==0 || num_col==0){
		plot(0, type="n", xlim=c(-1,1), ylim=c(-1,1), xaxt="n", yaxt="n", bty="n", xlab="", ylab="",
			main=title);
		text(0,0, "No data to plot...");
		return();
	}

	# Flips the rows, so becuase origin is bottom left
        mat=mat[rev(1:num_row),, drop=F];

	# Generate a column scheme
        num_colors=50;
        color_arr=rainbow(num_colors, start=0, end=4/6);
        if(high_is_hot){
                color_arr=rev(color_arr);
        }

	# Provide a means to map values to an (color) index 
        remap=function(in_val, in_range, out_range){
                in_prop=(in_val-in_range[1])/(in_range[2]-in_range[1])
                out_val=in_prop*(out_range[2]-out_range[1])+out_range[1];
                return(out_val);
        }

	# If range is not specified, find it based on the data
        if(is.na(plot_min)){
                plot_min=min(mat);
        }
        if(is.na(plot_max)){
                plot_max=max(mat);
        }
        cat("Plot min/max: ", plot_min, "/", plot_max, "\n");

	# Get Label lengths
	row_max_nchar=max(nchar(row_names));
	col_max_nchar=max(nchar(col_names));
	cat("Max Row Names Length: ", row_max_nchar, "\n");
	cat("Max Col Names Length: ", col_max_nchar, "\n");

	##################################################################################################
	
	get_dendrogram=function(in_mat, type){
		if(type=="row"){
			dendist=dist(in_mat);
		}else{
			dendist=dist(t(in_mat));
		}
		
		get_clstrd_leaf_names=function(den){
		# Get a list of the leaf names, from left to right
			den_info=attributes(den);
			if(!is.null(den_info$leaf) && den_info$leaf==T){
				return(den_info$label);
			}else{
				lf_names=character();
				for(i in 1:2){
					lf_names=c(lf_names, get_clstrd_leaf_names(den[[i]]));
				}
				return(lf_names);
			}
		}

		hcl=hclust(dendist, method="ward.D2");
		dend=list();
		dend[["tree"]]=as.dendrogram(hcl);
		dend[["names"]]=get_clstrd_leaf_names(dend[["tree"]]);
		return(dend);
	}


	##################################################################################################
	# Comput Layouts
	col_dend_height=ceiling(num_row*.1);
	row_dend_width=ceiling(num_col*.2);
	
	heatmap_height=num_row;
	heatmap_width=num_col;

	if(num_row==1){
		plot_row_dendr=F;
	}
	if(num_col==1){
		plot_col_dendr=F;
	}

	# Don't plot dendrogram if there are any NAs in the matrix
	if(any(is.na(mat))){
		plot_col_dendr=F;
		plot_row_dendr=F;
	}

	if(plot_col_dendr && plot_row_dendr){
		layoutmat=matrix(
			c(
			rep(c(rep(4, row_dend_width), rep(3, heatmap_width)), col_dend_height),
			rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height)
			), byrow=T, ncol=row_dend_width+heatmap_width);

		col_dendr=get_dendrogram(mat, type="col");
		row_dendr=get_dendrogram(mat, type="row");

		mat=mat[row_dendr[["names"]], col_dendr[["names"]], drop=F];
		
	}else if(plot_col_dendr){
		layoutmat=matrix(
			c(
			rep(rep(2, heatmap_width), col_dend_height),
			rep(rep(1, heatmap_width), heatmap_height)
			), byrow=T, ncol=heatmap_width); 

		col_dendr=get_dendrogram(mat, type="col");
		mat=mat[, col_dendr[["names"]], drop=F];
		
	}else if(plot_row_dendr){
		layoutmat=matrix(
			rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height),
			byrow=T, ncol=row_dend_width+heatmap_width);

		row_dendr=get_dendrogram(mat, type="row");
		mat=mat[row_dendr[["names"]],,drop=F];
	}else{
		layoutmat=matrix(
			rep(1, heatmap_height*heatmap_width), 
			byrow=T, ncol=heatmap_width);
	}

	#print(layoutmat);
	layout(layoutmat);

	##################################################################################################
	
	par(oma=c(col_max_nchar*.60, 0, 3, row_max_nchar*.60));
	par(mar=c(0,0,0,0));
        plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), xaxt="n", yaxt="n", bty="n", xlab="", ylab="");
	mtext(title, side=3, line=0, outer=T, font=2);

        # x-axis
        axis(side=1, at=seq(.5, num_col-.5, 1), labels=colnames(mat), las=2, line=-1.75);
        axis(side=4, at=seq(.5, num_row-.5, 1), labels=rownames(mat), las=2, line=-1.75);

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

                        if(mat[y,x]!=0 || label_zeros){
                                if(counts){
                                        text_lab=sprintf("%i", mat[y,x]);
                                }else{
                                        text_lab=sprintf(paste("%0.", deci_pts, "f", sep=""), mat[y,x]);
                                }
                                text(x-.5, y-.5, text_lab, srt=atan(num_col/num_row)/pi*180, cex=value.cex, font=2);
                        }
                }
        }

	##################################################################################################

	par(mar=c(0, 0, 0, 0));

	if(plot_row_dendr && plot_col_dendr){
		rdh=attributes(row_dendr[["tree"]])$height;
		cdh=attributes(col_dendr[["tree"]])$height;
		plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", bty="n", xlim=c(rdh, 0));
		plot(col_dendr[["tree"]], leaflab="none",xaxt="n", yaxt="n", bty="n", ylim=c(0, cdh));
		plot(0,0, type="n", bty="n", xaxt="n", yaxt="n");
		#text(0,0, "Placeholder");
	}else if(plot_row_dendr){
		rdh=attributes(row_dendr[["tree"]])$height;
		plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", bty="n", xlim=c(rdh, 0));
		#text(0,0, "Row Dendrogram");
	}else if(plot_col_dendr){
		cdh=attributes(col_dendr[["tree"]])$height;
		plot(col_dendr[["tree"]], leaflab="none", xaxt="n", yaxt="n", bty="n", ylim=c(0, cdh));
		#text(0,0, "Col Dendrogram");
	}

	par(orig.par);

}


plot_histograms=function(table){
	num_cols=ncol(table);	
	orig.par=par(no.readonly=T);

	par(mfrow=c(5,3));
	par(mar=c(2,2,2,2));
	par(oma=c(2,2,2,2));
	colname=colnames(table);
	for(i in 1:num_cols){
		vals=table[,i];
		if(mode(vals)!="numeric" || is.factor(vals)){
			vals=as.factor(vals);
			barplot(prop.table(table(vals)), main=colname[i], col="white");
		}else{
			hist(vals, main=colname[i], xlab="values", ylab="");
		}
	}

	par(orig.par);
}

compute_correlations=function(mat){
        num_col=ncol(mat);
        cor_mat=matrix(0, nrow=num_col, ncol=num_col);
        pval_mat=matrix(0, nrow=num_col, ncol=num_col);
        rownames(cor_mat)=colnames(mat);
        colnames(cor_mat)=colnames(mat);
        rownames(pval_mat)=colnames(mat);
        colnames(pval_mat)=colnames(mat);
        for(i in 1:num_col){
                for(j in 1:i){
                        v1=mat[,i];
                        v2=mat[,j];
                        notna=!(is.na(v1) | is.na(v2));
                        #cor_mat[i,j]=cor(v1[notna], v2[notna]);
                        test=cor.test(v1[notna], v2[notna]);
                        pval_mat[i,j]=test$p.value;
                        pval_mat[j,i]=test$p.value;
                        cor_mat[i,j]=test$estimate;
                        cor_mat[j,i]=test$estimate;
                }
        }
        res=list();
        res[["val"]]=cor_mat;
        res[["pval"]]=pval_mat;;
        res[["dist"]]=as.dist(1-abs(cor_mat));
        return(res);
}

mask_matrix=function(val_mat, mask_mat, mask_thres, mask_val){
        masked_matrix=val_mat;
        masked_matrix[mask_mat>mask_thres]=mask_val;
        return(masked_matrix);
}

holm_bon_sym_mat_correction=function(pval_mat){
	tri_val=pval_mat[lower.tri(pval_mat)];
	ranks=rank(tri_val);
	num_val=length(tri_val);
	adj=numeric(num_val);

	for(i in 1:num_val){
		adj[i]=min(tri_val[i]*(num_val-ranks[i]+1), 1.0);
	}

	out_mat=matrix(NA, nrow=nrow(pval_mat), ncol=ncol(pval_mat));
	out_mat[lower.tri(out_mat)]=adj;
	out_mat=as.matrix(as.dist(out_mat));
	return(out_mat);
}

##############################################################################
##############################################################################

# Load summary file table counts 
cat("Loading summary table...\n");
counts=load_summary_file(SummaryFile);

# Remove zero count samples
cat("Looking for zero-count samples to remove.\n");
tot=apply(counts, 1, sum);
nonzero=tot>0;
if(!(all(nonzero))){
	cat("WARNING: Zero count samples found:\n");
	samp_names=rownames(counts);
	print(samp_names[!nonzero]);
	cat("\n");
	counts=counts[nonzero,,drop=F];
}else{
	cat("\tGood.  None found.\n");
}

num_categories=ncol(counts);
num_samples=nrow(counts);

# Shorten cateogry names
if(ShortenCategoryNames!=""){
	cat("Shortening category names...\n");
	full_names=colnames(counts);
	splits=strsplit(full_names, ShortenCategoryNames);
	short_names=character();
	for(i in 1:length(full_names)){
		short_names[i]=tail(splits[[i]], 1);

		short_names[i]=gsub("_unclassified$", "_uncl", short_names[i]);
		short_names[i]=gsub("_group", "_grp", short_names[i]);
	}
	colnames(counts)=short_names;
	
	cat("\tNames have been shortened.\n");
}

if(Prefix!=""){
	colnames(counts)=paste(Prefix, ".", colnames(counts), sep="");
}

# Normalize
cat("Normalizing counts to proportions.\n");
pure_counts=counts;
counts=counts+.5;
normalized=normalize(counts);
#print(normalized);

if(UseRemaining){
	category_names=colnames(counts);	
	uc_cat_names=toupper(category_names);
	remaining_ix=which(uc_cat_names=="REMAINDER" | uc_cat_names=="REMAINING");
	if(length(remaining_ix)!=1){
		cat("*******************************************************\n");
		cat("*  WARNING:  Could not identify remaining column.     *\n");
		cat("*******************************************************\n");
		UseRemaining=F;
	}else{
		cat("Remaining original column: ", remaining_ix, "\n");
		# Take out "remaining" column so it doesn't end up as a top column
		normalized_remaining_col_dat=normalized[,remaining_ix, drop=F];
		normalized=normalized[,-remaining_ix];
		cat("  Remaining column isolated.\n");
	}
}else{
	cat("Assuming no categories named Remaining or Remainder.\n");
}

# Reorder by abundance
cat("Reordering summary table categories by abundance...\n");
mean_abund=apply(normalized, 2, mean);
ix=order(mean_abund, decreasing=TRUE);
normalized=normalized[,ix];
mean_abund=mean_abund[ix];

if(UseRemaining){
	cat("Placing 'remaining' column at end of sorted mean abundance.\n");
	normalized=cbind(normalized, normalized_remaining_col_dat);
	mean_abund=c(mean_abund, mean(normalized_remaining_col_dat));
}

sorted_taxa_names=colnames(normalized);

num_top_categories=NumVariables;
num_top_categories=min(c(num_top_categories, num_categories));
prop_abundance_represented=sum(mean_abund[1:num_top_categories]);

cat("\nThe top ", num_top_categories, " categories are:\n", sep="");
for(i in 1:num_top_categories){
	cat("\t", sorted_taxa_names[i], "\t[", mean_abund[i], "]\n", sep="");
}
cat("\n");

cat("Accounting for ", prop_abundance_represented, " of taxa.\n", sep="");
cat("\n");

##############################################################################

# Assign 0's to values smaller than smallest abundance across entire dataset
min_assay=min(normalized[normalized!=0]);
#cat("Lowest non-zero value: ", min_assay, "\n", sep="");
#zero_replacment=min_assay/10;
#cat("\tSubstituting 0's with: ", zero_replacment, "\n", sep="");
#normalized[normalized==0]=zero_replacment;
#cat("Renormalizing to remove zero replacement artifacts.\n");
#normalized=normalize(normalized);

##############################################################################

if(num_top_categories>= num_categories){
	num_top_categories = (num_categories-1);
	cat("Number of taxa to work on was changed to: ", num_top_categories, "\n");
}

##############################################################################

pdf_fname=paste(OutputRoot, ".top", num_top_categories, ".alr_summary.pdf", sep="");
pdf(pdf_fname, height=11, width=9.5);

coverpage_info=c(
	paste("Summary Table Filename: ", SummaryFile),
	paste("Output Filename Root:", OutputRoot),
	"",
	paste("Number of Top ALR to include: ", NumVariables),
	paste("  Look for 'Remaining' category column? ", UseRemaining),
	"",
	paste("Lowest non-zero value: ", min_assay)
	#paste("Substituting 0's with: ", zero_replacment, ", i.e. 1/10th of Lowest non-zero = ", min_assay, "/10.", sep="")
);
plot_text(coverpage_info);

##############################################################################

# Load additional categories
additional_categories=c();

cat("\n");
cat("Extracting top categories.\n");
cat_abundances=extract_top_categories(normalized, num_top_categories, additional_cat=additional_categories);
resp_alr_struct=additive_log_rato(cat_abundances);
alr_categories_val=resp_alr_struct$transformed;
alr_cat_names=colnames(alr_categories_val);

plot_text(c(
	paste("ALR Categories (Top ", num_top_categories, ")", sep=""),
	capture.output(print(alr_cat_names))
));


##############################################################################

cat("\n");
s=summary(alr_categories_val);
plot_text(c(
	"ALR Categories Summary:",
	"\n",
	capture.output(print(s))
));
cat("Plotting ALR Category Histograms:\n");
#print(alr_categories_val);
plot_histograms(alr_categories_val);

##############################################################################

num_samples=nrow(normalized);
stat_names=c("CategoryName", "MeanAbnd", "NumSamp>0", "Preval%", "Rs95NumSamp>0", "Rs95Preval%");
num_stats=length(stat_names);

prevalence_matrix=matrix("", nrow=num_top_categories, ncol=num_stats);
colnames(prevalence_matrix)=stat_names;
rownames(prevalence_matrix)=sprintf("%3i.", 1:num_top_categories);


pure_totals=apply(pure_counts, 1, sum);
pure_normalized=normalize(pure_counts);
pure_num_samples=nrow(pure_counts);
pure_norm_means=apply(pure_normalized, 2, mean);
order=order(pure_norm_means, decreasing=T);
pure_normalized=pure_normalized[,order,drop=F];
pure_cat_names=colnames(pure_normalized);



for(i in 1:num_top_categories){

	abnd=pure_normalized[,i];
	abv=sum(abnd>0);

	prob95=integer(pure_num_samples);
	for(smp in 1:pure_num_samples){
		prob95[smp]=ifelse(pbinom(0, pure_totals[smp], pure_normalized[smp,i])<.05, 1, 0);
	}
	num_prob95=sum(prob95);

	prevalence_matrix[i,]=c(
		pure_cat_names[i], 
		sprintf("%3.4f", mean(abnd)),
		sprintf("%5i", abv), 
		sprintf("%5.1f", abv/num_samples*100.0),
		sprintf("%5i", num_prob95),
		sprintf("%5.1f", num_prob95/num_samples*100.0)
	);

}

print(prevalence_matrix, quote=F);

#------------------------------------------------------------------------------

prevalence_matrix_val=matrix(0, nrow=num_top_categories, ncol=num_stats-1);
colnames(prevalence_matrix_val)=stat_names[2:num_stats];
rownames(prevalence_matrix_val)=pure_cat_names[1:num_top_categories];

for(cix in 1:(num_stats-1)){
	prevalence_matrix_val[,cix]=as.numeric(prevalence_matrix[,cix+1]);
}

plot_prevalence_barplots=function(prev_mat, numsamp){
	orig_par=par(no.readonly=T);

	par(mfrow=c(3, 1));
	par(mar=c(10,5,2,5));

	cat_names=rownames(prev_mat);
	num_cat=length(cat_names);
	prop_rep=sum(prev_mat[,"MeanAbnd"]);
	annot=c(0,.25,.33,.5,.66,.75,1.);

	mids=barplot(prevalence_matrix_val[,"MeanAbnd"], xaxt="n", ylab="Mean Abundance", 
		col="lightblue");

	title(main=paste("Proportion of Abundance Represented in Top ", num_cat, ":", sep="")); 
	title(main=sprintf("%3.2f %%", prop_rep*100), line=-1);
		

	# Compute and label categories beneath barplot
        bar_width=mids[2]-mids[1];
        plot_range=par()$usr;
        plot_height=plot_range[4];
        label_size=min(c(1,.7*bar_width/par()$cxy[1]));

	text(mids, par()$cxy[2]/2, 1:num_top_categories, cex=label_size*.75, font=2);
        text(mids-par()$cxy[1]/2, rep(-par()$cxy[2]/2, num_top_categories), 
		cat_names, srt=-45, xpd=T, pos=4, cex=label_size);

	barplot(prevalence_matrix_val[,"NumSamp>0"], xaxt="n", ylim=c(0, numsamp),
		ylab="Num Samples w/ Taxon", col="green");
	axis(side=4, labels=annot, at=annot*numsamp, las=2, cex.axis=.7);
	abline(h=annot*numsamp, col="black", lwd=.25, lty="dashed");
	mtext("Proportion of Samples", side=4, line=4, cex=.7);

	text(mids, par()$cxy[2]/2, 1:num_top_categories, cex=label_size*.75, font=2);
        text(mids-par()$cxy[1]/2, rep(-par()$cxy[2]/2, num_top_categories), 
		cat_names, srt=-45, xpd=T, pos=4, cex=label_size);

	barplot(prevalence_matrix_val[,"Rs95NumSamp>0"], xaxt="n", ylim=c(0, numsamp),
		ylab="Num Samples w/ Taxon\n[Resampled]", col="yellow");
	axis(side=4, labels=annot, at=annot*numsamp, las=2, cex.axis=.7);
	abline(h=annot*numsamp, col="black", lwd=.25, lty="dashed");
	mtext("Proportion of Samples", side=4, line=4, cex=.7);
	text(mids, par()$cxy[2]/2, 1:num_top_categories, cex=label_size*.75, font=2);
        text(mids-par()$cxy[1]/2, rep(-par()$cxy[2]/2, num_top_categories), 
		cat_names, srt=-45, xpd=T, pos=4, cex=label_size);

	
	par(orig_par);	
}

plot_prevalence_barplots(prevalence_matrix_val, num_samples);

#------------------------------------------------------------------------------

plot_text(c(
	"The next page provides the following statistics for each category in a table:",	
	"", 
	"CategoryName:    Name of the category",
	"MeanAbnd:        Mean abundance/proportion",
	"NumSamp>0:       Number of samples with abundances > 0.0",
	"Preval%:         Percent of samples with abundances > 0.0",
	"Rs95NumSamp>0:   Number of samples with abundances > 0.0, 95% of the time if resampled", 
	"Rs95Preval%>0:   Percent of samples with abundances > 0.0, 95% of the time if resampled", 
	"",
	"Note that the Rsmp95 statistics assume the binomial distribution.",
	"For example, if a sample has a read depth of 3000, and the abundance of the category",
	"  was 0.001, then the probabilty of seeing 0 reads associated with that category",
	"  would be 0.0497.  This means that more than 1-0.0497 = 0.9503 or >95% of the time",
	"  at least 1 read will be associated with this category for this sample."
));
plot_text(c(
	paste("Total Number of Samples: ", num_samples),
	"",
	capture.output(print(prevalence_matrix, quote=F))
));

##############################################################################
# Write prevalence to file
prev_fname=paste(OutputRoot, ".top", num_top_categories, ".prevalence.csv", sep="");
cat("Writing CSV of prevalence values to file: ", prev_fname, "\n", sep="");
write.table(prevalence_matrix, file=prev_fname, quote=F, sep=",", row.names=T, col.names=NA);

# Write ALR to csv file
csv_fname=paste(OutputRoot, ".top", num_top_categories, ".csv", sep="");
cat("Writing CSV of ALR values to file: ", csv_fname, "\n", sep="");
write.table(alr_categories_val, file=csv_fname, quote=F, sep=",", row.names=T, col.names=NA);

# Write ALR to tsv file 
tsv_fname=paste(OutputRoot, ".top", num_top_categories, ".tsv", sep="");
cat("Writing TSV of ALR values to file: ", tsv_fname, "\n", sep="");
fh=file(tsv_fname, "w");
cat(file=fh, "SampleID");
close(fh);
write.table(alr_categories_val, file=tsv_fname, quote=F, sep="\t", row.names=T, col.names=NA, append=T);

# Write categories to grp file
grps_fname=paste(OutputRoot, ".top", num_top_categories, ".grps", sep="");
grp_mat=cbind(colnames(alr_categories_val), Prefix);
write.table(grp_mat, file=grps_fname, quote=F, sep="\t", row.names=F, col.names=F);

##############################################################################

cor_rec=compute_correlations(alr_categories_val);

print(cor_rec);

par(oma=c(1,1,1,1));

paint_matrix(cor_rec$val, title="Correlation: Coefficients", value.cex=.7, 
	high_is_hot=F, plot_min=-1, plot_max=1, deci_pts=2);

paint_matrix(cor_rec$pval, title="Correlation: P-values (Null Hypothesis, H0: correl=0)", 
	high_is_hot=F, value.cex=.7, plot_min=0, plot_max=1, deci_pts=2);


for(sig in c(.1, .05, .01, .001)){
	cor_at_sig=mask_matrix(cor_rec$val, cor_rec$pval, sig, 0);
	paint_matrix(cor_at_sig, title=paste("Signif Correl at: Individual p-value < ", sig, sep=""),
		value.cex=.7, plot_min=-1, plot_max=1, deci_pts=2, label_zeros=F, high_is_hot=F);
}

holm_correct_pval=holm_bon_sym_mat_correction(cor_rec$pval);

for(sig in c(.1, .05, .01, .001)){
	cor_at_sig=mask_matrix(cor_rec$val, holm_correct_pval, sig, 0);
	paint_matrix(cor_at_sig, title=paste("Signif Correl wth Holm-Bonferroni Corrected, p-value < ", sig, sep=""),
		value.cex=.7, plot_min=-1, plot_max=1, deci_pts=2, label_zeros=F, high_is_hot=F);
}

paint_matrix(as.matrix(cor_rec$dist), title="Distance: 1-|cor(x,y)|  0/1 is close/far", value.cex=.7, high_is_hot=F,
	plot_min=0, plot_max=1, deci_pts=2);


# Plot dendrogram

hcl=hclust(cor_rec$dist);
dend=as.dendrogram(hcl);
par(mar=c(20,4,4,1));
plot(dend, xlab="", ylab="Distance: 1-|cor(x,y)|", main="Associated Taxa");

# Plot correlations by

num_categories=length(alr_cat_names);

plot_cor_by_row=function(coef, pval, alpha=0.05){
	target_name=rownames(coef);
	compare_names=colnames(coef);

	self=which(compare_names==target_name)

	coef=coef[,-self];
	pval=pval[,-self];
	num_comp=length(coef);

	order=order(coef);
	coef=coef[order];
	pval=pval[order];
	compare_names=names(coef);

	coef_str=gsub("0\\.", ".", sprintf("%3.2f", coef));

	pval_ranks=rank(pval)
	adj_pval=pval;
	# Calculate Holm-Bonferroni adjustment
	for(i in 1:num_comp){
		adj_pval[i]=min(1, pval[i]*(num_comp-pval_ranks[i]+1));	
	}

	get_color=function(rho){
		prop=(rho+1)/2;
		val=rainbow(1, start=prop*4/6, end=prop*4/6+0.0001);
		return(val);
	}

	# If holm-bon signf, make bold
	# If uncorrect sigf, make regular
	# else make grey

	font_signf=rep(1, num_comp);
	col_signf=rep("grey50", num_comp);
	cex_signf=rep(.8, num_comp);
	box_lwd=rep(.5, num_comp);
	for(i in 1:num_comp){
		if(pval[i]<=alpha){
			col_signf[i]="black";
			box_lwd[i]=1;
			cex_signf[i]=1;
		}

		if(adj_pval[i]<=alpha){
			font_signf[i]=2;
			box_lwd[i]=2;
			cex_signf[i]=1.2;
		}
	}
	
	#print(coef_str);
	#print(pval);
	#print(adj_pval);

	signf=(adj_pval<=alpha);	
	cat(target_name, ":\n");
	print(coef[signf]);

	cat("Num Comparisions: ", num_comp, "\n");
	plot(0, type="n", xlim=c(0, num_comp), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n",
		main=target_name, bty="n");
		
	cxy=par()$cxy;
	for(i in 1:num_comp){

		# Draw color filled boxes
		rect(xleft=i-1, ybottom=0, xright=i, ytop=1, col=get_color(coef[i]), lwd=box_lwd[i]);

		# Label coeff inside of boxes
		text(i-.5, .5, coef_str[i], cex=(1.125+-0.0125*(num_comp+1)), font=font_signf[i], col=col_signf[i]);

		# Label taxa
		text(i-.5-cxy[1]/2, -cxy[2]/2, compare_names[i], xpd=T, pos=4, srt=-45,
			font=font_signf[i], col=col_signf[i], cex=cex_signf[i]);
	}
	
}


plot_text(c(
	"The following section contain the correlations for each taxa sorted by their",
	"correlation with another taxa.",
	"",
	"Red: Negatively correlated  (mnemonic: In bloody competition with)",
	"Green: Not correlated       (mnemonic: Sharing green pastures, minding their own business)",
	"Blue: Positively correlated (mnemonic: Together the sky is the limit)",
	"",
	"Bold Font / Thickly Boxed : significant association < alpha, even with Holm-Bonferroni correction",
	"Regular Font / Boxed : significant associations < alpha, uncorrected p-values",
	"Grey Font / Thinly Boxed : not significant",
	"",
	"Notes:",
	paste("  Alpha cutoff: ", Alpha),
	"",
	"  The number of tests assumed for the Holm-Bonferroni correction in these strip plots is: ", 
	paste("  (Number of Categories - 1) = ", num_categories - 1),
	"",
	"  The number of tests assumed in the full matrices is the half matrix area: ",
	paste("  (Number of Categories - 1) * (Numer of Categories) / 2 = ", (num_categories - 1)*num_categories/2),
	"",
	"There are fewer assumed tests in the strips because we are assuming that not all categories",
	"will be considered equally important."
));

par(mfrow=c(5,1));
par(mar=c(12, 2, 2, 6));
for(i in 1:num_categories){
	plot_cor_by_row(cor_rec$val[i,,drop=F], cor_rec$pval[i,,drop=F], Alpha);
	cat("\n\n");
}



##############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
