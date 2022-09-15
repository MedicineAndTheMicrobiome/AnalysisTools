#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');

options(useFancyQuotes=F);


params=c(
	"hypothetical_ref_summary_file", "y", 1, "character",
	"historical_ref_summary_file", "h", 1, "character",
	"qry_summary_file", "q", 1, "character",
	"output_root", "o", 1, "character",
	"shorten_category_names", "x", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-y <Hypothetical Reference summary file table>\n",
	"	-h <Historical Reference summary file table>\n",
	"	-q <Query summary file table>]\n",
	"	-o <Output filename root>\n",
	"\n",
	"	[-x <shorten category names, with separator in double quotes (default=\"\")>]\n",
	"\n",
	"This script will use the hypothatical reference summary table (1-sample) as a reference for\n",
	"what categories to focus on, as well as an 'ideal' composition.\n",
	"The historical reference a should a be a relatively large set of historical data\n",
	"so that actual performance can be compared against.\n",
	"\n");


if(
	!length(opt$hypothetical_ref_summary_file) || 
	!length(opt$historical_ref_summary_file) || 
	!length(opt$qry_summary_file) || 
	!length(opt$output_root) 
){
	cat(usage);
	q(status=-1);
}


if(length(opt$shorten_category_names)){
	ShortenCategoryNames=opt$shorten_category_names;
}else{
	ShortenCategoryNames="";
}

HypotheticalRefSummaryFile=opt$hypothetical_ref_summary_file;
HistoricalRefSummaryFile=opt$historical_ref_summary_file;
QrySummaryFile=opt$qry_summary_file;
OutputRoot=opt$output_root;

param_msg=capture.output({
	cat("\n");
	cat("Hypothetical Reference Summary File: ", HypotheticalRefSummaryFile, "\n", sep="");
	cat("Historical Reference Summary File: ", HistoricalRefSummaryFile, "\n", sep="");
	cat("Query     Summary File: ", QrySummaryFile, "\n", sep="");
	cat("Output File: ", OutputRoot, "\n", sep="");
	cat("\n");
	cat("Shorten Category Names: '", ShortenCategoryNames, "'\n", sep="");
	cat("\n");
});
print(param_msg, quote=F);

if(ShortenCategoryNames==TRUE){
	cat("Error:  You need to specify a delimitor to split the category names.\n");
	cat("        i.e., this -x option is not a flag, it requires a parameter.\n");
	quit(status=-1);
}

options(width=100);
cat("Text Line Width: ", options()$width, "\n", sep="");

##############################################################################
##############################################################################

load_summary_file=function(fname){
	inmat=as.matrix(
		read.table(fname, sep="\t", header=TRUE, check.names=FALSE, 
			comment.char="", quote="", row.names=1))

	counts_mat=inmat[,2:(ncol(inmat)), drop=F];

	# Clean category names a little
	cat_names=colnames(counts_mat);
	cat_names=gsub("-", "_", cat_names);
	colnames(counts_mat)=cat_names;
	
	cat("Num Categories in Summary Table: ", ncol(counts_mat), "\n", sep="");
	return(counts_mat);
}

merge_summary_tables=function(st1, st2){
	st1_cat_names=colnames(st1);
	st2_cat_names=colnames(st2);
	st1_samp=rownames(st1);
	st2_samp=rownames(st2);

	samp_names=sort(c(st1_samp, st2_samp));
	num_samp=length(samp_names);

	cat_names=sort(unique(c(st1_cat_names, st2_cat_names)));
	num_cat=length(cat_names);
	
	# Allocate
	merged_st=matrix(0, nrow=num_samp, ncol=num_cat);
	colnames(merged_st)=cat_names;
	rownames(merged_st)=samp_names;

	# Copy over
	merged_st[st1_samp, st1_cat_names]=st1[st1_samp, st1_cat_names];
	merged_st[st2_samp, st2_cat_names]=st2[st2_samp, st2_cat_names];
	return(merged_st);
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
}

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, deci_pts=4, 
	label_zeros=T, counts=F, value.cex=1, 
	plot_col_dendr=F,
	plot_row_dendr=F
){

	cat("Working on: ", title, "\n");

        num_row=nrow(mat);
        num_col=ncol(mat);

	if(num_row==0 || num_col==0){
		cat("Nothing to plot.\n");
		return();
	}

	any_nas=any(is.na(mat));

	if(num_row==1 || any_nas){
		plot_row_dendr=F;
	}
	if(num_col==1 || any_nas){
		plot_col_dendr=F;
	}

	row_names=rownames(mat);
	col_names=colnames(mat);

	orig.par=par(no.readonly=T);

        cat("Num Rows: ", num_row, "\n");
        cat("Num Cols: ", num_col, "\n");

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
                plot_min=min(mat, na.rm=T);
        }
        if(is.na(plot_max)){
                plot_max=max(mat, na.rm=T);
        }

	if(plot_min>=-1 && plot_max<=1){
		fractions_only=T;	
	}else{
		fractions_only=F;
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
		mat=mat[row_dendr[["names"]],, drop=F];
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

                        if(is.na(mat[y,x]) || mat[y,x]!=0 || label_zeros){
                                if(counts){
                                        text_lab=sprintf("%i", mat[y,x]);
                                }else{
                                        text_lab=sprintf(paste("%0.", deci_pts, "f", sep=""), mat[y,x]);
					if(fractions_only){
						if(!is.na(mat[y,x])){
							if(mat[y,x]==-1 || mat[y,x]==1){
								text_lab=as.integer(mat[y,x]);	
							}else{
								text_lab=gsub("0\\.","\\.", text_lab);
							}
						}
					}
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

##############################################################################
##############################################################################

# Open main output file
pdf(paste(OutputRoot, ".ref_cmp.pdf", sep=""), height=8.5, width=14);

plot_text(param_msg);

# Load summary file table counts 
cat("\n");
cat("Loading summary table...\n");
hypo_ref_st=load_summary_file(HypotheticalRefSummaryFile);
hist_ref_st=load_summary_file(HistoricalRefSummaryFile);
query_st=load_summary_file(QrySummaryFile);

dim_msg=capture.output({
	cat("\n");
	cat("Hypothetical Reference Dimensions:\n", dim(hypo_ref_st), "\n");
	cat("Historical Reference Dimensions:\n", dim(hist_ref_st), "\n");
	cat("Query Dimensions:\n", dim(query_st), "\n");
}); 
print(dim_msg, quote=F);
plot_text(dim_msg);

hypo_ref_cat=colnames(hypo_ref_st);
hist_ref_cat=colnames(hist_ref_st);
qry_cat=colnames(query_st);

counts=merge_summary_tables(hypo_ref_st, hist_ref_st);
counts=merge_summary_tables(counts, query_st);

hypo_ref_samp_ids=rownames(hypo_ref_st);
hist_ref_samp_ids=rownames(hist_ref_st);
qry_samp_ids=rownames(query_st);
sample_ids=rownames(counts);
num_samples=nrow(counts);

ids_msg=capture.output({
	cat("\n\n");

	cat("Hypothetical Reference Sample ID:\n");
	print(hypo_ref_samp_ids);
	cat("\n");

	cat("Historical Reference Sample ID:\n");
	print(hist_ref_samp_ids);
	cat("\n");

	cat("Query Sample IDs:\n");
	print(qry_samp_ids);

	cat("\n\n");
});
print(ids_msg, quote=F);
plot_text(ids_msg);

##############################################################################

# Shorten cateogry names
if(ShortenCategoryNames!=""){
	full_names=colnames(counts);
	splits=strsplit(full_names, ShortenCategoryNames);
	short_names=character();
	for(i in 1:length(full_names)){
		short_names[i]=tail(splits[[i]], 1);
		short_names[i]=gsub("_unclassified$", "_uncl", short_names[i]);
		short_names[i]=gsub("\\[", "", short_names[i]);
		short_names[i]=gsub("\\]", "", short_names[i]);
	}
	colnames(counts)=short_names;
	cat("Names have been shortened.\n");
}else{
	cat("Keeping original category names...\n");
}

# Remove 0 count categories
num_categories_merged=ncol(counts);
cat("Number of Categories in merged (hyp ref + hist ref + qry) summary tables: ", num_categories_merged, "\n");
total_counts=apply(counts, 2, sum);

cat("Removing categories if they are 0's across all samples.\n");
counts=counts[,total_counts>0];
num_categories_nz=ncol(counts);
cat("Number of Categories with non-zero counts: ", num_categories_nz, "\n");

cat_names=colnames(counts);

# Collapse everything down to what is in reference and other
hypo_ref_counts=counts[hypo_ref_samp_ids,, drop=F];
in_hypo_ref_cat_ix=hypo_ref_counts>0;

in_hypo_ref_cat=cat_names[in_hypo_ref_cat_ix];
rem_cat=setdiff(cat_names, in_hypo_ref_cat);

cat_mesg=capture.output({
	cat("\nHypothetical Reference Categories:\n");
	print(in_hypo_ref_cat);

	cat("\nNon-Reference (remaining) Categories:\n");
	print(rem_cat);
});
print(cat_mesg, quote=F);
plot_text(cat_mesg);

###############################################################################

# Normalize
cat("Normalizing counts...\n");
#counts=counts+.05; # Use smaller adjustment to maintain values similar relative abundance
normalized=normalize(counts);

# Add centroids for hist and qry
hist_centroid=apply(normalized[hist_ref_samp_ids,,drop=F], 2, mean);
qry_centroid=apply(normalized[qry_samp_ids,,drop=F], 2, mean);
normalized=rbind(normalized, hist_centroid, qry_centroid);
#sample_ids=rownames(normalized);
sample_ids=rownames(normalized);


# Reorder by abundance
cat("Reordering summary table categories by abundance...\n");
mean_abund=apply(normalized, 2, mean);
ix=order(mean_abund, decreasing=TRUE);
normalized=normalized[,ix];
mean_abund=mean_abund[ix];

# Extract ref cat only 
refcat_only_mat=normalized[,in_hypo_ref_cat,drop=F];
refcat_only_renorm_mat=normalize(refcat_only_mat);

cat("Reference Categories only:\n");
remaining_abd=apply(refcat_only_mat, 1, function(x){1-sum(x)});
refcat_only_mat=cbind(refcat_only_mat, remaining_abd);
colnames(refcat_only_mat)=c(colnames(refcat_only_mat)[1:(ncol(refcat_only_mat)-1)], "Remaining");
print(refcat_only_mat);

cat("Reference Categories renormalized:\n");
print(refcat_only_renorm_mat);

###############################################################################

refcat_only_dist=dist(refcat_only_mat, "manhattan");
refcat_only_mds=cmdscale(refcat_only_dist);

refcat_only_renorm_dist=dist(refcat_only_renorm_mat, "manhattan");
refcat_only_renorm_mds=cmdscale(refcat_only_renorm_dist);

par(mfrow=c(1,2));
par(mar=c(2,2,3,1));

num_mds_points=num_samples+2;

colors=rep("blue", num_mds_points);
names(colors)=sample_ids;
colors[hypo_ref_samp_ids]="black";
colors[qry_samp_ids]="green";
colors["hist_centroid"]="blue";
colors["qry_centroid"]="green";

sizes=rep(1, num_mds_points);
names(sizes)=sample_ids;
sizes[hypo_ref_samp_ids]=3;
sizes[qry_samp_ids]=2;
sizes["hist_centroid"]=2.5;
sizes["qry_centroid"]=2.5;

lwds=rep(1, num_mds_points);
names(lwds)=sample_ids;
lwds[hypo_ref_samp_ids]=2;
lwds[qry_samp_ids]=3;
lwds["hist_centroid"]=1.5;
lwds["qry_centroid"]=1.5;

cexs=rep(.5, num_mds_points);
names(cexs)=sample_ids;
cexs[hypo_ref_samp_ids]=.75;
cexs[qry_samp_ids]=1;
cexs["hist_centroid"]=.7;
cexs["qry_centroid"]=.7;

pchs=rep(1, num_mds_points);
names(pchs)=sample_ids;
pchs[hypo_ref_samp_ids]=1;
pchs[qry_samp_ids]=1;
pchs["hist_centroid"]=3;
pchs["qry_centroid"]=3;

sample_ids_no_centroids=sample_ids;
names(sample_ids_no_centroids)=sample_ids;
sample_ids_no_centroids["hist_centroid"]="";
sample_ids_no_centroids["qry_centroid"]="";


calc_exp_lim=function(vals, mar){
	r=range(vals);
	span=r[2]-r[1];
	return(c(r[1]-mar*span, r[2]+mar*span));
}


# Plot points
plot(refcat_only_mds[,1], refcat_only_mds[,2], xlab="", ylab="", 
	xlim=calc_exp_lim(refcat_only_mds[,1], mar=.05),
	main="MDS Reference Categories Only", 
	col=colors, cex=sizes, lwd=lwds, pch=pchs);

plot(refcat_only_renorm_mds[,1], refcat_only_renorm_mds[,2], xlab="", ylab="",
	xlim=calc_exp_lim(refcat_only_renorm_mds[,1], mar=.05),
	main="MDS Reference Categories Only, with Renormal.", 
	col=colors, cex=sizes, lwd=lwds, pch=pchs);


# Plot with labels
plot(refcat_only_mds[,1], refcat_only_mds[,2], xlab="", ylab="", 
	xlim=calc_exp_lim(refcat_only_mds[,1], .25),
	main="MDS Reference Categories Only", 
	col="grey", cex=sizes/2);
text(refcat_only_mds[,1], refcat_only_mds[,2], sample_ids_no_centroids, cex=cexs, col=colors);

plot(refcat_only_renorm_mds[,1], refcat_only_renorm_mds[,2], xlab="", ylab="",
	xlim=calc_exp_lim(refcat_only_renorm_mds[,1], .25),
	main="MDS Reference Categories Only, with Renormal.", 
	col="grey", cex=sizes/2);
text(refcat_only_renorm_mds[,1], refcat_only_renorm_mds[,2], sample_ids_no_centroids, cex=cexs, col=colors);


###############################################################################

plot_against_reference=function(norm_mat, hypo_ids, hist_ids, qry_ids, title){

	# Sort decreasing by hypothetical abundance
	hypo_ref_abund=norm_mat[hypo_ids,,drop=F];	
	order_ix=order(hypo_ref_abund, decreasing=T);
	norm_mat=norm_mat[,order_ix];

	hypo_abund=norm_mat[hypo_ids,,drop=F];
	hist_abund=norm_mat[hist_ids,,drop=F];
	qry_abund=norm_mat[qry_ids,,drop=F];

	median_hist_abund=apply(hist_abund, 2, median);
	median_qry_abund=apply(qry_abund, 2, median);
	
	num_cat=ncol(norm_mat);
	num_hist_ids=length(hist_ids);
	num_qry_ids=length(qry_ids);

	max_abund=max(norm_mat);

	par(mar=c(13, 3, 4, 1));
	mids=barplot(hypo_abund,
		ylim=c(0, max_abund*1.2), col="white",
		las=2,
		main=title
	);

	scat_var=(mids[2]-mids[1])/16;
	scat=rnorm(num_qry_ids+num_hist_ids, 0, scat_var);
	names(scat)=c(qry_ids, hist_ids);

	for(catix in 1:num_cat){
		for(qix in hist_ids){
			points(mids[catix]+scat[qix], norm_mat[qix, catix], col="blue");
		}

		for(qix in qry_ids){
			points(mids[catix]+scat[qix], norm_mat[qix, catix], col="green", lwd=1.1);
		}

		points(mids[catix], median_hist_abund[catix], pch="-", cex=3, col="blue");
		points(mids[catix], median_qry_abund[catix], pch="-", cex=3, col="green", lwd=1.1);
	}

}

plot_against_reference(refcat_only_mat, 
	hypo_ref_samp_ids, hist_ref_samp_ids, qry_samp_ids,
	"Reference Categories Only");

plot_against_reference(refcat_only_renorm_mat,
	hypo_ref_samp_ids, hist_ref_samp_ids, qry_samp_ids,
	"Reference Categories Only, with Renormal.");

###############################################################################

plot_overlapping_histograms=function(norm_mat, hypo_ids, hist_ids, qry_ids){
	
	# Sort decreasing by hypothetical abundance
        hypo_ref_abund=norm_mat[hypo_ids,,drop=F];
        order_ix=order(hypo_ref_abund, decreasing=T);
        norm_mat=norm_mat[,order_ix];

        hypo_abund=norm_mat[hypo_ids,,drop=F];
        hist_abund=norm_mat[hist_ids,,drop=F];
        qry_abund=norm_mat[qry_ids,,drop=F];

        median_hist_abund=apply(hist_abund, 2, median);
        median_qry_abund=apply(qry_abund, 2, median);

        num_cat=ncol(norm_mat);
        num_hist_ids=length(hist_ids);
        num_qry_ids=length(qry_ids);

	par(mfrow=c(3,4));
        par(mar=c(3,3,4,1));

	catnames=colnames(norm_mat);
	max_abund=max(norm_mat);

        for(catix in 1:num_cat){

		hypo_val=norm_mat[hypo_ids,catix];
		hist_val=norm_mat[hist_ids,catix];
		qry_val=norm_mat[qry_ids,catix];

		test.res=wilcox.test(hist_val, qry_val);
		pval=test.res$p.value;

		if(pval<0.05){
			sigcol="red";
		}else{
			sigcol="darkgreen";
		}

		h=hist(hist_val, main=catnames[catix], 
			xlab="Relative Abundance", ylab="Frequency", 
			col="blue",
			breaks=seq(0, max_abund, length.out=20),
		);

		title(main=paste("p-value: ", round(pval,4), sep=""),
			line=.9, cex.main=.85, col.main=sigcol);

		abline(v=hypo_val, col="black", lwd=3);
		abline(v=qry_val, col="green", lwd=2);
		abline(v=qry_val, col="black", lwd=.5);


        }

}

plot_overlapping_histograms(refcat_only_mat, hypo_ref_samp_ids, hist_ref_samp_ids, qry_samp_ids);

###############################################################################

generate_qry_ref_stats_table=function(norm_mat, hypo_ids, hist_ids, qry_ids){
	
	# Sort decreasing by hypothetical abundance
        hypo_ref_abund=norm_mat[hypo_ids,,drop=F];
        order_ix=order(hypo_ref_abund, decreasing=T);
        norm_mat=norm_mat[,order_ix];

        hypo_abund=norm_mat[hypo_ids,,drop=F];
        hist_abund=norm_mat[hist_ids,,drop=F];
        qry_abund=norm_mat[qry_ids,,drop=F];

        median_hist_abund=apply(hist_abund, 2, median);
        median_qry_abund=apply(qry_abund, 2, median);

        num_cat=ncol(norm_mat);
        num_hist_ids=length(hist_ids);
        num_qry_ids=length(qry_ids);

	catnames=colnames(norm_mat);

	info_hdr=c(
		"Category",
		"Hypo_Abd",
		"Hist_Abd",
		"Qry_Abd",
		"Hist_Stdev",
		"Qry_Stdev",
		"HypoQry_Diff",
		"HistQry_Diff",
		"HistQry_Pval",
		"Signf"
	);

	tab=matrix("", nrow=num_cat, ncol=length(info_hdr));
	colnames(tab)=info_hdr

	signed_str=function(x){
		if(x>0){
			s="+";
		}else if(x<0){
			s="-";
		}else{
			s=" ";
		}

		return(
			paste(s, round(abs(x), 4), sep="")
		);
	}

        for(catix in 1:num_cat){

		hypo_val=norm_mat[hypo_ids,catix];
		hist_val=norm_mat[hist_ids,catix];
		qry_val=norm_mat[qry_ids,catix];

		test.res=wilcox.test(hist_val, qry_val);
		pval=test.res$p.value;

		tab[catix, "Category"]=catnames[catix];
		tab[catix, "Hypo_Abd"]=round(hypo_val, 4);
		tab[catix, "Hist_Abd"]=round(mean(hist_val), 4);
		tab[catix, "Qry_Abd"]=round(mean(qry_val), 4);
		tab[catix, "Hist_Stdev"]=round(sd(hist_val), 4);
		tab[catix, "Qry_Stdev"]=round(sd(qry_val), 4);

		tab[catix, "HypoQry_Diff"]=signed_str(mean(qry_val)-hypo_val);
		tab[catix, "HistQry_Diff"]=signed_str(mean(qry_val)-mean(hist_val));
		tab[catix, "HistQry_Pval"]=round(pval, 4);

		s="";
		if(pval<0.0001){
			s="***";
		}else if(pval<0.001){
			s="**";
		}else if(pval<0.01){
			s="*";
		}else if(pval<0.05){
			s="-";
		}else if(pval<0.1){
			s=".";	
		}

		tab[catix, "Signf"]=s;
        }
	
	return(tab);

}

table=generate_qry_ref_stats_table(refcat_only_mat, hypo_ref_samp_ids, hist_ref_samp_ids, qry_samp_ids);

options(width=200);
tab=capture.output({
	print(table, quote=F, width=200);
});

print(tab);

par(mfrow=c(1,1));
plot_text(
	tab
);

###############################################################################

# Extract only "remaining" categories not in Hypothetical
all_nz_cat=colnames(normalized);
all_nz_remaining_cat=setdiff(all_nz_cat, in_hypo_ref_cat);
remaining_normalized=normalized[,all_nz_remaining_cat];

# Sort by query 
qry_remaining_abd=remaining_normalized[qry_samp_ids,,drop=F];
mean_qry_remaining_abd=apply(qry_remaining_abd, 2, mean);
qry_rem_sort_ix=order(mean_qry_remaining_abd, decreasing=T, method="shell");
qry_remain_ordered_abd=remaining_normalized[,qry_rem_sort_ix,drop=F];

# Sort by hist ref
hist_remaining_abd=remaining_normalized[hist_ref_samp_ids,,drop=F];
mean_hist_remaining_abd=apply(hist_remaining_abd, 2, mean);
hist_rem_sort_ix=order(mean_hist_remaining_abd, decreasing=T, method="shell");
hist_remain_ordered_abd=remaining_normalized[,hist_rem_sort_ix,drop=F];

# Sort by combined
mean_comb_remaining_abd=(mean_qry_remaining_abd+mean_hist_remaining_abd)/2;
comb_rem_sort_ix=order(mean_comb_remaining_abd, decreasing=T, method="shell");
comb_remain_ordered_abd=remaining_normalized[,comb_rem_sort_ix,drop=F];


plot_qry_vs_ref=function(norm_mat, hist_ids, qry_ids, num_top, title){

	num_cat=ncol(norm_mat);
	num_top=min(num_top, num_cat);
	
	norm_mat=norm_mat[,1:num_top];

	hist_abund=norm_mat[hist_ids,,drop=F];
	qry_abund=norm_mat[qry_ids,,drop=F];

	median_hist_abund=apply(hist_abund, 2, median);
	median_qry_abund=apply(qry_abund, 2, median);
	
	num_hist_ids=length(hist_ids);
	num_qry_ids=length(qry_ids);

	max_abund=max(norm_mat);

	par(mar=c(13, 3, 4, 1));
	placeholder=rep(0, num_top);
	names(placeholder)=colnames(norm_mat);
	mids=barplot(placeholder,
		ylim=c(0, max_abund*1.2), col="white",
		las=2,
		main=title
	);

	scat_var=(mids[2]-mids[1])/16;
	scat=rnorm(num_qry_ids+num_hist_ids, 0, scat_var);
	names(scat)=c(qry_ids, hist_ids);

	for(catix in 1:num_top){
		for(qix in hist_ids){
			if(norm_mat[qix, catix]>0){
				points(mids[catix]+scat[qix], norm_mat[qix, catix], col="blue");
			}
		}

		for(qix in qry_ids){
			if(norm_mat[qix, catix]>0){
				points(mids[catix]+scat[qix], norm_mat[qix, catix], col="green", lwd=1.1);
			}
		}

		points(mids[catix], median_hist_abund[catix], pch="-", cex=3, col="blue");
		points(mids[catix], median_qry_abund[catix], pch="-", cex=3, col="green", lwd=1.1);
	}

}

plot_qry_vs_ref(comb_remain_ordered_abd, hist_ref_samp_ids, qry_samp_ids, num_top=5, "Top 5 Combined Remainders");
plot_qry_vs_ref(comb_remain_ordered_abd, hist_ref_samp_ids, qry_samp_ids, num_top=15, "Top 15 Combined Remainders");
plot_qry_vs_ref(comb_remain_ordered_abd, hist_ref_samp_ids, qry_samp_ids, num_top=50, "Top 50 Combined Remainders");

plot_qry_vs_ref(qry_remain_ordered_abd, hist_ref_samp_ids, qry_samp_ids, num_top=5, "Top 5 Query Remainders");
plot_qry_vs_ref(qry_remain_ordered_abd, hist_ref_samp_ids, qry_samp_ids, num_top=15, "Top 15 Query Remainder");
plot_qry_vs_ref(qry_remain_ordered_abd, hist_ref_samp_ids, qry_samp_ids, num_top=50, "Top 50 Query Remainder");

plot_qry_vs_ref(hist_remain_ordered_abd, hist_ref_samp_ids, qry_samp_ids, num_top=5, "Top 5 Historical Remainders");
plot_qry_vs_ref(hist_remain_ordered_abd, hist_ref_samp_ids, qry_samp_ids, num_top=15, "Top 15 Historical Remainders");
plot_qry_vs_ref(hist_remain_ordered_abd, hist_ref_samp_ids, qry_samp_ids, num_top=50, "Top 50 Historical Remainders");


###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
