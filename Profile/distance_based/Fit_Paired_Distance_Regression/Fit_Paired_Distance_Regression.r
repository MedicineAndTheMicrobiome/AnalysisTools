#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
library(car);

source('~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r');

options(useFancyQuotes=F);

DEF_DISTTYPE="euc";

params=c(
	"summary_file", "s", 1, "character",

	"factors", "f", 1, "character",
	"factor_samp_id_name", "F", 1, "character",
	"model_var", "M", 1, "character",
	"required", "q", 2, "character",

	"pairings", "p", 1, "character",
	"B_minuend", "B", 1, "character",
	"A_subtrahend", "A", 1, "character",

	"dist_type", "d", 2, "character",
	"outputroot", "o", 2, "character",

	"reference_levels", "c", 2, "character",
	"shorten_category_names", "x", 2, "character",

	"tag_name", "t", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

NUM_TOP_RESP_CAT=35;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table>\n",
	"\n",
	"	-f <factors file, contains covariates and factors>\n",
	"	-F <column name of sample ids in factor file>\n",
	"	-M <list of covariate X's names to include in the model from the factor file>\n",
	"	[-q <required list of variables to include after NA removal>]\n",
	"\n",
	"	-p <pairings map, pairing sample IDs from two groups. Must have header/column names>\n",
	"	-B <Sample Group B, (column name in pairings file)\n",
	"	-A <Sample Group A, (column name in pairings file)\n",
	"\n",
	"	[-d <euc/wrd/man/bray/horn/bin/gow/tyc/minkp5/minkp3, default =", DEF_DISTTYPE, ">]\n",
	"\n",
	"	[-o <output filename root>]\n",
	"\n",
	"	[-c <reference levels file for Y's in factor file>]\n",
        "       [-x <shorten category names, with separator in double quotes (default=\"\")>]\n",
	"\n",
	"	[-t <tag name>]\n",
	"\n",
	"This script will fit the following model for diversity:\n",
	"\n",
	"	1.) Distance(A[i], B[i]) = covariates\n",
	"\n",
	"The intercept of the model represents the difference between B and A.\n",
	"\n", sep="");

if(
	!length(opt$summary_file) || 
	!length(opt$factors) || 
	!length(opt$model_var) || 
	!length(opt$A_subtrahend) || 
	!length(opt$B_minuend) || 
	!length(opt$pairings) ||
	!length(opt$factor_samp_id_name)
){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputRoot=gsub(".summary_table.xls", "", opt$summary_file);
	OutputRoot=gsub(".summary_table.tsv", "", OutputRoot);
}else{
	OutputRoot=opt$outputroot;
}

if(!length(opt$reference_levels)){
        ReferenceLevelsFile="";
}else{
        ReferenceLevelsFile=opt$reference_levels;
}

if(length(opt$required)){
	RequiredFile=opt$required;
}else{
	RequiredFile="";
}

if(length(opt$dist_type)){
	DistType=opt$dist_type;
}else{
	DistType=DEF_DISTTYPE;

}

if(length(opt$shorten_category_names)){
        ShortenCategoryNames=opt$shorten_category_names;
}else{
        ShortenCategoryNames="";
}

if(length(opt$factor_samp_id_name)){
	FactorSampleIDName=opt$factor_samp_id_name;
}else{
	FactorSampleIDName=1;
}

if(length(opt$tag_name)){
        TagName=opt$tag_name;
        cat("Setting TagName Hook: ", TagName, "\n");
        setHook("plot.new",
                function(){
                        #cat("Hook called.\n");
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

SummaryFile=opt$summary_file;
FactorsFile=opt$factors;
ModelVarFile=opt$model_var;
PairingsFile=opt$pairings;
A_subtrahend=opt$A_subtrahend;
B_minuend=opt$B_minuend;

FactorSampleIDName=opt$factor_samp_id_name;

OutputRoot=paste(OutputRoot, ".a_", A_subtrahend, ".b_", B_minuend, sep="");

cat("\n");
cat("         Summary File: ", SummaryFile, "\n", sep="");
cat("         Factors File: ", FactorsFile, "\n", sep="");
cat("Factor Sample ID Name: ", FactorSampleIDName, "\n", sep="");
cat(" Model Variables File: ", ModelVarFile, "\n", sep="");
cat("        Pairings File: ", PairingsFile, "\n", sep="");
cat("                   A : ", A_subtrahend, "\n", sep="");
cat("                   B : ", B_minuend, "\n", sep="");
cat("          Output File: ", OutputRoot, "\n", sep="");
cat("\n");
cat("        Distance Type: ", DistType, "\n", sep="");
cat("\n");
cat("Reference Levels File: ", ReferenceLevelsFile, "\n", sep="");
cat("\n");

options(width=100);
cat("Text Line Width: ", options()$width, "\n", sep="");

##############################################################################
##############################################################################

mask_matrix=function(val_mat, mask_mat, mask_thres, mask_val){
	masked_matrix=val_mat;
	masked_matrix[mask_mat>mask_thres]=mask_val;
	return(masked_matrix);
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
        axis(side=1, at=seq(.5, num_col-.5, 1), labels=colnames(mat), las=2, line=-1.75, cex.axis=value.cex);
        axis(side=4, at=seq(.5, num_row-.5, 1), labels=rownames(mat), las=2, line=-1.75, cex.axis=value.cex);

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

add_sign_col=function(coeff){
	cnames=colnames(coeff);
	pval_ix=which(cnames=="Pr(>|t|)");
	pval=coeff[,pval_ix];

	sig_char=function(val){
		if(val == 0){ return("***");}
		if(val <= .001){ return("** ");}
		if(val <= .01){ return("*  ");}
		if(val <= .05){ return(".  ");}
		return(" ");
	}

	sig_arr=sapply(pval, sig_char);
	fdr=round(p.adjust(pval, method="fdr"), 4);
	fdr_sig_char=sapply(fdr, sig_char);
	#out_mat=cbind(coeff, fdr);
	out_mat=cbind(coeff, sig_arr, fdr, fdr_sig_char);
	colnames(out_mat)=c(cnames, "Signf", "FDR", "Signf");
	
	return(formatC(out_mat, format="f", digits=5,));
	
}

mask_matrix=function(val_mat, mask_mat, mask_thres, mask_val){
        masked_matrix=val_mat;
        masked_matrix[mask_mat>mask_thres]=mask_val;
        return(masked_matrix);
}

compute_dist=function(norm_st, type){

        if(type=="euc"){
                dist_mat=dist(norm_st);
        }else if (type=="wrd"){
                dist_mat=weight_rank_dist_opt(norm_st, deg=4);
        }else if (type=="man"){
                dist_mat=vegdist(norm_st, method="manhattan");
        }else if (type=="bray"){
                dist_mat=vegdist(norm_st, method="bray");
        }else if (type=="horn"){
                dist_mat=vegdist(norm_st, method="horn");
        }else if (type=="bin"){
                dist_mat=vegdist(norm_st, method="bin");
        }else if (type=="gow"){
                dist_mat=vegdist(norm_st, method="gower");
        }else if (type=="tyc"){
                dist_mat=thetaYC(norm_st);
        }else if (type=="minkp3"){
                dist_mat=dist(norm_st, method="minkowski", p=1/3);
        }else if (type=="minkp5"){
                dist_mat=dist(norm_st, method="minkowski", p=1/2);
        }

        dist_mat[dist_mat==0]=1e-323;

        return(dist_mat);
}


##############################################################################
##############################################################################

# Open main output file
pdf(paste(OutputRoot, ".paird_dist.", DistType, ".pdf", sep=""), height=11, width=9.5);

input_files=load_and_reconcile_files(
                sumtab=list(fn=SummaryFile, shorten_cat_names_char=ShortenCategoryNames,
                        return_top=NULL, specific_cat_fn=NULL),
                factors=list(fn=FactorsFile),
                pairs=list(fn=PairingsFile, a_cname=A_subtrahend, b_cname=B_minuend),
                covariates=list(fn=ModelVarFile),
                grpvar=list(fn=""),
                reqvar=list(fn=RequiredFile)
        );

counts=input_files[["SummaryTable_counts"]];
factors=input_files[["Factors"]];
good_pairs_map=input_files[["PairsMap"]];
model_var_arr=input_files[["Covariates"]];
required_arr=input_files[["RequiredVariables"]];
normalized=input_files[["SummaryTable_normalized"]];


##############################################################################
##############################################################################

# Order/Split the data....
# Order the pairings map by response sample IDs
A_sample_ids=good_pairs_map[,A_subtrahend];
B_sample_ids=good_pairs_map[,B_minuend];
AB_sample_ids=c(A_sample_ids, B_sample_ids);

##############################################################################
# Compute diversity indices
cat("Computing distance matrix:\n");

distmat=as.matrix(compute_dist(normalized, DistType));

for(i in 1:length(distmat)){
        if(distmat[i]==0){
                distmat[i]=1e-323;
        }
}

paired_dist=diag(distmat[A_sample_ids,B_sample_ids]);
names(paired_dist)=A_sample_ids;
#print(paired_dist);

mds=cmdscale(distmat);
#print(mds);

cat("A Sample IDs:\n");
print(A_sample_ids);

cat("B Sample IDs:\n");
print(B_sample_ids);

A_centroid=c(mean(mds[A_sample_ids,1]), mean(mds[A_sample_ids,2]));
B_centroid=c(mean(mds[B_sample_ids,1]), mean(mds[B_sample_ids,2]));

cat("A centroid: \n");
print(A_centroid);
cat("B centroid: \n");
print(B_centroid);

num_shared_sample_ids=length(AB_sample_ids);
colors=numeric(num_shared_sample_ids)
names(colors)=AB_sample_ids;
colors[A_sample_ids]="green";
colors[B_sample_ids]="blue";

glyphs=character(num_shared_sample_ids);
names(glyphs)=AB_sample_ids;

# Find a short but unique abbreviation
shorter_name_len=min(c(nchar(A_subtrahend), nchar(B_minuend)));
abbrev_len=min(shorter_name_len, 5);
while(substr(A_subtrahend, 1, abbrev_len)==substr(B_minuend, 1, abbrev_len) && abbrev_len<=shorter_name_len){
	abbrev_len=abbrev_len+1;	
}
glyphs[A_sample_ids]=substr(A_subtrahend, 1, abbrev_len);
glyphs[B_sample_ids]=substr(B_minuend, 1, abbrev_len);

remap_val=function(inval, end_rng){
	in_rng=range(inval);
	norm=(inval-in_rng[1])/(in_rng[2]-in_rng[1]);
	outval=norm*(end_rng[2]-end_rng[1])+end_rng[1];
	return(outval);
}


mds_x_range=range(mds[,1]);
mds_y_range=range(mds[,2]);

mds_x_width=diff(mds_x_range);
mds_y_height=diff(mds_y_range);


draw_centroids=function(acent, bcent, acol="green", bcol="blue"){
	points(acent[1], acent[2], pch=19, col=acol, cex=4);	
	points(acent[1], acent[2], pch=19, col="white", cex=3);	

	points(bcent[1], bcent[2], pch=19, col=bcol, cex=4);	
	points(bcent[1], bcent[2], pch=19, col="white", cex=3);	
}

par(mar=c(4,4,4,4));

#------------------------------------------------------------------------------
# Plot points
plot(0,0, type="n", xlab="Dim 1", ylab="Dim 2", 
	xlim=c(mds_x_range[1]-mds_x_width*.1, mds_x_range[2]+mds_x_width*.1),
	ylim=c(mds_y_range[1]-mds_y_height*.1, mds_y_range[2]+mds_y_height*.1),
	main="Points Only");

draw_centroids(A_centroid, B_centroid);

points(mds[,1], mds[,2], col=colors)

#------------------------------------------------------------------------------
# Plot labels
plot(0,0, type="n", xlab="Dim 1", ylab="Dim 2", 
	xlim=c(mds_x_range[1]-mds_x_width*.1, mds_x_range[2]+mds_x_width*.1),
	ylim=c(mds_y_range[1]-mds_y_height*.1, mds_y_range[2]+mds_y_height*.1),
	main="Labelled by Sample Type");

draw_centroids(A_centroid, B_centroid);

text(mds[,1], mds[,2], glyphs, col=colors);

#------------------------------------------------------------------------------
# Plot with sample IDs 
plot(0,0, type="n", xlab="Dim 1", ylab="Dim 2", 
	xlim=c(mds_x_range[1]-mds_x_width*.1, mds_x_range[2]+mds_x_width*.1),
	ylim=c(mds_y_range[1]-mds_y_height*.1, mds_y_range[2]+mds_y_height*.1),
	main=paste("Labelled by ", A_subtrahend, sep=""));

draw_centroids(A_centroid, B_centroid);

label_by_A=glyphs;
label_by_A[A_sample_ids]=A_sample_ids;
text(mds[,1], mds[,2], label_by_A, col=colors, cex=.75);

# Plot with sample IDs 
plot(0,0, type="n", xlab="Dim 1", ylab="Dim 2", 
	xlim=c(mds_x_range[1]-mds_x_width*.1, mds_x_range[2]+mds_x_width*.1),
	ylim=c(mds_y_range[1]-mds_y_height*.1, mds_y_range[2]+mds_y_height*.1),
	main=paste("Labelled by ", B_minuend, sep=""));

draw_centroids(A_centroid, B_centroid);

label_by_B=glyphs;
label_by_B[B_sample_ids]=B_sample_ids;
text(mds[,1], mds[,2], label_by_B, col=colors, cex=.75);

#------------------------------------------------------------------------------
# Plot connected unweighted
plot(0,0, type="n", xlab="Dim 1", ylab="Dim 2", 
	xlim=c(mds_x_range[1]-mds_x_width*.1, mds_x_range[2]+mds_x_width*.1),
	ylim=c(mds_y_range[1]-mds_y_height*.1, mds_y_range[2]+mds_y_height*.1),
	main="Pairs Connected");

draw_centroids(A_centroid, B_centroid);

for(ix in 1:(num_shared_sample_ids/2)){
	aid=A_sample_ids[ix];
	bid=B_sample_ids[ix];
	points(
		c(mds[aid,1], mds[bid,1]),
		c(mds[aid,2], mds[bid,2]),
		type="l", lwd=.5, col="grey70");
}
text(mds[,1], mds[,2], glyphs, col=colors);

#------------------------------------------------------------------------------
# Plot connected weighted
plot(0,0, type="n", xlab="Dim 1", ylab="Dim 2", 
	xlim=c(mds_x_range[1]-mds_x_width*.1, mds_x_range[2]+mds_x_width*.1),
	ylim=c(mds_y_range[1]-mds_y_height*.1, mds_y_range[2]+mds_y_height*.1),
	main="Closer Pairs Connected with Thicker Lines");

draw_centroids(A_centroid, B_centroid);

line_wts=remap_val(paired_dist, c(4,.25));
for(ix in 1:(num_shared_sample_ids/2)){
	aid=A_sample_ids[ix];
	bid=B_sample_ids[ix];
	points(
		c(mds[aid,1], mds[bid,1]),
		c(mds[aid,2], mds[bid,2]),
		type="l", lwd=line_wts[ix], col="black");
}
text(mds[,1], mds[,2], glyphs, col=colors);

#------------------------------------------------------------------------------

plot_paired_mds_colored_by_group=function(a_ids, b_ids, mds_res, fact, labels){

	if(!all(names(a_ids)==rownames(fact))){
		cat("Error!  Subject ID if A IDs don't match the Subject IDs of Factor Matrix.\n");
		quit(-1);
	}

	orig_par=par(no.readonly=T);

	# Make space for plot and legend
	layout_mat=matrix(c(1,1,1,1,1,2), nrow=6, ncol=1);
	layout(layout_mat);

	num_factors=ncol(fact);
	if(num_factors<1){
		return;
	}

	cat("Num of Factors: ", num_factors, "\n");
	factor_names=colnames(fact);
	for(fi in 1:num_factors){

		cur_fact=as.character(fact[,fi,]);
		cur_fact_name=factor_names[fi];
		uniq_levels=sort(unique(cur_fact));
		num_levels=length(uniq_levels);
		if(num_levels>10){
			cat(cur_fact_name, ": Too many levels to color. (", num_levels, ")\n", sep="");
			next;
		}else{
			cat("Working on: ", cur_fact_name, "\n");
			cat("Num Levels: ", num_levels, "\n");
			print(cur_fact);

			color_map=rainbow(num_levels, start=0, end=4/6);
			names(color_map)=uniq_levels;

			cat("Color Map:\n");
			print(color_map);

			cat("Color Assignment:\n");
			colors_list=color_map[cur_fact];
			print(colors_list);

			# Generate MDS plot
			mds_x_range=range(mds_res[,1]);
			mds_y_range=range(mds_res[,2]);

			mds_x_width=diff(mds_x_range);
			mds_y_height=diff(mds_y_range);

			par(mar=c(4,4,4,4));
			plot(0,0, type="n", xlab="Dim 1", ylab="Dim 2", 
				xlim=c(mds_x_range[1]-mds_x_width*.05, mds_x_range[2]+mds_x_width*.05),
				ylim=c(mds_y_range[1]-mds_y_height*.05, mds_y_range[2]+mds_y_height*.05),
				main=paste("Colored by: ", cur_fact_name, sep=""));

			draw_centroids(A_centroid, B_centroid);

			for(ix in 1:length(a_ids)){
				aid=a_ids[ix];
				bid=b_ids[ix];
				points(
					c(mds_res[aid,1], mds_res[bid,1]),
					c(mds_res[aid,2], mds_res[bid,2]),
					type="l", lwd=1, col=colors_list[ix]);
			}
			text(mds_res[,1], mds_res[,2], labels, col="black", cex=1.3);

			# Generate Legend
			par(mar=c(0,0,0,0));
			plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n", xlab="", ylab="", bty="n", 
				xaxt="n", yaxt="n", main="");

			legend(0, 1, fill=color_map, legend=uniq_levels, bty="n");	
		}

		cat("\n\n");
	}

	par(orig_par);
}


plot_paired_mds_colored_by_group(A_sample_ids, B_sample_ids, mds, factors, glyphs);

#------------------------------------------------------------------------------

# Plot Histogram
par(mfrow=c(1,1));
par(oma=c(0,0,2,0));
par(mar=c(4,4,1,1));

cutoffs=quantile(paired_dist, c(0,.25,.75, 1));
plot_labl=c("Most Similar Quartile (0-25%)", "Median Quartile (25-50%)", "Most Dissimilar Quartile (75-100%)");

h=hist(paired_dist, breaks=20,
	xlab=paste(DistType, " Distances", sep=""), main="Paired Distances Distribution");
abline(v=cutoffs, col="blue", lty=2);
text(cutoffs, mean(range(h$counts)), c("0%", "25%", "75%", "100%"), adj=c(.5, -.5), srt=90, col="blue");

#------------------------------------------------------------------------------ 
layout_mat=matrix(
	c(2,2,
	  1,3), byrow=T, nrow=2);
layout(layout_mat);
par(mar=c(3,3,5,1));

for(i in 1:3){

	plot(0,0, type="n", xlab="Dim 1", ylab="Dim 2",
		xlim=c(mds_x_range[1]-mds_x_width*.1, mds_x_range[2]+mds_x_width*.1),
		ylim=c(mds_y_range[1]-mds_y_height*.1, mds_y_range[2]+mds_y_height*.1),
		main=plot_labl[i]);
	
	for(ix in 1:(num_shared_sample_ids/2)){
		if(paired_dist[ix]>=cutoffs[i] && paired_dist[ix]<=cutoffs[i+1]){
			aid=A_sample_ids[ix];
			bid=B_sample_ids[ix];
			points(
				c(mds[aid,1], mds[bid,1]),
				c(mds[aid,2], mds[bid,2]),
				type="l", lwd=.5, col="grey30");
			text( c(mds[aid,1], mds[bid,1]),
				c(mds[aid,2], mds[bid,2]),
				glyphs[c(aid,bid)], col=colors[c(aid,bid)]
			);
		}

	}
}

##############################################################################
print(paired_dist);

num_model_pred=length(model_var_arr);

cov_model_string=paste(model_var_arr, collapse="+");
mmat=model.matrix(as.formula(paste("~", cov_model_string, "-1")), data=as.data.frame(factors));
cov_coeff_names=c("(Intercept)", colnames(mmat));
num_cov_coeff_names=length(cov_coeff_names);

cat("Anticipated Coefficient Names:\n");
print(cov_coeff_names);
cat("\n");

covariates_coef_mat  =matrix(NA, nrow=1, ncol=num_cov_coeff_names, 
	dimnames=list("distance", cov_coeff_names));
covariates_pval_mat  =matrix(NA, nrow=1, ncol=num_cov_coeff_names, 
	dimnames=list("distance", cov_coeff_names));

rsqrd_mat	=matrix(NA, nrow=1, ncol=2, 
	dimnames=list("distance", c("R^2", "Adj-R^2")));

model_pval_mat	=matrix(NA, nrow=1, ncol=1, 
	dimnames=list("distance", c("F-statistic P-value")));

###############################################################################

model_data=cbind(factors, paired_dist);

if(length(model_var_arr)==0){
	model_str=paste("paired_dist ~ 1");
}else{
	model_str=paste("paired_dist ~ ", paste(model_var_arr, collapse=" + "), sep="");
}

NUM_BS=2000+length(model_var_arr)*350;

#NUM_BS=40;

if(NUM_BS<500){
	plot_text(c(
		"***************************************************************",
		paste("WARNING: Number of Bootstraps is <500.  (", NUM_BS, ")"),
		"***************************************************************"
		));
}

num_data_rows=nrow(model_data);

obs_lm_fit=lm(as.formula(model_str), data=model_data);
obs_lm_fit_sum=summary(obs_lm_fit);
obs_lm_fit_sum_text=capture.output(print(obs_lm_fit_sum));


cat("Running regression bootstraps on:\n");
cat("Num Bootstraps: ", NUM_BS, "\n");
cat(model_str, "\n");

bs_ix=1;
while(bs_ix <=NUM_BS){
	bs_data=model_data[sample(num_data_rows, replace=T),, drop=F];

	result=tryCatch({
		lm_fit=lm(as.formula(model_str), data=bs_data);
		}, error=function(e){ cat("Degenerate bootstrap:, ", as.character(e), "\n")}
	);

	if(is.null(result)){
		next;
	}else{
		lm_fit=result;
	}

	if(bs_ix==1){
		# Allocated matrix based on observed samples, to ensure we have
		# columns for all coefficients.
		# It's possible for a bootstrap outcome to be missing categories for a
		# categorical variable, so dummy variables (and thus coefficients) are missing.
		coefficients=matrix(NA, nrow=NUM_BS, ncol=length(obs_lm_fit$coefficients));
		colnames(coefficients)=names(obs_lm_fit$coefficients);
	}

	coefficients[bs_ix, names(lm_fit$coefficients)]=lm_fit$coefficients;
	bs_ix=bs_ix+1;
}

#print(coefficients);

num_coef=ncol(coefficients);
regression_table=matrix(NA, nrow=num_coef, ncol=4);
rownames(regression_table)=colnames(coefficients);
colnames(regression_table)=c(
	"Estimate",
	"LB 95%",
	"UB 95%",
	"P-value"
);
	

means=apply(coefficients, 2, function(x){mean(x, na.rm=T)});
lb95=apply(coefficients, 2, function(x){quantile(x, .025, na.rm=T)});
ub95=apply(coefficients, 2, function(x){quantile(x, .975, na.rm=T)});
prob_gt0=apply(coefficients, 2, function(x){ notna=!is.na(x); x=x[notna]; mean(x<=0) });

not0_pval=sapply(prob_gt0, function(x){ min(min(x, 1-x)*2,1)});

regression_table[,"Estimate"]=means;
regression_table[,"LB 95%"]=lb95;
regression_table[,"UB 95%"]=ub95;
regression_table[,"P-value"]=not0_pval;

reg_tab_char=apply(regression_table, 1:2, function(x){sprintf("%7.4f",x)});

Signf=sapply(regression_table[,"P-value"], function(x){
	if(is.na(x)){ return("");}
	else if(x<=.001){return("***");}
	else if(x<=.01){return("**");}
	else if(x<=.05){return("*");}
	else if(x<=.1){return(".");}
	return("");
});

regression_table_text=capture.output(print(cbind(reg_tab_char, Signf), quote=F));

print(regression_table_text, quote=F);
print(obs_lm_fit_sum);

par(mfrow=c(1,1));
plot_text(c(
"Bootstrapped Stats:",
regression_table_text));


################################################################################

outfn=paste(OutputRoot, ".paired_dist.pval.tsv", sep="");
if(TagName!=""){
	fh=file(outfn, "w");
	cat(file=fh, TagName, "\n\t");
}
write.table(round(regression_table[,c("Estimate","P-value")], 4), 
	file=outfn, append=T,
        sep="\t", row.names=T, col.names=T);

################################################################################

plot_text(c(
"Standard Stats:",
obs_lm_fit_sum_text));

###############################################################################

rename_intercept=function(mat, old_n, new_n){
	cname=colnames(mat);
	old_ix=which(cname==old_n);
	cname[old_ix]=new_n;
	colnames(mat)=cname;
	return(mat);
}

coef_mat=regression_table[,"Estimate", drop=F];
pval_mat=regression_table[,"P-value", drop=F];

coef_mat=rename_intercept(coef_mat, "(Intercept)", "\"Difference\"");
pval_mat=rename_intercept(pval_mat, "(Intercept)", "\"Difference\"");

paint_matrix(coef_mat, 
        title="Regression Model Bootstrapped Coefficient Values",
        high_is_hot=T, deci_pts=2, value.cex=1);

paint_matrix(pval_mat, 
        title="Regression Model Bootstrapped Coeff. P-values",
	plot_min=0, plot_max=1,
        high_is_hot=F, deci_pts=2, value.cex=1);

signf_coef_mat=mask_matrix(coef_mat, pval_mat, .1, 0);
paint_matrix(signf_coef_mat,
        title="Regression Model Bootstrapped Significant Coefficients", 
        high_is_hot=T, deci_pts=2, value.cex=1, label_zeros=F);
mtext(text="(P-Values < 0.10 Shown)", side=3, cex=.9, font=3, line=2, outer=T);

################################################################################

plot_abund_wCI=function(abundances, lb, ub, num_top_categories=10, category_colors, ymax, title, ylab){
        # This draws a single Rank Abundance Plot

        if(!is.null(dim(abundances))){
                abundances=abundances[1,];
        }

	cat_names=names(abundances);
        mids=barplot(abundances, col=category_colors, names.arg="", ylim=c(0, ymax*1.05), ylab=ylab);

	dist_bt_mids=mids[2]-mids[1];
	error_bar_half_wid=dist_bt_mids/8;

	num_top_categories=min(num_top_categories, length(abundances));

	for(catix in 1:num_top_categories){
		if(ub[catix]!=lb[catix]){
			points(
				c(mids[catix]-error_bar_half_wid, mids[catix]+error_bar_half_wid),
				rep(ub[catix],2), lwd=.8, type="l");
			points(
				c(mids[catix]-error_bar_half_wid, mids[catix]+error_bar_half_wid),
				rep(lb[catix],2), lwd=.8, type="l");
			points(
				rep(mids[catix], 2),
				c(ub[catix], lb[catix]), lwd=.8, type="l");
		}

	}

        title(main=title, font.main=2, cex.main=1.75, line=-1);

        # Compute and label categories beneath barplot
        bar_width=mids[2]-mids[1];
        plot_range=par()$usr;
        plot_height=plot_range[4];
        label_size=min(c(1,.7*bar_width/par()$cxy[1]));
        text(mids-par()$cxy[1]/2, rep(-par()$cxy[2]/2, num_top_categories), 
		cat_names, srt=-45, xpd=T, pos=4, cex=label_size);

}

plot_lograt_wCI=function(lrat, lb, ub, num_top_categories=10, ymax, category_colors, title, ylab){
        # This draws a single Rank Abundance Plot

	cat_names=names(lrat);
        plot(0,0, type="n", ylim=c(-ymax*1.05, ymax*1.05), xlim=c(0, num_top_categories+2), ylab=ylab,
		bty="n", xaxt="n", yaxt="n");
	abline(h=c(-1, 1), col="grey75", lty="dotted", lwd=.7);
	abline(h=c(-2, 2), col="grey50", lty="dashed", lwd=.8);
        mids=barplot(lrat, col=category_colors, names.arg="", add=T);

	dist_bt_mids=mids[2]-mids[1];
	error_bar_half_wid=dist_bt_mids/8;

	for(catix in 1:num_top_categories){
		if(ub[catix]!=lb[catix]){
			points(
				c(mids[catix]-error_bar_half_wid, mids[catix]+error_bar_half_wid),
				rep(ub[catix],2), lwd=.8, type="l");
			points(
				c(mids[catix]-error_bar_half_wid, mids[catix]+error_bar_half_wid),
				rep(lb[catix],2), lwd=.8, type="l");
			points(
				rep(mids[catix], 2),
				c(ub[catix], lb[catix]), lwd=.8, type="l");
		}

	}

        title(main=title, font.main=2, cex.main=1.75, line=0);

        # Compute and label categories beneath barplot
        bar_width=mids[2]-mids[1];
        plot_range=par()$usr;
        plot_height=plot_range[4];
        label_size=min(c(1,.7*bar_width/par()$cxy[1]));
        text(mids-par()$cxy[1]/2, -ymax+rep(-par()$cxy[2]/2, num_top_categories), 
		cat_names, srt=-45, xpd=T, pos=4, cex=label_size);

}

plot_comparisons=function(a_profs, b_profs, global_colormap, title, a_name, b_name, topN=15, ylab, max_abs_lr){
	cat("Plotting: ", ylab, "\n");


	comb_prof=rbind(a_profs, b_profs);
	avg_prof=apply(comb_prof, 2, mean);
	
	sort_ix=order(avg_prof, decreasing=T);

	topN=min(topN, ncol(comb_prof));

	avg_prof=avg_prof[sort_ix[1:topN]];
	a_profs=a_profs[,sort_ix[1:topN], drop=F];
	b_profs=b_profs[,sort_ix[1:topN], drop=F];

	top_categories=names(avg_prof);
	bar_col=global_colormap[top_categories];
	bar_col[is.na(bar_col)]="grey";


	minabund=min(a_profs[a_profs>0], b_profs[b_profs>0])/10;
	a_mean_prof=apply(a_profs, 2, mean);
	b_mean_prof=apply(b_profs, 2, mean);
	lrab_prof=log10((a_mean_prof+minabund)/(b_mean_prof+minabund));

	# Compute 95% CI
	num_samp=nrow(a_profs);
	if(num_samp>=38){

		NUMBS=1000;
		a_bs_means=matrix(NA, nrow=NUMBS, ncol=topN);
		b_bs_means=matrix(NA, nrow=NUMBS, ncol=topN);
		lrab_bs=matrix(NA, nrow=NUMBS, ncol=topN);

		for(bsix in 1:NUMBS){
			samples=sample(num_samp, replace=T);
			a_bs_means[bsix,]=apply(a_profs[samples,,drop=F], 2, mean);
			b_bs_means[bsix,]=apply(b_profs[samples,,drop=F], 2, mean);
			lrab_bs[bsix,]=log10((a_bs_means[bsix,]+minabund)/(b_bs_means[bsix,]+minabund));
		}

		a_95ci_ub=apply(a_bs_means, 2, function(x){ quantile(x, .975, na.rm=T);});
		a_95ci_lb=apply(a_bs_means, 2, function(x){ quantile(x, .025, na.rm=T);});

		b_95ci_ub=apply(b_bs_means, 2, function(x){ quantile(x, .975, na.rm=T);});
		b_95ci_lb=apply(b_bs_means, 2, function(x){ quantile(x, .025, na.rm=T);});

		lrab_95ci_ub=apply(lrab_bs, 2, function(x){ quantile(x, .975, na.rm=T);});
		lrab_95ci_lb=apply(lrab_bs, 2, function(x){ quantile(x, .025, na.rm=T);});
	}else{
		a_95ci_ub=a_mean_prof;
		a_95ci_lb=a_mean_prof;
		b_95ci_ub=b_mean_prof;
		b_95ci_lb=b_mean_prof;
		lrab_95ci_ub=lrab_prof;
		lrab_95ci_lb=lrab_prof;
	}

	ymax=max(a_95ci_ub, b_95ci_ub);
	lr_ymax=max(abs(c(lrab_95ci_ub, lrab_95ci_lb, max_abs_lr)));

	par(mar=c(6,6,2,1));
	plot_abund_wCI(a_mean_prof, a_95ci_lb, a_95ci_ub, num_top_categories=topN, 
		ymax=ymax, bar_col, title=a_name, ylab=ylab);

	par(mar=c(6,2,2,1));
	plot_abund_wCI(b_mean_prof, b_95ci_lb, b_95ci_ub, num_top_categories=topN, 
		ymax=ymax, bar_col, title=b_name, ylab="");

	par(mar=c(6,2,2,3));
	plot_lograt_wCI(lrab_prof, lb=lrab_95ci_lb, ub=lrab_95ci_ub, num_top_categories=topN, 
		ymax=lr_ymax, bar_col, title=paste("Log10(", a_name, "/", b_name,")", sep=""), ylab="");

	stats=list();
	stats[["max_abs_lrab"]]=max(abs(c(lrab_95ci_ub, lrab_95ci_lb)));
	stats[["max_ub_abund"]]=max(a_95ci_ub, b_95ci_ub);
	return(stats);

}

################################################################################

rownames(good_pairs_map)=good_pairs_map[,1];
sorted_paired_dist=sort(paired_dist, decreasing=T);
sorted_A_samp_ids=names(sorted_paired_dist);
sorted_good_pairs_map=good_pairs_map[sorted_A_samp_ids,];

print(sorted_paired_dist)
print(sorted_good_pairs_map);

num_dist=length(sorted_paired_dist);
splits_arr=1:(min(6,num_dist-1));
num_splits=length(splits_arr);

num_cat_to_plot=15;
num_colors_to_asgn=round(num_cat_to_plot*1.1)
num_categories=ncol(normalized);
color_map=rainbow(num_colors_to_asgn, start=0, end=2/3);

all_means=apply(normalized, 2, mean);
all_means_sorted=sort(all_means, decreasing=T);
names(color_map)=names(all_means_sorted[1:num_colors_to_asgn]);

glob_max_abs_lr=0;
for(spix in 1:num_splits){

	splits=splits_arr[spix];
	cat("Splitting data to: ", splits, "\n");
	split_pts=unique(c(as.integer(seq(1, num_dist, num_dist/splits)), num_dist));

	cat("Ranges:\n");
	print(split_pts);

	par(mfrow=c(splits,3));

	for(spl_pts in 1:splits){
		sidx=split_pts[spl_pts];
		eidx=split_pts[spl_pts+1];
		
		quantile_lb=round(100*(sidx-1)/(num_dist-1));
		quantile_ub=round(100*(eidx-1)/(num_dist-1));

		if(spl_pts==splits){
			eidx=eidx+1;
		}
	
		cat("Grabbing data from ", sidx, "<= i <", eidx, "\n");	

		mean_dist=mean(sorted_paired_dist[sidx:(eidx-1)]);
		stats=plot_comparisons(
			a_profs=normalized[sorted_good_pairs_map[sidx:(eidx-1), 1],, drop=F],
			b_profs=normalized[sorted_good_pairs_map[sidx:(eidx-1), 2],, drop=F],
			global_colormap=color_map,
			a_name=A_subtrahend,
			b_name=B_minuend,
			topN=num_cat_to_plot,
			ylab=paste("Mean Dist = ", round(mean_dist,3), 
				"\nQuantile : ", quantile_lb, "% - ", quantile_ub, "%", 
				"\nN = ", (eidx-sidx), sep=""),
			max_abs_lr=glob_max_abs_lr
			);

		glob_max_abs_lr=max(glob_max_abs_lr, stats[["max_abs_lrab"]]);
	}

}

################################################################################
# Perform statistics on dispersion

plot_mds=function(mds, title="", label=F, col="black"){

	x=mds[,1];
	y=mds[,2];
	samp_ids=rownames(mds);

	xcentroid=mean(x);
	ycentroid=mean(y);

	xrange=range(x);
	yrange=range(y);

	xpad=diff(xrange)*.2;
	ypad=diff(yrange)*.1;

	plot(0, type="n",
		xlim=c(xrange[1]-xpad, xrange[2]+xpad),
		ylim=c(yrange[1]-ypad, yrange[2]+ypad),
		main=title
	)

	points(xcentroid, ycentroid, type="p", bg=col, col="black", cex=3, pch=21);

	if(label==F){
		points(x,y, type="p", col=col);
	}else{
		text(x,y, samp_ids, cex=.7, col=col);
	}

}

plot_dist_bars=function(Adist_arr, Bdist_arr, AtoBMap, title="", acol="black", bcol="black"){

	Bnames_map=AtoBMap[,2];
	names(Bnames_map)=AtoBMap[,1];

	par(mfrow=c(2,1));
	

	max_dist=max(c(Adist_arr, Bdist_arr));

	a_order_ix=order(Adist_arr, decreasing=T);
	Adist_arr_sorted=Adist_arr[a_order_ix];
	Anames=names(Adist_arr_sorted);
	Bnames=Bnames_map[Anames];

	num_samp=length(Adist_arr);

	cexadj=min(1,60/num_samp);
	cat("cex adj: ", cexadj, "\n");

	par(mar=c(13,4,4,1));
	mids=barplot(Adist_arr_sorted, ylim=c(0, max_dist*1.05), 
		names.arg=Anames, las=2, main=title, col=acol, cex.names=cexadj);
	points(mids, Bdist_arr[a_order_ix], type="p", pch=3, cex=cexadj*.9, col=bcol);

	par(mar=c(13,4,0,1));
	mids=barplot(Bdist_arr[Bnames], ylim=c(0, max_dist*1.05),
		names.arg=Bnames, las=2, col=bcol, cex.names=cexadj);
	points(mids, Adist_arr[a_order_ix], type="p", pch=3, cex=cexadj*.9, col=acol);

}

plot_indiv_dist_diff=function(Adist_arr, Bdist_arr, AtoBMap, Aname, Bname, acol, bcol){
	
	diff=Bdist_arr-Adist_arr;
	cor_pear=cor.test(Adist_arr, Bdist_arr, method="pearson");
	cor_spear=cor.test(Adist_arr, Bdist_arr, method="spearman");

	par(mfrow=c(1,1));
	par(mar=c(4,13,4,13));
	rngs=range(diff);
	max_diff=max(abs(rngs));

	diff_sorted=sort(diff, decreasing=F);
	

	mids=barplot(diff_sorted, horiz=T, xlim=c(-max_diff*1.03, max_diff*1.03),  names.arg="",
		main=paste("Ordered by Differences from Centroid: ", Bname, "-", Aname));

	title(main=paste("Pearson's Cor: ", round(cor_pear$estimate, 3), 
		", p-value: ", round(cor_pear$p.value,3), sep=""), line=0, cex.main=.8, font.main=1);
	title(main=paste("Spearman's Cor: ", round(cor_spear$estimate, 3), 
		", p-value: ", round(cor_spear$p.value,3), sep=""), line=-.8, cex.main=.8, font.main=1);

	sorted_names=names(diff_sorted);

	a_to_b=AtoBMap[,1];
	names(a_to_b)=AtoBMap[,2];

	num_samp=length(Adist_arr);
	cexadj=min(1,60/num_samp);

	axis(4, mids, sorted_names, las=2, cex.axis=cexadj);
	axis(2, mids, a_to_b[sorted_names], las=2, cex.axis=cexadj);

	axis(4, -3, paste(Bname, "\nfurther from centroid"), 
		cex.axis=1, font.axis=2, las=2, line=-2, tick=F, col.axis=bcol);

	axis(2, -3, paste(Aname, "\nfurther from centroid"), 
		cex.axis=1, font.axis=2, las=2, line=-2, tick=F, col.axis=acol);

}

plot_grp_diff=function(Adist_arr, Bdist_arr, Aname, Bname, acol, bcol){

	diff=Bdist_arr-Adist_arr;

	wt_res=wilcox.test(Adist_arr, Bdist_arr, alternative=c("two.sided"));
	mean_diff=mean(diff);
	med_diff=median(diff);
	diff_range=range(diff);
	max_diff=max(abs(diff(diff_range)));

	par(mfrow=c(1,1));
	par(mar=c(4.1,4.1,6,1));
	hist(diff, breaks=2*nclass.Sturges(diff),
		xlim=c(-max_diff*1.1, max_diff*1.1),
		xlab="Differences in Distance from Centroid", 
		main=""
		);

	cur_par=par();
	plot_limits=cur_par$usr;
	max_yplot=plot_limits[4];

	# Annotate with vertical lines
	abline(v=0, col="grey", lwd=3);
	abline(v=med_diff, col="blue", lwd=2);
	abline(v=mean_diff, col="purple", lwd=2);

	# Label vertical lines
	if(med_diff>mean_diff){
		tpos=c(4,2);
	}else{
		tpos=c(2,4);
	}
	text(med_diff, max_yplot, "\nMedian", col="blue", pos=tpos[1]);
	text(mean_diff, max_yplot, "\nMean", col="purple", pos=tpos[2]);
	
	# Label statistics
	title(main=paste("Difference of Dispersion: ", Bname, " - ", Aname, sep=""),
		cex.main=1.5, font.main=2, line=4);
	title(main=paste("Mean: ", round(mean_diff, 4), "  Median: ", round(med_diff,4), sep=""), 
		cex.main=.8, font.main=1, line=2);
	title(main=paste("Wilcoxon Difference in Dispersion: p-value = ", round(wt_res$p.value, 4), sep=""),
		cex.main=.8, font.main=1, line=1);
	
	
}

plot_dist_hist=function(a_dist_arr, b_dist_arr, a_name, b_name, acol, bcol){

	comb_arr=c(a_dist_arr, b_dist_arr);
	comb_range=range(comb_arr);
	comb_breaks=2*nclass.Sturges(comb_arr);
	hist_res=hist(comb_arr, comb_breaks, plot=F);

	a_hist_res=hist(a_dist_arr, hist_res$breaks, plot=F);
	b_hist_res=hist(b_dist_arr, hist_res$breaks, plot=F);
	print(hist_res);

	a_mean=mean(a_dist_arr);
	b_mean=mean(b_dist_arr);

	par(mfrow=c(2,1));
	par(mar=c(4.1,4,4,2));

	# Plot and Label for A
	hist(a_dist_arr, breaks=hist_res$breaks,
		xlim=c(0, comb_range[2]),
		xlab=paste("Distances from ", a_name, "'s Centroid", sep=""),
		main=a_name, col=acol
		);
	title(main=paste("Avg Distance from Centroid: ", round(a_mean,3), sep=""), 
		font.main=1, cex.main=.8, line=1);
	abline(v=a_mean, col="black", lwd=2);
	abline(v=b_mean, col=bcol, lwd=1, lty=3);

	# Plot and Label for B
	hist(b_dist_arr, breaks=hist_res$breaks,
		xlim=c(0, comb_range[2]),
		xlab=paste("Distances from ", b_name, "'s Centroid", sep=""),
		main=b_name, col=bcol
		);
	title(main=paste("Avg Distance from Centroid: ", round(b_mean,3), sep=""), 
		font.main=1, cex.main=.8, line=1);
	abline(v=b_mean, col="black", lwd=2);
	abline(v=a_mean, col=acol, lwd=1, lty=3);

}

regress_dispersion=function(Adist_arr, Bdist_arr, Aname, Bname, model_var, factors){

	A_samp_ids=names(Adist_arr);
	
	num_model_var=length(model_var);
	cat("RegressDispersion: Number of Predictors: ", num_model_var, "\n");

	pred_string=paste(model_var, collapse=" + ");
	A_disp_model_string=paste("Adist_arr ~ ", pred_string, sep="");
	B_disp_model_string=paste("Bdist_arr ~ ", pred_string, sep="");

	Adisp_fit=lm(as.formula(A_disp_model_string), data=factors);
	print(Adisp_fit);
	print(summary(Adisp_fit));

	Bdisp_fit=lm(as.formula(B_disp_model_string), data=factors);
	print(Bdisp_fit);
	print(summary(Bdisp_fit));

	skip=function(arr, n){
		coef_ix=min(which(arr=="Coefficients:"));
		len=length(arr);
		return(arr[coef_ix:len]);
	}

	par(mfrow=c(1,1));
	plot_text(c(
		"Associations on Dispersion (Assuming Normally Distributed Residuals)",
		"",
		paste("[", Aname, "]", sep=""),
		skip(capture.output({print(summary(Adisp_fit))}), 8),
		"",
		"",
		paste("[", Bname, "]", sep=""),
		skip(capture.output({print(summary(Bdisp_fit))}), 8)
	));

}

regress_diff_dispersion=function(Adist_arr, Bdist_arr, Aname, Bname, model_var, factors){

	centr_diff=Bdist_arr-Adist_arr;
	A_samp_ids=names(centr_diff);
	
	num_model_var=length(model_var);
	cat("RegressDifferenceDispersion: Number of Predictors: ", num_model_var, "\n");

	pred_string=paste(model_var, collapse=" + ");
	diff_disp_model_string=paste("centr_diff ~ ", pred_string, sep="");

	diff_fit=lm(as.formula(diff_disp_model_string), data=factors);
	print(diff_fit);
	print(summary(diff_fit));

	skip=function(arr, n){
		coef_ix=min(which(arr=="Coefficients:"));
		len=length(arr);
		return(arr[coef_ix:len]);
	}

	par(mfrow=c(1,1));
	plot_text(c(
		paste("Associations on Difference of Dispersion: ", Bname, " - ", Aname, sep=""),
		paste("  (Assuming normally distributed residuals)"),
		"",
		skip(capture.output({print(summary(diff_fit))}), 8)
	));
}


bootstrap_regression_dispersion=function(dist_arr, name, model_var, factors, num_bs){

	bs_ix=1;
	num_samp=length(dist_arr);
	while(bs_ix <=num_bs){

		bs_samp_ix=sample(num_samp, replace=T)

		bs_dists=dist_arr[bs_samp_ix];
		bs_factors=factors[bs_samp_ix,, drop=F];
		bs_data=cbind(bs_dists, bs_factors);

		model_str=paste("bs_dists ~ ", paste(model_var, collapse=" + ", sep=""), sep="");

		#cat("Model String: ", model_str, "\n");
		#cat("BS Idx: ", bs_samp_ix, "\n");

		result=tryCatch({
			lm_fit=lm(as.formula(model_str), data=bs_data);
			}, error=function(e){ cat("Degenerate bootstrap:, ", as.character(e), "\n")}
		);

		if(is.null(result)){
			next;
		}else{
			lm_fit=result;
		}

		if(bs_ix==1){
			# Allocated matrix based on observed samples, to ensure we have
			# columns for all coefficients.
			# It's possible for a bootstrap outcome to be missing categories for a
			# categorical variable, so dummy variables (and thus coefficients) are missing.
			coefficients=matrix(NA, nrow=num_bs, ncol=length(obs_lm_fit$coefficients));
			colnames(coefficients)=names(obs_lm_fit$coefficients);
		}

		coefficients[bs_ix, names(lm_fit$coefficients)]=lm_fit$coefficients;
		bs_ix=bs_ix+1;
	}

	means=apply(coefficients, 2, function(x){mean(x, na.rm=T)});
	lb95=apply(coefficients, 2, function(x){quantile(x, .025, na.rm=T)});
	ub95=apply(coefficients, 2, function(x){quantile(x, .975, na.rm=T)});
	prob_gt0=apply(coefficients, 2, function(x){ notna=!is.na(x); x=x[notna]; mean(x<=0) });

	not0_pval=sapply(prob_gt0, function(x){ min(min(x, 1-x)*2,1)});

	regression_table[,"Estimate"]=means;
	regression_table[,"LB 95%"]=lb95;
	regression_table[,"UB 95%"]=ub95;
	regression_table[,"P-value"]=not0_pval;

	reg_tab_char=apply(regression_table, 1:2, function(x){sprintf("%7.4f",x)});

	Signf=sapply(regression_table[,"P-value"], function(x){
		if(is.na(x)){ return("");}
		else if(x<=.001){return("***");}
		else if(x<=.01){return("**");}
		else if(x<=.05){return("*");}
		else if(x<=.1){return(".");}
		return("");
	});

	regression_table_text=capture.output(print(cbind(reg_tab_char, Signf), quote=F));

	return(regression_table_text);
}


################################################################################

par(mfrow=c(1,1));
plot(0,0, xlim=c(-1,1), ylim=c(-1,1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="", type="n");
text(0,0, "Dispersion Analysis", cex=3, font=2);


Anorm=normalized[A_sample_ids,];
Bnorm=normalized[B_sample_ids,];

Amean=apply(Anorm, 2, mean);
Bmean=apply(Bnorm, 2, mean);

Anorm=rbind(Anorm, Amean);
Bnorm=rbind(Bnorm, Bmean);

Adistmat=as.matrix(compute_dist(Anorm, DistType));
Bdistmat=as.matrix(compute_dist(Bnorm, DistType));

numAsamp=length(A_sample_ids);
numBsamp=length(B_sample_ids);

for(i in 1:length(Adistmat)){
        if(Adistmat[i]==0){
                Adistmat[i]=1e-323;
        }
}

for(i in 1:length(Bdistmat)){
        if(Bdistmat[i]==0){
                Bdistmat[i]=1e-323;
        }
}

Amds=cmdscale(Adistmat);
Bmds=cmdscale(Bdistmat);

par(mfrow=c(2,2));
plot_mds(Amds[1:numAsamp,], title=A_subtrahend, label=F, col="green");
plot_mds(Amds[1:numAsamp,], title=A_subtrahend, label=T, col="green");
plot_mds(Bmds[1:numBsamp,], title=B_minuend, label=F, col="blue");
plot_mds(Bmds[1:numBsamp,], title=B_minuend, label=T, col="blue");

A_dist_fr_centr=Adistmat["Amean", 1:numAsamp];
B_dist_fr_centr=Bdistmat["Bmean", 1:numBsamp];

plot_dist_bars(A_dist_fr_centr, B_dist_fr_centr, good_pairs_map[,c(1,2)], 
	paste("Ordered By", A_subtrahend), acol="green", bcol="blue");

plot_dist_bars(B_dist_fr_centr, A_dist_fr_centr, good_pairs_map[,c(2,1)], 
	paste("Ordered By", B_minuend), acol="blue", bcol="green");

plot_dist_hist(A_dist_fr_centr, B_dist_fr_centr, A_subtrahend, B_minuend, acol="blue", bcol="green");

regress_dispersion(A_dist_fr_centr, B_dist_fr_centr, A_subtrahend, B_minuend, model_var_arr, factors);

A_bs_reg_tab=bootstrap_regression_dispersion(A_dist_fr_centr, A_subtrahend, model_var_arr, 
	factors, NUM_BS);
B_bs_reg_tab=bootstrap_regression_dispersion(B_dist_fr_centr, B_minuend, model_var_arr, 
	factors, NUM_BS);

plot_text(c(
	paste("Associations on Dispersion (Bootstrapped Regression, num bootstraps: ", NUM_BS, ")", sep=""),
	"",
	paste("[", A_subtrahend, "]", sep=""),
	"Coefficients:",
	A_bs_reg_tab,
	"",
	"",
	paste("[", B_minuend, "]", sep=""),
	"Coefficients:",
	B_bs_reg_tab
));

par(mfrow=c(1,1));
plot(0,0, xlim=c(-1,1), ylim=c(-1,1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="", type="n");
text(0,0, "Differences in\nDispersion Analysis", cex=3, font=2);

plot_indiv_dist_diff(A_dist_fr_centr, B_dist_fr_centr, good_pairs_map[,c(1,2)],
	A_subtrahend, B_minuend, acol="blue", bcol="green");

plot_grp_diff(A_dist_fr_centr, B_dist_fr_centr,
	A_subtrahend, B_minuend, acol="blue", bcol="green");

regress_diff_dispersion(A_dist_fr_centr, B_dist_fr_centr, A_subtrahend, B_minuend, model_var_arr, factors);
diff_bs_reg_tab=bootstrap_regression_dispersion(
	B_dist_fr_centr-A_dist_fr_centr, 
	paste(B_minuend, "-", A_subtrahend), model_var_arr, factors, NUM_BS);

par(mfrow=c(1,1));
plot_text(c(
	paste("Associations on Difference of Dispersion: ", B_minuend, " - ", A_subtrahend, sep=""),
	paste("  (Bootstrapped Regression, num bootstraps: ", NUM_BS, ")", sep=""),
	"",
	"Coefficients:",
	diff_bs_reg_tab
));


################################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
