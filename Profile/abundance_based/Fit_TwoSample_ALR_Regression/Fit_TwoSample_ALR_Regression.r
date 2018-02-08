#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
library(car);

source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');

options(useFancyQuotes=F);


params=c(
	"summary_file", "s", 1, "character",
	"summary_file2", "S", 2, "character",
	"pairings", "p", 1, "character",
	"factors", "f", 1, "character",
	"factor_samp_id_name", "F", 2, "character",
	"model_var", "M", 1, "character",
	"required", "q", 2, "character",
	"response", "e", 1, "character",
	"predictor", "g", 1, "character",

	"num_resp_var", "u", 2, "numeric",
	"num_pred_var", "v", 2, "numeric",

	"reference_levels", "c", 2, "character",
	"outputroot", "o", 2, "character",

	"contains_remaining", "R", 2, "logical",
	"shorten_category_names", "x", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

NUM_TOP_RESP_CAT=35;
NUM_TOP_PRED_CAT=10;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table>\n",
	"	[-S <second summary file table, in case Resp and Pred were in different files.>]\n",
	"\n",
	"	-p <pairings map, pairing Resp and Pred sample IDs>\n",
	"	-f <factors file, contains covariates and factors>\n",
	"	-F <column name of sample ids in factors, should be either response or predictor ALR name>\n",
	"	-M <list of covariate X's names to include in the model from the factor file>\n",
	"	-e <response ALR name, (column name in pairings file)\n",
	"	-g <predictor ALR name, (column name in pairings file)\n",
	"	[-q <required list of variables to include after NA removal>]\n",
	"\n",
	"	[-u <number of top response categories to analyze, default=", NUM_TOP_RESP_CAT, ">]\n",
	"	[-v <number of top predictor (as ALR) categories to include, default=", NUM_TOP_PRED_CAT, ">]\n",
	"\n",
	"	[-o <output filename root>]\n",
	"\n",
	"	[-c <reference levels file for Y's in factor file>]\n",
	"	[-R (pay attention to 'remaining' category)]\n",
	"	[-x <shorten category names, with separator in double quotes (default=\"\")>]\n",
	"\n",
	"This script will fit the following 2 models for each category, i = 1 to ", NUM_TOP_RESP_CAT, ":\n",
	"  where the top P=", NUM_TOP_PRED_CAT, " categories are included.\n",
	"\n",
	"	1.) ALR_Resp[i] = ALR_Pred[i] + ALR_Pred[top-p] + covariates\n",
	"\n",
	" 	2.) (ALR_Resp[i]-ALR_Pred[i]) = ALR_Pred[top-p] + covariates\n",
	"\n",
	"If the -R flag is set, a 'remaining' category will be be included in the denominator\n",
	"	independent of how large it is.  I.e., do not use it as one of the response variables.\n",
	"\n", sep="");

if(!length(opt$summary_file) || !length(opt$factors) || !length(opt$model_var) || 
	 !length(opt$response) || !length(opt$pairings) || !length(opt$factor_samp_id_name)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputRoot=gsub(".summary_table.xls", "", opt$summary_file);
	OutputRoot=gsub(".summary_table.tsv", "", OutputRoot);
}else{
	OutputRoot=opt$outputroot;
}

if(!length(opt$summary_file2)){
	SecondSummaryTable="";
}else{
	SecondSummaryTable=opt$summary_file2;
}

if(!length(opt$num_resp_var)){
	NumRespVariables=NUM_TOP_RESP_CAT;
}else{
	NumRespVariables=opt$num_resp_var;
}

if(!length(opt$num_pred_var)){
	NumPredVariables=NUM_TOP_PRED_CAT;
}else{
	NumPredVariables=opt$num_pred_var;
}

NumMaxALRVariables=max(NumRespVariables, NumPredVariables);

if(!length(opt$reference_levels)){
        ReferenceLevelsFile="";
}else{
        ReferenceLevelsFile=opt$reference_levels;
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

if(length(opt$required)){
	RequiredFile=opt$required;
}else{
	RequiredFile="";
}

SummaryFile=opt$summary_file;
FactorsFile=opt$factors;
ModelVarFile=opt$model_var;
PairingsFile=opt$pairings;
ResponseName=opt$response;
PredictorName=opt$predictor;
FactorSampleIDName=opt$factor_samp_id_name;

cat("\n");
cat("         Summary File: ", SummaryFile, "\n", sep="");
if(SecondSummaryTable!=""){
	cat("     2nd Summary File: ", SecondSummaryTable, "\n", sep="");
}
cat("         Factors File: ", FactorsFile, "\n", sep="");
cat("Factor Sample ID Name: ", FactorSampleIDName, "\n", sep="");
cat(" Model Variables File: ", ModelVarFile, "\n", sep="");
cat("        Pairings File: ", PairingsFile, "\n", sep="");
cat("        Response Name: ", ResponseName, "\n", sep="");
cat("       Predictor Name: ", PredictorName, "\n", sep="");
cat("          Output File: ", OutputRoot, "\n", sep="");
cat("\n");
cat("Number of Predictor Variables: ", NumPredVariables, "\n", sep="");
cat(" Number of Response Variables: ", NumRespVariables, "\n", sep="");
cat("           Max ALR Var to Fit: ", NumMaxALRVariables, "\n", sep="");
cat("\n");
cat("Reference Levels File: ", ReferenceLevelsFile, "\n", sep="");
cat("Use Remaining? ", UseRemaining, "\n");
cat("Shorten Category Names: '", ShortenCategoryNames, "'\n", sep="");
cat("\n");

if(ShortenCategoryNames==TRUE){
	cat("Error:  You need to specify a delimitor to split the category names.\n");
	cat("        i.e., this -x option is not a flag, it requires a parameter.\n");
	quit(status=-1);
}

options(width=100);
cat("Text Line Width: ", options()$width, "\n", sep="");

##############################################################################
##############################################################################

load_factors=function(fname){
	factors=data.frame(read.table(fname,  sep="\t", header=TRUE, row.names=1, check.names=FALSE));
	factor_names=colnames(factors);

	ignore_idx=grep("^IGNORE\\.", factor_names);

	if(length(ignore_idx)!=0){
		return(factors[-ignore_idx]);
	}else{
		return(factors);
	}
}

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

load_reference_levels_file=function(fname){
        inmat=as.matrix(read.table(fname, sep="\t", header=F, check.names=FALSE, comment.char="#", row.names=1))
        colnames(inmat)=c("ReferenceLevel");
        print(inmat);
        cat("\n");
        if(ncol(inmat)!=1){
                cat("Error reading in reference level file: ", fname, "\n");
                quit(status=-1);
        }
        return(inmat);
}

relevel_factors=function(factors, ref_lev_mat){
        num_factors_to_relevel=nrow(ref_lev_mat);
        relevel_names=rownames(ref_lev_mat);
	factor_names=colnames(factors);
        for(i in 1:num_factors_to_relevel){
		
		target_relev_name=relevel_names[i];
		if(any(target_relev_name==factor_names)){
			tmp=factors[,target_relev_name];
			#print(tmp);
			tmp=relevel(tmp, ref_lev_mat[i, 1]);
			#print(tmp);
			factors[,target_relev_name]=tmp;
		}else{
			cat("Note: ", target_relev_name, " not in model.  Ignoring reference releveling.\n\n", sep="");
		}
        }
        return(factors);
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

load_list=function(filename){
	val=scan(filename, what=character(), comment.char="#");
	return(val);
}

load_mapping=function(filename, src, dst){
	mapping=read.table(filename, sep="\t", header=T, comment.char="#", quote="", row.names=NULL);
	column_names=colnames(mapping);
	if(all(column_names!=src)){
		cat("Error: Could not find ", src, " in header of map file.\n");
		quit(status=-1);
	}
	if(all(column_names!=dst)){
		cat("Error: Could not find ", dst, " in header of map file.\n");
		quit(status=-1);
	}

	map=cbind(as.character(mapping[,src]), as.character(mapping[,dst]));
	colnames(map)=c(src, dst);
	return(map);
}

extract_top_categories=function(ordered_normalized, top){

	num_samples=nrow(ordered_normalized);
	num_categories=ncol(ordered_normalized);

	cat("Samples: ", num_samples, "\n");
	cat("Categories: ", num_categories, "\n");
	
	num_saved=min(c(num_categories, top+1));

	cat("Top Requested to Extract: ", top, "\n");
	cat("Columns to Extract: ", num_saved, "\n");

	top_cat=matrix(0, nrow=num_samples, ncol=num_saved);
	top=num_saved-1;

	# Extract top categories requested
	top_cat[,1:top]=ordered_normalized[,1:top];

	# Included remaineder as sum of remaining categories
	top_cat[,(top+1)]=apply(
		ordered_normalized[,(top+1):num_categories, drop=F],
		1, sum);

	rownames(top_cat)=rownames(ordered_normalized);
	colnames(top_cat)=c(colnames(ordered_normalized)[1:top], "Remaining");

	return(top_cat);
			
}

additive_log_rato=function(ordered_matrix){
# Assumes last column will be the denominator

	num_cat=ncol(ordered_matrix);
	num_samp=nrow(ordered_matrix);

	denominator=ordered_matrix[,num_cat];
	alr_mat=matrix(0, nrow=num_samp, ncol=(num_cat-1));
	
	for(i in 1:num_samp){
		alr_mat[i,]=log(ordered_matrix[i,1:(num_cat-1)]/denominator[i]);
		#print(alr_mat[i,]);
	}

	rownames(alr_mat)=rownames(ordered_matrix)
	colnames(alr_mat)=head(colnames(ordered_matrix), num_cat-1);

	alr_struct=list();
	alr_struct[["transformed"]]=alr_mat;
	alr_struct[["denominator"]]=denominator;

	return(alr_struct);
}

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

        num_row=nrow(mat);
        num_col=ncol(mat);

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

		mat=mat[row_dendr[["names"]], col_dendr[["names"]]];
		
	}else if(plot_col_dendr){
		layoutmat=matrix(
			c(
			rep(rep(2, heatmap_width), col_dend_height),
			rep(rep(1, heatmap_width), heatmap_height)
			), byrow=T, ncol=heatmap_width); 

		col_dendr=get_dendrogram(mat, type="col");
		mat=mat[, col_dendr[["names"]]];
		
	}else if(plot_row_dendr){
		layoutmat=matrix(
			rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height),
			byrow=T, ncol=row_dend_width+heatmap_width);

		row_dendr=get_dendrogram(mat, type="row");
		mat=mat[row_dendr[["names"]],];
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

##############################################################################
##############################################################################

# Open main output file
pdf(paste(OutputRoot, ".2smp_alr.pdf", sep=""), height=11, width=9.5);

# Load summary file table counts 
cat("\n");
cat("Loading summary table...\n");
counts1=load_summary_file(SummaryFile);

counts2=NULL;
if(SecondSummaryTable!=""){
	cat("Loading second summary table...\n");
	counts2=load_summary_file(SecondSummaryTable);	
	cat("Merging second summary table..\n");
	counts=merge_summary_tables(counts1, counts2);
	cat("Merged Summary Table Samples: ", nrow(counts), "\n", sep="");
	cat("Merged Summary Table Categories: ", ncol(counts), "\n", sep="");
}else{
	counts=counts1;
}

# Remove zero count samples
tot=apply(counts, 1, sum);
nonzero=tot>0;
if(!(all(nonzero))){
	cat("WARNING: Zero count samples found:\n");
	samp_names=rownames(counts);
	print(samp_names[!nonzero]);
	cat("\n");
	counts=counts[nonzero,,drop=F];
}

num_st_categories=ncol(counts);
num_st_samples=nrow(counts);

# Load resp/pred sample mappings
all_pairings_map=load_mapping(PairingsFile, ResponseName, PredictorName);
st_samples=rownames(counts);
num_pairings_loaded=nrow(all_pairings_map);
cat("\n");
cat("Pairing entries loaded: ", num_pairings_loaded, "\n");

available_pairings=matrix(F, nrow=nrow(all_pairings_map), ncol=2, dimnames=dimnames(all_pairings_map));
rownames(available_pairings)=all_pairings_map[,1];
available_pairings[intersect(st_samples, all_pairings_map[,1]),1]=T;
rownames(available_pairings)=all_pairings_map[,2];
available_pairings[intersect(st_samples, all_pairings_map[,2]),2]=T;

print(available_pairings);
num_resp_avail=sum(available_pairings[,ResponseName]);
num_pred_avail=sum(available_pairings[,PredictorName]);

# Extract complete pairs
pairs_complete=apply(available_pairings, 1, function(x){all(x)});
pairs_complete_samples=names(pairs_complete[pairs_complete]);
pairings_map=all_pairings_map[pairs_complete,];
num_complete_pairings=sum(pairs_complete);

pairs_incomplete=apply(available_pairings, 1, function(x){any(!x)});
pairs_incomplete_samples=names(pairs_complete[pairs_incomplete]);
pairings_incomplete_map=all_pairings_map[pairs_incomplete,];
num_incomplete_pairings=sum(pairs_incomplete);

cat("Number of complete pairings: ", num_complete_pairings, "\n");
#print(pairs_complete_samples);
paired_responses=pairings_map[,1];
paired_predictors=pairings_map[,2];
paired_samples=c(paired_responses, paired_predictors);
counts=counts[paired_samples,];
num_samples=nrow(counts);
num_categories=ncol(counts);
cat("Count Matrix:\n");
cat("  Num Samples: ", num_samples, "\n");
cat("  Num Categories: ", num_categories, "\n");
cat("\n");

loaded_sample_info=c(
	"Summary Table Info:",
	paste(" 1st Summary Table Name: ", SummaryFile, sep=""),
	paste("    Samples: ", nrow(counts1), " Categories: ", ncol(counts1), sep=""),
	paste(" 2nd Summary Table Name: ", SecondSummaryTable, sep=""),
	paste("    Samples: ", nrow(counts2), " Categories: ", ncol(counts2), sep=""),
	"",
	paste("  Total Number Samples Loaded: ", num_st_samples, sep=""),
	paste("  Total Number Categories Loaded: ", num_st_categories, sep=""),
	"",
	"Sample Pairing Info:",
	paste("  Mapping Name: ", PairingsFile, sep=""),
	paste("  Number of Possible Pairings Loaded: ", num_pairings_loaded, sep=""),
	paste("  Number of ", ResponseName, " Samples: ", num_resp_avail, sep=""),
	paste("  Number of ", PredictorName, " Samples: ", num_pred_avail, sep=""),
	"",
	paste("Number of Complete/Matched Pairings: ", num_complete_pairings, sep=""),
	paste("Number of InComplete/UnMatched Pairings: ", num_incomplete_pairings, sep="")
);

incomplete_pairing_info=c(
	"Incomplete Pairings: ",
	"",
	"Missing Sample Info Availability:",
	capture.output(print(available_pairings[pairs_incomplete,])),
	"",
	"Missing Sample Pairing:",
	capture.output(print(pairings_incomplete_map))
);

plot_text(loaded_sample_info);
plot_text(incomplete_pairing_info);

##############################################################################

# Shorten cateogry names
if(ShortenCategoryNames!=""){
	full_names=colnames(counts);
	splits=strsplit(full_names, ShortenCategoryNames);
	short_names=character();
	for(i in 1:length(full_names)){
		short_names[i]=tail(splits[[i]], 1);
		short_names[i]=gsub("_unclassified$", "_uncl", short_names[i]);
	}
	colnames(counts)=short_names;
	cat("Names have been shortened.\n");
}else{
	cat("Keeping original category names...\n");
}

# Normalize
cat("Normalizing counts...\n");
normalized=normalize(counts);

# Assign 0's to values smaller than smallest abundance across entire dataset
min_assay=min(normalized[normalized!=0]);
cat("Lowest non-zero value: ", min_assay, "\n", sep="");
zero_replacment=min_assay/10;
cat("Substituting 0's with: ", zero_replacment, "\n", sep="");
normalized[normalized==0]=zero_replacment;

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
	}
}else{
	cat("Assuming not categories called 'remainder' or 'remaining'\n");
}

# Reorder by abundance
cat("Reordering summary table categories by abundance...\n");
mean_abund=apply(normalized, 2, mean);
ix=order(mean_abund, decreasing=TRUE);
normalized=normalized[,ix];
mean_abund=mean_abund[ix];

if(UseRemaining){
	normalized=cbind(normalized, normalized_remaining_col_dat);
	mean_abund=c(mean_abund, mean(normalized_remaining_col_dat));
}

sorted_taxa_names=colnames(normalized);

num_top_taxa=NumMaxALRVariables;
num_top_taxa=min(c(num_top_taxa, num_st_categories));
prop_abundance_represented=sum(mean_abund[1:num_top_taxa]);

cat("\nThe top ", num_top_taxa, " taxa are:\n", sep="");
for(i in 1:num_top_taxa){
	cat("\t", sorted_taxa_names[i], "\t[", mean_abund[i], "]\n", sep="");
}
cat("\n");

cat("Accounting for ", prop_abundance_represented, " of taxa.\n", sep="");
cat("\n");

##############################################################################

# Load factors
cat("Loading Factors...\n");
factors=load_factors(FactorsFile);
factor_names=colnames(factors);
num_factors=ncol(factors);
factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);

num_loaded_factors=num_factors;
num_loaded_factor_samp=num_factor_samples;

cat("\n");
cat(num_factors, " Factor(s) Loaded:\n", sep="");
print(factor_names);
cat("\n");

# Load predictorss to include in model
if(ModelVarFile!=""){
	model_var_arr=load_list(ModelVarFile);
	cat("Model Variables:\n");
	print(model_var_arr);
	cat("\n");
}else if(ModelString!=""){
	cat("Error: Model String option features not implemented yet.\n");
	quit(status=-1);
}

# Load variables to require after NA removal
required_arr=NULL;
if(""!=RequiredFile){
	required_arr=load_list(RequiredFile);
	cat("Required Variables:\n");
	print(required_arr);
	cat("\n");
	missing_var=setdiff(required_arr, factor_names);
	if(length(missing_var)>0){
		cat("Error: Missing required variables from factor file:\n");
		print(missing_var);
	}
}else{
	cat("No Required Variables specified...\n");
}

plot_text(c(
	"Variables Targeted:",
	"",
	"Predictors:",
	capture.output(print(model_var_arr)),
	"",
	"Required Variables:",
	 capture.output(print(required_arr))
));

# Confirm we can find all the factors
missing_fact=setdiff(model_var_arr, factor_names);
if(length(missing_fact)>0){
	cat("Error: Factors in model, missing Factor File.\n");
	print(missing_fact);
	quit(status=-1);
}else{
	cat("All model variables found in factor file...\n");
}
factors=factors[,model_var_arr, drop=F];
factor_names=colnames(factors);
num_factors=ncol(factors);

# Relevel factor levels
if(ReferenceLevelsFile!=""){
        ref_lev_mat=load_reference_levels_file(ReferenceLevelsFile)
        factors=relevel_factors(kept_factors, ref_lev_mat);
}else{
        cat("No Reference Levels File specified.\n");
}

##############################################################################
# Reconcile factors with samples

cat("\nReconciling samples between factor file and (paired responses) summary table...\n");
factor_sample_ids=rownames(factors);

counts_sample_ids=pairings_map[,FactorSampleIDName];

shared_sample_ids=sort(intersect(factor_sample_ids, counts_sample_ids));

num_shared_sample_ids=length(shared_sample_ids);
num_factor_sample_ids=length(factor_sample_ids);
num_counts_sample_ids=length(counts_sample_ids);

cat("Num counts (paired responses) sample IDs: ", num_counts_sample_ids, "\n");
cat("Num factor sample IDs: ", num_factor_sample_ids, "\n");
cat("Num shared sample IDs: ", num_shared_sample_ids, "\n");
cat("\n");

cat("Samples missing from factor information:\n");
samples_missing_factor_info=(setdiff(counts_sample_ids, factor_sample_ids));
num_samp_missing_fact_info=length(samples_missing_factor_info);
cat("\n");
cat("Total samples shared: ", num_shared_sample_ids, "\n");

# Remove samples not in summary table
rownames(pairings_map)=pairings_map[,FactorSampleIDName];
pairings_map=pairings_map[shared_sample_ids,];

paired_responses=pairings_map[,1];
paired_predictors=pairings_map[,2];


paired_samples=sort(c(paired_responses, paired_predictors));

# Reorder data by sample id
normalized=normalized[paired_samples,];
num_samples=nrow(normalized);
recon_factors=factors[shared_sample_ids,,drop=F];

factor_file_info=c(
	paste("Factor File Name: ", FactorsFile, sep=""),
	"",
	paste("Num Loaded Factors/Variables: ", num_loaded_factors, sep=""),
	paste("Num Samples in Factor File: ", num_loaded_factor_samp, sep=""),
	"",
	paste("Num Samples Shared between Factors and Pairable Samples: ", num_shared_sample_ids, sep=""),
	"",
	paste("Num Samples Missing Factor Information: ", num_samp_missing_fact_info, sep=""),
	"",
	"Samples missing info: ",
	capture.output(print(samples_missing_factor_info))
);

plot_text(factor_file_info);

##############################################################################
# Remove samples with NAs

cat("Identifying samples/factors to keep with NAs...\n");
num_samples_recon=nrow(recon_factors);
num_factors_recon=ncol(recon_factors);
num_samples_before_na_removal=num_samples_recon;
num_factors_before_na_removal=num_factors_recon;
factors_wo_nas_res=remove_sample_or_factors_wNA_parallel(recon_factors, 
	required=required_arr, num_trials=64000, num_cores=64, outfile=paste(OutputRoot, ".noNAs", sep=""));

factors_wo_nas=factors_wo_nas_res$factors;
factor_names_wo_nas=colnames(factors_wo_nas);
factor_sample_ids_wo_nas=rownames(factors_wo_nas);
model_var_arr=intersect(model_var_arr, factor_names_wo_nas);

# Subset pairing map based on factor sample IDs
rownames(pairings_map)=pairings_map[,FactorSampleIDName];
pairings_map=pairings_map[factor_sample_ids_wo_nas,];
paired_responses=pairings_map[,1];
paired_predictors=pairings_map[,2];
paired_samples=sort(c(paired_responses, paired_predictors));

# Subset the normalized counts based on pairing map
normalized=normalized[paired_samples,, drop=F];
num_samples_wo_nas=nrow(factors_wo_nas);
num_factors_wo_nas=ncol(factors_wo_nas);

#cat("Num Samples w/o NAs: ", num_samples_wo_nas, "\n");
#cat("Num Factors w/o NAs: ", num_factors_wo_nas, "\n");
#cat("\n");

##############################################################################
# Prepping for ALR calculations

if(NumMaxALRVariables >= num_st_categories){
	NumMaxALRVariables = (num_st_categories-1);
	cat("Number of taxa to work on was changed to: ", NumMaxALRVariables, "\n");
	
	NumRespVariables=min(NumMaxALRVariables, NumRespVariables);
	NumPredVariables=min(NumMaxALRVariables, NumPredVariables);

	cat("Number of ALR Predictors to include: ", NumPredVariables, "\n");
	cat("Number of ALR Responses to analyze: ", NumRespVariables, "\n");
}

##############################################################################

cat("\n");
cat("Extracting Top categories: ", NumMaxALRVariables, " from amongst ", ncol(normalized), "\n", sep="");

cat_abundances=extract_top_categories(normalized, NumMaxALRVariables);
resp_alr_struct=additive_log_rato(cat_abundances);
alr_categories_val=resp_alr_struct$transformed;
alr_cat_names=colnames(alr_categories_val);

plot_text(c(
	paste("Num (Reconciled) Samples before NA removal: ", num_samples_before_na_removal, sep=""),
	paste("Num Factors before NA removal: ", num_factors_before_na_removal, sep=""),
	"",
	"Acceptable Variables after NA Removal:",
	"",
	capture.output(print(factor_names_wo_nas)),
	"",
	paste("Num Samples w/o NAs: ", num_samples_wo_nas, sep=""),
	paste("Num Factors w/o NAs: ", num_factors_wo_nas, sep="")
));

plot_text(c(
	paste("ALR Categories (Top ", NumMaxALRVariables, ")", sep=""),
	capture.output(print(alr_cat_names))
));


cat("\n");
cat("ALR Category Summary:\n");
s=summary(alr_categories_val);
print(s);
plot_text(c(
	"ALR Categories Summary:",
	"\n",
	capture.output(print(s))
));
plot_histograms(alr_categories_val);

##############################################################################

# Order/Split the data....

# Order the pairings map by response sample IDs
rownames(pairings_map)=pairings_map[,ResponseName];
sorted_response_ids=sort(rownames(pairings_map));
pairings_map=pairings_map[sorted_response_ids,];

# Get the ordered response and predictor sample IDs
response_sample_ids=pairings_map[,ResponseName];
predictor_sample_ids=pairings_map[,PredictorName];

# Extract the response/predictor ALR and factors values in the right order
response_alr=alr_categories_val[response_sample_ids,,drop=F];
predictor_alr=alr_categories_val[predictor_sample_ids,,drop=F];
factors=factors_wo_nas[pairings_map[,FactorSampleIDName],];

##############################################################################

alr_names=colnames(alr_categories_val);
plots_per_page=6
par(mfrow=c(plots_per_page,2));
par(oma=c(0,0,0,0));
par(mar=c(4,4,2,2));
median_delta=numeric(NumRespVariables);
names(median_delta)=alr_names;
for(cat_ix in 1:NumRespVariables){

	cur_alr_resp=response_alr[,cat_ix];
	cur_alr_pred=predictor_alr[,cat_ix];


	if((cat_ix %% plots_per_page)==0){
		bottom_label=PredictorName;
	}else{
		bottom_label="";
	}

	# Plot matching pred/resp off reference line
	plot(cur_alr_pred, cur_alr_resp, type="n", main=alr_names[cat_ix],
		xlab=bottom_label, ylab=ResponseName);
	abline(0, 1, col="blue");
	points(cur_alr_pred, cur_alr_resp, main=alr_names[cat_ix]);

	deltas=cur_alr_resp-cur_alr_pred;
	median_delta[cat_ix]=median(deltas);
	hist(deltas, main=paste(alr_names[cat_ix], " ALR: ", ResponseName, " - ", PredictorName, sep=""), 
		xlab="");
	
}

par(mfrow=c(2,1));
par(mar=c(20,4,3,1));
color_arr=rainbow(NumRespVariables, start=0, end=4/6);
barplot(median_delta, horiz=F, las=2, 
	main=paste("Median Difference between ", ResponseName, " and ", PredictorName, sep=""),
	col=color_arr);

dec_del_ix=order(median_delta, decreasing=F);
barplot(median_delta[dec_del_ix], horiz=F, las=2, 
	main=paste("Median Difference between ", ResponseName, " and ", PredictorName, sep=""),
	col=color_arr[dec_del_ix]
	);

par(mfrow=c(1,1));
#print(response_alr);
#print(predictor_alr);
#print(factors);


# Store the results
num_model_pred=length(model_var_arr);

cov_model_string=paste(model_var_arr, collapse="+");
mmat=model.matrix(as.formula(paste("~", cov_model_string, "-1")), data=as.data.frame(factors));
cov_coeff_names=colnames(mmat);
num_cov_coeff_names=length(cov_coeff_names);

cat("Anticipated Coefficient Names:\n");
print(cov_coeff_names);
cat("\n");

category_alr_coef_mat=matrix(NA, nrow=NumRespVariables, ncol=NumRespVariables,
	dimnames=list(alr_cat_names, alr_cat_names));
category_alr_pval_mat=matrix(NA, nrow=NumRespVariables, ncol=NumRespVariables,
	dimnames=list(alr_cat_names, alr_cat_names));

covariates_coef_mat  =matrix(NA, nrow=NumRespVariables, ncol=num_cov_coeff_names, 
	dimnames=list(alr_cat_names, cov_coeff_names));
covariates_pval_mat  =matrix(NA, nrow=NumRespVariables, ncol=num_cov_coeff_names, 
	dimnames=list(alr_cat_names, cov_coeff_names));

rsqrd_mat             =matrix(NA, nrow=NumRespVariables, ncol=4, 
	dimnames=list(alr_cat_names[1:NumRespVariables], c("R^2", "Adj-R^2", "Reduced Adj-R^2", "ALR Contrib")));

model_pval_mat             =matrix(NA, nrow=NumRespVariables, ncol=1, 
	dimnames=list(alr_cat_names[1:NumRespVariables], c("F-statistic P-value")));

# Fit the regression model

for(resp_ix in 1:NumRespVariables){
	alr_resp=response_alr[,resp_ix,drop=F];
	resp_cat_name=colnames(alr_resp);
	alr_resp=as.vector(alr_resp);
	
	cat("Fitting: ", resp_ix, ".) ", resp_cat_name, "\n");

	if(resp_ix<=NumPredVariables){
		# Response category is already in predictors
		alr_pred=predictor_alr[,1:NumPredVariables];
	}else{
		# Include response category in predictor
		alr_pred=cbind(predictor_alr[,resp_cat_name, drop=F], predictor_alr[,1:NumPredVariables]);
	}
	alr_pred_names=colnames(alr_pred);
	
	model_pred_df=as.data.frame(cbind(alr_pred, factors));

	#print(resp);
	#print(model_pred_df);
	model_str=paste("alr_resp ~ ", paste(c(alr_pred_names,model_var_arr), collapse="+"), sep="");
	model_reduced_str=paste("alr_resp ~ ", paste(model_var_arr, collapse="+"), sep="");

	cat("Model String: \n");
	cat("Full:\n");
	print(model_str);
	cat("Reduced:\n");
	print(model_reduced_str);

	lm_fit=lm(as.formula(model_str), data=model_pred_df);
	lm_reduced_fit=lm(as.formula(model_reduced_str), data=model_pred_df);

	sum_fit=summary(lm_fit);
	sum_reduced_fit=summary(lm_reduced_fit);

	anova_res=anova(lm_fit);

	plot_text(c(
		paste(resp_ix, ".) ", resp_cat_name, ":", sep=""),
		"",
		capture.output(print(sum_fit))
		)
	);

	plot_text(c(
		paste(resp_ix, ".) ", resp_cat_name, ":", sep=""),
                "",
		capture.output(print(anova_res))
		)
	);

	mmps(lm_fit);

	model_coef_names=setdiff(rownames(sum_fit$coefficients), "(Intercept)");

	cat_names=intersect(model_coef_names, alr_pred_names);
	category_alr_coef_mat[resp_cat_name, cat_names]=sum_fit$coefficients[cat_names,"Estimate"];
	category_alr_pval_mat[resp_cat_name, cat_names]=sum_fit$coefficients[cat_names,"Pr(>|t|)"];

	cat_names=intersect(model_coef_names, cov_coeff_names);
	covariates_coef_mat[resp_cat_name, cat_names]=sum_fit$coefficients[cat_names,"Estimate"];
	covariates_pval_mat[resp_cat_name, cat_names]=sum_fit$coefficients[cat_names,"Pr(>|t|)"];

	rsqrd_mat[resp_cat_name, "R^2"]=sum_fit$r.squared;
	rsqrd_mat[resp_cat_name, "Adj-R^2"]=sum_fit$adj.r.squared;
	rsqrd_mat[resp_cat_name, "Reduced Adj-R^2"]=sum_reduced_fit$adj.r.squared;
	rsqrd_mat[resp_cat_name, "ALR Contrib"]=sum_fit$adj.r.squared-sum_reduced_fit$adj.r.squared;

	model_pval_mat[resp_cat_name, "F-statistic P-value"]=1-pf(sum_fit$fstatistic[1], sum_fit$fstatistic[2], sum_fit$fstatistic[3]);

	cat("*************************************************\n");

}

#print(category_alr_coef_mat);
#print(category_alr_pval_mat);

#print(covariates_coef_mat);
#print(covariates_pval_mat);

all.nas=apply(covariates_coef_mat, 2, function(x){all(is.na(x))});
covariates_coef_mat=covariates_coef_mat[,!all.nas,drop=F];

all.nas=apply(covariates_pval_mat, 2, function(x){all(is.na(x))});
covariates_pval_mat=covariates_pval_mat[,!all.nas,drop=F];

print(rsqrd_mat);

#paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, deci_pts=4,
#       label_zeros=T, counts=F, value.cex=1,
#        plot_col_dendr=F,
#        plot_row_dendr=F

par(oma=c(2,1,5,2));

# ALR Coefficients
paint_matrix(category_alr_coef_mat, 
	title=paste("Top ", NumPredVariables, " ", PredictorName," Predictor ALR Coefficients for Top ", NumRespVariables, " ", ResponseName, " Responses ALR", sep=""), 
	deci_pts=2, value.cex=.8);
mtext(PredictorName, side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# ALR Coefficients Clustered
paint_matrix(category_alr_coef_mat[,1:NumPredVariables], 
	title=paste("Top ", NumPredVariables, " ", PredictorName," Predictor ALR Coefficients for Top ", NumRespVariables, " ", ResponseName, " Responses ALR", sep=""), 
	 plot_col_dendr=T, plot_row_dendr=T,
	deci_pts=2, value.cex=.8);
mtext(PredictorName, side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# ALR P-Values
paint_matrix(category_alr_pval_mat, plot_min=0, plot_max=1,
	title=paste("Top ", NumPredVariables, " ", PredictorName," Predictor ALR P-Values for Top ", NumRespVariables, " ", ResponseName, " Responses ALR", sep=""), 
	high_is_hot=F, deci_pts=2, value.cex=.8);
mtext(PredictorName, side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# ALR Coefficients w/ P-value masked
limit=max(abs(range(category_alr_coef_mat)));
category_alr_coef_masked_mat=mask_matrix(
	val_mat=category_alr_coef_mat, 
	mask_mat=category_alr_pval_mat, 
	mask_thres=0.05, 
	mask_val=0.0);
paint_matrix(category_alr_coef_masked_mat,
	title=paste("Top ", NumPredVariables, " ", 
		PredictorName," Predictor ALR Coeff for Top ", NumRespVariables, " ", ResponseName,
		" Responses ALR P-Val(<.05) Maskd", sep=""), 
	label_zeros=F, high_is_hot=F, deci_pts=2, value.cex=.8);
mtext(PredictorName, side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# ALR P-Values Clustered
paint_matrix(category_alr_pval_mat[,1:NumPredVariables], plot_min=0, plot_max=1,
	title=paste("Top ", NumPredVariables, " ", PredictorName," Predictor ALR P-Values for Top ", NumRespVariables, " ", ResponseName, " Responses ALR", sep=""), 
	plot_col_dendr=T, plot_row_dendr=T,
	high_is_hot=F, deci_pts=2, value.cex=.8);
mtext(PredictorName, side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# Covariate Coefficients
paint_matrix(covariates_coef_mat, 
	title=paste("Covariates Coefficients for Top ", NumRespVariables, " ", ResponseName, " Categories", sep=""),
	deci_pts=2, value.cex=.8);
mtext("Covariates", side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# Covariate Coefficients w/ Clustering
paint_matrix(covariates_coef_mat, 
	title=paste("Covariates Coefficients for Top ", NumRespVariables, " ", ResponseName, " Categories", sep=""),
	plot_col_dendr=T, plot_row_dendr=T,
	deci_pts=2, value.cex=.8);
mtext("Covariates", side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# Covariate P-Values
paint_matrix(covariates_pval_mat, plot_min=0, plot_max=1, 
	title=paste("Covariates P-Values for Top ", NumRespVariables, " ", ResponseName, " Categories", sep=""),
	high_is_hot=F, deci_pts=2, value.cex=.8);
mtext("Covariates", side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# Covariate Coefficients w/ P-value masked
covariates_coef_masked_mat=mask_matrix(
        val_mat=covariates_coef_mat,
        mask_mat=covariates_pval_mat,
        mask_thres=0.05,
        mask_val=0.0);
paint_matrix(covariates_coef_masked_mat,,
        title=paste("Top ", NumPredVariables, " ",
                PredictorName," Covariates Coeff for Top ", NumRespVariables, " ", ResponseName,
                " Responses ALR P-Val(<.05) Maskd", sep=""),
        label_zeros=F, high_is_hot=F, deci_pts=2, value.cex=.8);
mtext("Covariates", side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);


# Covariate P-Values Clustered
paint_matrix(covariates_pval_mat, plot_min=0, plot_max=1, 
	title=paste("Covariates P-Values for Top ", NumRespVariables, " ", ResponseName, " Categories", sep=""),
	plot_col_dendr=T, plot_row_dendr=T,
	high_is_hot=F, deci_pts=2, value.cex=.8);
mtext("Covariates", side=1, cex=2, font=2, line=.75);
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# R^2
paint_matrix(rsqrd_mat, plot_min=0, plot_max=1, title=paste("Explained Variation for Top ", NumRespVariables, " Responses: R^2", sep=""));
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# Model pval
paint_matrix(model_pval_mat, plot_min=0, plot_max=1, high_is_hot=F, 
	title=paste("Model F-Statistic P-Values for Top ", NumRespVariables, sep=""));
mtext(ResponseName, side=4, cex=2, font=2, line=.75);

# Plot Pred->Resp for ALR
plot_pred_resp_bar=function(coef_mat, pval_mat, title=""){
	val=diag(coef_mat);
	pval=diag(pval_mat);
	cat_names=rownames(coef_mat);

	ord=order(val, decreasing=T);

	val=val[ord];
	pval=pval[ord];
	cat_names=cat_names[ord];

        num_colors=50;
        color_arr=rainbow(num_colors, start=0, end=4/6);

        # Provide a means to map values to an (color) index
        remap=function(in_val, in_range, out_range){
                in_prop=(in_val-in_range[1])/(in_range[2]-in_range[1])
                out_val=in_prop*(out_range[2]-out_range[1])+out_range[1];
                return(out_val);
        }

	colors=ceiling(remap(pval, c(0,1), c(1,num_colors)));

	barplot(val, names.arg=cat_names, las=2, col=color_arr[colors],
		main=title
		);
}

par(oma=c(20, 5, 5, 1));
plot_pred_resp_bar(category_alr_coef_mat, category_alr_pval_mat);
mtext(paste("Predictability of ", ResponseName, " ALR based on ", PredictorName, " ALR", sep=""),
	side=3, font=2, line=2
	);
mtext("After Controlling for Covariates", side=3, font=2, line=1);
mtext("Regression Coefficient", side=2, font=1, line=3.75);


cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
