#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
library(car);

source('~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r');

options(useFancyQuotes=F);


params=c(
	"summary_file", "s", 1, "character",

	"factors", "f", 1, "character",
	"factor_samp_id_name", "F", 2, "character",
	"factor_subj_id_name", "S", 2, "character",
	"model_var", "M", 1, "character",
	"required", "q", 2, "character",

	"pairings", "p", 1, "character",
	"B_minuend", "B", 1, "character",
	"A_subtrahend", "A", 1, "character",

	"num_top", "u", 2, "numeric",
	"alr_list_file", "a", 2, "character",

	"outputroot", "o", 2, "character",

	"reference_levels", "c", 2, "character",
	"contains_remaining", "R", 2, "logical",
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
	"	  -F <Sample ID Column Name>\n",
	"	  -S <Subject ID Column Name>\n",
	"	-M <list of covariate X's names to include in the model from the factor file>\n",
	"	[-q <required list of variables to include after NA removal>]\n",
	"\n",
	"	-p <pairings map, pairing sample IDs from two groups. Must have header/column names>\n",
	"	-B <Sample Group B (minuend: B-A=diff), (column name in pairings file)\n",
	"	-A <Sample Group A (subtrahend: B-A=diff), (column name in pairings file)\n",
	"\n",
	"	ALR Options:\n",
	"	[-u <number of top ALR categories to analyze, default=", NUM_TOP_RESP_CAT, "]\n",
	"	[-a <list of ALR categories to use in additon to top>]\n",
	"\n",
	"	[-o <output filename root>]\n",
	"\n",
	"	[-c <reference levels file for Y's in factor file>]\n",
	"	[-x <shorten category names, with separator in double quotes (default=\"\")>]\n",
	"\n",
	"	[-t <tag name>]\n",
	"\n",
	"This script will fit the following model for each ALR category from the top ", NUM_TOP_RESP_CAT, " taxa\n",
	"  and/or the specified list.\n",
	"\n",
	"	1.) (ALR_B_minuend[i] - ALR_A_subtrahend[i]) = covariates\n",
	"\n",
	"If the -R flag is set, a 'remaining' category will be be included in the denominator\n",
	"	independent of how large it is.  I.e., do not use it as one of the response variables.\n",
	"\n", sep="");

if(
	!length(opt$summary_file) || 
	!length(opt$factors) || 
	!length(opt$model_var) || 
	!length(opt$A_subtrahend) || 
	!length(opt$B_minuend) || 
	!length(opt$pairings)
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

if(!length(opt$num_top)){
	NumTopALR=NUM_TOP_RESP_CAT;
}else{
	NumTopALR=opt$num_top;
}

if(!length(opt$reference_levels)){
        ReferenceLevelsFile="";
}else{
        ReferenceLevelsFile=opt$reference_levels;
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

if(length(opt$alr_list_file)){
	ALRCategListFile=opt$alr_list_file;
}else{
	ALRCategListFile="";
}

if(length(opt$factor_samp_id_name)){
	FactorSampleIDName=opt$factor_samp_id_name;
	if(FactorSampleIDName==TRUE){
		FactorSampleIDName="";
	}
}else{
	FactorSampleIDName="";
}

if(length(opt$factor_subj_id_name)){
	FactorSubjectIDName=opt$factor_subj_id_name;
	if(FactorSubjectIDName==TRUE){
		FactorSubjectIDName="";
	}
}else{
	FactorSubjectIDName="";
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


OutputRoot=paste(OutputRoot, ".a_", A_subtrahend, ".b_", B_minuend, sep="");

cat("\n");
cat("         Summary File: ", SummaryFile, "\n", sep="");
cat("         Factors File: ", FactorsFile, "\n", sep="");
cat("Factor Sample ID Name: ", FactorSampleIDName, "\n", sep="");
cat(" Model Variables File: ", ModelVarFile, "\n", sep="");
cat("        Pairings File: ", PairingsFile, "\n", sep="");
cat("            A Minuend: ", A_subtrahend, "\n", sep="");
cat("         B Subtrahend: ", B_minuend, "\n", sep="");
cat("          Output File: ", OutputRoot, "\n", sep="");
cat("\n");
cat("List of ALR Categories to Include (instead of using Top): ", ALRCategListFile, "\n", sep="");
cat("\n");
cat("Reference Levels File: ", ReferenceLevelsFile, "\n", sep="");
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

plot_text=function(strings, size_multi=1){
	par(mfrow=c(1,1));
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
		text(0, top-i, strings[i], pos=4, cex=text_size*size_multi); 
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

signf_as_table=function(coef_mat, pval_mat){

        num_rows=nrow(coef_mat);
        num_cols=ncol(coef_mat);
        num_entries=num_rows*num_cols;

        cnames=colnames(coef_mat);
        rnames=rownames(coef_mat);
        tab_headers=c("Row", "Column", "Coefficient", "P-value", "Formatted");
        comb_tab=matrix(NA, nrow=num_entries, ncol=length(tab_headers));
        colnames(comb_tab)=tab_headers;

        pval_val=numeric(num_entries);

        line_ix=1;
        for(i in 1:num_rows){
                for(j in 1:num_cols){

                        pval=pval_mat[i,j];
                        if(!is.na(pval) && !is.nan(pval) && pval < 0.10){
                                comb_tab[line_ix, "Row"]=rnames[i];
                                comb_tab[line_ix, "Column"]=cnames[j];
                                comb_tab[line_ix, "Coefficient"]=sprintf("%.5g", coef_mat[i,j]);
                                comb_tab[line_ix, "P-value"]=sprintf("%.5g", pval);

                                comb_tab[line_ix, "Formatted"]=
                                        sprintf("(coef = %.4f, p-val = %.3g)", coef_mat[i,j], pval);
                                pval_val[line_ix]=pval;
                                line_ix=line_ix+1;
                        }
                }
        }

        num_signf=line_ix-1;

        if(num_signf>=1){

                comb_tab=comb_tab[1:num_signf,,drop=F];
                pval_val=pval_val[1:num_signf];

                sorted_ix=order(pval_val);

                comb_tab=comb_tab[sorted_ix,,drop=F];
                rownames(comb_tab)=1:num_signf;
        }

        return(comb_tab);

}

##############################################################################
##############################################################################

# Open main output file
pdf(paste(OutputRoot, ".paird_diff_alr.pdf", sep=""), height=11, width=9.5);

input_files=load_and_reconcile_files(
                sumtab=list(fn=SummaryFile, shorten_cat_names_char=ShortenCategoryNames,
                        return_top=NumTopALR, specific_cat_fn=ALRCategListFile),
		factors=list(fn=FactorsFile, sbj_cname=FactorSubjectIDName,
                        smp_cname=FactorSampleIDName),
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

keyed=F;

if(FactorSubjectIDName!=""){

        subject_ids=as.character(factors[,FactorSubjectIDName]);

        dup_sbj_id=duplicated(subject_ids);
        if(any(dup_sbj_id)){
                cat("Duplicated Subject IDs found.  Can not key off of FactorSubjectIDName.\n");
                print(subject_ids[dup_sbj_id]);
        }else{
                rownames(factors)=factors[,FactorSubjectIDName,];
                keyed=T;
        }

}

# Assume using factors from A_subtrahend

if(!keyed){
        if(FactorSampleIDName!=""){
                rownames(factors)=factors[,FactorSampleIDName];
                factors=factors[good_pairs_map[,A_subtrahend],];
        }
}

write_file_report(input_files[["Report"]]);

##############################################################################
##############################################################################

# Normalize
cat("Normalizing counts...\n");
counts=counts+.5;
normalized=normalize(counts);

resp_alr_struct=additive_log_rato(normalized);
alr_categories_val=resp_alr_struct$transformed;
alr_cat_names=colnames(alr_categories_val);
num_used_alr_cat=length(alr_cat_names);

NumRespVariables=ncol(alr_categories_val);

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
A_sample_ids=good_pairs_map[,A_subtrahend];
B_sample_ids=good_pairs_map[,B_minuend];
subject_ids=rownames(good_pairs_map);

# Extract the predictor ALR and factors values in the right order
A_alr_values=alr_categories_val[A_sample_ids,,drop=F];
B_alr_values=alr_categories_val[B_sample_ids,,drop=F];

##############################################################################

num_model_pred=length(model_var_arr);

cov_model_string=paste(model_var_arr, collapse="+");
mmat=model.matrix(as.formula(paste("~", cov_model_string, "-1")), data=as.data.frame(factors));
cov_coeff_names=c("(Intercept)", colnames(mmat));
num_cov_coeff_names=length(cov_coeff_names);

cat("Anticipated Coefficient Names:\n");
print(cov_coeff_names);
cat("\n");

wilcoxon_pval_mat=matrix(NA, nrow=num_used_alr_cat, ncol=1,
	dimnames=list(alr_cat_names, "Wilcoxon Paired"));

covariates_coef_mat  =matrix(NA, nrow=num_used_alr_cat, ncol=num_cov_coeff_names, 
	dimnames=list(alr_cat_names, cov_coeff_names));
covariates_pval_mat  =matrix(NA, nrow=num_used_alr_cat, ncol=num_cov_coeff_names, 
	dimnames=list(alr_cat_names, cov_coeff_names));

rsqrd_mat             =matrix(NA, nrow=num_used_alr_cat, ncol=2, 
	dimnames=list(alr_cat_names[1:num_used_alr_cat], c("R^2", "Adj-R^2")));

model_pval_mat             =matrix(NA, nrow=num_used_alr_cat, ncol=1, 
	dimnames=list(alr_cat_names[1:num_used_alr_cat], c("F-statistic P-value")));

###############################################################################

shorten_paired_names=function(namea, nameb){

	all_names=c(namea, nameb);
	cat("Original:\n");
	print(all_names);
	name_len=nchar(all_names);
	shortest_name_len=min(name_len);

	# Trim Shared Beginning
	for(i in 1:shortest_name_len){
		prefix=substr(all_names, 1, i);
		if(length(unique(prefix))>1){
			break;	
		}
	}
	
	all_names=substr(all_names, i, name_len);
	cat("Front Trimmed:\n");
	print(all_names);

	# Trim Shared Ending
	name_len=nchar(all_names);
	shortest_name_len=min(name_len);
	for(i in 1:shortest_name_len){
		postfix=substr(all_names, name_len-i+1, name_len);
		if(length(unique(postfix))>1){
			break;	
		}
	}

	all_names=substr(all_names, 1, name_len-i+1);
	cat("End Trimmed:\n");
	print(all_names);

	res=list();
	res[["a"]]=all_names[1:length(namea)];
	res[["b"]]=all_names[(length(namea)+1):length(all_names)];
	return(res);
}


adjust_label_pos=function(x, y){
	
	gp=par();

	xlimits=gp[["usr"]][c(1,2)];
	ylimits=gp[["usr"]][c(3,4)];

	text_h=strheight("X", units="user");
	text_w=strwidth("XXX", units="user");

	num_samp=length(x);

	cat("Min Y Spacing: ", text_h, "\n");

	#----------------------------------------------------------------------

	y_orig_names=rownames(y);
	y=y[order(y),];

	mid=round(num_samp/2);

	y[1:mid]=y[1:mid]-text_h;
	y[(mid+1):num_samp]=y[(mid+1):num_samp]+text_h;
	
	max_iterations=num_samp^2;

	for(it in 1:max_iterations){

		adj=F;
		# Adjust down
		for(i in 1:(num_samp-1)){
			if((y[i+1]-y[i]) < text_h){
				y[i]=y[i]-(text_h/2);
				#cat("adj down\n");

				# If getting pushed off bottom of plot, come back a little
				if(y[i]<ylimits[1]){
					y[i]=y[i]+(text_h/4);
					cat("pushed off bottom\n");
				}
				adj=T;
			}
		}

		# Adjust up
		for(i in 1:(num_samp-1)){
			if((y[i+1]-y[i]) < text_h){
				y[i+1]=y[i+1]+(text_h/2);
				#cat("adj up\n");
				
				# If getting pushed off top of plot, come back a little
				if(y[i+1]>ylimits[2]){
					y[i+1]=y[i+1]-(text_h/4);
					cat("pushed off top\n");
				}
				adj=T;
			}
		}
		
		if(adj==F){
			cat("Label adjustments converged.\n");
			break;
		}

	}
	if(it==max_iterations){
		cat("Label adjustments not converged yet.\n");
	}else{
		cat("Num Iterations taken: ", it, "\n");
	}

	y=y[y_orig_names];

	#----------------------------------------------------------------------

	x_orig_names=rownames(x);
	x=x[order(x),];
	
	for(i in 1:num_samp){
		x[i]=x[i]+(i-num_samp/2)/(num_samp/2)*text_w;
	}

	x=x[x_orig_names];

	#----------------------------------------------------------------------

	res=list();
	res[["x"]]=x;
	res[["y"]]=y;
	return(res);

}

plot_ab_comparisons=function(a, b, aname, bname, pval, title){
# Plot histogram of differences and paired samples
	layout_mat=matrix(c(
		1,2,
		3,3,
		4,4
	), ncol=2, byrow=T);
	layout(layout_mat);
	par(oma=c(1,1,3,1));
	par(mar=c(6,6,8,.25));
	par(family="serif");

	ranges=range(c(a,b));
	ab=c(a,b);
	ab_hrec=hist(ab, plot=F);

	# Plot individual A/B
	ahrec=hist(A_alr_val, breaks=ab_hrec$breaks, plot=F);
	bhrec=hist(B_alr_val, breaks=ab_hrec$breaks, plot=F);
	max_count=max(c(ahrec$counts, bhrec$counts));
	ahrec=hist(A_alr_val, breaks=ab_hrec$breaks, main=aname, ylim=c(0, max_count), 
		cex.axis=1.5, cex.lab=1.5, cex.main=1.5, xlab="ALR", las=1);
	bhrec=hist(B_alr_val, breaks=ab_hrec$breaks, main=bname, ylim=c(0, max_count), 
		cex.axis=1.5, cex.lab=1.5, cex.main=1.5, xlab="ALR", las=1);

	# Plot difference
	alr_dif=b-a;
	magn=max(abs(alr_dif));
	diff_plot_range=c(-magn*1.2, magn*1.2);
	print(diff_plot_range);
	hist(alr_dif, xlim=diff_plot_range, yaxt="n", xaxt="n", main="", 
		xlab="", ylab="", las=1, breaks=30);

	yaxis_info=par()$yaxp;
	y_ats=as.integer(seq(yaxis_info[1], yaxis_info[2], 
		length.out=max(min(yaxis_info[2], yaxis_info[3]), 3)));
	axis(2, at=y_ats, labels=y_ats, las=1, cex.axis=2.5, padj=.5);

	xaxis_info=par()$xaxp;
	x_ats=as.integer(seq(xaxis_info[1], xaxis_info[2], 
		length.out=max(min(xaxis_info[2], xaxis_info[3]), 5)));
	axis(1, at=x_ats, labels=x_ats, las=1, cex.axis=2.5, padj=1);

	title(main=gsub("^\\d+\\.\\)", "", title), line=5, cex.main=5);
	title(main=paste("Wilcoxon paired p-value: ", round(pval, 4), sep=""), line=2.5, cex.main=2);
	title(xlab=paste("ALR difference = ", bname, " - ", aname, sep=""), cex.lab=2.5, line=5);
	title(ylab="Frequency", cex.lab=2.5, line=4);

	abline(v=0, col="blue", lty=2, lwd=2);

	#----------------------------------------------------------------------
	# Plot scatter
	plot(a,b, xlim=ranges, ylim=ranges, xlab=aname, ylab=bname, 
		col="red",
		cex=1.2, cex.axis=1.5, cex.lab=1.5, las=1);
	abline(a=0, b=1, col="blue", lty=2);
	mtext(title, side=3, line=0, outer=T, font=2, cex=2);

	
	#----------------------------------------------------------------------
	# Plot scatter with labels
	par(mfrow=c(1,1));
	short_names=shorten_paired_names(rownames(a), rownames(b));

	expand_range=function(r, prop){
		span=(r[2]-r[1])
		pad=span*prop;
		return(c(r[1]-pad, r[2]+pad));
	}

	adj_lab_pos=adjust_label_pos(a,b);
	
	exp_rng=expand_range(ranges, .15);
	plot(a,b, xlim=exp_rng, ylim=exp_rng, 
		xlab=aname, ylab=bname, 
		cex=1, cex.axis=1.5, cex.lab=1.5, las=1, type="p", 
		pch=20, col="red");

	arr_width=strwidth(".", units="user")*.15;

	# Draw line from orig pos to adjusted text
	for(i in 1:length(a)){
		points(
			c(a[i], adj_lab_pos[["x"]][i]+arr_width),
			c(b[i], adj_lab_pos[["y"]][i]),
			type="l", col="grey25");
		points(
			c(a[i], adj_lab_pos[["x"]][i]-arr_width),
			c(b[i], adj_lab_pos[["y"]][i]),
			type="l", col="grey25");
		points(
			c(a[i], adj_lab_pos[["x"]][i]),
			c(b[i], adj_lab_pos[["y"]][i]),
			type="l", col="pink4", lwd=.8);
	}

	# Draw adjusted text
	text(adj_lab_pos[["x"]],adj_lab_pos[["y"]], short_names[["a"]]);
	#text(a,b, short_names[["b"]]);

	abline(a=0, b=1, col="blue", lty=2);
	mtext(title, side=3, line=0, outer=T, font=2, cex=2);

}

###############################################################################

alr_cat_names=colnames(A_alr_values);
for(cat_ix in 1:num_used_alr_cat){

	A_alr_val=A_alr_values[,cat_ix,drop=F];
	B_alr_val=B_alr_values[,cat_ix,drop=F];

	cur_cat_name=alr_cat_names[cat_ix];
	
	cat("Fitting: ", cat_ix, ".) ", cur_cat_name, "\n");

	alr_dif=(B_alr_val-A_alr_val);
	names(alr_dif)=subject_ids;

	model_pred=cbind(factors, alr_dif);

	if(length(model_var_arr)==0){
		model_str=paste("alr_dif ~ 1");
	}else{
		model_str=paste("alr_dif ~ ", paste(model_var_arr, collapse=" + "), sep="");
	}

	lm_fit=lm(as.formula(model_str), data=model_pred);
	sum_fit=summary(lm_fit);

	print(sum_fit);
	ab_wilcox_res=wilcox.test(A_alr_val, B_alr_val, paired=T);
	wilcoxon_pval_mat[cur_cat_name, 1]=ab_wilcox_res$p.value;


	plot_ab_comparisons(A_alr_val, B_alr_val, A_subtrahend, B_minuend, 
		pval=ab_wilcox_res$p.value,
		title=paste(cat_ix, ".) ", cur_cat_name, sep=""));


	par(mfrow=c(1,1));

	# ANOVA on full model
	anova_res=anova(lm_fit);

	plot_text(c(
		paste(cat_ix, ".) ", cur_cat_name, ":", sep=""),
		"",
		capture.output(print(sum_fit))
		)
	);

	plot_text(c(
		paste(cat_ix, ".) ", cur_cat_name, ":", sep=""),
                "",
		capture.output(print(anova_res))
		)
	);

	mmps(lm_fit);

	#model_coef_names=setdiff(rownames(sum_fit$coefficients), "(Intercept)");
	# The intercept will be zero if there is no difference, so we should report this.
	model_coef_names=rownames(sum_fit$coefficients);

	print(model_coef_names);

	# Save the covariate result stats
	cat_names=intersect(model_coef_names, cov_coeff_names);
	covariates_coef_mat[cur_cat_name, cat_names]=sum_fit$coefficients[cat_names,"Estimate"];
	covariates_pval_mat[cur_cat_name, cat_names]=sum_fit$coefficients[cat_names,"Pr(>|t|)"];

	rsqrd_mat[cur_cat_name, "R^2"]=sum_fit$r.squared;
	rsqrd_mat[cur_cat_name, "Adj-R^2"]=sum_fit$adj.r.squared;

	if(!is.null(sum_fit$fstatistic)){
		model_pval_mat[cur_cat_name, "F-statistic P-value"]=
			1-pf(sum_fit$fstatistic[1], sum_fit$fstatistic[2], sum_fit$fstatistic[3]);
	}else{
		model_pval_mat[cur_cat_name, "F-statistic P-value"]=NA;
	}

	cat("\n\n*************************************************\n\n");

}

###############################################################################

all.nas=apply(covariates_coef_mat, 2, function(x){all(is.na(x))});
covariates_coef_mat=covariates_coef_mat[,!all.nas,drop=F];

#all.nas=apply(covariates_pval_mat, 2, function(x){all(is.na(x))});
# Keep NAs for placeholder if coefficients were calculable
covariates_pval_mat=covariates_pval_mat[,!all.nas,drop=F];

print(covariates_coef_mat);
print(covariates_pval_mat);
print(rsqrd_mat);
print(model_pval_mat);

par(oma=c(2,1,5,2));

rename_intercept=function(mat, old_n, new_n){
	cname=colnames(mat);
	old_ix=which(cname==old_n);
	cname[old_ix]=new_n;
	colnames(mat)=cname;
	return(mat);
}

covariates_coef_mat=rename_intercept(covariates_coef_mat, "(Intercept)", "\"Difference\"");
covariates_pval_mat=rename_intercept(covariates_pval_mat, "(Intercept)", "\"Difference\"");

print(covariates_coef_mat);

paint_matrix(wilcoxon_pval_mat, plot_min=0, plot_max=1,
	title="Wilcoxon Difference in ALR: P-values",
	high_is_hot=F, deci_pts=4, value.cex=.8); 
mtext(text="(No controlling for covariates)", side=3, cex=.8, font=3, line=2, outer=T);

paint_matrix(covariates_coef_mat, 
        title="Regression Model Coefficient Values",
        high_is_hot=T, deci_pts=3, value.cex=.8);

paint_matrix(covariates_pval_mat, plot_min=0, plot_max=1,
        title="Regression Model Coeff P-Values", 
        high_is_hot=F, deci_pts=4, value.cex=.8);
mtext(text="(H0: Coefficients equal zero, H1: Non-zero Coefficients)", side=3, cex=.9, font=3, line=2, outer=T);

signf_coef_mat=mask_matrix(covariates_coef_mat, covariates_pval_mat, .1, 0);
paint_matrix(signf_coef_mat,
        title="Regression Model Significant Coefficients", 
        high_is_hot=T, deci_pts=3, value.cex=.8, label_zeros=F);
mtext(text="(P-Values < 0.10 Shown)", side=3, cex=.9, font=3, line=2, outer=T);

signf_coef_mat=mask_matrix(covariates_coef_mat, covariates_pval_mat, .05, 0);
paint_matrix(signf_coef_mat,
        title="Regression Model Significant Coefficients", 
        high_is_hot=T, deci_pts=3, value.cex=.8, label_zeros=F);
mtext(text="(P-Values < 0.05 Shown)", side=3, cex=.9, font=3, line=2, outer=T);

signf_coef_mat=mask_matrix(covariates_coef_mat, covariates_pval_mat, .01, 0);
paint_matrix(signf_coef_mat,
        title="Regression Model Significant Coefficients", 
        high_is_hot=T, deci_pts=3, value.cex=.8, label_zeros=F);
mtext(text="(P-Values < 0.01 Shown)", side=3, cex=.9, font=3, line=2, outer=T);


# Plot table of significant associations
stab=signf_as_table(covariates_coef_mat, covariates_pval_mat);
options(width=1000);
plot_text(capture.output(print(stab, quote=F)), .8);


paint_matrix(rsqrd_mat, plot_min=0, plot_max=1,
        title="Regression R^2's", 
        high_is_hot=T, deci_pts=3, value.cex=.8);

paint_matrix(model_pval_mat, plot_min=0, plot_max=1,
        title="Regression Model Fit P-values", 
        high_is_hot=F, deci_pts=4, value.cex=.8);
mtext(text="(H0: Predictors have no contribution to model fit)", side=3, cex=.8, font=3, line=2, outer=T);

###############################################################################

write.table(covariates_pval_mat, file=paste(OutputRoot, ".pval.tsv", sep=""), sep="\t", row.names=T, col.names=T);
write.table(covariates_coef_mat, file=paste(OutputRoot, ".coef.tsv", sep=""), sep="\t", row.names=T, col.names=T);

###############################################################################

outfn=paste(OutputRoot, ".alr_diff.anova.tsv", sep="");
if(TagName!=""){
        fh=file(outfn, "w");
        cat(file=fh, TagName, "\n\t");
}
write.table(round(model_pval_mat[,c("F-statistic P-value"), drop=F], 4),
        file=outfn, append=T,
        sep="\t", row.names=T, col.names=T);

###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
