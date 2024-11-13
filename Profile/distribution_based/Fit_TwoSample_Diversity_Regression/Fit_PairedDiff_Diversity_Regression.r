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

	"outputroot", "o", 2, "character",

	"reference_levels", "c", 2, "character",

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
	"	  -F <column name of sample ids in factor file>\n",
	"	  -S <column name of subject ids in factor file>\n",
	"	-M <list of covariate X's names to include in the model from the factor file>\n",
	"	[-q <required list of variables to include after NA removal>]\n",
	"\n",
	"	-p <pairings map, pairing sample IDs from two groups. Must have header/column names>\n",
	"	-B <Sample Group B (minuend: B-A=diff), (column name in pairings file)\n",
	"	-A <Sample Group A (subtrahend: B-A=diff), (column name in pairings file)\n",
	"\n",
	"	[-o <output filename root>]\n",
	"\n",
	"	[-c <reference levels file for Y's in factor file>]\n",
	"\n",
	"	[-t <tag name>]\n",
	"\n",
	"This script will fit the following model for diversity:\n",
	"\n",
	"	1.) (Diversity_B_minuend[i] - Diversity_A_subtrahend[i]) = covariates\n",
	"\n",
	"The intercept of the model represents the difference between B and A.\n",
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

if(length(opt$factor_samp_id_name)){
	FactorSampleIDName=opt$factor_samp_id_name;
}else{
	FactorSampleIDName="";
}

if(length(opt$factor_subj_id_name)){
	FactorSubjectIDName=opt$factor_subj_id_name;
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
cat("           Sample ID Name: ", FactorSampleIDName, "\n", sep="");
cat("           Subject ID Name: ", FactorSubjectIDName, "\n", sep="");
cat(" Model Variables File: ", ModelVarFile, "\n", sep="");
cat("        Pairings File: ", PairingsFile, "\n", sep="");
cat("            A Minuend: ", A_subtrahend, "\n", sep="");
cat("         B Subtrahend: ", B_minuend, "\n", sep="");
cat("          Output File: ", OutputRoot, "\n", sep="");
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
pdf(paste(OutputRoot, ".paird_diff_div.pdf", sep=""), height=11, width=9.5);

input_files=load_and_reconcile_files(
                sumtab=list(fn=SummaryFile, shorten_cat_names_char=NULL,
                        return_top=NULL, specific_cat_fn=NULL),
		factors=list(fn=FactorsFile, sbj_cname=FactorSubjectIDName,
                        smp_cname=FactorSampleIDName),
                pairs=list(fn=PairingsFile, a_cname=A_subtrahend, b_cname=B_minuend),
                covariates=list(fn=ModelVarFile),
                grpvar=list(fn=""),
                reqvar=list(fn=RequiredFile)
        );

counts=input_files[["SummaryTable_counts"]];
normalized=input_files[["SummaryTable_normalized"]];
factors=input_files[["Factors"]];
good_pairs_map=input_files[["PairsMap"]];
model_var_arr=input_files[["Covariates"]];
required_arr=input_files[["RequiredVariables"]];

# Make the primary key for the factor mat the Sample ID for A
if(FactorSampleIDName!=""){
	rownames(factors)=factors[,FactorSampleIDName];
	factors=factors[good_pairs_map[,A_subtrahend],,drop=F];
}

if(FactorSubjectIDName!=""){
	rownames(factors)=factors[,FactorSubjectIDName];
}

write_file_report(input_files[["Report"]]);

##############################################################################
##############################################################################

# Order/Split the data....
# Order the pairings map by response sample IDs
A_sample_ids=good_pairs_map[,A_subtrahend];
B_sample_ids=good_pairs_map[,B_minuend];
AB_sample_ids=c(A_sample_ids, B_sample_ids);

#------------------------------------------------------------------------------

# Export the sample ids that were actually used
fh=file(paste(OutputRoot, ".used_pairs.meta.tsv", sep=""), "w");
cat(file=fh, "SampleID\tType\n");
for(i in 1:nrow(good_pairs_map)){
	cat(file=fh, A_sample_ids[i], "\t", A_subtrahend, "\n", sep="");
	cat(file=fh, B_sample_ids[i], "\t", B_minuend, "\n", sep="");
}
close(fh);

#------------------------------------------------------------------------------

fh=file(paste(OutputRoot, ".used_pairs.paired.tsv", sep=""), "w");
cat(file=fh, c(A_subtrahend, "\t", B_minuend, "\n"), sep="");
for(i in 1:nrow(good_pairs_map)){
	cat(file=fh, A_sample_ids[i], "\t", B_sample_ids[i], "\n");
}
close(fh);

##############################################################################
# Compute diversity indices
cat("Computing diversity indices:\n");

div_names=c("Tail", "Shannon", "Simpson", "Evenness", "SimpsonsRecip");
num_div_idx=length(div_names);

num_shared_sample_ids=length(AB_sample_ids);
div_mat=matrix(0, nrow=num_shared_sample_ids, ncol=num_div_idx);
colnames(div_mat)=div_names;
rownames(div_mat)=rownames(normalized);

cat("Computing diversity indices across samples.\n");
cat("Num Samples Paired Samples: ", num_shared_sample_ids, "\n");
for(i in 1:num_shared_sample_ids){
        curNorm=normalized[i,];
        zeroFreeNorm=curNorm[curNorm>0]
        div_mat[i,"Tail"]=tail_statistic(zeroFreeNorm);
        div_mat[i,"Shannon"]=-sum(zeroFreeNorm*log(zeroFreeNorm));
        div_mat[i,"Simpson"]=1-sum(curNorm^2);
        div_mat[i,"Evenness"]=div_mat[i,"Shannon"]/log(length(zeroFreeNorm));
        div_mat[i,"SimpsonsRecip"]=1/sum(curNorm^2);
}

print(div_mat);

# Set evenness to 0 if there is only 1 category.
evenness=div_mat[,"Evenness"];
div_mat[is.na(evenness),"Evenness"]=0;

A_div_val=div_mat[A_sample_ids,,drop=F];
B_div_val=div_mat[B_sample_ids,,drop=F];

##############################################################################

num_model_pred=length(model_var_arr);

cov_model_string=paste(model_var_arr, collapse="+");
mmat=model.matrix(as.formula(paste("~", cov_model_string, "-1")), data=as.data.frame(factors));
cov_coeff_names=c("(Intercept)", colnames(mmat));
num_cov_coeff_names=length(cov_coeff_names);

cat("Anticipated Coefficient Names:\n");
print(cov_coeff_names);
cat("\n");

wilcoxon_pval_mat=matrix(NA, nrow=num_div_idx, ncol=1,
        dimnames=list(div_names, "Wilcoxon Paired"));

covariates_coef_mat  =matrix(NA, nrow=num_div_idx, ncol=num_cov_coeff_names, 
	dimnames=list(div_names, cov_coeff_names));
covariates_pval_mat  =matrix(NA, nrow=num_div_idx, ncol=num_cov_coeff_names, 
	dimnames=list(div_names, cov_coeff_names));

rsqrd_mat	=matrix(NA, nrow=num_div_idx, ncol=2, 
	dimnames=list(div_names[1:num_div_idx], c("R^2", "Adj-R^2")));

model_pval_mat	=matrix(NA, nrow=num_div_idx, ncol=1, 
	dimnames=list(div_names[1:num_div_idx], c("F-statistic P-value")));

###############################################################################

plot_ab_comparisons=function(a, b, aname, bname, pval, title){
# Plot histogram of differences and paired samples
	layout_mat=matrix(c(
		1,2,
		3,3,
		4,5
	), ncol=2, byrow=T);
	layout(layout_mat);
	par(oma=c(1,1,3,1));
	par(mar=c(4,4,3,1));

	ranges=range(c(a,b));
	ab=c(a,b);
	ab_hrec=hist(ab, plot=F);

	# Plot individual A/B
	ahrec=hist(a, breaks=ab_hrec$breaks, plot=F);
	bhrec=hist(b, breaks=ab_hrec$breaks, plot=F);
	max_count=max(c(ahrec$counts, bhrec$counts));
	ahrec=hist(a, breaks=ab_hrec$breaks, main=aname, ylim=c(0, max_count), xlab="diversity");
	bhrec=hist(b, breaks=ab_hrec$breaks, main=bname, ylim=c(0, max_count), xlab="diversity");

	# Plot difference
	dif=b-a;
	magn=max(abs(dif));
	diff_plot_range=c(-magn*1.2, magn*1.2);
	print(diff_plot_range);
	hist(dif, xlim=diff_plot_range,
		xlab="",
		main=paste("diff = ", bname, "-", aname, "\n", 
			"Wilcoxon Paired P-value:", round(pval, 4), sep=""), breaks=30);
	abline(v=0, col="blue", lty=2);

	# Plot scatter
	plot(a,b, xlim=ranges, ylim=ranges, xlab=aname, ylab=bname, cex=1.2);
	abline(a=0, b=1, col="blue", lty=2);

	# Plot barplot
	amean=mean(a);
	bmean=mean(b);

	bs_mean=function(x){
		if(length(x)<40){
			return(rep(mean(x),2));
		}
		num_bs=160;
		means=numeric(num_bs);
		for(i in 1:num_bs){
			means[i]=mean(sample(x, replace=T));
		}
		return(quantile(means, c(.025, .975)));
	}
	
	a95ci=bs_mean(a);
	b95ci=bs_mean(b);

	val_max=max(c(a, b));
	mids=barplot(c(amean, bmean), ylim=c(0, val_max*1.2), col="white");
	title(ylab="Bootstrapped 95% CI around Mean", line=2, cex.lab=.8, col.lab="blue");
	axis(1, at=mids, labels=c(aname, bname));
	barsep=diff(mids);
	
	draw95ci=function(ci95, x, width){
		if(ci95[1]==ci95[2]){return};
		points(c(x,x), ci95, type="l", col="blue", lwd=.7);
		points(c(x-width/2,x+width/2), rep(ci95[1],2), type="l", col="blue", lwd=.8);
		points(c(x-width/2,x+width/2), rep(ci95[2],2), type="l", col="blue", lwd=.8);
	}

	drawsamples=function(val, x, width){
		points(x+rnorm(length(val),0,width/5), val, cex=.5, col="grey");
	}

	drawsigbar=function(xs, height, pval){
		if(pval<=.1){
			if(pval<=.001){
				sigch="***";
			}else if(pval<=.01){
				sigch="**";
			}else if(pval<=.05){
				sigch="*";
			}else if(pval<=.1){
				sigch="";
			}
	
			points(xs, rep(height,2), type="l", col="black", lwd=.9);
			points(rep(xs[1],2), c(height,height*.97), type="l", col="black", lwd=.9);
			points(rep(xs[2],2), c(height,height*.97), type="l", col="black", lwd=.9);
			text(mean(xs), height, sigch, adj=c(.5,-1.2));
	
		}
	}

	drawsamples(a, mids[1], barsep/2);
	drawsamples(b, mids[2], barsep/2);

	draw95ci(a95ci, mids[1], barsep/2);
	draw95ci(b95ci, mids[2], barsep/2);

	#label mean and CI
	text(mids[1], a95ci[2], sprintf("mean = %5.3f", amean), col="blue", adj=c(.5,-2.4));

	if(length(a)>=40){
		text(mids[1], a95ci[2], sprintf("95%% CI = (%5.3f, %5.3f)", a95ci[1], a95ci[2]), 
			col="blue", adj=c(.5,-.8), cex=.8);
	}

	text(mids[2], b95ci[2], sprintf("mean = %5.3f", bmean), col="blue", adj=c(.5,-2.4));

	if(length(b)>=40){
		text(mids[2], b95ci[2], sprintf("95%% CI = (%5.3f, %4.3f)", b95ci[1], b95ci[2]), 
			col="blue", adj=c(.5,-.8), cex=.8);
	}

	drawsigbar(mids, val_max*1.1, pval);
	
	mtext(title, side=3, line=0, outer=T, font=2, cex=2);
}

###############################################################################

for(div_ix in 1:num_div_idx){

	A_div=A_div_val[,div_ix,drop=F];
	B_div=B_div_val[,div_ix,drop=F];

	cur_div_name=div_names[div_ix];
	
	cat("Fitting: ", div_ix, ".) ", cur_div_name, "\n");

	div_dif=(B_div-A_div);
	rownames(div_dif)=rownames(good_pairs_map);

	if(!all(rownames(factors) == rownames(div_dif))){
		message("Error: Rownames of Factors and div_dif Matrix aren't the same.");
		quit();
	}

	model_pred=cbind(factors, div_dif);

	if(length(model_var_arr)==0){
		model_str=paste("div_dif ~ 1");
	}else{
		model_str=paste("div_dif ~ ", paste(model_var_arr, collapse=" + "), sep="");
	}

	lm_fit=lm(as.formula(model_str), data=model_pred);
	sum_fit=summary(lm_fit);

	print(sum_fit);
	ab_wilcox_res=wilcox.test(A_div, B_div, paired=T);
	wilcoxon_pval_mat[cur_div_name, 1]=ab_wilcox_res$p.value;

	plot_ab_comparisons(A_div, B_div, A_subtrahend, B_minuend, 
		pval=ab_wilcox_res$p.value,
		title=paste(div_ix, ".) ", cur_div_name, sep=""));


	par(mfrow=c(1,1));

	# ANOVA on full model
	anova_res=anova(lm_fit);

	plot_text(c(
		paste(div_ix, ".) ", cur_div_name, ":", sep=""),
		"",
		capture.output(print(sum_fit))
		)
	);

	plot_text(c(
		paste(div_ix, ".) ", cur_div_name, ":", sep=""),
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
	coef_names=intersect(model_coef_names, cov_coeff_names);
	covariates_coef_mat[cur_div_name, coef_names]=sum_fit$coefficients[coef_names,"Estimate"];
	covariates_pval_mat[cur_div_name, coef_names]=sum_fit$coefficients[coef_names,"Pr(>|t|)"];

	rsqrd_mat[cur_div_name, "R^2"]=sum_fit$r.squared;
	rsqrd_mat[cur_div_name, "Adj-R^2"]=sum_fit$adj.r.squared;

	if(is.null(sum_fit$fstatistic)){
		model_pval_mat[cur_div_name, "F-statistic P-value"]=NA;
	}else{
		model_pval_mat[cur_div_name, "F-statistic P-value"]=
			1-pf(sum_fit$fstatistic[1], sum_fit$fstatistic[2], sum_fit$fstatistic[3]);
	}

	cat("\n\n*************************************************\n\n");


}

###############################################################################

all.nas=apply(covariates_coef_mat, 2, function(x){all(is.na(x) | is.nan(x))});
covariates_coef_mat=covariates_coef_mat[,!all.nas,drop=F];

# If covariates can be estimated, go nan for pvalues
#all.nas=apply(covariates_pval_mat, 2, function(x){all(is.na(x) || is.nan(x))});
covariates_pval_mat=covariates_pval_mat[,!all.nas,drop=F];

cat("Covariates Coef:\n");
print(covariates_coef_mat);
cat("Covariates Pval:\n"); 
print(covariates_pval_mat);
cat("Rsqrd:\n");
print(rsqrd_mat)
cat("Model Pval:\n");
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
	high_is_hot=F, deci_pts=4, value.cex=1); 
mtext(text="(No controlling for covariates)", side=3, cex=.8, font=3, line=2, outer=T);

paint_matrix(covariates_coef_mat, 
        title="Regression Model Coefficient Values",
        high_is_hot=T, deci_pts=2, value.cex=1);

paint_matrix(covariates_pval_mat, plot_min=0, plot_max=1,
        title="Regression Model Coeff P-Values", 
        high_is_hot=F, deci_pts=4, value.cex=1);
mtext(text="(H0: Coefficients equal zero, H1: Non-zero Coefficients)", side=3, cex=.9, font=3, line=2, outer=T);

signf_coef_mat=mask_matrix(covariates_coef_mat, covariates_pval_mat, .1, 0);
paint_matrix(signf_coef_mat,
        title="Regression Model Significant Coefficients", 
        high_is_hot=T, deci_pts=2, value.cex=1, label_zeros=F);
mtext(text="(P-Values < 0.10 Shown)", side=3, cex=.9, font=3, line=2, outer=T);

signf_coef_mat=mask_matrix(covariates_coef_mat, covariates_pval_mat, .05, 0);
paint_matrix(signf_coef_mat,
        title="Regression Model Significant Coefficients", 
        high_is_hot=T, deci_pts=2, value.cex=1, label_zeros=F);
mtext(text="(P-Values < 0.05 Shown)", side=3, cex=.9, font=3, line=2, outer=T);

signf_coef_mat=mask_matrix(covariates_coef_mat, covariates_pval_mat, .01, 0);
paint_matrix(signf_coef_mat,
        title="Regression Model Significant Coefficients", 
        high_is_hot=T, deci_pts=2, value.cex=1, label_zeros=F);
mtext(text="(P-Values < 0.01 Shown)", side=3, cex=.9, font=3, line=2, outer=T);


# Plot table of significant associations
stab=signf_as_table(covariates_coef_mat, covariates_pval_mat);
options(width=1000);
plot_text(capture.output(print(stab, quote=F)), .8);


paint_matrix(rsqrd_mat, plot_min=0, plot_max=1,
        title="Regression R^2's", 
        high_is_hot=T, deci_pts=2, value.cex=1);

paint_matrix(model_pval_mat, plot_min=0, plot_max=1,
        title="Regression Model Fit P-values", 
        high_is_hot=F, deci_pts=2, value.cex=1);
mtext(text="(H0: Predictors have no contribution to model fit)", side=3, cex=.8, font=3, line=2, outer=T);

###############################################################################

outfn=paste(OutputRoot, ".paired_div_diff.model.anova.tsv", sep="");
if(TagName!=""){
        fh=file(outfn, "w");
        cat(file=fh, TagName, "\n\t");
}
write.table(round(model_pval_mat, 4),
	file=outfn, append=T,
	sep="\t", row.names=T, col.names=T);

###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
