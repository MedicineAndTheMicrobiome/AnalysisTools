#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);

params=c(
	"pred_file", "x", 1, "character",
	"resp_file", "y", 1, "character",
	"output_file", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];
script_path=paste(head(strsplit(script_name, "/")[[1]], -1), collapse="/");

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-x <as_pred p-value matrix input file>\n",
	"	-y <as_resp p-value matrix input file>\n",
	"	-o <output file name root>\n",
	"\n",
	"This script will read in two matrices:\n",
	" 1.) Predictor p-values\n",
	" 2.) Response p-values\n",
	"\n",
	"The matrices should have the format:\n",
	"\n",
	"Rows: categories (taxa)\n",
	"Cols: factors (variables)\n",
	"\n",
	"The numbers factors and categories should be the same\n",
	"(as well as the covariates) so that the values can\n",
	"be comparable.\n",
	"\n", sep="");

if(!length(opt$pred_file) || !length(opt$resp_file) || !length(opt$output_file)){
	cat(usage);
	q(status=-1);
}


PredFile=opt$pred_file;
RespFile=opt$resp_file;
OutputRoot=opt$output_file;

input_summary=capture.output({
	cat("\n");
	cat("Input As Predictor File: ", PredFile, "\n");
	cat("Input As Response File: ", RespFile, "\n");
	cat("Output File Root: ", OutputRoot, "\n", sep="");
	cat("\n");
});
cat(input_summary, sep="\n");

##############################################################################

plot_text=function(strings){
	par(family="Courier");
	par(oma=rep(.5,4));
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


##############################################################################
# Open PDF output

pdf(paste(OutputRoot,".pred_vs_resp.pdf", sep=""), height=8.5, width=8.5);

plot_text(input_summary);

#############################################################################	

as_pred=read.table(PredFile);
as_resp=read.table(RespFile);

num_pred_fact=ncol(as_pred);
num_pred_cat=nrow(as_pred);
pred_fact_names=colnames(as_pred);
pred_cat_names=rownames(as_pred);

num_resp_fact=ncol(as_resp);
num_resp_cat=nrow(as_resp);
resp_fact_names=colnames(as_resp);
resp_cat_names=rownames(as_resp);

shrd_fact_names=intersect(pred_fact_names, resp_fact_names);
shrd_cat_names=intersect(pred_cat_names, resp_cat_names);

var_summary=capture.output({
	cat("As Predictor:\n");
	cat("Categories:\n");
	print(pred_cat_names);
	cat("Factors:\n");
	print(pred_fact_names);
	cat("\n");

	cat("As Response:\n");
	cat("Categories:\n");
	print(resp_cat_names);
	cat("Factors:\n");
	print(resp_fact_names);
	cat("\n");
});

shrd_var_summary=capture.output({
	cat("Shared:\n");
	cat("Categories:\n");
	print(shrd_cat_names);
	cat("Factors:\n");
	print(shrd_fact_names);
});


cat(var_summary, sep="\n");
cat(shrd_var_summary, sep="\n");
plot_text(var_summary);
plot_text(shrd_var_summary);

#############################################################################	

plot_resp_pred_scatter=function(factor_name, categories, as_pred_mat, as_resp_mat){


	par(family="");
	x_val=-log10(as_pred_mat[categories, factor_name]);
	y_val=-log10(as_resp_mat[categories, factor_name]);

	# Specify reference significance plots
	siglines_val=c(.1, .05, .01, .001);
	siglines_pos=-log10(siglines_val);

	xrange=range(c(x_val, siglines_pos));
	yrange=range(c(y_val, siglines_pos));

	plot(x_val, y_val, xlab="As Predictor", ylab="As Response",
		main=factor_name,
		xlim=c(xrange[1]-.25, xrange[2]+1), 
		ylim=c(yrange[1]-.25, yrange[2]+1),
		type="n",
		cex.main=2
	);

	# Above untransformed p-values
	axis(side=3, at=siglines_pos, labels=siglines_val, cex.axis=.7);
	axis(side=4, at=siglines_pos, labels=siglines_val, cex.axis=.7);

	abline(a=0, b=1, lwd=2, col="#8888FF");
	abline(h=siglines_pos, col="grey");
	abline(v=siglines_pos, col="grey");

	# Size based on log pvalues
	labelsize=sqrt(x_val^2 + y_val^2)/2;

	# plots points
	points(x_val, y_val, col="red", pch=16, cex=labelsize);

	# Grow label size if it is more significant
	num_cat=length(categories);
	for(i in 1:num_cat){
		labelsize[i]=max(labelsize[i], .4);
		labelsize[i]=min(labelsize[i], 2);
	}

	# Label categories/taxa
	text(x_val, y_val, categories, pos=3, cex=labelsize);

}

#############################################################################

sgn_chr=function(x){

	num_val=length(x);
	result=character(num_val);

	for(i in 1:num_val){
		p=x[i];

		if(p<=.000001){	
			char="*****";
		}else if(p<=.00001){
			char="****";
		}else if(p<=.0001){
			char="*** ";
		}else if(p<=.001){
			char="**  ";
		}else if(p<=.01){
			char="*   ";
		}else if(p<=.05){
			char=":   ";
		}else if(p<=.10){
			char=".   ";
		}else{
			char="    ";
		}

		result[i]=char;
	}
	return(result);
}

#############################################################################

summarize_comparisons=function(factor_name, categories, as_pred_mat, as_resp_mat, cutoff){
	
	# Extract factors
	kept_pred_arr=as_pred_mat[categories, factor_name, drop=F];
	kept_resp_arr=as_resp_mat[categories, factor_name, drop=F];

	#print(kept_pred_arr);
	#print(kept_resp_arr);
	
	# Extract categories above cutoff
	above_cutoff_ix=(kept_pred_arr<=cutoff) | (kept_resp_arr<=cutoff);
	pred_abv=kept_pred_arr[above_cutoff_ix,];
	resp_abv=kept_resp_arr[above_cutoff_ix,];
	cat_abv=categories[above_cutoff_ix];
	#print(pred_abv);
	#print(resp_abv);

	comb_score=exp(-sqrt(log(pred_abv)^2+log(resp_abv)^2));
	#print(comb_score);
	num_kept=sum(above_cutoff_ix);

	# Sort by increasing score
	sort_ix=order(comb_score, decreasing=F);
	pred_abv=pred_abv[sort_ix];
	resp_abv=resp_abv[sort_ix];
	comb_score=comb_score[sort_ix];
	cat_abv=cat_abv[sort_ix];
	log_pval_ratio=log2(pred_abv/resp_abv);

	# Place into value matrix
	val_mat=matrix(numeric(), nrow=num_kept, ncol=4);
	rownames(val_mat)=cat_abv;
	colnames(val_mat)=c(
		"As_Response",
		"As_Predictor",
		"Combined_Score",
		"LogPvalRat"
		);

	val_mat[,"As_Response"]=resp_abv;
	val_mat[,"As_Predictor"]=pred_abv;
	val_mat[,"Combined_Score"]=comb_score;
	val_mat[,"LogPvalRat"]=log_pval_ratio;


	# Place into a formatted text matrix
	text_mat=matrix(character(), nrow=num_kept, ncol=7);
	rownames(text_mat)=cat_abv;
	colnames(text_mat)=c(
		"As_Response", "rsgn", 
		"As_Predictor", "psgn",
		"Combined_Score", "csgn",
		"LogPvalRat"
		);

	text_mat[,"As_Response"]=sprintf("%6.4f", resp_abv);
	text_mat[,"rsgn"]=sgn_chr(resp_abv);
	text_mat[,"As_Predictor"]=sprintf("%6.4f", pred_abv);
	text_mat[,"psgn"]=sgn_chr(pred_abv);
	text_mat[,"Combined_Score"]=sprintf("%6.4f", comb_score);
	text_mat[,"csgn"]=sgn_chr(comb_score);
	text_mat[,"LogPvalRat"]=sprintf("%+2.4f", log_pval_ratio);
	

	# Prepare return values
	result=list();
	result[["values"]]=val_mat;
	result[["formatted"]]=text_mat;
	
	if(num_kept==0){
		return("No categories of significance.");
	}else{
		return(result);
	}

}

#############################################################################	

plot_ratio_comparisons=function(comparison_mat, title=""){

	#print(comparison_mat);

	x=comparison_mat[,"LogPvalRat"];
	y=-log10(comparison_mat[,"Combined_Score"]);
	cat_names=rownames(comparison_mat);

	x_mag=max(abs(range(x)));
	x_range=c(-x_mag*1.2, x_mag*1.2);

	y_max=max(y)*1.1;
	y_range=c(0, y_max);

	num_pred=sum(x<0);
	num_resp=sum(x>0);


	plot(x, y,
		xlim=x_range, ylim=y_range,
		xlab="Log P-value Ratio", ylab="Combined P-value Score",
		type="n",
		main=title,
		cex.main=2
	);

	abline(v=0, col="#BBBBFF", lwd=5);
	abline(v=0, col="#8888FF", lwd=2);

	sig_cutoffs=c(.1, .05, .1, .01, .001);
	log_sig_cut=-log10(sig_cutoffs);

	# Top
	axis(side=3, at=c(-(x_mag*1.2)/2, (x_mag*1.2)/2), 
		labels=c(num_pred, num_resp), 
		cex.axis=1.4, font=2,
		tick=F, line=0);

	pred_txt=ifelse(num_pred==1, "Predictor", "Predictors");
	resp_txt=ifelse(num_resp==1, "Responder", "Responders");
	axis(side=3, at=c(-(x_mag*1.2)/2, (x_mag*1.2)/2), 
		labels=c(pred_txt, resp_txt), 
		cex.axis=1.4, font=2,
		tick=F, line=-1);

	# Right
	abline(h=log_sig_cut, col="black", lwd=1, lty="dotted");
	axis(side=4, at=log_sig_cut, labels=sig_cutoffs, cex.axis=.7);

	# 
	points(x,y, col="red", pch=16);
	text(x, y, label=cat_names, pos=3);

}

#############################################################################	

options(width=200);

for(cur_fact in shrd_fact_names){

	# Plot resp vs pred
	par(mar=c(5,5,6,3));
	plot_resp_pred_scatter(cur_fact, shrd_cat_names, as_pred, as_resp);
	
	# Summary comparisons
	result=summarize_comparisons(cur_fact, shrd_cat_names, as_pred, as_resp, .1);

	# Plot summary (transformed stats)
	par(mar=c(5,5,6,3));
	plot_ratio_comparisons(result[["values"]], title=cur_fact);

	values_mat=result[["values"]];
	formatted_mat=result[["formatted"]];

	# Plot tables in PDF
	# Ordered by combined score
	sorted_ix=order(values_mat[,"Combined_Score"]);
	summary_txt_byScore=capture.output(noquote(formatted_mat[sorted_ix,,drop=F]));

	# Ordered by Log Pvalue Ratio
	sorted_ix=order(values_mat[,"LogPvalRat"]);
	summary_txt_byRatio=capture.output(noquote(formatted_mat[sorted_ix,,drop=F]));
	
	plot_text(c(
		"Sorted By Combined Significance Score:",
		"",
		summary_txt_byScore,
		"",
		"",
		"Sorted By Log P-Value Ratio:",
		"",
		summary_txt_byRatio));
}


#############################################################################	
# Close PDF output

dev.off();

#############################################################################

cat("Done.\n");
print(warnings());
q(status=0);
