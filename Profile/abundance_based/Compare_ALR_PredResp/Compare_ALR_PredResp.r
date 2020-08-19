#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);

params=c(
	"pred_pval_file", "x", 1, "character",
	"resp_pval_file", "y", 1, "character",
	"pred_coef_file", "u", 1, "character",
	"resp_coef_file", "v", 1, "character",
	"output_file", "o", 1, "character",
	"sig_cutoff", "p", 2, "numeric",
	"ratio_cutoff", "r", 2, "numeric",
	"highlight_diag", "d", 2, "logical",
	"as_pred_A_name", "a", 2, "character",
	"as_resp_B_name", "b", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];
script_path=paste(head(strsplit(script_name, "/")[[1]], -1), collapse="/");


DEF_SIGNF_CUTOFF=0.01;
DEF_RATIO_CUTOFF=1; # Log2(2)=1

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-x <as_pred p-value matrix input file (rows=alr/pred, cols=var/resp)>\n",
	"	-y <as_resp p-value matrix input file (rows=alr/resp, cols=var/pred)>\n",
	"	-u <as_pred coefficients matrix input file>\n",
	"	-v <as_resp coefficients matrix input file>\n",
	"	-o <output file name root>\n",
	"	-p <significance cutoff, default=", DEF_SIGNF_CUTOFF, ">\n",
	"	[-r <ratio cutoff for reporting P/R, default=", DEF_RATIO_CUTOFF, ">]\n",
	"	[-d (highlight diagonal for symmetric P/R analyses)]\n",
	"	[-a <'as_pred', sample type>]\n", 
	"	[-b <'as_resp', sample type>]\n", 
	"\n",
	"This script will read in two pairs of matrices:\n",
	" 1.) Predictor p-values and coefficients\n",
	" 2.) Response p-values and coefficients\n",
	"\n",
	"The matrices should have the format:\n",
	"\n",
	"Rows: categories (taxa)\n",
	"Cols: factors (variables)\n",
	"\n",
	"The numbers factors and categories should be the same\n",
	"(as well as the covariates) so that the values can\n",
	"be comparable.\n",
	"\n",
	"If the predictor and response are two samples with the same taxa names, then\n",
	"make sure you specify the -a and -b sample types\n",
	"\n", sep="");

if(!length(opt$pred_pval_file) || !length(opt$resp_pval_file) ||
   !length(opt$pred_coef_file) || !length(opt$resp_coef_file) || !length(opt$output_file)){
	cat(usage);
	q(status=-1);
}


PredPvalFile=opt$pred_pval_file;
RespPvalFile=opt$resp_pval_file;
PredCoefFile=opt$pred_coef_file;
RespCoefFile=opt$resp_coef_file;
OutputRoot=opt$output_file;

if(length(opt$sig_cutoff)){
	SignifCutoff=opt$sig_cutoff;
}else{
	SignifCutoff=DEF_SIGNF_CUTOFF;
}

if(length(opt$ratio_cutoff)){
        RatioCutoff=opt$ratio_cutoff;
}else{
        RatioCutoff=DEF_RATIO_CUTOFF;
}

if(length(opt$highlight_diag)){
        HighlightDiag=T;
}else{
	HighlightDiag=F;
}

AsPredAName="[ALR Categories]";
AsRespBName="[Factors]";
if(length(opt$as_pred_A_name)){
	AsPredAName=opt$as_pred_A_name;
}
if(length(opt$as_resp_B_name)){
	AsRespBName=opt$as_resp_B_name;
}

signf_ext=sprintf("p%02g", SignifCutoff*100);

OutputRoot=paste(OutputRoot, ".pred_vs_resp.", signf_ext, sep="");

input_summary=capture.output({
	cat("\n");
	cat("Pvals:\n");
	cat("Input As Predictor File: ", PredPvalFile, "\n");
	cat("Input As Response File: ", RespPvalFile, "\n");
	cat("\n");
	cat("Coefs:\n");
	cat("Input As Predictor File: ", PredCoefFile, "\n");
	cat("Input As Response File: ", RespCoefFile, "\n");
	cat("\n");
	cat("Output File Root: ", OutputRoot, "\n", sep="");
	cat("\n");
	cat("Significant Cutoff: ", SignifCutoff, "\n", sep="");
	cat("Highlight Diagonal: ", HighlightDiag, "\n", sep="");
	cat("\n");
	cat("AsPred A Name: ", AsPredAName, "\n", sep="");
	cat("AsResp B Name: ", AsRespBName, "\n", sep="");
});
cat(input_summary, sep="\n");

##############################################################################

plot_text=function(strings){

	orig_par=par(no.readonly=T);
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

	par(orig_par);
}


##############################################################################
# Open PDF output

pdf(paste(OutputRoot,".pdf", sep=""), height=7, width=14);

plot_text(input_summary);

#############################################################################	

as_pred_pval=read.table(PredPvalFile);
as_resp_pval=read.table(RespPvalFile);
as_pred_coef=read.table(PredCoefFile);
as_resp_coef=read.table(RespCoefFile);

if(all(is.na(as_pred_pval))){
	plot_text(c(
		"All 'As Predictor' coefficient p-values were NA.",
		"Perhaps sample size was too small for number of covariates + predictors fit."
	));	
	quit(status=0);
}

if(all(is.na(as_resp_pval))){
	plot_text(c(
		"All 'As Response' coefficients p-values were NA.",
		"Perhaps sample size was too small for number of covariates fit."
	));	
	quit(status=0);
}


num_pred_fact=ncol(as_pred_pval);
num_pred_cat=nrow(as_pred_pval);
pred_fact_names=colnames(as_pred_pval);
pred_cat_names=rownames(as_pred_pval);

num_resp_fact=ncol(as_resp_pval);
num_resp_cat=nrow(as_resp_pval);
resp_fact_names=colnames(as_resp_pval);
resp_cat_names=rownames(as_resp_pval);

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
	cat("\n");
	cat("Categories:\n");
	print(shrd_cat_names);
	cat("\n");
	cat("Factors:\n");
	print(shrd_fact_names);
});


cat(var_summary, sep="\n");
cat(shrd_var_summary, sep="\n");
plot_text(var_summary);
plot_text(shrd_var_summary);

#############################################################################	

plot_resp_pred_scatter=function(factor_name, categories, as_pred_mat, as_resp_mat, a_name, b_name){


	par(family="");
	x_val=-log10(as_pred_mat[categories, factor_name]);
	y_val=-log10(as_resp_mat[categories, factor_name]);

	# Specify reference significance plots
	siglines_val=c(.1, .05, .01, .001);
	siglines_pos=-log10(siglines_val);

	xrange=range(c(x_val[!is.na(x_val)], siglines_pos));
	yrange=range(c(y_val[!is.na(y_val)], siglines_pos));

	plot(x_val, y_val, xlab="As Predictor", ylab="As Response",
		main=factor_name,
		xlim=c(xrange[1]-.25, xrange[2]+1), 
		ylim=c(yrange[1]-.25, yrange[2]+1),
		type="n",
		cex.main=2
	);

        mtext(a_name, side=1, line=3.7, cex=1, font=2, col="#FF000088");
        mtext(a_name, side=2, line=3.7, cex=1, font=2, col="#FF000088");
        mtext(b_name, side=3, line=3.7, cex=1, font=2, col="#0000FF88");


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

summarize_comparisons=function(factor_name, categories, 
	as_pred_pval_mat, as_resp_pval_mat, 
	as_pred_coef_mat, as_resp_coef_mat, 
	cutoff){
	
	# Extract factors
	kept_pred_pval_arr=as_pred_pval_mat[categories, factor_name, drop=F];
	kept_resp_pval_arr=as_resp_pval_mat[categories, factor_name, drop=F];
	kept_pred_coef_arr=as_pred_coef_mat[categories, factor_name, drop=F];
	kept_resp_coef_arr=as_resp_coef_mat[categories, factor_name, drop=F];

	print(kept_pred_pval_arr);
	print(kept_resp_pval_arr);
	print(kept_pred_coef_arr);
	print(kept_resp_coef_arr);
	
	# Extract categories above cutoff
	above_cutoff_ix=!is.na(kept_pred_pval_arr) & ((kept_pred_pval_arr<=cutoff) | (kept_resp_pval_arr<=cutoff));
	num_kept=sum(above_cutoff_ix);
	if(is.na(num_kept)){
		num_kept=0;
	}

	# Allocate matrices
	val_mat=matrix(numeric(), nrow=num_kept, ncol=6);
	colnames(val_mat)=c(
		"As_Response",
		"As_Predictor",
		"Combined_Score",
		"LogPvalRat",
		"As_Resp_Dir",
		"As_Pred_Dir"
		);

	text_mat=matrix(character(), nrow=num_kept, ncol=11);
	colnames(text_mat)=c(
		"As_Response", "rsgn", 
		"As_Predictor", "psgn",
		"Combined_Score", "csgn",
		"LogPvalRat",
		"As_Resp_Dir", "rdir",
		"As_Pred_Dir", "pdir"
		);

	if(num_kept>0){

		pred_pval_abv=kept_pred_pval_arr[above_cutoff_ix,];
		resp_pval_abv=kept_resp_pval_arr[above_cutoff_ix,];
		pred_coef_abv=kept_pred_coef_arr[above_cutoff_ix,];
		resp_coef_abv=kept_resp_coef_arr[above_cutoff_ix,];
		cat_abv=categories[above_cutoff_ix];
		#print(pred_abv);
		#print(resp_abv);

		comb_score=exp(-sqrt(log(pred_pval_abv)^2+log(resp_pval_abv)^2));
		#print(comb_score);

		# Sort by increasing score
		sort_ix=order(comb_score, decreasing=F);
		pred_pval_abv=pred_pval_abv[sort_ix];
		resp_pval_abv=resp_pval_abv[sort_ix];
		comb_score=comb_score[sort_ix];
		cat_abv=cat_abv[sort_ix];
		log_pval_ratio=log2(pred_pval_abv/resp_pval_abv);
		pred_coef_abv=pred_coef_abv[sort_ix];
		resp_coef_abv=resp_coef_abv[sort_ix];

		pred_dir_score=-log10(pred_pval_abv)*ifelse(pred_coef_abv>0, 1, -1);
		resp_dir_score=-log10(resp_pval_abv)*ifelse(resp_coef_abv>0, 1, -1);

		# Place into value matrix
		rownames(val_mat)=cat_abv;

		val_mat[,"As_Response"]=resp_pval_abv;
		val_mat[,"As_Predictor"]=pred_pval_abv;
		val_mat[,"Combined_Score"]=comb_score;
		val_mat[,"LogPvalRat"]=log_pval_ratio;
		val_mat[,"As_Resp_Dir"]=resp_dir_score;
		val_mat[,"As_Pred_Dir"]=pred_dir_score;
		
		# Place into a formatted text matrix
		rownames(text_mat)=cat_abv;

		text_mat[,"As_Response"]=sprintf("%6.4f", resp_pval_abv);
		text_mat[,"rsgn"]=sgn_chr(resp_pval_abv);
		text_mat[,"As_Predictor"]=sprintf("%6.4f", pred_pval_abv);
		text_mat[,"psgn"]=sgn_chr(pred_pval_abv);
		text_mat[,"Combined_Score"]=sprintf("%6.4f", comb_score);
		text_mat[,"csgn"]=sgn_chr(comb_score);
		text_mat[,"LogPvalRat"]=sprintf("%+2.4f", log_pval_ratio);
		text_mat[,"As_Resp_Dir"]=sprintf("%+2.2f", resp_dir_score);
		text_mat[,"rdir"]=ifelse(resp_dir_score>0, "+", "-");
		text_mat[,"As_Pred_Dir"]=sprintf("%+2.2f", pred_dir_score);
		text_mat[,"pdir"]=ifelse(pred_dir_score>0, "+", "-");

	}

	# Prepare return values
	result=list();
	result[["values"]]=val_mat;
	result[["formatted"]]=text_mat;
	
	return(result);

}

#############################################################################	

plot_ratio_comparisons=function(comparison_mat, title=""){

	#print(comparison_mat);
	if(nrow(comparison_mat)==0){

		plot(0, 0,
			xlim=c(-2,2), ylim=c(0,1),
			xlab="Log P-value Ratio", ylab="Combined P-value Score",
			type="n",
			main=title,
			cex.main=2
		);
		text(0,.5, "No categories above cutoff");
		
		return();
	}


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

	# Plot points
	chw=par()$cxy[1];

	col=ifelse(comparison_mat[,"As_Pred_Dir"]>0, "green", "red");
	size=sapply(abs(comparison_mat[,"As_Pred_Dir"]), function(x){ min(x, 2.5)});
	shape=ifelse((comparison_mat[,"As_Pred_Dir"])>0, 24, 25);
	points(x-chw/2,y, col="black", bg=col, pch=shape, cex=size);

	col=ifelse(comparison_mat[,"As_Resp_Dir"]>0, "green", "red");
	size=sapply(abs(comparison_mat[,"As_Resp_Dir"]), function(x){ min(x, 2.5)});
	shape=ifelse((comparison_mat[,"As_Resp_Dir"])>0, 24, 25);
	points(x+chw/2,y, col="black", bg=col, pch=shape, cex=size);

	text(x, y, label=cat_names, pos=3);

}

#############################################################################

plot_combined_ratio_comparisons=function(
	comparison_list,
	title="",
	abbrev=0,
	draw_symbol_centroid=T,
	label_factor_centroid=T,
	label_categories=T,
	draw_lines=T
){

	factor_names=names(comparison_list);
	num_factors=length(factor_names);


	if(num_factors==0){
		plot(0,0, xlab="", ylab="", main=title, xaxt="n", yaxt="n");
		text(0,0, "No factors to plot.");	
		return();
	}

	# Pre transform the combined scores
	for(cur_fact in factor_names){
		comparison_list[[cur_fact]][,"Combined_Score"]=
			-log10(comparison_list[[cur_fact]][,"Combined_Score"]);
	}
	

	# Store centroids
	means=matrix(0, nrow=num_factors, ncol=2);
	rownames(means)=factor_names;
	colnames(means)=c("Combined_Score", "LogPvalRat");

	x_mag=0;
	y_max=0;
	y_min=1;

	# Calculate means and identify ranges
	for(cur_fact in factor_names){
		means[cur_fact, "Combined_Score"]=mean(comparison_list[[cur_fact]][,"Combined_Score"]);
		means[cur_fact, "LogPvalRat"]=mean(comparison_list[[cur_fact]][,"LogPvalRat"]);

		x_mag=max(x_mag, abs(comparison_list[[cur_fact]][,"LogPvalRat"]));
		y_max=max(y_max, comparison_list[[cur_fact]][,"Combined_Score"]);
		y_min=min(y_min, comparison_list[[cur_fact]][,"Combined_Score"]);

	}

	cat("Centroids:\n");
	print(means);

	x_range=c(-(x_mag*1.2), x_mag*1.2);
	y_range=c(y_min*.95, y_max*1.05);
	cat("X Range:\n");
	print(x_range);
	cat("Y Range:\n");
	print(y_range);

	par(family="");
	plot(0,0,
		xlim=x_range, ylim=y_range,
		type="n",
		xlab="Log P-value Ratio", ylab="Combined P-value Score",
		main=title,
		cex.main=2
		);

	abline(v=0, col="#BBBBFF", lwd=5);
	abline(v=0, col="#8888FF", lwd=2);

	sig_cutoffs=c(.1, .05, .1, .01, .001);
	log_sig_cut=-log10(sig_cutoffs);

	axis(side=3, at=c(-(x_mag*1.2)/2, (x_mag*1.2)/2), 
		labels=c("Predictors", "Responders"), 
		cex.axis=1.4, font=2,
		tick=F, line=-1);

	# Right
	abline(h=log_sig_cut, col="black", lwd=1, lty="dotted");
	axis(side=4, at=log_sig_cut, labels=sig_cutoffs, cex.axis=.7);

	#----------------------------------------------------------------------
	# Optional Section:

	if(draw_symbol_centroid){
		for(i in 1:num_factors){
			cur_fact=factor_names[i];
			if(is.nan(means[cur_fact, "LogPvalRat"])){
				next;
			}
			points(
				means[cur_fact, "LogPvalRat"],
				means[cur_fact, "Combined_Score"],
				bg=i, cex=2, col="black", pch=21
			)
		}
	}

	if(draw_lines){
		for(i in 1:num_factors){
			cur_fact=factor_names[i];
			cat_names=rownames(comparison_list[[cur_fact]]);
			num_rows=length(cat_names);
			if(num_rows==0){
				next;
			}
			for(cat_ix in 1:num_rows){
				points(
					c(means[cur_fact, "LogPvalRat"],
					comparison_list[[cur_fact]][cat_ix,"LogPvalRat"]),
					c(means[cur_fact, "Combined_Score"],
					comparison_list[[cur_fact]][cat_ix,"Combined_Score"]),
					col=i,
					type="l"
				);

				points(
					comparison_list[[cur_fact]][cat_ix,"LogPvalRat"],
					comparison_list[[cur_fact]][cat_ix,"Combined_Score"],
					col=i,
					type="p",
					pch=15,
					cex=.25
				);
			}
		}


	}

	
	if(label_categories){
		for(i in 1:num_factors){
			cur_fact=factor_names[i];
			cat_names=rownames(comparison_list[[cur_fact]]);
			num_rows=length(cat_names);
			if(num_rows==0){
				next;
			}
			for(cat_ix in 1:num_rows){

				xpos=comparison_list[[cur_fact]][cat_ix,"LogPvalRat"];
				ypos=comparison_list[[cur_fact]][cat_ix,"Combined_Score"];

				adj_x=.5;
				adj_y=.5;
				label_size=1;
			
				if(draw_lines){
					# Move labels away from lines

					if(xpos>means[cur_fact, "LogPvalRat"]){
						adj_x=-.25
					}else if(xpos<means[cur_fact, "LogPvalRat"]){
						adj_x=1
					}

					if(ypos>means[cur_fact, "Combined_Score"]){
						adj_y=-.5
					}else if(ypos<means[cur_fact, "Combined_Score"]){
						adj_y=1
					}
					
					label_size=.75;
				}
			
				label=cat_names[cat_ix];
				if(abbrev>0){
					label=substr(label,1,abbrev);
				}

				text(
					xpos, ypos,
					label=label,
					cex=label_size,
					col=i,
					adj=c(adj_x, adj_y)
				);
			}
		}
	}
	

	if(label_factor_centroid){
	# Label centroids with factor names
		for(i in 1:num_factors){
			cur_fact=factor_names[i];
			if(is.nan(means[cur_fact, "LogPvalRat"])){
				next;
			}
			legend(
				means[cur_fact, "LogPvalRat"],
				means[cur_fact, "Combined_Score"],
				legend=cur_fact,
				xjust=.5, yjust=.5,
				x.intersp=-.4, y.intersp=.2,
				text.col=i, bg="white"
			);
		}
	}
	
}

#############################################################################	

plot_legend_for_combined_ratio_comparisons=function(comparison_list, abbrev=0){

	factor_names=names(comparison_list);
	num_factors=length(factor_names);

	if(num_factors==0){
		plot(0,0, ylim=c(0,1), xlim=c(0,1),
			xlab="", ylab="", main="",
			type="n", bty="n", xaxt="n", yaxt="n");
		text(.5,.5, "No colors assigned.");
		return();
	}
		
	orig_par=par(no.readonly=T);
	par(mar=c(1,1,.5,1));
	plot(0,0, ylim=c(0,1), xlim=c(0,1),
		xlab="", ylab="", main="",
		type="n", bty="n", xaxt="n", yaxt="n");

	#legend(0,1, bty="n", legend=factor_names, fill=1:num_factors, y.intersp=0);

	splits=max(num_factors, 20);

	divisions=1/splits;
	label_size=divisions*10;	

	for(i in 1:num_factors){
		# Plot glyph
		points(0,1-i*divisions, col="black", bg=i, pch=22, cex=label_size*4); 

		# Plot Grouping name
		text(.03, 1-i*divisions, factor_names[i], cex=label_size*1.25, pos=4);

		# Plot Group members
		members=rownames(comparison_list[[factor_names[i]]]);
		if(abbrev>0){
			members=substr(members, 1, abbrev);
		}

		members_str=paste(members, collapse=", ");
		text(.06, 1-(i*divisions+divisions/2), members_str, cex=label_size*.85, pos=4, font=3);
	}
	par(orig_par);

}

#############################################################################	

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

#############################################################################	

invert_records=function(comparison_list){
	factor_names=names(comparison_list);
	num_factors=length(factor_names);

	categories=c();
	for(i in 1:num_factors){
		cur_factor=factor_names[i];
		cur_mat=comparison_list[[cur_factor]];
		categories=c(categories, rownames(cur_mat));
	}

	categories=unique(categories);

	cat("Unique Categories in Comparison List:\n");
	print(categories);

	inverted_list=list();

	for(i in 1:num_factors){
		cur_factor=factor_names[i];
                cur_mat=comparison_list[[cur_factor]];

		num_rows=nrow(cur_mat);
		if(num_rows>0){
			row_names=rownames(cur_mat);

			for(row_ix in 1:num_rows){
				row_name=row_names[row_ix];
				if(is.null(inverted_list[[row_name]])){
			
					cat("Allocating for: ", row_name, "\n");
					val_mat=matrix(numeric(), nrow=0, ncol=6);
					colnames(val_mat)=c(
						"As_Response",
						"As_Predictor",
						"Combined_Score",
						"LogPvalRat",
						"As_Resp_Dir",
						"As_Pred_Dir"
						);
				
					inverted_list[[row_name]]=val_mat;
				}

				cat("Saving for: ", cur_factor, "\n");

				prev_rownames=rownames(inverted_list[[row_name]]);

				inverted_list[[row_name]]=
					rbind(inverted_list[[row_name]], cur_mat[row_ix,]);

				
				rownames(inverted_list[[row_name]])=c(
					prev_rownames,
					cur_factor
				);

			}
		}
	}		
	return(inverted_list);

}

#############################################################################	

output_text_summary=function(comparison_list, root_filename){
	factor_names=names(comparison_list);
        num_factors=length(factor_names);

	fh=file(paste(root_filename, ".stats.tsv", sep=""), "w");

	outstr=paste(
		"Group", 
		"Factor",
		"Category",
		"As_Resp_Pval",
		"As_Resp_Signf",
		"As_Pred_Pval",
		"As_Pred_Signf",
		"Combined_Score",
		"Combined_Signf",
		"LogPvalRatio",
		sep="\t");

	cat(file=fh, "#", outstr, "\n", sep="");

	for(i in 1:num_factors){

		fact_name=factor_names[i];

		mat=comparison_list[[fact_name]];
		cat_names=rownames(mat);
		num_cat=nrow(mat);

		if(num_cat>0){
			for(cat_ix in 1:num_cat){
				outstr=paste(
					root_filename,
					fact_name,
					cat_names[cat_ix],
					mat[cat_ix, "As_Response"],		
					mat[cat_ix, "rsgn"],		
					mat[cat_ix, "As_Predictor"],		
					mat[cat_ix, "psgn"],		
					mat[cat_ix, "Combined_Score"],		
					mat[cat_ix, "csgn"],		
					mat[cat_ix, "LogPvalRat"],
					sep="\t"
					);

			cat(file=fh, outstr, "\n");
			}
		}else{
			cat("No records for: ", fact_name, "\n", sep="");
		}


	}


}


#############################################################################	

summarize_to_matrix=function(combined_ratio_coef, shrd_fact_names, shrd_cat_names, ratio_thres){

	cat("Plot summary matrices with ratio > ", ratio_thres, ".\n");
	#print(combined_ratio_coef);
	#print(shrd_fact_names);
	#print(shrd_cat_names);

	absthres=abs(ratio_thres);

	shrd_fact_names=sort(shrd_fact_names);
	shrd_cat_names=sort(shrd_cat_names);

	num_shrd_fact_names=length(shrd_fact_names);
	num_shrd_cat_names=length(shrd_cat_names);

	fact_cat_mat_coefdir=matrix(0, nrow=num_shrd_cat_names, ncol=num_shrd_fact_names);
	rownames(fact_cat_mat_coefdir)=shrd_cat_names;
	colnames(fact_cat_mat_coefdir)=shrd_fact_names;

	fact_cat_mat_predresp=matrix("", nrow=num_shrd_cat_names, ncol=num_shrd_fact_names);
	rownames(fact_cat_mat_predresp)=shrd_cat_names;
	colnames(fact_cat_mat_predresp)=shrd_fact_names;


	#print(fact_cat_mat_coefdir);
	#print(fact_cat_mat_predresp);

	#print(combined_ratio_coef);
	signf_fact=names(combined_ratio_coef);

	for(fact in signf_fact){
		cat_table=combined_ratio_coef[[fact]];
		categories=rownames(cat_table);
	
		for(cat in categories){
			directionality=0;

			# Assign to predictor or responder
			if(cat_table[cat, "LogPvalRat"]>absthres){

				# Responder
				fact_cat_mat_predresp[cat, fact]="R";

				if(cat_table[cat, "As_Resp_Dir"]>0){
					directionality=1;
				}else if(cat_table[cat, "As_Resp_Dir"]<0){
					directionality=-1;
				}

			}else if(cat_table[cat, "LogPvalRat"]<(-absthres)){

				# Predictor	
				fact_cat_mat_predresp[cat, fact]="P";

				if(cat_table[cat, "As_Pred_Dir"]>0){
					directionality=1;
				}else if(cat_table[cat, "As_Pred_Dir"]<0){
					directionality=-1;
				}

			}else{

				# If directionality is ambiguous, just report
				# directionality if it is consistent
				if(
					cat_table[cat, "As_Pred_Dir"]>0 &&
					cat_table[cat, "As_Resp_Dir"]>0){
						directionality=1;
				}else if(
					cat_table[cat, "As_Pred_Dir"]<0 &&
					cat_table[cat, "As_Resp_Dir"]<0){
						directionality=-1;
				}else{
					directionality=0;
				}

				fact_cat_mat_predresp[cat, fact]="";
			}

			fact_cat_mat_coefdir[cat, fact]=directionality;
		}
	}

	#print(fact_cat_mat_predresp);
	#print(fact_cat_mat_coefdir);
	res=list();
	res[["pred.resp"]]=fact_cat_mat_predresp;
	res[["coeff.dir"]]=fact_cat_mat_coefdir;
	return(res);
}

#############################################################################

plot_predresp_matrix=function(matrices, highlight_diag=F, ratio_thres, signf_thres, a_name="", b_name=""){
	
	pr_mat=matrices[["pred.resp"]];
	dir_mat=matrices[["coeff.dir"]];

	#print(pr_mat);

	num_facts=ncol(pr_mat);
	num_cat=nrow(pr_mat);

	cat_names=rownames(pr_mat);
	fact_names=colnames(pr_mat);


	# Look for rows/columns with associations
	hl=function(x){
		return(any(x!=""));
	}	

	col_highlight=apply(pr_mat, 2, hl);
	row_highlight=apply(pr_mat, 1, hl);
	col_hl_ix=which(col_highlight);
	row_hl_ix=which(row_highlight);


	# Start plot
	par(mar=c(4,14,14,1));
	off=.25;
	plot(0,0, xlim=c(1, num_facts+1), ylim=c(1, num_cat+1), xaxt="n", yaxt="n",
		xlab="", ylab="", type="n");
	title(xlab=paste("p-values < ", signf_thres, "      |log2(p-value ratio)| > ", ratio_thres, sep=""),
		line=1);

	# Help text at top left
	a_name_placeholder=paste(rep(" ", nchar(a_name)), collapse="");
	b_name_placeholder=paste(rep(" ", nchar(b_name)), collapse="");
	mtext("Interpretations:", 
		3, adj=0, line=11, at=0, cex=1, font=2, col="black", family="mono");

	mtext(paste(a_name_placeholder, " [R]esponds to ", b_name_placeholder, sep=""), 
		3, adj=0, line=10, at=0, cex=1, font=2, col="orange", family="mono");
	mtext(paste(a_name_placeholder, " [P]redicts ", b_name_placeholder, sep=""), 
		3, adj=0, line= 9, at=0, cex=1, font=2, col="blue", family="mono");

	mtext(paste(a_name, "               ", b_name_placeholder, sep=""), 
		3, adj=0, line=10, at=0, cex=1, font=2, col="#FF000088", family="mono");
	mtext(paste(a_name, "            ", b_name_placeholder, sep=""), 
		3, adj=0, line= 9, at=0, cex=1, font=2, col="#FF000088", family="mono");

	mtext(paste(a_name_placeholder, "               ", b_name, sep=""), 
		3, adj=0, line=10, at=0, cex=1, font=2, col="#0000FF88", family="mono");
	mtext(paste(a_name_placeholder, "            ", b_name, sep=""), 
		3, adj=0, line= 9, at=0, cex=1, font=2, col="#0000FF88", family="mono");



	if(highlight_diag){
		abline(a=0, b=1, lwd=12, col="grey80");
	}

	# Comput locations for glyphs
	horiz=(1:num_cat)+off;
	vert=(1:num_facts)+off;

	# Draw guide lines for empty associations
	abline(h=horiz, lwd=2, col="grey97");
	abline(v=vert, lwd=3, col="grey97");
	# Draw guid lines for associations
	abline(h=horiz[row_highlight], lwd=1, col="grey80");
	abline(v=vert[col_highlight], lwd=1, col="grey80");

	glyph_size=1.2;
	for(x in 1:num_facts){
		for(y in 1:num_cat){

			if(pr_mat[y,x]=="P"){
				text(x+off,y+off, "P", col="blue", cex=glyph_size);
			}else if(pr_mat[y,x]=="R"){
				text(x+off,y+off, "R", col="orange", cex=glyph_size);
			}

			if(dir_mat[y,x]==1){
				points(x+1.5*off,y+off, pch=24, col="green", bg="green", cex=glyph_size);
			}else if(dir_mat[y,x]==-1){
				points(x+1.5*off,y+off, pch=25, col="red", bg="red", cex=glyph_size);
			}

		}

	}

	# Draw axes labels with associations darker
	axis(side=2, at=horiz[row_highlight], labels=cat_names[row_highlight], 
		las=2, col.axis="black", col="black");
	axis(side=3, at=vert[col_highlight], labels=fact_names[col_highlight], 
		las=2, col.axis="black", col="black");
	# Draw axes labels without associations lighter
	axis(side=2, at=horiz[!row_highlight], labels=cat_names[!row_highlight], 
		las=2, col.axis="grey80", col="grey80");
	axis(side=3, at=vert[!col_highlight], labels=fact_names[!col_highlight], 
		las=2, col.axis="grey80", col="grey80");

	# Label direction type
	cat("a/b: ", a_name, "/", b_name, "\n");
	mtext(a_name, side=2, line=-1.5, cex=1.75, font=2, col="#FF000044");
	mtext(b_name, side=3, line=-1.5, cex=1.75, font=2, col="#0000FF44");

}

#############################################################################

plot_predresp_line_diagram=function(
		matrices, ratio_thres, signf_thres, a_name="", b_name="", 
		subset="both", title=""
	){
	
	cat("Starting plot_predresp_line_diagram plot:\n");

	dir_mat=matrices[["coeff.dir"]];
	pr_mat=matrices[["pred.resp"]];

	if(subset=="both"){
	}else if(subset=="responders"){
		pr_mat=apply(pr_mat, 1:2, function(x){ ifelse(x=="P", "", x);});
	}else if(subset=="predictors"){
		pr_mat=apply(pr_mat, 1:2, function(x){ ifelse(x=="R", "", x);});
	}

	#print(pr_mat);
	#print(dir_mat);

	nrows_A_left=nrow(pr_mat);
	ncols_B_right=ncol(pr_mat);

	rows_A_left_names=rownames(pr_mat);
	cols_B_right_names=colnames(pr_mat);
	
	# Calc spacing between variables on each side
	spacing_A=seq(0,1, length.out=nrows_A_left+2)[2:(nrows_A_left+1)];
	spacing_B=seq(0,1, length.out=ncols_B_right+2)[2:(ncols_B_right+1)];

	par(mar=c(5,20,3,20));
	line_margin=.05;
	plot(0, type="n", 
		xlim=c(0-line_margin,1+line_margin), 
		ylim=c(0+.05,1-.05), xlab="", ylab="", xaxt="n", yaxt="n");

	title(main=title);

	mtext(a_name, side=2, line=15, cex=4, col="#FF000088");
	mtext(b_name, side=4, line=15, cex=4, col="#0000FF88");

	# Label variables on left/right axis
	for(rowix in 1:nrows_A_left){
		if(all(pr_mat[rowix,]=="")){
			label_col="grey";
		}else{
			label_col="#FF0000";
		}
		axis(side=2, at=spacing_A[rowix], labels=rows_A_left_names[rowix], las=2, col.axis=label_col);
	}
	for(colix in 1:ncols_B_right){
		if(all(pr_mat[,colix]=="")){
			label_col="grey";
		}else{
			label_col="#0000FF";
		}
		axis(side=4, at=spacing_B[colix], labels=cols_B_right_names[colix], las=2, col.axis=label_col);
	}

	# Draw arrows when they are responders
	for(rowix in 1:nrows_A_left){
		if(any(pr_mat[rowix,]=="R")){
			arrows(0, spacing_A[rowix], 0-line_margin, spacing_A[rowix], 
				length=.1, angle=22.5, col="orange", lwd=1.5);
		}
	}
	for(colix in 1:ncols_B_right){
		if(any(pr_mat[,colix]=="P")){
			arrows(1, spacing_B[colix], 1+line_margin, spacing_B[colix], 
				length=.1, angle=22.5, col="blue", lwd=1.5);
		}
	}

	# Function to place up/down green/red arrows base on slop of the line
	label_direction=function(xs, ys, side, dir, outline_col, lift, dist_from_end=.05){
		# y=mx+b;
		slope=(ys[2]-ys[1])/(xs[2]-xs[1]);
		yintc=ys[1]-slope*xs[1];

		xloc=ifelse(side==0, dist_from_end, 1-dist_from_end);
		yloc=slope*xloc+yintc;
	
		redgreen=ifelse(dir==1, "green", "red");
		updown=ifelse(dir==1, 24, 25);
		points(xloc, yloc+lift, pch=updown, col=outline_col, bg=redgreen, cex=1.5);
	}

	# Draw lines before drawing arrows
	for(rowix in 1:nrows_A_left){
		for(colix in 1:ncols_B_right){
			predresp=pr_mat[rows_A_left_names[rowix], cols_B_right_names[colix]];
			if(predresp==""){
				#noop
			}else if(predresp=="P"){
				points(c(0,1), c(spacing_A[rowix], spacing_B[colix]), type="l", lwd=1.5,
					 col="blue");
			}else if(predresp=="R"){
				points(c(0,1), c(spacing_A[rowix], spacing_B[colix]), type="l", lwd=1.5,
					col="orange");
			}else{
				cat("ERROR: bad pred/resp character.\n");
				quit(status=-1);
			}
		}
	}

	# How far above line to draw coefficient arrows
	lift_a=(spacing_A[2]-spacing_A[1])*.1;
	lift_b=(spacing_B[2]-spacing_B[1])*.1;

	# Draw direction arrow over the lines
	for(rowix in 1:nrows_A_left){
		for(colix in 1:ncols_B_right){
			predresp=pr_mat[rows_A_left_names[rowix], cols_B_right_names[colix]];
			if(predresp==""){
				#noop
			}else if(predresp=="P"){
				label_direction(
					c(0,1), c(spacing_A[rowix], spacing_B[colix]), side=0,
					dir_mat[rowix,colix], outline_col="blue", lift=lift_a
					);
			}else if(predresp=="R"){
				label_direction(
					c(0,1), c(spacing_A[rowix], spacing_B[colix]), side=1,
					dir_mat[rowix,colix], outline_col="orange", lift=lift_b
					);
			}else{
				cat("ERROR: bad pred/resp character.\n");
				quit(status=-1);
			}
		}
	}

	# Help text at top 
	a_name_placeholder=paste(rep(" ", nchar(a_name)), collapse="");
	b_name_placeholder=paste(rep(" ", nchar(b_name)), collapse="");
	xpos=0;
	mtext("Interpretations:", 
		1, adj=0, line=1, at=xpos, cex=1, font=2, col="black", family="mono");

	mtext(paste(a_name_placeholder, " <-Responds to ", b_name_placeholder, sep=""), 
		1, adj=0, line=2, at=xpos, cex=1, font=2, col="orange", family="mono");
	mtext(paste(a_name_placeholder, " Predicts-> ", b_name_placeholder, sep=""), 
		1, adj=0, line=3, at=xpos, cex=1, font=2, col="blue", family="mono");

	mtext(paste(a_name, "               ", b_name_placeholder, sep=""), 
		1, adj=0, line=2, at=xpos, cex=1, font=2, col="#FF000088", family="mono");
	mtext(paste(a_name, "            ", b_name_placeholder, sep=""), 
		1, adj=0, line=3, at=xpos, cex=1, font=2, col="#FF000088", family="mono");

	mtext(paste(a_name_placeholder, "               ", b_name, sep=""), 
		1, adj=0, line=2, at=xpos, cex=1, font=2, col="#0000FF88", family="mono");
	mtext(paste(a_name_placeholder, "            ", b_name, sep=""), 
		1, adj=0, line=3, at=xpos, cex=1, font=2, col="#0000FF88", family="mono");
	

	cat("Done with PredResp Line Diagram.\n");
}

#############################################################################	

plot_venn=function(matrices, a_name="", b_name=""){

	pr_mat=matrices[["pred.resp"]];
	dir_mat=matrices[["coeff.dir"]];

	fact_names=colnames(pr_mat);
	cat_names=rownames(pr_mat);

	num_factors=length(fact_names);
	num_categories=length(cat_names);

	vennize=function(x){
		preds=any(x=="P");
		resps=any(x=="R");
		if(preds && resps){
			return("Both");
		}else if(preds){
			return("Predictor");
		}else if(resps){
			return("Responder");
		}else{
			return("Neither");
		}
	}


	plot_venn_table=function(x, title, lab_col, flip){
		
		width=4;

		pred_ix=(x=="Predictor");
		resp_ix=(x=="Responder");
		both_ix=(x=="Both");
		neit_ix=(x=="Neither");

		num_pred=sum(pred_ix);
		num_resp=sum(resp_ix);
		num_both=sum(both_ix);
		num_neit=sum(neit_ix);

		var_names=names(x);

		max_list_length=max(c(num_pred, num_resp, num_both, num_neit, 15));

		plot(0,0, type="n", xaxt="n", yaxt="n", xlab="", ylab="",
			xlim=c(-1,4), ylim=c(0, max_list_length));

		if(flip==F){
			axis(3, at=c(0,1,2,3), labels=c("Predictors", "Both", "Responders", "Neither"),
				font.axis=2, cex.axis=1.5
			);
		}else{
			axis(3, at=c(0,1,2,3), labels=c("Responders", "Both", "Predictors", "Neither"),
				font.axis=2, cex.axis=1.5
			);
		}

		mtext(title, side=3, line=3, cex=1.5, font=2, col=lab_col);

		plot_list=function(varlist, x, y){
			llen=length(varlist);
			if(llen>0){
				for(i in 1:llen){
					text(x, y-i, varlist[i]); 
				}
			}
		}

		plot_list(var_names[pred_ix], 0, max_list_length);
		plot_list(var_names[both_ix], 1, max_list_length);
		plot_list(var_names[resp_ix], 2, max_list_length);
		plot_list(var_names[neit_ix], 3, max_list_length);


	}

	#print(pr_mat);

	cat_pr_venn=apply(pr_mat, 1, vennize);
	fact_pr_venn=apply(pr_mat, 2, vennize);

	
	par(mar=c(1,1,5,1));

	plot_venn_table(cat_pr_venn, a_name, "#FF000044", flip=F);
	plot_venn_table(fact_pr_venn, b_name, "#0000FF44", flip=T);

	#print(cat_pr_venn);
	#print(fact_pr_venn);

	
}

#############################################################################	

options(width=200);

combined_records=list();
combined_records_txt=list();

num_shrd_factors=length(shrd_fact_names);
factor_colors=get_colors(max(2,num_shrd_factors), alpha=1);
palette(factor_colors);

layout_mat=matrix(c(1,2), nrow=1, ncol=2);


for(cur_fact in shrd_fact_names){

	# Plot resp vs pred
	par(mar=c(5,5,6,3));
	layout(layout_mat);
	plot_resp_pred_scatter(cur_fact, shrd_cat_names, as_pred_pval, as_resp_pval, AsPredAName, AsRespBName);
	
	# Summary comparisons
	result=summarize_comparisons(cur_fact, shrd_cat_names, 
		as_pred_pval, as_resp_pval, 
		as_pred_coef, as_resp_coef, 
		SignifCutoff);

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
	
	par(mfrow=c(1,1));
	plot_text(c(
		"Factor:",
		cur_fact,
		"",
		"Sorted By Combined Significance Score:",
		"",
		summary_txt_byScore,
		"",
		"",
		"Sorted By Log P-Value Ratio:",
		"",
		summary_txt_byRatio));

	# Store summarized for combined
	combined_records[[cur_fact]]=values_mat;
	combined_records_txt[[cur_fact]]=formatted_mat;
}


#############################################################################	

print(combined_records);

par(mfrow=c(1,1));
plot_text(capture.output(combined_records));

par(mar=c(5,5,6,3));

layout(layout_mat);
plot_combined_ratio_comparisons(combined_records, title="Combined: Scattered Categories", abbrev=4,
	label_factor_centroid=F, draw_symbol_centroid=F, label_categories=T, draw_lines=F);
plot_legend_for_combined_ratio_comparisons(combined_records, abbrev=6);


layout(layout_mat);
plot_combined_ratio_comparisons(combined_records, title="Combined: Connected Categories", abbrev=4,
	label_factor_centroid=F, draw_symbol_centroid=F, label_categories=T, draw_lines=T);
plot_legend_for_combined_ratio_comparisons(combined_records, abbrev=6);


layout(layout_mat);
plot_combined_ratio_comparisons(combined_records, title="Combined: Scattered Categories", abbrev=4,
	label_factor_centroid=F, draw_symbol_centroid=T, label_categories=T, draw_lines=T);
plot_legend_for_combined_ratio_comparisons(combined_records, abbrev=6);


layout(layout_mat);
plot_combined_ratio_comparisons(combined_records, title="Combined: Factors Labeled", abbrev=4,
	label_factor_centroid=T, draw_symbol_centroid=F, label_categories=T, draw_lines=T);
plot_legend_for_combined_ratio_comparisons(combined_records, abbrev=6);



# Merge records
inverted_records=invert_records(combined_records);
print(inverted_records);
num_categories=length(names(inverted_records));
category_colors=get_colors(max(2,num_categories), alpha=1);
palette(category_colors);

par(mfrow=c(1,1));
plot_text(capture.output(inverted_records));

layout(layout_mat);
plot_combined_ratio_comparisons(inverted_records, title="Combined: Scattered Factors", abbrev=0,
	label_factor_centroid=F, draw_symbol_centroid=F, label_categories=T, draw_lines=F);
plot_legend_for_combined_ratio_comparisons(inverted_records);


layout(layout_mat);
plot_combined_ratio_comparisons(inverted_records, title="Combined: Connected Factors", abbrev=0,
	label_factor_centroid=F, draw_symbol_centroid=F, label_categories=T, draw_lines=T);
plot_legend_for_combined_ratio_comparisons(inverted_records);


layout(layout_mat);
plot_combined_ratio_comparisons(inverted_records, title="Combined: Scattered Factors", abbrev=0,
	label_factor_centroid=F, draw_symbol_centroid=T, label_categories=T, draw_lines=T);
plot_legend_for_combined_ratio_comparisons(inverted_records);


layout(layout_mat);
plot_combined_ratio_comparisons(inverted_records, title="Combined: Categories Labeled", abbrev=0,
	label_factor_centroid=T, draw_symbol_centroid=F, label_categories=T, draw_lines=T);
plot_legend_for_combined_ratio_comparisons(inverted_records);


#############################################################################	
print(combined_records);
par(mfrow=c(1,1));

# Loose cutoffs
pred_resp_mat=summarize_to_matrix(combined_records, shrd_fact_names, shrd_cat_names, 
	ratio_thres=0);
plot_predresp_matrix(pred_resp_mat, highlight_diag=HighlightDiag, 
	ratio_thres=0, signf_thres=SignifCutoff, a_name=AsPredAName, b_name=AsRespBName);
plot_venn(pred_resp_mat, a_name=AsPredAName, b_name=AsRespBName);

plot_predresp_line_diagram(pred_resp_mat,
        ratio_thres=0, signf_thres=SignifCutoff, a_name=AsPredAName, b_name=AsRespBName,
	subset="both", title=paste("Predictors and Responders"));
plot_predresp_line_diagram(pred_resp_mat,
        ratio_thres=0, signf_thres=SignifCutoff, a_name=AsPredAName, b_name=AsRespBName,
	subset="predictors", title=paste("Predictors Only"));
plot_predresp_line_diagram(pred_resp_mat,
        ratio_thres=0, signf_thres=SignifCutoff, a_name=AsPredAName, b_name=AsRespBName,
	subset="responders", title=paste("Responders Only"));


# More stringent/user defined cutoff
pred_resp_mat=summarize_to_matrix(combined_records, shrd_fact_names, shrd_cat_names, 
	ratio_thres=RatioCutoff);
plot_predresp_matrix(pred_resp_mat, highlight_diag=HighlightDiag, 
	ratio_thres=RatioCutoff, signf_thres=SignifCutoff, a_name=AsPredAName, b_name=AsRespBName);
plot_venn(pred_resp_mat, a_name=AsPredAName, b_name=AsRespBName);

plot_predresp_line_diagram(pred_resp_mat,
        ratio_thres=0, signf_thres=SignifCutoff, a_name=AsPredAName, b_name=AsRespBName,
	subset="both", title=paste("Predictors and Responders: (Ratio Cutoff>", RatioCutoff, ")", sep=""));
plot_predresp_line_diagram(pred_resp_mat,
        ratio_thres=0, signf_thres=SignifCutoff, a_name=AsPredAName, b_name=AsRespBName,
	subset="predictors", title=paste("Predictors Only: (Ratio Cutoff>", RatioCutoff, ")", sep=""));
plot_predresp_line_diagram(pred_resp_mat,
        ratio_thres=0, signf_thres=SignifCutoff, a_name=AsPredAName, b_name=AsRespBName,
	subset="responders", title=paste("Responders Only: (Ratio Cutoff>", RatioCutoff, ")", sep=""));


#############################################################################	
# Close PDF output

dev.off();

cat("Writing statistics to text file...\n");
output_text_summary(combined_records_txt, OutputRoot);
cat("Done...\n");


#############################################################################

cat("Done.\n");
print(warnings());
q(status=0);
