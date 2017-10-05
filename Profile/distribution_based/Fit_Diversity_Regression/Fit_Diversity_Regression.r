#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
library(car);
options(useFancyQuotes=F);

params=c(
	"summary_file", "s", 1, "character",
	"factors", "f", 1, "character",
	"reference_levels", "r", 2, "character",
	"outputroot", "o", 2, "character",
	"model_formula", "m", 2, "character",
	"testing_flag", "t", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table>\n",
	"	-f <factors>\n",
	"	[-r <reference levels file>]\n",
	"	[-o <output filename root>]\n",
	"	[-m \"<model formula string>\"]\n",
	"\n",
	"	[-t (testing flag)]\n",
	"\n",
	"Computes the diversity measuares for each sample, performs a box-cox\n",
	"transformation and then performs a linear regression based on the\n",
	"factors in the factor file.\n",
	"\n",
	"For testing purposes, you can append the factor name with 'IGNORE.' and\n",
	"the factor will be ignored.\n",
	"\n",
	"The format of the reference levels file is: \n",
	"	<factor name1>\\t<reference level1>\\n\n",
	"	<factor name2>\\t<reference level2>\\n\n",
	"	...\n",
	"	<factor nameN>\\t<reference levelN>\\n\n",
	"\n");

if(!length(opt$summary_file) || !length(opt$factors)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputRoot=gsub(".summary_table.xls", "", opt$summary_file);
}else{
	OutputRoot=opt$outputroot;
}

if(!length(opt$reference_levels)){
	ReferenceLevelsFile="";
}else{
	ReferenceLevelsFile=opt$reference_levels;
}

if(!length(opt$model_formula)){
	ModelFormula="";
}else{
	ModelFormula=opt$model_formula;
}

SummaryFile=opt$summary_file;
FactorsFile=opt$factors;

cat("\n");
cat("Summary File : ", SummaryFile, "\n", sep="");
cat("Factors File: ", FactorsFile, "\n", sep="");
cat("Output File: ", OutputRoot, "\n", sep="");
cat("Reference Levels File: ", ReferenceLevelsFile, "\n", sep="");
cat("\n");

TestingMode=ifelse(length(opt$testing_flag)>0, T, F);

if(ModelFormula!=""){
	cat("Model Formula: ", ModelFormula, "\n", sep="");

	RegressionFormula=gsub("\\+\\s*Error\\(.*\\)", "", ModelFormula);
	cat("Regression Formula: ", RegressionFormula, "\n", sep="");
}

##############################################################################

load_factors=function(fname){
	cat("Loading Factors: ", fname, "\n");
	factors=data.frame(read.table(fname,  header=TRUE, row.names=1, check.names=FALSE, comment.char="", quote=""));
	factor_names=colnames(factors);

	ignore_idx=grep("^IGNORE\\.", factor_names);

	if(length(ignore_idx)!=0){
		return(factors[-ignore_idx]);
	}else{
		return(factors);
	}
}

load_summary_file=function(fname){
	cat("Loading Summary Table: ", fname, "\n");
	inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", row.names=1))
	counts_mat=inmat[,2:(ncol(inmat))];
	return(counts_mat);
}

load_reference_levels_file=function(fname){
	cat("Loading Reference Levels: ", fname, "\n");
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
	for(i in 1:num_factors_to_relevel){
		tmp=factors[,relevel_names[i]];
		#print(tmp);
		tmp=relevel(tmp, ref_lev_mat[i, 1]);
		#print(tmp);
		factors[,relevel_names[i]]=tmp;
	}
	return(factors);
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

plot_correl_heatmap=function(mat, title="", noPrintZeros=F, guideLines=F){

        if(is.null(dim(mat))){
                cat(title, " Matrix is NULL.  No heatmap generated.\n");
                return();
        }

        cat("Plotting: ", title, "\n");

        par(family="Courier");
        par(oma=c(10, 10, 1, .5));
        par(mar=c(5.1, 4.1, .5, .5));

        # Generate colors from red to blue
        colors=(rainbow(2^16, start=0, end=0.65));

        # Remember that rows and columsn are reversed in the image
        image(1:nrow(mat),1:ncol(mat), mat,
                xaxt="n", yaxt="n",
                xlab="", ylab="",
                col=colors
        );

        # Pad strings
        cnames=paste(colnames(mat), " ", sep="");
        rnames=paste(rownames(mat), " ", sep="");

        # Get longest length of each column or row label
        cname_max_len=max(nchar(cnames));
        rname_max_len=max(nchar(rnames));

        # Get the number of rows and columns
        ncols=ncol(mat);
        nrows=nrow(mat);

        cscale=min(c(25/cname_max_len, 25/ncols));
        rscale=min(c(25/rname_max_len, 25/nrows));

        max_width=max(nchar(sprintf("%.2f",mat)));
        cell_cex=(3.5/max_width)*sqrt(min(c(cscale, rscale))^2);

        for(i in 1:nrow(mat)){
                for(j in 1:ncol(mat)){

                        if(!is.na(mat[i,j]) && (noPrintZeros && mat[i,j]==0)){
                                # Skip
                        }else{
                                str=sprintf("%.2f",mat[i,j]);
                                str=gsub("^0\\.",".", str);
                                text(i,j,labels=str, cex=cell_cex, srt=45);
                        }
                }
        }

        # Plot guidelines
        if(guideLines){

                splits=c(4,5,6,7);

                if(nrows>8){
                        h_remainder=nrows %% splits;
                        best_h_split=min(which(min(h_remainder)==h_remainder));
                        num_hor_splits=nrows %/% best_h_split;
                        h_line_pos=seq(5, 5*(num_hor_splits-1), 5)+.5;
                        abline(h=h_line_pos, col="black", lty="dashed");
                        abline(h=h_line_pos, col="white", lty="dotted");
                }

                if(ncols>8){
                        v_remainder=ncols %% splits;
                        best_v_split=min(which(min(v_remainder)==v_remainder));
                        num_ver_splits=ncols %/% best_v_split;
                        v_line_pos=seq(5, 5*(num_ver_splits-1), 5)+.5;
                        abline(v=v_line_pos, col="black", lty="dashed");
                        abline(v=v_line_pos, col="white", lty="dotted");
                }

        }

        # Plot the labels
        mtext(cnames, at=1:ncols, side=2, las=2, cex=cscale);
        mtext(rnames, at=1:nrows, side=1, las=2, cex=rscale);

        # Plot the title
        mtext(title, line=0, at=nrows*.5, side=3, font=2);

}

##############################################################################

# Load matrix
counts=load_summary_file(SummaryFile);
num_taxa=ncol(counts);
num_samples=nrow(counts);
#print(counts);

# Shorten names:
full_names=colnames(counts);
splits=strsplit(full_names, " ");
short_names=character();
for(i in 1:length(full_names)){
	short_names[i]=tail(splits[[i]], 1);
}
colnames(counts)=short_names;

# Normalize
normalized=normalize(counts);
#print(normalized);

##############################################################################

# Load factors
factors=load_factors(FactorsFile);
factor_names=colnames(factors);
num_factors=ncol(factors);
factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);

cat("\n");
cat(num_factors, " Factor(s) Loaded:\n", sep="");
#print(factor_names);
cat("\n");

#print(factors);
#summary(factors);

# Relevel factor levels
if(ReferenceLevelsFile!=""){
	ref_lev_mat=load_reference_levels_file(ReferenceLevelsFile)
	factors=relevel_factors(factors, ref_lev_mat);
}else{
	cat("No Reference Levels File specified.\n");
}

continuous_factors=factors;
is_continous_factor=logical(num_factors);

for(f in 1:num_factors){
	level_info=levels(factors[,f]);
	is_continous_factor[f]=is.null(level_info);

	if(is_continous_factor[f]){
		# do nothing
	}else if(length(level_info)==2){
		# Convert two level factors to numeric
		is_continous_factor[f]=TRUE;
		continuous_factors[,f]=as.integer(continuous_factors[,f]);
	}else{
		is_continous_factor[f]=FALSE;
	}
}

continuous_factors=continuous_factors[,is_continous_factor, drop=F];
#print(continuous_factors);
factor_correlations=cor(continuous_factors);

##############################################################################

# Reconcile factors with samples
factor_sample_ids=rownames(factors);
counts_sample_ids=rownames(counts);

#print(factor_sample_id);
#print(counts_sample_id);

shared_sample_ids=intersect(factor_sample_ids, counts_sample_ids);
num_shared_sample_ids=length(shared_sample_ids);
num_factor_sample_ids=length(factor_sample_ids);
num_counts_sample_ids=length(counts_sample_ids);

cat("Num counts sample IDs: ", num_counts_sample_ids, "\n");
cat("Num factor sample IDs: ", num_factor_sample_ids, "\n");
cat("Num shared sample IDs: ", num_shared_sample_ids, "\n");
cat("\n");

cat("Samples missing from count information:\n");
print(setdiff(factor_sample_ids, counts_sample_ids));
cat("\n");
cat("Samples missing from factor information:\n");
print(setdiff(counts_sample_ids, factor_sample_ids));
cat("\n");
cat("Total samples shared: ", num_shared_sample_ids, "\n");

shared_sample_ids=sort(shared_sample_ids);

# Reorder data by sample id
normalized=normalized[shared_sample_ids,];
num_samples=nrow(normalized);
factors=factors[shared_sample_ids,, drop=F];

##############################################################################

plot_text=function(strings){
	par(mfrow=c(1,1));
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

##############################################################################

if(TestingMode){
	rand=sprintf(".%04i", sample(1000,1));
}else{
	rand="";
}
pdf(paste(OutputRoot, rand, ".div_reg.pdf", sep=""), height=14, width=8.5);

# Output the factor correlations
if(nrow(factor_correlations)>0){
	plot_correl_heatmap(factor_correlations, title="Factor Correlations");
}else{
	plot_text("No factor correlation heatmap generated, because the number of ordinal factors was 0.\n");
}


text=character();

##############################################################################
# Compute diversity indices

div_names=c("Tail", "Shannon", "Simpson", "Evenness", "SimpsonsRecip");
num_div_idx=length(div_names);

div_mat=matrix(0, nrow=num_samples, ncol=num_div_idx);
colnames(div_mat)=div_names;
rownames(div_mat)=rownames(normalized);

cat("Computing diversity indices across samples.\n");
for(i in 1:num_samples){
	curNorm=normalized[i,];
	zeroFreeNorm=curNorm[curNorm>0]
	div_mat[i,"Tail"]=tail_statistic(zeroFreeNorm);
	div_mat[i,"Shannon"]=-sum(zeroFreeNorm*log(zeroFreeNorm));
	div_mat[i,"Simpson"]=1-sum(curNorm^2);
	div_mat[i,"Evenness"]=div_mat[i,"Shannon"]/log(length(zeroFreeNorm));
	div_mat[i,"SimpsonsRecip"]=1/sum(curNorm^2);
}

cat("Plotting histograms of raw diversity indices.\n");
par(mfrow=c(3,2));
par(oma=c(1, 1, 1, 1));
par(mar=c(5,4,4,2));
for(i in 1:num_div_idx){
	hist(div_mat[, div_names[i]], main=div_names[i], xlab=div_names[i], 
		breaks=15);
}
mtext("Sample Diversity Indices Distributions", outer=T);

##############################################################################
# Perform Box-Cox transformation

lambda=numeric(num_div_idx);
transformed=matrix(0, nrow=num_samples, ncol=num_div_idx);
colnames(transformed)=div_names;

cat("Computing and plotting Box-Cox transformations.\n\n");
BOX_COX_SEARCH_RANGE=2;
SEARCH_FREQUENCY=10;
par(mfrow=c(3,2));
for(i in 1:num_div_idx){

	cat("Computing lambda for: ", div_names[i], "\n", sep="");

	raw=div_mat[, div_names[i]];

	if(ModelFormula==""){
		model_string= paste("raw ~", 
			paste(factor_names, collapse=" + "));
	}else{
		model_string= paste("raw ~", RegressionFormula);
	}

	#cat("Model used in Box-Cox transformation: ", model_string, "\n", sep="");

	# Look for lambda approximately
	lambda_found=FALSE;
	lambda_start=-BOX_COX_SEARCH_RANGE;
	lambda_end=BOX_COX_SEARCH_RANGE;
	search_trial=0;
	while(!lambda_found){

		search_trial=search_trial+1;

		# Perform trial boxcox search
		cat("[", search_trial, "] Search range: ", lambda_start, " to ", lambda_end, "\n", sep="");
		cat("Model: ", model_string, "\n");
		bc=boxcox(as.formula(model_string), data=factors, 
			lambda=seq(lambda_start, lambda_end, length.out=SEARCH_FREQUENCY),
			plotit=FALSE
		);

		# Determine where to expand search depending if peak is on left or right
		max_idx=which(bc$y==max(bc$y));
		if(max_idx==1){
			#lambda_end=lambda_start+1;
			lambda_start=lambda_start-(BOX_COX_SEARCH_RANGE^search_trial);
		}else if(max_idx==length(bc$y)){
			#lambda_start=lambda_end-1;
			lambda_end=lambda_end+BOX_COX_SEARCH_RANGE^search_trial;
		}else{
			search_tolerance=(lambda_end-lambda_start)/SEARCH_FREQUENCY;
			approx_lambda=bc$x[max_idx];
			lambda_found=TRUE;
			cat("Lambda found around: ", approx_lambda, "\n");	
			cat("   Search tolerance: ", search_tolerance, "\n");	
		}
	}

	# Rerun and plot the lambda search with a smaller increments
	par(mar=c(5,4,4,2));
	cat("Refining search around: ", approx_lambda-search_tolerance, " - ", approx_lambda+search_tolerance, "\n");
	bc=boxcox(as.formula(model_string), data=factors, 
		lambda=seq(approx_lambda-search_tolerance, approx_lambda+search_tolerance, length.out=40));
	title(main=div_names[i]);

	# Store find grain results
	max_idx=which(bc$y==max(bc$y));
	lambda[i]=bc$x[max_idx];
	
	cat(div_names[i], ": Box-Cox Transformation Lambda = ", 
		lambda[i], "\n\n", sep="");

	# Apply transform to raw data
	if(lambda[i]==0){
		transformed[, div_names[i]]=log(raw);
	}else{
		transformed[, div_names[i]]=((raw^lambda[i])-1)/lambda[i];
	}
}
mtext("Box-Cox Lambda Intervals", outer=T);

par(mfrow=c(3,2));
par(oma=c(1, 1, 1, 1));
par(mar=c(5,4,4,2));
for(i in 1:num_div_idx){
	hist(transformed[, div_names[i]], 
		main=paste(div_names[i], " (lambda=", sprintf("%3.2f", lambda[i]), ")", sep=""),
		xlab=div_names[i], 
		breaks=15);
}
mtext("Transformed Sample Diversity Indices Distributions", outer=T);

##############################################################################
# Output reference factor levels

text=character();
text[1]="Reference factor levels:";
text[2]="";

for(i in 1:num_factors){
	fact_levels=levels(factors[,i]);
	if(!is.null(fact_levels)){
		fact_info=paste(factor_names[i], ": ", fact_levels[1], sep="");	
	}else{
		fact_info=paste(factor_names[i], ": None (ordered factor)", sep="");	
	}
	text=c(text, fact_info);
}

text=c(text, "");
text=c(text, paste("Number of Samples: ", num_samples, sep=""));
text=c(text, "");
text=c(text, "Description of Factor Levels and Samples:");
text=c(text, capture.output(summary(factors)));

plot_text(text);

##############################################################################

plot_diversity_with_factors=function(raw, factors, model_string){

	palette(c(
		"red",
		"green",
		"blue",
		"orange",
		"purple",
		"black",
		"brown",
		"darkgoldenrod3"
	));
	
	# Extract out predictors
	predictor_string=strsplit(model_string, "~")[[1]][2];
	pred_arr=strsplit(predictor_string, "\\+")[[1]];
	pred_arr=gsub(" ", "", pred_arr);
	num_pred=length(pred_arr);
	num_values=length(raw);
	raw_range=range(raw);

	# Sort records by magnitude of raw values
	sort_ix=order(raw);
	raw=raw[sort_ix];
	factors=factors[sort_ix, , drop=F];

	par(mar=c(5,5,5,5));
	zeros=rep(0, num_values);
	
	min_spacing=(raw_range[2]-raw_range[1])/60;
	cat("Min spacing: ", min_spacing, "\n");
	#abline(h=seq(raw_range[1], raw_range[2], length.out=50));

	adj=min_spacing/10;
	cat("Adjustment Increment: ", adj, "\n");
	adj_pos=raw;
	
	for(adj_ix in 1:10000){
		adjusted=F
		for(forw in 1:(num_values-1)){
			if(abs(adj_pos[forw]-adj_pos[forw+1])<min_spacing){
				adj_pos[forw]=adj_pos[forw]-adj;
				adjusted=T;
				break;
			}	
		}
		for(rev in (num_values:2)){
			if(abs(adj_pos[rev]-adj_pos[rev-1])<min_spacing){
				adj_pos[rev]=adj_pos[rev]+adj;
				adjusted=T;
				break;
			}	
		}
		if(!adjusted){
			cat("Done adjusting at iteration: ", adj, "\n");
			break;	
		}
	}
	
	#cat("Raw:\n");
	#print(raw);
	#cat("Adj:\n");
	#print(adj_pos);

	extra_sample_space=2;

	raw_adj_range=range(c(raw, adj_pos));
	plot(0, xlab="", ylab=div_names[i], type="n", xaxt="n", main=div_names[i],
		ylim=c(raw_adj_range[1], raw_adj_range[2]),
		xlim=c(-1, num_pred+1+extra_sample_space));

	for(i in 1:num_values){
		lines(x=c(-1, -.75), y=c(raw[i], raw[i]));
		lines(x=c(-.75, -.25), y=c(raw[i], adj_pos[i]));

	}
	# Label samples
	text(zeros-.20, adj_pos, label=names(raw), pos=4, cex=.5);

	# Label predictors
	tot_levels=0;	
	for(j in 1:num_pred){
		if(!length(grep(":", pred_arr[j]))){
			factor=as.factor(factors[, pred_arr[j]]);

			abbreviate=substr(factors[,pred_arr[j]], 1, 10);

			text(zeros+j+extra_sample_space, adj_pos, label=abbreviate, 
				col=tot_levels+as.numeric(factor),
				cex=.5
				);
			tot_levels=tot_levels+length(levels(factor));
		}
	}
}

###############################################################################

plot_overlapping_histograms=function(raw, factors, model_string, title){

	orig.par=par(no.readonly=T);
	palette(c(
		"red",
		"green",
		"blue",
		"orange",
		"purple",
		"black",
		"brown",
		"darkgoldenrod3"
	));

	# Extract out predictors
	predictor_string=strsplit(model_string, "~")[[1]][2];
	pred_arr=strsplit(predictor_string, "\\+")[[1]];
	pred_arr=gsub(" ", "", pred_arr);
	num_pred=length(pred_arr);
	num_values=length(raw);
	raw_range=range(raw);

	cat("Num Predictors:", num_pred, "\n");
	cat("Range: ", raw_range[1], " - ", raw_range[2], "\n");

	factor_names=colnames(factors);

	par(mfrow=c(4,1));

	for(pix in 1:num_pred){
		cur_pred=pred_arr[pix];
		cat("  Working on: ", cur_pred, "\n");

		if(any(cur_pred==factor_names)){

			# Get levels for each value
			cur_fact_val=factors[,cur_pred];
			is_factor=is.factor(cur_fact_val);
			if(is_factor){

				#num_bins=nclass.Sturges(raw)*2;			
				#overall_hist=hist(raw, breaks=num_bins, plot=F);
				#print(overall_hist);
				#cat("  Num bins: ", num_bins, "\n");

				levels=levels(cur_fact_val);
				num_levels=length(levels);
				
				# Compute the density for each level
				dens_list=list();
				max_dens=0;
				for(lix in 1:num_levels){
					level_val=raw[cur_fact_val==levels[lix]];
					#hist_list[[levels[lix]]]=hist(level_val, breaks=overall_hist$breaks, plot=F);
					dens_list[[lix]]=density(level_val);
					max_dens=max(max_dens, dens_list[[lix]]$y);	
				}

				# Plot the legend for each factor
				plot(0,0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", bty="n");
				legend(0,1, legend=levels, fill=1:num_levels, bty="n");

				# Open a blank plot
				plot(0,0, type="n", xlim=raw_range, ylim=c(0,max_dens), 
					main=cur_pred,
					xlab=title, ylab="Density");

				# Draw the curves for each factor level
				for(lix in 1:num_levels){
					points(dens_list[[lix]], type="l", col=lix);
				}

			}else{
				cat("Not plotting continuous metadata.\n");
			}
		}else{
			cat("Not ploting interactions...\n");
		}
	}

	par(orig.par);
	
}

##############################################################################

# Matrix to store R^2's
rsqrd_mat=matrix(0, nrow=num_div_idx, ncol=2);
colnames(rsqrd_mat)=c("R^2", "Adj. R^2");
rownames(rsqrd_mat)=div_names;

# Build model matrix so we know how many coefficients to expect.
if(ModelFormula==""){
	model_string= paste("trans ~", 
	paste(factor_names, collapse=" + "));
}else{
	model_string= paste("trans ~", 
		RegressionFormula);
}
cat("\nFitting this regression model: ", model_string, "\n");

trans=rep(0, nrow(factors));
model_matrix=model.matrix(as.formula(model_string), data=factors);
num_coeff=ncol(model_matrix);

# Matrices to store coefficents and p-values
coeff_matrix=matrix(NA, nrow=num_coeff, ncol=num_div_idx);
pval_matrix=matrix(NA, nrow=num_coeff, ncol=num_div_idx);

colnames(coeff_matrix)=div_names;
rownames(coeff_matrix)=colnames(model_matrix);

colnames(pval_matrix)=div_names;
rownames(pval_matrix)=colnames(model_matrix);


for(i in 1:num_div_idx){

	cat("Working on: ", div_names[i], "\n", sep="");

	raw=div_mat[, div_names[i]];
	plot_diversity_with_factors(raw, factors, model_string);
	plot_overlapping_histograms(raw, factors, model_string, title=div_names[i]);

	trans=transformed[, div_names[i]];
	fit=lm(as.formula(model_string), data=factors);

	summ=summary(fit);

	rsqrd_mat[i, 1]=summ$r.squared;
	rsqrd_mat[i, 2]=summ$adj.r.squared;

	# Output ANOVA
	if(ModelFormula!=""){
		aov_model=paste("trans ~", ModelFormula);
	}else{
		aov_model=model_string;
	}

	aov_sumtext=capture.output(summary(aov(as.formula(aov_model), data=factors)));

	plot_text(c(
		paste("ANOVA for: ", div_names[i]),
		"",
		aov_model,
		"", aov_sumtext));

	# Output regression coefficients and significance
	sumtext=(capture.output(summary(fit)));
	plot_text(c(
		paste("Multiple Regression for: ", div_names[i]), 
		"",
		model_string,
		"", sumtext));

	# Generate marginal model plots
	mmps(fit);



	sum_fit=summary(fit);

	computable_coef_names=names(sum_fit$coefficients[,"Estimate"]);

	coeff_matrix[computable_coef_names,i]=sum_fit$coefficients[computable_coef_names,"Estimate"];
	pval_matrix[computable_coef_names,i]=sum_fit$coefficients[computable_coef_names,"Pr(>|t|)"];
}

# R-squred summary
print(rsqrd_mat);
plot_correl_heatmap(t(rsqrd_mat), title="R^2 Values");

# Remove intercept from matrix
intercept_ix=which(rownames(coeff_matrix)=="(Intercept)");
coeff_matrix=coeff_matrix[-intercept_ix, , drop=F];
pval_matrix=pval_matrix[-intercept_ix, , drop=F];

plot_correl_heatmap(coeff_matrix, title="Coefficients");
plot_correl_heatmap(pval_matrix, title="P-values");

significant_coeff=(pval_matrix<0.05)*coeff_matrix;
plot_correl_heatmap(significant_coeff, title="Significant Coefficients", noPrintZeros=T);


##############################################################################

fh=file(paste(OutputRoot, ".coeff.tsv", sep=""), "w");

# Output Header
cat(file=fh, paste(c("Coefficients", "Estimates:", div_names, "p-values:", div_names), collapse="\t"));
cat(file=fh, "\n");

# Output values
coeff_names=rownames(pval_matrix);
for(i in 1:length(coeff_names)){
	cat(file=fh, 
		coeff_names[i], 
		"",
		paste(coeff_matrix[i,], collapse="\t"),
		"",
		paste(pval_matrix[i,], collapse="\t"),
		sep="\t");
	cat(file=fh, "\n");
}

close(fh);

##############################################################################
	
dev.off();

##############################################################################

cat("Done.\n");
print(warnings());
q(status=0);
