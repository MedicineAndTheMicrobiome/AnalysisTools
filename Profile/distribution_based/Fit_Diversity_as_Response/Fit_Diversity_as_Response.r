#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
library(car);
options(useFancyQuotes=F);

source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');

params=c(
	"summary_file", "s", 1, "character",
	"factors", "f", 1, "character",
	"model_filename", "M", 2, "character",
	"required_var", "q", 2, "character",
	"reference_levels", "r", 2, "character",
	"outputroot", "o", 2, "character",
	"model_formula", "m", 2, "character",
	"testing_flag", "T", 2, "logical",
	"tag_name", "t", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table>\n",
	"	-f <factors>\n",
	"	[-M <model variables in a list file>]\n",
	"	[-q <required variables>]\n",
	"	[-r <reference levels file>]\n",
	"	[-o <output filename root>]\n",
	"	[-m \"<model formula string>\"]\n",
	"\n",
	"	[-T (testing flag)]\n",
	"	[-t <tag name>]\n",
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
	OutputRoot=gsub(".summary_table.tsv", "", OutputRoot);
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

ModelFilename="";
if(length(opt$model_filename)){
        ModelFilename=opt$model_filename;
}

RequiredFile="";
if(length(opt$required_var)){
        RequiredFile=opt$required_var;
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

cat("\n");
cat("Summary File : ", SummaryFile, "\n", sep="");
cat("Factors File: ", FactorsFile, "\n", sep="");
cat("Output File: ", OutputRoot, "\n", sep="");
cat("Reference Levels File: ", ReferenceLevelsFile, "\n", sep="");
cat("Required Variables File: ", RequiredFile, "\n", sep="");
cat("\n");

TestingMode=ifelse(length(opt$testing_flag)>0, T, F);

##############################################################################

load_factors=function(fname){
	cat("Loading Factors: ", fname, "\n");
	factors=data.frame(read.table(fname, sep="\t",  header=TRUE, row.names=1, check.names=FALSE, comment.char="", quote=""));
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
	factor_names=colnames(factors);

	for(i in 1:num_factors_to_relevel){
		relevel_target=relevel_names[i];

		if(length(intersect(relevel_target, factor_names))){
			target_level=ref_lev_mat[i, 1];
			tmp=factors[,relevel_target];
			if(length(intersect(target_level, tmp))){
				tmp=relevel(tmp, target_level);
    				factors[,relevel_target]=tmp;
			}else{
				cat("WARNING: Target level '", target_level,
					"' not found in '", relevel_target, "'!!!\n", sep="");
			}
		}else{
			cat("WARNING: Relevel Target Not Found: '", relevel_target, "'!!!\n", sep="");
		}
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

sig_char=function(val){
        if(!is.null(val) && !is.nan(val) && !is.na(val)){
                if(val <= .0001){ return("***");}
                if(val <= .001 ){ return("** ");}
                if(val <= .01  ){ return("*  ");}
                if(val <= .05  ){ return(":  ");}
                if(val <= .1   ){ return(".  ");}
        }
        return(" ");
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

# Load variables to require after NA removal
required_arr=NULL;
if(""!=RequiredFile){
	required_arr=scan(RequiredFile, what=character(), comment.char="#");
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

# Build model and select variables
if(ModelFormula!=""){
	cat("Model Formula: ", ModelFormula, "\n", sep="");
	model_string=paste("raw ~ ", ModelFormula);
	
}else if(ModelFilename!=""){
	cat("Model Variables Filename: ", ModelFilename, "\n");
	variables=scan(ModelFilename, what=character(), comment.char="#");
	shared=intersect(variables, factor_names);
	if(length(shared)!=length(variables)){
		cat("Missing variables specified in Model File:\n");
		print(setdiff(variables, shared));
	}
	model_string=paste("raw ~", paste(shared,  collapse=" + "), sep="");
}else{
	model_string= paste("raw ~", 
		paste(factor_names, collapse=" + "));
}

cat("Model String used for Regression: \n");
print(model_string);
model_var=get_var_from_modelstring(model_string);
print(model_var);
factors=factors[,model_var, drop=F];

# Handle NAs
#factors=remove_sample_or_factors_wNA_parallel(factors, required=NULL, num_trials=6400, num_cores=64, outfile=paste(OutputRoot, ".noNAs", sep=""));
remove_na_res=remove_sample_or_factors_wNA_parallel(factors, required=required_arr, num_trials=640000, num_cores=64, outfile=paste(OutputRoot, ".noNAs", sep=""));
factors=remove_na_res$factors;
samp_wo_nas=rownames(factors);
factor_names=colnames(factors);
num_factors=length(factor_names);
normalized=normalized[samp_wo_nas,];
num_samples=length(samp_wo_nas);

print(model_string);
print(factor_names);
model_string=rem_missing_var_from_modelstring(model_string, factor_names);
cat("New Model String with Factors with NAs removed:\n");
print(model_string);
num_model_variables=length(strsplit(model_string, "\\+")[[1]]);


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
pdf(paste(OutputRoot, rand, ".div_as_resp.pdf", sep=""), height=14, width=8.5);

plot_text(c(
	paste("Summary File : ", SummaryFile, sep=""),
	paste("Factors File: ", FactorsFile, sep=""),
	paste("Output File: ", OutputRoot, sep=""),
	paste("Reference Levels File: ", ReferenceLevelsFile, sep=""),
	paste("Required Variables File: ", RequiredFile, sep="")
));


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

# Set evenness for NAs to 0.  Assume evenness is very low, but it's degenerate
# because of insufficient sequencing depth.
evenness=div_mat[,"Evenness"];
div_mat[is.na(evenness),"Evenness"]=0;

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


		degenerates=which(is.nan(raw));
		if(length(degenerates)){
			cat("*************************************\n");
			cat("*  Degenerate diversities: \n");
			print(raw[degenerates]);
			cat("*************************************\n");
			raw[degenerates]=0;
		}

		zero_ix=(raw==0);
		if(any(zero_ix)){
			cat("Zeros found in diversity:\n");
			print(raw[zero_ix]);
			min_nonzero=min(raw[!zero_ix]);
			min_subst=min_nonzero/10;
			cat("Adding ", min_nonzero, " / 10 = ", min_subst, " to all responses...\n", sep="");
			raw=raw+min_subst;
			print(sort(raw));
		}
		
		bc=tryCatch({
			boxcox(as.formula(model_string), data=factors, 
				lambda=seq(lambda_start, lambda_end, length.out=SEARCH_FREQUENCY),
				plotit=FALSE);
			}, 
			error=function(cond){
				cat("Error finding Box-Cox transform lambda value.\n");
				print(cond);
				return(NULL);
			});

		if(is.null(bc)){
			approx_lambda=1;
			lambda_found=TRUE;
		}else{

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
	}


	if(approx_lambda!=1){
		# Rerun and plot the lambda search with a smaller increments
		par(mar=c(5,4,4,2));
		cat("Refining search around: ", 
			approx_lambda-BOX_COX_SEARCH_RANGE, " - ", approx_lambda+BOX_COX_SEARCH_RANGE, "\n");

		bc=boxcox(as.formula(model_string), data=factors, 
			lambda=seq(approx_lambda-BOX_COX_SEARCH_RANGE, approx_lambda+BOX_COX_SEARCH_RANGE, 
			length.out=20));

		title(main=div_names[i]);

		# Store fine grain results
		max_idx=which(bc$y==max(bc$y));
		lambda[i]=bc$x[max_idx];
	}else{
		cat("Warning: Box-Cox Lambda value not found.  Going with 1 (i.e. no transform)\n");
		cat("\n");
		cat("If the sample size was close to the number of variables in the model, it is likely\n");
		cat("that the model was overfit, so it was not possible to find lambda that minimizes the residuals.\n");
		cat("\n");
		cat("Num Model Variables: ", num_model_variables, "\n");
		cat("Num Samples: ", num_samples, "\n");
		cat("\n");
		lambda[i]=1;

	}
	
	cat(div_names[i], ": Box-Cox Transformation Lambda = ", lambda[i], "\n\n", sep="");

	# Apply transform to raw data
	if(lambda[i]==0){
		transformed[, div_names[i]]=log(raw);
	}else{
		transformed[, div_names[i]]=((raw^lambda[i])-1)/lambda[i];
	}
}

rownames(transformed)=rownames(div_mat);

mtext("Box-Cox Lambda Intervals", outer=T);

cat("Plotting transformed histograms...\n");
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

cat("Outputing factor levels...\n");
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

bin_continuous_values=function(values, num_bins=10){
        minv=min(values);
        maxv=max(values);
        range=maxv-minv;
        # Map values between 0 and 1
        prop=(values-minv)/range;
        # Scale value up to bin, and round, to quantize
        closest=round(prop*(num_bins-1),0);
        log10range=log10(range);
        trunc=signif(closest/(num_bins-1)*range+minv, 5)
	bin_range=mean(diff(sort(unique(trunc))))/2;

	lb=trunc-bin_range;
	ub=trunc+bin_range;
        # Remap values to original range and location
        return(paste("(", lb, ", ", ub, ")", sep=""));
}

plot_diversity_with_factors=function(raw, factors, model_string, stat_name, bin_cont=5){

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
	
	# Remove interaction terms
	interact_ix=grep(":", pred_arr);
	if(length(interact_ix)){
		pred_arr=pred_arr[-interact_ix]
	}

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
	predictors_per_plot=5;
	pred_names=colnames(factors);
	raw_adj_range=range(c(raw, adj_pos));

	num_plots=num_pred %/% predictors_per_plot + 1;

	for(plot_ix in 1:num_plots){
		# Plot the samples on the y axis
		plot(0, xlab="", ylab=stat_name, type="n", xaxt="n", main=stat_name,
			ylim=c(raw_adj_range[1], raw_adj_range[2]),
			xlim=c(-1, predictors_per_plot+1+extra_sample_space));

		# Plot 
		for(i in 1:num_values){
			lines(x=c(-1, -.75), y=c(raw[i], raw[i]));
			lines(x=c(-.75, -.25), y=c(raw[i], adj_pos[i]));

		}
		# Label samples
		text(zeros-.20, adj_pos, label=names(raw), pos=4, cex=.5);

		# Label predictors
		predictors_per_plot=min(predictors_per_plot, num_pred);

		if(predictors_per_plot>0){
			for(j in 1:predictors_per_plot){
				pred_ix=((plot_ix-1)*predictors_per_plot)+(j-1)+1;
			
				if(!is.na(pred_arr[pred_ix]) && !length(grep(":", pred_arr[pred_ix]))){

					fact_val=factors[, pred_arr[pred_ix]];
					uniq_val=unique(fact_val);

					if(length(uniq_val)>=bin_cont && !is.factor(fact_val)){
						factor=signif(fact_val, 4);
						coloring=as.numeric(as.factor(bin_continuous_values(fact_val, num_bins=bin_cont)));
					}else{
						factor=as.factor(fact_val);
						coloring=as.numeric(factor);
					}

					if(is.factor(factor)){
						factor=substr(factor, 1, 10);
					}

					text(zeros+j+extra_sample_space, adj_pos, label=factor, 
						col=coloring,
						cex=.5
						);
				}
				# label predictor/factor name
				lab_cex=min(1, 8/nchar(pred_arr[pred_ix]));
				text(j+extra_sample_space, raw_adj_range[2], 
					pred_arr[pred_ix], cex=lab_cex*.75, col="black", family="", font=2, pos=3);
			}
		}
	}

}

###############################################################################

plot_overlapping_histograms=function(raw, factors, model_string, title, bin_cont=5){

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

	#par(mfrow=c(4,1));
	layout_mat=matrix(c(1,1,2,3,3,4,5,5,6), ncol=1);
	layout(layout_mat);

	for(pix in 1:num_pred){
		cur_pred=pred_arr[pix];
		cat("  Working on: ", cur_pred, "\n");

		if(any(cur_pred==factor_names)){

			# Get levels for each value
			cur_fact_val=factors[,cur_pred];
			is_factor=is.factor(cur_fact_val);

			num_uniq=length(unique(cur_fact_val));
			if(num_uniq>=bin_cont && !is_factor){
				cur_fact_val=as.factor(bin_continuous_values(cur_fact_val, num_bins=bin_cont));
                        }else{
				cur_fact_val=as.factor(cur_fact_val);
			}

			#num_bins=nclass.Sturges(raw)*2;			
			#overall_hist=hist(raw, breaks=num_bins, plot=F);
			#print(overall_hist);
			#cat("  Num bins: ", num_bins, "\n");

			levels=levels(cur_fact_val);
			num_levels=length(levels);
			print(cur_fact_val);
			cat("Num Levels: ", num_levels, "\n");
			
			# Compute the density for each level
			dens_list=list(num_levels);
			level_samp_size=numeric();
			max_dens=0;
			for(lix in 1:num_levels){
				level_val=raw[cur_fact_val==levels[lix]];
				#hist_list[[levels[lix]]]=hist(level_val, breaks=overall_hist$breaks, plot=F);
				num_samp_at_level=length(level_val);
				level_samp_size[lix]=num_samp_at_level;
				if(num_samp_at_level==0){
					dens_list[[lix]]=NA;
				}else{
					if(num_samp_at_level==1){
						level_val=c(level_val, level_val);
					}

					dens_list[[lix]]=density(level_val);
					max_dens=max(max_dens, dens_list[[lix]]$y);	
				}
			}

			# Open a blank plot
			par(mar=c(4,3,2,1));
			plot(0,0, type="n", xlim=raw_range, ylim=c(0,max_dens), 
				main=cur_pred,
				xlab=title, ylab="Density",
				bty="l");

			# Draw the curves for each factor level
			for(lix in 1:num_levels){
				dens=dens_list[[lix]];
				if(!any(is.na(dens))){
					points(dens, type="l", col=lix);
				}
			}

			# Plot the legend for each factor
			legend_labels=character();
			for(lix in 1:num_levels){
				legend_labels[lix]=paste(levels[lix], " [n=", level_samp_size[lix], "] ", sep="");
			}

			par(mar=c(0,0,0,0));
			plot(0,0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", bty="n");
			legend(0,1, legend=legend_labels, fill=1:num_levels, bty="n");

		}else{
			cat("Not ploting interactions...\n");
		}
	}

	par(orig.par);
	
}

##############################################################################

# Matrix to store R^2's
cat("Allocating R^2 Matrix...\n");
rsqrd_mat=matrix(0, nrow=num_div_idx, ncol=2);
colnames(rsqrd_mat)=c("R^2", "Adj. R^2");
rownames(rsqrd_mat)=div_names;

# Build model matrix so we know how many coefficients to expect.
if(ModelFormula==""){
	model_string= paste("trans ~", 
	paste(factor_names, collapse=" + "));
}

cat("\nFitting this regression model: ", model_string, "\n");

trans=rep(0, nrow(factors));
model_matrix=model.matrix(as.formula(model_string), data=factors);
num_coeff=ncol(model_matrix);

# Matrices to store coefficents and p-values
cat("Allocating coeff and pval matrices...\n");
coeff_matrix=matrix(NA, nrow=num_coeff, ncol=num_div_idx);
pval_matrix=matrix(NA, nrow=num_coeff, ncol=num_div_idx);

colnames(coeff_matrix)=div_names;
rownames(coeff_matrix)=colnames(model_matrix);

colnames(pval_matrix)=div_names;
rownames(pval_matrix)=colnames(model_matrix);


# Fit regression model and ANOVA analysis
for(i in 1:num_div_idx){

	cat("Working on: ", div_names[i], "\n", sep="");

	raw=div_mat[, div_names[i]];

	plot_diversity_with_factors(raw, factors, model_string, div_names[i], 6);
	plot_overlapping_histograms(raw, factors, model_string, title=div_names[i], 6);

	trans=transformed[, div_names[i]];
	print(model.matrix(as.formula(model_string), data=factors));
	fit=lm(as.formula(model_string), data=factors);
	print(fit);

	summ=summary(fit);
	print(summ);

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
		"",
		paste("Diversity Mean: ", mean(raw)),
		paste("Diversity Stderr: ", sd(raw)/sqrt(length(raw))),
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


dev.off();

##############################################################################

fh=file(paste(OutputRoot, ".div_as_resp.regr_stats.tsv", sep=""), "w");

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
# Write coefficient p-values to file

write.table(t(pval_matrix), file=paste(OutputRoot, ".div_as_resp.pvals.tsv", sep=""),
        sep="\t", quote=F, col.names=NA, row.names=T);

write.table(t(coeff_matrix), file=paste(OutputRoot, ".div_as_resp.coefs.tsv", sep=""),
        sep="\t", quote=F, col.names=NA, row.names=T);
	
##############################################################################
# Write diversity matrix to file

write.table(div_mat, file=paste(OutputRoot, ".diversity_indices.tsv", sep=""),
	sep="\t", quote=F, col.names=NA, row.names=T);

##############################################################################
# Compute multivariate regression for a MANOVA

#print(transformed);
#print(model_string);
#print(factors);
#print(model_string);

trans=transformed;
mv_fit=lm(as.formula(model_string), data=factors);

try_res=try({manova_res=anova(mv_fit)});

if(class(try_res)=="try-error"){
	manova_res=NULL;
	cat("MANOVA Failed.\n");
}else{
	cat("MANOVA:\n");
	print(manova_res);
}

if(length(manova_res)){
        num_variables=nrow(manova_res)-1;
        outmat=matrix("", nrow=num_variables, ncol=3);
        colnames(outmat)=c(TagName, "Pr(>F)", "Signf");
        varnames=unlist(rownames(manova_res));
        pvals=unlist(manova_res["Pr(>F)"]);
        outmat[,TagName]=varnames[1:num_variables];
        outmat[,"Pr(>F)"]=sprintf("%4.4f", pvals[1:num_variables]);
	outmat[,"Signf"]=sapply(pvals[1:num_variables], sig_char);
}else{
        outmat=matrix("-", nrow=1, ncol=2);
        colnames(outmat)=c(TagName, "Pr(>F)");
}

write.table(outmat, file=paste(OutputRoot, ".div_as_resp.anova.summary.tsv", sep=""),
        sep="\t", quote=F, col.names=T, row.names=F);


##############################################################################

cat("Done.\n");
print(warnings());
q(status=0);
