#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
options(useFancyQuotes=F);

source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');

params=c(
	"summary_file", "s", 1, "character",
	"factors", "f", 1, "character",
	"outputroot", "o", 2, "character",

	"model_formula", "m", 2, "character",

	"response_var", "y", 2, "character",
	"covariates_var", "c", 2, "character",
	"required_var", "q", 2, "character",

	"reference_levels", "l", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary table file (.tsv)>\n",
	"	-f <factors/metadata file>\n",
	"	[-o <output filename root>]\n",
	"\n",
	"	Model Specification:\n",
	"	[-m \"<model formula string>\"]\n",
	"\n",
	"	Variable Specification:\n",
	"	[-y <response variables>]\n",
	"	[-c <covariates variables>]\n",
	"	[-q <required variables>]\n",
	"	[-l <reference levels file>]\n",
	"\n",
	"This script will fit the following types of models:\n",
	"	<responses> = <covariates> + <microbiome diversity>\n",	
	"\n",
	"The analysis will cycle through the various diversity\n",
	"indices.\n",
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
ResponseFile="";
if(length(opt$response_var)){
        ResponseFile=opt$response_var;
}
CovariatesFile="";
if(length(opt$covariates_var)){
        CovariatesFile=opt$covariates_var;
}

SummaryFile=opt$summary_file;
FactorsFile=opt$factors;


summary_text=c();

out=capture.output({
		cat("\n");
		cat("Summary File : ", SummaryFile, "\n", sep="");
		cat("Factors File: ", FactorsFile, "\n", sep="");
		cat("Output File: ", OutputRoot, "\n", sep="");
		cat("\n");
		cat("Reference Levels File: ", ReferenceLevelsFile, "\n", sep="");
		cat("Required Variables File: ", RequiredFile, "\n", sep="");
		cat("Response Variables File: ", ResponseFile, "\n", sep="");
		cat("Covariates Variables File: ", CovariatesFile, "\n", sep="");
		cat("\n");
});
cat(out, sep="\n");
summary_text=c(summary_text, "\n", out);

##############################################################################

load_factors=function(fname){
	cat("Loading Factors: ", fname, "\n");
	factors=data.frame(read.table(fname,  sep="\t", header=TRUE, row.names=1, 
		check.names=FALSE, comment.char="", quote=""));
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
	inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, 
		comment.char="", row.names=1))
	counts_mat=inmat[,2:(ncol(inmat))];
	return(counts_mat);
}

load_reference_levels_file=function(fname){
	cat("Loading Reference Levels: ", fname, "\n");
	inmat=as.matrix(read.table(fname, sep="\t", header=F, check.names=FALSE, 
		comment.char="#", row.names=1))
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

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, deci_pts=4,
        label_zeros=T, counts=F, value.cex=1,
        plot_col_dendr=F,
        plot_row_dendr=F
){
	# Remove any rows with NAs in them
	any_nas=apply(mat, 1, function(x){all(is.na(x))});
	mat=mat[!any_nas,,drop=F];
	any_nas=apply(mat, 2, function(x){all(is.na(x))});
	mat=mat[,!any_nas,drop=F];

        num_row=nrow(mat);
        num_col=ncol(mat);

	if(num_row==0 || num_col==0){
		plot(0, type="n", xlim=c(-1,1), ylim=c(-1,1), xaxt="n", yaxt="n", bty="n", xlab="", ylab="",
                        main=title);
                text(0,0, "No data to plot...");
		return();
	}

	if(num_row==1 || num_col==1){
		plot_row_dendr=F;
		plot_col_dendr=F;
	}

        row_names=rownames(mat);
        col_names=colnames(mat);

        orig.par=par(no.readonly=T);

	cat("Painting Matrix for: ", title, "\n"); 
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
                        #noop;
                }else{
                        in_mat=t(in_mat);
                }
		dendist=dist(in_mat);

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
        plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), xaxt="n", yaxt="n", 
		bty="n", xlab="", ylab="");
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
                                text(x-.5, y-.5, text_lab, srt=atan(num_col/num_row)/pi*180, 
					cex=value.cex, font=2);
                        }
                }
        }

        ##################################################################################################

        par(mar=c(0, 0, 0, 0));

        if(plot_row_dendr && plot_col_dendr){
                rdh=attributes(row_dendr[["tree"]])$height;
                cdh=attributes(col_dendr[["tree"]])$height;
                plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", bty="n", 
			xlim=c(rdh, 0));
                plot(col_dendr[["tree"]], leaflab="none",xaxt="n", yaxt="n", bty="n", ylim=c(0, cdh));
                plot(0,0, type="n", bty="n", xaxt="n", yaxt="n");
                #text(0,0, "Placeholder");
        }else if(plot_row_dendr){
                rdh=attributes(row_dendr[["tree"]])$height;
                plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", bty="n", 
			xlim=c(rdh, 0));
                #text(0,0, "Row Dendrogram");
        }else if(plot_col_dendr){
                cdh=attributes(col_dendr[["tree"]])$height;
                plot(col_dendr[["tree"]], leaflab="none", xaxt="n", yaxt="n", bty="n", ylim=c(0, cdh));
                #text(0,0, "Col Dendrogram");
        }

        par(orig.par);

}

plot_pred_vs_obs=function(lmfit_ful, lmfit_red, title=""){

	par.orig=par(no.readonly=T);

	#print(lmfit$fitted.values);
	#print(lmfit$y);

	obs=lmfit_ful$y;

	pred_red=lmfit_red$fitted.values;
	pred_ful=lmfit_ful$fitted.values;

	sum_ful=summary(lmfit_ful);
	sum_red=summary(lmfit_red);

	if(is.null(ncol(obs))){
		obs=matrix(obs, nrow=length(obs), ncol=1, dimnames=list(names(pred_red), "obs"));
		pred_red=matrix(pred_red, nrow=length(pred_red), ncol=1, 
			dimnames=list(names(pred_red), "y"));
		pred_ful=matrix(pred_ful, nrow=length(pred_ful), ncol=1, 
			dimnames=list(names(pred_ful), "y"));

		sum_ful=list(sum_ful);
		sum_red=list(sum_red);
	}

	num_mv_resp=ncol(pred_ful);
	response_names=colnames(pred_ful);

	plot_row=2;
	plot_col=4;
	plots_per_page=(plot_row*plot_col)/2;
	par(mfcol=c(plot_row, plot_col));
	
	for(resp_ix in 1:num_mv_resp){

		obs_cur=obs[,resp_ix];
		rngs=range(c(obs_cur, pred_red[,resp_ix], pred_ful[,resp_ix]));	

		# P-values
		red_pval=1-pf(sum_red[[resp_ix]]$fstatistic["value"], 
			sum_red[[resp_ix]]$fstatistic["numdf"], sum_red[[resp_ix]]$fstatistic["dendf"]);
		ful_pval=1-pf(sum_ful[[resp_ix]]$fstatistic["value"], 
			sum_ful[[resp_ix]]$fstatistic["numdf"], sum_ful[[resp_ix]]$fstatistic["dendf"]);

		red_pval=round(red_pval,4);
		ful_pval=round(ful_pval,4);
	
		# R-sqrds
		red_adjsqrd=sum_red[[resp_ix]]$adj.r.squared;
		ful_adjsqrd=sum_ful[[resp_ix]]$adj.r.squared;

		red_adjsqrd=round(red_adjsqrd,3);
		ful_adjsqrd=round(ful_adjsqrd,3);

		# Plot Reduced
		cat("Plotting Reduced:\n");
		shrd=intersect(names(obs_cur), names(pred_red[,resp_ix]));

		plot(obs_cur[shrd], pred_red[shrd,resp_ix], main=response_names[resp_ix], 
			xlim=rngs, ylim=rngs,
			xlab="", ylab="Reduced Predicted");
		mtext(paste("adj. R^2=",red_adjsqrd, "  p-val=", red_pval, sep=""), line=0, cex=.7)
		abline(a=0, b=1, col="blue");

		# Plot Full
		cat("Plotting Full:\n");
		plot(obs_cur[shrd], pred_ful[shrd,resp_ix], main="", 
			xlim=rngs, ylim=rngs,
			xlab="Observed", ylab="Full Predicted (w/ Diversity)");
		mtext(paste("adj. R^2=",ful_adjsqrd, "  p-val=", ful_pval, sep=""), line=0, cex=.7)
		abline(a=0, b=1, col="blue");

		# Label page
		if((resp_ix-1)%%(plots_per_page)==0){
			mtext(title, outer=T);
		}
	}

	par(par.orig);
}

load_list=function(filename){
	arr=scan(filename, what=character(), comment.char="#");
	return(arr);
}

plot_text=function(strings, title=""){
	par.orig=par(no.readonly=T);

	par(mfrow=c(1,1));
        par(family="Courier");
        par(oma=c(.5, .5, 1.5, .5));
        par(mar=rep(0,4));

        num_lines=length(strings);

        top=max(as.integer(num_lines), 52);

        plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                );
	if(title!=""){
		mtext(title, outer=T, font=2);
	}

        text_size=max(.01, min(.8, .8 - .003*(num_lines-52)));
        #print(text_size);

        for(i in 1:num_lines){
                #cat(strings[i], "\n", sep="");
                strings[i]=gsub("\t", "", strings[i]);
                text(0, top-i, strings[i], pos=4, cex=text_size);
        }
	par(par.orig);
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

sig_char=function(val){
	if(!is.null(val) && !is.nan(val) && !is.na(val)){
		if(val <= .0001){ return("***");}
		if(val <= .001 ){ return("** ");}
		if(val <= .01  ){ return("*  ");}
		if(val <= .05  ){ return(":  ");}
		if(val <= .1   ){ return(".  ");}
	}
	return("   ");
}


##############################################################################

pdf(paste(OutputRoot, ".div_as_pred.pdf", sep=""), height=8.5, width=11);

##############################################################################
# Load matrix

counts=load_summary_file(SummaryFile);
num_categories=ncol(counts);
num_samples=nrow(counts);
#print(counts);

# Normalize
normalized=normalize(counts);
#print(normalized);

out=capture.output({
	cat("Summary Table:\n");
	cat("  Original Num Samples: ", num_samples, "\n");
	cat("  Original Num Categories: ", num_categories, "\n");
});
cat(out, sep="\n");
summary_text=c(summary_text, "\n", out);

##############################################################################
# Load factors

factors=load_factors(FactorsFile);
factor_names=colnames(factors);
num_factors=ncol(factors);
factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);

out=capture.output({
	cat("\n");
	cat("Factors:\n");
	cat("  Original Num Samples: ", num_factor_samples, "\n");
	cat("  Original Num Factors: ", num_factors, "\n");
	cat("\n");
});
cat(out, sep="\n");
summary_text=c(summary_text, "\n", out);

##############################################################################
# Relevel factor levels

if(ReferenceLevelsFile!=""){
	ref_lev_mat=load_reference_levels_file(ReferenceLevelsFile)
	factors=relevel_factors(factors, ref_lev_mat);
}else{
	cat("No Reference Levels File specified.\n");
}

##############################################################################
# Reconcile factors with samples

cat("Reconciling samples between factor file and summary table:\n");
factor_sample_ids=rownames(factors);
counts_sample_ids=rownames(counts);

shared_sample_ids=intersect(factor_sample_ids, counts_sample_ids);
num_shared_sample_ids=length(shared_sample_ids);
num_factor_sample_ids=length(factor_sample_ids);
num_counts_sample_ids=length(counts_sample_ids);

out=capture.output({
	cat("Num counts sample IDs: ", num_counts_sample_ids, "\n");
	cat("Num factor sample IDs: ", num_factor_sample_ids, "\n");
	cat("Num shared sample IDs: ", num_shared_sample_ids, "\n");
	cat("\n");
});
cat(out, sep="\n");
summary_text=c(summary_text, "\n", out);

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

##############################################################################
# Build model and select variables

if(ModelFormula!=""){
	cat("Model Formula: ", ModelFormula, "\n", sep="");
	model_string=ModelFormula;
	
}else if(CovariatesFile!=""){
	cat("Model Covariates Filename: ", CovariatesFile, "\n");
	variables=load_list(CovariatesFile);
	shared=intersect(variables, factor_names);
	if(length(shared)!=length(variables)){
		cat("Missing variables specified in Covariates File:\n");
		print(setdiff(variables, shared));
	}
	model_string=paste(shared,  collapse=" + ");
}else{
	cat("No Model Specified, using all variables.\n");
	model_string= paste(factor_names, collapse=" + ");
}

cat("Model String used for Regression: \n");
print(model_string);
model_var=get_var_from_modelstring(model_string);
cat("Predictors:\n");
print(model_var);

cat("\n");
cat("Response Variables:\n");
responses_arr=load_list(ResponseFile);
print(responses_arr);
cat("\n");

cat("Extracting predictors+responses from available factors...\n");
all_var=c(model_var, responses_arr);

factors=factors[,all_var, drop=F];
cat("\n");

##############################################################################
# Handle NAs

cat("Working on NA Removal...\n");
remove_na_res=remove_sample_or_factors_wNA_parallel(factors, required=required_arr, 
	num_trials=640000, num_cores=64, outfile=OutputRoot);
#remove_na_res=remove_sample_or_factors_wNA_parallel(factors, required=required_arr, 
#	num_trials=1000, num_cores=64, outfile=OutputRoot);

factors=remove_na_res$factors;
samp_wo_nas=rownames(factors);
factor_names=colnames(factors);
num_factors=length(factor_names);
normalized=normalized[samp_wo_nas,];
num_samples=length(samp_wo_nas);

model_string=rem_missing_var_from_modelstring(model_string, factor_names);
model_var=intersect(model_var, factor_names);

cat("\n");
cat("New Model String with Factors with NAs removed:\n");
print(model_string);
cat("\n");
responses_arr=intersect(responses_arr, factor_names);

summary_text=c(summary_text, "\n", remove_na_res$summary_text);

##############################################################################

plot_text(summary_text);

##############################################################################
# Compute diversity indices
cat("Computing diversity indices:\n");

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

# Set evenness to 0 if there is only 1 category.
evenness=div_mat[,"Evenness"];
div_mat[is.na(evenness),"Evenness"]=0;

cat("Plotting histograms of raw diversity indices.\n");
par(mfrow=c(3,2));
par(oma=c(1, 1, 1.5, 1));
par(mar=c(5,4,4,2));

for(i in 1:num_div_idx){
	hist(div_mat[, div_names[i]], main=div_names[i], xlab=div_names[i], 
		breaks=15);
}

mtext("Sample Diversity Indices Distributions", outer=T);

##############################################################################
# Output reference factor levels

predictors_mat=factors[samp_wo_nas,model_var,drop=F];
num_pred_var=ncol(predictors_mat);
responses_mat=as.matrix(factors[samp_wo_nas,responses_arr,drop=F]);
num_resp_var=ncol(responses_mat);

##############################################################################
# Plot Response Correlations
cat("Plotting correlations among responses...\n");
paint_matrix(cor(responses_mat), title="Response Correlations", plot_min=-1, plot_max=1, deci_pts=2);

##############################################################################

# Allocation matrices 
resp=rep(1, length(samp_wo_nas));
model_matrix=model.matrix(as.formula(paste("resp ~", model_string)),data=predictors_mat);
regression_variables=setdiff(colnames(model_matrix), "(Intercept)");
cat("Expected Regression Variables:\n");
print(regression_variables);
num_regression_variables=length(regression_variables);
cat("\n");

pval_list=list();
coeff_list=list();
adjrsqrd_list=list();

diversity_coef=matrix(NA, nrow=num_resp_var, ncol=num_div_idx,
		dimnames=list(responses_arr, div_names));
diversity_pval=matrix(NA, nrow=num_resp_var, ncol=num_div_idx,
		dimnames=list(responses_arr, div_names));

diversity_adj.rsqrd=matrix(NA, nrow=num_resp_var, ncol=num_div_idx,
                dimnames=list(responses_arr, div_names));
diversity_adj.rsqrd_delta=matrix(NA, nrow=num_resp_var, ncol=num_div_idx,
                dimnames=list(responses_arr, div_names));


anova_pval=matrix(NA, nrow=num_pred_var+1, ncol=num_div_idx,
		dimnames=list(c("diversity",model_var), div_names));


for(i in 1:num_div_idx){

	cat("Working on: ", div_names[i], "\n", sep="");

	diversity=div_mat[samp_wo_nas, div_names[i], drop=F];
	all_predictors_mat=cbind(diversity[samp_wo_nas,,drop=F], predictors_mat[samp_wo_nas,,drop=F]);

	# Full
	full_model_string=paste("responses_mat ~ diversity + ", model_string, sep="");
	print(full_model_string);

	mv_fit=lm(as.formula(full_model_string), data=as.data.frame(all_predictors_mat), y=T);
	sum_fit=summary(mv_fit);

	# Compute ANOVA full model
	num_resid_df=mv_fit$df.residual;
	num_responses=length(mv_fit);
	if(num_resid_df<num_responses){
		msg=paste("There are not enough residual degrees of freedom (",
			num_resid_df, ") for the number of responses (", num_responses, 
			") to perform MANOVA.", sep="");
		plot_text(msg);
		mv_anova=NULL;
	}else{

		try_res=try({mv_anova=anova(mv_fit)});
		
		if(class(try_res)=="try-error"){
			mv_anova=NULL;
			
			msg=c(
				"Error occured while computed MANOVA: ",
				paste(try_res)
			);
			plot_text(msg);

		}else{

			plot_text(capture.output(print(mv_anova)), title=div_names[i]);
			anova_pval[c("diversity",model_var),div_names[i]]=mv_anova[c("diversity",model_var),"Pr(>F)"];
		}
	}

	# Reduced
	red_model_string=paste("responses_mat ~ ", model_string, sep="");
	print(red_model_string);

	red_mv_fit=lm(as.formula(red_model_string), data=as.data.frame(all_predictors_mat), y=T);
	red_sum_fit=summary(red_mv_fit);

	# Allocate data structures
	estimates_matrix=matrix(NA, nrow=num_regression_variables, ncol=num_resp_var,
			dimnames=list(regression_variables, responses_arr));
	pvalues_matrix=matrix(NA, nrow=num_regression_variables, ncol=num_resp_var,
			dimnames=list(regression_variables, responses_arr));

	sum_resp_names=paste("Response ", responses_arr, sep="");

	if(num_resp_var==1){
		cat("Putting single variate response into list...\n");
		tmp=list();
		tmp[[sum_resp_names[1]]]=sum_fit;
		sum_fit=tmp;
		class(sum_fit)="listof";

		tmp=list();
		tmp[[sum_resp_names[1]]]=red_sum_fit;
		red_sum_fit=tmp;
		class(red_sum_fit)="listof";
	}


	for(resp_ix in 1:num_resp_var){


		missing=setdiff(regression_variables, 
			rownames(sum_fit[[sum_resp_names[resp_ix]]]$coefficients));

		avail_regression_variables=regression_variables;
		if(length(missing)>0){
			cat("***************************************************\n");
			cat("Warning: Not all regression coefficient calculable:\n");
			print(missing);
			cat("***************************************************\n");
			avail_regression_variables=setdiff(regression_variables, missing);
		}
		

		estimates_matrix[avail_regression_variables, resp_ix]=
			sum_fit[[sum_resp_names[resp_ix]]]$coefficients[avail_regression_variables, "Estimate"];
		pvalues_matrix[avail_regression_variables, resp_ix]=
			sum_fit[[sum_resp_names[resp_ix]]]$coefficients[avail_regression_variables, "Pr(>|t|)"];

		# Store diversity
		diversity_coef[resp_ix, div_names[i]]=
			sum_fit[[sum_resp_names[resp_ix]]]$coefficients["diversity", "Estimate"];
		diversity_pval[resp_ix, div_names[i]]=
			sum_fit[[sum_resp_names[resp_ix]]]$coefficients["diversity", "Pr(>|t|)"];


		# Store Adj-Rsqrd
		diversity_adj.rsqrd[resp_ix, div_names[i]]=
			sum_fit[[sum_resp_names[resp_ix]]]$adj.r.squared;

		diversity_adj.rsqrd_delta[resp_ix, div_names[i]]=
                        sum_fit[[sum_resp_names[resp_ix]]]$adj.r.squared-
                        red_sum_fit[[sum_resp_names[resp_ix]]]$adj.r.squared;
	}

	coeff_list[[div_names[i]]]=estimates_matrix;
	pval_list[[div_names[i]]]=pvalues_matrix;
	
	#paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, deci_pts=4,
	paint_matrix(estimates_matrix, title=paste("Coefficients: ", div_names[i], " Covariates", sep=""),
		plot_col_dendr=T, plot_row_dendr=T);
	paint_matrix(pvalues_matrix, title=paste("P-values: ", div_names[i], " Covariates", sep=""), 
		plot_col_dendr=T, plot_row_dendr=T,
		high_is_hot=F, plot_min=0, plot_max=1);


	# Plot predicted vs observed
	orig.par=par(no.readonly=T);
	par(mfrow=c(1,2));
	plot_pred_vs_obs(red_mv_fit, mv_fit, title=div_names[i]);
	par(orig.par);


}

#print(coeff_list);
#print(pval_list);

#print(diversity_coef);
#print(diversity_pval);

# With dendrograms
paint_matrix(diversity_coef, title="Diversity Coefficients", plot_col_dendr=T, plot_row_dendr=T);
paint_matrix(diversity_pval, title="Diversity P-values", high_is_hot=F, plot_min=0, plot_max=1,
	plot_col_dendr=T, plot_row_dendr=T);
paint_matrix(diversity_adj.rsqrd, title="Diversity Adjusted R^2 (Full Model)", 
	high_is_hot=F, plot_min=0, plot_max=1,
	plot_col_dendr=T, plot_row_dendr=T);
paint_matrix(diversity_adj.rsqrd_delta, title="Diversity Delta Adjusted R^2 (Full-Reduced Model)", 
	high_is_hot=F,
	plot_col_dendr=T, plot_row_dendr=T);

if(!is.null(mv_anova)){
	paint_matrix(anova_pval, title="Full Model MANOVAs (Diversity+Covariates) P-values", 
		high_is_hot=F, plot_min=0, plot_max=1);
}else{
	plot_text("Full Model MANOVAs (Diversity+Covariates) P-values are not available.");
}

dev.off();

##############################################################################

# Output
# MANOVA
div_ix="Tail";
sigdig=5;

fh=file(paste(OutputRoot, ".div_as_pred.model_stats.", div_ix, ".tsv", sep=""), "w");

cat(file=fh, paste("Name:", OutputRoot, "", "", sep="\t"), "\n", sep="");
cat(file=fh, paste("Diversity:", div_ix, "", "", sep="\t"), "\n", sep="");

# Output MANOVA (1 row per predictor)
anova_var=rownames(anova_pval);
cat(file=fh, "\t\t\t\n");
cat(file=fh, "Predictors:\t\tMANOVA p-values:\tsignf:\n");
for(i in 1:nrow(anova_pval)){
	cat(file=fh, paste(
		anova_var[i], 
		"",
		signif(anova_pval[i, "Tail"], sigdig), 
		sig_char(anova_pval[i, "Tail"]), 
		sep="\t"), "\n", sep="");
}
cat(file=fh, "\t\t\t\n");

cat(file=fh, "Responses:\t", div_ix, "univar coeff:\t", div_ix, "univar p-values:\tsignf:\n");
response_names=rownames(diversity_coef);
for(i in 1:nrow(diversity_coef)){ 
	cat(file=fh, paste(
		response_names[i], 
		signif(diversity_coef[i, "Tail"], sigdig), 
		signif(diversity_pval[i, "Tail"], sigdig),
		sig_char(diversity_pval[i, "Tail"]),
		sep="\t"), "\n", sep="");
}

close(fh);

##############################################################################

# Export coefficience and p-values for pred/resp analysis

write.table(t(diversity_pval), file=paste(OutputRoot, ".div_as_pred.pvals.tsv", sep=""),
        sep="\t", quote=F, col.names=NA, row.names=T);

write.table(t(diversity_coef), file=paste(OutputRoot, ".div_as_pred.coefs.tsv", sep=""),
        sep="\t", quote=F, col.names=NA, row.names=T);


##############################################################################
	
##############################################################################

cat("Done.\n");
print(warnings());
q(status=0);
