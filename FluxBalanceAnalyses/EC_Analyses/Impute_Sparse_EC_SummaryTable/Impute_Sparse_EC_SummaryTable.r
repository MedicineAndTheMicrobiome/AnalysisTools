#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
library(parallel);

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 1, "character",
	"num_processors", "n", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	-o <output normalized file, imputed>\n",
	"	-n <num processors, default=1>\n",
	"\n",
	"This script will read in the summary table and then\n",
	"for each of the categories impute abundances for the\n",
	"values that had zero counts.\n",
	"\n",
	"\n", sep="");

if(
	!length(opt$input_file) || 
	!length(opt$output_file)){
	cat(usage);
	q(status=-1);
}

InputFilename=opt$input_file;
OutputFilename=opt$output_file;

NumProcessors=1;
if(!length(opt$num_processors)){
	NumProcessors=opt$num_processors;
}



###############################################################################

load_factors=function(fname){
        factors=data.frame(read.table(fname,  header=TRUE, check.names=FALSE, row.names=1, 
		comment.char="", quote="", sep="\t", stringsAsFactors=TRUE));
        dimen=dim(factors);
        cat("Rows Loaded: ", dimen[1], "\n");
        cat("Cols Loaded: ", dimen[2], "\n");
        return(factors);
}

load_summary_file=function(fname){
        cat("Loading Summary Table: ", fname, "\n");
        inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", row.names=1))
        counts_mat=inmat[,2:(ncol(inmat))];
        return(counts_mat);
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

	orig.par=par(no.readonly=T);
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

	par(orig.par);
}

log_trans=function(normalized_mat){
	# This function will take in a normalized matrix,
	# then return a log10 transformed matrix.
	# Since log10 of 0 is -Inf, the log10 of the
	# nonzero values are first transformed, and then
	# the mean is computed in transformed space.
	# The zero replacement value is then set to 1 SD lower
	# then the lowest nonzero value.

	num_cat=ncol(normalized_mat);
	num_samples=nrow(normalized_mat);

	out_mat=matrix(NA, nrow=num_samples, ncol=num_cat);
	colnames(out_mat)=colnames(normalized_mat);
	rownames(out_mat)=rownames(normalized_mat);

	for(cix in 1:num_cat){

		values=normalized_mat[,cix];

		# Calculate mean and sd of log(non-zero values)
		nonz_ix=(values>0);
		nz_values=values[nonz_ix];
		log_nz_values=log10(nz_values);

		min_log_nz=min(log_nz_values);
		sd_log_nz=sd(log_nz_values);
		zero_replace=10^(min_log_nz - sd_log_nz);

		# Replace zeros with new value
		values[!nonz_ix]=zero_replace;
		out_mat[,cix]=log10(values);
	}
	return(out_mat);
}

###############################################################################

pdf(paste(OutputFilename, ".pdf", sep=""), height=11, width=8.5);

cat("Reading: ", InputFilename, "\n");
counts_mat=load_summary_file(InputFilename);

# Remove zero count categories
cat_counts=apply(counts_mat, 2, sum);
zero_counts_ix=(cat_counts==0);
if(any(zero_counts_ix)){
	cat("Zero count categories detected: ", sum(zero_counts_ix), "\n");
	counts_mat=counts_mat[, !zero_counts_ix];
}

num_categories=ncol(counts_mat);
num_samples=nrow(counts_mat);

category_names=colnames(counts_mat);

cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");

sample_depths=apply(counts_mat, 1, sum);
quantiles=quantile(sample_depths, c(.025, .975));

par(mfrow=c(4,1));

hist(sample_depths, main="Distribution of Depths", xlab="Sample Depth", breaks=40);
abline(v=quantiles, col=c("blue", "black", "blue"), lwd=c(1,2,1), lty=c("dashed", "solid", "dashed"));

hist(log10(sample_depths+1), main="Distribution of Log10(Depths+1)", xlab="Log10(Sample Depth + 1)", breaks=40);
abline(v=log10(quantiles), col=c("blue", "black", "blue"), lwd=c(1,2,1), lty=c("dashed", "solid", "dashed"));

zero_counts=apply(counts_mat, 1, function(x){ sum(x==0)});
hist(zero_counts, main="Distribution of Num Zeros (acrossed samples) in Categories", xlab="Num Zeros in Category",
	 breaks=40);
hist(log10(zero_counts+1), main="Distribution of Log10(Num Zero in Categories + 1) across Samples",
	xlab="Log10(Num Zero + 1)", breaks=40);

hist(zero_counts/num_categories, main="Distribution of Proportion of Zeros Categories across Samples", 
	xlab="Proportion of Zeros in Category",
	xlim=c(0,1), breaks=40);

cat("Normalizing Counts...\n");
normalized_mat=normalize(counts_mat);

cat("Log10 Transforming Normalized Data...\n");
log_normz_mat=log_trans(normalized_mat);

if(0){
	par(mfrow=c(5,2));
	for(i in 1:ncol(log_norm_mat)){
		cat("Plotting: ", i, "/", category_names[i], "\n");
		hist(normalized_mat[,i], breaks=40);
		hist(log_norm_mat[,i], breaks=40);
	}
}

###############################################################################

simplify_model=function(model, resp, data, verbose=T, max_keep=10, min_keep=3, adj_pval_cutoff=0.01){

	coefficients=summary(model)$coefficients;
	no_intercept=setdiff(rownames(coefficients), "(Intercept)");

	if(length(no_intercept)==0){
		return(model);
	}

	coefficients=coefficients[no_intercept,,drop=F];
	pvalue_sort_ix=order(coefficients[,"Pr(>|t|)"]);
	sorted_coefficients=coefficients[pvalue_sort_ix,,drop=F];
	# print(sorted_coefficients);

	# Keep predictions within between min and max keep, depending on FDR pvalue
	fdr_adj_pval=p.adjust(sorted_coefficients[,"Pr(>|t|)"]);
	gt_pval_ix=fdr_adj_pval <= adj_pval_cutoff;
	keep_ix=gt_pval_ix;

	if(sum(keep_ix)<min_keep){
		keep_ix[1 : min(min_keep, length(keep_ix))]=T;
	}

	if(sum(keep_ix) >= max_keep){
		keep_ix[(max_keep+1) : (length(keep_ix))]=F;
	}

	
	# Rebuild model with only the kept predictors
	keep_predictors=rownames(sorted_coefficients[keep_ix,,drop=F]);
	new_model_str=paste("resp ~ ", paste(keep_predictors, collapse=" + "));
	new_model=lm(as.formula(new_model_str), data=data);	

	if(verbose){
		cat("-------------------------------------------------------\n");
		print(keep_ix);
		cat("-------------------------------------------------------\n");
		cat("Variable Selection Model:\n");
		print(summary(model)$coefficients);
		cat("\n");
		cat("Simplified Model:\n");
		print(summary(new_model)$coefficients);
		cat("\n");
		cat("-------------------------------------------------------\n");
		cat("\n");
	}
	
	return(new_model);

}

impute_zeros=function(id, qry_name, normalized_mat, lognormlzd_mat, max_predictors=Inf, 
	do_simplify_model=T, verbose=T){
	
	target_resp=normalized_mat[,qry_name];
	
	# Split off zeros from training set
	zeros_ix=which(target_resp==0);
	num_zero_set=length(zeros_ix);
	prop_zero=num_zero_set/length(target_resp);
	nonzero_ix=which(target_resp!=0);
	num_nonzero_set=length(nonzero_ix);

	if(num_nonzero_set<10){
		results=list();
		results[["id"]]=id;
		results[["name"]]=qry_name;
		results[["success"]]=F;
		results[["num_zeros"]]=num_zero_set;
		results[["prop_zeros"]]=prop_zero;
		results[["num_pred_to_select_from"]]=NA;
		results[["num_variables_selected"]]=NA;
		results[["fit"]]=NA;
		results[["sumfit"]]=NA;
		results[["imputed.lognrmzd"]]=1;
		results[["imputed.normalized"]]=0;
		results[["values"]]=NA;
		results[["obs_vs_pred.lognrmzd"]]=0;
		return(results);
	}

	# Build Training set and remove qry from predictors
	resp_values=lognormlzd_mat[nonzero_ix, qry_name];
	predictors_names=setdiff(colnames(lognormlzd_mat), qry_name);
	predictors_lognormlzd_mat=lognormlzd_mat[nonzero_ix, predictors_names];

	if(verbose){
		cat("Target: ", qry_name, "\n");
		cat("Num Zeros: ", num_zero_set, "\n");
		cat("Num NonZeros: ", num_nonzero_set, "\n");
	}
	
	# Compute correlations to sort predictors:
	correls=cor(resp_values, predictors_lognormlzd_mat);
	names(correls)=predictors_names;	
	correl_magn_sorted_ix=order(abs(correls), decreasing=T);
	correls_sorted=correls[correl_magn_sorted_ix];

	# Remove invalid correls
	valid_cor=which(!is.na(correls_sorted));
	correls_sorted=correls_sorted[valid_cor];
	num_valid_cor=length(correls_sorted);

	if(verbose){
		cat("Top correlated: \n");
		print(correls_sorted[1:10]);
	}

	# Sort data.frame
	sorted_pred_lognormlzd_df=as.data.frame(predictors_lognormlzd_mat[,names(correls_sorted)]);
	
	# Define stepwise variables it can see
	num_pred_to_select_from=min(c(num_valid_cor, num_nonzero_set-2, max_predictors)); 

	null_model=lm(resp_values ~ 1 , data=sorted_pred_lognormlzd_df[,1:num_pred_to_select_from]);
	full_model=lm(resp_values ~ . , data=sorted_pred_lognormlzd_df[,1:num_pred_to_select_from]);

	# Perform stepwise selection
	var_sel=stepAIC(null_model, direction="both", 
		scope=list(upper=full_model, lower=null_model), trace=ifelse(verbose, 1, 0));
	
	num_var_stepwise_selected=nrow(summary(var_sel)$coefficients)-1;

	if(do_simplify_model){
		var_sel=simplify_model(var_sel, resp_values, 
			data=sorted_pred_lognormlzd_df[,1:num_pred_to_select_from],
			verbose=verbose);
	}

	obs_vs_pred=cbind(resp_values, var_sel$fitted.values);
	colnames(obs_vs_pred)=c("observed", "predicted");

	sumfit=summary(var_sel);
	num_variables_selected=nrow(sumfit$coefficients)-1;
	if(verbose){
		cat("Num variables allowed for selection: ", num_pred_to_select_from, "\n");
		cat("Number of Variables Selected: ", num_variables_selected, "\n");
		print(sumfit);
	}



	# Predict the zero set
	zero_lognormlzd_df=as.data.frame(lognormlzd_mat[zeros_ix,,drop=F]);
	imputed_lognormlzd=predict(var_sel, newdata=zero_lognormlzd_df);
	imputed_normalized=10^imputed_lognormlzd;

	imputed_normalized[imputed_normalized==0]=1e-323;

	if(verbose){
		cat("Imputed Log10(Normalized):\n");
		print(imputed_lognormlzd);
		cat("\nImputed Normalized:\n");
		print(imputed_normalized);
	}

	imputed_sample_ids=rownames(zero_lognormlzd_df);
	combined=target_resp;
	combined[imputed_sample_ids]=imputed_normalized[imputed_sample_ids];
	values=cbind(target_resp, combined);
	colnames(values)=c("Original", "Imputed");


	# Package results to return
	results=list();
	results[["id"]]=id;
	results[["name"]]=qry_name;
	results[["success"]]=T;
	results[["num_zeros"]]=num_zero_set;
	results[["prop_zeros"]]=prop_zero;
	results[["num_pred_to_select_from"]]=num_pred_to_select_from;
	results[["num_variables_stepwise_selected"]]=num_var_stepwise_selected;
	results[["num_variables_selected"]]=num_variables_selected;
	results[["fit"]]=var_sel;
	results[["sumfit"]]=sumfit;
	results[["imputed.lognrmzd"]]=imputed_lognormlzd;
	results[["imputed.normalized"]]=imputed_normalized;
	results[["values"]]=values;
	results[["obs_vs_pred.lognrmzd"]]=obs_vs_pred;

	return(results);

}


#cross_validate=function(id, qry_name, normalized_mat, lognormlzd_mat){
	# Remove zeros
	# if number non-0 remaining > 40
	#impute_zeros(id, qry_name, normalized_mat, lognormlzd_mat, max_predictors=Inf, verbose=T);
#}


###############################################################################
# Paralle#l execute

NumProcessors=180;
num_categories_to_process=min(num_categories, Inf);
max_allowed_predictors=min(Inf);

cat("Launching in parallel...\n");
process_list=1:num_categories_to_process;
#process_list=c(46);
#process_list=1:100;


go_verbose=length(process_list)==1;

all_results=mclapply(process_list,
	function(id){
		cat(".");
		tryCatch({
			res=impute_zeros(id, category_names[id], 
				normalized_mat, log_normz_mat, max_predictors=max_allowed_predictors, 
				verbose=go_verbose);
		}, error=function(e){
                        cat("Error occurred on ", id, "th column during impute_zeros.\n");
                        print(e);
                });
		return(res);
	},
	mc.cores=NumProcessors
	);

###############################################################################

add_samples_to_hist=function(hist_res, samples){
	hist_out=hist(samples, breaks=hist_res$breaks, plot=F);
	num_counts=length(hist_out$counts);
	for(i in 1:num_counts){
		points(
			rep(hist_out$mids[i], 2),
			c(0, hist_out$counts[i]),
			type="l", col="green4",
			lwd=4, lend="butt"
		);
	}
	
	par=par();
	top=par$usr[3]+(par$usr[4]-par$usr[3])*.85;
	right=par$usr[1]+(par$usr[2]-par$usr[1])*.75;
	legend(right, top, fill="green4", legend="Imputed");

}

plot_diagnostics=function(imp_res, sample_depths, counts_mat, norm_mat, lognormz_mat){

	par(mfrow=c(3,2));
	cat_name=imp_res[["name"]];

	mar_title=paste("[", imp_res[["id"]], "]: ", cat_name, sep="");

	if(!imp_res[["success"]]){
		plot(0,0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",
			xaxt="n", yaxt="n", bty="n");

		text(.5, .5, 
			paste("Too many zeros to impute.\n", 
				"Num Zeros: ", imp_res[["num_zeros"]], "\n",
				"Prop Zeros: ", round(imp_res[["prop_zeros"]], 4)),	
			col="red", cex=2);

		hist(log10(counts_mat[,cat_name]+1),
			xlab="Counts",
			main="Log10(Target Counts)");
		hist(norm_mat[,cat_name],
			xlab="Abundance",
			main="Abundance"
			);
		hist(lognormz_mat[,cat_name],
			xlab="Log10(Normalized)",
			main="Log10(Normalized)"
			);
		mtext(mar_title, side=3, outer=T, cex=1.5, font=2);
		return();	
	}
	
	
	imputed_norml=imp_res[["imputed.normalized"]];
	imputed_lognrmz=imp_res[["imputed.lognrmzd"]];

	norm_vals=norm_mat[,cat_name];
	lognorm_vals=lognormz_mat[,cat_name];

	norm_range=range(c(imputed_norml, norm_vals));
	lognrm_range=range(c(imputed_lognrmz, lognorm_vals));

	zero_ix=(counts_mat[,cat_name]==0);
	num_zeros=sum(zero_ix);
	prop_zeros=num_zeros/nrow(counts_mat);

	# Plot relative abundance
	ra_hst=hist(norm_mat[,cat_name], breaks=seq(norm_range[1], norm_range[2], length.out=40),
		main="Distribution of\nOriginal Relative Abundances", 
		xlab="Relative Abundance"
		);
	add_samples_to_hist(ra_hst, imputed_norml);
	mtext(text=paste("Num Zeros: ", num_zeros, "\nProp Zeros: ", round(prop_zeros, 4), sep=""),
		line=-3, cex=.75);

	# Plot log_normlz 
	ln_hst=hist(lognormz_mat[,cat_name], breaks=seq(lognrm_range[1], lognrm_range[2], length.out=40),
		main="Distribution of\nOriginal Log10(Relative Abundances)",
		xlab="Log10(Relative Abundance)"
		);
	add_samples_to_hist(ln_hst, imputed_lognrmz);

	# Plot obs vs predict
	fit=imp_res[["obs_vs_pred.lognrmzd"]]
	sumfit=imp_res[["sumfit"]];
	rsqrd=sumfit$r.squared;
	rsqrd.adj=sumfit$adj.r.squared;

	fit_range=range(10^c(fit[,"observed"], fit[,"predicted"]));
	plot(10^fit[,"observed"], 10^fit[,"predicted"],
		main="Training Set: Observed vs. Predicted\n(Abundance)",
		xlab="Observed (abundance)", ylab="Predicted (abundance)",
		xlim=fit_range, ylim=fit_range
		);
	abline(0,1, col="blue");

	fit_range=range(c(fit[,"observed"], fit[,"predicted"]));
	plot(fit[,"observed"], fit[,"predicted"],
		main="Training Set: Observed vs. Predicted\n(Log(Normalized))",
		xlab="Observed (lognormlz)", ylab="Predicted (lognormlz)",
		xlim=fit_range, ylim=fit_range
		);
	abline(0,1, col="blue");
	mtext(text=paste(
		" R^2: ", round(rsqrd,4), 
		"\n Adj-R^2: ", round(rsqrd.adj, 4), 
		"\n Max Pred Allowed: ", imp_res[["num_pred_to_select_from"]],
		"\n Num Pred Stepwise Selected: ", imp_res[["num_variables_stepwise_selected"]],
		"\n Num Pred Finally Selected: ", imp_res[["num_variables_selected"]],
		sep=""),
		line=-4, adj=0, cex=.55, family="mono", col="blue");

	# Plot depth vs counts
	log10p1_counts=log10(counts_mat[,cat_name]+1);
	log10p1_depths=log10(sample_depths);

	imputed_norm=imp_res[["imputed.normalized"]];
	imputed_names=names(imputed_norm);
	log10p1_imputed_sample_depths=log10(sample_depths[imputed_names]);
	log10p1_imputed_counts=log10(imputed_norm*sample_depths[imputed_names]);

	combined_logdepth_range=range(c(log10p1_depths, log10p1_imputed_sample_depths));
	combined_logcount_range=range(c(log10p1_counts, log10p1_imputed_counts));

	combined_logcounts=log10(imp_res[["values"]][,"Imputed"]*sample_depths);
	imputed_and_observed_fit=lm(combined_logcounts~log10p1_depths);

	plot(log10p1_depths, log10p1_counts,
		xlim=combined_logdepth_range,
		ylim=combined_logcount_range,
		xlab="Log10(Sample Depths)",
		ylab="Log10(Target Counts)",
		main="Target Counts vs Sample Depths"
	);
	abline(h=0, col="pink");
	points(log10p1_depths[zero_ix], rep(0, num_zeros), col="red");

	abline(imputed_and_observed_fit, col="blue");
	points(log10p1_imputed_sample_depths, log10p1_imputed_counts, col="green");

	# Write coefficients and selected predictors to screen
	old_par=par(no.readonly=T);
	par(mar=c(0,0,0,0));
	plot(0,0, type="n", xlab="", ylab="", main="", bty="n", xaxt="n", yaxt="n",
		xlim=c(0,1), ylim=c(0,1));
	char_height=par()$cxy[2];
	coef_table=imp_res[["sumfit"]]$coefficients[-1,c("Estimate", "Pr(>|t|)"), drop=F];
	coef_order_ix=order(coef_table[,"Pr(>|t|)"]);
	coef_text=capture.output({print(coef_table[coef_order_ix,,drop=F])});
	num_lines=length(coef_text);
	text_scale=min(1, num_lines/char_height);

	if(text_scale==1){
		max_lines=1/char_height;
		coef_text=c(coef_text, rep("", max(0, max_lines-num_lines)));
		num_lines=length(coef_text);
	}

	for(i in num_lines:1){
		text(0,1-(i/num_lines*text_scale), adj=0, labels=coef_text[i], cex=text_scale*.95, family="mono");
	}
	par(old_par);

	# Plot sample depth vs target depth		
	mtext(mar_title, side=3, outer=T, cex=1.5, font=2);
}


#------------------------------------------------------------------------------

cat("Imputations completed.\n")

par(oma=c(0,0,2,0));
# Output Diagnostic Plots
cat("Length of results: ", length(all_results), "\n");
for(i in 1:length(all_results)){
	tryCatch({
		plot_diagnostics(all_results[[i]], sample_depths, counts_mat, normalized_mat, log_normz_mat);
		}, error=function(e){
			cat("Error occurred on ", i, "th column during plot_diagnostics.\n");
			print(e);
		}
	);
}

###############################################################################

# Output Imputed Abundances
cat("Collecting Imputed Abundances into Output Matrix...\n");
out_imputed_norm=normalized_mat;
nsamples=nrow(normalized_mat);
ncategories=ncol(normalized_mat);
for(i in 1:length(all_results)){

	tryCatch({

			if(is.null(all_results[[i]])){
				cat("Unexpected failure for column: ", i, 
					", (null) exporting original values.\n", sep="");
				next;
			}

			if(all_results[[i]][["success"]]){
				out_imputed_norm[,i]=all_results[[i]][["values"]][,"Imputed"];
			}else{
				cat("Failed imputation for column: ", i, 
					", (Too many zeros) exporting original values.\n", sep="");
			}

		}, error=function(e){
			cat("Error occured on ", i, "th column during export to output matrix.\n");
			print(e);
		}
		);
}

#------------------------------------------------------------------------------

write_summary_file=function(out_mat, fname){
        fc=file(fname, "w");
        cat(file=fc, paste("sample_id\ttotal", paste(colnames(out_mat), collapse="\t"), sep="\t"));
        cat(file=fc, "\n");
        sample_names=rownames(out_mat);
        num_samples=nrow(out_mat);
        for(samp_idx in 1:num_samples){
                total=sum(out_mat[samp_idx,]);
                outline=paste(sample_names[samp_idx], total,
                        paste(out_mat[samp_idx,], collapse="\t"), sep="\t");
                cat(file=fc, outline);
                cat(file=fc, "\n");
        }
        close(fc);
}

cat("Writing Imputed Abundances to Summary Table...\n");
write_summary_file(out_imputed_norm, OutputFilename);

###############################################################################

dev.off();

###############################################################################

cat("Done.\n")
warn=warnings();
if(length(warn)){
	print(warn);
}
q(status=0)
