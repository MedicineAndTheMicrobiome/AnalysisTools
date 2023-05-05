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

par(mfrow=c(2,1));

hist(sample_depths, main="Distribution of Depths");
abline(v=quantiles, col=c("blue", "black", "blue"), lwd=c(1,2,1), lty=c("dashed", "solid", "dashed"));
hist(log10(sample_depths+1), main="Distribution of Log10(Depths+1)");
abline(v=log10(quantiles), col=c("blue", "black", "blue"), lwd=c(1,2,1), lty=c("dashed", "solid", "dashed"));

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

impute_zeros=function(qry_name, normalized_mat, lognormlzd_mat, max_predictors=Inf, verbose=T){
	
	target_resp=normalized_mat[,qry_name];
	
	# Split off zeros from training set
	zeros_ix=which(target_resp==0);
	num_zero_set=length(zeros_ix);
	nonzero_ix=which(target_resp!=0);
	num_nonzero_set=length(nonzero_ix);

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
	num_pred_to_select_from=min(c(num_valid_cor, max_predictors)); 

	null_model=lm(resp_values ~ 1 , data=sorted_pred_lognormlzd_df[,1:num_pred_to_select_from]);
	full_model=lm(resp_values ~ . , data=sorted_pred_lognormlzd_df[,1:num_pred_to_select_from]);

	# Perform stepwise selection
	var_sel=stepAIC(null_model, direction="both", 
		scope=list(upper=full_model, lower=null_model), trace=0);

	obs_vs_pred=cbind(resp_values, var_sel$fitted.values);
	colnames(obs_vs_pred)=c("observed", "predicted");


	sumfit=summary(var_sel);
	if(verbose){
		cat("Num variables allowed for selection: ", num_pred_to_select_from, "\n");
		print(sumfit);
	}

	# Predict the zero set
	zero_lognormlzd_df=as.data.frame(lognormlzd_mat[zeros_ix,,drop=F]);
	imputed_lognormlzd=predict(var_sel, newdata=zero_lognormlzd_df);
	imputed_normalized=10^imputed_lognormlzd;

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
	results[["name"]]=qry_name;
	results[["num_pred_to_select_from"]]=num_pred_to_select_from;
	results[["fit"]]=var_sel;
	results[["sumfit"]]=sumfit;
	results[["imputed.lognrmzd"]]=imputed_lognormlzd;
	results[["imputed.normalized"]]=imputed_normalized;
	results[["values"]]=values;
	results[["obs_vs_pred.lognrmzd"]]=obs_vs_pred;

	return(results);

}


###############################################################################
# Paralle#l execute

NumProcessors=10;
num_categories_to_process=min(num_categories, 10);

cat("Launching in parallel...\n");
all_results=mclapply(1:num_categories_to_process,
#all_results=mclapply(1:num_categories,
	function(id){
		cat(".");
		res=impute_zeros(category_names[id], 
			normalized_mat, log_normz_mat, max_predictors=50, verbose=F);
		return(res);
	},
	mc.cores=NumProcessors
	);


#for(i in 1:num_categories){
#	target_name=category_names[i];
#	cat("Working on: ", target_name, " [", i, "/", num_categories, "]\n", sep="");
#	res=impute_zeros(target_name, normalized_mat, log_normz_mat, max_predictors=80, verbose=T);
#	print(res);
#}

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

plot_diagnostics=function(imp_res, depths_arr, counts_mat, norm_mat, lognormz_mat){
	
	par(mfrow=c(3,2));
	
	#print(imp_res);
	cat_name=imp_res[["name"]];
	
	imputed_norml=imp_res[["imputed.normalized"]];
	imputed_lognrmz=imp_res[["imputed.lognrmzd"]];

	norm_vals=norm_mat[,cat_name];
	lognorm_vals=lognormz_mat[,cat_name];

	norm_range=range(c(imputed_norml, norm_vals));
	lognrm_range=range(c(imputed_lognrmz, lognorm_vals));

	num_zeros=sum(counts_mat[,cat_name]==0);
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

	plot(10^fit[,"observed"], 10^fit[,"predicted"],
		main="Training Set: Observed vs. Predicted\n(Abundance)",
		xlab="Observed (abundance)", ylab="Predicted (abundance)"
		);
	abline(0,1, col="blue");

	plot(fit[,"observed"], fit[,"predicted"],
		main="Training Set: Observed vs. Predicted\n(Log(Normalized))",
		xlab="Observed (lognormlz)", ylab="Predicted (lognormlz)"
		);
	abline(0,1, col="blue");
	mtext(text=paste("R^2: ", round(rsqrd,4), "\nAdj-R^2: ", round(rsqrd.adj, 4), sep=""),
		line=-3, cex=.75);

	# Plot depth vs counts
	#plot(depths_arr, 

	# Plot sample depth vs target depth		

	mtext(cat_name, side=3, outer=T, cex=1.5, font=2);
}


#------------------------------------------------------------------------------

par(oma=c(0,0,2,0));
# Output Diagnostic Plots
for(i in 1:num_categories_to_process){
	plot_diagnostics(all_results[[i]], depths, counts_mat, normalized_mat, log_normz_mat);
}

# Output Imputed summary table

quit();



reassess_zeros=function(target, depths, counts, normalized, lognorm, diag_plots=T){

	if(diag_plots){
		par(mfrow=c(3,2));
	}

	num_samples=nrow(counts);

	target_counts=counts[,target];
	target_prop=normalized[,target];
	target_lnorm=lognorm[,target];

	samples_w_zeros=(target_counts==0);
	num_samples_w_zeros=sum(samples_w_zeros);
	prop_samples_w_zeros=num_samples_w_zeros/num_samples;

	cat("Num Samples w/ zero: ", num_samples_w_zeros, "\n");
	cat("Prop Samples w/ zero: ", prop_samples_w_zeros, "\n");

	if(prop_samples_w_zeros>.95){
		cat("Too many zeros\n");
		return();
	}

	if(diag_plots){
		# Plot normalized
		hist(target_prop, breaks=40, 
			main=paste("Distribution of Relative Abundances:\npr(", target, "))",sep=""),
			xlab="Abundance"
			);
		title(main=paste(
			"Mean Abund: ", signif(mean(target_prop), 4), 
			" / Perc 0's: ", round(100*num_samples_w_zeros/num_samples,2), "%", sep=""),
			line=-1, font.main=1, cex.main=.9
			);


		# Plot lognorm
		hist(target_lnorm, breaks=40, 
			main=paste("Distribution of Log10(Rel. Abund.)\nw/ 0-replacement",sep=""),
			xlab="Log10(Abundance)"
			);
	}

	#----------------------------------------------------------------------
	# Predict abundance based on depth

	nonzero_set=!samples_w_zeros;
	nonzero_set_idx=which(nonzero_set);
	zero_set=samples_w_zeros;
	zero_set_idx=which(zero_set);

	nz_counts=counts[nonzero_set_idx,target];
	nz_depths=depths[nonzero_set_idx];

	nz_log10_depths=log10(nz_depths+1);
	nz_log10_counts=log10(nz_counts+1);

	# Fit
	nz_fit=lm(nz_log10_counts~nz_log10_depths);

	# Predict
	zero_set_depths=depths[zero_set_idx]
	zero_set_log10_depths=log10(zero_set_depths+1);
	zero_set_log10_counts=log10(counts[zero_set_idx, target]+1);
	df=as.data.frame(zero_set_log10_depths);
	colnames(df)="nz_log10_depths";
	log10_expected_counts=predict(nz_fit, newdata=df);

	#----------------------------------------------------------------------
	# Build training set around samples without zeros in target

	all_sample_ids=names(depths);
	training_set=all_sample_ids[!zero_set];
	impute_set=all_sample_ids[zero_set];

	selected_model=find_best_predictors(target, lognorm[training_set,]);

	zero_set_lognorm_df=as.data.frame(lognorm[zero_set_idx,,drop=F]);

	#print(selected$model);
	#print(selected$coefficients);
	#print(zero_set_normalized);

	imputed_values_lognorm=predict.lm(selected_model, newdata=zero_set_lognorm_df);
	imputed_values_norm=10^(imputed_values_lognorm);
	imputed_values_counts=imputed_values_norm*depths[zero_set];

	#cat("Imputed Normalized:\n");
	#print(imputed_values_norm)
	#cat("Imputed Counts:\n");
	#print(imputed_values_counts);

	if(diag_plots){
		# Combine Points so we can draw a new green line
		log10_combined_counts=log10(c(counts[training_set, target], imputed_values_counts));
		log10_combined_depths=log10(c(depths[training_set], depths[zero_set]));
		log10_combined_fit=lm(log10_combined_counts~log10_combined_depths);

		# Plot the regression line, and red/blue for zeros/recomputed
		plot(log10(depths+1), log10(target_counts+1), 
			ylim=range(c(log10(target_counts+1), log10_expected_counts, log10(imputed_values_counts))),
			main="Expected Target Count based on Sample Depth",
			xlab="Log10(Sample Depth)", ylab="Log10(Targets+1)");
		points(zero_set_log10_depths, zero_set_log10_counts, col="red");
		points(zero_set_log10_depths, log10_expected_counts, col="blue");
		abline(nz_fit, col="blue");
		abline(h=0, col="pink");
		abline(log10_combined_fit, col="green");
		points(zero_set_log10_depths, log10(imputed_values_counts), col="green");


		# Plot the distribution of imputed
		hist(c(imputed_values_norm, target_prop[nonzero_set]), 
			breaks=40,
			xlab="Abundance",
			main="Distribution of Relative Abundances\nw/ Imputed Samples");

		hist(c(imputed_values_lognorm, lognorm[nonzero_set, target]), 
			breaks=40,
			xlab="Log10(Abund)",
			main="Distribution of Log10(Relative Abundances)\nw/ Imputed Samples");
	}

}

###############################################################################

dev.off();

###############################################################################

cat("Done.\n")
warn=warnings();
if(length(warn)){
	print(warn);
}
q(status=0)
