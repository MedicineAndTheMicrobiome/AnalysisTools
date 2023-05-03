#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"groupings_file", "g", 1, "character",
	"output_file", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	-g <groupings file (pathways, etc.)>\n",
	"	-o <output file, imputed>\n",
	"\n",
	"This script will read in the summary table and groupings and\n",
	"then try to 'recompute' the 0's based on the using the other\n",
	"reactions from the same pathway.\n",
	"\n",
	"\n", sep="");

if(
	!length(opt$input_file) || 
	!length(opt$groupings_file) || 
	!length(opt$output_file)){
	cat(usage);
	q(status=-1);
}

InputFilename=opt$input_file;
GroupingsFilename=opt$groupings_file;
OutputFilename=opt$output_file;



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
quantiles=quantile(sample_depths, c(.025, .5, .975));
log10_quantiles=log10(quantiles);

par(mfrow=c(2,1));

hist(sample_depths, main="Distribution of Depths");
abline(v=quantiles, col=c("blue", "black", "blue"), lwd=c(1,2,1), lty=c("dashed", "solid", "dashed"));
hist(log10(sample_depths+1), main="Distribution of Log10(Depths+1)");
abline(v=log10_quantiles, col=c("blue", "black", "blue"), lwd=c(1,2,1), lty=c("dashed", "solid", "dashed"));

normalized_mat=normalize(counts_mat);
log_norm_mat=log_trans(normalized_mat);

if(0){
	par(mfrow=c(5,2));
	for(i in 1:ncol(log_norm_mat)){
		cat("Plotting: ", i, "/", category_names[i], "\n");
		hist(normalized_mat[,i], breaks=40);
		hist(log_norm_mat[,i], breaks=40);
	}
}


num_samples_below_cutoff=sum(sample_depths<quantiles[1]);
cat("Number of Samples with depth less than 2.5% percentile: ", num_samples_below_cutoff, "\n");

###############################################################################

find_best_predictors=function(qry_name, ln_db_vals, testn=20){
	
	# Function should only receive the rows it is allowed to look at.
	# ln_db_vals: log normal database values

	resp=ln_db_vals[,qry_name];	
	num_samples_avail=nrow(ln_db_vals);
	
	# Remove qry from db_vals, then compute correlations so
	# we can order by them in decreasing magnitude
	predictors_noqry=setdiff(colnames(ln_db_vals), qry_name);
	correls=cor(resp, ln_db_vals[,predictors_noqry]);
	names(correls)=predictors_noqry;
	correls_mag_sort_ix=order(abs(correls), decreasing=T);
	correls_sorted=correls[correls_mag_sort_ix];

	# Remove cor that are 0, (i.e. giving NAs)
	valid_cor=which(!is.na(correls_sorted));
	correls_sorted=correls_sorted[valid_cor];
	num_valid_cor=length(correls_sorted);
	#print(correls_sorted);

	db_vals_cormag_sorted=as.data.frame(ln_db_vals[,names(correls_sorted)]);
	num_max_var=min(c(num_samples_avail-1, num_valid_cor, testn));
	cat("Num Available Samples: ", num_samples_avail, "\n");
	cat("Num Valid Correlations: ", num_valid_cor, "\n");
	cat("Number of Predictors to Include: ", num_max_var, "\n");

	# Define stepwise variables it can see
	null_model=lm(resp ~ 1 , data=db_vals_cormag_sorted[,1:num_max_var]);
	full_model=lm(resp ~ . , data=db_vals_cormag_sorted[,1:num_max_var]);

	var_sel=stepAIC(null_model, direction="both", 
		scope=list(upper=full_model, lower=null_model), trace=0);

	sumfit=summary(var_sel);
	print(sumfit);

	plot(resp, var_sel$fitted.values, 
		main="Predict vs. Observed LogNorm\nfor Counts > 0",
		xlab="LogNorm(Observed)", ylab="LogNorm(Predicted)");
	abline(a=0, b=1, col="blue", lty="dashed");
	title(main=paste(
		"R-squared: ", round(sumfit$r.squared, 4),
		"\n",
		"Adj-R-squared: ", round(sumfit$adj.r.squared, 4),
		sep=""),
		line=-2, font.main=1, cex.main=.9
		);
	
	return(var_sel);	

}

reassess_zeros=function(target, depths, counts, normalized, lognorm){

	par(mfrow=c(3,2));

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

	# Plot depth vs counts
	#plot(log10(depths+1), log10(target_counts+1), 
	#	main=paste("Log10(", target, " Counts) vs. Log10(Sample Depth)", sep=""), 
	#	xlab="Log10(Sample Depth + 1)", ylab="Log10(Targets + 1)");
	#abline(h=0, col="pink");

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

	# Plot the regression line, and red/blue for zeros/recomputed
	#plot(log10(depths+1), log10(target_counts+1), 
	#	ylim=range(c(log10(target_counts+1), log10_expected_counts)),
	#	main="Expected Target Count based on Sample Depth",
	#	xlab="Log10(Sample Depth + 1)", ylab="Log10(Targets + 1)");
	#points(zero_set_log10_depths, zero_set_log10_counts, col="red");
	#points(zero_set_log10_depths, log10_expected_counts, col="blue");
	#abline(nz_fit, col="blue");
	#abline(h=0, col="pink");

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

	# Plot the regression line, and red/blue for zeros/recomputed
	plot(log10(depths+1), log10(target_counts+1), 
		ylim=range(c(log10(target_counts+1), log10_expected_counts, log10(imputed_values_counts))),
		main="Expected Target Count based on Sample Depth",
		xlab="Log10(Sample Depth)", ylab="Log10(Targets+1)");
	points(zero_set_log10_depths, zero_set_log10_counts, col="red");
	points(zero_set_log10_depths, log10_expected_counts, col="blue");
	abline(nz_fit, col="blue");
	abline(h=0, col="pink");
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

	#quit();

}

###############################################################################

#correlation_mat=cor(impute_trainingset_mat);


#for(i in which(category_names=="1.1.1.125")){
for(i in 1:num_categories){

	# select predictors
	target_name=category_names[i];
	cat("Working on: ", target_name, " [", i, "/", num_categories, "]\n", sep="");

	reassess_zeros(target_name, sample_depths, counts_mat, normalized_mat, log_norm_mat);

	cat("\n");
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
