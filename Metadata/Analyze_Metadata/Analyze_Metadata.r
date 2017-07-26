#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
library(plotrix);

options(useFancyQuotes=F);
options(digits=5)

params=c(
	"factors", "f", 1, "character",
	"outputroot", "o", 2, "character",
	"exclude", "x", 2, "character",
	"include", "n", 2, "character",
	"response", "y", 2, "character"
);

DEF_MIN_NONNA_PROP=.65;

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <factors>\n",
	"\n",
	"	User specify which subset of variables to LASSO:\n",
	"	[--include <subset variables list filename>]\n",
	"	[--exclude <file with variables to exclude>]\n",
	"\n",
	"	NA Handling:\n",
	"	[-p <min nonNA proportion, default=", DEF_MIN_NONNA_PROP, ">]\n",
	"	Note: if you allow NAs, they will be replaced with average across nonNA\n",
	"\n",
	"This script will read in a metadata file and generate some\n",
	"analyses on the samples underlying them.\n",
	"\n");

if(!length(opt$factors)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputFnameRoot=gsub(".tsv", "", opt$factors);
}else{
	OutputFnameRoot=opt$outputroot;
}
OutputFnameRoot=paste(OutputFnameRoot, ".fact_info", sep="");

if(!length(opt$include)){
	VariableIncludeListFname="";
}else{
	VariableIncludeListFname=opt$include;
}

if(!length(opt$exclude)){
	VariableExcludeListFname="";
}else{
	VariableExcludeListFname=opt$exclude;
}

if(!length(opt$min_nonna_pro)){
	MinNonNAProp=DEF_MIN_NONNA_PROP;
}else{
	MinNonNAProp=opt$min_nonna_pro;
}

FactorsFname=opt$factors;

cat("Factors Filename: ", FactorsFname, "\n", sep="");
cat("Output Filename Root: ", OutputFnameRoot, "\n", sep="");
cat("Minimum Non-NA Proportion: ", MinNonNAProp, "\n", sep="");
cat("\n");

if(VariableIncludeListFname!=""){
	cat("Using subset of variables from: ", VariableIncludeListFname, " (Inclusion List)\n");
}

if(VariableExcludeListFname!=""){
	cat("Excluding subset of variables from: ", VariableExcludeListFname, " (Exclusion List)\n");
}

###############################################################################

load_factors=function(fname){
	factors=as.data.frame(read.delim(fname,  header=TRUE, row.names=1, check.names=FALSE, sep="\t"));
	return(factors);
}

load_factors_as_text=function(fname){
	factors=read.delim(fname,  header=TRUE, row.names=1, check.names=FALSE, sep="\t", colClasses="character");
	return(factors);
}

##############################################################################

plot_text=function(strings, max_lines=75){

	plot_page=function(strings){
		orig_par=par(no.readonly=T);

		par(mfrow=c(1,1));
		par(family="Courier");
		par(oma=rep(.5,4));
		par(mar=rep(0,4));

		num_lines=length(strings);

		top=max_lines;

		plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
			xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
			);
		for(i in 1:num_lines){
			#cat(strings[i], "\n", sep="");
			text(0, top-i, strings[i], pos=4, cex=.8);
		}

		par(orig_par);
	}

	num_lines=length(strings);
	num_pages=ceiling(num_lines / max_lines);
	#cat("Num Pages: ", num_pages, "\n");
	for(page_ix in 1:num_pages){
		start=(page_ix-1)*max_lines+1;
		end=start+max_lines-1;
		end=min(end, num_lines);
		##print(c(start,end));
		plot_page(strings[start:end]);
	}
}

##############################################################################

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, counts=F){

	orig_par=par(no.readonly=T);
	par(mfrow=c(1,1));

        num_row=nrow(mat);
        num_col=ncol(mat);

        cat("Num Rows: ", num_row, "\n");
        cat("Num Cols: ", num_col, "\n");

        mat=mat[rev(1:num_row),, drop=F];

        num_colors=50;
        color_arr=rainbow(num_colors, start=0, end=4/6);
        if(high_is_hot){
                color_arr=rev(color_arr);
        }

        remap=function(in_val, in_range, out_range){
                in_prop=(in_val-in_range[1])/(in_range[2]-in_range[1])
                out_val=in_prop*(out_range[2]-out_range[1])+out_range[1];
                return(out_val);
        }

        if(is.na(plot_min)){
                plot_min=min(mat);
        }
        if(is.na(plot_max)){
                plot_max=max(mat);
        }
        cat("Plot min/max: ", plot_min, "/", plot_max, "\n");
	par(mar=c(10,10,1,1));
        plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), xaxt="n", yaxt="n", bty="n", xlab="", ylab="", main=title);

        # x-axis
        axis(side=1, at=seq(.5, num_col-.5, 1), labels=colnames(mat), las=2, cex.axis=.7);
        axis(side=2, at=seq(.5, num_row-.5, 1), labels=rownames(mat), las=2, cex.axis=.7);

        if(log_col){
                plot_min=log10(plot_min+.0125);
                plot_max=log10(plot_max+.0125);
        }

	text_size=min(1, 17.5/num_row);
	cat("Text Size: ", text_size, "\n");

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

                        if(counts){
                                text_lab=sprintf("%i", mat[y,x]);
                        }else{
                                text_lab=mat[y,x];
                        }
                        text(x-.5, y-.5, text_lab, srt=45, cex=text_size, font=2);
                }
        }

	par(orig_par);

}

##############################################################################

process_factor_NAs=function(factors, min_non_NA_prop=.95){
	num_factors=ncol(factors);
	num_samples=nrow(factors);

	max_accetable_na=1-min_non_NA_prop;
	keep=rep(F, num_factors);
	factor_names=colnames(factors);
	cat("Processing factors for NAs:\n");
	cat("  Max percentage of NAs allowed: ", round(max_accetable_na*100, 1) , "%\n\n", sep="");
	for(fact_ix in 1:num_factors){
		vals=factors[,fact_ix];
		na_ix=is.na(vals);	
		numNAs=sum(na_ix);
		propNA=numNAs/num_samples;
		if(propNA<max_accetable_na){
			keep[fact_ix]=TRUE;
			if(numNAs>0){
				cat(factor_names[fact_ix], ": \n     ", sep="");
				non_nas_val=vals[!na_ix];
				num_unique=length(unique(non_nas_val));
		
				if(is.numeric(vals) && num_unique>2){
					# If values appear to be continuous
					med_val=median(non_nas_val);
					factors[na_ix, fact_ix]=med_val;
					cat(numNAs, " NA values replaced with: ", med_val, "\n", sep="");
				}else{
					# If values appear to be categorical
					resampled=sample(non_nas_val, numNAs, replace=TRUE);
					factors[na_ix, fact_ix]=resampled;
					cat(numNAs, " NA values replaced from: ", paste(head(resampled), collapse=", "), "\n", sep="");
				}
			}
		}else{
			cat(factor_names[fact_ix], ":\n     Removed for too many NAs (", round(propNA*100,1) , "% NA)\n", sep="");
		}

		# Check to see if filling in NAs lead to single value
		num_unique=length(unique(factors[, fact_ix]));
		if(num_unique==1){
			keep[fact_ix]=F;	
			cat(factor_names[fact_ix], ":\n     Removed for having no variance (identical values).\n", sep="");
		}
	}

	num_kept=sum(keep);
	cat("\n");
	cat(num_kept, "/", num_factors, " Factors Kept.\n\n", sep="");

	return(factors[, keep]);
}

##############################################################################

compute_correl=function(factors, pval_cutoff=0.05, abs_correl_cutoff=0.5){
# This function will calculate the correlation and pvalue between all
# factors and then recommend interaction terms 

	num_factors=ncol(factors);
	factor_names=colnames(factors);
	pvalue_mat=matrix(NA, nrow=num_factors, ncol=num_factors, dimnames=list(factor_names,factor_names));
	correl_mat=matrix(NA, nrow=num_factors, ncol=num_factors, dimnames=list(factor_names,factor_names));
	corrltd_mat=matrix("", nrow=num_factors, ncol=num_factors, dimnames=list(factor_names,factor_names));

	for(i in 1:num_factors){

		for(j in 1:num_factors){
			if(i<j){

				# Remove NAs and then compute correl 
				not_na_ix=!(is.na(factors[,i]) | is.na(factors[,j]));

				if(sum(not_na_ix)>2){
					test_res=cor.test(factors[not_na_ix,i], factors[not_na_ix,j]);
					correl_mat[i,j]=test_res$estimate;
					pvalue_mat[i,j]=test_res$p.value;
					is_cor=abs(test_res$estimate)>abs_correl_cutoff && test_res$p.value<=pval_cutoff;
					corrltd_mat[i,j]=ifelse(is_cor, "X", ".");
				}else{
					correl_mat[i,j]=0;
					pvalue_mat[i,j]=1;
					corrltd_mat[i,j]=".";
				}

				# Copy over symmetric values
				correl_mat[j,i]=correl_mat[i,j];
				pvalue_mat[j,i]=pvalue_mat[i,j];
				corrltd_mat[j,i]=corrltd_mat[i,j];
			}
		}
	}	

	results=list();
	results$pval_cutoff=pval_cutoff;
	results$correl_cutoff=abs_correl_cutoff;
	results$corrltd_mat=corrltd_mat;
	results$correl_mat=correl_mat;
	results$pvalue_mat=pvalue_mat;

	return(results);

}

##############################################################################

recode_non_numeric_factors=function(factors, var_dependencies){
	num_factors=ncol(factors);
	num_samples=nrow(factors);
	fact_names=colnames(factors);
		
	out_factors=matrix(nrow=num_samples, ncol=0);
	for(i in 1:num_factors){

		fact_val=factors[,i];
		nonNAix=!is.na(fact_val);

		if(!is.numeric(fact_val)){
			
			cat("\n", fact_names[i], ": Recoding.\n", sep="");

			unique_types=sort(unique(fact_val[nonNAix]));
			num_unique=length(unique_types);

			if(num_unique==1){
				cat("No information, all available types are the same.\n");
				next;
			}else if(num_unique==0){
				cat("No information, no available types that are not NA.\n");
				next;
			}

			cat("Unique types:\n");
			print(unique_types);

			recoded_mat=numeric();
			# Put reference and type in new variable name
			for(type_ix in 2:num_unique){
				recoded_mat=cbind(recoded_mat, unique_types[type_ix]==fact_val);	
			}
			new_label=paste(fact_names[i], "_ref", unique_types[1], "_is", unique_types[2:num_unique], sep="");
			cat("Recoded Factors:\n");
			print(new_label);
			cat("\n");
			colnames(recoded_mat)=new_label;

			for(new_name in new_label){
				var_dependencies[[new_name]]=fact_names[i];
			}
			out_factors=cbind(out_factors, recoded_mat);

		}else{
			cat(fact_names[i], ": Leaving as is.\n", sep="");
			out_factors=cbind(out_factors, factors[,i, drop=F]);
		}
	}
	results=list();
	results[["factors"]]=out_factors;
	results[["var_dep"]]=var_dependencies;
	return(results);
}

##############################################################################

remove_outliers=function(factors, outlier_sd=10){

	num_factors=ncol(factors);
	num_samples=nrow(factors);
	new_fact=factors;

	factor_names=colnames(factors);
	sample_names=rownames(factors);

	cat("Outlier Analysis:\n");
	cat("Locating measurements >", outlier_sd, "x standard deviations away from mean.\n");
	cat("\n");

	for(i in 1:num_factors){
		val=factors[,i];
		mean=mean(val, na.rm=T);
		sd=sd(val, na.rm=T);
		
		dev=abs(val-mean);
		outliers=dev>(outlier_sd*sd);
		outlier_ix=which(outliers==T);
		num_outliers=length(outlier_ix);

		if(num_outliers>0){
			cat(num_outliers, " outlier(s) removed for: ", factor_names[i], "   Mean: ", mean, " Std: ", sd, "\n", sep="");

			mat=matrix(val[outlier_ix], nrow=num_outliers, ncol=1, dimnames=list(sample_names[outlier_ix], ""));
			print(mat);
			cat("\n");
			new_fact[outlier_ix,i]=NA;
		}
	}

	return(new_fact);	
}


##############################################################################
# Main Program Starts Here!
##############################################################################

pdf(paste(OutputFnameRoot, ".pdf", sep=""), height=11, width=8.5);

# Load factors
cat("Loading Factors...\n");
factors=load_factors(FactorsFname);
factor_names=colnames(factors);
num_factors=ncol(factors);
num_factor_orig_samples=nrow(factors);
cat(num_factors, " Factor(s) Loaded:\n", sep="");
print(factor_names);
cat("\n");

# Subset factors
if(VariableIncludeListFname!=""){
	variable_subset=scan(VariableIncludeListFname, what=character());
	cat("Variable Inclusion List:\n");
	print(variable_subset);
	cat("\n");
	shared_variables=intersect(factor_names, variable_subset);
	cat("Identified:\n");
	print(shared_variables);

	factors=factors[,shared_variables];
	factor_names=colnames(factors);
	num_factors=ncol(factors);
}

if(VariableExcludeListFname!=""){
	variable_subset=scan(VariableExcludeListFname, what=character());
	cat("Variable Exclusion List:\n");
	print(variable_subset);
	cat("\n");
	remaining_variables=setdiff(factor_names, variable_subset);
	cat("Remaining:\n");
	print(remaining_variables);

	factors=factors[,remaining_variables];
	factor_names=colnames(factors);
	num_factors=ncol(factors);
}

factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);
num_samples=num_factor_samples;

###############################################################################

var_dep=list();

results=recode_non_numeric_factors(factors, var_dep);
factors=results[["factors"]];
var_dep=results[["var_dep"]];
factor_names=colnames(factors);
num_factors=ncol(factors);

for(i in 1:num_factors){
	categories=unique(unique(factors[,i]));
	num_cat= length(categories);
	numNAs=sum(is.na(factors[,i]));
	percNA=round(numNAs/num_samples*100, 2);
	
	cat("\n");
	cat("'", factor_names[i], "' has ", num_cat, " unique values, ", numNAs, " (", percNA, "%) NAs\n", sep="");
	cat("\tUnique: ", paste(head(categories, n=10), collapse=", "), sep="");
	if(num_cat>10){
		cat(", ...");
	}
	cat("\tExample: ", paste(head(factors[,i], n=10), collapse=", "), sep="");
	if(num_cat>10){
		cat(", ...");
	}
	cat("\n");
}

cat("\nProcessing NAs in factors...\n");
factors_preNAproc=factors; # Keep track of values before filling in values


##############################################################################

plot_factor_histograms=function(factors, main=""){
	num_factors=ncol(factors);
	factor_names=colnames(factors);
	orig_par=par(no.readonly=T);
	par(oma=c(0,0,2,0));
	par(mfrow=c(6,5));
	par(mar=c(3,3,1.5,1));
	for(i in 1:num_factors){
		hist(factors[,factor_names[i]],
			xlab="", ylab="", main=factor_names[i], cex.main=.75);	
	}
	par(orig_par);
}

plot_factor_histograms(factors);
outlier_info=capture.output(invisible({factors=remove_outliers(factors)}));
plot_text(outlier_info);
plot_factor_histograms(factors);

##############################################################################

NA_info=capture.output(invisible({factors=process_factor_NAs(factors, MinNonNAProp)}));
plot_text(NA_info);

##############################################################################

CorEffCutoff=.5;
CorPvalCutoff=.05;

correl_res=compute_correl(factors, CorPvalCutoff, CorEffCutoff);
correl=apply(correl_res$correl_mat,c(1,2), function(x){ ifelse(is.na(x),1,x)});
correl=apply(correl, c(1,2), function(x){ round(x,2)});

abs_correl=abs(correl);
avg_correl=apply(abs_correl, 1, mean);
ord_correl=order(avg_correl);
paint_matrix(correl, title="Sorted by Factor Name, Alphabetically");
paint_matrix(correl[ord_correl, rev(ord_correl)], title="Sorted by Absolute Correlation");

##############################################################################

target_variance=.95;

num_factors=ncol(factors);
res=princomp(factors, cor=T, scores=T);
var_per_component=(res$sdev)^2;
tot_var=sum(var_per_component);
prop_var=var_per_component/tot_var;

cat("PCA Components:\n");
cat("Proportion of Variance Explained:\n");
print(prop_var);
cum_var_explained=cumsum(prop_var);
cat("Cumulative Variance Explained:\n");
print(cum_var_explained);
top_components=cum_var_explained<=target_variance;
num_components_to_keep=sum(top_components)+1;

cat("Num Components Kept to Keep ", target_variance*100, "% of Variance: ", num_components_to_keep, "\n", sep="");

bar_colors=rep("grey", num_factors);
bar_colors[num_components_to_keep]="pink";
barplot(cum_var_explained, xlab="Number of Principal Components", ylab="Proportion of Variance Explained", 
	names.arg=1:num_factors, cex.names=.5, col=bar_colors);	
abline(h=target_variance, col="blue");

# pr$sdev, eigenvalues, but you need to square to get variance explained
# pr$loadings, eigenvectors (based on correlation)
# pr$centers (means)
# pr$scale (close to stdev?)
# pr$scores, transformed (correl * eigenvectors)

top2_var_explained=sum(prop_var[1:2]);
plot(res$scores[,1], -res$scores[,2], 
	main=paste("Top 2 Principal Components (", round(top2_var_explained*100,2),"%)"),
	xlab=paste("Component 1 (", round(prop_var[1]*100, 2), "%)"),
	ylab=paste("Component 2 (", round(prop_var[2]*100, 2), "%)")
);

# Compute distance using only top components
pca_dist=dist(res$scores[,1:num_components_to_keep]);
is_zero=(pca_dist==0);
pca_dist[is_zero]=1e-323;

cmds=cmdscale(pca_dist);
plot(cmds, type="n", xlab="Dim 1", ylab="Dim 2", main="NonMetric MDS on Tranformed Coefficients");
points(cmds);

isomds=isoMDS(pca_dist);
plot(isomds$points[,1], isomds$points[,2], type="n", xlab="Dim 1", ylab="Dim 2", main="NonMetric MDS on Tranformed Coefficients");
points(isomds$points[,1], isomds$points[,2]);


##############################################################################

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
