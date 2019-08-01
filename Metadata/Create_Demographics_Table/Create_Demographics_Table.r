#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

options(useFancyQuotes=F);
options(digits=5)

params=c(
	"factors", "f", 1, "character",
	"outputroot", "o", 2, "character",
	"exclude", "x", 2, "character",
	"include", "n", 2, "character",
	"split_var", "s", 2, "character"
);

DEF_MIN_NONNA_PROP=.65;

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <factor/metadata tsv file>\n",
	"	[-o <output root>]\n",
	"\n",
	"	Specify which subset of variables to keep/remove:\n",
	"	[--include <subset variables list filename>]\n",
	"	[--exclude <file with variables to exclude>]\n",
	"\n",
	"	[-s <variable name for splitting/comparing>]\n",
	"\n",
	"This script will read in a metadata file and generate\n",
	"a summary demographics table for each of the columns\n",
	"in the input factor/metadata file.\n",
	"\n",
	"For Boolean variables, the count and percent will be reported.\n",
	"For categorical variables, the categories will be reported separately.\n",
	"For integer/real variables, the mean and 95%CI (assuming normal dist) will be reported.\n",
	"\n",
	"If the split variable name is specified, then ANOVA p-value will be calculated\n",
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

if(!length(opt$split)){
	SplitVarName="";
}else{
	SplitVarName=opt$split;
}

if(SplitVarName!=""){
	OutputFnameRoot=paste(OutputFnameRoot, ".", SplitVarName, sep="");
}


FactorsFname=opt$factors;

cat("Factors Filename/Metadata: ", FactorsFname, "\n", sep="");
cat("Output Root Filename: ", OutputFnameRoot, "\n", sep="");

if(VariableIncludeListFname!=""){
	cat("Using subset of variables from: ", VariableIncludeListFname, " (Inclusion List)\n");
}

if(VariableExcludeListFname!=""){
	cat("Excluding subset of variables from: ", VariableExcludeListFname, " (Exclusion List)\n");
}

if(SplitVarName!=""){
	cat("Splitting data by: ", SplitVarName, "\n");
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

		num_lines=length(strings);

		top=max_lines;

		plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
			xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
			);
		for(i in 1:num_lines){
			#cat(strings[i], "\n", sep="");
			text(0, top-i*6, strings[i], pos=4, cex=1);
		}

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
# Main Program Starts Here!
##############################################################################

# Load factors
cat("Loading Factors...\n");
factors=load_factors(FactorsFname);
factor_names=colnames(factors);
num_factors=ncol(factors);
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

# Create matrix of sample membership based on split variable 
if(SplitVarName==""){
	sp_ix=matrix(rep(T, num_samples), nrow=num_samples, ncol=1);
	rownames(sp_ix)=rownames(factors);
	colnames(sp_ix)=SplitVarName;
}else{
	split_vals=factors[,SplitVarName];
	splt_tab=table(split_vals);
	print(splt_tab);

	num_uniq_vals=length(splt_tab);
	uniq_vals=names(splt_tab);
	sp_ix=matrix(NA, nrow=num_samples, ncol=num_uniq_vals);
	rownames(sp_ix)=rownames(factors);
	colnames(sp_ix)=uniq_vals;

	for(val in uniq_vals){
		sp_ix[,val]=val==(split_vals)
	}

	cat("Removing split variable from variables to place into table.\n");
	new_fact_names=setdiff(colnames(factors), SplitVarName);
	factors=factors[, new_fact_names];
}


###############################################################################

calc_bin_conf=function(counts){

	dimen=dim(counts);
	ub_mat=matrix(NA, nrow=dimen[1], ncol=dimen[2]);
	lb_mat=matrix(NA, nrow=dimen[1], ncol=dimen[2]);

	rownames(ub_mat)=rownames(counts);
	rownames(lb_mat)=rownames(counts);
	colnames(ub_mat)=colnames(counts);
	colnames(lb_mat)=colnames(counts);

	colsums=apply(counts, 2, sum);

	for(c_ix in 1:dimen[2]){
		for(r_ix in 1:dimen[1]){
			bt_res=binom.test(counts[r_ix, c_ix], colsums[c_ix]);
			lb_mat[r_ix, c_ix]=bt_res$conf.int[1];
			ub_mat[r_ix, c_ix]=bt_res$conf.int[2];
		}
	}

	res=list();
	res[["lb"]]=lb_mat;
	res[["ub"]]=ub_mat;
	return(res);
}

###############################################################################

num_factors=ncol(factors);
num_splits=ncol(sp_ix);

split_categories=colnames(sp_ix);

cat("Split Categories:\n");
print(split_categories);
cat("\n");

factor_names=colnames(factors);
MAX_LEVELS=10;

demo_list=list();

# If continous variable:
#   calc mean, 95%CI of mean, ANOVA
#
# if categorical or boolean:
#   count and percent, 95% CI of mean, chi-squared

cat_stats=c("count", "percent", "lb", "ub");
num_cat_stats=length(cat_stats);
num_stats=c("mean", "lb", "ub")

for(f_ix in 1:num_factors){

	cur_factor=factors[,f_ix];
	cur_factname=factor_names[f_ix];
	uniq_vals=sort(unique(cur_factor[!is.na(cur_factor)]));
	num_uniq=length(uniq_vals);

	cat("Working on: ", cur_factname, "\n");
	isfactor=is.factor(cur_factor);
	
	if(isfactor && num_uniq>MAX_LEVELS){
		cat("  Factor has too many levels (", num_uniq, ").\n  Skipping...\n", sep="");
		next;
	}

	if(isfactor || num_uniq==2){
		cat("Treating as categorical.\n");

		counts_mat=matrix(0, ncol=num_splits, nrow=num_uniq);
		colnames(counts_mat)=split_categories;
		rownames(counts_mat)=uniq_vals;

		for(splcat in split_categories){
			cat_ix=sp_ix[,splcat];
			cat_fact_vals=cur_factor[cat_ix];
			cat_fact_vals=cat_fact_vals[!is.na(cat_fact_vals)];
			tab=table(cat_fact_vals);
			counts_mat[names(tab), splcat]=tab;
		}

		#print(counts_mat);
		colcounts=apply(counts_mat, 2, sum);

		binom_ci=calc_bin_conf(counts_mat);

		prop_mat=counts_mat;
		for(i in 1:num_splits){
			prop_mat[,i]=counts_mat[,i]/colcounts[i];
		}
		#print(prop_mat);
		chisq_res=chisq.test(counts_mat);
		#print(chisq_res);

		res=list();
		res[["type"]]="categorical";
		res[["counts"]]=counts_mat;
		res[["prop"]]=prop_mat;
		res[["pval"]]=chisq_res$p.value;
		res[["lb"]]=binom_ci[["lb"]];
		res[["ub"]]=binom_ci[["ub"]];

		demo_list[[cur_factname]]=res;

	}else{
		cat("Treating as continuous\n");

		means=numeric(num_splits);
		counts=numeric(num_splits);
		lb=numeric(num_splits);
		ub=numeric(num_splits);
	
		names(means)=split_categories;
		names(counts)=split_categories;
		names(lb)=split_categories;
		names(ub)=split_categories;
	
		for(splcat in split_categories){
			cat_ix=sp_ix[,splcat];
			cat_cont_vals=cur_factor[cat_ix];
			cat_cont_vals=cat_cont_vals[!is.na(cat_cont_vals)];
			means[splcat]=mean(cat_cont_vals);
			counts[splcat]=length(cat_cont_vals);
			
			tres=tryCatch({
				t.test(cat_cont_vals);
			}, error=function(e){
			});

			if(!is.null(tres)){
				lb[splcat]=tres[["conf.int"]][1];
				ub[splcat]=tres[["conf.int"]][2];
			}else{
				lb[splcat]=NA;
				ub[splcat]=NA;
			}
		}

		# overall anova
		nona=!is.na(cur_factor) & !is.na(split_vals);
		xval=split_vals[nona];
		yval=cur_factor[nona];
		fit=lm(yval~xval);
		anova_res=anova(fit);

		res=list();
		res[["type"]]="continuous";
		res[["means"]]=means;
		res[["counts"]]=counts;
		res[["lb"]]=lb;
		res[["ub"]]=ub;
		res[["pval"]]=anova_res[["Pr(>F)"]][1];

		demo_list[[cur_factname]]=res;

	}

}

#print(demo_list);

sigchar=function(x){
	if(is.na(x) || is.nan(x)){
		return("");
	}

	if(x<=.001){
		return("***");
	}else if(x<=.010){
		return("**");
	}else if(x<=.050){
		return("*");
	}else if(x<=.1){
		return(".");
	}else{
		return("");
	}
}


output_demo_list=function(dl, splts, fh){

	list_names=names(dl);

	num_splts=length(splts);

	cat(file=fh, "Variable\tCategories\t", 
		paste(splts,collapse="\t"), "\t",
		paste("N", "p-value", "Signf", "Method", "DummyName", sep="\t"),
		"\n", sep="");

	for(nm in list_names){

		info=dl[[nm]];
		col_str=character(num_splts);
		names(col_str)=splts;

		if(info[["type"]]=="categorical"){

			mn=info[["prop"]];
			ct=info[["counts"]];
			lb=info[["lb"]];
			ub=info[["ub"]];
			pv=info[["pval"]];

			num_rows=nrow(mn);
			cat_names=rownames(mn);

			cat(file=fh, nm, ":\t", sep="");

			for(cat in cat_names){
				totaln=sum(ct[cat,]);
				for(sp in splts){
					col_str[sp]=paste(
						round(mn[cat,sp], 3), 
						" (", round(lb[cat,sp],3), ", ", 
							round(ub[cat,sp],3), ") [",	
							ct[cat,sp], "]", sep="");
				}

				if(cat!=cat_names[1]){
					cat(file=fh, "\t");
				}
				cat(file=fh, cat, col_str, totaln, pv, sigchar(pv), "ChiSq", 
					paste(nm,".",cat,sep=""), sep="\t"); 
				cat(file=fh, "\n");
			}
		
		} else if(info[["type"]]=="continuous"){
	
			mn=info[["means"]];
			ct=info[["counts"]];
			lb=info[["lb"]];
			ub=info[["ub"]];
			pv=info[["pval"]];

			totaln=sum(ct);

			for(sp in splts){
				col_str[sp]=paste(
					round(mn[sp], 3), 
					" (", round(lb[sp],3), ", ", round(ub[sp],3), ") [", ct[sp], "]", sep="");
			}
			cat(file=fh, paste(nm, ":\t", sep=""), col_str, totaln, pv, sigchar(pv), "ANOVA", 
					nm, sep="\t"); 
			cat(file=fh, "\n");
		}


	}
}

outfh=file(paste(OutputFnameRoot, ".demo_tab.tsv", sep=""), "w");

output_demo_list(demo_list, split_categories, outfh);

close(outfh);





##############################################################################

cat("Done.\n");

print(warnings());
q(status=0);
