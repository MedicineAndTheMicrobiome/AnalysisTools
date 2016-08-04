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
	"num_variables", "p", 2, "character",
	"outputroot", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table>\n",
	"	-f <factors>\n",
	"	[-p <number of variables>]\n",
	"	[-o <output filename root>]\n",
	"\n",
	"This script will read in the summary file table, then perform\n",
	"a multivariate logistic regression on the the top categories\n",
	"using the factors/predictors in the factor file.\n",
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

if(!length(opt$num_variables)){
	NumVariables=20;
}else{
	NumVariables=opt$num_variables;
}

SummaryFile=opt$summary_file;
FactorsFile=opt$factors;

cat("\n");
cat("Summary File : ", SummaryFile, "\n", sep="");
cat("Factors File: ", FactorsFile, "\n", sep="");
cat("Number of Response Variables: ", NumVariables, "\n", sep="");
cat("Output File: ", OutputRoot, "\n", sep="");
cat("\n");

##############################################################################

load_factors=function(fname){
	factors=data.frame(read.table(fname,  header=TRUE, row.names=1, check.names=FALSE));
	return(factors);
}

load_summary_file=function(fname){
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

# Reorder by abundance
mean_abund=apply(normalized, 2, mean);
ix=order(mean_abund, decreasing=TRUE);
normalized=normalized[,ix];
mean_abund=mean_abund[ix];
sorted_taxa_names=colnames(normalized);

num_top_taxa=NumVariables;
num_top_taxa=min(c(num_top_taxa, num_taxa));
prop_abundance_represented=sum(mean_abund[1:num_top_taxa]);

cat("\nThe top ", num_top_taxa, " taxa are:\n", sep="");
for(i in 1:num_top_taxa){
	cat("\t", sorted_taxa_names[i], "\t[", mean_abund[i], "]\n", sep="");
}
cat("\n");

cat("Accounting for ", prop_abundance_represented, " of taxa.\n", sep="");
cat("\n");

##############################################################################

# Load factors
factors=load_factors(FactorsFile);
factor_names=colnames(factors);
num_factors=ncol(factors);
factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);

cat("\n");
cat(num_factors, " Factor(s) Loaded:\n", sep="");
print(factor_names);
cat("\n");

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
factors=factors[shared_sample_ids,];

##############################################################################

extract_top_categories=function(ordered_normalized, top){
	num_samples=nrow(ordered_normalized);
	num_categories=ncol(ordered_normalized);
	num_saved=min(c(num_categories, top+1));

	top_cat=matrix(0, nrow=num_samples, ncol=num_saved);

	top_cat[,1:top]=ordered_normalized[,1:top];
	top_cat[,(top+1)]=apply(
		ordered_normalized[,(top+1):num_categories],
		1, sum);

	rownames(top_cat)=rownames(ordered_normalized);
	colnames(top_cat)=c(colnames(ordered_normalized)[1:top], "Remaining");
	return(top_cat);
			
}

additive_log_rato=function(ordered_matrix){
	num_cat=ncol(ordered_matrix);
	num_samp=nrow(ordered_matrix);

	denominator=ordered_matrix[,num_cat];
	alr_mat=matrix(0, nrow=num_samp, ncol=(num_cat-1));
	
	for(i in 1:num_samp){
		alr_mat[i,]=log(ordered_matrix[i,1:(num_cat-1)]/denominator[i]);
		#print(alr_mat[i,]);
	}

	rownames(alr_mat)=rownames(ordered_matrix)
	colnames(alr_mat)=head(colnames(ordered_matrix), num_cat-1);

	alr_struct=list();
	alr_struct[["transformed"]]=alr_mat;
	alr_struct[["denominator"]]=denominator;

	return(alr_struct);
}

##############################################################################

# Assign 0's to values smaller than smallest abundance across entire dataset
min_assay=min(normalized[normalized!=0]);
cat("Lowest non-zero value: ", min_assay, "\n\n", sep="");
zero_replacment=min_assay/10;
normalized[normalized==0]=zero_replacment;

##############################################################################

# Output histograms of tranformed data
pdf(paste(OutputRoot, ".mlr.hist.pdf", sep=""), height=11, width=8.5);
responses=extract_top_categories(normalized, num_top_taxa);
resp_alr=additive_log_rato(responses)$transformed;

par(mfrow=c(3,3));
for(i in 1:num_top_taxa){
	hist(resp_alr[,i], breaks=20, xlab="ALR Transformed Abundance", main=sorted_taxa_names[i]);
}

plot_matrix_scatter_plot=T;
if(plot_matrix_scatter_plot){
	# Plot matrix scatter plots
	num_scatter_plot_taxa=min(c(num_top_taxa, 10));
	par(mfrow=c(num_scatter_plot_taxa, num_scatter_plot_taxa));
	par(mar=c(0,0,0,0));

	for(factor_ix in 1:num_factors){
		colors=as.numeric(as.factor(as.vector(factors[,factor_ix])));
		cat("Factors: ", factor_names[factor_ix], "\n", sep="");
		print(colors);
		for(i in 1:num_scatter_plot_taxa){
			for(j in 1:num_scatter_plot_taxa){
				if(i==j){
					plot(0,0, type="n",  xaxt="n", yaxt="n");
					text(0,0, sorted_taxa_names[i], cex=1, srt=45);
				}else{
					plot(resp_alr[,i], resp_alr[,j], cex=.3, xaxt="n", yaxt="n", col=colors);
				}
			}
		}
		cat("ok.\n");
	}
}

dev.off();

##############################################################################

plot_text=function(strings){
	par(family="Courier");
	par(oma=rep(.5,4));
	par(mar=rep(0,4));
	plot(0,0, xlim=c(0,70), ylim=c(0,52), type="n",  xaxt="n", yaxt="n",
		xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
		);
	num_lines=length(strings);
	for(i in 1:num_lines){
		cat(strings[i], "\n", sep="");
		text(0, 53-i, strings[i], pos=4, cex=.8); 
	}
}
	
##############################################################################

pdf(paste(OutputRoot, ".mlr.mmp.pdf", sep=""), height=11, width=8.5);

cat("Performing regression.\n");
cat("Extracting: ", num_top_taxa, " + 1 (remaining) categories.\n", sep="");

responses=extract_top_categories(normalized, num_top_taxa);
resp_alr_struct=additive_log_rato(responses);
transformed=resp_alr_struct$transformed;

model_string= paste("transformed ~", paste(factor_names, collapse=" + "));
cat("\nFitting this multivariate model: ", model_string, "\n");
mv_fit=lm(as.formula(model_string), data=factors);

text=character();
manova_txt=capture.output(anova(mv_fit));

text[1]=paste("Multivariate Regression with ", num_top_taxa, " top taxa", sep="");
text[2]=paste("Proportion of overall mean abundance represented: ", prop_abundance_represented, sep="");
text[3]="";
text=c(text, manova_txt);
plot_text(text);

#plot_text(capture.output(mv_fit));

##############################################################################

plot_correl_heatmap=function(mat, title=""){

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
			str=sprintf("%.2f",mat[i,j]);
			str=gsub("0\\.",".", str);
                        text(i,j,labels=str, cex=cell_cex, srt=45);
                }
        }

        # Plot the labels
        mtext(cnames, at=1:ncols, side=2, las=2, cex=cscale);
        mtext(rnames, at=1:nrows, side=1, las=2, cex=rscale);

	# Plot the title
	mtext(title, line=0, at=nrows*.5, side=3, font=2);

}

cor_mat=cor(transformed);
plot_correl_heatmap(cor_mat, title="Taxonomic Correlations");

##############################################################################

uv_fit=list();

#model_string=paste("transformed[,1] ~", paste(factor_names, collapse=" + "));
#model_matrix=model.matrix(as.formula(model_string), data=factors);

#model_variables=attributes(model_matrix)$dimnames[[2]];
#num_model_variables=length(model_variables);

#print(model_variables);
#print(num_model_variables);

uv_pval_mat=matrix(NA, nrow=num_factors, ncol=num_top_taxa,
		dimnames=list(factor_names, sorted_taxa_names[1:num_top_taxa]));

rsqrd=numeric(num_top_taxa);
adj_rsqrd=numeric(num_top_taxa);

for(var_ix in 1:num_top_taxa){
	summary_txt=character();

	cat("##########################################################################\n");
	cat("#                                                                        #\n");
	cat("# ", sorted_taxa_names[var_ix], "\n");
	cat("#                                                                        #\n");
	cat("##########################################################################\n");

	ALR_Abundance=transformed[,var_ix];
	model_string= paste("ALR_Abundance ~", paste(factor_names, collapse=" + "));
	uv_fit[[var_ix]]=lm(model_string, data=factors);

	uv_summ=summary(uv_fit[[var_ix]]);
	rsqrd[var_ix]=uv_summ$r.squared;
	adj_rsqrd[var_ix]=uv_summ$adj.r.squared;
	
	uv_anova=anova(uv_fit[[var_ix]]);
	pval=uv_anova[["Pr(>F)"]][1:num_factors];
	
	#print(uv_anova);
	#print(pval);

	uv_pval_mat[,var_ix]=pval;
}

# Plot pvalues
plot_correl_heatmap(uv_pval_mat, title="Univariate Pr(>F)");

# Plot R^2
rsqrd_mat=rbind(rsqrd, adj_rsqrd);
rownames(rsqrd_mat)=c("R^2", "Adjusted R^2");
colnames(rsqrd_mat)=sorted_taxa_names[1:num_top_taxa];
plot_correl_heatmap(rsqrd_mat, title="Univariate R Squared");

##############################################################################
# Plot univariate analyses

for(var_ix in 1:num_top_taxa){

	summary_txt=c();

	summary_txt[1]="Univariate Regression:";
	summary_txt[2]="";
	summary_txt[3]=paste(var_ix, ".) ", sorted_taxa_names[var_ix], sep="");
	summary_txt[4]="";
	summary_txt[5]=paste("Mean abundance: ",  sprintf("%3.1f%%",mean_abund[var_ix]*100), sep="");
	summary_txt[6]="";
	summary_txt[7]=paste("R^2: ", sprintf("%3.4f", rsqrd_mat[1,var_ix]), sep="");
	summary_txt[8]=paste("Adjusted R^2: ", sprintf("%3.4f", rsqrd_mat[2, var_ix]), sep="");
	summary_txt[9]="";
	summary_txt=c(summary_txt, capture.output(anova(uv_fit[[var_ix]])));
	plot_text(summary_txt);	

	mmps(uv_fit[[var_ix]], main=paste(var_ix, ".) ", sorted_taxa_names[var_ix], 
		sprintf(" [%3.1f%%]",mean_abund[var_ix]*100), sep=""));
}

##############################################################################
	
dev.off();

##############################################################################

sink(paste(OutputRoot, ".mlr.log.txt", sep=""));

cat("\nFactor information:\n\n");
summary(factors);

sink();

##############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
