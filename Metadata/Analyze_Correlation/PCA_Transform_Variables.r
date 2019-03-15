#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

options(useFancyQuotes=F);
options(digits=5)

params=c(
	"factors", "f", 1, "character",
	"outputroot", "o", 1, "character",
	"predictor", "p", 1, "character",
	"response", "r", 2, "character",
	"pc_coverage", "c", 2, "numeric",
	"export_orig", "e", 2, "character",
	"donnot_transform", "t", 2, "character"
);

NORM_PVAL_CUTOFF=.2;
PCA_COVERAGE=.95;
CURATED_PREFIX="crtd";
NO_CHANGE="orig"

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <factors/metadata file name>\n",
	"	-o <output filename root.\n",
	"\n",
	"	-p <targeted predictor variable list>\n",
	"	[-r <targeted response variable to plot against]\n",
	"\n",
	"	[-c PC Coverage, default=", PCA_COVERAGE, "\n",
	"	[-e (export original predictor variable list, default: remove from output)\n",
	"	[-t (do not transform/autocurate predictors, default: transform if necessary)\n",
	"\n",
	"This script will take in the specified predictor variable\n",
	"list and perform a PCA, keeping the top PCs.\n",
	"\n",
	"If the -t option is not selected, the variables will be\n",
	"automatically checked for non-normality and log'd for improvement.\n",
	"\n");

if(
	!length(opt$factors) || 
	!length(opt$outputroot) || 
	!length(opt$predictor)
){
	cat(usage);
	q(status=-1);
}

FactorsFname=opt$factors;
OutputFnameRoot=opt$outputroot;
ResponseListName=opt$response;

PredictorListName="";
PCCoverage=PCA_COVERAGE;
DonnotTransform=F;
ExportOrig=F;

if(length(opt$predictor)){
	PredictorListName=opt$predictor;
}

if(length(opt$pc_coverage)){
	PCCoverage=opt$pc_coverage;
}

if(length(opt$donnot_transform)){
	DonnotTransform=T;
}

if(length(opt$export_orig)){
	ExportOrig=T;
}


param_text=capture.output({
	cat("\n");
	cat("Factor File Name: ", FactorsFname, "\n");
	cat("Output File Name Root: ", OutputFnameRoot, "\n");
	cat("Response List Name: ", ResponseListName, "\n");
	cat("Predictor List Name: ", PredictorListName, "\n");
	cat("PC Min Coverage: ", PCCoverage, "\n");
	cat("\n");
	cat("Donnot Transform Variables: ", DonnotTransform, "\n");
	cat("Export original variables: ", ExportOrig, "\n");
	cat("\n");
});

cat(param_text);

###############################################################################

load_factors=function(fname){
	factors=as.data.frame(read.delim(fname,  header=TRUE, row.names=1, check.names=FALSE, sep="\t"));
	return(factors);
}

load_list=function(fname){
	lst=scan(fname, what=character());
	return(lst);	
}

test_and_apply_log_transform=function(mat_val, pval_cutoff=.2){
	nrows=nrow(mat_val);
	ncols=ncol(mat_val);

	trans_mat=mat_val;
	orig_colnames=colnames(mat_val);
	new_colnames=character();

	for(var in orig_colnames){
		values=mat_val[,var];

		num_unique_val=length(setdiff(unique(values), NA));
		cat("\n", var, ": Num Unique Values: ", num_unique_val, "\n");

		test_res=shapiro.test(values);

		if(test_res$p.value<=pval_cutoff && num_unique_val>2){
			cat(" Not normal: ", test_res$p.value, "\n");

			log_values=log(values+1);
			test_log_res=shapiro.test(log_values);

			if(test_log_res$p.value < test_res$p.value){
				# Keep original
				cat("  No Improvement: ", test_log_res$p.value, "\n");
				new_colnames=c(new_colnames, paste("orig_", var, sep=""));
			}else{
				# Keep log transformed
				cat("  Transformation Effective: ", test_log_res$p.value, "\n");
				new_colnames=c(new_colnames, paste("log_", var, sep=""));
				trans_mat[, var]=log_values;
			}
		}else{
			cat(" Normal enough. ", test_res$p.value, "\n");
			new_colnames=c(new_colnames, var);
		}
	}
	colnames(trans_mat)=new_colnames;
	return(trans_mat);
}

plot_text=function(strings){

	orig.par=par(no.readonly=T);
	
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

	par(orig.par);
}

compute_correlations=function(mat){
	num_col=ncol(mat);
	cor_mat=matrix(0, nrow=num_col, ncol=num_col);
	pval_mat=matrix(0, nrow=num_col, ncol=num_col);
	rownames(cor_mat)=colnames(mat);
	colnames(cor_mat)=colnames(mat);
	rownames(pval_mat)=colnames(mat);
	colnames(pval_mat)=colnames(mat);
	for(i in 1:num_col){
		for(j in 1:i){
			v1=mat[,i];
			v2=mat[,j];
			notna=!(is.na(v1) | is.na(v2));
			#cor_mat[i,j]=cor(v1[notna], v2[notna]);
			test=cor.test(v1[notna], v2[notna]);
			pval_mat[i,j]=test$p.value;
			pval_mat[j,i]=test$p.value;
			cor_mat[i,j]=test$estimate;
			cor_mat[j,i]=test$estimate;
		}
	}
	res=list();
	res[["val"]]=cor_mat;
	res[["pval"]]=pval_mat;;
	res[["dist"]]=as.dist(1-abs(cor_mat));
	return(res);
}


##############################################################################
# Main Program Starts Here!
##############################################################################

pdf(paste(OutputFnameRoot, ".pdf", sep=""), height=11, width=8.5);

plot_text(param_text);

# Load factors
cat("Loading Factors...\n");
loaded_factors=load_factors(FactorsFname);
loaded_factor_names=colnames(loaded_factors);
loaded_sample_names=rownames(loaded_factors);

cat("Loaded factors:\n");
print(loaded_factor_names);
cat("\n");
cat("Loaded sample ids:\n");
print(loaded_sample_names);
cat("\n");

# Subset factors
predictors_arr=load_list(PredictorListName);

responses_arr=c();
if(ResponseListName!=""){
	responses_arr=load_list(ResponseListName);
}else{
	cat("Response variable list not specified.\n");
}


cat("\n");
cat("Targeted Responses:\n");
print(responses_arr);
missing=setdiff(responses_arr, loaded_factor_names);
if(length(missing)>0){
	cat("Missing:\n");
	print(missing);
	quit(status=-1);
}
cat("\n");

cat("Targeted Predictors:\n");
print(predictors_arr);
missing=setdiff(predictors_arr, loaded_factor_names);
if(length(missing)>0){
	cat("Missing:\n");
	print(missing);
	quit(status=-1);
}
cat("\n");

num_pred=length(predictors_arr);
num_resp=length(responses_arr);

pred_mat=loaded_factors[, predictors_arr, drop=F];
resp_mat=loaded_factors[, responses_arr, drop=F];

orig_pred_names=predictors_arr;
orig_resp_names=responses_arr;

cat("Testing Predictor Variables for normality.\n");
curated_pred_mat=test_and_apply_log_transform(pred_mat, NORM_PVAL_CUTOFF);

cat("Testing Response Variables for normality.\n");
curated_resp_mat=test_and_apply_log_transform(resp_mat, NORM_PVAL_CUTOFF);

##############################################################################

curated_predictors_arr=colnames(curated_pred_mat);
curated_responses_arr=colnames(curated_resp_mat);

##############################################################################

plot_text(c(
	"Relationship between predictors (x-axis) and responses (y-axis):",
	"",
	"Plots are sorted by decreasing statistical significance (increasing p-value).",
	"",
	"To determine if the Log transform of the predictors or responses were necessary,",
	"  the following algorithm was applied:",
	"",
	"  If the Shapiro-Wilks (SW) Test for normality rejects the distribution as normal",
	paste("  at alpha<", NORM_PVAL_CUTOFF, ", then the log transform is attempted.", sep=""),
	"  If the SW Test p-value on the transformed distribution is greater than that of",
	"  the untransformed distribution, then the transformed distribution is retained,",
	"  else the untransformed distribution is retained."
));

par(oma=c(1,1,4,1));

pred_ix=1;
for(pred_name in curated_predictors_arr){
	par(mfrow=c(4,3));

	pval_arr=numeric(num_resp);
	coef_arr=numeric(num_resp);
	rsqd_arr=numeric(num_resp);
	fit_arr=list();

	pred_val=curated_pred_mat[,pred_name];

	# Fit against each response variable
	i=1;
	for(resp_name in curated_responses_arr){

		cat("Working on: ", resp_name, " vs ", pred_name, "\n");
		resp_val=curated_resp_mat[,resp_name];

		# Fit and summarize
		fit=lm(resp_val~pred_val);
		sumfit=summary(fit);

		# Pull regression stats
		pval_arr[i]=sumfit$coefficients["pred_val", "Pr(>|t|)"];
		coef_arr[i]=sumfit$coefficients["pred_val", "Estimate"];
		rsqd_arr[i]=sumfit$r.squared;
		fit_arr[[i]]=fit;
	
		i=i+1;
	}

	# Sort by increasing pvalue 
	sort_ix=order(pval_arr, decreasing=F);

	# Plot each fit by decreasing significance
	for(i in sort_ix){
		resp_val=curated_resp_mat[,i];
		resp_name=curated_responses_arr[i];

		signf_char="";
		if(pval_arr[i]<=.001){
			signf_char=" ***";
		}else if(pval_arr[i]<=.01){
                        signf_char=" **";
		}else if(pval_arr[i]<=.05){
                        signf_char=" *";
		}else if(pval_arr[i]<=.1){
                        signf_char=" '";
		}

		plot(pred_val, resp_val, xlab=pred_name, ylab=resp_name, 
			main=paste(orig_resp_names[i], signf_char, sep=""));

		# Draw regression line
		abline(fit_arr[[i]], col="blue");
		stat_info=paste(
			"coeff=", signif(coef_arr[i],2), 
			"  p-val=", sprintf("%3.3f", pval_arr[i]),
			"  R^2=", signif(rsqd_arr[i],2), 
			sep="");

		mtext(stat_info, side=3, outer=F, cex=.5, col="blue");
		mtext(orig_pred_names[pred_ix], side=3, outer=T, cex=2);
	}

	pred_ix=pred_ix+1;
}

#############################################################################


impute_cell=function(target_predictors, responses, predictors, verbose=F){

	num_samples=nrow(responses);
	avail_predictors=ncol(predictors);

	num_pred_to_use=min(num_samples, avail_predictors)-2;
	
	cat("Num Samples Available for Imputation: ", num_samples, "\n");
	cat("Num Available Predictors: ", avail_predictors, "\n");
	cat("Num Predictors to Use: ", num_pred_to_use, "\n");
	
	# Compute correlation between response and predictors
	corr=apply(predictors, 2, function(x){ 
		non_na=!is.na(x); 
		cor(x[non_na], responses[non_na,]);}
		);

	# Order predictors by decreasing correlation
	corr_order=order(abs(corr), decreasing=T);
	corr_ordered=(corr[corr_order]);
	#print(corr_ordered);
	top_pred=names(corr_ordered)[1:num_pred_to_use];


	# Replace NAs with median value from the same column/predictor
	predictors_nona=predictors;
	for(j in 1:ncol(predictors_nona)){
		cur_median=median(predictors[,j], na.rm=T);	
		predictors_nona[is.na(predictors[,j]),j]=cur_median;
	}

	# Create model with the predictors with greatest correlaton with response first
	form_str=paste(colnames(responses), "~", paste(top_pred, collapse="+"));
	print(form_str);
	fit=lm(formula(form_str), data=cbind(responses, predictors_nona));

	# Predict NA with the values that we have
	imputed_val=predict(fit, new=target_predictors[top_pred]);

	#print(fit);
	if(verbose){
		print(summary(fit));
	}
	#eg=rbind(predictors, target_predictors);
	#print(eg);
	#x=predict(fit, new=eg);
	obs_resp_range=range(responses);
	cat("--------------------------------------------------------\n");
	cat("Name: ", colnames(responses), "\n");
	cat("Response Range: ", obs_resp_range[1], " - ", obs_resp_range[2], "\n");
	cat("Inputed Value: ", imputed_val, "\n");
	if(imputed_val<obs_resp_range[1] || imputed_val>obs_resp_range[2]){
		cat("WARNING: Inputed Value Outside Range of Observed Values\n");
	}
	cat("--------------------------------------------------------\n");
	return(imputed_val);
} 

impute_matrix=function(mat_wna){

	num_rows=nrow(mat_wna);
	num_cols=ncol(mat_wna);

	cat("Original Matrix Dimensions: ", num_rows, " x ", num_cols, 
		" matrix. (", num_rows*num_cols, ")\n", sep="");

	# remove rows with all NAs
	usable_rows=apply(mat_wna, 1, function(x){!all(is.na(x))});
	usable_cols=apply(mat_wna, 2, function(x){!all(is.na(x))});
	usable_mat_wna=mat_wna[usable_rows,usable_cols];

	usable_num_rows=nrow(usable_mat_wna);
	usable_num_cols=ncol(usable_mat_wna);

	cat("Removed rows/cols with all NAs: ", usable_num_rows, " x ", usable_num_cols, 
		" matrix. (", usable_num_rows * usable_num_cols,")\n", sep="");

	#print(mat_wna);

	cat("Looking for NAs...\n");
	na_pos=numeric();

	if(1){
		for(rix in 1:usable_num_rows){
			for(cix in 1:usable_num_cols){
				if(is.na(usable_mat_wna[rix, cix])){
					na_pos=rbind(na_pos, c(rix, cix));
				}
			}
		}
	}else{
		# For validation
		for(rix in 1:usable_num_rows){
			for(cix in 1:usable_num_cols){
				na_pos=rbind(na_pos, c(rix, cix));
			}
		}
	}
	
	num_nas_to_impute=nrow(na_pos);
	cat("Num NAs to try to impute:",  num_nas_to_impute, "\n");

	filled_matrix=usable_mat_wna;
	for(na_ix in 1:num_nas_to_impute){

		target_row=na_pos[na_ix,1];
		target_column=na_pos[na_ix,2];

		cell=usable_mat_wna[target_row, target_column, drop=F];

		if(!is.na(cell)){
			cat("Error:  Trying to input cell not NA.\n");
			quit();
		}
	
		cat("(", na_ix, "/", num_nas_to_impute, ") Imputing: ", rownames(cell), " / ", colnames(cell), "\n");

		non_na_row=!is.na(usable_mat_wna[,target_column,drop=F]);
		non_na_col=!is.na(usable_mat_wna[target_row,,drop=F]);

		imputed_val=impute_cell(
			target_predictors=usable_mat_wna[target_row, non_na_col, drop=F], 
			responses=usable_mat_wna[non_na_row, target_column, drop=F],
			predictors=usable_mat_wna[non_na_row, non_na_col, drop=F]
			);

		filled_matrix[target_row, target_column]=imputed_val;

	}

	repl_rows=rownames(filled_matrix);
	mat_wna[repl_rows,]=filled_matrix[repl_rows,];
	return(mat_wna);
	
}

imputed_mat=impute_matrix(curated_pred_mat);
print(imputed_mat);

##############################################################################

# Compute the dendrogram, clustering based on the correlation between responses and predictors

plot_text(c(
	"The following dendrogram illustrates the relationship between the",
	"  response and predictor variables based on their correlation.",
	"",
	"The greater the magnitude of the correlation, i.e. abs(cor),",
	"  the shorter the distance between variables.",
	"",
	"The Euclidean distance is still used to compare the correlation profiles",
	"  between variables to satisfy the 'triangle inequality'.",
	"",
	"The predictor variables are colored red.",
	"Note that the relationship between variables is actually closer to R^2",
	"  than the estimated coefficient magnitudes because the scale/unit for",
	"  each of the variables may not be the same."
));

correl=compute_correlations(cbind(curated_pred_mat, curated_resp_mat));

par(mfrow=c(1,1));
par(mar=c(15,2,1,2));
hcl=hclust(correl$dist, method="ward.D2");
dend=as.dendrogram(hcl);

highlight_predictors=function(x){
	if(is.leaf(x)){
		leaf_attr=attributes(x);
		label=leaf_attr$label;
		print(label);
		if(any(label==curated_predictors_arr)){
			color="black";
			font=2;
		}else{
			color="red";
			font=1;
		}
		attr(x, "nodePar")=c(leaf_attr$nodePar, list(lab.font=font, lab.col=color, cex=0));
	}
	return(x);
}

dend=dendrapply(dend, highlight_predictors);

plot(dend, main="Ward's Minimum Variance: dist(1-abs(cor))");


##############################################################################

pred_correl=compute_correlations(curated_pred_mat);

# Note that the princomp R function takes the resp values directly and the Standard Deviations
# are the squareroot of the eigenvalues

# Component Loading: num_pred x num_pred matrix, i.e. correlation between num_pred and PCs
#	Loadings and eigenvectors contain similiar information (and may be equal, but are not the same)
# Use varimax to realign loadings, since there are multiple solutions to component loadings
#
# Component Score: num_obs x num_pred matrix
#	scale(pred_mat, center, scale) %*% loadings

# Compute PCA
eigen=eigen(pred_correl$val);

# Compute variance contribution of each PC
pca_propvar=eigen$values/sum(eigen$values);
pca_propcumsum=cumsum(pca_propvar);
num_pc_at_cutoff=sum(pca_propcumsum<PCCoverage);

# Compute per sample scores
scores=(scale(imputed_mat, center=T, scale=T) %*% eigen$vectors);

plot_text(c(
	"Principal Components Analysis:",
	"(Eigenvalues/vectors on Response Correlation Matrix)",
	"",
	"Proportion of Variance in each PC:",
	capture.output(print(pca_propvar)),
	"",
	"Cumulative Variance:",
	capture.output(print(pca_propcumsum)),
	"",
	paste("Number of PCs to retain >", (PCCoverage*100.0), "% of Variance:"),
	num_pc_at_cutoff
));

# Plot bar plots of PC variance explanation
par(mfrow=c(2,1));
par(mar=c(7,4,2,2));
colors=rep("grey",num_pred);
colors[1:num_pc_at_cutoff]="darkcyan";

mids=barplot(pca_propvar, las=2, names.arg=1:num_pred, xlab="PCs", 
	col=colors,
	ylab="Proportion", main="PCA Proportion of Variance");

mids=barplot(pca_propcumsum, las=2, names.arg=1:num_pred, xlab="PCs", 
	col=colors,
	ylab="Proportion", main="PCA Cumulative Variance");
abline(h=PCCoverage, col="blue", lty=2);

# Plot sample ordination based on scores
sample_ids=rownames(resp_mat);
par(mfrow=c(3,2));
par(mar=c(3,3,1,1));

for(i in seq(1,num_pc_at_cutoff+1,2)){
	xpos=scores[,i];
	ypos=scores[,i+1];
	xrange=range(xpos, na.rm=T);
	yrange=range(ypos, na.rm=T);

	xspan=diff(xrange);
	yspan=diff(yrange);

	# Plot labelled samples
	plot(xpos, ypos, type="n", 
		xlim=c(xrange[1]-xspan/10, xrange[2]+xspan/10),
		ylim=c(yrange[1]-yspan/10, yrange[2]+yspan/10),
		xlab="", ylab="", main=""
	)

	title(
		xlab=paste("PC",i,sep=""), 
		ylab=paste("PC",i+1,sep=""), 
		line=2
	);
	text(xpos, ypos, sample_ids, cex=.7);

	# Plot points
	plot(xpos, ypos, type="p", 
		xlim=c(xrange[1]-xspan/10, xrange[2]+xspan/10),
		ylim=c(yrange[1]-yspan/10, yrange[2]+yspan/10),
		xlab="", ylab="", main=""
	)

	title(
		xlab=paste("PC",i,sep=""), 
		ylab=paste("PC",i+1,sep=""), 
		line=2
	);

}

##############################################################################







quit();
##############################################################################

clean_resp_mat=resp_mat;
clean_resp_names=gsub("\\(","_", colnames(resp_mat));
clean_resp_names=gsub("\\)","", clean_resp_names);
colnames(clean_resp_mat)=clean_resp_names;

clean_pred_mat=pred_mat;
clean_pred_names=gsub("\\(","_", colnames(pred_mat));
clean_pred_names=gsub("\\)","", clean_pred_names);
colnames(clean_pred_mat)=clean_pred_names;

##############################################################################

cat("Outputing Factor file with PCs included and Resp Variables removed...\n");

orig_factor_variables=colnames(loaded_factors);
new_factor_variables=setdiff(orig_factor_variables, orig_resp_names);

colnames(scores)=paste("PC",1:num_resp, sep="_");
out_factors=cbind(loaded_factors[,new_factor_variables], scores[,1:num_pc_at_cutoff]);

fname=paste(OutputFnameRoot, ".pca.tsv", sep="");
fh=file(fname, "w");
cat(file=fh, "SampleID");
close(fh);
write.table(out_factors, file=fname, col.names=NA, append=T, quote=F, sep="\t");

##############################################################################

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
