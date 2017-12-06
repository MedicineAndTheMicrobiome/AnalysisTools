#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

options(useFancyQuotes=F);
options(digits=5)

params=c(
	"factors", "f", 1, "character",
	"outputroot", "o", 1, "character",
	"response", "r", 1, "character",
	"predictor", "p", 2, "character",
	"pc_coverage", "c", 2, "numeric"
);

PCA_COVERAGE=.95;

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <factors/metadata file name>\n",
	"	-o <output filename root.\n",
	"\n",
	"	-r <response variable list>\n",
	"	[-p <predictor variable list>]\n",
	"\n",
	"	[-c PC Coverage, default=", PCA_COVERAGE, "\n",
	"\n",
	"This script will perform a PCA analysis\n",
	"across the response variables and\n",
	"recommend a set of principal coordinates.\n",
	"\n",
	"If the predictor variable list is specified\n",
	"then plots between each predictor (x) and\n",
	"all the responses (y) will be made.\n",
	"\n",
	"The predictors will be clustered based on\n",
	"their correlation magnitude.\n",
	"\n");

if(
	!length(opt$factors) || 
	!length(opt$outputroot) || 
	!length(opt$response)
){
	cat(usage);
	q(status=-1);
}

FactorsFname=opt$factors;
OutputFnameRoot=opt$outputroot;
ResponseListName=opt$response;

if(length(opt$predictor)){
	PredictorListName=opt$predictor;
}else{
	PredictorListName="";
}

if(length(opt$pc_coverage)){
	PCCoverage=opt$pc_coverage;
}else{
	PCCoverage=PCA_COVERAGE;
}

cat("\n");
cat("Factor File Name: ", FactorsFname, "\n");
cat("Output File Name Root: ", OutputFnameRoot, "\n");
cat("Response List Name: ", ResponseListName, "\n");
cat("Predictor List Name: ", PredictorListName, "\n");
cat("PC Min Coverage: ", PCCoverage, "\n");
cat("\n");

###############################################################################

load_factors=function(fname){
	factors=as.data.frame(read.delim(fname,  header=TRUE, row.names=1, check.names=FALSE, sep="\t"));
	return(factors);
}

load_list=function(fname){
	lst=scan(fname, what=character());
	return(lst);	
}

test_log=function(mat_val){
	nrows=nrow(mat_val);
	ncols=ncol(mat_val);

	trans_mat=mat_val;
	orig_colnames=colnames(mat_val);
	new_colnames=character();

	for(var in orig_colnames){
		values=mat_val[,var];

		num_unique_val=length(setdiff(unique(values), NA));
		cat(var, ": Num Unique Values: ", num_unique_val, "\n");

		test_res=shapiro.test(values);

		if(test_res$p.value<.20 && num_unique_val>2){
			log_values=log(values+1);
			test_log_res=shapiro.test(log_values);

			if(test_log_res$p.value < test_res$p.value){
				# Keep original
				new_colnames=c(new_colnames, var);
			}else{
				# Keep log transformed
				new_colnames=c(new_colnames, paste("log(", var, ")", sep=""));
				trans_mat[, var]=log_values;
			}
		}else{
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

##############################################################################
# Main Program Starts Here!
##############################################################################

pdf(paste(OutputFnameRoot, ".pdf", sep=""), height=11, width=8.5);

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
responses_arr=load_list(ResponseListName);
predictors_arr=load_list(PredictorListName);

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

pred_mat=test_log(pred_mat);
resp_mat=test_log(resp_mat);

predictors_arr=colnames(pred_mat);
responses_arr=colnames(resp_mat);

plot_text(c(
	paste("Factor File Name: ", FactorsFname),
	paste("Output File Name Root: ", OutputFnameRoot),
	paste("Response List Name: ", ResponseListName),
	paste("Predictor List Name: ", PredictorListName),
	paste("PC Min Coverage: ", PCCoverage),
	"",
	"Targeted Responses:",
	capture.output(responses_arr),
	"",
	"Targeted Predictors:",
	capture.output(predictors_arr)
));

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
	"  at alpha<.2, then the log transform is attempted.",
	"  If the SW Test p-value on the transformed distribution is greater than that of",
	"  the untransformed distribution, then the transformed distribution is retained,",
	"  else the untransformed distribution is retained."
));

par(oma=c(1,1,4,1));

pred_ix=1;
for(pred in predictors_arr){
	par(mfrow=c(4,3));

	pval_arr=numeric(num_resp);
	coef_arr=numeric(num_resp);
	rsqd_arr=numeric(num_resp);
	fit_arr=list();

	i=1;
	for(resp in responses_arr){

		cat("Working on: ", resp, " vs ", pred, "\n");
		pred_val=pred_mat[,pred];
		resp_val=resp_mat[,resp];

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

	# Plot
	for(i in sort_ix){
		resp_val=resp_mat[,i];
		resp_name=responses_arr[i];

		plot(pred_val, resp_val, xlab=pred, ylab=resp_name, main=orig_resp_names[i]);

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

##############################################################################

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

correl=compute_correlations(cbind(pred_mat, resp_mat));

par(mfrow=c(1,1));
par(mar=c(15,2,1,2));
hcl=hclust(correl$dist, method="ward.D2");
dend=as.dendrogram(hcl);

highlight_predictors=function(x){
	if(is.leaf(x)){
		leaf_attr=attributes(x);
		label=leaf_attr$label;
		print(label);
		if(any(label==predictors_arr)){
			color="red";
			font=2;
		}else{
			color="black";
			font=1;
		}
		attr(x, "nodePar")=c(leaf_attr$nodePar, list(lab.font=font, lab.col=color, cex=0));
	}
	return(x);
}

dend=dendrapply(dend, highlight_predictors);

plot(dend, main="Ward's Minimum Variance: dist(1-abs(cor))");


##############################################################################

resp_correl=compute_correlations(resp_mat);

# Note that the princomp R function takes the resp values directly and the Standard Deviations
# are the squareroot of the eigenvalues

# Component Loading: num_pred x num_pred matrix, i.e. correlation between num_pred and PCs
#	Loadings and eigenvectors contain similiar information (and may be equal, but are not the same)
# Use varimax to realign loadings, since there are multiple solutions to component loadings
#
# Component Score: num_obs x num_pred matrix
#	scale(pred_mat, center, scale) %*% loadings

# Compute PCA
eigen=eigen(resp_correl$val);

# Compute variance contribution of each PC
pca_propvar=eigen$values/sum(eigen$values);
pca_propcumsum=cumsum(pca_propvar);
num_pc_at_cutoff=sum(pca_propcumsum<PCCoverage)+1;

# Compute per sample scores
scores=(scale(resp_mat, center=T, scale=T) %*% eigen$vectors);

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
colors=rep("grey",num_resp);
colors[1:num_pc_at_cutoff]="darkcyan";

mids=barplot(pca_propvar, las=2, names.arg=1:num_resp, xlab="PCs", 
	col=colors,
	ylab="Proportion", main="PCA Proportion of Variance");

mids=barplot(pca_propcumsum, las=2, names.arg=1:num_resp, xlab="PCs", 
	col=colors,
	ylab="Proportion", main="PCA Cumulative Variance");
abline(h=PCCoverage, col="blue", lty=2);

# Plot sample ordination based on scores
sample_ids=rownames(resp_mat);
par(mfrow=c(1,1));

for(i in 1:(num_pc_at_cutoff-1)){
	xpos=scores[,i];
	ypos=scores[,i+1];
	xrange=range(xpos, na.rm=T);
	yrange=range(ypos, na.rm=T);

	xspan=diff(xrange);
	yspan=diff(yrange);

	plot(xpos, ypos, type="n", 
		xlim=c(xrange[1]-xspan/10, xrange[2]+xspan/10),
		ylim=c(yrange[1]-yspan/10, yrange[2]+yspan/10),
		xlab=paste("PC",i,sep=""), 
		ylab=paste("PC",i+1,sep=""), 
		main="PCA Ordination Plot");
	text(xpos, ypos, sample_ids);

}

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
