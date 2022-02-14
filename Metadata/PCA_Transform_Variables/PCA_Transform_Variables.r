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
	"donnot_transform", "t", 2, "logical",
	"pc_coverage", "c", 2, "numeric",
	"pc_indiv_coverage", "i", 2, "numeric",
	"export_orig", "O", 2, "logical",
	"export_curated", "C", 2, "logical",
	"export_imputed", "I", 2, "logical",
	"export_PC", "P", 2, "logical"
);

NORM_PVAL_CUTOFF=.2;
PCA_COVERAGE=.95;
PCA_INDIV_COVERAGE=0.05;
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
	"	[-t (do not transform/autocurate predictors, default: transform if necessary)\n",
	"	[-c PC Coverage, default=", PCA_COVERAGE, "\n",
	"	[-i PC Individual coverage mininum, default=", PCA_INDIV_COVERAGE, "\n",
	"\n",
	"Output Options:\n",
	"	[-O (export Original predictor variable list)\n",
	"	[-C (export 'Curated', i.e. log and orig values depending on Shapiro-Wilkes)\n",
	"	[-I (export Imputed values)\n",
	"	[-P (export Principal Components)\n",
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
PredictorListName=opt$predictor;

ResponseListName="";
PCCoverage=PCA_COVERAGE;
PCIndivCoverage=PCA_INDIV_COVERAGE;
DonnotTransform=F;

ExportOrig=F;
ExportCurated=F;
ExportImputed=F;
ExportPC=F;

if(length(opt$response)){
	ResponseListName=opt$response;
}

if(length(opt$pc_coverage)){
	PCCoverage=opt$pc_coverage;
}

if(length(opt$pc_indiv_coverage)){
	PCIndivCoverage=opt$pc_indiv_coverage;
}

if(length(opt$donnot_transform)){
	DonnotTransform=T;
}

if(length(opt$export_orig)){
	ExportOrig=T;
}
if(length(opt$export_curated)){
	ExportCurated=T;
}
if(length(opt$export_imputed)){
	ExportImputed=T;
}
if(length(opt$export_PC)){
	ExportPC=T;
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

print(param_text, quote=F);

###############################################################################

load_factors=function(fname){
	factors=as.data.frame(read.delim(fname,  header=TRUE, row.names=1, check.names=FALSE, sep="\t"));
	return(factors);
}

load_list=function(fname){
	cat("Loading: ", fname, "\n");
	lst=read.delim(fname, header=F, check.names=F, comment.char="#", as.is=T);
	return(lst[,1]);	
}

test_and_apply_log_transform=function(mat_val, pval_cutoff=.2, plot_before_after=T){
	nrows=nrow(mat_val);
	ncols=ncol(mat_val);

	trans_mat=mat_val;
	orig_colnames=colnames(mat_val);
	new_colnames=character();

	if(plot_before_after){
		orig_par=par(no.readonly=T);
		par(mfrow=c(5,2));
	}

	delete_list=c();
	for(var in orig_colnames){
		values=mat_val[,var];


		log_transformed=F;
		sqrt_transformed=F;

		num_unique_val=length(setdiff(unique(values), NA));
		values_nona=values[!is.na(values)];
		num_nona=length(values_nona);

		if(!is.numeric(values_nona)){
			cat("Error: Values not numeric for: ", var, "\n", sep="");
			print(values_nona);
		}

		if(num_nona<=3){
			cat("Not enough non NA values to measure normality.\n");
			new_colnames=c(new_colnames, var);
			next;
		}

		if(any(values_nona<0)){
			cat("Negative values.  Skipping transformation.\n");
			new_colnames=c(new_colnames, var);
			next;
		}

		cat("\n", var, ": Num Unique Values: ", num_unique_val, "\n");

		if(num_unique_val>1){

			test_res=shapiro.test(values);
			test_log_res=NULL;

			if(test_res$p.value<=pval_cutoff && num_unique_val>2){
				cat(" Not normal: ", test_res$p.value, "\n");

				if(any(values_nona==0)){
					log_values=log(values+1);
				}else{
					log_values=log(values);
				}
				sqrt_values=sqrt(values);

				test_log_res=shapiro.test(log_values);
				test_sqrt_res=shapiro.test(sqrt_values);

				if(test_log_res$p.value < test_res$p.value && 
					test_sqrt_res$p.value < test_res$p.value){
					# Keep original
					cat("  No Improvement: ", test_log_res$p.value, "\n");
					new_colnames=c(new_colnames, paste("orig_", var, sep=""));
				}else{

					cat("     Log p-val : ", test_log_res$p.value, "\n");
					cat("    Sqrt p-val : ", test_sqrt_res$p.value, "\n");
					
					if(test_log_res$p.value > test_sqrt_res$p.value){
						# Keep log transformed
						cat("  Log Transformation Effective: ", test_log_res$p.value, "\n");
						new_colnames=c(new_colnames, paste("log_", var, sep=""));
						trans_mat[, var]=log_values;
						log_transformed=T;		
					}else{
						# Keep sqrt transformed
						cat("  Sqrt Transformation Effective: ", test_sqrt_res$p.value, "\n");
						new_colnames=c(new_colnames, paste("sqrt_", var, sep=""));
						trans_mat[, var]=sqrt_values;
						sqrt_transformed=T;		
					}
				}
			}else{
				cat(" Normal enough. ", test_res$p.value, "\n");
				new_colnames=c(new_colnames, var);
			}

		}else{
			cat("  All values identical, removing...\n");
			new_name=paste("all_ident_", var, sep="");
			new_colnames=c(new_colnames, new_name);
			delete_list=c(delete_list, new_name);
		}

		if(plot_before_after){
			nclass=nclass.Sturges(values)*4;

			hist(values, main=var, breaks=nclass);
			title(main=sprintf("p-value: %4.4f", test_res$p.value), cex.main=.8, line=.5);
			
			if(log_transformed){
				hist(log_values, breaks=nclass, main=paste("log(", var,")", sep=""));
				title(main=sprintf("p-value: %4.4f", test_log_res$p.value), cex.main=.8, line=.5);
			}else if(sqrt_transformed){
				hist(sqrt_values, breaks=nclass, main=paste("sqrt(", var,")", sep=""));
				title(main=sprintf("p-value: %4.4f", test_sqrt_res$p.value), cex.main=.8, line=.5);
			}else{
				plot(0,0, xlab="", ylab="", main="", xaxt="n", yaxt="n", bty="n", type="n");

				if(test_res$p.value>pval_cutoff){
					text(0,0, "Transform not necessary");
				}else{
					text(0,0, "Transform not beneficial");
				}
			}
		}

	}

	colnames(trans_mat)=new_colnames;

	trans_mat=trans_mat[,setdiff(new_colnames, delete_list),drop=F];

	if(plot_before_after){
		par(orig_par);
	}

	

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

	cat("Computing Correlations: \n");
	num_col=ncol(mat);
	cat("Num Columns: ", num_col, "\n");

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

			if(!is.na(test$estimate)){
				pval_mat[i,j]=test$p.value;
				pval_mat[j,i]=test$p.value;
				cor_mat[i,j]=test$estimate;
				cor_mat[j,i]=test$estimate;
			}else{
				pval_mat[i,j]=1;
				pval_mat[j,i]=1;
				cor_mat[i,j]=0;
				cor_mat[j,i]=0;
			}
		}
	}
	res=list();
	res[["val"]]=cor_mat;
	res[["pval"]]=pval_mat;;
	res[["dist"]]=as.dist(1-abs(cor_mat));

	return(res);
}

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, deci_pts=4,
        label_zeros=T, counts=F, value.cex=1,
        plot_col_dendr=F,
        plot_row_dendr=F
){

        num_row=nrow(mat);
        num_col=ncol(mat);

        row_names=rownames(mat);
        col_names=colnames(mat);

        orig.par=par(no.readonly=T);

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
                        dendist=dist(in_mat);
                }else{
                        dendist=dist(t(in_mat));
                }

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
                        c(
                        rep(rep(2, heatmap_width), col_dend_height),
                        rep(rep(1, heatmap_width), heatmap_height)
                        ), byrow=T, ncol=heatmap_width);

                col_dendr=get_dendrogram(mat, type="col");
                mat=mat[, col_dendr[["names"]]];

        }else if(plot_row_dendr){
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
        plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), xaxt="n", yaxt="n", bty="n", xlab="", ylab="");
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
                                text(x-.5, y-.5, text_lab, srt=atan(num_col/num_row)/pi*180, cex=value.cex, font=2);
                        }
                }
        }

        ##################################################################################################

        par(mar=c(0, 0, 0, 0));

        if(plot_row_dendr && plot_col_dendr){
                rdh=attributes(row_dendr[["tree"]])$height;
                cdh=attributes(col_dendr[["tree"]])$height;
                plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", bty="n", xlim=c(rdh, 0));
                plot(col_dendr[["tree"]], leaflab="none",xaxt="n", yaxt="n", bty="n", ylim=c(0, cdh));
                plot(0,0, type="n", bty="n", xaxt="n", yaxt="n");
                #text(0,0, "Placeholder");
        }else if(plot_row_dendr){
                rdh=attributes(row_dendr[["tree"]])$height;
                plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", bty="n", xlim=c(rdh, 0));
                #text(0,0, "Row Dendrogram");
        }else if(plot_col_dendr){
                cdh=attributes(col_dendr[["tree"]])$height;
                plot(col_dendr[["tree"]], leaflab="none", xaxt="n", yaxt="n", bty="n", ylim=c(0, cdh));
                #text(0,0, "Col Dendrogram");
        }

        par(orig.par);

}

##############################################################################
# Main Program Starts Here!
##############################################################################

pdf(paste(OutputFnameRoot, ".pca.pdf", sep=""), height=11, width=8.5);

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

if(DonnotTransform){
	curated_pred_mat=pred_mat;
	curated_resp_mat=resp_mat;
}else{
	cat("Testing Predictor Variables for normality.\n");
	curated_pred_mat=test_and_apply_log_transform(pred_mat, NORM_PVAL_CUTOFF);

	cat("Testing Response Variables for normality.\n");
	curated_resp_mat=test_and_apply_log_transform(resp_mat, NORM_PVAL_CUTOFF);
}

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
		if(length(intersect("pred_val", rownames(sumfit$coefficients)))){
			pval_arr[i]=sumfit$coefficients["pred_val", "Pr(>|t|)"];
			coef_arr[i]=sumfit$coefficients["pred_val", "Estimate"];
			rsqd_arr[i]=sumfit$r.squared;
			fit_arr[[i]]=fit;
		}else{
			pval_arr[i]=1;
			coef_arr[i]=0;
			rsqd_arr[i]=0;
			fit_arr[[i]]=NA;
		}
	
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
		if(!is.na(fit_arr[[i]])){
			abline(fit_arr[[i]], col="blue");
		}
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

mask_matrix=function(val_mat, mask_mat, mask_thres, mask_val){
        masked_matrix=val_mat;
        masked_matrix[mask_mat>mask_thres]=mask_val;
        return(masked_matrix);
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

	if(length(na_pos)==0){
		num_nas_to_impute=0;
	}else{
		num_nas_to_impute=nrow(na_pos);
	}

	cat("Num NAs to try to impute:",  num_nas_to_impute, "\n");

	filled_matrix=usable_mat_wna;

	if(!is.null(num_nas_to_impute) && num_nas_to_impute>0){
		for(na_ix in 1:num_nas_to_impute){

			target_row=na_pos[na_ix,1];
			target_column=na_pos[na_ix,2];

			cell=usable_mat_wna[target_row, target_column, drop=F];

			if(!is.na(cell)){
				cat("Error:  Trying to input cell not NA.\n");
				quit();
			}
		
			cat("(", na_ix, "/", num_nas_to_impute, ") Imputing: ", 
				rownames(cell), " / ", colnames(cell), "\n");

			non_na_row=!is.na(usable_mat_wna[,target_column,drop=F]);
			non_na_col=!is.na(usable_mat_wna[target_row,,drop=F]);

			imputed_val=impute_cell(
				target_predictors=usable_mat_wna[target_row, non_na_col, drop=F], 
				responses=usable_mat_wna[non_na_row, target_column, drop=F],
				predictors=usable_mat_wna[non_na_row, non_na_col, drop=F]
				);

			filled_matrix[target_row, target_column]=imputed_val;

		}
	}

	repl_rows=rownames(filled_matrix);
	mat_wna[repl_rows,]=filled_matrix[repl_rows,];
	return(mat_wna);
	
}

imputed_mat=impute_matrix(curated_pred_mat);
cat("Inputed Matrix Dimensions: ", nrow(imputed_mat), " x ", ncol(imputed_mat), "\n");

#print(imputed_mat);

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
			font=1;
		}else{
			color="red";
			font=2;
		}
		attr(x, "nodePar")=c(leaf_attr$nodePar, list(lab.font=font, lab.col=color, cex=0));
	}
	return(x);
}

dend=dendrapply(dend, highlight_predictors);

plot(dend, main="Ward's Minimum Variance: dist(1-abs(cor))");

print(correl$val);
paint_matrix(correl$val, deci_pts=2, title="All Correlations");

signf_10=mask_matrix(correl$val, correl$pval, 0.1, mask_val=0);
paint_matrix(signf_10, deci_pts=2, label_zeros=F, title="Correlations (p-val<0.10)");

signf_05=mask_matrix(correl$val, correl$pval, 0.05, mask_val=0);
paint_matrix(signf_05, deci_pts=2, label_zeros=F, title="Correlations (p-val<0.05)");

signf_01=mask_matrix(correl$val, correl$pval, 0.01, mask_val=0);
paint_matrix(signf_01, deci_pts=2, label_zeros=F, title="Correlations (p-val<0.01)");

num_comparisons=nrow(correl$val)*(nrow(correl$val)-1)/2;
cat("Num comparisons to correct for: ", num_comparisons, "\n");

signf_05_bonf=mask_matrix(correl$val, correl$pval, 0.05/num_comparisons, mask_val=0);
paint_matrix(signf_05_bonf, deci_pts=2, label_zeros=F, title="Correlations (Bonferroni Corrected: p-val<0.05)");

signf_01_bonf=mask_matrix(correl$val, correl$pval, 0.01/num_comparisons, mask_val=0);
paint_matrix(signf_01_bonf, deci_pts=2, label_zeros=F, title="Correlations (Bonferroni Corrected: p-val<0.01)");

signf_001_bonf=mask_matrix(correl$val, correl$pval, 0.001/num_comparisons, mask_val=0);
paint_matrix(signf_001_bonf, deci_pts=2, label_zeros=F, title="Correlations (Bonferroni Corrected: p-val<0.001)");

##############################################################################

pred_correl=compute_correlations(curated_pred_mat);
cat("Predictor Correlation Matrix Dimensions: ", nrow(pred_correl$val), " x ", ncol(pred_correl$val), "\n");

# Note that the princomp R function takes the resp values directly and the Standard Deviations
# are the squareroot of the eigenvalues

# Component Loading: num_pred x num_pred matrix, i.e. correlation between num_pred and PCs
#	Loadings and eigenvectors contain similiar information (and may be equal, but are not the same)
# Use varimax to realign loadings, since there are multiple solutions to component loadings
#
# Component Score: num_obs x num_pred matrix
#	scale(pred_mat, center, scale) %*% loadings

# Compute PCA
cat("Calculating Eigen Values...\n");
eigen=eigen(pred_correl$val);

# Compute variance contribution of each PC
pca_propvar=eigen$values/sum(eigen$values);
pca_propcumsum=cumsum(pca_propvar);
num_pc_at_cutoff=sum(pca_propcumsum<PCCoverage)+1;

# Compute per sample scores
scores=(scale(imputed_mat, center=T, scale=T) %*% eigen$vectors);

# Top PC with coverate >5%
top_PCs=pca_propvar>PCIndivCoverage;
num_pc_above_indiv_cutoff=sum(top_PCs);

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
	num_pc_at_cutoff,
	"",
	paste("Number of PCs with each with >", (PCIndivCoverage*100.0), "% of Variance:"),
	num_pc_above_indiv_cutoff

));

# Plot bar plots of PC variance explanation
num_kept_pred=length(pca_propvar);
par(mfrow=c(2,1));
par(mar=c(7,4,2,2));

###############################################################################

cat("Generating PCA variance barplots...\n");

# Individual
colors=rep("grey",num_kept_pred);
colors[1:num_pc_above_indiv_cutoff]="darkgreen";
mids=barplot(pca_propvar, las=2, names.arg=1:num_kept_pred, xlab="PCs", 
	col=colors,
	ylab="Proportion", main="PCA Proportion of Variance");
abline(h=PCIndivCoverage, col="darkgreen", lty=2);

# Cumulative
colors=rep("grey",num_kept_pred);
colors[1:num_pc_at_cutoff]="darkcyan";
mids=barplot(pca_propcumsum, las=2, names.arg=1:num_kept_pred, xlab="PCs", 
	col=colors,
	ylab="Proportion", main="PCA Cumulative Variance");
abline(h=PCCoverage, col="darkcyan", lty=2);

###############################################################################

cat("Generating PCA variance barplots (zoomed in)...\n");

# Individual
colors=rep("grey",num_kept_pred);
colors[1:num_pc_above_indiv_cutoff]="darkgreen";
mids=barplot(pca_propvar, las=2, names.arg=1:num_kept_pred, xlab="PCs", 
	col=colors,
	ylab="Proportion", main="PCA Proportion of Variance\n(zoomed)",
	xlim=c(0, num_pc_at_cutoff+2)
	);
abline(h=PCIndivCoverage, col="darkgreen", lty=2);

# Cumulative
colors=rep("grey",num_kept_pred);
colors[1:num_pc_at_cutoff]="darkcyan";
mids=barplot(pca_propcumsum, las=2, names.arg=1:num_kept_pred, xlab="PCs", 
	col=colors,
	ylab="Proportion", main="PCA Cumulative Variance\n(zoomed)",
	xlim=c(0, num_pc_at_cutoff+2)
	);
abline(h=PCCoverage, col="darkcyan", lty=2);

###############################################################################

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
		xlab="", ylab="", main="", cex=2
	)

	title(
		xlab=paste("PC",i,sep=""), 
		ylab=paste("PC",i+1,sep=""), 
		line=2
	);

}

##############################################################################

# Calculate correlation between each PC and the predictors
cat("Calculating correlation between each PCA and the predictors...\n");

par(mfrow=c(2,2));
par(mar=c(12,3,2,1));
positive_scores=scores;

cat("Num PCs at Cutoff: ", num_pc_at_cutoff, "\n");
cat("Num Scores: ", ncol(scores), "\n");

out_pc_cormat=matrix("", ncol=num_pc_at_cutoff*3, nrow=ncol(scores));
out_pc_cormat_header=rep("", num_pc_at_cutoff*3);

pc_name=paste("PC", sprintf("%02g", 1:num_kept_pred), sep="");
for(i in 1:num_pc_at_cutoff){

	cat("Working on: PC ", i, "\n");
	
	pc=scores[,i];
	nonna_pc=!is.na(pc);

	pc_pred_cor=numeric(num_kept_pred);
	names(pc_pred_cor)=colnames(curated_pred_mat)
	for(pix in 1:num_kept_pred){

		prd=curated_pred_mat[,pix];
		nonna_prd=!is.na(prd);
		both_nonna=(nonna_pc & nonna_prd);

		pc_pred_cor[pix]=cor(pc[both_nonna], prd[both_nonna]);

	}

	
	mag_order=order(abs(pc_pred_cor), decreasing=T);
	pc_pred_cor_ordered=pc_pred_cor[mag_order];

	# If the most correlated is negative, flip the PC values
	flipped="";
	if(pc_pred_cor_ordered[1]<0){
		positive_scores[,i]=positive_scores[,i]*-1;
		pc_pred_cor_ordered=pc_pred_cor_ordered*-1
		flipped=" (flipped)";
	}
	
	ordered_names=names(pc_pred_cor_ordered);
	
	proxyname=paste(
		pc_name[i], "_",
		round(pc_pred_cor_ordered[1]*100, 0), "_", 
		ordered_names[1], sep="");

	# Max bars to plot
	num_corr_bars_to_plot=min(20, length(pc_pred_cor_ordered));

	barplot(pc_pred_cor_ordered[1:num_corr_bars_to_plot], 
		names.arg=names(pc_pred_cor_ordered)[1:num_corr_bars_to_plot], 
		col=ifelse(i<=num_pc_above_indiv_cutoff, "darkgreen", "grey"),
		ylim=c(-1,1),
		ylab="Correlation",
		las=2, cex.names=.7,
		main=paste(proxyname, flipped, sep="")
		);

	pc_name[i]=proxyname;

	# Save to output matrix for export
	out_pc_cormat[, ((i-1)*3)+1]=names(pc_pred_cor_ordered);
	out_pc_cormat[, ((i-1)*3)+2]=round(pc_pred_cor_ordered, 4);
	out_pc_cormat_header[((i-1)*3)+1]=proxyname;
	out_pc_cormat_header[((i-1)*3)+2]="Correlation";
}

colnames(out_pc_cormat)=out_pc_cormat_header;
rownames(out_pc_cormat)=paste(1:ncol(scores), ".)", sep="");

cat("Outputing PC to Variable Correlations:\n");
fname=paste(OutputFnameRoot, ".pc_var_cor.tsv", sep="");
write.table(out_pc_cormat, file=fname, col.names=NA, append=T, quote=F, sep="\t");

##############################################################################

colnames(positive_scores)=pc_name;

correl_wpca=compute_correlations(cbind(positive_scores[,1:num_pc_at_cutoff], curated_pred_mat, curated_resp_mat));
hcl=hclust(correl_wpca$dist, method="ward.D2");
dend=as.dendrogram(hcl);

highlight_pcs=function(x){
	if(is.leaf(x)){
		leaf_attr=attributes(x);
		label=leaf_attr$label;
		print(label);
		if(any(label==pc_name)){
			if(any(label==pc_name[1:num_pc_above_indiv_cutoff])){
				color="darkgreen";
			}else{
				color="darkcyan";
			}
			font=2;
			cex=1.05;
		}else{
			color="grey33";
			font=1;
			cex=1/1.05;
		}
		attr(x, "nodePar")=c(leaf_attr$nodePar, list(lab.font=font, lab.col=color, lab.cex=cex, cex=0));
	}
	return(x);
}

dend=dendrapply(dend, highlight_pcs);

par(mfrow=c(1,1));
par(mar=c(2,1,1,20));
plot(dend, horiz=T, main="Ward's Minimum Variance: dist(1-abs(cor)) with PCs");


##############################################################################

out_factors=loaded_factors;

if(!ExportOrig){
	cn=colnames(out_factors);
	orig_names=c(colnames(curated_pred_mat), colnames(curated_resp_mat));
	kept=setdiff(cn, orig_names);
	out_factors=out_factors[,kept];
}

append_columns=function(original_mat, additional_mat){

	origmat_dim=dim(original_mat);
	addmat_dim=dim(additional_mat);

	cat("Inserting: ", addmat_dim[1], "x", addmat_dim[2], " into ", origmat_dim[1], "x", origmat_dim[2], "\n");

	orig_cnames=colnames(original_mat);
	add_cnames=colnames(additional_mat);	
	samp_ids=rownames(original_mat);
	add_ids=rownames(additional_mat);

	comb_mat=as.data.frame(matrix(NA, nrow=origmat_dim[1], ncol=origmat_dim[2]+addmat_dim[2]));

	rownames(comb_mat)=samp_ids;
	colnames(comb_mat)=c(orig_cnames, add_cnames);

	#comb_mat[samp_ids, orig_cnames]=original_mat[samp_ids, orig_cnames];

	# Copy original mat over
	for(cnames in orig_cnames){
		comb_mat[,cnames]=original_mat[,cnames];
	}

	# Copy addition mat over
	avail_ids=intersect(samp_ids, add_ids);
	for(cnames in add_cnames){
		comb_mat[avail_ids, cnames]=additional_mat[avail_ids, cnames];
	}

	dim_return=dim(comb_mat);
	cat("Returning Matrix: ", dim_return[1], " x ", dim_return[2], "\n"); 
	return(comb_mat);
}


if(ExportCurated){
	cat("Appending Curated Predictors...\n");
	out_factors=append_columns(out_factors, curated_pred_mat)	
	cat("Appending Curated Responders...\n");
	out_factors=append_columns(out_factors, curated_resp_mat)	
}

if(ExportImputed){
	cat("Appending Imputed Variables\n");
	colnames(imputed_mat)=paste("imp.", colnames(imputed_mat), sep="");
	out_factors=append_columns(out_factors, imputed_mat);
}

if(ExportPC){
	cat("Appending Computed PCs\n");
	out_factors=append_columns(out_factors, positive_scores[,1:num_pc_at_cutoff, drop=F]);
}

##############################################################################

cat("Outputing New Factor File Values:\n");
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
