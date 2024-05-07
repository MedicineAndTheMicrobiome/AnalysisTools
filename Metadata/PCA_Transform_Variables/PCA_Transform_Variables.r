#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

options(useFancyQuotes=F);
options(digits=5)

params=c(
	"factors", "f", 1, "character",
	"outputroot", "o", 1, "character",
	"predictor", "p", 2, "character",
	"response", "r", 2, "character",
	"subsample", "s", 2, "numeric",
	"donnot_transform", "t", 2, "logical",
	"pc_coverage", "c", 2, "numeric",
	"pc_indiv_coverage", "i", 2, "numeric",
	"pc_var_correlation", "v", 2, "numeric",
	"export_orig", "O", 2, "logical",
	"export_curated", "C", 2, "logical",
	"export_imputed", "I", 2, "logical",
	"export_PC", "P", 2, "logical"
);

NORM_PVAL_CUTOFF=.2;
PCA_COVERAGE=.95;
PCA_INDIV_COVERAGE=0.025;
PCA_VAR_CORREL_CUTOFF=0.5;

CURATED_PREFIX="crtd";
NO_CHANGE="orig"

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

script_path=paste(head(strsplit(script_name, "/")[[1]], -1), collapse="/");
source(paste(script_path, "/Impute_Matrix.r", sep=""));

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <factors/metadata file name>\n",
	"	-o <output filename root.\n",
	"\n",
	"	[-p <targeted predictor variable list, default is all variables>]\n",
	"	[-s <subsample from predictors, default=use all samples>]\n",
	"	[-r <targeted response variable to plot against]\n",
	"\n",
	"	[-t (do not transform/autocurate predictors, default: transform if necessary)\n",
	"	[-c PC Coverage, default=", PCA_COVERAGE, "\n",
	"	[-i PC Individual coverage mininum, default=", PCA_INDIV_COVERAGE, "\n",
	"	[-v Min correlation at cutoff, default=", PCA_VAR_CORREL_CUTOFF, "\n",
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
	!length(opt$outputroot)
){
	cat(usage);
	q(status=-1);
}

FactorsFname=opt$factors;
OutputFnameRoot=opt$outputroot;

ResponseListName="";
PCCumulCoverage=PCA_COVERAGE;
PCIndivCoverage=PCA_INDIV_COVERAGE;
PCVariableCorrelationCutoff=PCA_VAR_CORREL_CUTOFF;
DonnotTransform=F;
SubsampleTargets=-1;

ExportOrig=F;
ExportCurated=F;
ExportImputed=F;
ExportPC=F;


if(length(opt$predictor)){
	PredictorListName=opt$predictor;
}else{
	PredictorListName="";
}	

if(length(opt$response)){
	ResponseListName=opt$response;
}

if(length(opt$pc_coverage)){
	PCCumulCoverage=opt$pc_coverage;
}

if(length(opt$pc_indiv_coverage)){
	PCIndivCoverage=opt$pc_indiv_coverage;
}

if(length(opt$pc_var_correlation)){
	PCVariableCorrelationCutoff=opt$pc_var_correlation;
}

if(length(opt$donnot_transform)){
	DonnotTransform=T;
}

if(length(opt$subsample)){
	SubsampleTargets=opt$subsample;
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
	cat("PC Min Cumulative Coverage: ", PCCumulCoverage, "\n");
	cat("PC Min Individual Cutoff: ", PCIndivCoverage, "\n");
	cat("PC Variable Correlation Cutoff: ", PCVariableCorrelationCutoff, "\n");
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

		if(heatmap_height*heatmap_width < 10000){
			layoutmat=matrix(
				rep(1, heatmap_height*heatmap_width),
				byrow=T, ncol=heatmap_width);
		}else{
			larger_dim=max(heatmap_height, heatmap_width);

			hm_height_norm=heatmap_height/larger_dim;	
			hm_width_norm=heatmap_width/larger_dim;

			lo_height=100*hm_height_norm;
			lo_width=100*hm_width_norm;

			layoutmat=matrix(
				rep(1, lo_height*lo_width),
				byrow=T, ncol=lo_width);

		}
			
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

	calc_guidelines=function(num_cells, rev=F){
		if(num_cells<8){
			return(NA);
		}else{
			grps=c(4, 5, 6, 7);
			rem=num_cells%%grps
			if(any(rem==0)){
				# Take largest of 0 remainders
				use_grp=grps[max(which(rem==0))];
			}else{
				# Take largest remainder
				grp_ix=max(which(rem==max(rem)));
				use_grp=grps[grp_ix]; 
			}

			if(rev){
				guide_pos=seq(0,num_cells,use_grp);
			}else{
				guide_pos=seq(num_cells,0,-use_grp);
			}
			guide_pos=setdiff(guide_pos, c(0, num_cells));
			return(guide_pos);
		}
	}

	abline(h=1:(num_row-1), lwd=.125, lty="dotted", col="grey95");
	abline(h=calc_guidelines(num_row), lty="dotted", lwd=2, col="white");
	abline(h=calc_guidelines(num_row), lty="dotted", lwd=1, col="black");

	abline(v=1:(num_col-1), lwd=.125, lty="dotted", col="grey95");
	abline(v=calc_guidelines(num_col, rev=T), lty="dotted", lwd=2, col="white");
	abline(v=calc_guidelines(num_col, rev=T), lty="dotted", lwd=1, col="black");

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

plot_title_page=function(title, subtitle="", title_cex=3){

        orig.par=par(no.readonly=T);
        par(family="serif");
        par(mfrow=c(1,1));

        plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n",  xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                );

        # Title
        title_line=1;
        text(0.5, title_line, title, cex=title_cex, font=2, adj=c(.5,1));

        # Subtitle
        num_subt_lines=length(subtitle);
        cxy=par()$cxy;
        for(i in 1:num_subt_lines){
                text(.5, title_line -title_cex*cxy[2] -i*cxy[2], subtitle[i], adj=.5);
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
if(PredictorListName!=""){
	predictors_arr=load_list(PredictorListName);
}else{
	predictors_arr=loaded_factor_names;
}

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

if(SubsampleTargets!=-1){
	if(SubsampleTargets<num_pred){
		predictors_arr=sample(predictors_arr, SubsampleTargets, replace=F);
		num_pred=SubsampleTargets;
		cat("Targets have been subsampled.\n");
	}
}

pred_mat=loaded_factors[, predictors_arr, drop=F];
resp_mat=loaded_factors[, responses_arr, drop=F];

orig_pred_names=predictors_arr;
orig_resp_names=responses_arr;

if(DonnotTransform){
	curated_pred_mat=pred_mat;
	curated_resp_mat=resp_mat;
}else{

	plot_title_page("Checking Variable Normality", c(
		"Variables are tested for normality with the Shapiro-Wilks test.",
		paste("If the p-value is < ", NORM_PVAL_CUTOFF, ", then both a log and sqrt", sep=""),
		"transform are applied and normality are retested.  If either transformation",
		"increases the normality (and increases the p-value), then the transformation",
		"that has the greater p-value is retained.",
		"",
		"In the following pages, the left column contains the variable's original",
		"distribution before any transformation is attempted.  The right column",
		"contains the accepted transformation if it is closer to a normal distribution.",
		"The p-values annotated below the variable names (titles) are calculated with",
		"the Shapiro-Wilks test.",
		"",
		"If a transformation was not applied, then the right column will hold the indication:",
		"\"Transform not necessary\", if the original distribution was sufficiently normal.",
		"\"Transform not beneficial\", if the transformation did not make the variable more normal."
	));

	cat("Testing Predictor Variables for normality.\n");
	curated_pred_mat=test_and_apply_log_transform(pred_mat, NORM_PVAL_CUTOFF);

	cat("Testing Response Variables for normality.\n");
	curated_resp_mat=test_and_apply_log_transform(resp_mat, NORM_PVAL_CUTOFF);
}

##############################################################################

curated_predictors_arr=colnames(curated_pred_mat);
curated_responses_arr=colnames(curated_resp_mat);

##############################################################################

if(length(responses_arr)){
	plot_text(c(
		"Relationship between predictors (x-axis) and responses (y-axis):",
		"",
		"Plots are sorted by decreasing statistical significance (increasing p-value)."
	));
}

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

##############################################################################

imputed_mat=impute_matrix(curated_pred_mat);
cat("Inputed Matrix Dimensions: ", nrow(imputed_mat), " x ", ncol(imputed_mat), "\n");

#print(imputed_mat);

##############################################################################

# Compute the dendrogram, clustering based on the correlation between responses and predictors

plot_title_page("Dendrogram of Variable Similarity", title_cex=2, c(
	"Variables with a higher degree of correlation (abs(cor)) are clustered",
	"together.  If 'predictor' variables were specified in this analyis, then",
	"they will be colored red."
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

##############################################################################

plot_title_page("Correlation Heatmaps", c(
	"The following heatmaps illustrate the correlation between variables at various",
	"cutoffs of statistical significance (both unadjusted and adjusted for multiple testing).",
	"",
	"Red is more positively correlated, Blue is more negatively correlated."
));

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
num_pc_at_cutoff=sum(pca_propcumsum<PCCumulCoverage)+1;

# Compute per sample scores
scores=(scale(imputed_mat, center=T, scale=T) %*% eigen$vectors);

# Top PC with coverate >5%
top_PCs=pca_propvar>PCIndivCoverage;
num_pc_above_indiv_cutoff=sum(top_PCs);

plot_title_page("PCA Results", c(
	"The following pages of results were based on running a PCA.",
	"",
	"Two criteria were provided as a means to select the number of PCs to",
	"represent a dataset.",
	"",
	paste("1.) Number of PCs to cumulatively acquire >", (PCCumulCoverage*100.0), 
		"% of the variance in dataset.", sep=""),
	"This criteria tries to capture a minimum amount of 'information' in a dataset.",
	"",
	paste("2.) Number of PCs individually with >", (PCIndivCoverage*100.0),
		"% of coverage.", sep=""),
	"This criteria tries to select only the PCs with a minimum amount of 'information'",
	"in them, with the assumption that PCs with low coverage may be representing 'noise'.",
	"",
	"",
	"Bar plots are also generated:",
	"1.) PCA Proportion of Variance",
	"The number of PCs that exceed the individual cutoff threshold are colored green.",
	"",
	"2.) PCA Cumulative Variance",
	"The number of PCs necessary to cumulative exceed the cumulative cutoff are colored teal."
));

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
	paste("Number of PCs to retain >", (PCCumulCoverage*100.0), "% of Variance:"),
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
abline(h=PCCumulCoverage, col="darkcyan", lty=2);

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
abline(h=PCCumulCoverage, col="darkcyan", lty=2);

###############################################################################

plot_title_page("PC Ordination", c(
	"These scatter plots illustrate the separation of samples across the top PCs.",
	"PC[i] and PC[i+1] are represented along the x and y axes, respectively.",
	"",
	"Left figures are labeled with sample IDs to help identify outliers, and the right",
	"figures only have glyphs drawn to more clearly see the distribution of samples",
	"across both axes/PCs."
));

# Plot sample ordination based on scores
sample_ids=rownames(resp_mat);
par(mfrow=c(3,2));
par(mar=c(3,3,1,1));

last_pc_to_plot=min(num_pc_at_cutoff+1, ncol(scores)-1);
for(i in seq(1, last_pc_to_plot, 2)){
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

top_pc_var_correl_header=c("PCID", "Correlation", "IndivCoverage", "CumulCoverage" );
top_pc_var_correl_mat=matrix(NA, nrow=num_pc_at_cutoff, ncol=length(top_pc_var_correl_header));
colnames(top_pc_var_correl_mat)=top_pc_var_correl_header;
top_pc_var_name=character();

par(mfrow=c(2,2));
par(mar=c(12,3,2,1));
positive_scores=scores;

cat("Num PCs at Cutoff: ", num_pc_at_cutoff, "\n");
cat("Num Scores: ", ncol(scores), "\n");

out_pc_cormat=matrix("", ncol=num_pc_at_cutoff*3, nrow=ncol(scores));
out_pc_cormat_header=rep("", num_pc_at_cutoff*3);

plot_title_page("Identifying PC Proxies", c(
	"Since PCs are composed of characteristics of the underlying observed variables in a",
	"dataset, often times it is more useful to identify a set of these observed variables",
	"to represent these PCs that can be more easily interpreted with domain-specific knowledge,",
	"rather than an abstract PC.",
	"",
	"In the following bar plots, a correlation is calculated between each PC and all the",
	"observed variables and the variables with the greatest magnitude of correlation",
	"with the PC are listed.  The PC is then annotated/named with the percent correlation",
	"and name of the observed variable.  If the relationship between the PC and the observed",
	"variable with the most correlation is negative, then the PC is flipped (multiplied by",
	"-1) so the correlation between the modified PC and the assigned observed variable is positive.",
	"",
	"The number of PCs annotated/plotted is the number of PCs to cumulatively cover the",
	"cumulative PC threshold.  The PCs with individual PCs exceeding the individual PC",
	"cutoff have their bars colored in green",
	"",
	"We have also observed that when a PC is not at least 50% correlated with one of the",
	"observed variables, then its ensuing utility is also potentially low (i.e. noise)."
));

fname=paste(OutputFnameRoot, ".proxy_top_correl.tsv", sep="");
fh=file(fname, "w");

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

	# Save top variable correlated to PC
	top_pc_var_name=c(top_pc_var_name, names(pc_pred_cor_ordered)[1]) ;
	top_pc_var_correl_mat[i,]=c(i, pc_pred_cor_ordered[1], pca_propvar[i], pca_propcumsum[i]);
	
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

	# Write proxy top correlates to file
	cat(file=fh, "ProxyName\tTopCorrelates\tCorrel_wProxy\n");
	cat(file=fh, proxyname, "\t\n", sep="");
	cat(file=fh, paste("",
		names(pc_pred_cor_ordered)[1:num_corr_bars_to_plot], 
		sprintf("%3.4f", pc_pred_cor_ordered[1:num_corr_bars_to_plot]),
		sep="\t"), sep="\n");
	cat(file=fh, "\t\n");
}

close(fh);

rownames(top_pc_var_correl_mat)=top_pc_var_name;

###############################################################################
# Select variables based on the difference cutoff options

plot_title_page("Selected Variables", c(
	"The following pages report the PCA-based proxy observed variable selection",
	"results based on the 3 previously described criteria.",
	"",
	"The next page reports the selected observed variables as well as a table of the",
	"criteria used to select them including: correlations, individual and cumulative",
	"coverage."
));	

selected_variables=list();

abv_cor_max_ix=max(which(top_pc_var_correl_mat[,"Correlation"]>=PCVariableCorrelationCutoff));
abv_ind_max_ix=max(which(top_pc_var_correl_mat[,"IndivCoverage"]>=PCIndivCoverage));

selected_variables[["by_correl"]]=unique(top_pc_var_name[1:abv_cor_max_ix]);
selected_variables[["by_indiv"]]=unique(top_pc_var_name[1:abv_ind_max_ix]);
selected_variables[["by_cumul"]]=unique(top_pc_var_name);

plot_text(c(
	"Selected Variables By Cutoff:",
	"",
	capture.output(print(top_pc_var_correl_mat)),
	"",
	paste("Correlation Cutoff:", PCVariableCorrelationCutoff),
	selected_variables[["by_correl"]],
	"",
	paste("Individual Coverage Cutoff:", PCIndivCoverage),
	selected_variables[["by_indiv"]],
	"",
	paste("Cumulative Coverage Cutoff:", PCCumulCoverage),
	selected_variables[["by_cumul"]]
));

colnames(out_pc_cormat)=out_pc_cormat_header;
rownames(out_pc_cormat)=paste(1:ncol(scores), ".)", sep="");

cat("Outputing PC to Variable Correlations:\n");
fname=paste(OutputFnameRoot, ".pc_var_cor.tsv", sep="");
write.table(out_pc_cormat, file=fname, col.names=NA, append=T, quote=F, sep="\t");

##############################################################################

# Plot dendrogram with selected variable

find_height_at_k=function(hclust, k){
# Computes the height on the dendrogram for a particular k

	heights=hclust$height;
	num_heights=length(heights);
	num_clust=numeric(num_heights);
	for(i in 1:num_heights){
		num_clust[i]=length(unique(cutree(hclust, h=heights[i])));
	}
	height_idx=which(num_clust==k);
	midpoint=(heights[height_idx+1]+heights[height_idx])/2;
	return(midpoint);
}

var_cor=compute_correlations(curated_pred_mat);
hcl=hclust(var_cor$dist, method="ward.D2");
dend=as.dendrogram(hcl);

plot_title_page("Selected Variable Dendrograms", c(
	"The following 3 dendrograms illustrate the relative location of selected",
	"observed variables in the context of all observed variables.",
	"",
	"Selected variables are colored teal.",
	"",
	"A dotted line is drawn based on the number of variables selected.",
	"The cutree algorithm will identify the height at which to cut a hierarchically",
	"cluster dendrogram so that a number of clusters (m) can be specified.  The number of",
	"clusters (m) specified in each figure is based on the number of variables selected.",
	"Ideally, if each of the clusters has exactly one selected variable in it,",
	"then the selected variables were able to represent the entire dataset reasonably",
	"efficiently."
));

par(mar=c(14, 2, 4, 1));
par(mfrow=c(1,1));
for(sel_type in c("by_correl", "by_indiv", "by_cumul")){

	cat("Generating Dendrogram for: ", sel_type, "\n");
	cur_sel_var=selected_variables[[sel_type]];

	highlight_vars=function(x){
		if(is.leaf(x)){
			leaf_attr=attributes(x);
			label=leaf_attr$label;
			#print(label);
			if(any(label==cur_sel_var)){
				cat("Coloring: ", label, "\n");
				color="darkcyan";
				font=2;
				cex=1.05;
			}else{
				color="black";
				font=1;
				cex=.95;
			}
			attr(x, "nodePar")=c(leaf_attr$nodePar, 
				list(lab.font=font, lab.col=color, lab.cex=cex, cex=0));
		}
		return(x);
	}

	
	cur_dend=dendrapply(dend, highlight_vars);

	if(sel_type == "by_correl"){
		msg=paste("Selected Variables with Correlation with PCs  > ", 
			PCVariableCorrelationCutoff, sep="");
	}else if(sel_type == "by_indiv"){
		msg=paste("Selected Variables Most Similar with PCs Individually Covering > ", 
			PCIndivCoverage, sep="");
	}else if(sel_type == "by_cumul"){
		msg=paste("Selected Variables Most Similar to PCs Cumulatively Covering > ", 
			PCCumulCoverage, sep="");
	}

	plot(cur_dend, , main=paste(msg, "\nWard's Minimum Variance:\ndist(1-abs(cor))", sep=""));


	effect_cut=find_height_at_k(hcl, k=length(cur_sel_var));
	abline(h=effect_cut, col="blue", lty="dashed");


}

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

plot_title_page("Dendrogram w/ PCs", c(
	"In this dendrogram, the PCs (that have been annotated/renamed with the observed variable",
	"they are most similar (correlated) with have been included in the dendrogram with all",
	"observed variables.  The green PCs are the variables that have been identified as capturing",
	"more than the minimum individual variance threshold.  The remaining teal PCs cumulatively exceed",
	"the min cumulative variance threshold."
));

par(mfrow=c(1,1));
par(mar=c(2,1,4,20));
plot(dend, horiz=T, main="Relationship of PCs with Original Variables\nWard's Minimum Variance:\ndist(1-abs(cor)) with PCs");


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
# Export selected variables

ext_list=list();
ext_list[["by_correl"]]=sprintf("var_corr%02.1f", PCVariableCorrelationCutoff*100);
ext_list[["by_indiv"]]= sprintf("pc_indiv%02.1f", PCIndivCoverage*100);
ext_list[["by_cumul"]]= sprintf("pc_cumul%02.1f", PCCumulCoverage*100);

for(criteria in c("by_correl", "by_indiv", "by_cumul")){
	sel_mat=curated_pred_mat[,selected_variables[[criteria]]];

	fname=paste(OutputFnameRoot, ".selected_transf_var.", ext_list[[criteria]], ".tsv", sep="");
	fh=file(fname, "w");
	cat(file=fh, "SampleID");
	close(fh);
	write.table(sel_mat, file=fname, col.names=NA, append=T, quote=F, sep="\t");	
}


##############################################################################

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
