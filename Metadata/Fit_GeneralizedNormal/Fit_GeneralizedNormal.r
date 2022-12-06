#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
library('gnorm'); # for gnrom functions
library('stats4'); # For mle

options(useFancyQuotes=F);
options(digits=5)

INDIV_PCA_CUTOFF=0.05;
CUMUL_PCA_CUTOFF=0.80;
NUM_PC_GLOBAL=5;

params=c(
	"factors", "f", 1, "character",
	"outputroot", "o", 1, "character",
	"beta_export_cutoff", "b", 2, "numeric",
	"debug_targets", "d", 2, "character",
	"fast_fit", "F", 2, "logical",
	"indiv_pca_cutoff", "i", 2, "numeric",
	"cumul_pca_cutoff", "c", 2, "numeric",
	"num_pcs_glob_cor", "p", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

#script_path=paste(head(strsplit(script_name, "/")[[1]], -1), collapse="/");
#source(paste(script_path, "/Impute_Matrix.r", sep=""));

NORM_PVAL_CUTOFF=.2;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <factors/metadata file name>\n",
	"	-o <output filename root.\n",
	"	[-b <Beta export cutoff for variable export>]\n",
	"\n",
	"	[-d <debug list of variable targets>]\n",
	"	[-F (fast fit)>]\n",
	"\n",
	"	[-i <Individual PC Cutoff, default=", INDIV_PCA_CUTOFF, ">]\n",
	"	[-c <Cumulative PC Cutoff, default=", CUMUL_PCA_CUTOFF, ">]\n",
	"\n",
	"	[-p <Num PCs for global correlation, default=", NUM_PC_GLOBAL, ">]\n",
	"\n",
	"This script will try to sort the variables based on\n",
	"how much 'information' the variables have.  Since\n",
	"one can't strictly compare entropy between continuous\n",
	"distributions, and the variance of a variable is\n",
	"dependent on scale, we fit the generalized normal distribution\n",
	"and then preferentially analyze variables with a larger\n",
	"shape parameters (flatter).\n",
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

DebugVarList="";
if(length(opt$debug_targets)){
	DebugVarList=opt$debug_targets;	
}

BetaExportCutoff=NA;
if(length(opt$BetaExportCutoff)){
	BetaExportCutoff=opt$BetaExportCutoff;
	OutputFnameRoot=paste(OutputFnameRoot, ".b", BetaExportCutoff, sep="");
}

FastFit=F;
if(length(opt$fast_fit)){
	FastFit=T;
}

if(FastFit){
	OutputFnameRoot=paste(OutputFnameRoot, ".fast", sep="");
}

Indiv_PC_Cutoff=INDIV_PCA_CUTOFF;
Cumul_PC_Cutoff=CUMUL_PCA_CUTOFF;

if(length(opt$indiv_pca_cutoff)){
	Indiv_PC_Cutoff=opt$indiv_pca_cutoff;
}

if(length(opt$cumul_pca_cutoff)){
	Cumul_PC_Cutoff=opt$cumul_pca_cutoff;
}

Num_PC_Global_Correl=NUM_PC_GLOBAL;
if(length(opt$num_pcs_glob_cor)){
	Num_PC_Global_Correl=opt$num_pcs_glob_cor;
}

param_text=capture.output({
	cat("\n");
	cat("Factor File Name: ", FactorsFname, "\n");
	cat("Output File Name Root: ", OutputFnameRoot, "\n");
	cat("Fast Fit? ", FastFit, "\n");
	cat("\n");
});

print(param_text, quote=F);

###############################################################################

load_factors=function(fname){
	factors=as.data.frame(read.delim(fname,  header=TRUE, row.names=1, check.names=FALSE, sep="\t"));
	factors=factors[,-1,drop=F];

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

		if(heatmap_height < 200 &&  heatmap_width < 200){
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

pdf(paste(OutputFnameRoot, ".fit_gnd.raw.pdf", sep=""), height=11, width=8.5);

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

targeted_mat=loaded_factors;

cat("Testing Predictor Variables for normality.\n");
curated_targets_mat=test_and_apply_log_transform(targeted_mat, NORM_PVAL_CUTOFF);

if(DebugVarList!=""){
	debug_tar_list=scan(DebugVarList, what=character(), comment.char="#");
	cat("Debug target variable list:\n");
	print(debug_tar_list);
	curated_targets_mat=curated_targets_mat[,debug_tar_list,drop=F];
}


curated_target_names=colnames(curated_targets_mat);


##############################################################################

num_variables=length(curated_target_names);

params_matrix=matrix(NA, nrow=num_variables, ncol=3);
colnames(params_matrix)=c("mu", "alpha", "beta");
rownames(params_matrix)=curated_target_names;

gn_density=list();

normal_estim_matrix=matrix(NA, nrow=num_variables, ncol=2);
colnames(normal_estim_matrix)=c("mu", "sd");
rownames(normal_estim_matrix)=curated_target_names;

norm_density=list();


failed_variables=c();
par(mfrow=c(2,2));
counter=1;



if(!FastFit){
	# Complete
	mu_pert=3;
	alpha_pert=5;
	beta_pert=7;

	pert_range=4;
	max_successes=20;
}else{
	# Fast
	mu_pert=3;
	alpha_pert=3;
	beta_pert=3;

	pert_range=3;
	max_successes=3
}

# Compute range and steps for perturbation
mu_adj=2^(seq(-pert_range, pert_range, length.out=mu_pert));
alpha_adj=2^(seq(-pert_range, pert_range, length.out=alpha_pert));
beta_adj=2^(seq(-pert_range, pert_range, length.out=beta_pert));

total_perts=mu_pert*alpha_pert*beta_pert;

# Compute how far each perturbation is from 'no perturbation', and order.
pert_matrix=matrix(0, nrow=total_perts, ncol=3);
colnames(pert_matrix)=c("mu", "alpha", "beta");
pert_distance=numeric(total_perts);
pert_ix=1;
for(i in 1:mu_pert){
	for(j in 1:alpha_pert){
		for(k in 1:beta_pert){

			pert_matrix[pert_ix,]=c(mu_adj[i], alpha_adj[j], beta_adj[k]);

			pert_distance[pert_ix]=
				sqrt(log2(mu_adj[i])^2+log2(alpha_adj[j])^2+log2(beta_adj[k])^2);

			pert_ix=pert_ix+1;
		}
	}
}

pert_dist_sort_ix=order(pert_distance);
pert_matrix=pert_matrix[pert_dist_sort_ix,];

cat("Pertubration Matrix (Sorted):\n");
print(pert_matrix);

##############################################################################

MIN_REAL=1e-323;

for(ctn in curated_target_names){
	data=curated_targets_mat[,ctn];

	cat("********************************************************************\n");
	cat("Working on: ", ctn, "\n");

	mu0=mean(data);
	variance=var(data);
	alpha0=sqrt(variance*2);
	beta0=2;

	normal_estim_matrix[ctn, "mu"]=mu0;
	normal_estim_matrix[ctn, "sd"]=sqrt(variance);

	nLL=function(mu, alpha, beta){
		if(alpha<MIN_REAL){alpha=MIN_REAL;};
		if(beta<MIN_REAL){beta=MIN_REAL};
		return(-sum(dgnorm(data, mu, alpha, beta, log=T)));
	};

	cat("Initial Values: mean=",mu0, " alpha=", alpha0, 
		" (var=", variance, ") beta=", beta0, "\n", sep="");
	muf=NA; alphaf=NA; betaf=NA;

	pert_ix=1;

	pert_res_matrix=matrix(NA, nrow=total_perts, ncol=4);
	colnames(pert_res_matrix)=c("nLL", "mu", "alpha", "beta");
	
	num_success=0;

	while(num_success<max_successes && pert_ix<=total_perts){
		success=F;

		#cat("Perturbations: ", pert_ix, "  Trials: ", tries, "\n");

		muTry=mu0*pert_matrix[pert_ix, "mu"];
		alphaTry=alpha0*pert_matrix[pert_ix, "alpha"];
		betaTry=beta0*pert_matrix[pert_ix, "beta"];

	
		for(method in c("Nelder-Mead", "SANN")){
			method_success=F;
			tryCatch(
				{

					cat(method, ": Adj Init: (",
						muTry, ", ",
						alphaTry, ", ",
						betaTry,"): \n",sep="");

					names(muTry)=NULL;
					names(alphaTry)=NULL;
					names(betaTry)=NULL;

					start_param=list(mu=muTry, alpha=alphaTry, beta=betaTry);
					fit=mle(nLL, start=start_param, method);

					#print(fit);
					fit_res=attributes(fit);
					#print(fit_res);
					muf=fit_res$coef[1];	
					alphaf=max(fit_res$coef[2], MIN_REAL);
					betaf=max(fit_res$coef[3], MIN_REAL);
					num_success=num_success+1;
					method_success=T;
				}, 
				error=function(e){
					cat("Error:\n");
					print(e);
				}
			);
			if(method_success){
				cat(method, ": Success!\n");
				break;
			}else{
				cat(method, ": Failure.\n");
			}
		}

		if(method_success){
			nll_pert=nLL(mu=muf, alpha=alphaf, beta=betaf);
			pert_res_matrix[pert_ix,]=c(nll_pert, muf, alphaf, betaf);
		}

		pert_ix=pert_ix+1;		

		cat("\n");
	}

	# Find perturbation with min nll_pert
	nona_ix=!is.na(pert_res_matrix[,"nLL"]);
	nona_pert_res_matrix=pert_res_matrix[nona_ix,,drop=F];

	
	print(nona_pert_res_matrix);

	min_nLL=min(nona_pert_res_matrix[,"nLL"]);
	min_ix=min(which(nona_pert_res_matrix[,"nLL"]==min_nLL));

	mu_best=nona_pert_res_matrix[min_ix,"mu"];
	alpha_best=nona_pert_res_matrix[min_ix,"alpha"];
	beta_best=nona_pert_res_matrix[min_ix,"beta"];

	cat("Best: (", mu_best, ", ", alpha_best, ", ", beta_best, ")\n", sep="");

	hist(data, breaks=30, main=paste(counter, ".) ", ctn, sep=""), freq=F, xlab="");

	if(num_success>0){
		cat("Overall Success.\n")
		title(main=paste("Mean =", round(mu_best,4)), line=-1);
		title(main=paste("Alpha =", round(alpha_best,4)), line=-2);
		title(main=paste("Beta =", round(beta_best,4)), line=-3);

		range=range(data);
		span=diff(range);
		margin=span*.15;
		x=seq(range[1]-margin, range[2]+margin, length.out=100);
		#print(c(muf, alphaf, betaf));
		y=dgnorm(x, mu_best, alpha_best, beta_best);
		#print(x);
		#print(y);
		points(x,y, col="blue");
	
		params_matrix[ctn,]=c(mu_best, alpha_best, beta_best);
		gn_density[[ctn]]=list();
		gn_density[[ctn]][["x"]]=x;
		gn_density[[ctn]][["y"]]=y;

		y=dnorm(x, normal_estim_matrix[ctn, "mu"], normal_estim_matrix[ctn, "sd"]);
		points(x,y, col="red", cex=.5);
		norm_density[[ctn]]=list();
		norm_density[[ctn]][["x"]]=x;
		norm_density[[ctn]][["y"]]=y;
		

	}else{
		cat("Overall Failure.\n");
		title(main="Unsuccessful Fit.", line=-1);
		failed_variables=c(failed_variables, ctn);	
	}


	cat("\n\n");
	counter=counter+1;
}

##############################################################################

dev.off();

##############################################################################

failed_fn=paste(OutputFnameRoot, ".failed.meta.tsv", sep="");
fh=file(failed_fn, "w");
cat(file=fh, "SampleID\t");
close(fh);

write.table(
	curated_targets_mat[,failed_variables,drop=F], 
	file=failed_fn, append=T,
	quote=F, sep="\t");

##############################################################################

failed_fn=paste(OutputFnameRoot, ".failed.lst", sep="");
fh=file(failed_fn, "w");
cat(file=fh, paste(failed_variables, collapse="\n"), "\n", sep="");
close(fh);

##############################################################################
##############################################################################
# Write out the variable histograms sorted by decreasing beta

pdf(paste(OutputFnameRoot, ".sorted.pdf", sep=""), height=11, width=8.5);
par(mfrow=c(3,2));

beta_arr=params_matrix[,"beta"];
beta_arr_sort_ix=order(beta_arr, decreasing=F);

params_matrix_beta_sort=params_matrix[beta_arr_sort_ix,, drop=F];
beta_arr_sorted=params_matrix_beta_sort[, "beta"];
mean_arr_sorted=params_matrix_beta_sort[, "mu"];
var_sorted=rownames(params_matrix_beta_sort);

counter=1;
for(varname in var_sorted){
	
	data=curated_targets_mat[,varname];
	mu=params_matrix[varname,"mu"];
	alpha=params_matrix[varname,"alpha"];
	beta=params_matrix[varname,"beta"];

	hist(data, main=paste(counter, ".) ", varname, sep=""), xlab="", freq=F, breaks=30);
	title(main=paste("Mean =", round(mu,4)), line=-1);
	title(main=paste("Alpha =", round(alpha,4)), line=-2);
	title(main=paste("Beta =", round(beta,4)), line=-3);

	points(gn_density[[varname]][["x"]], gn_density[[varname]][["y"]], col="blue");
	points(norm_density[[varname]][["x"]], norm_density[[varname]][["y"]], col="red", cex=.5);

	counter=counter+1;
}

dev.off();

##############################################################################
##############################################################################

pdf(paste(OutputFnameRoot, ".summaries.pdf", sep=""), height=11, width=8.5);

par(mfrow=c(4,1));
hist(beta_arr, main="Beta Distribution (All)", xlab="");
hist(beta_arr[beta_arr<25], main="Beta Distribution (0-25)", xlab="");
hist(beta_arr[beta_arr<10], main="Beta Distribution (0-10)", xlab="");
hist(beta_arr[beta_arr<5], main="Beta Distribution (0-5)", xlab="");

##############################################################################
# Calculate correlation across variables from different beta values

max_group_size=50;
num_var=length(beta_arr);

low_beta_vnames=var_sorted[1:max_group_size];
high_beta_vnames=var_sorted[(num_var-max_group_size):num_var];

gt2_beta_ix=max(which(beta_arr_sorted>2));
ub=max(1, gt2_beta_ix-max_group_size/2);
lb=min(num_var, gt2_beta_ix+max_group_size/2);
mid_beta_vnames=var_sorted[lb:ub];

cat("Low Beta:\n");
print(beta_arr_sorted[low_beta_vnames]);
cat("Medium Beta:\n");
print(beta_arr_sorted[mid_beta_vnames]);
cat("High Beta:\n");
print(beta_arr_sorted[high_beta_vnames]);

low_beta_values_mat=curated_targets_mat[,low_beta_vnames,drop=F];
mid_beta_values_mat=curated_targets_mat[,mid_beta_vnames,drop=F];
high_beta_values_mat=curated_targets_mat[,high_beta_vnames,drop=F];

par(mfrow=c(1,1));
paint_matrix(cor(low_beta_values_mat), plot_min=-1, plot_max=1, title="Low Beta Variables", 
	value.cex=.5, deci_pts=2);
paint_matrix(cor(mid_beta_values_mat), plot_min=-1, plot_max=1, title="Mid Beta Variables",
	value.cex=.5, deci_pts=2);
paint_matrix(cor(high_beta_values_mat), plot_min=-1, plot_max=1, title="High Beta Variables",
	value.cex=.5, deci_pts=2);


##############################################################################
# compute sliding window

abs_group_mean_cor=function(x){
	mean(abs(as.dist(cor(x))));	
}

slide_window_cor=function(mat, beta_arr, window_size, steps){
	num_var=ncol(mat);

	window_size=min(window_size, round(num_var/8, 0));
	steps=min(num_var, steps);
	
	indices=as.integer(seq(1, num_var-window_size, length.out=steps));
	separation=indices[2]-indices[1];

	if(any(indices<0)){
		return(rep(NA, num_var));
	}
	
	mean_abs_cor=numeric(steps);

	i=1;
	for(offset in indices){
		mean_abs_cor[i]=abs_group_mean_cor(mat[,offset:(offset+window_size)]);
		#cat(offset, "-", (offset+window_size), "\n");
		i=i+1;
	}	

	res=list();
	res$x=beta_arr[indices];
	res$y=mean_abs_cor;
	res$steps=steps;
	res$window_size=window_size;
	res$indices=indices;
	res$separation=separation;

	return(res);
}

mac=slide_window_cor(
	curated_targets_mat[,var_sorted],
	beta_arr_sorted,
	window_size=100,
	steps=300);

cat("Betas:\n");
print(mac$x);

cat("Mean(Abs(Correl)):\n");
print(mac$y);

##############################################################################
# windowed PCA

pca_group_counts=function(x, indiv_pc_cutoff=0.03, cumulative_pc_cutoff=.80){

	correl_mat=cor(x);
	eigen=eigen(correl_mat);

	pca_propvar=eigen$values/sum(eigen$values);
	pca_propcumsum=cumsum(pca_propvar);

	num_pc_to_exceed_cutoff=sum(pca_propcumsum<cumulative_pc_cutoff)+1;
	num_pc_above_indiv_cutoff=sum(pca_propvar>indiv_pc_cutoff);

	return(c(num_pc_to_exceed_cutoff, num_pc_above_indiv_cutoff));
}

slide_window_pca=function(mat, beta_arr, indiv_pc_cutoff, cumul_pc_cutoff, window_size, steps){
	num_var=ncol(mat);

	window_size=min(window_size, round(num_var/8, 0));
	steps=min(num_var, steps);
	
	indices=as.integer(seq(1, num_var-window_size, length.out=steps));
	separation=indices[2]-indices[1];

	if(any(indices<0)){
		return(rep(NA, num_var));
	}
	
	num_pc_to_exceed_cutoff=numeric(steps);
	num_pc_above_indiv_cutoff=numeric(steps);

	i=1;
	for(offset in indices){
		counts=pca_group_counts(mat[,offset:(offset+window_size)], indiv_pc_cutoff, cumul_pc_cutoff);
		num_pc_to_exceed_cutoff[i]=counts[1];
		num_pc_above_indiv_cutoff[i]=counts[2];
		i=i+1;
	}	

	res=list();
	res$x=beta_arr[indices];
	res$y_num_pc_to_exceed_cutoff=num_pc_to_exceed_cutoff;
	res$y_num_pc_above_indiv_cutoff=num_pc_above_indiv_cutoff;
	res$steps=steps;
	res$window_size=window_size;
	res$indices=indices;
	res$separation=separation;

	return(res);
}

pcam=slide_window_pca(
	curated_targets_mat[,var_sorted],
	beta_arr_sorted,
	Indiv_PC_Cutoff,
	Cumul_PC_Cutoff,
	window_size=100,
	steps=300);

print(pcam$indices);

cat("Betas:\n");
print(pcam$x);

cat("PCA Measures:\n");
print(pcam$y_num_pc_to_exceed_cutoff);
print(pcam$y_num_pc_above_indiv_cutoff);

###############################################################################
# Calculate Global Correl

global_abs_correl=function(x, reference){
	cormat=abs(cor(x, reference));
	ref_cor=apply(cormat, 2, mean);
	return(ref_cor);
}

slide_window_global_cor=function(mat, beta_arr, num_glob_PCs, window_size, steps){

	cat("Calculating sliding window global correlation.\n");
	cat("Num Global PCs to include: ", num_glob_PCs, "\n");
	num_var=ncol(mat);

	window_size=min(window_size, round(num_var/8, 0));
	steps=min(num_var-window_size, steps);
	
	indices=as.integer(seq(1, num_var-window_size, length.out=steps));
	separation=indices[2]-indices[1];

	cat("Num Variables: ", num_var, "\n");
	cat("Window Size: ", window_size, "\n");
	cat("Steps: ", steps, "\n");
	#print(indices);

	subsample_size=min(num_var/5, 200);
	cat("Subsample size for Glocal PCA: ", subsample_size, "\n");
	mat_subsamp=mat[,sample(num_var, subsample_size)];

	sample_cor_mat=cor(mat_subsamp);	
	eigen=eigen(sample_cor_mat);
	#print(eigen);
	scores=(scale(mat_subsamp, center=T, scale=T) %*% eigen$vectors);
	target_scores=scores[,1:num_glob_PCs];
	#print(scores);

	glob_abs_cor_mat=matrix(NA, nrow=steps, ncol=num_glob_PCs);

	i=1;
	for(offset in indices){
		glob_abs_cor_mat[i,]=global_abs_correl(
			mat[,offset:(offset+window_size)], 
			target_scores);
		i=i+1;
	}	

	res=list();
	res$x=beta_arr[indices];
	res$y_mat=glob_abs_cor_mat;
	res$steps=steps;
	res$window_size=window_size;
	res$indices=indices;
	res$perc_var=eigen$values/sum(eigen$values);
	res$separation=separation;

	return(res);
}

magc=slide_window_global_cor(
	curated_targets_mat[,var_sorted],
	beta_arr_sorted,
	Num_PC_Global_Correl,
	window_size=100,
	steps=300);

cat("Betas:\n");
print(magc$x);

cat("Global Mean(Abs(Correl)):\n");
print(magc$y_mat);

###############################################################################
# Generate Plots

beta_landmarks=sort(unique(c(0:5, 1.5, 10, 25, 50, 100, 200, 400)));
par(mfrow=c(4,1));

# Plot Beta vs. Metric
range_x=range(beta_arr);
cat("Beta Range:\n");
print(range_x);

cat("Plotting Beta vs. AutoCor...\n");
plot(mac$x, mac$y, xlab="Beta", ylab="Mean Abs(Cor)", main="Abs Auto-Correlation", 
	xlim=range_x,
	ylim=c(0,1));
title(main=paste(
	"Window Size=", mac$window_size, "  Steps=", mac$steps, "  Sep=", mac$separation, sep=""), 
	line=.25, cex.main=.8);
abline(v=beta_landmarks, col="blue", lty="dotted");

cat("Plotting Beta vs. Num PCS to Exceed Cutoff\n");
plot(pcam$x, pcam$y_num_pc_to_exceed_cutoff, 
	xlim=range_x,
	ylim=c(0, max(pcam$y_num_pc_to_exceed_cutoff)+2),
	xlab="Beta", ylab="Num PCs", main=paste("To Exceed Cumulative Cutoff:", Cumul_PC_Cutoff));
abline(v=beta_landmarks, col="blue", lty="dotted");

cat("Plotting Beta vs. Num PCS to Exceed Cutoff\n");
plot(pcam$x, pcam$y_num_pc_above_indiv_cutoff, 
	xlim=range_x,
	ylim=c(0, max(pcam$y_num_pc_above_indiv_cutoff)+2),
	xlab="Beta", ylab="Num PCs", main=paste("To Exceed Individual Cutoff:", Indiv_PC_Cutoff));
abline(v=beta_landmarks, col="blue", lty="dotted");

hist(beta_arr, main="Beta Distribution", xlab="", xlim=range_x, breaks=100);
abline(v=beta_landmarks, col="blue", lty="dotted");

# Plot Log10(Beta) vs. Metric
range_logx=range(log10(beta_arr));


cat("Plotting Log10(Beta) vs. AutoCor...\n");
plot(log10(mac$x), mac$y, xlab="Log10(Beta)", ylab="Mean Abs(Cor)", main="Abs Auto-Correlation", 
	xlim=c(range_logx),
	ylim=c(0,1));
title(main=paste(
	"Window Size=", mac$window_size, "  Steps=", mac$steps, "  Sep=", mac$separation, sep=""), 
	line=.25, cex.main=.8);
abline(v=log10(beta_landmarks), col="blue", lty="dotted");
text(log10(beta_landmarks), 0, beta_landmarks, font=2);

cat("Plotting Log10(Beta) vs. Num PCS to Exceed Cutoff\n");
plot(log10(pcam$x), pcam$y_num_pc_to_exceed_cutoff,
	xlim=c(range_logx),
	ylim=c(0, max(pcam$y_num_pc_to_exceed_cutoff)+2),
	xlab="Log10(Beta)", ylab="Num PCs", main=paste("To Exceed Cumulative Cutoff:", Cumul_PC_Cutoff));
abline(v=log10(beta_landmarks), col="blue", lty="dotted");
text(log10(beta_landmarks), 0, beta_landmarks, font=2);

cat("Plotting Log10(Beta) vs. Num PCS to Exceed Cutoff\n");
plot(log10(pcam$x), pcam$y_num_pc_above_indiv_cutoff,
	xlim=c(range_logx),
	ylim=c(0, max(pcam$y_num_pc_above_indiv_cutoff)+2),
	xlab="Log10(Beta)", ylab="Num PCs", main=paste("To Exceed Individual Cutoff:", Indiv_PC_Cutoff));
abline(v=log10(beta_landmarks), col="blue", lty="dotted");
text(log10(beta_landmarks), 0, beta_landmarks, font=2);

hist(log10(beta_arr), main="Log10(Beta) Distribution", xlab="", 
	xlim=c(range_logx),
	breaks=seq(range_logx[1], range_logx[2], length.out=100));
abline(v=log10(beta_landmarks), col="blue", lty="dotted");

##############################################################################

par(mar=c(2,4,1,1));
par(mfrow=c(Num_PC_Global_Correl+1,1));

for(i in 1:Num_PC_Global_Correl){
	plot(magc$x, magc$y_mat[,i],
		ylab="Abs(Correlation)",
		xlim=range_x,
		ylim=c(0,1),
		main=paste("PC[",i,"]: ", paste(round(magc$perc_var[i]*100,2)), "%", sep="")
	);
	abline(v=beta_landmarks, col="blue", lty="dotted");
}
hist(beta_arr, main="Beta Distribution", xlab="", xlim=range_x, breaks=100);
abline(v=beta_landmarks, col="blue", lty="dotted");

for(i in 1:Num_PC_Global_Correl){
	#plot(log10(magc$x), magc$y_mat[,i],
	plot(0, type="n",
		ylab="Abs(Correlation)",
		xlim=range_logx,
		ylim=c(0,1),
		main=paste("PC[",i,"]: ", paste(round(magc$perc_var[i]*100,2)), "%", sep="")
	);
	abline(v=log10(beta_landmarks), col="blue", lty="dotted");
	text(log10(beta_landmarks), 0, beta_landmarks, font=2);

	range_y=range(magc$y_mat[,i]);
	abline(h=range_y[1], col="blue", lwd=2);
	abline(h=range_y[2], col="red", lwd=2);

	points(log10(magc$x), magc$y_mat[,i]);
	
}
hist(log10(beta_arr), main="Log10(Beta) Distribution", xlab="", 
	xlim=c(range_logx),
	breaks=seq(range_logx[1], range_logx[2], length.out=100));
abline(v=log10(beta_landmarks), col="blue", lty="dotted");


##############################################################################
mean_zero_cutoff=-.025;

par(mar=c(4,4,4,1));
par(mfrow=c(3,1));

beta_arr_sorted=params_matrix_beta_sort[beta_arr_sort_ix, "beta", drop=F];
mean_arr_sorted=params_matrix_beta_sort[beta_arr_sort_ix, "mu", drop=F];

hist(mean_arr_sorted, main="Mean Distribution", xlab="means", breaks=100);

plot(log10(beta_arr_sorted), mean_arr_sorted,
	xlab="Log10(Beta)",
	xlim=range_logx,
	ylab="Mean",
	main="Mean vs. Beta"
);
abline(v=log10(beta_landmarks), col="blue", lty="dotted");
abline(h=mean_zero_cutoff, col="red", lty="dotted");
text(log10(beta_landmarks), 0, beta_landmarks, font=2);

hist(log10(beta_arr), main="Log10(Beta) Distribution", xlab="", 
	xlim=c(range_logx),
	breaks=seq(range_logx[1], range_logx[2], length.out=100));
abline(v=log10(beta_landmarks), col="blue", lty="dotted");

##############################################################################

# Look at means greater than >-.025

mean_zeros_ix=(mean_arr_sorted>mean_zero_cutoff);
max_mean=max(mean_arr_sorted);

mean_zeros_mean=mean_arr_sorted[mean_zeros_ix];
mean_zeros_logbeta=log10(beta_arr_sorted[mean_zeros_ix]);


hist(mean_zeros_mean, 
	main=paste("Mean Distribution (Cutoff=", mean_zero_cutoff, ")"), 
	xlab="means", breaks=100);

plot(mean_zeros_logbeta, mean_zeros_mean,
	xlab="Log10(Beta)",
	xlim=range_logx,
	ylab="Mean",
	main="Means at Cutoff vs. Beta"
);
abline(v=log10(beta_landmarks), col="blue", lty="dotted");
text(log10(beta_landmarks), 0, beta_landmarks, font=2);

hist(mean_zeros_logbeta, main="Log10(Beta) Distribution of Means at Cutoff", xlab="",
	xlim=c(range_logx),
	breaks=seq(range_logx[1], range_logx[2], length.out=100));
abline(v=log10(beta_landmarks), col="blue", lty="dotted");

##############################################################################

#

par(mfrow=c(2,2));

log10_beta=log10(params_matrix[,"beta"]);
beta_log10_range=range(log10_beta);

#print(params_matrix[params_matrix[,"alpha"]<=0,"alpha"]);

#quit();


log10_alpha=log10(params_matrix[,"alpha"]);
alpha_log10_range=range(log10_alpha);

mu_range=range(params_matrix[,"mu"]);



plot(log10_beta, log10_alpha,
	xlab="log10(beta), shape", ylab="log10(alpha), scale",
	xlim=beta_log10_range, ylim=alpha_log10_range
);

plot(0,type="n", xlab="", ylab="", main="", xaxt="n", yaxt="n", bty="n");

plot(log10_beta, params_matrix[,"mu"],
	xlab="log10(beta), shape", ylab="mean (mu), location",
	xlim=beta_log10_range, ylim=mu_range
);

plot(log10_alpha, params_matrix[,"mu"],
	xlab="log10(alpha), scale", ylab="mean (mu), location",
	xlim=alpha_log10_range, ylim=mu_range
);

mtext("Scatter Plots for Generalized Normal Parameters", side=3, line=0, outer=T, font=2);

##############################################################################

BetaExportCutoff

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
