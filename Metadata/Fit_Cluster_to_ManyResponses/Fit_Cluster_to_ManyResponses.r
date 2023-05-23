#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

options(useFancyQuotes=F);

params=c(
	"factor_fn", "f", 1, "character",
	"responses_fn", "r", 1, "character",
	"covariates_fn", "c" , 1, "character",
	"cluster_names_fn", "k", 1, "character",
	"outputroot", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"\n",
	"	-f <Factor Filename>\n",
	"	-r <Response Variables List, filename>\n",
	"	-c <Covariates List, and other predictors, filename>\n",
	"	-k <cluster cut (categorical variable) name list, filename>\n",
	"	-o <Output Filename Root>\n",
	"\n",
	"This script will fit the following regression model:\n",
	"\n",
	"	[response_n] = [covariates] + [cluster-cut]\n",
	"\n",
	"\n", sep="");

if(!(length(opt$factor_fn) && 
	length(opt$responses_fn) &&
	length(opt$covariates_fn) && 
	length(opt$cluster_names_fn) && 
	length(opt$outputroot)
)){
	cat(usage);
	q(status=-1);
}

FactorsFile=opt$factor_fn;
ResponsesFile=opt$responses_fn;
CovariatesFile=opt$covariates_fn;
ClusterNamesFile=opt$cluster_names_fn;
OutputRoot=opt$outputroot;

params=capture.output({
cat("\n");
cat("      Factors File: ", FactorsFile, "\n", sep="");
cat("    Responses File: ", ResponsesFile, "\n", sep="");
cat("   Covariates File: ", CovariatesFile, "\n", sep="");
cat("Cluster Names File: ", ClusterNamesFile, "\n", sep="");
cat("       Output Root: ", OutputRoot, "\n", sep="");
cat("\n");
});

print(params);
options(width=120);

##############################################################################

load_factors=function(fname){
	cat("Loading Factors/Metadata: ", fname, "\n", sep="");
	factors=data.frame(read.table(fname,  sep="\t", header=TRUE, row.names=1, check.names=FALSE));
	return(factors);
}

load_map=function(fname){
	# Variable name / group
	map=read.delim(fname, header=F, sep="\t", comment.char="#", stringsAsFactors=F);

	if(ncol(map)<2){
		cat("Error:  Map file needs two columns.\n");
		quit(status=-1);
	}

	map_list=list();
	grps=sort(unique(map[,2]));
	for(g in grps){
		map_list[[g]]=map[map[,2]==g, 1];
	}

	return(map_list);
}

load_list=function(filename){
	cat("Loading List: ", filename, "\n", sep="");
	val=scan(filename, what=character(), comment.char="#");
	return(val);
}

plot_text=function(strings, max_lines_pp=Inf, oma_tag=""){

        orig.par=par(no.readonly=T);

        par(mfrow=c(1,1));
        par(family="Courier");

	if(oma_tag==""){
		par(oma=rep(.5,4));
	}else{
		par(oma=c(1, .5, .5, .5));
	}

        par(mar=rep(0,4));

        num_lines=length(strings);
        num_pages=max(1, ceiling(num_lines/max_lines_pp));

        cat("Num Pages for ", num_lines, " lines: ", num_pages, "\n", sep="");

        lines_pp=min(num_lines, max_lines_pp);
        for(p in 1:num_pages){

                top=max(as.integer(lines_pp), 52);
	
                plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
                        xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                        );

		if(oma_tag!=""){
			mtext(paste("[", oma_tag, "]", sep=""), side=1, col="grey25");
		}

                text_size=max(.01, min(.7, .7 - .003*(lines_pp-52)));
                #print(text_size);

                start=(p-1)*lines_pp+1;
                end=start+lines_pp-1;
                end=min(end, num_lines);
                line=1;
                for(i in start:end){
                        #cat(strings[i], "\n", sep="");
                        strings[i]=gsub("\t", "", strings[i]);
                        text(0, top-line, strings[i], pos=4, cex=text_size);
                        line=line+1;
                }

        }

        par(orig.par);
}

plot_title_page=function(title, subtitle=""){

        orig.par=par(no.readonly=T);
        par(family="serif");
        par(mfrow=c(1,1));

        plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n",  xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                );

        # Title
        title_cex=3;
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


plot_page_separator=function(title, subtitle="", notes=""){

	par(mfrow=c(1,1));
	plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n",  xaxt="n", yaxt="n",
		xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
		);

	text(.5, .8, title, adj=c(.5, -1), cex=4, font=2); 
	text(.5, .5, subtitle, adj=c(.5, 0), cex=2, font=1); 
	text(.5, .4, notes, adj=c(.5, 1), cex=1, font=3); 
}

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, deci_pts=4, 
	label_zeros=T, counts=F, value.cex=2, 
	plot_col_dendr=F,
	plot_row_dendr=F,
	v_guide_lines=c(), h_guide_lines=c()
){

        num_row=nrow(mat);
        num_col=ncol(mat);

	row_names=rownames(mat);
	col_names=colnames(mat);

	orig.par=par(no.readonly=T);

	cat("Painting Matrix: ", title, "\n");
        cat("Num Rows: ", num_row, "\n");
        cat("Num Cols: ", num_col, "\n");


	if(num_row==0 || num_col==0){
		plot(0, type="n", xlim=c(-1,1), ylim=c(-1,1), xaxt="n", yaxt="n", bty="n", xlab="", ylab="",
			main=title);
		text(0,0, "No data to plot...");
		return();
	}

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
                plot_min=min(c(mat,Inf), na.rm=T);
        }
        if(is.na(plot_max)){
                plot_max=max(c(mat,-Inf), na.rm=T);
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

	if(num_row==1){
		plot_row_dendr=F;
	}
	if(num_col==1){
		plot_col_dendr=F;
	}

	# Don't plot dendrogram if there are any NAs in the matrix
	#if(any(is.na(mat))){
	#	plot_col_dendr=F;
	#	plot_row_dendr=F;
	#}

	if(plot_col_dendr && plot_row_dendr){
		layoutmat=matrix(
			c(
			rep(c(rep(4, row_dend_width), rep(3, heatmap_width)), col_dend_height),
			rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height)
			), byrow=T, ncol=row_dend_width+heatmap_width);

		col_dendr=get_dendrogram(mat, type="col");
		row_dendr=get_dendrogram(mat, type="row");

		mat=mat[row_dendr[["names"]], col_dendr[["names"]], drop=F];
		
	}else if(plot_col_dendr){
		layoutmat=matrix(
			c(
			rep(rep(2, heatmap_width), col_dend_height),
			rep(rep(1, heatmap_width), heatmap_height)
			), byrow=T, ncol=heatmap_width); 

		col_dendr=get_dendrogram(mat, type="col");
		mat=mat[, col_dendr[["names"]], drop=F];
		
	}else if(plot_row_dendr){
		layoutmat=matrix(
			rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height),
			byrow=T, ncol=row_dend_width+heatmap_width);

		row_dendr=get_dendrogram(mat, type="row");
		mat=mat[row_dendr[["names"]],,drop=F];
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

                        if(mat[y,x]!=0 || label_zeros){
                                if(counts){
                                        text_lab=sprintf("%i", mat[y,x]);
                                }else{
                                        text_lab=sprintf(paste("%0.", deci_pts, "f", sep=""), mat[y,x]);
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

	if(length(v_guide_lines)){
		for(i in 1:length(v_guide_lines)){
			points(rep(v_guide_lines[i],2), c(0, num_row), col="grey50", type="l", lwd=.75);
		}
	}

	if(length(h_guide_lines)){
		for(i in 1:length(h_guide_lines)){
			points(c(0, num_col), rep(h_guide_lines[i],2), col="grey50", type="l", lwd=.75);
		}
	}

	par(orig.par);

}

calc_guide_lines=function(num_cells, min_cuts=5, max_cuts=10){

	if(num_cells<max_cuts){return(c())}

	cuts=min_cuts:max_cuts;
	cut_remainders=num_cells %% cuts;
	names(cut_remainders)=cuts;
	sorted_remainders=sort(cut_remainders, method="shell");
	recommended_cut=as.numeric(names(sorted_remainders[1]));
	guides=seq(0, num_cells, recommended_cut);
	num_guides=length(guides);
	return(guides[2:(num_guides-1)]);
}

sig_char=function(val){
	if(!is.null(val) && !is.nan(val) && !is.na(val)){
		if(val <= .0001){ return("***");}
		if(val <= .001 ){ return("** ");}
		if(val <= .01  ){ return("*  ");}
		if(val <= .05  ){ return(":  ");}
		if(val <= .1   ){ return(".  ");}
	}
	return(" ");
}

plot_text_mini=function(text, main=""){
	plot(0,0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", main=main, bty="n",
		xaxt="n", yaxt="n"
		);
	#print(text);

	nlines=length(text);
	cwidth=par()$cxy[1];
	cheight=par()$cxy[2];

	max_lines_wo_scaling=1/cheight;
	scale=1;
	if(nlines>max_lines_wo_scaling){
		scale=max_lines_wo_scaling/nlines/.65;
	}
	text(0, 1, paste(text, collapse="\n"), family="mono", cex=scale, adj=c(0,1)); 

}

##############################################################################
##############################################################################

factors_loaded=load_factors(FactorsFile);
cat("Factor File:\n");
factors_num_samples=nrow(factors_loaded);
factors_num_variables=ncol(factors_loaded);
cat("Num Samples: ", factors_num_samples, " x Num Variables: ", factors_num_variables, "\n");

available_variables=colnames(factors_loaded);
covariates_arr=load_list(CovariatesFile);
cluster_names_arr=load_list(ClusterNamesFile);
response_var_arr=load_list(ResponsesFile);

num_response_variables=length(response_var_arr);

num_cluster_var=length(cluster_names_arr);

pdf(paste(OutputRoot, ".mresp_clst.pdf", sep=""), height=11, width=8.5);

var_info=capture.output({
	cat("Covariates:\n");
	print(covariates_arr);
	cat("\n");
	cat("Clusters:\n");
	print(cluster_names_arr);
	cat("\n");
	cat("Number of Response Variables: ", num_response_variables, "\n");
	cat("Example Response Variables:\n");
	print(response_var_arr[1:min(15, num_response_variables)]);
	});

print(var_info, quote=F);
plot_text(var_info);

##############################################################################

#response_var_arr=response_var_arr[1:20];

library('digest');
response_matrix=factors_loaded[,response_var_arr,drop=F];
checksums=apply(response_matrix, 2, digest);
num_checksums=length(checksums);
unique_resp_cs=unique(checksums);
num_unique_response_cs=length(unique_resp_cs);

cat("Unique Responses: ", num_unique_response_cs, " of ", num_checksums, "\n");

unique_responses=list();
for(ucs in unique_resp_cs){
	keep_ix = which(checksums==ucs);
	checksum_group=names(keep_ix);
	representative=checksum_group[1];
	unique_responses[[representative]]=checksum_group;
}

unique_responses_var_arr=names(unique_responses);
num_unique_responses=length(unique_responses_var_arr);


##############################################################################

if(length(covariates_arr)){
	covariates_string=paste(" + ", paste(covariates_arr, collapse=" + "));
}else{
	covariates_string="";
}

best_cluster_fit=list();
i=1;
for(resp_nm in unique_responses_var_arr){

	cat("Working on Response: ", resp_nm, " [", i, "/", num_unique_responses, "]\n");
	resp=factors_loaded[,resp_nm];

	cluster_pval=numeric(num_cluster_var);
	names(cluster_pval)=cluster_names_arr;
	cluster_results=list();
	for(clust_nm in cluster_names_arr){

		cat("Working on Cluster: ", clust_nm, "\n");

		# Build model
		model_string=paste("resp ~ ", clust_nm, covariates_string, sep="");
		#print(model_string);

		fit=lm(formula(model_string), data=as.data.frame(factors_loaded));

		# Store fstat pvalue and fit
		fit_summary=summary(fit);
		fs=fit_summary$fstatistic;

		fstat_pval=1-pf( fs[1], fs[2], fs[3]);

		cluster_results[[clust_nm]]=list();
		cluster_results[[clust_nm]][["fit"]]=fit;
		cluster_results[[clust_nm]][["sum_fit"]]=fit_summary;
		cluster_pval[[clust_nm]]=fstat_pval;

	}

	# File best cluster model
	min_pval=min(cluster_pval);
	best_model_ix=min(which(cluster_pval==min_pval));
	most_signf_model=cluster_names_arr[best_model_ix];

	# Store best model	
	fit_details=list();
	fit_details[["clust_id"]]=most_signf_model;	
	fit_details[["fit"]]=cluster_results[[most_signf_model]];
	fit_details[["pval"]]=min_pval;
	fit_details[["all_clust_pvals"]]=unlist(cluster_pval);

	best_cluster_fit[[resp_nm]]=fit_details;

	cat("Best Cluster Model: ", most_signf_model, " / ", min_pval, "\n");
	cat("\n");

	i=i+1;
}

##############################################################################

resp_pval=numeric(num_unique_responses);
names(resp_pval)=unique_responses_var_arr;

for(resp_nm in unique_responses_var_arr){
	resp_pval[resp_nm]=best_cluster_fit[[resp_nm]][["pval"]];
}

resp_pval_order=order(resp_pval);
resp_pval_sorted=resp_pval[resp_pval_order];
fdr_adj=p.adjust(resp_pval_sorted, method="fdr");

cat("Top Unadjusted P-values:\n");
print(head(resp_pval_sorted, 10));
cat("\n");
cat("Top FDR Adjusted P-values:\n");
print(head(fdr_adj, 10));
cat("\n");

min_pval=resp_pval_sorted[1];
nl10_min_pval=-log10(min_pval);

plot(resp_pval_sorted, fdr_adj, main="FDR Adjusted P-values",
	xlab="Original Unadjusted P-values", ylab="FDR Adjusted P-values");

##############################################################################

signf_char=function(x){

	sc=lapply(x, function(x){
		if(x<.001){
			return("***");
		}else if(x<.01){
			return("**");
		}else if(x<.05){
			return("*");
		}else if(x<.1){
			return(".");
		}
		return("ns");
	});

	return(sc);

}

layout_mat=matrix(c(
	1,1,3,3,
	4,4,4,4,
	2,2,2,2,
	2,2,2,2,
	2,2,2,2), byrow=T, ncol=4);
layout(layout_mat);

resp_sorted_names=names(resp_pval_sorted);
par(oma=c(0,0,4,0));

pval_ticks=c(.1, .05, .01, .001, .0001);
nl10_pval_ticks=-log10(pval_ticks);
for(resp_name in resp_sorted_names){

	best_cl_fit=best_cluster_fit[[resp_name]];	

	#----------------------------------------------------------------------
	# Plot pvalues across clusters
	nl10_pval=-log10(best_cl_fit[["all_clust_pvals"]]);
	par(mar=c(6.1, 4.1, 4.1, 4.1));
	plot(1:num_cluster_var, nl10_pval, type="b",
		xaxt="n",
		main="Cluster Config vs. Unadjusted P-values",
		xlab="", ylab="-Log10(P-values)", ylim=c(0, nl10_min_pval),
		las=2
		);
	max_nl10pval=max(nl10_pval);
	abline(h=max_nl10pval, col="blue", lty="dotted");
	
	cl_ix=which(best_cl_fit[["clust_id"]]==cluster_names_arr);
	abline(v=cl_ix, col="blue");
	for(i in 1:num_cluster_var){
		axis(side=1, at=i, labels=cluster_names_arr[i], las=2, font=ifelse(i==cl_ix, 2, 1));
	}

	axis(side=4, at=nl10_pval_ticks, labels=pval_ticks, las=2);
	mtext(line=.25, sprintf("Best p-val = %5.4g", best_cl_fit[["pval"]]), 
		cex=.75, font=3, col="blue");

	#----------------------------------------------------------------------

	# Split coefficients into clusters and covariates 
	coefficients=best_cl_fit[["fit"]][["sum_fit"]]$coefficients;
	num_coefficients=nrow(coefficients);
	coef_rownames=rownames(coefficients);
	
	coef_names=setdiff(rownames(coefficients), "(Intercept)");
		
	cluster_categories=setdiff(sort(unique(as.character(factors_loaded[, best_cl_fit[["clust_id"]]]))), NA);
	num_categories=length(cluster_categories);

	# Allocate matrix to store cluster coefficients
	cluster_category_coefficients=matrix(c(0, 1), byrow=T, nrow=num_categories, ncol=2);
	rownames(cluster_category_coefficients)=cluster_categories;
	colnames(cluster_category_coefficients)=c("Coefficient", "P-value");

	# Allocate matrix to store covariates coefficients
	covariates_coefficients=matrix(c(0, 1), byrow=T, nrow=0, ncol=2);
	colnames(covariates_coefficients)=c("Coefficient", "P-value");
	cov_ix=1;
	covariates_names=c();

	# Copy cluster coefficients to matrix
	clid_prefix=paste("^", best_cl_fit[["clust_id"]], sep="");
	for(i in 1:num_coefficients){
		if(length(grep(clid_prefix, coef_rownames[i]))){
			# Clusters
			category_name=gsub(clid_prefix, "", coef_rownames[i]);
			cluster_category_coefficients[category_name,"Coefficient"]=coefficients[i, "Estimate"];
			cluster_category_coefficients[category_name,"P-value"]=coefficients[i, "Pr(>|t|)"];
		}else{
			# Covariates
			if(coef_rownames[i]!="(Intercept)"){
				covariates_names=c(covariates_names, coef_rownames[i]);
				covariates_coefficients=rbind(
					covariates_coefficients,
					c(coefficients[i, "Estimate"], coefficients[i, "Pr(>|t|)"]));
				cov_ix=cov_ix+1;
			}
		}
	}

	# ---------------------------------------------------------------------
	# Bar plot of cluster coefficients
	par(mar=c(30, 5.1, 4.1, 2.1));
	signif_chars=signf_char(cluster_category_coefficients[,"P-value"]);

	range=c(-1,1) * max(abs(cluster_category_coefficients[,"Coefficient"]));
	mids=barplot(cluster_category_coefficients[,"Coefficient"], las=2,
		main=paste(best_cl_fit[["clust_id"]], ":\nCluster Category Coefficients", sep=""),
		ylim=range,
		col=ifelse(cluster_category_coefficients[,"Coefficient"]>0, "blue", "red")
		);

	text(mids, 0, labels=signif_chars, font=3, 
		pos=ifelse(cluster_category_coefficients[,"Coefficient"]>0, 1, 3),
		cex=ifelse(cluster_category_coefficients[,"P-value"]>.1, .75, 3)
		);

	# ---------------------------------------------------------------------
	# Covariate Coefficient/Pvalues
	par(mar=c(0, .5, 4.1,.5));
	out_table=cbind(
		sprintf("%10.4g", covariates_coefficients[,"Coefficient"]), 
		sprintf("%10.4g", covariates_coefficients[,"P-value"]), 
		covariates_names);


	colnames(out_table)=c("Coefficients", "P-values", "Variable Names");
	options(width=2000);
	out_text=capture.output(print(as.data.frame(out_table), row.names=F));
	plot_text_mini(out_text, main=paste(best_cl_fit[["clust_id"]], ":\nCovariates", sep="") );
	
	# ---------------------------------------------------------------------
	# Cluster Coefficient/Pvalues
	par(mar=c(0, .5, 4.1,.5));
	out_table=cbind(
		sprintf("%10.4g", cluster_category_coefficients[,"Coefficient"]), 
		sprintf("%10.4g", cluster_category_coefficients[,"P-value"]), 
		rownames(cluster_category_coefficients));
	colnames(out_table)=c("Coefficients", "P-values", "Variable Names");
	options(width=2000);
	out_text=capture.output(print(as.data.frame(out_table), row.names=F));
	plot_text_mini(out_text, main=paste(best_cl_fit[["clust_id"]], ":\nCluster Categories", sep="") );


	
	# ---------------------------------------------------------------------
	# Label response in outer margins
	mtext(resp_name, outer=T, line=2, font=2);	
	altnames=unique_responses[[resp_name]];
	if(length(altnames)>1){
		mtext(paste("aliases: ",  paste(altnames, collapse=", "), sep=""), 
			outer=T, line=.5, cex=.75);
	}
}









quit();

###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
