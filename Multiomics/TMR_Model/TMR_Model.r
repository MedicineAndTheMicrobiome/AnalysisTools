#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

options(useFancyQuotes=F);
options(digits=5);
options(width=300);

params=c(
	"factors", "f", 1, "character",
	"groupings", "g", 1, "character",
	"outputroot", "o", 1, "character",
	"covtrt_grp_list_fn", "c", 1, "character",
	"measured_grp_list_fn", "m", 1, "character",
	"response_grp_list_fn", "r", 2, "character",
	"repeated_measure_analysis_id", "a", 2, "character"
);

NO_CHANGE="orig";
NORM_PVAL_CUTOFF=0.20;
MIN_PC_PROP_CUTOFF=0.10;

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

script_path=paste(head(strsplit(script_name, "/")[[1]], -1), collapse="/");

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <factors/metadata file name, (samples x var_names)>\n",
	"	-g <variable groupings, (var_name, grp_name)>\n",
	"	-o <output filename root>\n",
	"\n",
	"	-c <treatment/covariate group list file>\n",
	"	-m <measured variable group list file>\n",
	"	[-r <response group list file>]\n",
	"\n",
	"	[-a <time offset, or repeated measure identifer to insert into output table>]\n",
	"            (This does not affect computations.)\n",
	"\n",
	"This script will automatically generate an analysis for the following\n",
	"groups of variables.\n",
	"\n",
	"Covariates: These are variables that are expected to affect the response\n",
	"	through the measured variables, that are either treatments or covariates.\n",
	"Measured: These are variables which we have measurements, but we aren't\n",
	"	sure of how they interact together, so their associations are latent.\n",
	"Response: These are variables that are the final outcome of the experiment\n",
	"	we are interested in.  They are result of the mechanisms that the\n",
	"	Measured variables are illuminating, but triggered by the Treatments\n",
	"\n",
	"\n", sep="");

if(
	!length(opt$factors) || 
	!length(opt$groupings) || 
	!length(opt$outputroot) || 
	!length(opt$covtrt_grp_list_fn) || 
	!length(opt$measured_grp_list_fn) 
){
	cat(usage);
	q(status=-1);
}

FactorsFname=opt$factors;
VariableGroupingsFname=opt$groupings;
OutputFnameRoot=opt$outputroot;
CovTrtGrpFname=opt$covtrt_grp_list_fn;
MeasuredGrpFname=opt$measured_grp_list_fn;


if(length(opt$response_grp_list_fn)){
	ResponseGrpFname=opt$response_grp_list_fn;
}else{
	ResponseGrpFname="";
}

if(length(opt$repeated_measure_analysis_id)){
	RepeatedMeasureAnalysisID=opt$repeated_measure_analysis_id;
}else{
	RepeatedMeasureAnalysisID=NULL;
}

param_text=capture.output({
	cat("\n");
	cat("Factor/Metadata Filename:            ", FactorsFname, "\n");
	cat("Variable Groupings Filename:         ", VariableGroupingsFname, "\n");
	cat("Output Filename Root:                ", OutputFnameRoot, "\n");
	cat("Covariate/Treatment Groups Filename: ", CovTrtGrpFname, "\n");
	cat("Measured Groups Filename:            ", MeasuredGrpFname, "\n");
	cat("Response Groups Filename:            ", ResponseGrpFname, "\n");
	cat("Repeated Measure Analysis ID:        ", RepeatedMeasureAnalysisID, "\n");
	cat("\n");
});

cat(paste(param_text, collapse="\n"), "\n");

###############################################################################

load_factors=function(fname){
	factors=as.data.frame(read.delim(fname,  header=TRUE, row.names=1, check.names=FALSE, sep="\t"));
	return(factors);
}

#-----------------------------------------------------------------------------#

load_list=function(fname){
	cat("Loading: ", fname, "\n");

	if(fname==""){
		cat("Could not load list: File not specified\n");
		return(c());
	}

	lst=read.delim(fname, header=F, check.names=F, comment.char="#", as.is=T);
	return(lst[,1]);	
}

#-----------------------------------------------------------------------------#

load_groupings=function(fname, var_col=1, grp_col=2){

	cat("Loading Grouping: ", fname, "\n", sep="");
	data=read.table(fname, header=F, as.is=T, comment.char="#");
	#print(data);
	grps=data[,grp_col];
	uniq_grps=sort(unique(grps));

	uniq_vars=sort(unique(data[,var_col]));

	grp_list=list();
	for(grp in uniq_grps){
		tar_var_ix=(grps==grp);
		grp_list[[grp]]=as.character(data[tar_var_ix, var_col]);
	}

	group_map_rec=list();
	group_map_rec[["Groups"]]=uniq_grps;
	group_map_rec[["NumGroups"]]=length(uniq_grps);
	group_map_rec[["GrpVarMap"]]=grp_list;
	group_map_rec[["Variables"]]=uniq_vars;
	group_map_rec[["NumUniqVar"]]=length(uniq_vars);
	
	return(group_map_rec);

}

#-----------------------------------------------------------------------------#

plot_text=function(strings, max_lines_pp=Inf){

	orig.par=par(no.readonly=T);

	par(mfrow=c(1,1));
	par(family="Courier");
	par(oma=rep(.5,4));
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

		text_size=max(.01, min(.8, .8 - .003*(lines_pp-52)));
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

#-----------------------------------------------------------------------------#

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

#-----------------------------------------------------------------------------#

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, deci_pts=4,
        label_zeros=T, counts=F, value.cex=2,
        plot_col_dendr=F,
        plot_row_dendr=F
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
                plot_min=min(mat);
        }
        if(is.na(plot_max)){
                plot_max=max(mat);
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
        if(any(is.na(mat))){
                plot_col_dendr=F;
                plot_row_dendr=F;
        }

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

        par(orig.par);

}

##############################################################################
# Main Program Starts Here!
##############################################################################

pdf(paste(OutputFnameRoot, ".tmr.pdf", sep=""), height=8.5, width=11);
plot_text(param_text);

plot_title_page("TMR", c(
	"The TMR is a framework for integrating various datasets",
	"so that inter-dataset relationships can be systematically identified.",
	"",
	"Datasets can be measurements from different biological compartments,",
	"-omics types, panels, etc.",
	"",
	"In TMR, a dataset is represented by a group of variables that has been previously",
	"curated down (eg. with PCA) so that a majority of its variance can be captured",
	"in fewer variables and also so that the members are not highly co-linear.",
	"It is important not to perform a single PCA across all groups, since an important",
	"goal of the TMR framework is to identify relationships between datasets",
	"through building and testing specific but numerous models, in order to identify an",
	"overall 'network' of interactions among the biological compartments measured.",
	"",
	"Variable groups are assigned into one of 3 groups:",
	"1.) Treatment/Covariates: Variables that are determined by the investigator (Treatments)",
	"or are variables that could affect the experiment, but are essentially static, e.g.",
	"age, sex, smoking status.",
	"2.) Measured: These variables are quantitatively measured and may be",
	"highly dimensional.  These may be affected by the treatment/covariates, or they",
	"may be responsible for affecting another Measured variables.  These were measured",
	"for in the experiment, because they were believe to have an affect on the Response.",
	"3.) Response: These variables are typically measured or used to define the severity",
	"or manifestation of a disease, or the ultimate outcome of the experiment.  The essential",
	"characteristic of the Response is that it is assume to not affect a Measured variable."
));

# Load factors
cat("Loading Factors...\n");
loaded_factors=load_factors(FactorsFname);
loaded_factor_names=colnames(loaded_factors);
loaded_sample_names=rownames(loaded_factors);

options(width=80);
cat("Loaded factors:\n");
print(loaded_factor_names);
cat("\n");
cat("Loaded sample ids:\n");
print(loaded_sample_names);
cat("\n");

# Load Groupings
groupings_rec=load_groupings(VariableGroupingsFname);
print(groupings_rec);

# Load variable types
covariates_list=load_list(CovTrtGrpFname);
measured_list=load_list(MeasuredGrpFname);
response_list=load_list(ResponseGrpFname);

varlists_text=capture.output({
	cat("Number of Samples:", length(loaded_sample_names), "\n");
	cat("\n");
	cat("Variable Group Type Assigments\n");
        cat("\n");
	cat("Covariates:\n");
	print(covariates_list);
	cat("\n");
	cat("Measured:\n");
	print(measured_list);
	cat("\n");
	cat("Response:\n");
	print(response_list);
        cat("\n");
});

cat(paste(varlists_text, collapse="\n"),"\n");
plot_text(varlists_text);

##############################################################################

split_factors_to_groups=function(data_matrix, group_list, grp_arr){

	matrix_vars=colnames(data_matrix);

	data_grp=list();
	for(grp_id in grp_arr){
		cat("Placing '", grp_id, "' into data structures.\n", sep="");

		grp_vars=group_list[["GrpVarMap"]][[grp_id]];

		missing_var=setdiff(grp_vars, matrix_vars);
		if(length(missing_var)==0){
			data_grp[[grp_id]]=data_matrix[,grp_vars, drop=F];
		}else{
			cat("Error:  Could not find variables in data matrix.\n");
			print(missing_var);
			quit(status=-1);
		}
	}
	return(data_grp);

}

#print(loaded_factors);

variables_rec=list();
variables_rec[["Covariates"]]=split_factors_to_groups(loaded_factors, groupings_rec, covariates_list);
variables_rec[["Measured"]]=split_factors_to_groups(loaded_factors, groupings_rec, measured_list);
variables_rec[["Response"]]=split_factors_to_groups(loaded_factors, groupings_rec, response_list);

#print(variables_rec);

# Generate list of how variables will be used
model_variables_summary=capture.output({
	for(type in names(variables_rec)){
		cat(type, ":\n\n");
		for(grp in names(variables_rec[[type]])){
			cat("   ", grp, "\n");
			vars=colnames(variables_rec[[type]][[grp]]);
			for(var in vars){
				cat("      ", var, "\n");
			}
			cat("\n");
		}
		cat("\n\n");
	}
});


plot_title_page("TMR Framework Model Variables", c(
	"The following page lists the variable types and variable groups that",
	"will be used in the TMR Framework."
));

#print(model_variables_summary);
plot_text(c(
	"Model Variables:",
	"",	
	model_variables_summary
), max_lines_pp=80);

cat(paste(model_variables_summary, collapse="\n"));

##############################################################################

standardize_matrix=function(x){
	# Standardize (mean=0, std=1) all columns

	sd=apply(x, 2, sd);
	mean=apply(x, 2, mean);
	
	nrow=nrow(x);
	ncol=ncol(x);

	st_mat=matrix(NA, nrow=nrow, ncol=ncol);
	colnames(st_mat)=colnames(x);
	rownames(st_mat)=rownames(x);

	for(var_ix in 1:ncol){
		st_mat[,var_ix]=(x[,var_ix]-mean[var_ix])/sd[var_ix];
	}

	return(st_mat);

}

calc_mds_matrix=function(x){
	dist=dist(x);
	mds=cmdscale(dist);
	return(mds);
}

apply_fun_to_var_rec=function(var_rec, mat_funct){
	std_var_rec=list();
	for(type in names(var_rec)){
		tmp_list=list();
		for(grp in names(var_rec[[type]])){
			data=var_rec[[type]][[grp]];
			tmp_list[[grp]]=(mat_funct)(data);
		}
		std_var_rec[[type]]=tmp_list;
	}
	return(std_var_rec);
}

cat("Standardizing Matrix...\n");
standardized_variables_rec=apply_fun_to_var_rec(variables_rec, standardize_matrix);
#print(standardized_variables_rec);

cat("Calculate Distances...\n");
dist_rec=apply_fun_to_var_rec(standardized_variables_rec, dist);
#print(dist_rec);

cat("Calculating MDS...\n");
mds_rec=apply_fun_to_var_rec(standardized_variables_rec, calc_mds_matrix);
#print(mds_rec);

##############################################################################

calculate_intergroup_correlation=function(d_rec){

	flatten=list();

	types=names(d_rec);
	for(t in types){
		group_names=names(d_rec[[t]]);
		for(g in group_names){
			name=paste(t, ":", g, sep="");
			flatten[[name]]=d_rec[[t]][[g]];
		}

	}


	#print(flatten);
	num_distmats=length(flatten);
	distmat_names=names(flatten);
	distmat_correl=matrix(NA, nrow=num_distmats, ncol=num_distmats);
	rownames(distmat_correl)=distmat_names;
	colnames(distmat_correl)=distmat_names;

	diag(distmat_correl)=1;

	for(i in 1:num_distmats){
		for(j in 1:num_distmats){
			if(i<j){
				distmat_correl[i,j]=
					cor(flatten[[distmat_names[i]]], flatten[[distmat_names[j]]],
						method="spearman"
						);
			}else{
				distmat_correl[i,j]=distmat_correl[j,i];
			}
		}
	}

	return(distmat_correl);

}

plot_distances_on_line=function(dist_rec, plots_per_page=4){

	dist_mat=as.matrix(dist_rec);
	#print(dist_mat);
	nvar=ncol(dist_mat);
	varnames=rownames(dist_mat);

	orig_par=par(no.readonly=T);
	par(mfrow=c(plots_per_page, 1));
	par(mar=c(3,2,4,2));
	
	for(i in 1:nvar){
		
		dist_arr=dist_mat[i,];
		names(dist_arr)=varnames;

		cat(varnames[i], "\n");
		#print(dist_arr)

		plot(0, type="n", main=paste("1-|cor(dist)| from:  ", varnames[i], sep=""),
			xlim=c(0,1), ylim=c(-.05,3), bty="n", yaxt="n", xlab="", ylab="");

		axis(side=1, at=0, "More Similar", tick=F, line=1, cex.axis=1, font.axis=2);
		axis(side=1, at=1, "More Different", tick=F, line=1, cex.axis=1, font.axis=2);

		for(j in 1:nvar){
			if(j!=i){
				points(dist_arr[j], -0.1, pch=15);
				text(dist_arr[j], 0.1, varnames[j], pos=4, 
					offset=0, srt=90, cex=.75);
			}
		}

	}

	par(orig_par);
}

plot_title_page("Inter-group Similarity Dendrogram", c(
	"The following dendrogram illustrates the similarity of variable groups",
	"based on the measurements collected for each sample.",
	"",
	"Intersample distance matrices were calculated for each of the variable",
	"groups using dist=1-abs(cor) as input towards generating the Ward hierarchical",
	"clusters.",
	"",
	"When two groups cluster closely together, their similarity may",
	"suggest signally between compartments, similarity (redundancy) in",
	"the mechanisms/biomarkers being assayed, etc."  
));

cat("Calculating intergroup correlations...\n");
dist_cor=calculate_intergroup_correlation(dist_rec);
cor_as_dist=1-abs(dist_cor);

hcl=hclust(as.dist(cor_as_dist), method="ward.D2");
plot(hcl, xlab="", ylab="Distance", main="Group Dendrogram: 1-|cor(distances)|");

plot_title_page("Intergroup Distance Heat Map", c(
	"The following heat map illustrates the correlation between the inter-sample distance matrices",
	"computed between the variable groups." 
));

paint_matrix(dist_cor, title="Correlation Among Distances",
	plot_min=-1, plot_max=1, deci_pts=2, value.cex=1);

plot_title_page("Intergroup Distances", c(
	"The following distance plots illustrate the distance between groups",
	"on a line, using each variable group as a reference (0). Groups that",
	"are more similar to the reference variable group will be located towards",
	"the left of the line, whereas variable groups more different from the",
	"reference will be located towards the right."
));

plot_distances_on_line(cor_as_dist);

##############################################################################

covariates_data=matrix(NA, nrow=nrow(loaded_factors), ncol=0);
for(grp in names(variables_rec[["Covariates"]])){
	covariates_data=cbind(covariates_data, variables_rec[["Covariates"]][[grp]]);
}

covariate_variable_names=colnames(covariates_data);

##############################################################################

colors=matrix(NA, nrow=nrow(covariates_data), ncol=ncol(covariates_data));
colnames(colors)=colnames(covariates_data);
rownames(colors)=rownames(covariates_data);

quantize=function(x, steps){
	# Assign value to bin by quantizig
	min=min(x);
	max=max(x);
	range=max-min;
	norm=(x-min)/range;
	return(floor(norm*(steps-1))+1);
}

centers=function(x, steps){
	# Calculate the bin centers in 
	min=min(x);
	max=max(x);
	breaks=seq(min, max, length.out=steps+1);
	breaks=round(breaks+(breaks[2]-breaks[1])/2, 3);
	return(head(breaks, steps));
}

max_cat=5;
legends_values=list();
for(var_ix in covariate_variable_names){
	data=covariates_data[,var_ix];
	colors[,var_ix]=quantize(data, max_cat);
	legends_values[[var_ix]]=centers(data, max_cat);
}

print(colors);
print(legends_values);

palette(rainbow(max_cat, end=4/6));

##############################################################################
# Generate MDS Plots for each group of measurements and response

plot_title_page("MDS Plots", c(
	"The following Multi-Dimensional Scaling plots are based on the",
	"intersample distances calculated on by variable group.",
	"",
	"Samples that are more spatially separated are more different from",
	"each other based on the information captured by the data in each",
	"specific group.",
	"",
	"Although only Measured and Response group variable MDS plots have",
	"been generated, the samples have been colored by the values in the",
	"Treatment/Covariate group."
));

for(type in c("Measured", "Response")){
	for(grp in names(mds_rec[[type]])){
		mds_coord=mds_rec[[type]][[grp]];		
	
		for(var_ix in covariate_variable_names){

			plot(mds_coord[,1], mds_coord[,2], col=colors[, var_ix],
				xlab="Dim 1", ylab="Dim 2", type="n",
				main=paste("Group: ", grp, sep=""),
				);
			title(main=paste("(Colored by: ", var_ix, ")", sep=""), line=.66, cex.main=.75);
			points(mds_coord[,1], mds_coord[,2], col=colors[,var_ix], cex=1.5, lwd=3);
			points(mds_coord[,1], mds_coord[,2], col="black", cex=1.5, lwd=.25);
			text(mds_coord[,1], mds_coord[,2], rownames(mds_coord), cex=.5, pos=3);

			# Determine which part of plot has the least points to place legend
			plot_ranges=par()$usr; # left, right, bottom, top
			xmid=(plot_ranges[1]+plot_ranges[2])/2;
			ymid=(plot_ranges[3]+plot_ranges[4])/2;
			left=sum(mds_coord[,1]<xmid)<sum(mds_coord[,1]>xmid);
			bottom=sum(mds_coord[,2]<ymid)<sum(mds_coord[,2]>ymid);

			xrange=plot_ranges[2]-plot_ranges[1];
			yrange=plot_ranges[4]-plot_ranges[3];

			legend(
				ifelse(left, plot_ranges[1]+xrange/8, plot_ranges[2]-xrange/4),
				ifelse(bottom, plot_ranges[3]+yrange/4, plot_ranges[4]-yrange/8),
				0,0, 
				fill=1:max_cat, border="black",
				title=var_ix,
				legend=legends_values[[var_ix]]);
		}

	}
}

##############################################################################
# Fit Covariates to Predict Measured

model_results=list();

model_results[["Cov_to_Msd"]]=list();
model_results[["Msd_to_Msd"]]=list();
model_results[["Msd_to_Rsp"]]=list();

num_covariates=ncol(covariates_data);

##############################################################################

# Get mapping from variable name to covariates grouping
covtrt_to_group_map=list();
for(gix in covariates_list){
	for(vix in groupings_rec[["GrpVarMap"]][[gix]]){
		covtrt_to_group_map[[vix]]=gix;
	}
}

for(msd_ix in measured_list){

	cat("Fitting: Covariates as Predictor to :", msd_ix, "\n");
	msd_resp=as.matrix(variables_rec[["Measured"]][[msd_ix]]);
	num_resp=ncol(msd_resp);
	msd_varnames=colnames(msd_resp);

	model_string=paste("msd_resp ~ ", paste(covariate_variable_names, collapse=" + "));

	cat("Model: \n");
	print(model_string);	

	fit=lm(as.formula(model_string), data=covariates_data);
	sum_fit=summary(fit);

	# Matrices for storing coef and pvalues
	pval_mat=matrix(NA, nrow=num_covariates, ncol=num_resp);
	rownames(pval_mat)=covariate_variable_names;
	colnames(pval_mat)=msd_varnames;

	coef_mat=matrix(NA, nrow=num_covariates, ncol=num_resp);
	rownames(coef_mat)=covariate_variable_names;
	colnames(coef_mat)=msd_varnames;
	
	# Copy values from summary to matrices
	for(var_ix in names(sum_fit)){
		varname=gsub("Response ", "", var_ix);

		# Sometimes NAs in fit drop the predictor name from coefficients table
		avail_pred=intersect(
			rownames(sum_fit[[var_ix]][["coefficients"]]),
			covariate_variable_names
		);

		pval_mat[avail_pred, varname]=
			sum_fit[[var_ix]][["coefficients"]][avail_pred,"Pr(>|t|)"];

		coef_mat[avail_pred, varname]=
			sum_fit[[var_ix]][["coefficients"]][avail_pred,"Estimate"];

	}

	# Store in record
	model_results[["Cov_to_Msd"]][[msd_ix]][["pval"]]=pval_mat;
	model_results[["Cov_to_Msd"]][[msd_ix]][["coef"]]=coef_mat;

}

#print(model_results[["Cov_to_Msd"]]);

cat("\n\n");

##############################################################################
# Fit (Measured + Covariates) to Predict Measured

for(pred_msd_ix in measured_list){
	for(resp_msd_ix in measured_list){

		if(pred_msd_ix==resp_msd_ix){
			next;
		}

		cat("\n\nFitting: Measured to Measured: (Pred)", pred_msd_ix, " (Resp)", resp_msd_ix, "\n");
		analysis_string=paste(pred_msd_ix, "->", resp_msd_ix, sep="");

		msd_resp=as.matrix(variables_rec[["Measured"]][[resp_msd_ix]]);
		num_resp=ncol(msd_resp);
		msd_resp_varnames=colnames(msd_resp);

		msd_pred=as.matrix(variables_rec[["Measured"]][[pred_msd_ix]]);
		num_pred=ncol(msd_pred);
		msd_pred_varnames=colnames(msd_pred);

		model_string=paste("msd_resp ~ ", 
			paste(c(covariate_variable_names, msd_pred_varnames), collapse=" + "));

		cat("Model: \n");
		print(model_string);	

		fit=lm(as.formula(model_string), data=cbind(covariates_data, msd_pred));
		sum_fit=summary(fit);

		cov_and_pred_names=c(covariate_variable_names, msd_pred_varnames);
		num_cov_pred_var=num_covariates+num_pred;

		# Matrices for storing coef and pvalues
		pval_mat=matrix(NA, nrow=num_cov_pred_var, ncol=num_resp);
		rownames(pval_mat)=cov_and_pred_names;
		colnames(pval_mat)=msd_resp_varnames;

		coef_mat=matrix(NA, nrow=num_cov_pred_var, ncol=num_resp);
		rownames(coef_mat)=cov_and_pred_names;
		colnames(coef_mat)=msd_resp_varnames;
		
		# Copy values from summary to matrices
		for(var_ix in names(sum_fit)){

			varname=gsub("Response ", "", var_ix);

			# Sometimes NAs in fit drop the predictor name from coefficients table
			avail_pred=intersect(
				rownames(sum_fit[[var_ix]][["coefficients"]]),
				cov_and_pred_names
			);

			pval_mat[avail_pred, varname]=
				sum_fit[[var_ix]][["coefficients"]][avail_pred,"Pr(>|t|)"];


			coef_mat[avail_pred, varname]=
				sum_fit[[var_ix]][["coefficients"]][avail_pred,"Estimate"];

		}

		# Store in record
		model_results[["Msd_to_Msd"]][[analysis_string]][["pval"]]=pval_mat;
		model_results[["Msd_to_Msd"]][[analysis_string]][["coef"]]=coef_mat;

	}

}

#print(model_results[["Msd_to_Msd"]]);

cat("\n\n");

##############################################################################
# Fit (Covariates + Measured) to Predict Response

for(pred_msd_ix in measured_list){
	for(resp_ix in response_list){

		cat("Fitting: Measured to (as Predictor)", pred_msd_ix, " to ", resp_msd_ix, " (as Response)\n");

		analysis_string=paste(pred_msd_ix, "->", resp_ix, sep="");

		resp=as.matrix(variables_rec[["Response"]][[resp_ix]]);
		num_resp=ncol(resp);
		resp_varnames=colnames(resp);

		msd_pred=as.matrix(variables_rec[["Measured"]][[pred_msd_ix]]);
		num_pred=ncol(msd_pred);
		msd_pred_varnames=colnames(msd_pred);

		model_string=paste("resp ~ ", 
			paste(c(covariate_variable_names, msd_pred_varnames), collapse=" + "));

		cat("Model: \n");
		print(model_string);	

		fit=lm(as.formula(model_string), data=cbind(covariates_data, msd_pred));
		sum_fit=summary(fit);
		#print(sum_fit);

		cov_and_pred_names=c(covariate_variable_names, msd_pred_varnames);
		num_cov_pred_var=num_covariates+num_pred;

		# Matrices for storing coef and pvalues
		pval_mat=matrix(NA, nrow=num_cov_pred_var, ncol=num_resp);
		rownames(pval_mat)=cov_and_pred_names;
		colnames(pval_mat)=resp_varnames;

		coef_mat=matrix(NA, nrow=num_cov_pred_var, ncol=num_resp);
		rownames(coef_mat)=cov_and_pred_names;
		colnames(coef_mat)=resp_varnames;
		
		# Copy values from summary to matrices
		for(var_ix in names(sum_fit)){
			varname=gsub("Response ", "", var_ix);

			# Sometimes NAs in fit drop the predictor name from coefficients table
			avail_pred=intersect(
				rownames(sum_fit[[var_ix]][["coefficients"]]),
				cov_and_pred_names
			);


			pval_mat[avail_pred, varname]=
				sum_fit[[var_ix]][["coefficients"]][avail_pred,"Pr(>|t|)"];

			coef_mat[avail_pred, varname]=
				sum_fit[[var_ix]][["coefficients"]][avail_pred,"Estimate"];

		}

		# Store in record
		model_results[["Msd_to_Rsp"]][[analysis_string]][["pval"]]=pval_mat;
		model_results[["Msd_to_Rsp"]][[analysis_string]][["coef"]]=coef_mat;

	}
}

#print(model_results[["Msd_to_Rsp"]]);

##############################################################################

matrix_to_tables=function(results, pval_cutoff){
	
	tables=list();
	model_types=names(results);

	model_type=character();
	model_name=character();	
	predictor=character();
	response=character();
	pval=numeric();
	coef=numeric();

	for(mt in model_types){
		cat("Traversing: ", mt, "\n");
		
		tables[[mt]]=list();

		model_names=names(results[[mt]]);
		for(mn in model_names){
			cat("\t",mn, "\n");

			pval_mat=results[[mt]][[mn]][["pval"]];
			coef_mat=results[[mt]][[mn]][["coef"]];

			num_pred=nrow(pval_mat);
			num_resp=ncol(pval_mat);

			pred_names=rownames(pval_mat);
			resp_names=colnames(pval_mat);

			for(pix in 1:num_pred){
				for(rix in 1:num_resp){
					if(!is.na(pval_mat[pix, rix]) &&
					   pval_mat[pix, rix]<=pval_cutoff){
						model_type=c(model_type, mt);
						model_name=c(model_name, mn);
						predictor=c(predictor, pred_names[pix]);
						response=c(response, resp_names[rix]);
						pval=c(pval, pval_mat[pix,rix]);
						coef=c(coef, coef_mat[pix,rix]);
					}
				}
			}

		}

	}

	table=cbind(
		as.data.frame(cbind(model_type, model_name, predictor, response)),
		as.data.frame(cbind(coef, pval))
		);

	return(table);
}

#print(model_results);

denorm_results=list();
for(pvco in rev(c(1.0000, 0.1000, 0.050, 0.010, 0.005, 0.001, 0.0005, 0.0001))){
	denorm_results[[sprintf("%3.4f", pvco)]]=matrix_to_tables(model_results, pvco);
}

#print(denorm_results);

##############################################################################

calc_vertical_spacing=function(num_rows, max_rows_before_squeeze=10, start, end){
	positions=seq(start, end, length.out=max(num_rows, max_rows_before_squeeze));
	if(num_rows<max_rows_before_squeeze){
		offset=ceiling((max_rows_before_squeeze-num_rows)/2);
		return(positions[offset:(offset+num_rows-1)]);
	}else{
		return(positions);
	}
	
}

draw_squares_centered=function(xpos, ypos, height, width, grp_name, variables, text_align){

	#print(c(xpos, ypos, height, width, grp_name, variables, text_align));
	
	points(c(
		xpos-width/2, # tl 
		xpos+width/2, # tr
		xpos+width/2, # br
		xpos-width/2, # bl
		xpos-width/2  # tl
		),
		c(
		ypos+height/2,
		ypos+height/2,
		ypos-height/2,
		ypos-height/2,
		ypos+height/2
		),
		type="l");

	title_cex=.7;
	var_cex=.4;

	text(xpos, ypos+height/2, grp_name, font=2, cex=title_cex, pos=1);
	#text(xpos, ypos+height/2, paste(c(rep("",3), variables), collapse="\n"), cex=.4, pos=1);
	
	title_spc=par()$cxy[2]*title_cex*1.7;
	var_spc=par()$cxy[2]*var_cex;
	
	num_var=length(variables);

	if(num_var>0){
		text_pos=calc_vertical_spacing(num_var, start=ypos+height/2-title_spc, end=ypos-height/2+var_spc);


		#xpad=par()$cxy[1]*var_cex;
		xpad=1/32;
		if(text_align=="left"){
			posv=4;
			text_xpos=xpos-width/2-xpad;
			#aln_adj=c(-0.2,.3);
		}else if(text_align=="right"){
			posv=2;
			text_xpos=xpos+width/2+xpad;
			#aln_adj=c(1.2, .3);
		}else if(text_align=="center"){
			posv=NULL;
			text_xpos=xpos;
		}



		for(i in 1:num_var){
			text(text_xpos, text_pos[i], variables[i], cex=var_cex, pos=posv);
			#text(text_xpos, text_pos[i], variables[i], cex=var_cex, adj=aln_adj);
		}
		return(text_pos);
	}else{
		return(NA);
	}

}


plot_TMR_diagram=function(
	result_rec, title, subtitle="", cvtrt_to_grp_map, 
	grp_links=0,
	cvtrt_grps, msd_grps, rsp_grps){

	cat("Plotting TMR diagram: ", title, "\n");
	cat("Subtitle: ", subtitle, "\n");

	num_cvtrt_grps=length(cvtrt_grps);
	num_msd_grps=length(msd_grps);
	num_rsp_grps=length(rsp_grps);

	#num_trtcov_var=nrow(result_rec[["Cov_to_Msd"]][[1]][["pval"]]);

	#cat("Num Treatment/Covariates variables: ", num_trtcov_var, "\n");

	par(mar=c(0,0,0,0));

	covtrt_xpos=0;
	msd_1_resp_xpos=1;
	msd_1_pred_xpos=1.5;
	msd_2_resp_xpos=2.5;
	msd_2_pred_xpos=3;
	rsp_xpos=4;

	fig_xmar=.125;
	fig_ymar=.0;
	
	plot(0,0,type="n", ylim=c(0,1), xlim=c(covtrt_xpos-fig_xmar, rsp_xpos+fig_xmar),
		bty="n", xaxt="n", yaxt="n");

	#abline(h=c(0,1));

	#abline(v=c(covtrt_xpos, msd_1_xpos, msd_2_xpos, rsp_xpos), lwd=2, col="grey");
	text((rsp_xpos-covtrt_xpos)/2, 1+.0075, title, font=2, cex=2);
	text((rsp_xpos-covtrt_xpos)/2, 1-.025, subtitle, font=2, cex=1.25);

	text(covtrt_xpos, 0, "Covariates &\nTreatments", font=2, cex=1.5);
	text((msd_1_resp_xpos+msd_1_pred_xpos)/2, 0, "Measured", font=2, cex=1.5);
	text((msd_2_resp_xpos+msd_2_pred_xpos)/2, 0, "Measured", font=2, cex=1.5);
	text(rsp_xpos, 0, "Response", font=2, cex=1.5);

	#draw_square_centered(x, y, h, w);
	covtrt_ypos=head(tail(seq(0,1, length.out=(2+num_cvtrt_grps)), -1), -1);
	msd_ypos=head(tail(seq(0,1, length.out=(2+num_msd_grps)), -1), -1);
	rsp_ypos=head(tail(seq(0,1, length.out=(2+num_rsp_grps)), -1), -1);
	

	height_multiplier=.95;
	covtrt_height=(covtrt_ypos[2]-covtrt_ypos[1])*height_multiplier;
	msd_height=(msd_ypos[2]-msd_ypos[1])*height_multiplier;
	rsp_height=(rsp_ypos[2]-rsp_ypos[1])*height_multiplier;

	covtrt_height=ifelse(is.na(covtrt_height), .5, covtrt_height);
	msd_height=ifelse(is.na(msd_height), .5, msd_height);
	rsp_height=ifelse(is.na(rsp_height), .5, rsp_height);

	all_widths=1*.5;
	edge=all_widths/2;

	#----------------------------------------------------------------------

	get_group_linkages=function(rr, ct_grp_map, covtrt_g, msd_g, rsp_g){

		num_result_rows=nrow(rr);

		init_list=function(names){
			outlist=matrix(character(), nrow=0, ncol=5);
			colnames(outlist)=c("pred", "resp", "pred_var", "resp_var", "direction");
			return(outlist);
		}

		covtrt=init_list(covtrt_g);
		msd=init_list(msd_g);
		resp=init_list(rsp_g);

		for(i in 1:num_result_rows){
			if(num_result_rows==0) { break;}

			type=as.character(rr[i, "model_type"]);
			name=as.character(rr[i, "model_name"]);
			pred_var=as.character(rr[i, "predictor"]);
			resp_var=as.character(rr[i, "response"]);
			direction=ifelse(rr[i, "coef"]>0, "+", "-");
	
			grplink=strsplit(name, "->")[[1]];
			pred_grp=grplink[1];
			resp_grp=grplink[2];

			cat("Model Type: ", type, "\n");
		
			if(type=="Cov_to_Msd"){
				pred_grp=ct_grp_map[[pred_var]];
				resp_grp=grplink[1];

				pred_ix=which(pred_grp==covtrt_g);
				resp_ix=which(resp_grp==msd_g);

				covtrt=rbind(covtrt, c(pred_ix, resp_ix, pred_var, resp_var, direction));

			}else if(type=="Msd_to_Msd"){

				pred_ix=which(pred_grp==msd_g);
				resp_ix=which(resp_grp==msd_g);

				msd=rbind(msd, c(pred_ix, resp_ix, pred_var, resp_var, direction));

			}else if(type=="Msd_to_Rsp"){

				pred_ix=which(pred_grp==msd_g);
				resp_ix=which(resp_grp==rsp_g);
				
				resp=rbind(resp, c(pred_ix, resp_ix, pred_var, resp_var, direction));

			}else{
				cat("Type error.\n");
				quit(-1);
			}
		
		}

		add_group_offset=function(em){

			cat("Adding Group Offsets:\n");
			print(em);
			num_entries=nrow(em);

			pred_members=c();
			resp_members=c();
			if(num_entries>0){

				uniq_pred_grp_ix=sort(unique(em[,"pred"]));
				uniq_resp_grp_ix=sort(unique(em[,"resp"]));

				cat("Predictors:\n");
				print(uniq_pred_grp_ix);
				cat("Responders:\n");
				print(uniq_resp_grp_ix);

				pred_members=list();
				for(pred_grp_ix in uniq_pred_grp_ix){
					ingrp=(pred_grp_ix==em[,"pred"]);
					pred_members[[pred_grp_ix]]=sort(unique(em[ingrp,"pred_var"]));
				}
				resp_members=list();
				for(resp_grp_ix in uniq_resp_grp_ix){
					ingrp=(resp_grp_ix==em[,"resp"]);
					resp_members[[resp_grp_ix]]=sort(unique(em[ingrp,"resp_var"]));
				}

				cat("Pred Members:\n");
				print(pred_members);
				cat("Resp Members:\n");
				print(resp_members);			

				offsets=matrix(NA, nrow=0, ncol=2);
				colnames(offsets)=c("pred_off", "resp_off");
				for(i in 1:num_entries){
					pred_grp_ix=em[i, "pred"];
					resp_grp_ix=em[i, "resp"];

					pred_var=em[i, "pred_var"];
					resp_var=em[i, "resp_var"];

					pred_off=which(pred_members[[pred_grp_ix]]==pred_var);
					resp_off=which(resp_members[[resp_grp_ix]]==resp_var);

					offsets=rbind(offsets, c(pred_off, resp_off))

				}

				em=cbind(em, offsets);
			}

			res=list();
			res[["grp_links"]]=em;
			res[["pred_grp_members"]]=pred_members;
			res[["resp_grp_members"]]=resp_members;

			return(res);
		}

		grp_link_info=list();
		grp_lists=list();
		grp_lists[["covtrt"]]=covtrt;
		grp_lists[["msd"]]=msd;
		grp_lists[["resp"]]=resp;

		for(grp_type in names(grp_lists)){
			grp_link_info[[grp_type]]=add_group_offset(grp_lists[[grp_type]]);
		}

		return(grp_link_info);

	}

	#cat("Results:\n");
	#print(result_rec);

	extr_links=get_group_linkages(result_rec, cvtrt_to_grp_map, cvtrt_grps, msd_grps, rsp_grps);

	cat("Extracted Links:\n");
	print(extr_links);

	# Draw squares and calculate where variables are plotted

	cat("Initializing link locations struct.\n");
	extr_links[["covtrt"]][["pred_grp_members_loc"]]=list();
	extr_links[["covtrt"]][["resp_grp_members_loc"]]=list();
	extr_links[["msd"]][["pred_grp_members_loc"]]=list();
	extr_links[["msd"]][["resp_grp_members_loc"]]=list();
	extr_links[["resp"]][["pred_grp_members_loc"]]=list();
	extr_links[["resp"]][["resp_grp_members_loc"]]=list();

	cat("Drawing squares and calculating link x-positions.\n");
	cat("  Drawing covariates/treatments squares...\n");
	for(i in 1:num_cvtrt_grps){
		if(num_cvtrt_grps==0){break;};

		ichar=sprintf("%i", i);
		var_labels=extr_links[["covtrt"]][["pred_grp_members"]][[ichar]];

		loc=draw_squares_centered(covtrt_xpos, covtrt_ypos[i],
			height=covtrt_height, width=all_widths,
			grp_name=cvtrt_grps[i], 
			variables=var_labels, text_align=ifelse(grp_links, "center", "right"));

		extr_links[["covtrt"]][["pred_grp_members_loc"]][[ichar]]=loc;
	}


	for(i in 1:num_msd_grps){
	cat("  Drawing Measured-to-Measured squares...\n");
		ichar=sprintf("%i", i);

		pred_var_labels=extr_links[["msd"]][["pred_grp_members"]][[ichar]];
		resp_var_labels=extr_links[["msd"]][["resp_grp_members"]][[ichar]];

		covtrt_resp_var_labels=extr_links[["covtrt"]][["resp_grp_members"]][[ichar]];
		resp_pred_var_labels=extr_links[["resp"]][["pred_grp_members"]][[ichar]];

		loc=draw_squares_centered(msd_1_resp_xpos, msd_ypos[i],
			height=msd_height, width=all_widths,
			grp_name=msd_grps[i], 
			variables=covtrt_resp_var_labels, text_align=ifelse(grp_links, "center", "left"));
		extr_links[["covtrt"]][["resp_grp_members_loc"]][[ichar]]=loc;

		loc=draw_squares_centered(msd_1_pred_xpos, msd_ypos[i],
			height=msd_height, width=all_widths,
			grp_name="", 
			variables=pred_var_labels, text_align=ifelse(grp_links, "center", "right"));
		extr_links[["msd"]][["pred_grp_members_loc"]][[ichar]]=loc;

		loc=draw_squares_centered(msd_2_resp_xpos, msd_ypos[i],
			height=msd_height, width=all_widths,
			grp_name=msd_grps[i], 
			variables=resp_var_labels, text_align=ifelse(grp_links, "center", "left"));
		extr_links[["msd"]][["resp_grp_members_loc"]][[ichar]]=loc;

		loc=draw_squares_centered(msd_2_pred_xpos, msd_ypos[i],
			height=msd_height, width=all_widths,
			grp_name="",
			variables=resp_pred_var_labels, text_align=ifelse(grp_links, "center", "right"));
		extr_links[["resp"]][["pred_grp_members_loc"]][[ichar]]=loc;
	}

	cat("  Drawing Response squares...\n");
	for(i in 1:num_rsp_grps){
		if(num_rsp_grps==0){break;};

		ichar=sprintf("%i", i);
		var_labels=extr_links[["resp"]][["resp_grp_members"]][[ichar]];
		loc=draw_squares_centered(rsp_xpos, rsp_ypos[i],
			height=rsp_height, width=all_widths,
			grp_name=rsp_grps[i], 
			variables=var_labels, text_align=ifelse(grp_links, "center", "left"));
		extr_links[["resp"]][["resp_grp_members_loc"]][[ichar]]=loc;
	}

	cat("  Drawing links...\n");
	if(grp_links){
		for(type in names(extr_links)){

			cat("Drawing Group links for: ", type, "\n");
			link_tab=extr_links[[type]][["grp_links"]];
			num_links=nrow(link_tab);
			if(num_links>0){
				for(i in 1:num_links){
					# From covtrt to msd1_resp

					b_grp_off=as.numeric(link_tab[i,"pred"]);
					e_grp_off=as.numeric(link_tab[i,"resp"]);
					b_var_off=as.numeric(link_tab[i,"pred_off"]);
					e_var_off=as.numeric(link_tab[i,"resp_off"]);


					if(type=="covtrt"){
						x_pos=c(covtrt_xpos+edge, msd_1_resp_xpos-edge);
						y_pos=c(covtrt_ypos[b_grp_off], msd_ypos[e_grp_off]);
					}else if(type=="msd"){
						x_pos=c(msd_1_pred_xpos+edge, msd_2_resp_xpos-edge);
						y_pos=c(msd_ypos[b_grp_off], msd_ypos[e_grp_off]);
					}else if(type=="resp"){
						x_pos=c(msd_2_pred_xpos+edge, rsp_xpos-edge);
						y_pos=c(msd_ypos[b_grp_off], rsp_ypos[e_grp_off]);
					}

					points(x=x_pos, y=y_pos, type="l", col="black", lwd=1);

				}
			}
		}
	}else{
		for(type in names(extr_links)){

			cat("Drawing Group links for: ", type, "\n");
			link_tab=extr_links[[type]][["grp_links"]];
			num_links=nrow(link_tab);

			pred_var_grp_loc=extr_links[[type]][["pred_grp_members_loc"]];
			resp_var_grp_loc=extr_links[[type]][["resp_grp_members_loc"]];

			if(num_links>0){
				for(i in 1:num_links){
					# From covtrt to msd1_resp

					pred_grp_off=link_tab[i,"pred"];
					resp_grp_off=link_tab[i,"resp"];
					pred_var_off=link_tab[i,"pred_off"];
					resp_var_off=link_tab[i,"resp_off"];

					link_col=ifelse(link_tab[i,"direction"]=="+", "green", "red");

					pred_grp_off_num=as.numeric(pred_grp_off);
					resp_grp_off_num=as.numeric(resp_grp_off);
					pred_var_off_num=as.numeric(pred_var_off);
					resp_var_off_num=as.numeric(resp_var_off);

					pred_var_loc=pred_var_grp_loc[[pred_grp_off]][pred_var_off_num];
					resp_var_loc=resp_var_grp_loc[[resp_grp_off]][resp_var_off_num];
				

					if(type=="covtrt"){
						x_pos=c(covtrt_xpos+edge, msd_1_resp_xpos-edge);
						y_pos=c(
							pred_var_loc, 
							resp_var_loc);
							#covtrt_ypos[pred_grp_off_num]+pred_var_loc, 
							#msd_ypos[resp_grp_off_num]+resp_var_loc);
					}else if(type=="msd"){
						x_pos=c(msd_1_pred_xpos+edge, msd_2_resp_xpos-edge);
						y_pos=c(
							pred_var_loc, 
							resp_var_loc);
							#msd_ypos[pred_grp_off_num]+pred_var_loc, 
							#msd_ypos[resp_grp_off_num]+resp_var_loc);
					}else if(type=="resp"){
						x_pos=c( msd_2_pred_xpos+edge, rsp_xpos-edge);
						y_pos=c(
							pred_var_loc, 
							resp_var_loc);
							#msd_ypos[pred_grp_off_num]+pred_var_loc, 
							#rsp_ypos[resp_grp_off_num]+resp_var_loc);
					}
					
					points(x=x_pos, y=y_pos, type="l");
					points(x=x_pos, y=y_pos, type="l", col=link_col, lwd=1);
					points(x=x_pos, y=y_pos, type="l", col="black", lwd=.125);

				}
			}
		}
	}

	#----------------------------------------------------------------------

	cat("End of Plot TMR Diagram.\n");
}

remove_weaker_bidirectional_links=function(links_rec, log10_diff_thres=1){

	# Identify MSD to MSD links
	msd_to_msd_ix=links_rec[,"model_type"]=="Msd_to_Msd";

	# Extract out MSD to MSD links, and save other links for later
	msd_to_msd_rec=links_rec[msd_to_msd_ix,,drop=F];
	other_recs=links_rec[!msd_to_msd_ix,,drop=F];

	num_m2m_links=nrow(msd_to_msd_rec);
	if(num_m2m_links>1){

		link_hash=list();
		forward_dir=character();
		opposite_dir=character();

		# Build a hash so we can quickly determine existence of opposite link
		for(i in 1:num_m2m_links){

			grps=strsplit(as.character(msd_to_msd_rec[i, "model_name"]), "->")[[1]];
			pred_grp=grps[1];
			resp_grp=grps[2];

			# Generate group#variable key for pred/resp
			pred_str=paste(pred_grp, "#", msd_to_msd_rec[i, "predictor"], sep="");
			resp_str=paste(resp_grp, "#", msd_to_msd_rec[i, "response"], sep="");

			# Generate pred/resp key
			for_pair_str=paste(pred_str, resp_str, sep="|");
			opp_pair_str=paste(resp_str, pred_str, sep="|");

			# Save the pred/resp pval as store value in hash
			link_hash[[for_pair_str]]=msd_to_msd_rec[i, "pval"];

			# Keep track of pred/resp keys and resp/red keys
			forward_dir=c(forward_dir, for_pair_str);
			opposite_dir=c(opposite_dir, opp_pair_str);

		}

		remove_list=c();
		for(i in 1:num_m2m_links){
		
			cur_for=forward_dir[i];
			cur_opp=opposite_dir[i];	

			if(is.null(link_hash[[cur_opp]])){
				# If there is no, link in opposite direction, do nothing.
				next;
			}else{
				# If there is a link in the opposite direction, compare
				# the p-values.  If one is more significant than the 
				# other by more than the threshold, keep the more significant
				# one.

				cat("F/R link found:", cur_for, "\n");

				for_pval=link_hash[[cur_for]];
				opp_pval=link_hash[[cur_opp]];

				# log_diff < 0, if for_pval more signf than rev_pval
				log10_diff=log10(for_pval/opp_pval);

				if(log10_diff < (-log10_diff_thres)){
					# forward is more significant then opposite
					link_hash[[cur_opp]]=1;
				}else if(log10_diff > (log10_diff_thres)){
					# Opposite is more significant than forward
					# so mark it for removal
					link_hash[[cur_for]]=1;
					msd_to_msd_rec[i, "pval"]=1;
					remove_list=c(remove_list, i);
				}else{
					# Keep both
				}
			}
		}

		# Remove weaker links from table
		msd_to_msd_rec=msd_to_msd_rec[setdiff(1:num_m2m_links, remove_list),, drop=F];
	}

	# Combine filtered MSD to MSD records with other records
	out=rbind(other_recs, msd_to_msd_rec);

	if(nrow(out)>0){
		rownames(out)=1:nrow(out);
	}

	return(out);	

}

plot_title_page("TMR Diagrams", c(
	"The follow diagrams illustrate the group-to-group relationships that have been identified",
	"by the pair-wise regression models.",
	"",
	"The variable groups are organized as squares from left to right (with variables indicated",
	"internally) and lines connecting them, if significant associations between the variables in",
	"the groups had been identified:",
	"",
	"Covariates and Treatments (CT):",
	"Measured (M1):",
	"Measured (M2):",
	"Response (R):",
	"",
	"Note: Measured Groups (M1 and M2) are the same group, but are duplicated in the figure for ease",
	"of interpretation.",
	"When variables from a group on the left (predictors, fit with multiple regression) significantly",
	"predict a variable in a group to the right (response), a line is drawn between them.",
	"",
	"TMR diagrams are drawn at various p-value cutoffs.  For each cutoff, 3 TMR diagrans are drawn:",
	"1.) All Group: To emphasize inter-group relationship, a single line is drawn between groups, even",
	"if multiple variables associations have been identified between groups.",
	"2.) All Variables:  To show all possible relationships, lines are draw between variables of",
	"different groups.",
	"3.) Uni-directional Variables: To show most relevant relationships, when two variables in different",
	"groups in (M1 to M2) predict each other significantly, the less significant of the two links are removed.",
	"",
	"After the TMR diagrams (2) and (3), a table containing a list of the model type, model name, predictor,",
	"response, coefficients and p-values that were represented by the preceding TMR diagram is reported."
));

#for(cutoffs in c("0.0010", "0.1000")){
options(width=300);
for(cutoffs in names(denorm_results)){
	if(cutoffs=="1.0000"){break;}

	plot_TMR_diagram(denorm_results[[cutoffs]],
		paste("P-value Cutoff: ", cutoffs, sep=""),
		"All GROUP Links more significant than cutoff",
		covtrt_to_group_map,
		grp_links=1,
		covariates_list, measured_list, response_list);

	plot_TMR_diagram(denorm_results[[cutoffs]],
		paste("P-value Cutoff: ", cutoffs, sep=""),
		"All VARIABLE Links more significant than cutoff)",
		covtrt_to_group_map,
		grp_links=0,
		covariates_list, measured_list, response_list);

	plot_text(c(
		paste("P-value Cutoff: ", cutoffs, sep=""),
		paste("(All links above cutoff, ", nrow(denorm_results[[cutoffs]]), " links.)", sep=""),
		capture.output(print(denorm_results[[cutoffs]], quotes=""))
	), max_lines_pp=70);

	# Remove weaker of bi-directional links

	unidir_links=remove_weaker_bidirectional_links(denorm_results[[cutoffs]], log10_diff_thres=1);
	plot_TMR_diagram(unidir_links,
		paste("P-value Cutoff: ", cutoffs, sep=""),
		"Uni-directional VARIABLE Links with stronger associations",
		covtrt_to_group_map,
		grp_links=0,
		covariates_list, measured_list, response_list);

	plot_text(c(
		paste("P-value Cutoff: ", cutoffs, sep=""),
		paste("(Excluding weaker of bi-directional links, ", nrow(unidir_links), " links.)", sep=""),
		capture.output(print(unidir_links, quotes=""))
	), max_lines_pp=70);
		
}

##############################################################################
# Export full link table
cat("Exporting link table...\n");

outtable=denorm_results[["1.0000"]];
if(!is.null(RepeatedMeasureAnalysisID)){
	outtable=cbind(RepeatedMeasureAnalysisID, outtable);
}

write.table(
	outtable,
	file=paste(OutputFnameRoot, ".tmr.all_links.tsv", sep=""),
	quote=F, sep="\t", row.names=F, col.names=T
);

##############################################################################

fh=file(paste(OutputFnameRoot, ".tmr.done", sep=""));
cat(file=fh, "Completed.\n");
close(fh);

##############################################################################

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
