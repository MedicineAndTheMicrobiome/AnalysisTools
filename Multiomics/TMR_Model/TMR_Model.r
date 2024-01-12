#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
library(vegan);

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
	"repeated_measure_analysis_id", "a", 2, "character",
	"verbose", "v", 2, "logical"
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
	"	[-v (verbose flag.  Eg. Printing more intermediate tables.)]\n",
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
Verbose=opt$verbose;

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

if(length(opt$verbose)){
	Verbose=T;
}else{
	Verbose=F;
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
	cat("Verbose:                             ", Verbose, "\n");
	cat("\n");
});

cat(paste(param_text, collapse="\n"), "\n");

fh=file(paste(OutputFnameRoot, ".tmr.started", sep=""), "w");
close(fh);

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

	# Remove duplicated variables names
	unique_var=unique(data[,var_col]);
	num_unique=length(unique_var);
	if(num_unique!=nrow(data)){
		cat("Duplicates Found:\n");
		tmp_data=matrix("", nrow=num_unique, ncol=2);
		pos=1;
		for(i in 1:nrow(data)){
			if(all(data[i,var_col]!=tmp_data[,var_col])){
				tmp_data[pos, var_col]=data[i, var_col];
				tmp_data[pos, grp_col]=data[i, grp_col];
				pos=pos+1;
			}
		}
		data=tmp_data;
	}

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
        plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), xaxt="n", yaxt="n", 
		bty="n", xlab="", ylab="");
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
                                text(x-.5, y-.5, text_lab, srt=atan(num_col/num_row)/pi*180, 
					cex=value.cex, font=2);
                        }
                }
        }

        ##################################################################################################

        par(mar=c(0, 0, 0, 0));

        if(plot_row_dendr && plot_col_dendr){
                rdh=attributes(row_dendr[["tree"]])$height;
                cdh=attributes(col_dendr[["tree"]])$height;
                plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", 
			bty="n", xlim=c(rdh, 0));
                plot(col_dendr[["tree"]], leaflab="none",xaxt="n", yaxt="n", 
			bty="n", ylim=c(0, cdh));
                plot(0,0, type="n", bty="n", xaxt="n", yaxt="n");
                #text(0,0, "Placeholder");
        }else if(plot_row_dendr){
                rdh=attributes(row_dendr[["tree"]])$height;
                plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", 
			bty="n", xlim=c(rdh, 0));
                #text(0,0, "Row Dendrogram");
        }else if(plot_col_dendr){
                cdh=attributes(col_dendr[["tree"]])$height;
                plot(col_dendr[["tree"]], leaflab="none", xaxt="n", yaxt="n", 
			bty="n", ylim=c(0, cdh));
                #text(0,0, "Col Dendrogram");
        }

        par(orig.par);

}

###############################################################################

signf_char=function(x_arr){

	xlen=length(x_arr);
	sc=character(xlen);

	for(i in 1:xlen){
		sch="";
		x=x_arr[i];
		if(is.na(x)){
			sch=NA;
		}else{
			if(x<0.001){
				sch="***";
			}else if(x<0.01){
				sch="**";
			}else if(x<0.05){
				sch="*";
			}else if(x<0.10){
				sch="'";
			}
		}
		sc[i]=sch;
	}
	return(sc);
}

###############################################################################

fitsummary_to_list=function(summary, resp_colnames){

	# When the summary of a lm fit is univariate, the summary is not a list.
	# This function places the summary into a list, so univariate responses are not
	# a special case.

	list_names=names(summary);
	num_names=length(summary);

	if(length(grep("^Response ", list_names))==num_names){
		cat("Multivariate Response: \n");
		print(list_names);
		return(summary);
	}else{
		res=list();
		resp_name=resp_colnames[1];
		res[[paste("Response ", resp_name, sep="")]]=summary;
		cat("Univariate Reponse: \n");
		print(names(res));
		return(res);
	}
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

#------------------------------------------------------------------------------

remove_subjects_with_NAs=function(data_matrix, grp_rec){
	all_var=grp_rec$Variables;
	data_matrix=data_matrix[,all_var,drop=F];
	subj_wNA=apply(data_matrix, 1, function(x){any(is.na(x));});
	data_matrix=data_matrix[!subj_wNA,,drop=F];
	rem_rec=list();
	rem_rec[["noNA_matrix"]]=data_matrix;
	rem_rec[["RemovedSubjects"]]=names(subj_wNA[subj_wNA==1]);
	return(rem_rec);
}

#------------------------------------------------------------------------------
# Remove subjects with NAs 

na_removal_rec=remove_subjects_with_NAs(loaded_factors, groupings_rec);
loaded_factors=na_removal_rec[["noNA_matrix"]];

cat("\n");
cat("Subjects Removed Due to NAs:\n");
print(na_removal_rec[["RemovedSubjects"]]);

if(length(na_removal_rec[["RemovedSubjects"]])){
	plot_text(c(
		"Notice:  The TMR Analysis cannot proceed with NAs in the dataset.",
		"The following subjects have NAs in the variables requested for analysis:",
		"",
		capture.output(na_removal_rec[["RemovedSubjects"]]),
		"",
		"They have been removed.  Please confirm you anticipated this."
	));
}

num_samples=nrow(loaded_factors);

variables_rec=list();
variables_rec[["Covariates"]]=split_factors_to_groups(loaded_factors, groupings_rec, covariates_list);
variables_rec[["Measured"]]=split_factors_to_groups(loaded_factors, groupings_rec, measured_list);
variables_rec[["Response"]]=split_factors_to_groups(loaded_factors, groupings_rec, response_list);
#print(variables_rec);

##############################################################################
# Analyze degrees of freedom available based on covariates/treatments and
# measuread group sizes

# Calculate sums for each group and type
count_variables=function(var_rec){

	group_counts_rec=list();
	type_counts_rec=list();

	for(type in names(var_rec)){

		type_counts_rec[[type]]=0;
		group_counts_rec[[type]]=list();

		type_count=0;

		for(grp in names(var_rec[[type]])){

			vars=colnames(var_rec[[type]][[grp]]);

			grp_count=0;
			for(var in vars){
				grp_count=grp_count+1;
			}

			group_counts_rec[[type]][[grp]]=grp_count;
			type_count=type_count+grp_count;
		}
		type_counts_rec[[type]]=type_count;
	}

	counts_rec=list();
	counts_rec[["type"]]=type_counts_rec;
	counts_rec[["groups"]]=group_counts_rec;
	return(counts_rec);
}

counts_rec=count_variables(variables_rec);
#cat("Type Counts:\n");
#print(counts_rec[["type"]]);
#cat("Group Counts:\n");
#print(counts_rec[["groups"]]);

# Generate list of how variables will be used
model_variables_summary=capture.output({
	cat("\n");
	for(type in names(variables_rec)){
		cat(type, " [ ", counts_rec[["type"]][[type]], " ]:\n\n", sep="");
		for(grp in names(variables_rec[[type]])){
			cat("   ", grp, " [ ", counts_rec[["groups"]][[type]][[grp]], " ]:\n", sep="");
			vars=colnames(variables_rec[[type]][[grp]]);
			for(var in vars){
				cat("      ", var, "\n");
			}
			cat("\n");
		}
		cat("\n\n");
	}
});


df_rec=list();
df_min=Inf;
num_df_for_msd=num_samples - counts_rec[["type"]][["Covariates"]] - 1;
for(msd in names(counts_rec[["groups"]][["Measured"]])){
	df_rec[[msd]]=num_df_for_msd-counts_rec[["groups"]][["Measured"]][[msd]];
	df_min=min(df_min, df_rec[[msd]]);
}

df_summary=capture.output({
	cat("Num Samples: ", num_samples, "\n");
	cat("Num Covariates/Treatments: ", counts_rec[["type"]][["Covariates"]], "\n");
	cat("Num DF for Measured Variables (Num Samp - Num Cov/Trt - 1): ", num_df_for_msd, "\n");
	cat("\n");
	cat("Num DF remaining for each Measured Group:\n");
	for(msd in names(counts_rec[["groups"]][["Measured"]])){
		cat("  ", msd, " : ", df_rec[[msd]], "\n", sep="");
	}
	cat("\n");
	cat("Measured Group with Least Remaining DF: ", df_min, "\n");
	cat("\n\n");
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

plot_text(c(
	"Overfit Analysis:\n",
	df_summary
), max_lines_pp=80);

cat(paste(model_variables_summary, collapse="\n"));
cat(paste(df_summary, collapse="\n"));

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
		if(any(!is.finite(st_mat[,var_ix]))){
			st_mat[,var_ix]=rep(0, nrow);
		}
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

# Remove groups with all NAs
na_rows=apply(dist_cor, 1, function(x){ all(is.na(x) | x==1)});
na_cols=apply(dist_cor, 1, function(x){ all(is.na(x) | x==1)});
dist_cor=dist_cor[!na_rows, !na_cols];

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
	min=min(x, na.rm=T);
	max=max(x, na.rm=T);
	range=max-min;
	norm=(x-min)/range;
	return(floor(norm*(steps-1))+1);
}

centers=function(x, steps){
	# Calculate the bin centers in 
	min=min(x, na.rm=T);
	max=max(x, na.rm=T);
	breaks=seq(min, max, length.out=steps+1);
	breaks=round(breaks+(breaks[2]-breaks[1])/2, 3);
	return(head(breaks, steps));
}

max_cat=5;
legends_values=list();
for(var_ix in covariate_variable_names){
	cat("Standardizing: ", var_ix, "\n");
	data=covariates_data[,var_ix];
	colors[,var_ix]=quantize(data, max_cat);
	legends_values[[var_ix]]=centers(data, max_cat);
}

cat("Colors:\n");
print(colors);

cat("Legend Values:\n");
print(legends_values);

palette(rainbow(max_cat, end=4/6));

##############################################################################
# Compute PERMANOVA on measurements and response
doPermanova=F;
if(doPermanova){

plot_title_page("PERMANOVA", c(
	"PERMANOVA was run on the distance matrices for the Measured and",
	"Response groups, using the variables from the Treatment/Covariates",
	"groups as predictors.  R^2 and P-values are reported in tables",
	"and in heatmaps"
));

run_adonis=function(in_dist, in_factors, num_perm_per_var=1000){

	cat("\n");
	cat("Number of Variables: ", ncol(in_factors), "\n");
	cat("Number of Samples: ", nrow(in_factors), "\n");

	formula_str=paste("in_dist ~", paste(colnames(in_factors), collapse=" + "));

	cat("Running Adonis: \n", formula_str, "\n");

	num_perm_per_var=50;

	adon2_res=adonis2(
		as.formula(formula_str), 
		data=as.data.frame(in_factors),
		permutations=ncol(in_factors)*num_perm_per_var
		);
	anova_tab=as.matrix(adon2_res);
	return(anova_tab);
}

permanova_rec=list();
for(type in c("Measured", "Response")){
	for(grp in names(mds_rec[[type]])){
		cat("Working on: ", type, " / ", grp, "\n");
		distances=dist_rec[[type]][[grp]];
		perm_mod_name=paste(type, ": ", grp, sep="");

		permanova_rec[[perm_mod_name]]=run_adonis(distances, 
			covariates_data[,covariate_variable_names, drop=F]);

		print(permanova_rec[[perm_mod_name]]);
	}
}


# Consolidate PERMANOVA records into matrices
num_perm_models=length(permanova_rec);
perm_models=names(permanova_rec);

permanova_pval_mat=matrix(NA, nrow=length(covariate_variable_names), ncol=num_perm_models);
permanova_rsqr_mat=matrix(NA, nrow=length(covariate_variable_names), ncol=num_perm_models);

rownames(permanova_pval_mat)=covariate_variable_names;
rownames(permanova_rsqr_mat)=covariate_variable_names;
colnames(permanova_pval_mat)=perm_models;
colnames(permanova_rsqr_mat)=perm_models;

for(perm_grp in perm_models){

	perm_mat=permanova_rec[[perm_grp]];

	#print(perm_mat[,"Pr(>F)"]);
	avail_var=rownames(perm_mat);
	shared_var=intersect(avail_var, covariate_variable_names);
	num_shared_var=length(shared_var);

	if(num_shared_var>0){
		covariate_variable_names=shared_var;
	}else{
		next;
	}
	
	permanova_pval_mat[covariate_variable_names, perm_grp]=
		perm_mat[covariate_variable_names, "Pr(>F)"];
	permanova_rsqr_mat[covariate_variable_names, perm_grp]=
		perm_mat[covariate_variable_names, "R2"];
}

# Plot heatmaps for R^2 and P-values
paint_matrix(permanova_pval_mat, title="PERMANOVA P-Values", plot_min=0, plot_max=1, high_is_hot=F);
paint_matrix(permanova_rsqr_mat, title="PERMANOVA R^2", plot_min=0, plot_max=1, high_is_hot=T);

# Plot summary tables
plot_text(c(
	"PERMANOVA Results:",
	"",
	"P-values:",
	"",
	capture.output(print(permanova_pval_mat)),
	"",
	"",
	"R^2:",
	capture.output(print(permanova_rsqr_mat))
	), max_lines_pp=80);



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
	"Treatment/Covariate group.",
	"",
	"Plots have been annotated with PERMANOVA calculation from previous section."
));


par(mar=c(4,4,6,1));
for(type in c("Measured", "Response")){

	for(grp in names(mds_rec[[type]])){
		mds_coord=mds_rec[[type]][[grp]];		
		distances=dist_rec[[type]][[grp]];
	
		for(var_ix in covariate_variable_names){

			plot(mds_coord[,1], mds_coord[,2], col=colors[, var_ix],
				xlab="Dim 1", ylab="Dim 2", type="n",
				main=paste("Group: ", grp, sep=""),
				);
			title(main=paste("(Colored by: ", var_ix, ")", sep=""), line=1.5, cex.main=.95);

			# Label with permanova results
			perm_mod_name=paste(type, ": ", grp, sep="");
			perm_r2=round(permanova_rsqr_mat[var_ix, perm_mod_name], 4);
			perm_pv=paste(
				round(permanova_pval_mat[var_ix, perm_mod_name], 4),
				" ", signf_char(permanova_pval_mat[var_ix, perm_mod_name]),
				sep=""
				);
			title(main=paste("R^2 = ", perm_r2, " / P-val = ", perm_pv, sep=""), 
				line=.25, cex.main=.9);
			
			# Plot points and sample labels
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

}

#*******************************************

pick_and_fit_model=function(mod_var, pred_val, resp_val){

	if(0){
		cat("Picking and Fitting Models:\n");
		print(mod_str);
		cat("Num Samples:\n");
		print(nrow(pred_val));
		cat("Predictors:\n");
		print(colnames(pred_val));
		cat("Responses:\n");
		print(colnames(resp_val));
	}

	# Find number of responses, instead of using multivariate responses feature of lm
	num_responses=ncol(resp_val);
	resp_names=colnames(resp_val);
	cat("Number of Responses: ", num_responses, "\n");

	# Determine whether response is normal or binomial, or other?
	model_families=character(num_responses);
	names(model_families)=resp_names;
	for(i in 1:num_responses){
		cur_resp=resp_val[,i];
		uniq_resp=unique(cur_resp);
		num_uniq_resp=length(uniq_resp);
		min_resp=min(uniq_resp);
		max_resp=max(uniq_resp);

		if(num_uniq_resp==2 && min_resp==0 && max_resp==1){
			model_families[i]="binomial";	
		}else{
			model_families[i]="gaussian";
		}
	}

	cat("Model Families:\n");
	print(model_families);	
	
	# Old fitting implementation
	#fit=lm(as.formula(mod_str), data=pred_val);
	#sum_fit=fitsummary_to_list(summary(fit), resp_val);
	insufficient_residual_dfs=F;

	fit_list=list();
	sum_fit_list=list();
	ctl_lst=list(maxit=100);

	for(i in 1:num_responses){
		
		# Build formula from model string and current response
		resp_name=resp_names[i];
		y=resp_val[,i];
		mod_str=paste(mod_var, collapse="+");
		full_mod_str=paste("y ~ ", mod_str);
		full_mod_form=as.formula(full_mod_str);

	
		# Change the family based previously detected model family
		if(model_families[i]=="gaussian"){
			fit=glm(full_mod_form, data=pred_val, family=gaussian, control=ctl_lst);
		}else if(model_families[i]=="binomial"){
			fit=glm(full_mod_form, data=pred_val, family=binomial, control=ctl_lst)
		}else{
			cat("Error:  Unspecified Model Family.\n");
			quit(status=-1);
		}

	
		# Summarize (calculate p-values);
		sumfit=summary(fit);

		# Flag over parameterized models
		if(sumfit[["df.residual"]]==0){
			insufficient_residual_dfs=T;
			break;
		}

		# If gaussian, then Pr(>|t|), if binomial then Pr(>|z|)
		# Unify, so downstream column names are the same
		cnames=colnames(sumfit[["coefficients"]]);
		if(grep(cnames[4], "Pr(>|.|)")){
			cnames[4]="Pr(>|.|)";
		}
		colnames(sumfit[["coefficients"]])=cnames;

		# Save for export
		sum_fit_list[[resp_name]]=sumfit;
		fit_list[[resp_name]]=fit;

	}

	results=list();
	results[["predictors"]]=mod_var;
	results[["fit_list"]]=fit_list;
	results[["sumfit_list"]]=sum_fit_list;
	results[["sufficient_residuals"]]=!insufficient_residual_dfs;

	return(results);
}

#*******************************************

##############################################################################
# Fit Covariates to Predict Measured

model_results=list();

model_results[["Cov_to_Msd"]]=list();
model_results[["Msd_to_Msd"]]=list();
model_results[["Msd_to_Rsp"]]=list();

overfits=list();

overfits[["Cov_to_Msd"]]=list();
overfits[["Msd_to_Msd"]]=list();
overfits[["Msd_to_Rsp"]]=list();

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

	model_vars=covariate_variable_names;
	#fit=lm(as.formula(model_string), data=covariates_data);
	#sum_fit=fitsummary_to_list(summary(fit), msd_resp);
	model_fit=pick_and_fit_model(model_vars, covariates_data, msd_resp);
	sum_fit=model_fit[["sumfit_list"]];
	overfits[["Cov_to_Msd"]][[msd_ix]]=!model_fit[["sufficient_residuals"]];

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
			sum_fit[[var_ix]][["coefficients"]][avail_pred,"Pr(>|.|)"];

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
	
	any_overfits=F;

	for(resp_msd_ix in measured_list){

		if(pred_msd_ix==resp_msd_ix){
			next;
		}

		cat("\n\nFitting: Measured to Measured: (Pred)", pred_msd_ix, " (Resp)", 
			resp_msd_ix, "\n");

		analysis_string=paste(pred_msd_ix, "->", resp_msd_ix, sep="");

		msd_resp=as.matrix(variables_rec[["Measured"]][[resp_msd_ix]]);
		num_resp=ncol(msd_resp);
		msd_resp_varnames=colnames(msd_resp);

		msd_pred=as.matrix(variables_rec[["Measured"]][[pred_msd_ix]]);
		num_pred=ncol(msd_pred);
		msd_pred_varnames=colnames(msd_pred);

		#model_string=paste(paste(c(covariate_variable_names, msd_pred_varnames), collapse=" + "));
		#model_string=paste("msd_resp ~ ", 
		#	paste(c(covariate_variable_names, msd_pred_varnames), collapse=" + "));

		#cat("Model: \n");
		#print(model_string);	

		#fit=lm(as.formula(model_string), data=cbind(covariates_data, msd_pred));
		#sum_fit=fitsummary_to_list(summary(fit), msd_resp_varnames);
		model_vars=c(covariate_variable_names, msd_pred_varnames);

		model_fit=pick_and_fit_model(model_vars, cbind(covariates_data, msd_pred), 
			msd_resp);
		sum_fit=model_fit[["sumfit_list"]];

		any_overfits = any_overfits || !model_fit[["sufficient_residuals"]];

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
				sum_fit[[var_ix]][["coefficients"]][avail_pred,"Pr(>|.|)"];


			coef_mat[avail_pred, varname]=
				sum_fit[[var_ix]][["coefficients"]][avail_pred,"Estimate"];

		}

		# Store in record
		model_results[["Msd_to_Msd"]][[analysis_string]][["pval"]]=pval_mat;
		model_results[["Msd_to_Msd"]][[analysis_string]][["coef"]]=coef_mat;

	}

	overfits[["Msd_to_Msd"]][[pred_msd_ix]]=any_overfits;

}

#print(model_results[["Msd_to_Msd"]]);

cat("\n\n");

##############################################################################
# Fit (Covariates + Measured) to Predict Response

for(pred_msd_ix in measured_list){
	
	any_overfits=F;

	for(resp_ix in response_list){

		cat("Fitting: Measured to (as Predictor)", pred_msd_ix, " to ", 
			resp_ix, " (as Response)\n");

		analysis_string=paste(pred_msd_ix, "->", resp_ix, sep="");

		resp=as.matrix(variables_rec[["Response"]][[resp_ix]]);
		num_resp=ncol(resp);
		resp_varnames=colnames(resp);

		msd_pred=as.matrix(variables_rec[["Measured"]][[pred_msd_ix]]);
		num_pred=ncol(msd_pred);
		msd_pred_varnames=colnames(msd_pred);

		#model_string=paste(paste(c(covariate_variable_names, msd_pred_varnames), collapse=" + "));
		#model_string=paste("resp ~ ", 
		#	paste(c(covariate_variable_names, msd_pred_varnames), collapse=" + "));

		#cat("Model: \n");
		#print(model_string);	

		#fit=lm(as.formula(model_string), data=cbind(covariates_data, msd_pred));
		#sum_fit=fitsummary_to_list(summary(fit), resp_varnames);
		#print(sum_fit);

		model_vars=c(covariate_variable_names, msd_pred_varnames);
		model_fit=pick_and_fit_model(model_vars, cbind(covariates_data, msd_pred),
                        resp);
		sum_fit=model_fit[["sumfit_list"]];
		any_overfits = any_overfits || !model_fit[["sufficient_residuals"]];

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
			cat("Moving: ", var_ix, "\n");
			print(sum_fit[[var_ix]][["coefficients"]]);

			# Sometimes NAs in fit drop the predictor name from coefficients table
			avail_pred=intersect(
				rownames(sum_fit[[var_ix]][["coefficients"]]),
				cov_and_pred_names
			);

			pval_mat[avail_pred, varname]=
				sum_fit[[var_ix]][["coefficients"]][avail_pred,"Pr(>|.|)"];

			coef_mat[avail_pred, varname]=
				sum_fit[[var_ix]][["coefficients"]][avail_pred,"Estimate"];

		}

		# Store in record
		model_results[["Msd_to_Rsp"]][[analysis_string]][["pval"]]=pval_mat;
		model_results[["Msd_to_Rsp"]][[analysis_string]][["coef"]]=coef_mat;

	}

	overfits[["Msd_to_Rsp"]][[pred_msd_ix]]=any_overfits;
}

#print(model_results[["Msd_to_Rsp"]]);
#print(overfits);

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

#------------------------------------------------------------------------------

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

#------------------------------------------------------------------------------

remove_covariate_links_from_msd_pred=function(links_rec, covariates){
	#print(covariates);
	#print(links_rec);

	is_msdtomsd=("Msd_to_Msd"==links_rec[,"model_type"]);
	is_covar=links_rec[,"predictor"] %in% covariates;
	is_msdtomsd_and_covar=is_msdtomsd & is_covar;

	covar_removed=links_rec[!is_msdtomsd_and_covar,,drop=F];
	return(covar_removed);
}

##############################################################################

#print(model_results);

denorm_results=list();
unidir_results=list();

numerical_cutoffs=sort(c(1.0000, 0.1000, 0.050, 0.010, 0.005, 0.001, 0.0005, 0.0001), decreasing=T);
string_cutoffs=sprintf("%3.4f", numerical_cutoffs); 
num_cutoffs=length(numerical_cutoffs);

for(i in 1:num_cutoffs){

	denorm_results[[string_cutoffs[i]]]=
		matrix_to_tables(model_results, numerical_cutoffs[i]);

	unidir_results[[string_cutoffs[i]]]=
		remove_weaker_bidirectional_links(
			denorm_results[[string_cutoffs[i]]], 
			log10_diff_thres=1);
}

##############################################################################

extract_predictors=function(unidir_dnrm_table, model_type){
	#print(unidir_dnrm_table);
	#print(model_type);	
	cat("Extracting predictors for: ", model_type, "\n");

	model_rows_ix=unidir_dnrm_table[,"model_type"]==model_type;
	model_table=unidir_dnrm_table[model_rows_ix,,drop=F];

	predictors=sort(as.character(unique(model_table[,"predictor"])));
	print(predictors);

	return(predictors);
}

cat("=======================================================\n");
cat("=  Fitting Selected Predictors to Response Variables  =\n");
cat("=======================================================\n");

selected_predictors=list();
selected_fits=list();

for(pvco in string_cutoffs){

	selected_fits[[pvco]]=list();	
	selected_predictors[[pvco]]=list();

	cat("\n\n");
	cat("**** P-value Cutoff: ", pvco, "\n");	

	for(model_type in c("Cov_to_Msd", "Msd_to_Msd", "Msd_to_Rsp")){


		cat("**** Model Type: ", model_type, "\n");

		sel_pred=extract_predictors(unidir_results[[pvco]], model_type);
		selected_predictors[[pvco]][[model_type]]=sel_pred;

		if(length(sel_pred)){

			cat("**** Selected Predictors across Groups:\n");
			model_vars=sel_pred;

			for(resp_ix in response_list){

				# Fits multivariate response
				cat("**** Fitting Response Group: ", resp_ix, "\n");
				model_fit=pick_and_fit_model(
					model_vars,
					loaded_factors[,sel_pred,drop=F],
					variables_rec[["Response"]][[resp_ix]]
					);

				#cat("---------------------->>\n");
				#print(model_fit);
				#cat("<<----------------------\n");

				selected_fits[[pvco]][[model_type]]=model_fit;

			}
		}else{
			cat("No predictors found for: ", model_type, " at ", pvco, "\n");
			selected_fits[[pvco]][[model_type]]=list();
		}
			
	}
}

cat("============================================================\n");
cat("=  Done Fitting Selected Predictors to Response Variables  =\n");
cat("============================================================\n");

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

draw_squares_centered=function(xpos, ypos, height, width, grp_name, variables, text_align, 
	overfit_hash=NULL){

	print(c(xpos, ypos, height, width, grp_name, variables, text_align));
	
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

	if(!is.null(overfit_hash) && overfit_hash==1){
		text(xpos, ypos, "O", col="red", cex=2.5, adj=c(.5,.5), font=1);
		text(xpos, ypos, "overfit", col="red", cex=.5, adj=c(.5, 4), font=3);
	}

	title_cex=.7;
	var_cex=.4;

	text(xpos, ypos+height/2, grp_name, font=2, cex=title_cex, pos=1);
	#text(xpos, ypos+height/2, paste(c(rep("",3), variables), collapse="\n"), cex=.4, pos=1);
	
	title_spc=par()$cxy[2]*title_cex*1.7;
	var_spc=par()$cxy[2]*var_cex;
	
	num_var=length(variables);

	if(num_var>0){
		text_pos=calc_vertical_spacing(num_var, start=ypos+height/2-title_spc, 
			end=ypos-height/2+var_spc);


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
	cvtrt_grps, msd_grps, rsp_grps, overfit_list){

	cat("Plotting TMR diagram: ", title, "\n");
	cat("Subtitle: ", subtitle, "\n");

	num_cvtrt_grps=length(cvtrt_grps);
	num_msd_grps=length(msd_grps);
	num_rsp_grps=length(rsp_grps);


	# Reverse order so groups go from top to bottom
	cvtrt_grps=rev(cvtrt_grps);
	msd_grps=rev(msd_grps);
	rsp_grps=rev(rsp_grps);

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
			variables=var_labels, text_align=ifelse(grp_links, "center", "right"),
			overfit_hash=overfit_list[["Cov_to_Msd"]][[cvtrt_grps[i]]]);

		extr_links[["covtrt"]][["pred_grp_members_loc"]][[ichar]]=loc;
	}

	cat("  Drawing Measured-to-Measured squares...\n");
	for(i in 1:num_msd_grps){
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
			variables=pred_var_labels, text_align=ifelse(grp_links, "center", "right"),
			overfit_hash=overfit_list[["Msd_to_Msd"]][[msd_grps[i]]]);
		extr_links[["msd"]][["pred_grp_members_loc"]][[ichar]]=loc;

		loc=draw_squares_centered(msd_2_resp_xpos, msd_ypos[i],
			height=msd_height, width=all_widths,
			grp_name=msd_grps[i], 
			variables=resp_var_labels, text_align=ifelse(grp_links, "center", "left"));
		extr_links[["msd"]][["resp_grp_members_loc"]][[ichar]]=loc;

		loc=draw_squares_centered(msd_2_pred_xpos, msd_ypos[i],
			height=msd_height, width=all_widths,
			grp_name="",
			variables=resp_pred_var_labels, text_align=ifelse(grp_links, "center", "right"),
			overfit_hash=overfit_list[["Msd_to_Rsp"]][[msd_grps[i]]]);
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


extract_matrices_from_links=function(no_trtcov_unidir_links,
		covtrt_to_group_map, covariates_list, measured_list, response_list
){

	#print(covtrt_to_group_map);
	#print(covariates_list);
	#print(measured_list);
	#print(response_list);

	print(no_trtcov_unidir_links);

	cov_to_msd_tab=no_trtcov_unidir_links[
		no_trtcov_unidir_links[,"model_type"]=="Cov_to_Msd",,drop=F];
	msd_to_msd_tab=no_trtcov_unidir_links[
		no_trtcov_unidir_links[,"model_type"]=="Msd_to_Msd",,drop=F];
	msd_to_rsp_tab=no_trtcov_unidir_links[
		no_trtcov_unidir_links[,"model_type"]=="Msd_to_Rsp",,drop=F];


	# Put predictors on rows, responses on columns
	tab_to_mat=function(table){

		uniq_pred=unique(table[,"predictor"]);
		uniq_resp=unique(table[,"response"]);

		num_uniq_pred=length(uniq_pred);
		num_uniq_resp=length(uniq_resp);
		mat=matrix(0, nrow=num_uniq_pred, ncol=num_uniq_resp);
		rownames(mat)=uniq_pred;
		colnames(mat)=uniq_resp;

		nrow=nrow(table);
		cat("Num Rows: ", nrow, "\n");
		if(nrow>0){
			for(i in 1:nrow){
				rn=as.character(table[i,"predictor"]);
				cn=as.character(table[i,"response"]);
				coef=table[i,"coef"];

				if(!is.na(mat[rn,cn])){
					if(mat[rn,cn] >= 0  && coef > 0){
						mat[rn,cn] = 1;
					}else if(mat[rn,cn] <= 0  && coef < 0){
						mat[rn,cn] = -1;
					}else{
						mat[rn,cn] = NA;
					}
				}
			}
		}else{
			mat=matrix(NA, nrow=0, ncol=0);
		}

		return(mat);
	}

	results=list();
	results[["c2m"]]=tab_to_mat(cov_to_msd_tab);
	results[["m2m"]]=tab_to_mat(msd_to_msd_tab);
	results[["m2r"]]=tab_to_mat(msd_to_rsp_tab);

	return(results);
}

tmr_heatmap=function(mat, title="", subtitle="", pred_var_mat, resp_var_mat, value.cex=1){

	cat("Generating TMR Heatmap for: ", title, "\n", sep="");

        num_row=nrow(mat);
        num_col=ncol(mat);

	cat("Mat Dimensions:", num_row, " x ", num_col, "\n");
	print(mat);

        row_names=rownames(mat);
        col_names=colnames(mat);

	#----------------------------------------------------------------------

	subset_mat=function(targets, ordermat){
		# extract targets from the ordermat 

		if(length(targets)==0 || nrow(ordermat)==0){
			return(ordermat[c(),]);
		}

		ordermat=as.data.frame(ordermat);
		om_cname=colnames(ordermat);
		ordermat=cbind(ordermat,F);
		colnames(ordermat)=c(om_cname, "Keep");
		rownames(ordermat)=ordermat[,"Variable"];

		for(t in targets){
			ordermat[t,3]=T;
		}
		
		ordermat=ordermat[ordermat[,"Keep"],c(1,2),drop=F];
		ordermat=as.matrix(ordermat);

		rownames(ordermat)=c();
		return(ordermat);
	}

	#----------------------------------------------------------------------

	find_group_sizes=function(grp_lst){
		# Find group counts while preserving their order

		num_tot_items=length(grp_lst);

		item_lst=c();
		cur_item=grp_lst[1];
		item_lst=cur_item;
		if(num_tot_items>0){
			for(i in 1:num_tot_items){
				if(cur_item == grp_lst[i]){
					next;
				}else{
					cur_item=grp_lst[i];
					item_lst=c(item_lst, cur_item);
				}
			}
		}

		num_uniq_items=length(item_lst);

		item_counts=numeric(num_uniq_items);
		names(item_counts)=item_lst;
		if(num_uniq_items>0){
			for(i in 1:num_uniq_items){
				item_counts[i]=sum(grp_lst==item_lst[i]);
			}
		}
		
		return(item_counts);

	}

	#----------------------------------------------------------------------


        orig.par=par(no.readonly=T);

        cat("TMR Heatmap: ", title, "\n");
        cat("Num Rows: ", num_row, "\n");
        cat("Num Cols: ", num_col, "\n");

	pred_var_mat=subset_mat(row_names, pred_var_mat);
	resp_var_mat=subset_mat(col_names, resp_var_mat);

	cat("Predictor Mat:\n");
	print(pred_var_mat);
	cat("\n");
	
	cat("Response Mat:\n");
	print(resp_var_mat);
	cat("\n");

        if(nrow(pred_var_mat)==0 || nrow(resp_var_mat)==0){
		par(mar=c(5,5,5,5));
                plot(0, type="n", xlim=c(-1,1), ylim=c(-1,1), xaxt="n", yaxt="n", 
			bty="n", xlab="", ylab="");
                text(0,0, "No relationships to plot...");
		mtext(title, side=3, line=1.5, outer=F, font=2, cex=2);
		mtext(paste("P-value Cutoff: ", subtitle, sep=""), side=3, line=.5, 
			outer=F, font=1, cex=1);

                return();
        }


	row_names=pred_var_mat[,"Variable",drop=F];
	col_names=resp_var_mat[,"Variable",drop=F];

	pred_grp_sizes=find_group_sizes(pred_var_mat[,"Group"]);
	resp_grp_sizes=find_group_sizes(resp_var_mat[,"Group"]);

	data_mat=mat[row_names, col_names, drop=F];

        # Get Label lengths
        row_max_nchar=max(nchar(row_names));
        col_max_nchar=max(nchar(col_names));
        cat("Max Row Names Length: ", row_max_nchar, "\n");
        cat("Max Col Names Length: ", col_max_nchar, "\n");

	plot_max=max(mat);

        ########################################################################

	par(oma=c(2,2,0,0));
        par(mar=c(col_max_nchar*.40, row_max_nchar*.40, 5, col_max_nchar*.1));
        plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), xaxt="n", yaxt="n", 
		bty="n", xlab="", ylab="");
        mtext(title, side=3, line=1.5, outer=F, font=2, cex=2);
	mtext(paste("P-value Cutoff: ", subtitle, sep=""), side=3, line=.5, 
		outer=F, font=1, cex=1);

        for(x in 1:num_col){
                for(y in 1:num_row){

                        cell_val=data_mat[(num_row-y+1),x];

                       
			#remap_val=remap(cell_val, c(0, plot_max), c(1, num_colors));
                        #col_ix=ceiling(remap_val);
	
			if(is.na(cell_val)){
				cell_col="purple";
			}else if(cell_val>0){
				cell_col="green";
			}else if(cell_val<0){
				cell_col="red";
			}else{
				cell_col="grey85";
			}


			# Draw/Color the cell
                        rect(x-1, y-1, (x-1)+1, (y-1)+1, border=NA, col=cell_col);

			# Label the counts
                        if(cell_val>1){
                                text_lab=sprintf("%i", cell_val);
                                text(x-.5, y-.5, text_lab, srt=atan(num_col/num_row)/pi*180, 
					cex=value.cex, font=2);
                        }
                }
        }

	# Draw Grid Lines
	abline(h=0:num_row, lwd=.5, col="grey");
	abline(v=0:num_col, lwd=.5, col="grey");

	# Calculate the size of the labels
	txt_w=par()$cxy[1];
	txt_h=par()$cxy[2];
	label_size_x=min(c(1, 40/num_col));
	label_size_y=min(c(1, 25/num_row));

	# y-axis: predictor
        axis(side=2, at=seq(.5, num_row-.5, 1), labels=rev(row_names), las=2, line=-1.75,
		cex.axis=label_size_y, tick=F);

	# x-axis: response
        text((1:num_col)-1/2- label_size_x*txt_w*1.5, rep(-txt_h/2, num_col), col_names, 
		srt=-60, xpd=T, pos=4, cex=label_size_x);


	# Draw Group Lines
	pred_grp_lines=c(0, cumsum(rev(pred_grp_sizes)));
	resp_grp_lines=c(0, cumsum(resp_grp_sizes));

	for(i in 1:length(resp_grp_lines)){
		points(c(resp_grp_lines[i], resp_grp_lines[i]), c(0, num_row), type="l", 
			lwd=1.5, col="blue");
	}

	for(i in 1:length(pred_grp_lines)){
		points(c(0, num_col), c(pred_grp_lines[i], pred_grp_lines[i]), type="l",
			lwd=1.5, col="blue");
	}

        ########################################################################

	# Label Pred/Resp in margins
	mtext("Predictor Variables", side=2, outer=T, font=3, cex=2, col="grey");
	mtext("Response Variables", side=1, outer=T, font=3, cex=2, col="grey");

        par(orig.par);

}

################################################################################

tmr_heatmap_byGroup=function(mat, title="", subtitle="", pred_var_mat, resp_var_mat, value.cex=1){

	cat("Generating TMR Group Heatmap for: ", title, "\n", sep="");
        orig.par=par(no.readonly=T);

        num_row=nrow(mat);
        num_col=ncol(mat);

	cat("Mat Dimensions:", num_row, " x ", num_col, "\n");
	print(mat);

        row_names=rownames(mat);
        col_names=colnames(mat);

	#----------------------------------------------------------------------

	subset_mat=function(targets, ordermat){
		# Retrun the subset in the right order.
		# The targets may not be in the right order but the entries
		# in the ordermat will be.

		if(length(targets)==0 || nrow(ordermat)==0){
			return(ordermat[c(),]);
		}

		ordermat=as.data.frame(ordermat);
		om_cname=colnames(ordermat);
		ordermat=cbind(ordermat,F);
		colnames(ordermat)=c(om_cname, "Keep");
		rownames(ordermat)=ordermat[,"Variable"];

		for(t in targets){
			ordermat[t,3]=T;
		}
		
		ordermat=ordermat[ordermat[,"Keep"],c(1,2),drop=F];
		ordermat=as.matrix(ordermat);

		rownames(ordermat)=c();
		return(ordermat);
	}

	#----------------------------------------------------------------------

	pred_var_mat=subset_mat(row_names, pred_var_mat);
	resp_var_mat=subset_mat(col_names, resp_var_mat);

	cat("Predictor Mat:\n");
	print(pred_var_mat);
	cat("\n");
	
	cat("Response Mat:\n");
	print(resp_var_mat);
	cat("\n");

	pred_groups=unique(pred_var_mat[,"Group"]);
	resp_groups=unique(resp_var_mat[,"Group"]);

	num_unq_pred_grps=length(pred_groups);
	num_unq_resp_grps=length(resp_groups);

	grp_count_mat=matrix(0, nrow=num_unq_pred_grps, ncol=num_unq_resp_grps);
	rownames(grp_count_mat)=pred_groups;
	colnames(grp_count_mat)=resp_groups;

	if(num_unq_pred_grps>0 && num_unq_resp_grps>0){
		for(i in 1:num_unq_pred_grps){
			for(j in 1:num_unq_resp_grps){

				pred_grp=pred_groups[i];
				resp_grp=resp_groups[j];

				pred_vars=pred_var_mat[(pred_var_mat[,"Group"]==pred_grp), "Variable"];
				resp_vars=resp_var_mat[(resp_var_mat[,"Group"]==resp_grp), "Variable"];

				var_sub_mat=mat[pred_vars, resp_vars, drop=F];
			
				num_assoc=sum(abs(var_sub_mat));
				grp_count_mat[pred_grp, resp_grp]=num_assoc;

			}
		}
	}

	print(grp_count_mat);

        cat("TMR Group Heatmap: ", title, "\n");
        cat("Num Rows: ", num_unq_pred_grps, "\n");
        cat("Num Cols: ", num_unq_resp_grps, "\n");

        if(num_unq_pred_grps==0 || num_unq_resp_grps==0){
		par(mar=c(5,5,5,5));
                plot(0, type="n", xlim=c(-1,1), ylim=c(-1,1), xaxt="n", yaxt="n", 
			bty="n", xlab="", ylab="");
                text(0,0, "No relationships to plot...");
		mtext(title, side=3, line=1.5, outer=F, font=2, cex=2);

                return();
        }

        # Generate a color scheme
        num_colors=50;
        color_arr=rainbow(num_colors, start=0, end=4/6);
        color_arr=rev(color_arr);
	color_arr[1]="#FFFFFF";

        # Provide a means to map values to an (color) index
        remap=function(in_val, in_range, out_range){
                in_prop=(in_val-in_range[1])/(in_range[2]-in_range[1])
                out_val=in_prop*(out_range[2]-out_range[1])+out_range[1];
                return(out_val);
        }

        # Get Label lengths
        row_max_nchar=max(nchar(pred_groups));
        col_max_nchar=max(nchar(resp_groups));
        cat("Max Row Names Length: ", row_max_nchar, "\n");
        cat("Max Col Names Length: ", col_max_nchar, "\n");

	plot_max=max(grp_count_mat);

        ########################################################################

	par(oma=c(2,2,0,0));
        par(mar=c(col_max_nchar*.40, row_max_nchar*.40, 5, col_max_nchar*.1));
        plot(0, type="n", xlim=c(0,num_unq_resp_grps), ylim=c(0,num_unq_pred_grps), xaxt="n", yaxt="n", 
		bty="n", xlab="", ylab="");
        mtext(title, side=3, line=1.5, outer=F, font=2, cex=2);
        mtext(paste("P-value Cutoff: ", subtitle, sep=""), side=3, line=0.4, outer=F, font=1, cex=1);


	txt_w=par()$cxy[1];
	txt_h=par()$cxy[2];
	label_size_x=min(c(1, 40/num_unq_resp_grps));
	label_size_y=min(c(1, 25/num_unq_pred_grps));

	# y-axis: predictor
        axis(side=2, at=seq(.5, num_unq_pred_grps-.5, 1), labels=rev(pred_groups), las=2, line=-1.75,
		cex.axis=label_size_y, tick=F);

        # x-axis: response
        text((1:num_unq_resp_grps)-1/2- label_size_x*txt_w*1.5, rep(-txt_h/2, num_unq_resp_grps), 
		resp_groups, srt=-60, xpd=T, pos=4, cex=label_size_x);


        for(x in 1:num_unq_resp_grps){
                for(y in 1:num_unq_pred_grps){

                        cell_val=grp_count_mat[(num_unq_pred_grps-y+1),x];

                        remap_val=remap(cell_val, c(0, plot_max), c(1, num_colors));
                        col_ix=ceiling(remap_val);

			# Draw/Color the cell
                        rect(x-1, y-1, (x-1)+1, (y-1)+1, border=NA, col=color_arr[col_ix]);

			# Label the counts
                        if(cell_val>0){
                                text_lab=sprintf("%i", cell_val);
                                text(x-.5, y-.5, text_lab, srt=atan(num_col/num_row)/pi*180, 
					cex=value.cex, font=2);
                        }
                }
        }

	# Draw Grid Lines
	abline(h=0:num_unq_pred_grps, lwd=2, col="blue");
	abline(v=0:num_unq_resp_grps, lwd=2, col="blue");

	# Label Pred/Resp in margins
	mtext("Predictor Groups", side=2, outer=T, font=3, cex=2, col="grey");
	mtext("Response Groups", side=1, outer=T, font=3, cex=2, col="grey");

        ########################################################################

        par(orig.par);

	return(grp_count_mat);

}	

marginal_stacked_barplots_byGroup=function(t2m, m2m, m2r, grp_colors, cutoff){

	summag=function(x){
		x=x[!is.na(x)];
		sm=sum(abs(x));
		return(sm);
	}

	if(!is.null(t2m)){
		pred_t2m=apply(t2m, 1, summag);
		resp_t2m=apply(t2m, 2, summag);
	}else{
		pred_t2m=0;
		resp_t2m=0;
	}

	if(!is.null(m2m)){
		pred_m2m=apply(m2m, 1, summag);
		resp_m2m=apply(m2m, 2, summag);
	}else{
		pred_m2m=0;
		resp_m2m=0;
	}

	if(!is.null(m2r)){
		pred_m2r=apply(m2r, 1, summag);
		resp_m2r=apply(m2r, 2, summag);
	}else{
		pred_m2r=0;
		resp_m2r=0;
	}


	lengths=c(
		pred_t2m=length(pred_t2m),
		resp_t2m=length(resp_t2m),
		pred_m2m=length(pred_m2m),
		resp_m2m=length(resp_m2m),
		pred_m2r=length(pred_m2r),
		resp_m2r=length(resp_m2r)
		);

	sums=c(
		pred_t2m=sum(pred_t2m),
		resp_t2m=sum(resp_t2m),
		pred_m2m=sum(pred_m2m),
		resp_m2m=sum(resp_m2m),
		pred_m2r=sum(pred_m2r),
		resp_m2r=sum(resp_m2r)
		);
	
	max_len=max(lengths);
	max_sum=max(sums);

	calc_lab_pos=function(counts, parv){
		bottoms=c(0,cumsum(counts));
		halfs=c(counts,0)/2;
		midpos=bottoms+halfs;
		names(midpos)=names(counts);
		midpos=midpos[1:length(counts)];
		return(midpos);
	}

	sums_by_type=c(sums["pred_t2m"], sums["pred_m2m"], sums["pred_m2r"]);
	names(sums_by_type)=c("t2m", "m2m", "m2r");
	norm_by_max=sums_by_type/(max(sums_by_type));

	if(any(is.nan(norm_by_max))){
		norm_by_max[1:3]=c(0,0,0);
	}
	norm_by_max=sort(norm_by_max);
	target_min=c(2, 1.5, 1);

	scale=numeric(3);
	names(scale)=norm_by_max;
	scale[1]=max(target_min[1]*norm_by_max[1], norm_by_max[1]);
	scale[2]=max(target_min[2]*norm_by_max[2], norm_by_max[2]);
	scale[3]=max(target_min[3]*norm_by_max[3], norm_by_max[3]);
	names(scale)=names(norm_by_max);


	for(i in 1:3){
		scale[i]=ifelse(scale[i]==0, 1, scale[i]);
	}


	orig_par=par(no.readonly=T);
	par(mfrow=c(1,3));
	par(mar=c(5,5,8,1));
	par(oma=c(0,0,2,0));

	xlabs=c("Predictors", "Responses");

	#----------------------------------------------------------------------
	# Reverse these so labels go from top to bottom, instead of bottom to top

	rev_pred_t2m=rev(pred_t2m);
	rev_resp_t2m=rev(resp_t2m);
	rev_pred_m2m=rev(pred_m2m);
	rev_resp_m2m=rev(resp_m2m);
	rev_pred_m2r=rev(pred_m2r);
	rev_resp_m2r=rev(resp_m2r);

	#----------------------------------------------------------------------

	cat("Treatments and Covariates - Measured Barplots\n");
	m=matrix(0, nrow=max_len, ncol=2);
	m[1:lengths[1], 1]=rev_pred_t2m;
	m[1:lengths[2], 2]=rev_resp_t2m;
	mids=barplot(m, ylim=c(0, scale["t2m"]*max_sum), main="");
	axis(side=1, at=mids, 
		labels=c("\nPredictors:\nTreatments & Covariates", "\nResponses:\nMeasured"),
		tick=F, line=1);

	title(main="Treatments & Covariates\nto\nMeasured", cex.main=2);
	lab_pos_p=calc_lab_pos(rev_pred_t2m, par());
	lab_pos_r=calc_lab_pos(rev_resp_t2m, par());

	barplot(cbind(m[,1], 0), col=grp_colors[names(lab_pos_p)], add=T);
	barplot(cbind(0, m[,2]), col=grp_colors[names(lab_pos_r)], add=T);

	text(rep(mids[1], length(lab_pos_p)), lab_pos_p, labels=names(lab_pos_p));
	text(rep(mids[2], length(lab_pos_r)), lab_pos_r, labels=names(lab_pos_r));

	#----------------------------------------------------------------------

	cat("Measured - Measured Barplots\n");
	m=matrix(0, nrow=max_len, ncol=2);
	m[1:lengths[3], 1]=rev_pred_m2m;
	m[1:lengths[4], 2]=rev_resp_m2m;
	barplot(m, ylim=c(0, scale["m2m"]*max_sum), main="");
	axis(side=1, at=mids,
		labels=c("\nPredictors:\nMeasured,\n(Trt & Cov Excluded)", "Responses:\nMeasured"),
		tick=F, line=2);

	title(main="Measured\nto\nMeasured", cex.main=2);
	lab_pos_p=calc_lab_pos(rev_pred_m2m, par());
	lab_pos_r=calc_lab_pos(rev_resp_m2m, par());

	barplot(cbind(m[,1], 0), col=grp_colors[names(lab_pos_p)], add=T);
	barplot(cbind(0, m[,2]), col=grp_colors[names(lab_pos_r)], add=T);

	text(rep(mids[1], length(lab_pos_p)), lab_pos_p, names(lab_pos_p));
	text(rep(mids[2], length(lab_pos_r)), lab_pos_r, names(lab_pos_r));

	#----------------------------------------------------------------------

	cat("Measured - Response Barplots\n");
	m=matrix(0, nrow=max_len, ncol=2);
	m[1:lengths[5], 1]=rev_pred_m2r;
	m[1:lengths[6], 2]=rev_resp_m2r;
	barplot(m, ylim=c(0, scale["m2r"]*max_sum), main="");
	axis(side=1, at=mids,
		labels=c("\nPredictors:\nTrt & Cov,\nMeasured", "Responses:\nResponse"),
		tick=F, line=2);

	title(main="Measured\nto\nResponse", cex.main=2);
	lab_pos_p=calc_lab_pos(rev_pred_m2r, par());
	lab_pos_r=calc_lab_pos(rev_resp_m2r, par());

	barplot(cbind(m[,1], 0), col=grp_colors[names(lab_pos_p)], add=T);
	barplot(cbind(0, m[,2]), col=grp_colors[names(lab_pos_r)], add=T);

	text(rep(mids[1], length(lab_pos_p)), lab_pos_p, names(lab_pos_p));
	text(rep(mids[2], length(lab_pos_r)), lab_pos_r, names(lab_pos_r));

	#----------------------------------------------------------------------
	mtext(paste("P-value Cutoff: ", cutoff, sep=""), side=3, outer=T);

	par(orig_par);

	marginals=list();
	marginals[["pred_t2m"]]=pred_t2m;
	marginals[["resp_t2m"]]=resp_t2m;
	marginals[["pred_m2m"]]=pred_m2m;
	marginals[["resp_m2m"]]=resp_m2m;
	marginals[["pred_m2r"]]=pred_m2r;
	marginals[["resp_m2r"]]=resp_m2r;

	return(marginals);

}

##############################################################################

plot_tmr_matrices=function(lnk_rec, var_rec, grp_colors, cutoff){

	# Split records into trtcov-msd, msd-msd, and msd-resp matrices
	var_mat=matrix(NA, nrow=0, ncol=3);
	colnames(var_mat)=c("Type", "Group", "Variable");
	for(tmr_type in names(var_rec)){
		cat("TMR Type: ", tmr_type, "\n");
		for(grp in names(var_rec[[tmr_type]])){
			cat("  Group: ", grp, "\n");
			for(var in colnames(var_rec[[tmr_type]][[grp]])){
				cat("    Variable: ", var, "\n");
				var_mat=rbind(var_mat, c(tmr_type, grp, var));
			}
		}
	}
	
	trtcov_mat=var_mat[var_mat[,"Type"]=="Covariates", c("Group", "Variable"), drop=F];
	msd_mat=var_mat[var_mat[,"Type"]=="Measured", c("Group", "Variable"), drop=F];
	rsp_mat=var_mat[var_mat[,"Type"]=="Response", c("Group", "Variable"), drop=F];

	# +/- associations by variable
	tmr_heatmap(lnk_rec[["c2m"]], 
		title="Treatments/Covariates to Measured:  Assoc By Var", 
		subtitle=cutoff,
		trtcov_mat, msd_mat);

	tmr_heatmap(lnk_rec[["m2m"]], 
		title="Measured to Measured:  Assoc By Var", 
		subtitle=cutoff,
		msd_mat, msd_mat);

	tmr_heatmap(lnk_rec[["m2r"]], 
		title="Measured to Response:  Assoc By Var",
		subtitle=cutoff,
		rbind(trtcov_mat, msd_mat), rsp_mat);


	# Counts of associations by group
	trtcov_to_msd_mat=tmr_heatmap_byGroup(lnk_rec[["c2m"]], 
		title="Treatments/Covariates to Measured: Num Assoc By Group", 
		subtitle=cutoff,
		trtcov_mat, msd_mat);

	msd_to_msd_mat=tmr_heatmap_byGroup(lnk_rec[["m2m"]], 
		title="Measured to Measured:  Num Assoc By Group", 
		subtitle=cutoff,
		msd_mat, msd_mat);

	msd_to_rsp_mat=tmr_heatmap_byGroup(lnk_rec[["m2r"]], 
		title="Measured to Response:  Num Assoc By Group", 
		subtitle=cutoff,
		rbind(trtcov_mat, msd_mat), rsp_mat);

	marginals=marginal_stacked_barplots_byGroup(
		trtcov_to_msd_mat, msd_to_msd_mat, msd_to_rsp_mat, 
		grp_colors, cutoff);

	#quit();
}

##############################################################################

plot_combined_fits=function(fits, cutoff){	

	cat("Plotting combined predictors: ", cutoff, "\n");
	orig_par=par(no.readonly=T);
	par(mfrow=c(2,1));

	model_types=names(fits);

	for(mt in model_types){
		
		cat("Working on Model: ", mt, "\n");

		if(length(fits[[mt]])==1 && is.na(fits[[mt]])){
			plot_text(c(
				paste("Cutoff: ", cutoff),
				paste("Model: Using predictors from [", mt, 
					"] to predict response variables.", sep=""),
				"",
				"No significant predictors selected at this cutoff."
			));
			next;
		}
		
		num_resp=length(fits[[mt]][["sumfit_list"]]);
		res_names=names(fits[[mt]][["sumfit_list"]]);
		suff_resid=fits[[mt]][["sufficient_residuals"]];

		cat("Sufficient Residuals: ", suff_resid, "\n");
		cat("Num Responses: ", num_resp, "\n");
		cat("Responses: \n");
		print(res_names);

		if(num_resp && suff_resid){
			for(rix in 1:num_resp){

				cur_resp_name=res_names[rix];
				cat("Response Names: ", cur_resp_name, "\n");
				sumfit_rec=fits[[mt]][["sumfit_list"]][[cur_resp_name]];
				fit_rec=fits[[mt]][["fit_list"]][[cur_resp_name]];

				out_tab=sumfit_rec$coefficients[,c(1,4)];
				prednames=rownames(out_tab);
				out_tab=out_tab[setdiff(prednames, "(Intercept)"),,drop=F];

				# Generate text matrix
				out_formatted=matrix(character(), nrow=nrow(out_tab), ncol=2);
				rownames(out_formatted)=rownames(out_tab);
				out_formatted[,1]=sapply(out_tab[,1], function(x){sprintf("%12.4f", x)});
				out_formatted[,2]=sapply(out_tab[,2], function(x){sprintf("%8.4f", x)});

				# Append significance character to last column
				out_tab_char=cbind(out_formatted, signf_char(out_tab[,2]));
				colnames(out_tab_char)=c("Coefficients", "P-values", "Signif");

				plot_text(c(
					paste("Cutoff: ", cutoff),
					paste("Model: Using predictors from [", mt, 
						"] to predict [", cur_resp_name, "]", sep=""),
					"",
					capture.output(print(out_tab_char, quote=F))
				));	

				resp_range=range(c(fit_rec$y, fit_rec$fitted.values));
				par(mar=c(5,5,5,1));
				plot(fit_rec$y, fit_rec$fitted.values,
					ylim=resp_range, xlim=resp_range,
					xlab="Observed", ylab="Predicted",
					main=cur_resp_name
					);

				obs=fit_rec$y;
				prd=fit_rec$fitted.values;

				obs_pred_fit=lm(prd~obs);
				abline(obs_pred_fit, lty="dotted", lwd=2, col="blue");

			}	
		}
	}

	par(orig_par);

}

##############################################################################

assign_colors=function(name_arr){

	jumble_colors=function(num_colors){

		alloc_colors=rainbow(num_colors, start=0, end=max(1, num_colors-1)/num_colors);

		sqrt_num_col=sqrt(num_colors);
		num_sqr_cells=ceiling(sqrt_num_col)^2;

		holder=rep(NA, num_sqr_cells);
		holder[1:num_colors]=alloc_colors;

		mat=matrix(holder, ncol=ceiling(sqrt_num_col), byrow=T);

		jumbled=as.vector(mat);
		jumbled=jumbled[!is.na(jumbled)];
		return(jumbled);

	}
	
	grp_col=jumble_colors(length(name_arr));
	names(grp_col)=name_arr;

	return(grp_col);
}

###############################################################################

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

# Assign colors to each group
cat("Assigning Colors to: \n");
print(groupings_rec$Groups);
grp_colors=assign_colors(groupings_rec$Groups);
print(grp_colors);

options(width=300);
marginals=list();

#for(cutoffs in c("0.0010", "0.1000")){
#for(cutoffs in c("0.0050")){
for(cutoffs in names(denorm_results)){

	if(cutoffs=="1.0000"){next;}

	plot_TMR_diagram(denorm_results[[cutoffs]],
		paste("P-value Cutoff: ", cutoffs, sep=""),
		"All GROUP Links more significant than cutoff",
		covtrt_to_group_map,
		grp_links=1,
		covariates_list, measured_list, response_list, overfits);

	plot_TMR_diagram(denorm_results[[cutoffs]],
		paste("P-value Cutoff: ", cutoffs, sep=""),
		"All VARIABLE Links more significant than cutoff)",
		covtrt_to_group_map,
		grp_links=0,
		covariates_list, measured_list, response_list, overfits);

	if(Verbose){
		plot_text(c(
			paste("P-value Cutoff: ", cutoffs, sep=""),
			paste("(All links above cutoff, ", nrow(denorm_results[[cutoffs]]), 
				" links.)", sep=""),
			capture.output(print(denorm_results[[cutoffs]], quotes=""))
		), max_lines_pp=70);
	}

	# Remove weaker of bi-directional links

	unidir_links=remove_weaker_bidirectional_links(denorm_results[[cutoffs]], log10_diff_thres=1);
	plot_TMR_diagram(unidir_links,
		paste("P-value Cutoff: ", cutoffs, sep=""),
		"Uni-directional VARIABLE Links with stronger associations",
		covtrt_to_group_map,
		grp_links=0,
		covariates_list, measured_list, response_list, overfits);

	if(Verbose){
		plot_text(c(
			paste("P-value Cutoff: ", cutoffs, sep=""),
			paste("(Excluding weaker of bi-directional links, ", nrow(unidir_links), 
				" links.)", sep=""),
			capture.output(print(unidir_links, quotes=""))
		), max_lines_pp=70);
	}


	# Remove links from covariates/treatments

	no_trtcov_unidir_links=remove_covariate_links_from_msd_pred(
		unidir_links, covariate_variable_names);
	plot_TMR_diagram(no_trtcov_unidir_links,
		paste("P-value Cutoff: ", cutoffs, sep=""),
		"Uni-directional VARIABLE Links and w/o Treatments/Covariates",
		covtrt_to_group_map,
		grp_links=0,
		covariates_list, measured_list, response_list, overfits);

	plot_text(c(
		paste("P-value Cutoff: ", cutoffs, sep=""),
		paste("(Excluding weaker of bi-directional links, and Treatments/Covariates)", sep=""),
		capture.output(print(no_trtcov_unidir_links, quotes=""))
	), max_lines_pp=70);
		

	# Generate Heatmaps & Stacked Barplots
	link_rec=extract_matrices_from_links(no_trtcov_unidir_links);
	marginals[[cutoffs]]=plot_tmr_matrices(link_rec, variables_rec, grp_colors, cutoffs);


	# Show combined selected
	plot_combined_fits(selected_fits[[cutoffs]], cutoffs);
}

##############################################################################

plot_marginals_over_signif=function(marginals_list, grp_col){

	signf_arr=names(marginals_list);

	signf_numeric=as.numeric(signf_arr);
	# Reorder, so start with least signficant/loosest cutoff
	if(signf_numeric[1]<signf_numeric[2]){
		signf_arr=rev(signf_arr);
	}


	list_by_type=list();

	padded_bind=function(m, a){
		# This will add a line to the matrix with 0 padding
		if(is.null(m)){
			mat=matrix(a, nrow=1);
			colnames(mat)=names(a);
			return(mat)
		}else{
			mat=m;
		}
		
		mat=rbind(mat, 0);
		last_row=nrow(mat);
		mat[last_row, names(a)]=a;
		return(mat);
	}
	
	#----------------------------------------------------------------------

	for(signf in signf_arr){
		cat("Working on: ", signf, "\n");
		marginals=marginals_list[[signf]]

		for(type in names(marginals)){
			list_by_type[[type]]=padded_bind(list_by_type[[type]], marginals[[type]]);
		}
	}

	# Label signif and get all groups	
	uniq_grps_seen=c();
	for(type in names(list_by_type)){
		rownames(list_by_type[[type]])=signf_arr;
		uniq_grps_seen=unique(c(uniq_grps_seen, colnames(list_by_type[[type]])));
	}

	cat("Groups seen:\n");
	print(uniq_grps_seen);

	#----------------------------------------------------------------------

	plot_matrices=function(mat, colmap, title, m, al, ylab){

		cat("Working on: ", title, "\n");

		par(mar=m);
	
		num_xpts=nrow(mat);
		xlab=rownames(mat);
		max_counts=max(mat);
		num_var=ncol(mat);
		varnames=colnames(mat);

		if(sum(mat)==0){
			plot(NA, type="n", ylim=c(-1, 1), xlim=c(-1,1),
				main=title, xaxt="n", ylab=ylab, xlab="p-value");
			text(0,0, "No Associations Found.\n");
		
			return();
		}
	
		plot(NA, type="n", ylim=c(0, max_counts), xlim=c(1,num_xpts),
			main=title, xaxt="n", ylab=ylab, xlab="p-value");

		if(al){
			axis(side=1, at=1:num_xpts, labels=xlab);
		}else{
			axis(side=1, at=1:num_xpts, labels=rep("", num_xpts));
		}

		for(i in 1:num_var){
			cvn=varnames[i];
			col=colmap[cvn];
			points(1:num_xpts, mat[,i], type="l", col=col, lwd=1.25);
			points(1:num_xpts, mat[,i], type="l", col="black", lwd=.25);
		}
		
		pr=par()$usr;
		xr=pr[2]-pr[1];
		yr=pr[4]-pr[3];
		legend(x=pr[1]+xr*3/4, y=pr[3]+yr*7/8, fill=colmap[varnames], legend=varnames); 
			
	}

	#----------------------------------------------------------------------
	# Generate line plots with x=p-val / y=counts
	# Pair the Predictors with Responses

	par(mfrow=c(3,2));

	plot_matrices(list_by_type[["pred_t2m"]], grp_col, title="Predictors", 
		ylab="Covariates & Treatments", m=c(2,4,4,1), al=T);
	plot_matrices(list_by_type[["resp_t2m"]], grp_col, title="Responses", 
		ylab="Measured", m=c(2,4,4,1), al=T);
	plot_matrices(list_by_type[["pred_m2m"]], grp_col, title="", 
		ylab="Measured", m=c(2,4,1,1), al=T);
	plot_matrices(list_by_type[["resp_m2m"]], grp_col, title="",
		ylab="Measured", m=c(2,4,1,1), al=T);
	plot_matrices(list_by_type[["pred_m2r"]], grp_col, title="", 
		ylab="Measured", m=c(4,4,1,1), al=T);
	plot_matrices(list_by_type[["resp_m2r"]], grp_col, title="", 
		ylab="Responses", m=c(4,4,1,1), al=T);
}

plot_marginals_over_signif(marginals, grp_colors);

#-----------------------------------------------------------------------------

plot_signif_combined_predictors=function(fits_across_cutoffs){

	# The fits are grouped by model type: Cov_to_Msd, Msd_to_Msd, and Msd_to_Rsp

	orig_par=par(no.readonly=T);

	models=c();
	response_variables=c();
	max_predictors=0;

	signif_pred_rec=list();
	suff_res_rec=list();

	pval_cutoffs=names(fits_across_cutoffs);
	signf_pred_by_resp=list();

	# Extract the number of predictors, and those significant
	for(pvc in pval_cutoffs){
		cat("P-value: ", pvc, "\n");

		fits=fits_across_cutoffs[[pvc]];

		signif_pred_rec[[pvc]]=list();
		suff_res_rec[[pvc]]


		avail_fits=names(fits);
		for(mt in avail_fits){

			cat("Models: ", mt, "\n");
			
			sumfit_list=fits[[mt]][["sumfit_list"]];

			# Sufficient residuals apply to all response variables
			#   at the same cutoff and model type
			suf_res_list=fits[[mt]][["sufficient_residuals"]];
			predictors=fits[[mt]][["predictors"]];
			num_pred=length(predictors);
			#cat("SumFit List:\n");
			#print(sumfit_list);

			signif_pred_rec[[pvc]][[mt]]=list();
			models=c(models, mt);

			if(length(sumfit_list)){
				for(var in names(sumfit_list)){
		
					cat("Variables: ", var, "\n");

					sumfit=sumfit_list[[var]];

					coef_mat=sumfit[["coefficients"]];
					pred_names=setdiff(rownames(coef_mat), "(Intercept)");

					coef_mat=coef_mat[pred_names,,drop=F];
					#print(coef_mat);
					#num_pred=nrow(coef_mat)

					signf_ix=coef_mat[,"Pr(>|.|)"]<0.1;
					signif_pred=rownames(coef_mat)[signf_ix];
					num_signif_pred=length(signif_pred);

					cat("Num Significant: ", num_signif_pred, "\n");
					cat("\n\n");
				}
			}else{
				num_signif_pred=0;
				signif_pred=c();
			}

			arr=c(num_pred, num_signif_pred, 
				ifelse(is.null(suf_res_list), T, suf_res_list));

			names(arr)=c("Num Pred", "Num Signf Pred", "Suff Resid");
			signif_pred_rec[[pvc]][[mt]][[var]]=arr;

			signf_pred_by_resp[[var]]=
				unique(c(signf_pred_by_resp[[var]], signif_pred));

			response_variables=c(response_variables, var);
			max_predictors=max(c(max_predictors, num_pred));

		}
		cat("\n");
	}

	#cat("Extracted:\n");
	print(signif_pred_rec);

	#----------------------------------------------------------------------
	# Generate Barplots

	par(mfrow=c(3,1));
	par(mar=c(8,5,5,1));

	models=unique(models);
	num_models=3;
	#num_models=length(models);
	response_variables=unique(response_variables);

	num_pval_cutoffs=length(pval_cutoffs);

	cat("Max predictors: ", max_predictors, "\n");

	num_bars=num_pval_cutoffs*(num_models+1);

	for(resp_var in response_variables){
		cat("Working on: ", resp_var, "\n");
		
		mids=barplot(height=rep(NA, num_bars),
			border=NA,
			main=paste("Response Variable: ", resp_var, sep=""),
			xaxt="n",
			ylab="Number of Predictors Selected",
			ylim=c(0, max_predictors+1));

		axis(side=1, at=mids, labels=rep(c("TM","MM","MR", " "), num_pval_cutoffs), tick=F,
			 cex.axis=.95, line=0);

		pval_labels=rep("", num_bars);
		pval_labels[(0:(num_pval_cutoffs-1))*4+1]=pval_cutoffs;
		axis(side=1, at=mids, labels=pval_labels, tick=F, line=2, font.axis=2);

		x_buf=rep(0, (num_models+1)*num_pval_cutoffs);
		for(pvc_ix in 1:num_pval_cutoffs){

			pvc=pval_cutoffs[pvc_ix];

			for(m_ix in 1:num_models){

				mt=models[m_ix];

				counts=signif_pred_rec[[pvc]][[mt]][[resp_var]];
				if(is.null(counts) || length(counts)==0){
					counts=c(0,0,1);
					names(counts)=c("Num Pred", "Num Signf Pred", "Suff Resid");
				}

				column_ix=(pvc_ix-1)*(num_models+1)+ (m_ix-1) + 1;

				filled_buf=x_buf;
				filled_buf[column_ix]=counts["Num Pred"];

				if(counts["Suff Resid"]==0){
					barcol="grey";
				}else{
					barcol="green";
				}

				barplot(height=filled_buf, xaxt="n", yaxt="n", 
					border=NA,
					add=T, col=barcol, tick=F);

				filled_buf[column_ix]=counts["Num Signf Pred"];
				barplot(height=filled_buf, xaxt="n", yaxt="n", 
					border=NA,
					add=T, col="blue", tick=F);

				if(counts["Suff Resid"]==0){
					text(mids[column_ix], 0, labels="O", pos=3, 
						font=2, col="red");
				}
			}

			
		}

		#--------------------------------------------------------------

		# Label significant predictors
		mid_val=(max(mids)-min(mids))/2;
		signf_pred_banner=paste(sort(signf_pred_by_resp[[resp_var]]), collapse=", ");	
		axis(side=1, at=mid_val, labels="Significant Predictors:", tick=F, line=4, 
			font.axis=2, cex.axis=1, cex.axis=1.1);
		axis(side=1, at=mid_val, labels=signf_pred_banner, tick=F, line=5, 
			font.axis=3, cex.axis=1, cex.axis=1.1);

		# legend
		gr_sp=par()$usr;
		range=gr_sp[2]-gr_sp[1];
		
		legend(gr_sp[1]+range*2/3, gr_sp[4], 
			fill=c("green", "blue", "grey"),
			legend=c("Attempted", "Significant", "Overfit"),
			horiz=T
		);

		# Draw dividers between significance cutoffs
		abline(v=mids[((1:num_pval_cutoffs)-1)*4+4], lty="dotted");

		# Hack to remove 0's lines from between models
		barplot(height=x_buf, xaxt="n", yaxt="n", 
			border="white", lwd=1.25,
			add=T, col="white", tick=F);

	}

	par(orig_par);

}

plot_signif_combined_predictors(selected_fits);

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

file.rename(
	paste(OutputFnameRoot, ".tmr.started", sep=""),
	paste(OutputFnameRoot, ".tmr.done", sep="")
);

##############################################################################

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
