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
	"response_grp_list_fn", "r", 1, "character"
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
	"	-r <response group list file>\n",
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
	!length(opt$measured_grp_list_fn) || 
	!length(opt$response_grp_list_fn)
){
	cat(usage);
	q(status=-1);
}

FactorsFname=opt$factors;
VariableGroupingsFname=opt$groupings;
OutputFnameRoot=opt$outputroot;
CovTrtGrpFname=opt$covtrt_grp_list_fn;
MeasuredGrpFname=opt$measured_grp_list_fn;
ResponseGrpFname=opt$response_grp_list_fn;


param_text=capture.output({
	cat("\n");
	cat("Factor/Metadata Filename:            ", FactorsFname, "\n");
	cat("Variable Groupings Filename:         ", VariableGroupingsFname, "\n");
	cat("Output Filename Root:                ", OutputFnameRoot, "\n");
	cat("Covariate/Treatment Groups Filename: ", CovTrtGrpFname, "\n");
	cat("Measured Groups Filename:            ", MeasuredGrpFname, "\n");
	cat("Response Groups Filename:            ", ResponseGrpFname, "\n");
	cat("\n");
});

print(param_text, quote=F);

###############################################################################

load_factors=function(fname){
	factors=as.data.frame(read.delim(fname,  header=TRUE, row.names=1, check.names=FALSE, sep="\t"));
	return(factors);
}

#-----------------------------------------------------------------------------#

load_list=function(fname){
	cat("Loading: ", fname, "\n");
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

generic_distance=function(data){
	
}

##############################################################################
# Main Program Starts Here!
##############################################################################

pdf(paste(OutputFnameRoot, ".tmr.pdf", sep=""), height=8.5, width=11);
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

# Load Groupings
groupings_rec=load_groupings(VariableGroupingsFname);
#print(groupings_rec);

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

print(varlists_text, quote=F);
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

variables_rec=list();
variables_rec[["Covariates"]]=split_factors_to_groups(loaded_factors, groupings_rec, covariates_list);
variables_rec[["Measured"]]=split_factors_to_groups(loaded_factors, groupings_rec, measured_list);
variables_rec[["Response"]]=split_factors_to_groups(loaded_factors, groupings_rec, response_list);

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


print(model_variables_summary);
plot_text(c(
	"Model Variables:",
	"",	
	model_variables_summary
), max_lines_pp=80);


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

		pval_mat[covariate_variable_names, varname]=
			sum_fit[[var_ix]][["coefficients"]][covariate_variable_names,"Pr(>|t|)"];

		coef_mat[covariate_variable_names, varname]=
			sum_fit[[var_ix]][["coefficients"]][covariate_variable_names,"Estimate"];

	}

	# Store in record
	model_results[["Cov_to_Msd"]][[msd_ix]][["pval"]]=pval_mat;
	model_results[["Cov_to_Msd"]][[msd_ix]][["coef"]]=coef_mat;

}

#print(model_results[["Cov_to_Msd"]]);

##############################################################################
# Fit (Measured + Covariates) to Predict Measured

for(pred_msd_ix in measured_list){
	for(resp_msd_ix in measured_list){

		if(pred_msd_ix==resp_msd_ix){
			next;
		}

		cat("Fitting: Measured to Measured: (Pred)", pred_msd_ix, " (Resp)", resp_msd_ix, "\n");
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
		#print(sum_fit);

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

			pval_mat[cov_and_pred_names, varname]=
				sum_fit[[var_ix]][["coefficients"]][cov_and_pred_names,"Pr(>|t|)"];

			coef_mat[cov_and_pred_names, varname]=
				sum_fit[[var_ix]][["coefficients"]][cov_and_pred_names,"Estimate"];

		}

		# Store in record
		model_results[["Msd_to_Msd"]][[analysis_string]][["pval"]]=pval_mat;
		model_results[["Msd_to_Msd"]][[analysis_string]][["coef"]]=coef_mat;

	}

}

#print(model_results[["Msd_to_Msd"]]);

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

			pval_mat[cov_and_pred_names, varname]=
				sum_fit[[var_ix]][["coefficients"]][cov_and_pred_names,"Pr(>|t|)"];

			coef_mat[cov_and_pred_names, varname]=
				sum_fit[[var_ix]][["coefficients"]][cov_and_pred_names,"Estimate"];

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
					if(pval_mat[pix, rix]<=pval_cutoff){
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
for(pvco in rev(c(0.1000, 0.050, 0.010, 0.005, 0.001, 0.0005, 0.0001))){
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

			type=as.character(rr[i, "model_type"]);
			name=as.character(rr[i, "model_name"]);
			pred_var=as.character(rr[i, "predictor"]);
			resp_var=as.character(rr[i, "response"]);
			direction=ifelse(rr[i, "coef"]>0, "+", "-");
	
			grplink=strsplit(name, "->")[[1]];
			pred_grp=grplink[1];
			resp_grp=grplink[2];
		
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
	for(i in 1:num_cvtrt_grps){
		ichar=sprintf("%i", i);
		var_labels=extr_links[["covtrt"]][["pred_grp_members"]][[ichar]];
		loc=draw_squares_centered(covtrt_xpos, covtrt_ypos[i],
			height=covtrt_height, width=all_widths,
			grp_name=cvtrt_grps[i], 
			variables=var_labels, text_align=ifelse(grp_links, "center", "right"));

		extr_links[["covtrt"]][["pred_grp_members_loc"]][[ichar]]=loc;
	}


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

	for(i in 1:num_rsp_grps){
		ichar=sprintf("%i", i);
		var_labels=extr_links[["resp"]][["resp_grp_members"]][[ichar]];
		loc=draw_squares_centered(rsp_xpos, rsp_ypos[i],
			height=rsp_height, width=all_widths,
			grp_name=rsp_grps[i], 
			variables=var_labels, text_align=ifelse(grp_links, "center", "left"));
		extr_links[["resp"]][["resp_grp_members_loc"]][[ichar]]=loc;
	}

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
	rownames(out)=1:nrow(out);
	return(out);	

}

#for(cutoffs in c("0.0010", "0.1000")){
for(cutoffs in names(denorm_results)){

	plot_TMR_diagram(denorm_results[[cutoffs]],
		paste("P-value Cutoff: ", cutoffs, sep=""),
		"Group Links",
		covtrt_to_group_map,
		grp_links=1,
		covariates_list, measured_list, response_list);

	plot_TMR_diagram(denorm_results[[cutoffs]],
		paste("P-value Cutoff: ", cutoffs, sep=""),
		"(All links above cutoff)",
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
		"(Excluding weaker of bi-directional links)",
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

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
