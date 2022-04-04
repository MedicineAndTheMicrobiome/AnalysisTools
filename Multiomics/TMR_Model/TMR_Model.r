#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

options(useFancyQuotes=F);
options(digits=5);
options(width=120);

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

print(groupings_rec);

# Load variable types
covariates_list=load_list(CovTrtGrpFname);
measured_list=load_list(MeasuredGrpFname);
response_list=load_list(ResponseGrpFname);

varlists_text=capture.output({
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

##############################################################################

split_factors_to_groups=function(data_matrix, group_list, grp_arr){

	data_grp=list();
	for(grp_id in grp_arr){
		data_grp[[grp_id]]=data_matrix[,group_list[["GrpVarMap"]][[grp_id]], drop=F];
	}
	return(data_grp);

}

variables_rec=list();
variables_rec[["Covariates"]]=split_factors_to_groups(loaded_factors, groupings_rec, covariates_list);
variables_rec[["Measured"]]=split_factors_to_groups(loaded_factors, groupings_rec, measured_list);
variables_rec[["Response"]]=split_factors_to_groups(loaded_factors, groupings_rec, response_list);

print(variables_rec);

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

standardized_variables_rec=apply_fun_to_var_rec(variables_rec, standardize_matrix);
print(standardized_variables_rec);

mds_rec=apply_fun_to_var_rec(standardized_variables_rec, calc_mds_matrix);
print(mds_rec);

##############################################################################

covariates_data=matrix(NA, nrow=nrow(loaded_factors), ncol=0);
for(grp in names(variables_rec[["Covariates"]])){
print(variables_rec[["Covariates"]][[grp]]);
	covariates_data=cbind(covariates_data, variables_rec[["Covariates"]][[grp]]);
}

covariate_variable_names=colnames(covariates_data);

colors=matrix(NA, nrow=nrow(covariates_data), ncol=ncol(covariates_data));
colnames(colors)=colnames(covariates_data);
rownames(colors)=rownames(covariates_data);

quantize=function(x, steps){
	min=min(x);
	max=max(x);
	range=max-min;
	norm=(x-min)/range;
	return(floor(norm*(steps-1))+1);
}

centers=function(x, steps){
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

for(type in c("Measured", "Response")){
	for(grp in names(mds_rec[[type]])){
		mds_coord=mds_rec[[type]][[grp]];		
	
		for(var_ix in covariate_variable_names){

			plot(mds_coord[,1], mds_coord[,2], col=colors[, var_ix],
				xlab="Dim 1", ylab="Dim 2", type="n",
				main=paste("Group: ", grp, sep=""),
				);
			title(main=paste("(Colored by: ", var_ix, ")", sep=""), line=.66, cex.main=.75);

			plot_ranges=par()$usr; # left, right, bottom, top
			points(mds_coord[,1], mds_coord[,2], col=colors[,var_ix], cex=1.5, lwd=3);
			points(mds_coord[,1], mds_coord[,2], col="black", cex=1.5, lwd=.25);
			text(mds_coord[,1], mds_coord[,2], rownames(mds_coord), cex=.5, pos=3);

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

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
