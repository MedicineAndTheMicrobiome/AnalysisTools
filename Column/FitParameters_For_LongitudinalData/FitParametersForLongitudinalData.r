#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);
options(width=120);

params=c(
	"input", "i", 1, "character",
	"time_offset_cn", "t", 1, "character",
	"subject_id_cn", "s", 1, "character",
	"target_values_fn", "v", 1, "character",
	"output", "o", 1, "character",

	"find_line", "l", 2, "character",
	"find_limits", "m", 2, "character",
	"find_descriptive", "d", 2, "character",
	"find_ranges", "r", 2, "character",
	"find_timing", "n", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input tab-separated column data file>\n",
	"	--time_offset_cn=<column name of time offset>\n",
	"	--subject_id_cn=<column name of subject IDs>\n",
	"	--target_values_fn=<filename of list of target column names>\n",
	"\n",
	"	Specify at least one of these find/fit parameters:\n",
	"	[-l (fit line: intercept/slope)]\n",
	"	[-m (find limits: min/max/range)]\n",
	"	[-d (find descriptive stats: mean/stdev/median)]\n",
	"	[-r (find ranges: first/last/N)]\n",
	"	[-n (find timing: time_span/freq)]\n",
	"\n",
	"	-o <output filename>\n",
	"\n");

if(
	!length(opt$input) || 
	!length(opt$time_offset_cn) || 
	!length(opt$subject_id_cn) ||
	!length(opt$target_values_fn) ||
	!length(opt$output)

){
	cat(usage);
	q(status=-1);
}

InputFname=opt$input;
TimeOffsetColName=opt$time_offset_cn;
SubjectIDColName=opt$subject_id_cn;
TargetValuesFname=opt$target_values_fn;
OutputFname=opt$output;

Opt_FindLine=F;
Opt_FindLimits=F;
Opt_FindDescriptive=F;
Opt_FindRanges=F;
Opt_FindTiming=F;

if(length(opt$find_line)){
	Opt_FindLine=T;
}
if(length(opt$find_limits)){
	Opt_FindLimits=T;
}
if(length(opt$find_descriptive)){
	Opt_FindDescriptive=T;
}
if(length(opt$find_ranges)){
	Opt_FindRanges=T;
}
if(length(opt$find_timing)){
	Opt_FindTiming=T;
}

cat("\n");
cat("Input Filename:", InputFname , "\n");
cat("Time Offset Colname:", TimeOffsetColName, "\n");
cat("SubjectID Colname:", SubjectIDColName, "\n");
cat("Target Values List Filename:", TargetValuesFname, "\n");
cat("Output Filename:", OutputFname, "\n");
cat("\n");
cat("Find Line: ", Opt_FindLine, "\n");
cat("Find Limits: ", Opt_FindLimits, "\n");
cat("Find Descriptive: ", Opt_FindDescriptive, "\n");
cat("Find Ranges: ", Opt_FindRanges, "\n");
cat("Find Timing: ", Opt_FindTiming, "\n");

if(!(Opt_FindLine || Opt_FindLimits || Opt_FindDescriptive || Opt_FindRanges || Opt_findTiming)){
	cat("At least one of the find/fit parameters must be specified.\n");
}

##############################################################################

load_factors=function(fname){
	factors=data.frame(read.table(fname,  header=TRUE, check.names=FALSE, 
		as.is=T, comment.char="", quote="", sep="\t"));
	dimen=dim(factors);
	cat("Rows Loaded: ", dimen[1], "\n");
	cat("Cols Loaded: ", dimen[2], "\n");
	return(factors);
}

write_factors=function(fname, table){
	dimen=dim(table);
	cat("Rows Exporting: ", dimen[1], "\n");
	cat("Cols Exporting: ", dimen[2], "\n");
	write.table(table, fname, quote=F, row.names=F, sep="\t");
}

load_list=function(fname){
	arr=readLines(fname);	
	return(arr);
}

##############################################################################

# Load factors
all_factors=load_factors(InputFname);
num_samples=nrow(all_factors);

cat("Excerpt from: ", InputFname, "\n");
print(head(all_factors));

# Load targets
target_variable_list=load_list(TargetValuesFname);
num_target_variables=length(target_variable_list);
cat("Target Variables List:\n");
print(target_variable_list);


##############################################################################

# Extract out colums we need:
subject_ids_arr=all_factors[,SubjectIDColName];
time_offset_arr=all_factors[,TimeOffsetColName];
targeted_values_mat=all_factors[,target_variable_list, drop=F];

cat("\n");
cat("Subject IDs Excerpt:\n");
print(head(subject_ids_arr));
cat("\n");
cat("Time Offset Excerpt:\n");
print(head(time_offset_arr));
cat("\n");
cat("Targeted Values Excerpt:\n");
print(head(targeted_values_mat));
cat("\n");

##############################################################################

unique_subject_ids=sort(unique(subject_ids_arr));
num_unique_subject_ids=length(unique_subject_ids);
cat("Num Unique Subject IDs: ", num_unique_subject_ids, "\n");

##############################################################################

extensions=c();
if(Opt_FindLine){
	extensions=c(extensions, c("intercept", "slope", "sd_res"));
	calc_line=function(data){

		x=data[,1];
		y=data[,2];

		lmfit=lm(y~x);
		intercept=lmfit$coefficients["(Intercept)"];		
		slope=lmfit$coefficients["x"];		
		sd_res=sd(lmfit$residuals);

		return(c(intercept,slope, sd_res));
	}
}
if(Opt_FindLimits){
	extensions=c(extensions, c("min", "max", "range"));
	calc_limits=function(data){

		x=data[,1];
		y=data[,2];

		min=min(y);
		max=max(y);
		range=max-min;

		return(c(min, max, range));
	}
}
if(Opt_FindDescriptive){
	extensions=c(extensions, c("mean", "stdev", "median"));
	calc_descriptive=function(data){

		x=data[,1];
		y=data[,2];

		mean=mean(y);
		stdev=sd(y);
		median=median(y);

		return(c(mean, stdev, median));
	}
}
if(Opt_FindRanges){
	extensions=c(extensions, c("first", "last", "N"));
	calc_ranges=function(data){

		x=data[,1];
		y=data[,2];

		N=length(y);
		first=y[1];
		last=y[N];

		return(c(first, last, N));
	}
}
if(Opt_FindTiming){
	extensions=c(extensions, c("time_span", "freq"));
	calc_timing=function(data){

		x=data[,1];

		first_time=min(x);
		last_time=max(x);
		N=length(x);
		time_span=last_time-first_time;	
		freq=N/time_span;

		return(c(time_span, freq));
	}
}

##############################################################################

cat("Variable stats/extensions:\n");
print(extensions);
num_stats=length(extensions);
num_new_variables=num_stats*num_target_variables;
new_var_names=paste(target_variable_list, "_", extensions, sep="");
cat("\n");
cat("Num New Variables: ", num_new_variables, "\n");
print(new_var_names);

acc_matrix=matrix(NA, nrow=num_unique_subject_ids, ncol=num_stats)
colnames(acc_matrix)=new_var_names;
rownames(acc_matrix)=unique_subject_ids;

ranges_mat=matrix(NA, nrow=num_target_variables, ncol=2);
rownames(ranges_mat)=target_variable_list;
colnames(ranges_mat)=c("min", "max");

#ranges_mat[,"min"]=apply(targeted_values_mat, 2,  min);
#ranges_mat[,"max"]=apply(targeted_values_mat, 2,  max);

stdevs=apply(targeted_values_mat, 2,  sd);
means=apply(targeted_values_mat, 2,  mean);
ranges_mat[,"max"]=means+stdevs*3;
ranges_mat[,"min"]=means-stdevs*3;

print(ranges_mat);

multiplot_dim=ceiling(sqrt(num_target_variables));
cat("Multiplot Dimensions: ", multiplot_dim, "\n");

pdf(file=paste(OutputFname, ".longit.pdf", sep=""), height=8.5, width=11);


for(cur_subj in unique_subject_ids){

	subj_ix=(cur_subj==subject_ids_arr);
	subj_mat=all_factors[subj_ix, c(SubjectIDColName, TimeOffsetColName, target_variable_list)];

	sort_ix=order(subj_mat[,TimeOffsetColName], decreasing=F);
	subj_mat_sorted=subj_mat[sort_ix,,drop=F];

	#print(subj_mat_sorted);

	par(mfrow=c(multiplot_dim, multiplot_dim));
	for(targ_var_ix in target_variable_list){

		# Plot Variables
		plot(subj_mat_sorted[,TimeOffsetColName], subj_mat_sorted[,targ_var_ix],
			main=paste(cur_subj, " : ", targ_var_ix, sep=""),
			xlab=TimeOffsetColName,
			ylab=targ_var_ix,
			ylim=c(ranges_mat[targ_var_ix, "min"], ranges_mat[targ_var_ix, "max"])
		);

		# Calculate parameters
		if(Opt_FindLine){
			acc_matrix[cur_subj, paste(targ_var_ix, "_", c("intercept", "slope", "sd_res"), sep="")]=
				calc_line(subj_mat_sorted[,c(TimeOffsetColName, targ_var_ix),drop=F]);
		}

		if(Opt_FindLimits){
			acc_matrix[cur_subj, paste(targ_var_ix, "_", c("min", "max", "range"), sep="")]=
				calc_limits(subj_mat_sorted[,c(TimeOffsetColName, targ_var_ix),drop=F]);
		}

		if(Opt_FindDescriptive){
			acc_matrix[cur_subj, paste(targ_var_ix, "_", c("mean", "stdev", "median"), sep="")]=
				calc_descriptive(subj_mat_sorted[,c(TimeOffsetColName, targ_var_ix),drop=F]);
		}

		if(Opt_FindRanges){
			acc_matrix[cur_subj, paste(targ_var_ix, "_", c("first", "last", "N"), sep="")]=
				calc_ranges(subj_mat_sorted[,c(TimeOffsetColName, targ_var_ix),drop=F]);

		}
		if(Opt_FindTiming){
			acc_matrix[cur_subj, paste(targ_var_ix, "_", c("time_span", "freq"), sep="")]=
				calc_timing(subj_mat_sorted[,c(TimeOffsetColName, targ_var_ix),drop=F]);

		}

	}	

}

rounded_acc_mat=round(acc_matrix, 4);

out_hdr=c(SubjectIDColName, colnames(rounded_acc_mat));
outmat=cbind(unique_subject_ids, rounded_acc_mat);

colnames(outmat)=out_hdr;

write_factors(paste(OutputFname, ".longit.tsv", sep=""), outmat);

##############################################################################

cat("\nDone.\n");

print(warnings());
q(status=0);
