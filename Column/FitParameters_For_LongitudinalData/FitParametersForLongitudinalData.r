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
	"group_cn", "g", 2, "character",

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
	"	[--group_cn=<column name of subject grouping>]\n",
	"\n",
	"	Specify at least one of these find/fit parameters:\n",
	"	[-l (fit line: intercept/slope)]\n",
	"	[-m (find limits: min/max/range)]\n",
	"	[-d (find descriptive stats: mean/stdev/median)]\n",
	"	[-r (find ranges: first/last/N)]\n",
	"	[-n (find timing: timespan/freq)]\n",
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

GroupColName="";
if(length(opt$group_cn)){
	GroupColName=opt$group_cn;
}

cat("\n");
cat("Input Filename:", InputFname , "\n");
cat("Time Offset Colname:", TimeOffsetColName, "\n");
cat("SubjectID Colname:", SubjectIDColName, "\n");
cat("Target Values List Filename:", TargetValuesFname, "\n");
cat("Output Filename:", OutputFname, "\n");
cat("Group Colname: ", GroupColName, "\n");
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

		return(c(intercept, slope, sd_res));
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
	extensions=c(extensions, c("timespan", "freq"));
	calc_timing=function(data){

		x=data[,1];

		first_time=min(x);
		last_time=max(x);
		N=length(x);
		timespan=last_time-first_time;	
		freq=N/timespan;

		return(c(timespan, freq));
	}
}

##############################################################################

cat("Variable stats/extensions:\n");
print(extensions);
num_stats=length(extensions);
num_new_variables=num_stats*num_target_variables;

new_var_names=c();
target_to_newvar_map=list();
for(var_ix in target_variable_list){
	var_specific=character();
	for(ext_ix in extensions){
		name_w_ext=paste(var_ix, "_", ext_ix, sep="");
		new_var_names=c(new_var_names, name_w_ext);
		var_specific=c(var_specific, name_w_ext);
	}
	target_to_newvar_map[[var_ix]]=var_specific;
}

cat("\n");
cat("Num New Variables: ", num_new_variables, "\n");
print(new_var_names);
print(target_to_newvar_map);

acc_matrix=matrix(NA, nrow=num_unique_subject_ids, ncol=num_new_variables);
colnames(acc_matrix)=new_var_names;
rownames(acc_matrix)=unique_subject_ids;

ranges_mat=matrix(NA, nrow=num_target_variables, ncol=2);
rownames(ranges_mat)=target_variable_list;
colnames(ranges_mat)=c("min", "max");

#ranges_mat[,"min"]=apply(targeted_values_mat, 2,  min);
#ranges_mat[,"max"]=apply(targeted_values_mat, 2,  max);

stdevs=apply(targeted_values_mat, 2,  function(x){sd(x, na.rm=T)});
means=apply(targeted_values_mat, 2,  function(x){mean(x, na.rm=T)});
ranges_mat[,"max"]=means+stdevs*3;
ranges_mat[,"min"]=means-stdevs*3;

cat("\n");
cat("Plot Ranges (+/- 3*sd):\n");
print(ranges_mat);

multiplot_dim_r=ceiling(sqrt(num_target_variables));
multiplot_dim_c=multiplot_dim_r;

if(multiplot_dim_c*(multiplot_dim_r-1)>=num_target_variables){
	multiplot_dim_r=multiplot_dim_r-1;	
}
cat("Multiplot Dimensions: ", multiplot_dim_r, " x ", multiplot_dim_c, "\n");

plot_var=function(times, values, var_name, subject_id, val_lim){

	plot(times, values, type="n", 
		main=var_name,
		ylim=val_lim,
		xlab="", ylab="" 
	);

	fit=lm(values~times);
	slope=fit$coefficients["times"];
	intercept=fit$coefficients["(Intercept)"];

	lowest_ix=min(which(values==min(values, na.rm=T)));
	highest_ix=min(which(values==max(values, na.rm=T)));

	if(slope>=0){
		abline(a=intercept, b=slope, col="brown", lwd=.5, lty="dashed");
	}else{
		abline(a=intercept, b=slope, col="darkolivegreen", lwd=.5, lty="dashed");
	}

	points(times, values, type="b", col="grey", lwd=2);
	points(times, values, type="p", col="black");

	# Label min/max
	points(times[lowest_ix], values[lowest_ix], col="blue", pch=24, bg="blue", cex=1.5);
	points(times[highest_ix], values[highest_ix], col="red", pch=25, bg="red", cex=1.5);

}

var_by_subj=list();

pdf(file=paste(OutputFname, ".longit.pdf", sep=""), height=8.5, width=11);

par(oma=c(0,0,2,0));
par(mar=c(2,2,3,1));

time_ranges_min=min(all_factors[,TimeOffsetColName]);
time_ranges_max=max(all_factors[,TimeOffsetColName]);


for(cur_subj in unique_subject_ids){

	subj_ix=(cur_subj==subject_ids_arr);
	subj_mat=all_factors[subj_ix, c(SubjectIDColName, TimeOffsetColName, target_variable_list)];

	sort_ix=order(subj_mat[,TimeOffsetColName], decreasing=F);
	subj_mat_sorted=subj_mat[sort_ix,,drop=F];

	#print(subj_mat_sorted);
	var_by_subj[[cur_subj]]=list();

	par(mfrow=c(multiplot_dim_r, multiplot_dim_c));
	for(targ_var_ix in target_variable_list){

		# Plot Variables
		#plot(subj_mat_sorted[,TimeOffsetColName], subj_mat_sorted[,targ_var_ix],
		#	main=paste(cur_subj, " : ", targ_var_ix, sep=""),
		#	xlab=TimeOffsetColName,
		#	ylab=targ_var_ix,
		#	ylim=c(ranges_mat[targ_var_ix, "min"], ranges_mat[targ_var_ix, "max"])
		#);

		var_by_subj[[cur_subj]][[targ_var_ix]]=subj_mat_sorted[,c(TimeOffsetColName, targ_var_ix),drop=F];

		plot_var(
			times=subj_mat_sorted[,TimeOffsetColName], 
			values=subj_mat_sorted[,targ_var_ix],
			subject_id=cur_subj,
			var_name=targ_var_ix,
			val_lim=c(ranges_mat[targ_var_ix, "min"], ranges_mat[targ_var_ix, "max"])
		);
		

		# Calculate parameters
		if(Opt_FindLine){
			acc_matrix[cur_subj, paste(targ_var_ix, "_", 
				c("intercept", "slope", "sd_res"), sep="")]=
				calc_line(subj_mat_sorted[,c(TimeOffsetColName, targ_var_ix),drop=F]);
		}

		if(Opt_FindLimits){
			acc_matrix[cur_subj, paste(targ_var_ix, "_",
				c("min", "max", "range"), sep="")]=
				calc_limits(subj_mat_sorted[,c(TimeOffsetColName, targ_var_ix),drop=F]);
		}

		if(Opt_FindDescriptive){
			acc_matrix[cur_subj, paste(targ_var_ix, "_", 
				c("mean", "stdev", "median"), sep="")]=
				calc_descriptive(subj_mat_sorted[,c(TimeOffsetColName, targ_var_ix),drop=F]);
		}

		if(Opt_FindRanges){
			acc_matrix[cur_subj, paste(targ_var_ix, "_", 
				c("first", "last", "N"), sep="")]=
				calc_ranges(subj_mat_sorted[,c(TimeOffsetColName, targ_var_ix),drop=F]);

		}
		if(Opt_FindTiming){
			acc_matrix[cur_subj, paste(targ_var_ix, "_", 
				c("timespan", "freq"), sep="")]=
				calc_timing(subj_mat_sorted[,c(TimeOffsetColName, targ_var_ix),drop=F]);
		}

	}	

	mtext(text=cur_subj, side=3, line=0, outer=T, cex=2, font=2);
}


##############################################################################

generate_group_plots=function(group_info, stats_mat, var_name, grp_cols){

	cat("Generating Grouped Plots:\n");
	print(group_info);
	print(stats_mat);

	num_groups=length(group_info);
	group_ids=names(group_info);
	cat("Num Groups: ", num_groups, "\n");
	
	stats_names=colnames(stats_mat);
	for(stats_ix in stats_names){
		cat("\n\nStat Name: ", stats_ix, "\n");
		grouped_data=list();
		medians=list();
		pvalue=list();

		min=min(stats_mat[, stats_ix]);
		max=max(stats_mat[, stats_ix]);
		range=max-min;

		# Perform calculations
		for(grp_ix in group_ids){

			in_group=group_info[[grp_ix]];
			out_group=numeric();
			for(grp_ix_2 in setdiff(group_ids, grp_ix)){
				out_group=c(out_group, group_info[[grp_ix_2]]);
			}

			print(in_group);
			print(out_group);

			grp_data=stats_mat[in_group, stats_ix]; 
			outgrp_data=stats_mat[out_group, stats_ix];

			wilcox_res=wilcox.test(grp_data, outgrp_data);
			pvalue[[grp_ix]]=wilcox_res$p.value;
	
			grouped_data[[grp_ix]]=grp_data;
			medians[[grp_ix]]=median(grp_data);
		}
	
		# Generate Plots

		new_page=par()$page;

		plot(0,0, type="n", main=gsub(paste(var_name, "_", sep=""), "", stats_ix),
			xlim=c(1-1, num_groups+1), ylim=c(min-(range*.05), max+(range*.05)),
			xaxt="n"
			);

		if(new_page){
			mtext(var_name, outer=T, font=2);
		}

		axis(side=1, at=c(1:num_groups), labels=group_ids, cex.axis=.75);

		# 
		for(grp_ix in 1:num_groups){
			num_pts=length(grouped_data[[grp_ix]]);

			# Annotate median
			points(c(grp_ix-.20, grp_ix+.20), rep(medians[[grp_ix]], 2),
				type="l", lwd=2, col="black"
				);

			# Draw datapoints
			points(rep(grp_ix, num_pts), grouped_data[[grp_ix]], col=grp_cols[grp_ix], lwd=2);
			points(rep(grp_ix, num_pts), grouped_data[[grp_ix]], col="black", lwd=.5);

			adj.gly=c(-.05,-.2);
			adj.pvl=c(-.05,.4);	
			pts=.65;
			gly_mult=2;
		
			if(pvalue[[grp_ix]]<=0.001){

				text(grp_ix+.25, medians[[grp_ix]], "****",
					cex=pts*gly_mult, adj=adj.gly);
				text(grp_ix+.25, medians[[grp_ix]], sprintf("%0.4f", pvalue[[grp_ix]]), 
					cex=pts, adj=adj.pvl);

			}else if(pvalue[[grp_ix]]<=0.01){

				text(grp_ix+.25, medians[[grp_ix]], "*** ",
					cex=pts*gly_mult, adj=adj.gly);
				text(grp_ix+.25, medians[[grp_ix]], sprintf("%0.3f ", pvalue[[grp_ix]]), 
					cex=pts, adj=adj.pvl);

			}else if(pvalue[[grp_ix]]<=0.05){

				text(grp_ix+.25, medians[[grp_ix]], "**  ",
					cex=pts*gly_mult, adj=adj.gly);
				text(grp_ix+.25, medians[[grp_ix]], sprintf("%0.2f  ", pvalue[[grp_ix]]), 
					cex=pts, adj=adj.pvl);

			}else if(pvalue[[grp_ix]]<=0.1){

				text(grp_ix+.25, medians[[grp_ix]], "*   ",
					cex=pts*gly_mult, adj=adj.gly);
				text(grp_ix+.25, medians[[grp_ix]], sprintf("%0.1f   ", pvalue[[grp_ix]]), 
					cex=pts, adj=adj.pvl);

			}
		}

		
	}

	
}

if(GroupColName!=""){

	# Generate subject to group mapping
	group_mapping_matrix=all_factors[, c(GroupColName, SubjectIDColName)]; 
	group_map=list();
	for(subj_ix in unique_subject_ids){
		ix=min(which(subj_ix==group_mapping_matrix[,SubjectIDColName]));
		group_map[[subj_ix]]=group_mapping_matrix[ix, GroupColName];
	}
	print(group_map);
	
	# Get group info
	uniq_group_names=unique(group_mapping_matrix[,GroupColName]);
	num_uniq_groups=length(uniq_group_names);
	cat("Groups:\n");
	print(uniq_group_names);
	cat("Num Groups: ", num_uniq_groups, "\n");

	# Generate group members structure
	group_members=list();
	for(grp_ix in uniq_group_names){
		group_members[[grp_ix]]=list();
	}


	for(subj_ix in unique_subject_ids){
		subj_grp=group_map[[subj_ix]];
		grp_mem=group_members[[subj_grp]];

		if(length(grp_mem)){
			grp_mem=c(grp_mem, subj_ix);
		}else{
			grp_mem=subj_ix;
		}

		group_members[[subj_grp]]=grp_mem;

	}

	cat("Group Members:\n");
	print(group_members);

	#---------------------------------------------------------------------
	#var_by_subj[[cur_subj]][[targ_var_ix]]

	grp_plot_rows=ceiling(sqrt(num_uniq_groups+1));
	grp_plot_cols=grp_plot_rows;
	if(grp_plot_cols*(grp_plot_rows-1)>(num_uniq_groups+1)){
		grp_plot_rows=grp_plot_rows-1;
	}

	colors=c("blue", "red", "yellow", "purple", "orange", "green");
	num_colors=length(colors);
	if(num_uniq_groups>num_colors){
		colors=rainbow(num_uniq_groups);
	}
	colors=colors[1:num_uniq_groups];
	names(colors)=uniq_group_names;

	for(targ_var_ix in target_variable_list){
		cat("Current Target Variables:", targ_var_ix, "\n");

		par(mfrow=c(grp_plot_rows, grp_plot_cols));

		# plot by group
		grp_data=list();
		for(grp_ix in uniq_group_names){

			# Create empty plot
			plot(0,0, type="n", 
				main=grp_ix,
				xlim=c(time_ranges_min, time_ranges_max),
				ylim=c(ranges_mat[targ_var_ix, "min"], ranges_mat[targ_var_ix, "max"]));

			grp_members=group_members[[grp_ix]];

			grp_data[[grp_ix]]=numeric();

			# Draw lines with colors
			for(subj_ix in grp_members){
				sbj_dat=var_by_subj[[subj_ix]][[targ_var_ix]];
				grp_data[[grp_ix]]=rbind(grp_data[[grp_ix]], sbj_dat);
	
				points(sbj_dat[, TimeOffsetColName], sbj_dat[, targ_var_ix], 
					col=colors[grp_ix],
					type="l", lwd=2);
			}

			# Draw thin black lines
			for(subj_ix in grp_members){
				sbj_dat=var_by_subj[[subj_ix]][[targ_var_ix]];
				grp_data[[grp_ix]]=rbind(grp_data[[grp_ix]], sbj_dat);
	
				points(sbj_dat[, TimeOffsetColName], sbj_dat[, targ_var_ix], 
					col="black",
					type="l", lwd=.5);
			}
		}

		# Generate combined plots
		plot(0,0, type="n",
			main="[Combined]",
			xlim=c(time_ranges_min, time_ranges_max),
			ylim=c(ranges_mat[targ_var_ix, "min"], ranges_mat[targ_var_ix, "max"]));

	
		# Compute lowess for groups
		lowess_res=list();
		for(grp_ix in uniq_group_names){
			lowess_res[[grp_ix]]=lowess(
				grp_data[[grp_ix]][,TimeOffsetColName], 
				grp_data[[grp_ix]][,targ_var_ix]);
		}

		# Draw lines with colors
		for(grp_ix in uniq_group_names){
			points(lowess_res[[grp_ix]][["x"]],
				lowess_res[[grp_ix]][["y"]], 
				col=colors[grp_ix],
				lwd=2,
				type="l");

		}
		# Draw thin black lines
		for(grp_ix in uniq_group_names){
			points(lowess_res[[grp_ix]][["x"]],
				lowess_res[[grp_ix]][["y"]], 
				col="black",
				lwd=.5,
				type="l");
		}

		mtext(text=targ_var_ix, side=3, line=0, outer=T, cex=2, font=2);

		#--------------------------------------------------------------

		cur_varlist=target_to_newvar_map[[targ_var_ix]];

		par(mfrow=c(multiplot_dim_r, multiplot_dim_c));

		generate_group_plots(group_members, acc_matrix[,cur_varlist], targ_var_ix, colors);

	}


}

###############################################################################

rounded_acc_mat=round(acc_matrix, 4);

if(GroupColName!=""){
	collapsed_group_ids=character();
	for(sbj_ix in unique_subject_ids){
		collapsed_group_ids=c(collapsed_group_ids, group_map[[sbj_ix]]);
	}

	out_hdr=c(SubjectIDColName, GroupColName, colnames(rounded_acc_mat));
	outmat=cbind(unique_subject_ids, collapsed_group_ids,  rounded_acc_mat);
}else{
	out_hdr=c(SubjectIDColName, colnames(rounded_acc_mat));
	outmat=cbind(unique_subject_ids, rounded_acc_mat);
}

colnames(outmat)=out_hdr;
write_factors(paste(OutputFname, ".longit.tsv", sep=""), outmat);

##############################################################################

cat("\nDone.\n");

print(warnings());
q(status=0);
