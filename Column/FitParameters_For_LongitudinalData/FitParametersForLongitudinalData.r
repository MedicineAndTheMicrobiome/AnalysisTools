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
	"find_timing", "n", 2, "character",
	"find_time_of_limits", "T", 2, "character",
	"find_extrapolated_limits", "L", 2, "character",
	"find_recovery_rate_model", "R", 2, "character",

	"crop", "c", 2, "character",
	"prefix", "p", 2, "character"
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
	"	[-T (find timing of limits: min_time/max_time)\n",
	"	[-L (find extrapolated limits: extrp_start, extrp_end, extrp_slope)\n",
	"	[-R (fit recovery_rate_model: start, max, recovery_rate, decline_rate)\n",
	"\n",
	"	Specify focus on subset of timepoints:\n",
	"	[--crop=<start>,<end>]\n",
	"	[--prefix=<variable name prefix, e.g. mon1, mon1to6, mon6to12>]\n",
	"\n",
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
Opt_FindTimeOfLimits=F;
Opt_FindExtrapolatedLimits=F;
Opt_FindRecoveryRateModel=F;

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
if(length(opt$find_time_of_limits)){
	Opt_FindTimeOfLimits=T;
}
if(length(opt$find_extrapolated_limits)){
	Opt_FindExtrapolatedLimits=T;
}
if(length(opt$find_recovery_rate_model)){
	Opt_FindRecoveryRateModel=T;
}


GroupColName="";
if(length(opt$group_cn)){
	GroupColName=opt$group_cn;
}

Crop="";
CropLimits=c(-Inf, Inf);
if(length(opt$crop)){
	Crop=opt$crop;
	CropLimits=as.numeric(strsplit(Crop, ",")[[1]]);
}

Prefix="";
if(length(opt$prefix)){
	Prefix=opt$prefix;
}


cat("\n");
cat("Input Filename:", InputFname , "\n");
cat("Time Offset Colname:", TimeOffsetColName, "\n");
cat("SubjectID Colname:", SubjectIDColName, "\n");
cat("Target Values List Filename:", TargetValuesFname, "\n");
cat("Output Filename:", OutputFname, "\n");
cat("Group Colname: ", GroupColName, "\n");
cat("\n");
cat("Crop Start/End: (", paste(CropLimits, collapse="/"), ")\n");
cat("Variable Prefix: ", Prefix, "\n");
cat("\n");
cat("Find Line: ", Opt_FindLine, "\n");
cat("Find Limits: ", Opt_FindLimits, "\n");
cat("Find Descriptive: ", Opt_FindDescriptive, "\n");
cat("Find Ranges: ", Opt_FindRanges, "\n");
cat("Find Timing: ", Opt_FindTiming, "\n");
cat("Find Timing of Limits: ", Opt_FindTimeOfLimits, "\n");
cat("Find Extrapolated Limits: ", Opt_FindExtrapolatedLimits, "\n");
cat("Find Recovery Rate Model: ", Opt_FindRecoveryRateModel, "\n");

if(!(Opt_FindLine || Opt_FindLimits || Opt_FindDescriptive || 
	Opt_FindRanges || Opt_FindTiming || Opt_FindTimeOfLimits ||
	Opt_FindExtrapolatedLimits || Opt_FindRecoveryRateModel
)){
	cat("At least one of the find/fit parameters should be specified.\n");
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

plot_text=function(strings, no_touch_par=F, lines=52, cex=1){

        orig.par=par(no.readonly=T);
        par(family="Courier");

        if(!no_touch_par){
                par(oma=rep(.1,4));
                par(mar=rep(0,4));
        }

        num_lines=length(strings);

        top=max(as.integer(num_lines), lines);

        plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                );

        text_size=max(.01, min(.8, .8 - .003*(num_lines-52)))*cex;
        #print(text_size);

        for(i in 1:num_lines){
                #cat(strings[i], "\n", sep="");
                strings[i]=gsub("\t", "", strings[i]);
                text(0, top-i, strings[i], pos=4, cex=text_size);
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

		min=min(y, na.rm=T);
		max=max(y, na.rm=T);
		range=max-min;

		return(c(min, max, range));
	}
}
if(Opt_FindDescriptive){
	extensions=c(extensions, c("mean", "stdev", "median"));
	calc_descriptive=function(data){

		x=data[,1];
		y=data[,2];

		mean=mean(y, na.rm=T);
		stdev=sd(y, na.rm=T);
		median=median(y, na.rm=T);

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
if(Opt_FindTimeOfLimits){
	extensions=c(extensions, c("time_of_min", "time_of_max"));
	calc_timing_of_limits=function(data){
	
		minx=min(data[,2], na.rm=T);
		maxx=max(data[,2], na.rm=T);

		min_ix=min(which(data[,2]==minx));
		max_ix=min(which(data[,2]==maxx));
		
		min_time=data[min_ix,1];
		max_time=data[max_ix,1];

		return(c(min_time, max_time));

	}
}          

if(Opt_FindExtrapolatedLimits){
	extensions=c(extensions, c("extrp_start", "extrp_mid", "extrp_end", "extrp_slope"));

	calc_extrapolated_limits=function(data, beginx, endx){

		x=data[,1];
		y=data[,2];

		num_time_pts=nrow(data);

		if(num_time_pts<2){
			return(rep(NA, 4));
		}

		ord_ix=order(x);
		data=data[ord_ix,];

		earliest_y=data[1,2];
		last_y=data[num_time_pts, 2];

		end_weighted_data=rbind(
			c(beginx, earliest_y),
			data,
			c(endx, last_y)
		);

		ew_x=end_weighted_data[,1];
		ew_y=end_weighted_data[,2];

		lmfit=lm(ew_y~ew_x);
		intercept=lmfit$coefficients["(Intercept)"];		
		slope=lmfit$coefficients["ew_x"];		

		extrp_start=slope*beginx + intercept;
		extrp_end=slope*endx + intercept;
		extrp_mid=(extrp_end+extrp_start)/2;
		extrp_slope=slope;
		
		return(c(extrp_start, extrp_mid, extrp_end, extrp_slope));

	}
}

if(Opt_FindRecoveryRateModel){

	extensions=c(extensions, c("rrm_start", "rrm_max", "rrm_exp_rate", "rrm_lin_rate"));

	# This model consists of a exponential (short-term) and a linear (long-term)
	# recovery/deline rate.  The parameters allow for a non-zero starting and a maximum
	# recovery state/value.  This model was fit for lung transplant/FEV1 data, where
	# a subject would be expected to have a non-zero FEV1 after the surgery, then
	# improve with exponentially decreases of improvement to a maximum.  In many subjects,
	# there is a subsequent decline after the maximum recovery has been achieved.

	exp_lin_rate_model=function(t, start_val, max_val, exp_rate, lin_rate){
		y=start_val + (max_val-start_val)*2*(1/(1+exp(-t*exp_rate))-.5) + t*lin_rate;
		return(y);
	}

	#t=0:1000; plot(t, exp_lin_rate_model(t, 50, 100, .3, -.03), ylim=c(0, 200));

	calc_recovery_rate_params=function(data){
		time=data[,1];
		y=data[,2];
		num_time_pts=length(time);

		if(num_time_pts<5){
			return(rep(NA, 4));
		}

		# Order the data by time
		by_time=order(time);
		time=time[by_time];
		y=y[by_time];

		# Pre-calculations for estimating starting parameters
		max_y=max(y);
		min_y=min(y);
		time_of_max_y=min(which(max_y==y));

		last_half_time=time[num_time_pts]-median(time);
		
		init_start_val=y[1];
		init_max_val=max_y;
		init_lin_rate=sd(y)/last_half_time;
		init_exp_rate=0;

		#lin_fit=lm(y~time);
		#y_intercept=lin_fit$coefficients["(Intercept)"];

		# Define objective function
		model_ssd=function(p){
			obs=y;
			pred=exp_lin_rate_model(
				t=time, start_val=p[1], max_val=p[2], exp_rate=p[3], lin_rate=p[4]);
			ssd=sum(abs(obs-pred)^2);
			return(ssd);
		}

		cat("Starting Parameter Values:\n");
		print(c(init_start_val, init_max_val, init_exp_rate, init_lin_rate));

		# Calc normalization for each parameter
		# Set the scale for the exp_rate to be the same as for the lin_rate
		sum_param=sum(abs(c(init_start_val, init_max_val, init_lin_rate, init_lin_rate)));
		psc=c(init_start_val, init_max_val, init_lin_rate, init_lin_rate)/sum_param;

		cat("Parameter Scaling:\n");
		print(psc);
		
		# Run optimization
		opt_res=optim(
			par=c(init_start_val, init_max_val, 0, init_lin_rate),
			control=list(parscale=psc),
			lower=c(0, min_y/2, 0, -Inf),
			upper=c(y[1], max_y, Inf, Inf),
			fn=model_ssd, method="L-BFGS-B"
		);

		print(opt_res);

		return(opt_res$par);

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


trim_data=function(in_data, limits){
	x=in_data[,1];

	keep_ix=(x>=limits[1]) & (x<=limits[2]);
	num_kept=sum(keep_ix);

	if(num_kept==0){
		zmat=matrix(numeric(), nrow=0, ncol=2);
		colnames(zmat)=colnames(in_data);	
		return(zmat);
	}

	kept_data=in_data[keep_ix,,drop=F];
	return(kept_data);
}

pdf(file=paste(OutputFname, ".longit.pdf", sep=""), height=8.5, width=11);
	
plot_title_page("Individual Plots", c(
	"Longitudinal plots for each of the variables are",
	"generated for each subject.",
	"",
	"The subject ID is labeled in the top margin of each page.",
	"Red and blue glyphs demark the max and min values for",
	"each subject across time points.  The dotted line",
	"represents a line fit through the data points.",
	"The y and x axes represent each variable's value and",
	"time, respectively."
));

par(oma=c(0,0,2,0));
par(mar=c(2,2,3,1));

time_ranges_min=min(all_factors[,TimeOffsetColName]);
time_ranges_max=max(all_factors[,TimeOffsetColName]);

time_of_variable_arr=character();

for(cur_subj in unique_subject_ids){

	cat("Working on: ", cur_subj, "\n");

	subj_ix=(cur_subj==subject_ids_arr);
	subj_mat=all_factors[subj_ix, c(SubjectIDColName, TimeOffsetColName, target_variable_list)];

	sort_ix=order(subj_mat[,TimeOffsetColName], decreasing=F);
	subj_mat_sorted=subj_mat[sort_ix,,drop=F];

	#print(subj_mat_sorted);
	var_by_subj[[cur_subj]]=list();

	par(mfrow=c(multiplot_dim_r, multiplot_dim_c));
	for(targ_var_ix in target_variable_list){

		var_spec_tab=subj_mat_sorted[,c(TimeOffsetColName, targ_var_ix),drop=F];
		nona_ix=!is.na(var_spec_tab[,2]);
		
		var_by_subj[[cur_subj]][[targ_var_ix]]=var_spec_tab[nona_ix,];

		plot_var(
			times=subj_mat_sorted[,TimeOffsetColName], 
			values=subj_mat_sorted[,targ_var_ix],
			subject_id=cur_subj,
			var_name=targ_var_ix,
			val_lim=c(ranges_mat[targ_var_ix, "min"], ranges_mat[targ_var_ix, "max"])
		);
		
		target_data=subj_mat_sorted[,c(TimeOffsetColName, targ_var_ix), drop=F];
		times=target_data[,1];
		if(!is.null(CropLimits)){
			cat("Trimming data to: ", CropLimits, "\n");
			target_data=trim_data(target_data, CropLimits);
			print(target_data);
			if(is.finite(CropLimits[1])){
				abline(v=CropLimits[1], col="orange", lwd=1.5);
				mtext(paste(CropLimits[1], " >", sep=""), 
					at=CropLimits[1], side=1, line=-1, col="darkorange", adj=1);
			}
			if(is.finite(CropLimits[2])){
				abline(v=CropLimits[2], col="orange", lwd=1.5);
				mtext(paste("< ", CropLimits[2], sep=""),
					at=CropLimits[2], side=1, line=-1, col="darkorange", adj=0);
			}
		}

		#----------------------------------------------------------------------------
		# Calculate parameters
		if(Opt_FindLine){
			acc_matrix[cur_subj, paste(targ_var_ix, "_", 
				c("intercept", "slope", "sd_res"), sep="")]=
				calc_line(target_data);
		}

		if(Opt_FindLimits || Opt_FindTimeOfLimits){
			acc_matrix[cur_subj, paste(targ_var_ix, "_",
				c("min", "max", "range"), sep="")]=
				calc_limits(target_data);
		}

		if(Opt_FindDescriptive){
			acc_matrix[cur_subj, paste(targ_var_ix, "_", 
				c("mean", "stdev", "median"), sep="")]=
				calc_descriptive(target_data);
		}

		if(Opt_FindRanges){
			acc_matrix[cur_subj, paste(targ_var_ix, "_", 
				c("first", "last", "N"), sep="")]=
				calc_ranges(target_data);

		}
		if(Opt_FindTiming){
			acc_matrix[cur_subj, paste(targ_var_ix, "_", 
				c("timespan", "freq"), sep="")]=
				calc_timing(target_data);
		}
		if(Opt_FindTimeOfLimits){

			varnames=paste(targ_var_ix, "_", c("time_of_min", "time_of_max"), sep="");

			acc_matrix[cur_subj, varnames]=
				calc_timing_of_limits(target_data);

			time_of_variable_arr=unique(c(time_of_variable_arr, varnames));
		}

		if(Opt_FindExtrapolatedLimits){

			varnames=paste(targ_var_ix, "_", 
				c("extrp_start", "extrp_mid", "extrp_end", "extrp_slope"), sep="");

			params=calc_extrapolated_limits(target_data,
					CropLimits[1], CropLimits[2]
					);
			acc_matrix[cur_subj, varnames]=params;
		
			t=seq(CropLimits[1], CropLimits[2], length.out=100);
			intercept=params[1]-params[4]*CropLimits[1];

			points(t, params[4]*t+intercept, type="b", col="green", cex=.5);
			points(c(CropLimits[1], CropLimits[2]), c(params[1], params[2]), 
				type="p", cex=1.1, col="green");
			

		}

		if(Opt_FindRecoveryRateModel){

			varnames=paste(targ_var_ix, "_", 
				c("rrm_start", "rrm_max", "rrm_exp_rate", "rrm_lin_rate"),
				sep="");

			params=calc_recovery_rate_params(target_data);
			acc_matrix[cur_subj, varnames]=params;

			start_plot_time=max(CropLimits[1], 0);
			end_plot_time=min(CropLimits[2], max(times));

			t=sort(c(
				seq(start_plot_time, end_plot_time/16, length.out=100),
				seq(start_plot_time, end_plot_time/8, length.out=100),
				seq(start_plot_time, end_plot_time/4, length.out=100),
				seq(start_plot_time, end_plot_time, length.out=200)
			));

			points(t, exp_lin_rate_model(t, params[1], params[2], params[3], params[4]), 
				type="b", col="blue", cex=.25);
			
		}

		#----------------------------------------------------------------------------

		plot_par=par();	
		left=(plot_par$usr[2]-plot_par$usr[1])*.01 + plot_par$usr[1];
		top=(plot_par$usr[4]-plot_par$usr[3])*.9 + plot_par$usr[3];

		msg=paste(
			capture.output({print(t(acc_matrix[cur_subj, varnames, drop=F]), quote=F)}),
			collapse="\n");

		text(left, top, msg, pos=4, family="mono");


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

		min=min(stats_mat[, stats_ix], na.rm=T);
		max=max(stats_mat[, stats_ix], na.rm=T);
		range=max-min;

		# Perform calculations
		for(grp_ix in group_ids){

			in_group=group_info[[grp_ix]];
			out_group=numeric();
			for(grp_ix_2 in setdiff(group_ids, grp_ix)){
				out_group=c(out_group, group_info[[grp_ix_2]]);
			}

			grp_data=stats_mat[in_group, stats_ix]; 
			outgrp_data=stats_mat[out_group, stats_ix];

			wilcox_res=wilcox.test(grp_data, outgrp_data);

			pvalue[[grp_ix]]=wilcox_res$p.value;

			if(is.na(pvalue[[grp_ix]])){
				pvalue[[grp_ix]]=1.0;
			}
	
			grouped_data[[grp_ix]]=grp_data;
			medians[[grp_ix]]=median(grp_data, na.rm=T);
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

		jitter_list=list();
		for(grp_ix in 1:num_groups){
			jitter_list[[grp_ix]]=rnorm(length(grouped_data[[grp_ix]]), 0, .1);
		}
		 
		for(grp_ix in 1:num_groups){
			num_pts=length(grouped_data[[grp_ix]]);

			# Annotate median
			points(c(grp_ix-.20, grp_ix+.20), rep(medians[[grp_ix]], 2),
				type="l", lwd=2, col="black"
				);

			# Draw datapoints
			points(rep(grp_ix, num_pts)+jitter_list[[grp_ix]], 
				grouped_data[[grp_ix]], col=grp_cols[grp_ix], lwd=2);
			points(rep(grp_ix, num_pts)+jitter_list[[grp_ix]], 
				grouped_data[[grp_ix]], col="black", lwd=.5);

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

	plot_title_page("Group Plots", c(
		"The Group Plots (Longitudinal and Extracted Statistics) are generated in pairs for each variable.",
		"The variable name is labeled in the top margin of each page.",
		"",
		"Longitudinal Plots:", 
		"",
		"The longitudinal data for each subjected are overlaid together by their group membership",
		"and assigned a color.  Individual points in the plot are label with their subject ID.",
		"The last longitudinal plot contains the lowess curve representing each group overlaid together.",
		"",
		"",
		"Extracted Statistics Plots:",
		"",
		"These strip plots by group compare the statistics calculated for each subject's longitudinal data.",
		"Each glyph in the strip plot represents a different subject.  The horizontal bar",
		"represents the median value for each group.  If a group's median is significantly",
		"different from the median of the subjects not in the group, the p-value is labeled.",
		"P-values were calculated with the one-sample Wilcoxon test.",
		"",
		"intercept, slope, and sd_res (standard deviation of the residuals) were calculated with linear regression.",
		"min, max, range, mean, stdev, and median were calculated independent of time.",
		"timespan and freq (samples per time unit) consider the sampling time.",
		"time_of_min and time_of_max measure the earliest time point when a subject's",
		"    sample has reach a minimum or maximum in value, respectively."
	));

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

				# Draw colored lines
				points(sbj_dat[, TimeOffsetColName], sbj_dat[, targ_var_ix], 
					col=colors[grp_ix],
					type="l", lwd=2);

				# Draw points
				points(sbj_dat[, TimeOffsetColName], sbj_dat[, targ_var_ix], 
					col="black",
					type="p", cex=.7);

				# Label with subject id
				text(sbj_dat[, TimeOffsetColName], sbj_dat[, targ_var_ix],
					pos=3,
					subj_ix, col="grey20", cex=.4);
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

		num_curvar=length(cur_varlist);	
		var_stat_dim_c=ceiling(sqrt(num_curvar));	
		var_stat_dim_r=var_stat_dim_c;
		if(var_stat_dim_c*(var_stat_dim_r-1) > num_curvar){
			var_stat_dim_r=var_stat_dim_r-1;
		}
		
		par(mfrow=c(var_stat_dim_r, var_stat_dim_c));

		generate_group_plots(group_members, acc_matrix[,cur_varlist], targ_var_ix, colors);

	}


	###############################################################################
	###############################################################################

	if(length(time_of_variable_arr)){

		cat("Plotting Min/Max Time Lines:\n");

		plot_title_page("Timeline Plots", c(
			"Timeline plots place the variables along a timeline (x-axis) to provide context",
			"for when each variable's levels are peaking or troughing.",	
			"",
			"Plots are generated at various p-value cutoffs and for each group.",
			"For each group, a combined (red+blue), and then a mins (blue) and maxs (red)",
			"only plot are generated.  Since multiple variables may peak or trough at",
			"the sample time point, callouts are used to provide sufficient spacing",
			"between variables for legibility.",
			"",
			"The p-value cutoffs are used to exclude ploting a variable's min/max",
			"if there is not a significant difference between them.  The p-values",
			"use a Bonferroni correction, with then number of tests equal to the number",
			"of variables."
		));

		# Plot story lines for time_of variables
		#print(colors);
		#print(group_members);
		#print(time_of_variable_arr);
		#print(acc_matrix);
		#quit();

		tov_vnames=character();
		for(var_ix in target_variable_list){
			for(stat_ix in c("time_of_min", "time_of_max", "min", "max")){
				tov_vnames=c(tov_vnames, paste(var_ix, "_", stat_ix, sep=""));
			}
		}

		var2stat=function(varname, stat_name){
			paste(varname, "_", stat_name, sep="");
		}

		tov_matrix=acc_matrix[unique_subject_ids, tov_vnames]

		#print(tov_matrix);
		grp_medians_rec=list();
		groups=names(group_members);

		max_time=0;
		num_comparisons=num_uniq_groups*num_target_variables;

		# Compute means for the "time of" variables by group
		for(grp_ix in groups){

			grp_medians_rec[[grp_ix]]=list();

			for(var_ix in target_variable_list){

				grp_time_of_min=acc_matrix[
					group_members[[grp_ix]],
					var2stat(var_ix, "time_of_min"), drop=F];

				grp_time_of_max=acc_matrix[
					group_members[[grp_ix]],
					var2stat(var_ix, "time_of_max"), drop=F];

				grp_mins=acc_matrix[
					group_members[[grp_ix]],
					var2stat(var_ix, "min"), drop=F];

				grp_maxs=acc_matrix[
					group_members[[grp_ix]],
					var2stat(var_ix, "max"), drop=F];

				#cat("mins:\n");
				#print(grp_mins);
				#cat("maxs:\n");
				#print(grp_maxs);
				wc.res=wilcox.test(grp_mins, grp_maxs);			
				if(!is.finite(wc.res$p.value)){
					wc.res$p.value=1;
				}

				grp_medians_rec[[grp_ix]][[var2stat(var_ix, "time_of_min")]]=
					median(grp_time_of_min);
				grp_medians_rec[[grp_ix]][[var2stat(var_ix, "time_of_max")]]=
					median(grp_time_of_max);
				grp_medians_rec[[grp_ix]][[var2stat(var_ix, "range_signf")]]=
					wc.res$p.value;
				grp_medians_rec[[grp_ix]][[var2stat(var_ix, "range_signf_adj")]]=
					min(1, wc.res$p.value*num_comparisons);

				max_time=max(max_time, grp_time_of_max);
			}
		}


		par(mfrow=c(num_uniq_groups,1));

		#-----------------------------------------------------------------------------

		adjust_pos=function(orig_pos, range, char_h){

			num_pos=length(orig_pos);
			# Do not adjust if less than 2 positions
			if(num_pos<2){
				return(orig_pos);
			}
			
			orig_order=names(orig_pos);
			new_pos=sort(orig_pos);

			keep_adjusting=T;
			num_trials=0;
			max_trials=num_pos^2;
			while(keep_adjusting){

				keep_adjusting=F;

				# Move labels forward
				for(i in 1:(num_pos-1)){
					if(abs(new_pos[i]-new_pos[i+1])<char_h){
						adj=new_pos[i+1]+char_h/2;
						if(adj<=range[2]){
							new_pos[i+1]=adj;
							keep_adjusting=T;
						}
					}
				}

				# Move labels backward
				for(i in num_pos:2){
					if(abs(new_pos[i]-new_pos[i-1])<char_h){
						adj=new_pos[i-1]-char_h/2;
						if(adj>=range[1]){
							new_pos[i-1]=adj;
							keep_adjusting=T;
						}
					}
				}

				# If not enough space, make inter-label spacing smaller
				if(num_trials>num_pos){
					num_trials=0;
					char_h=char_h*.95;
					cat("Reducing spacing between labels: ", char_h, "\n");
				}

				num_trials=num_trials+1;

			}

			new_pos=new_pos[orig_order];
			return(new_pos);
		}

		#-----------------------------------------------------------------------------

		label=function(actual_pos, adjusted_pos, cxy_scale){

			arr_len=length(actual_pos);

			if(arr_len==0){return();};

			points(actual_pos, rep(0, arr_len), cex=.3);

			varnames=names(actual_pos);

			for(i in 1:arr_len){
				
				points(
					c(actual_pos[i], adjusted_pos[i]),
					c(0, .08),
					type="l", col="black", lwd=.5);

				label=gsub("_time_of", "", varnames[i]);
				label_length=nchar(label)*cxy_scale[2];
				label_size=2*min(1, 1/label_length);

				ismin=length(grep("_min$", label));

				label_col=ifelse(ismin, "blue", "red");

				text(adjusted_pos[i], .1, 
					label,
					adj=c(0,.25), srt=90, 
					col=label_col,
					cex=label_size);
			}
		}

		#-----------------------------------------------------------------------------
		print(grp_medians_rec);

		par(mfrow=c(3,1));
		for(cutoff_ix in c(1, 0.1, 0.05, 0.01, 0.001)){
			cat("Cutoff: ", cutoff_ix, "\n");

			for(grp_ix in groups){	
				cat("      Group: ", grp_ix, "\n");

				for(split_ix in c(F, T)){
					cat("   Split: ", split_ix, "\n");


					plot_title=paste(grp_ix, ", cutoff=", cutoff_ix, sep="");				
					# Start Plot
					#   We need to start the plot so we can get the character size calculations.
					#   The first plot is min when splitting, or both when combining.
					plot(NA, type="n", xlim=c(0, max_time), ylim=c(0,1), 
						yaxt="n",
						main="", xlab="", ylab="");
					title(main=
						paste(ifelse(split_ix, "Mins: ", "Mins & Maxs:"), plot_title), 
						cex.main=2, line=.5);
					cxy=par("cxy");

					# Get variables to plot based on cutoff
					plot_list=c();
					for(var_ix in target_variable_list){
						signf=grp_medians_rec[[grp_ix]][[var2stat(var_ix, "range_signf")]];
						if(signf<cutoff_ix){
							plot_list=c(plot_list, var_ix);
						}
					}

					var_min_time_arr=c();
					var_max_time_arr=c();
					var_both_time_arr=c();
					var_min_names_arr=c();
					var_max_names_arr=c();
					var_both_names_arr=c();
					for(var_ix in plot_list){
						tomin_name=var2stat(var_ix, "time_of_min");
						tomax_name=var2stat(var_ix, "time_of_max");

						var_min_names_arr=c(var_min_names_arr, tomin_name);
						var_max_names_arr=c(var_max_names_arr, tomax_name);
						var_both_names_arr=c(var_both_names_arr, c(tomin_name, tomax_name));

						tomin=grp_medians_rec[[grp_ix]][[tomin_name]];
						tomax=grp_medians_rec[[grp_ix]][[tomax_name]];
						var_min_time_arr=c(var_min_time_arr, tomin);
						var_max_time_arr=c(var_max_time_arr, tomax);
						var_both_time_arr=c(var_both_time_arr, c(tomin, tomax));
					}
					names(var_min_time_arr)=var_min_names_arr;
					names(var_max_time_arr)=var_max_names_arr;
					names(var_both_time_arr)=var_both_names_arr;

					# Adjust spacing of labels
					min_time_adjusted_pos=
						adjust_pos(var_min_time_arr, range=c(0, max_time), cxy[1]);
					max_time_adjusted_pos=
						adjust_pos(var_max_time_arr, range=c(0, max_time), cxy[1]);
					both_time_adjusted_pos=
						adjust_pos(var_both_time_arr, range=c(0, max_time), cxy[1]);

					if(split_ix){
						
						# Min
						label(var_min_time_arr, min_time_adjusted_pos, cxy);					
						# Max 
						plot(NA, type="n", xlim=c(0, max_time), ylim=c(0,1), 
							yaxt="n",
							main="", xlab="", ylab="");
								
						title(main=paste("Maxs: ", plot_title), 
							cex.main=2, line=.5);
						# Max	
						label(var_max_time_arr, max_time_adjusted_pos, cxy);					

					}else{
						# Both
						label(var_both_time_arr, both_time_adjusted_pos, cxy);					
					}
				}
			}

		}

	}

	#######################################################################
	#######################################################################

	cat("Plot Rankings...\n");

	plot_title_page("Extreme (Top/Bottom) Plots", c(
		"For selected statistics, groups with values that are significantly greater or less",
		"than the values from members of the other groups (out-groups) are marked on the matrix plot.",
		"",
		"If a group's values are greater or less than members from the out-groups, then a red up or",
		"blue down triangle glyph, respectively, is drawn at the intersection of the group name",
		"and the variable name.  It is possible for more than one group to be marked",
		"as higher or lower than members outside of their group.  When glyphs are not drawn",
		"across or down a group or variable, the group or variable is not significantly",
		"different from the values of other groups.",
		"",
		"Following each plot is a text list of the greater/less than relationships for ease",
		"of cutting/pasting into another document.",
		"",
		"",
		"slope: A slope suggests that over time, values are increasing or decreasing.  Remember",
		"that these plots only identify the greatest or least sloped groups, so a variable/group",
		"combination that is not marked could still have significant non-zero slopes.",
		"",
		"stdev: A greater stdev suggests a variable's values are more disrupted than the other groups.",
		"A group with consistently minimal stdev across all variables could be a characteristic",
		"of a control/no-treatment group."
	));

	find_signf_extremes=function(grp_mem, tar_var, acc_mat, tar_stat, pval_cutoff=0.1){

		grp_names=names(grp_mem);
		avail_subjects=rownames(acc_mat);	

		num_tar_var=length(tar_var);
		num_grps=length(grp_names);

		# Store results
		extreme_mat=matrix(0, nrow=num_grps, ncol=num_tar_var);
		extreme_pval_mat=matrix(numeric(), nrow=num_grps, ncol=num_tar_var);

		rownames(extreme_mat)=grp_names;
		colnames(extreme_mat)=tar_var;
		rownames(extreme_pval_mat)=grp_names;
		colnames(extreme_pval_mat)=tar_var;

		# Compare every variable for it's target statistic
		for(var_ix in tar_var){

			varstat_name=var2stat(var_ix, tar_stat);			
			values=acc_mat[,varstat_name];

			cat("Target Vars: ", varstat_name, "\n");
			
			# Compare each group against the rest
			for(grp_ix in grp_names){

				cat("Group: ", grp_ix, "\n");
				tar_subj=grp_mem[[grp_ix]];

				tar_val=values[tar_subj];
				oth_val=values[setdiff(avail_subjects, tar_subj)];

				# Compare
				wres=wilcox.test(tar_val, oth_val);

				if(median(tar_val) > median(oth_val)){
					rank=1;
				}else{
					rank=-1;
				}	

				extreme_pval_mat[grp_ix, var_ix]=wres$p.value;
				if(wres$p.value<pval_cutoff){
					extreme_mat[grp_ix, var_ix]=rank;
				}
			}
		}

		res=list();
		res[["rank"]]=extreme_mat;
		res[["pval"]]=extreme_pval_mat;
		return(res);

	}

	#----------------------------------------------------------------------

	plot_signf_extremes=function(extr_rec, title){

		ranks=extr_rec[["rank"]];
		num_vars=ncol(ranks);
		num_grps=nrow(ranks);


		par(mfrow=c(1,1));
		par(oma=c(0,1,2,0));
		par(mar=c(1,12,12,1));

		plot(0,0, type="n", xlim=c(0, num_vars+1), ylim=c(0, num_grps+1),
			xaxt="n", yaxt="n", xlab="", ylab="", main=""
			);

		mtext(title, side=3, outer=T, cex=3, font=2, line=-1, adj=0);

		# variables
		axis(side=3, at=1:num_vars, colnames(ranks), las=2);
		axis(side=2, at=1:num_grps, rownames(ranks), las=2);

		abline(v=1:num_vars, col="grey", lwd=.5);
		abline(h=1:num_grps, col="grey", lwd=.75);

		for(vix in 1:num_vars){
			for(gix in 1:num_grps){

				if(ranks[gix, vix]==1){
					rpch=24;
					rcol="red";
				}else if(ranks[gix, vix]==-1){
					rpch=25;
					rcol="blue"
				}
				
				if(ranks[gix, vix]!=0){
					points(vix, gix, pch=rpch, bg=rcol);
				}
			}
		}

		# Label 
		cxy=par()$cxy;

		points(0,.5, pch=24, bg="red");
		points(0,.5-cxy[2], pch=25, bg="blue");

		text(0,.5, pos=4, labels="Group's values greater than subjects in other groups.");
		text(0,.5-cxy[2], pos=4, labels="Group's values less than subjects in other groups.");

	}

	#----------------------------------------------------------------------

	list_signf_extremes=function(extr_rec, title){

		ranks=extr_rec[["rank"]];
                num_vars=ncol(ranks);
                num_grps=nrow(ranks);

		grps=rownames(ranks);
		vars=colnames(ranks);

		extr_list=list();

		for(grp_ix in grps){

			mins_arr=c();
			maxs_arr=c();
			for(var_ix in vars){
				r=ranks[grp_ix, var_ix];
				if(r==1){
					maxs_arr=c(maxs_arr, var_ix);
				}else if(r==-1){
					mins_arr=c(mins_arr, var_ix);
				}
			}

			extr_list[[grp_ix]]=list();
			extr_list[[grp_ix]][["Mins"]]=mins_arr;
			extr_list[[grp_ix]][["Maxs"]]=maxs_arr;

		}

		print(extr_list);

	}

	#----------------------------------------------------------------------

	if(Opt_FindLine){
	
		# Plot slope
		extr_rec=find_signf_extremes(group_members, target_variable_list, acc_matrix, "slope");
		plot_signf_extremes(extr_rec, "slope");
		slope_ext_list=list_signf_extremes(extr_rec, "slope");

		plot_text(c(
			"Slope:",
			"",
			capture.output(print(slope_ext_list, quote=F))
			));
	}

	if(Opt_FindDescriptive){

		# plot Stdev
		extr_rec=find_signf_extremes(group_members, target_variable_list, acc_matrix, "stdev");
		plot_signf_extremes(extr_rec, "stdev");
		stdev_ext_list=list_signf_extremes(extr_rec, "stdev");

		plot_text(c(
			"Stdev:",
			"",
			capture.output(print(stdev_ext_list, quote=F))
			));

	}

}

###############################################################################
# Output calculated variables into tsv

rounded_acc_mat=round(acc_matrix, 4);

out_colnames=colnames(rounded_acc_mat);
if(Prefix!=""){
	out_colnames=paste(Prefix, ".", out_colnames, sep="");
}

if(GroupColName!=""){
	# Add group id to output if specified
	collapsed_group_ids=character();
	for(sbj_ix in unique_subject_ids){
		collapsed_group_ids=c(collapsed_group_ids, group_map[[sbj_ix]]);
	}

	out_hdr=c(SubjectIDColName, GroupColName, out_colnames);
	outmat=cbind(unique_subject_ids, collapsed_group_ids,  rounded_acc_mat);
}else{
	out_hdr=c(SubjectIDColName, out_colnames);
	outmat=cbind(unique_subject_ids, rounded_acc_mat);
}

colnames(outmat)=out_hdr;
write_factors(paste(OutputFname, ".longit.tsv", sep=""), outmat);

##############################################################################

cat("\nDone.\n");

print(warnings());
q(status=0);
