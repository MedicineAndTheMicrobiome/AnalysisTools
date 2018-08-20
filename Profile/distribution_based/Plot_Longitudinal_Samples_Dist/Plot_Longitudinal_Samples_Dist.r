#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"offset_file", "t", 1, "character",
	"output_file", "o", 2, "character",
	"diversity_type", "d", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

DEF_DIVERSITY="tail";

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input summary_table.tsv file>\n",
	"	-t <offset file>\n",
	"	[-o <output file root name>]\n",
	"	[-d <diversity, default=", DEF_DIVERSITY, ".]\n",
	"\n",
	"	This script will read in the summary table\n",
	"	and a file describing the time from the first\n",
	"	sample.\n",
	"\n",
	"	For each sample a time series bar plot and\n",
	"	and diversity will be generated.\n",
	"\n",
	"	The format of the offset file is:\n",
	"\n",
	"	<sample id> \\t <subj id> \\t <time stamp> \\t <cohort/grp id> \\t <events> [\\t <add'l events>]\\n\n",
	"\n",
	"	Column Names are required.\n",
	"\n",
	"	Diversity types include:\n",
	"		shannon, tail, simpson, invsimpson\n",
	"\n");

if(!length(opt$input_file) || !length(opt$offset_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OffsetFileName=opt$offset_file;

DiversityType=DEF_DIVERSITY;
if(length(opt$diversity_type)){
	DiversityType=opt$diversity_type;
}

if(length(opt$output_file)>0){
	OutputFileRoot=opt$output_file;
}else{
	OutputFileRoot=InputFileName;
	OutputFileRoot=gsub("\\.summary_table\\.tsv$", "", OutputFileRoot);
	OutputFileRoot=gsub("\\.summary_table\\.xls$", "", OutputFileRoot);
	cat("No output file root specified.  Using input file name as root.\n");
}

###############################################################################

OutputFileRoot=paste(OutputFileRoot, ".", substr(DiversityType, 1, 4), sep="");

OutputPDF = paste(OutputFileRoot, ".div_ts.pdf", sep="");
cat("Output PDF file name: ", OutputPDF, "\n", sep="");
pdf(OutputPDF,width=8.5,height=11)

###############################################################################

load_offset=function(fname){
        cat("Loading Offsets: ", fname, "\n");
        offsets_mat=read.delim(fname,  header=TRUE, row.names=1, sep="\t", comment.char="#", quote="");

	num_col=ncol(offsets_mat);
	cat("Num Columns Found: ", num_col, "\n");
	
	extra_colnames=colnames(offsets_mat);
	print(extra_colnames);
	colnames(offsets_mat)=c("Indiv ID", "Offsets", "Group ID", extra_colnames[4:num_col])[1:num_col];

	# reset offsets
	if(is.numeric(offsets_mat[,"Indiv ID"])){
		offsets_mat[,"Indiv ID"]=paste("#", offsets_mat[,"Indiv ID"], sep="");
	}
	groups=unique(offsets_mat[,"Indiv ID"]);
	
	cat("Groups:\n");
	print(groups);
	cat("\n");

	# Reset offsets so they are relative to the first/smallest sample
	for(gid in groups){
		group_ix=(gid==offsets_mat[,"Indiv ID"]);
		offsets=offsets_mat[group_ix, "Offsets"];
		min_off=min(offsets);
		offsets_mat[group_ix, "Offsets"]=offsets-min_off;
	}

	offsets_data=list();
	offsets_data[["matrix"]]=offsets_mat;
	offsets_data[["IndivID"]]=extra_colnames[1];
	offsets_data[["Offsets"]]=extra_colnames[2];
	offsets_data[["GroupID"]]=extra_colnames[3];
	
	return(offsets_data);
}

load_summary_file=function(fname){
        cat("Loading Summary Table: ", fname, "\n");
        inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", row.names=1))
        counts_mat=inmat[,2:(ncol(inmat))];
        return(counts_mat);
}

normalize=function(counts){
        totals=apply(counts, 1, sum);
        num_samples=nrow(counts);
        normalized=matrix(0, nrow=nrow(counts), ncol=ncol(counts));

        for(i in 1:num_samples){
                normalized[i,]=counts[i,]/totals[i];
        }

        colnames(normalized)=colnames(counts);
        rownames(normalized)=rownames(counts);
        return(normalized);
}

simplify_matrix_categories=function(normalized_mat, top=4){
	keep_cat=character();
	num_samp=nrow(normalized_mat);
	samp_names=rownames(normalized_mat);
	
	for(i in 1:num_samp){
		#cat(samp_names[i], "\n");
		abund=sort(normalized_mat[i,], decreasing=T);
		top_cat=(names(abund)[1:top]);	
		#print(top_cat);
		#cat("\n");
		keep_cat=c(keep_cat, top_cat);
	}

	uniq_keep_cat=unique(keep_cat);
	cat("Top ", top, " across all samples:\n", sep="");
	print(uniq_keep_cat);

	# Keep categories across top categories
	keep_mat=normalized_mat[,uniq_keep_cat];

	# Sort categories in decreasing order
	avg_abund=apply(keep_mat, 2, mean);
	sort_ix=order(avg_abund, decreasing=T);
	keep_mat=keep_mat[, sort_ix];
	return(keep_mat);
	
}

###############################################################################

# Get color assignments
get_colors=function(num_col, alpha=1){
	colors=hsv(seq(0,1,length.out=num_col+1), c(1,.5), c(1,.75,.5), alpha=alpha);
	color_mat_dim=ceiling(sqrt(num_col));
	color_pad=rep("grey", color_mat_dim^2);
	color_pad[1:num_col]=colors[1:num_col];
	color_mat=matrix(color_pad, nrow=color_mat_dim, ncol=color_mat_dim);
	colors=as.vector(t(color_mat));
	colors=colors[colors!="grey"];
}

###############################################################################

plot_dist=function(x, y, width=20, abundances){
	
	rect(
		xleft=x-width/2,
		ybottom=0,
		xright=x+width/2,
		ytop=1,
		lwd=.01,
		col="grey"
	);

	num_abund=length(abundances);
	prev=0;
	for(i in 1:num_abund){
		rect(
			xleft=x-width/2,
			ybottom=prev,
			xright=x+width/2,
			ytop=prev+abundances[i],
			lwd=.01,
			col=i
		);	
		prev=prev+abundances[i];
	}
		
}

###############################################################################

plot_event=function(x, y, events, col="black"){
	nl=length(levels(events));
	nsamp=length(events);
	if(nl==0){
		pts=which(events==1);
		points(x[pts], y[pts], type="p", pch="*", font=2, cex=1.5, col=col);
		points(x[pts], y[pts], type="p", pch="*", cex=1, col="black");
	}else{
		text(x, y, as.character(events), col=col, cex=.7);
	}
}

###############################################################################

plot_sample_distributions_by_individual=function(diversity_arr, div_type, normalized_mat, 
	offsets_mat, col_assign, category_colors, ind_colors){

	sorted_sids=sort(rownames(offsets_mat));
	offsets_mat=offsets_mat[sorted_sids,];

	# Get Unique Groups
	groups=sort(unique(offsets_mat[,"Indiv ID"]));
	num_groups=length(groups);

	def_par=par(no.readonly=T);
	par(mfrow=c(3,1));

	# Get range of offsets
	offset_ranges=range(offsets_mat[,"Offsets"]);
	cat("Offset Range:\n");
	print(offset_ranges);

	diversity_ranges=range(diversity_arr);
	cat("Diversity Range:\n");
	print(diversity_ranges);

	# Get closest between offsets
	periods=numeric();
	for(i in 1:num_groups){
		grp_subset=which(offsets_mat[,"Indiv ID"]==groups[i]);
                offsets=offsets_mat[grp_subset, "Offsets"];
		periods=c(diff(sort(offsets)), periods);
	}
	periods_sorted=sort(unique(periods));
	if(periods_sorted[1]==0){
		min_period=periods_sorted[2];
	}else{
		min_period=periods_sorted[1];
	}
	cat("Minimum period: ", min_period, "\n");

	cat_names=colnames(normalized_mat);

	events_range=4:ncol(offsets_mat);
	num_event_types=length(events_range);;
	
	# Plot individual samples
	for(i in 1:num_groups){

		# Grab members from group
		cat("Plotting: ", groups[i], "\n");
		grp_subset=which(offsets_mat[,"Indiv ID"]==groups[i]);
		num_members=length(grp_subset);

		# Subset offsets, and sort by offset
		offset_info=offsets_mat[grp_subset,];
		sort_ix=order(offset_info[,"Offsets"]);
		offset_info=offset_info[sort_ix,];
		events_mat=offset_info[, events_range, drop=F];

		#print(offset_info);

		###############################################################

		# Subset diversity
		subset_samples=rownames(offset_info);
		subset_diversity=diversity_arr[subset_samples];
		#print(subset_diversity);


		# Plot Diversity lines
		palette(ind_colors);
		plot(offset_info[,"Offsets"], subset_diversity, main=groups[i],
			xlab="Time", ylab=paste("Diversity (", div_type, ")", sep=""), type="l", 
			col=col_assign[groups[i]], lwd=2.5,
			xlim=offset_ranges, 
			ylim=c(diversity_ranges[1], diversity_ranges[2]+(diff(diversity_ranges)*.2))
		);

		char_height=par()$cxy[2]*.75;

		points(offset_info[c(1,1, num_members),"Offsets"], subset_diversity[c(1,1, num_members)], 
			type="p", pch=c(17, 1, 15), cex=c(1, 2, 1.25));
		points(offset_info[,"Offsets"], subset_diversity, type="b", pch=16, cex=.5, lwd=.1);

		# Mark Events
		for(e_ix in 1:num_event_types){
			plot_event(offset_info[,"Offsets"], subset_diversity+char_height*e_ix, events_mat[, e_ix]);
		}

		###############################################################

		# Plot abundances bar plots
		palette(category_colors);
		subset_norm=normalized_mat[subset_samples,];
		plot(offset_info[,"Offsets"], subset_diversity, main=groups[i],
			 xlab="Time", ylab="", type="n", col=i, lwd=2,
			 yaxt="n",
			 xlim=offset_ranges, ylim=c(0,1+.25));

		abd_pos=c(0, .25, .33, .5, .66, .75, 1);
		axis(side=2, at=abd_pos, labels=sprintf("%2.2f", abd_pos), las=2, cex.axis=.75);
		abline(h=c(0, 1), lwd=.5, col="grey75");
		abline(h=c(.25, .5, .75), lwd=.4, col="grey85");
		abline(h=c(.33, .66), lwd=.4, col="grey95");

		char_height=par()$cxy[2]*.75;

		# Plot stacked abundances
		for(t in 1:num_members){
			plot_dist(offset_info[subset_samples[t],"Offsets"], y=0, 
				abundances=normalized_mat[subset_samples[t],], width=min_period);
		}

		# Mark Events
		for(e_ix in 1:num_event_types){
			plot_event(offset_info[subset_samples, "Offsets"], 
				rep(1.075+(e_ix-1)*char_height, num_members), events_mat[, e_ix]);
		}

		###############################################################

		# Plot color key/legend
		num_in_key=min(35, length(cat_names));
		plot(0, 0, main=groups[i],
			 xlab="", ylab="", type="n", col=i, lwd=2,
			 xlim=c(0,10), ylim=c(0,num_in_key), xaxt="n", yaxt="n");

		for(j in 1:num_in_key){
			rect(0, j-1, .9, j+.9-1, col=j, lwd=.1);
			text(1, j-1+.4, labels=cat_names[j], pos=4, cex=.8);
		}
		
	
	}
	par(def_par);
}

###############################################################################

plot_sample_diversity_by_group=function(diversity_arr, div_type, normalized_mat, 
	offsets_mat, col_assign, ind_colors, grp_name=""){

	sorted_sids=sort(rownames(offsets_mat));
	offsets_mat=offsets_mat[sorted_sids,];

	# Get Num Cohorts
	cohorts=sort(unique(offsets_mat[,"Group ID"]));	
	num_cohorts=length(cohorts);
	cat("Number of Cohorts: ", num_cohorts, "\n");
	print(cohorts);
	cat("\n");

	# Get range of offsets
	offset_ranges=range(offsets_mat[,"Offsets"]);
	cat("Offset Range:\n");
	print(offset_ranges);

	# Get range of diversity
	diversity_ranges=range(diversity_arr);
	cat("Diversity Range:\n");
	print(diversity_ranges);

	# Set up plots per page
	def_par=par(no.readonly=T);
	par(mfrow=c(num_cohorts,1));

	# Set palette for individuals
	palette(ind_colors);

	events_range=4:ncol(offsets_mat);
	num_event_types=length(events_range);

	for(g in 1:num_cohorts){

		cat("Plotting: ", cohorts[g], "\n");
		plot(0, 0, main=paste(grp_name, ": ", cohorts[g], sep=""),
			 xlab="Time", ylab=paste("Diversity (", div_type,")",sep=""), type="n", 
			 xlim=offset_ranges, ylim=diversity_ranges);

		coh_offset_mat=offsets_mat[ offsets_mat[,"Group ID"]==cohorts[g], ];
		events_mat=coh_offset_mat[, events_range, drop=F];
		#print(coh_offset_mat);

		# Get Unique Inidividuals
		indivs=sort(unique(coh_offset_mat[,"Indiv ID"]));
		num_indivs=length(indivs);
		cat("Number of Individuals: ", num_indivs, "\n");
		print(indivs);
		cat("\n");

		# Plot individual samples
		for(i in 1:num_indivs){

			# Grab from group
			cat("Plotting: ", indivs[i], "\n");
			ind_subset=which(coh_offset_mat[,"Indiv ID"]==indivs[i]);
			num_timepts=length(ind_subset);

			# Subset offsets, and sort by offset
			offset_info=coh_offset_mat[ind_subset,];
			sort_ix=order(offset_info[,"Offsets"]);
			offset_info=offset_info[sort_ix,];
			#print(offset_info);

			###############################################################

			# Subset diversity
			subset_samples=rownames(offset_info);
			subset_diversity=diversity_arr[subset_samples];
			#print(subset_diversity);

			# Plot Diversity

			points(offset_info[c(1,1, num_timepts),"Offsets"], subset_diversity[c(1,1, num_timepts)], 
				type="p", pch=c(17, 1, 15), cex=c(1, 2, 1.25), col=col_assign[indivs[i]]);
			points(offset_info[,"Offsets"], subset_diversity, type="l", lwd=2.5, col=col_assign[indivs[i]]);
			points(offset_info[,"Offsets"], subset_diversity, type="l", lwd=.1, col="black");

			# Mark Events
			#mark_event=rep(0, num_timepts);
			#mark_event[offset_info[, "Events"]==T]=2;
			#points(offset_info[,"Offsets"], subset_diversity, type="p", pch="*", 
		#		font=2, cex=mark_event*1.5, col=col_assign[indivs[i]]);
			#points(offset_info[,"Offsets"], subset_diversity, type="p", pch="*", 
			#	cex=mark_event);

		}
	}
	par(def_par);
}

###############################################################################

bootstrap_med=function(x, num_bs=2000){

	meds=numeric(num_bs);
	for(i in 1:num_bs){
		meds[i]=median(sample(x, replace=T));
	}
	return(meds);
}

plot_median_sample_diversity_by_group=function(diversity_arr, div_type, 
	normalized_mat, offsets_mat, col_assign, ind_colors, grp_name=""){

	sorted_sids=sort(rownames(offsets_mat));
	offsets_mat=offsets_mat[sorted_sids,];

	def_par=par(no.readonly=T);

	# Get Num Cohorts
	cohorts=sort(unique(offsets_mat[,"Group ID"]));	
	num_cohorts=length(cohorts);
	cat("Number of Cohorts: ", num_cohorts, "\n");
	print(cohorts);
	cat("\n");

	cohort_data=list();

	# Get range of offsets
	offset_pos=sort(unique(offsets_mat[,"Offsets"]));
	offset_ranges=range(offset_pos);
	num_offset_positions=length(offset_pos);
	cat("Offset Range:\n");
	print(offset_ranges);
	cat("Offset Positions:\n");
	print(offset_pos);
	cat("Num Offset Positions:\n");
	num_offset_positions=length(offset_pos);	

	# Get range of diversity
	diversity_ranges=range(diversity_arr);
	cat("Diversity Range:\n");
	print(diversity_ranges);
	cat("\n\n");


	# Build data structure
	for(g in 1:num_cohorts){

		cat("Extracting: ", cohorts[g], "\n");

		coh_offset_mat=offsets_mat[ offsets_mat[,"Group ID"]==cohorts[g], ];

		# Get Unique Inidividuals
		indivs=sort(unique(coh_offset_mat[,"Indiv ID"]));
		num_indivs=length(indivs);
		cat("Number of Subjects: ", num_indivs, "\n");
		cat("Subjects:\n");
		print(indivs);
		cat("\n");

		cohort_id=as.character(cohorts[g]);
		cohort_data[[cohort_id]]=list();

		# Plot individual samples
		for(i in 1:num_indivs){

			ind_subset=which(coh_offset_mat[,"Indiv ID"]==indivs[i]);

			# Subset offsets, and sort by offset
			offset_info=coh_offset_mat[ind_subset,];
			sort_ix=order(offset_info[,"Offsets"]);
			offset_info=offset_info[sort_ix,];

			indiv_id=as.character(indivs[i]);
			cohort_data[[g]][[indiv_id]]=offset_info;
		}
	}

	print(cohort_data);

	medians=list();
	loess=list();
	samps=list();
	
	# Compute Medians and compute loess per cohort
	for(g in 1:num_cohorts){

		cohort_id=as.character(cohorts[g]);
		indiv_ids=names(cohort_data[[cohort_id]]);
		num_indv=length(indiv_ids);

		# per cohort medians across time
		med_mat=matrix(NA, nrow=num_offset_positions, ncol=3);
		colnames(med_mat)=c("Median", "LB_95CI", "UB_95CI");
		rownames(med_mat)=as.character(offset_pos);

		# per cohort div vs time
		div_offs_mat=matrix(NA, nrow=0, ncol=2);
		colnames(div_offs_mat)=c("Offsets", "Diversity");
		

		for(offset in offset_pos){
			div_arr=numeric(num_indv)
			
			for(i in 1:num_indv){

				idv=indiv_ids[i];
				ind_off=cohort_data[[cohort_id]][[idv]]		
				samp_names=rownames(ind_off);	
				targ_off_ix=which(ind_off[,"Offsets"]==offset);
				cat(samp_names[targ_off_ix], "\n");

				div_arr[i]=diversity_arr[samp_names[targ_off_ix]];
				div_offs_mat=rbind(div_offs_mat, c(offset, div_arr[i]));
			}	

			cat("Offset: ", offset, "\n");


			bs_med=bootstrap_med(div_arr);
			qtl=quantile(bs_med, na.rm=T, probs=c(0.025, .5, 0.975));
			med_mat[as.character(offset), c("LB_95CI", "Median", "UB_95CI")]=qtl;
		}

		#print(div_offs_mat);
		samps[[cohort_id]]=div_offs_mat;

		# Fit loess	
		offset_tmp=div_offs_mat[,1];
		diversity_tmp=div_offs_mat[,2];

		loess_res=loess(diversity_tmp~offset_tmp);
		loess_fit_x=sort(unique(c(offset_pos, 
					seq(offset_ranges[1], offset_ranges[2], length.out=num_offset_positions*2))));
		loess_fit_y=predict(loess_res, loess_fit_x);
		loess[[cohort_id]]=cbind(loess_fit_x, loess_fit_y);

		medians[[cohort_id]]=med_mat;

	}

	print(medians);
	print(loess);
	print(samps);

	# Set up plots per page
	def_par=par(no.readonly=T);
	par(mfrow=c(3,1));

	colors_full  =rainbow(num_cohorts, start=0, end=4/6, s=1);
	colors_pastel=rainbow(num_cohorts, start=0, end=4/6, s=.3, alpha=.85);

	#----------------------------------------------------------------------
	# Plot Medians and PI's
	plot(0,0, type="n",
		xlim=c(offset_ranges[1], offset_ranges[2]),
		ylim=c(0, diversity_ranges[2]*1.2),
		xlab="Time",
		ylab="Median Diversity (w/ 95% CI)"
	);

	chwd=par()$cxy[1]/4;
	for(g in 1:num_cohorts){
		cohort_id=as.character(cohorts[g]);
		cur_med=medians[[cohort_id]]		

		# Prediction intervals
		x=as.numeric(rownames(cur_med));
		for(off_ix in 1:length(x)){
			y_lb=cur_med[off_ix,"LB_95CI"];
			y_ub=cur_med[off_ix,"UB_95CI"];

			xpos=x[off_ix];

			points(c(xpos-chwd, xpos+chwd), c(y_lb, y_lb), type="l", col=colors_pastel[g]);
			points(c(xpos-chwd, xpos+chwd), c(y_ub, y_ub), type="l", col=colors_pastel[g]);

			points(c(xpos,xpos), c(y_lb, y_ub), type="l", col=colors_pastel[g], lwd=2);
		}

		# Medians
		x=as.numeric(rownames(cur_med));
		y=cur_med[,"Median"];
		cat("[x,y]=[",x, ", ", y,"]\n");

		points(x, y, lwd=.5, type="l", col=colors_full[g]);
		points(x, y, type="p", col=colors_full[g]);
		
		

	}

	#----------------------------------------------------------------------
	# Plot loess fit
	plot(0,0, type="n",
		xlim=c(offset_ranges[1], offset_ranges[2]),
		ylim=c(0, diversity_ranges[2]*1.2),
		xlab="Time",
		ylab="Loess Fitted Diversity"
	);

	# Plot individuals points
	for(g in 1:num_cohorts){
		cohort_id=as.character(cohorts[g]);
		cur_samps=samps[[cohort_id]];

		x=cur_samps[,"Offsets"];
		y=cur_samps[,"Diversity"];
		points(x, y, type="p", cex=1, col=colors_pastel[g]);
	}
	
	# Draw loess line
	for(g in 1:num_cohorts){
		cohort_id=as.character(cohorts[g]);
		cur_loess=loess[[cohort_id]];

		x=cur_loess[,"loess_fit_x"];
		y=cur_loess[,"loess_fit_y"];
		points(x, y, type="l", col=colors_full[g], lwd=2);
	}

	#----------------------------------------------------------------------
	# Plot group legend
	plot(0,0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", bty="n", xaxt="n", yaxt="n");
	legend(0, 1, cohorts, fill=colors_full, bty="n", cex=2, title=grp_name);
	
	par(def_par);

}

###############################################################################

plot_sample_distributions_by_group=function(normalized_mat, offsets_mat, cat_colors, grp_name=""){

	sorted_sids=sort(rownames(offsets_mat));
	offsets_mat=offsets_mat[sorted_sids,];

	# Get Num Cohorts
	cohorts=sort(unique(offsets_mat[,"Group ID"]));	
	num_cohorts=length(cohorts);
	cat("\nNumber of Cohorts: ", num_cohorts, "\n");
	print(cohorts);
	cat("\n");

	# Get range of offsets
	offset_ranges=range(offsets_mat[,"Offsets"]);
	cat("Offset Range:\n");
	print(offset_ranges);

	# Get closest between offsets
	periods=numeric();
	indivs=unique(offsets_mat[,"Indiv ID"]);
	num_indivs=length(indivs);
	for(i in 1:num_indivs){
		grp_subset=which(offsets_mat[,"Indiv ID"]==indivs[i]);
                offsets=offsets_mat[grp_subset, "Offsets"];
		periods=c(diff(sort(offsets)), periods);
	}
	periods_sorted=sort(unique(periods));

	if(periods_sorted[1]==0){
		min_period=periods_sorted[2];
	}else{
		min_period=periods_sorted[1];
	}
	cat("Minimum period: ", min_period, "\n");

	palette(cat_colors);
	events_range=4:ncol(offsets_mat);
	num_event_types=length(events_range);

	# Set up plots per page
	def_par=par(no.readonly=T);
	for(g in 1:num_cohorts){

		coh_offset_mat=offsets_mat[ offsets_mat[,"Group ID"]==cohorts[g], ];
		#print(coh_offset_mat);

		# Get Unique Inidividuals
		indivs=sort(unique(coh_offset_mat[,"Indiv ID"]));
		num_indivs=length(indivs);
		cat("Number of Individuals: ", num_indivs, "\n");
		print(indivs);
		cat("\n");

		par(mfrow=c(num_indivs,1));
		par(oma=c(3,1,2,0));
		par(mar=c(0.1, 4.1, 0.1, 2.1)); # bot, left, top, right

		# Plot individual samples
		for(i in 1:num_indivs){

			# Grab from group
			cat("\n\nPlotting: ", as.character(indivs[i]), "\n");
			ind_subset=which(coh_offset_mat[,"Indiv ID"]==indivs[i]);
			num_timepts=length(ind_subset);
			cat("Num Time Pts: ", num_timepts, "\n");

			# Subset offsets, and sort by offset
			offset_info=coh_offset_mat[ind_subset,];
			sort_ix=order(offset_info[,"Offsets"]);
			offset_info=offset_info[sort_ix,];
			events_mat=offset_info[, events_range, drop=F];
			#print(offset_info);

			###############################################################
			if(i==num_indivs){
				xaxt_setting="s";	
			}else{
				xaxt_setting="n";
			}
		
			plot(0, 0, main="",
				xlab="", ylab=indivs[i], type="n", bty="n",
				xaxt=xaxt_setting, yaxt="n",
				xlim=offset_ranges, ylim=c(0, 1+.5*num_event_types));

			char_height=par()$cxy[2]*.5;

			subset_samples=rownames(offset_info);

			# Draw guide lines
			abline(h=.5, col="grey", lwd=20);
			xaxp=par()$xaxp;
			xguides=seq(xaxp[1], xaxp[2], (xaxp[2]-xaxp[1])/xaxp[3]);
			#print(xguides);
			abline(v=xguides, col="grey", lwd=.5);

			# Draw stacked barplot
			for(t in 1:num_timepts){
				plot_dist(offset_info[t,"Offsets"], y=0, 
					abundances=normalized_mat[subset_samples[t],], width=min_period);

			}

			# Mark Events
			for(e_ix in 1:num_event_types){
				plot_event(offset_info[subset_samples, "Offsets"], 
					rep(1.15+(e_ix-1)*char_height, num_timepts), events_mat[, e_ix]);
			}
		}

		mtext(paste(grp_name, ": ", cohorts[g], sep=""), side=3, outer=T, font=2);
	}
	par(def_par);
}

###############################################################################

plot_mean_distributions_by_group=function(normalized_mat, offsets_mat, cat_colors, grp_name=""){


	num_categories=ncol(normalized_mat);

	sorted_sids=sort(rownames(offsets_mat));
	offsets_mat=offsets_mat[sorted_sids,];

	# Get Num Cohorts
	cohorts=sort(unique(offsets_mat[,"Group ID"]));	
	num_cohorts=length(cohorts);
	cat("\nNumber of Cohorts: ", num_cohorts, "\n");
	print(cohorts);
	cat("\n");

	# Get range of offsets
	uniq_offsets=sort(unique(offsets_mat[,"Offsets"]));
	num_uniq_offsets=length(uniq_offsets);
	offset_ranges=range(uniq_offsets);
	cat("Offset Range:\n");
	print(offset_ranges);

	# Get closest between offsets
	periods=numeric();
	indivs=unique(offsets_mat[,"Indiv ID"]);
	num_indivs=length(indivs);
	for(i in 1:num_indivs){
		grp_subset=which(offsets_mat[,"Indiv ID"]==indivs[i]);
                offsets=offsets_mat[grp_subset, "Offsets"];
		periods=c(diff(sort(offsets)), periods);
	}
	periods_sorted=sort(unique(periods));

	if(periods_sorted[1]==0){
		min_period=periods_sorted[2];
	}else{
		min_period=periods_sorted[1];
	}
	cat("Minimum period: ", min_period, "\n");

	palette(cat_colors);
	events_range=4:ncol(offsets_mat);
	num_event_types=length(events_range);


	cat("Offsets:\n");
	print(uniq_offsets);
	cat("\n");

	# Set up plots per page
	def_par=par(no.readonly=T);
	par(oma=c(0,0,2,0));
	par(mfrow=c(num_cohorts+1, 1));

	for(g in 1:num_cohorts){

		# Extract offset for cohort
		coh_offset_mat=offsets_mat[offsets_mat[,"Group ID"]==cohorts[g],,drop=F];
		#print(coh_offset_mat);

		# Save mean abundances in matrix
		avg_off_comp=matrix(0, nrow=num_uniq_offsets, ncol=num_categories);
		rownames(avg_off_comp)=as.numeric(uniq_offsets);
		colnames(avg_off_comp)=colnames(normalized_mat);

		# Save num samples per offset in array
		samp_size_at_offset=numeric(num_uniq_offsets);
		names(samp_size_at_offset)=uniq_offsets;

		# Calculate means and sample size per offset
		for(cur_off in uniq_offsets){
			#cat("Offset Val: ", cur_off, "\n");
			off_mat_ix=coh_offset_mat[,"Offsets"]==cur_off;
			subj_id=rownames(coh_offset_mat[off_mat_ix,,drop=F]);

			offset_key=as.character(cur_off);
			avg_off_comp[offset_key, ]=apply(normalized_mat[subj_id,,drop=F], 2, mean);
			samp_size_at_offset[offset_key]=length(subj_id);
		
		}

		#--------------------------------------------------------------
		# Plot times series for cohort
		plot(0, 0, main="",
			xlab="", 
			ylab=paste(grp_name, ": ", cohorts[g]), 
			type="n", bty="n",
			xaxt="s", yaxt="n",
			xlim=offset_ranges, ylim=c(0, 1.1));

		# Draw guide lines
		abline(h=.5, col="grey", lwd=20);
		xaxp=par()$xaxp;
		cxy=par()$cxy;
		xguides=seq(xaxp[1], xaxp[2], (xaxp[2]-xaxp[1])/xaxp[3]);
		abline(v=xguides, col="grey", lwd=.5);

		# Draw stacked barplot
		for(t in 1:num_uniq_offsets){
			cur_off=uniq_offsets[t];
			offset_key=as.character(cur_off);
			plot_dist(cur_off, y=0, 
				abundances=avg_off_comp[offset_key,,drop=F], width=min_period);

			text(cur_off, 1+cxy[2]/2, 
				paste("n=",samp_size_at_offset[offset_key], sep=""), cex=.75, font=3);

		}

		#--------------------------------------------------------------

	}

	mtext("Mean Abundances by Group", side=3, outer=T, font=2);

	# Plot color key/legend
	cat_names=colnames(normalized_mat);
	num_in_key=min(35, length(cat_names));
	plot(0, 0, 
		 xlab="", ylab="", type="n", lwd=2,
		 xlim=c(0,10), ylim=c(0,num_in_key), xaxt="n", yaxt="n");

	for(j in 1:num_in_key){
		rect(0, j-1, .9, j+.9-1, col=j, lwd=.1);
		text(1, j-1+.4, labels=cat_names[j], pos=4, cex=.8);
	}

	par(def_par);

}

###############################################################################

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T){
	num_row=nrow(mat);
	num_col=ncol(mat);

	cat("Num Rows: ", num_row, "\n");
	cat("Num Cols: ", num_col, "\n");

	mat=mat[rev(1:num_row),];

	num_colors=50;
	color_arr=rainbow(num_colors, start=0, end=4/6);
	if(high_is_hot){
		color_arr=rev(color_arr);
	}
	
	remap=function(in_val, in_range, out_range){
		in_prop=(in_val-in_range[1])/(in_range[2]-in_range[1])	
		out_val=in_prop*(out_range[2]-out_range[1])+out_range[1];
		return(out_val);
	}

	if(is.na(plot_min)){
		plot_min=min(mat);	
	}
	if(is.na(plot_max)){
		plot_max=max(mat);
	}
	cat("Plot min/max: ", plot_min, "/", plot_max, "\n");

	par(oma=c(12,12,0,1));
	par(mar=c(0, 0, 2, 0));
	plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), xaxt="n", yaxt="n", bty="n", xlab="", ylab="", main=title);
	
	# x-axis
	axis(side=1, at=seq(.5, num_col-.5, 1), labels=colnames(mat), las=2);
	axis(side=2, at=seq(.5, num_row-.5, 1), labels=rownames(mat), las=2);

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
			text(x-.5, y-.5, sprintf("%0.4f", mat[y,x]), srt=45);
		}
	}
	
}

###############################################################################

plot_change_scatter=function(diversity_arr, offset_mat, grp_name=""){

	cat("Plotting Change Scatter:\n");
	#print(diversity_arr);
	#print(offset_mat);

	trt_levels=unique(offset_mat[,"Group ID"]);
	ind_ids=unique(offset_mat[,"Indiv ID"]);
	
	num_indiv=length(ind_ids);
	num_trt=length(trt_levels);

	cat("Treatment Levels: \n");
	print(trt_levels);
	cat("Individual IDs: \n");
	print(ind_ids);

	# Extract start/end diversity
	ends=matrix(NA, nrow=num_indiv, ncol=2, dimnames=list(ind_ids, c("start", "end")));
	for(ind_ix in ind_ids){
		offset_subset=offset_mat[offset_mat[,"Indiv ID"]==ind_ix,];
		offset_subset=offset_subset[order(offset_subset[,"Offsets"], decreasing=F),];
		samp_names=rownames(offset_subset)
		num_offsets=nrow(offset_subset);
		ends[ind_ix, "start"]=diversity_arr[samp_names[1]];
		ends[ind_ix, "end"]=diversity_arr[samp_names[num_offsets]];
	}

	end_ranges=range(ends);	
	axis_ticks=round(seq(end_ranges[1], end_ranges[2], length.out=10), 3);

	grp_col_transp=get_colors(num_trt, alpha=.2);
	grp_col_opaque=get_colors(num_trt, alpha=1);


	par(mar=c(5,6,3,2));
	# Plot by treatment group
	par(mfrow=c(num_trt+1, 1));
	for(trt_ix in 1:num_trt){
		trt=trt_levels[trt_ix];
		trt_ind_ids=unique(offset_mat[trt==offset_mat[,"Group ID"], "Indiv ID"]);
		trt_ends=ends[trt_ind_ids,, drop=F];

		#print(trt_ends);
		plot(0,0, type="n",
			main=paste(grp_name, ": ", trt, sep=""),
			xlab="", ylab="",
			xlim=end_ranges, ylim=end_ranges, 
			xaxt="n", yaxt="n");
		axis(side=1, at=axis_ticks, label=axis_ticks);
		axis(side=2, at=axis_ticks, label=axis_ticks, las=2);
		abline(a=0, b=1, col="grey", lwd=2);
		title(xlab="Start", line=3);
		title(ylab="End", line=4);

		median_end=median(trt_ends[,"end"]);
		median_start=median(trt_ends[,"start"]);

		abline(h=median_end, col=grp_col_opaque[trt_ix]);
		abline(v=median_start, col=grp_col_opaque[trt_ix]);

		axis(side=3, at=median_start, label=round(median_start, 3), 
			line=-1, tick=F, cex.axis=.75, col.axis=grp_col_opaque[trt_ix]);
		axis(side=4, at=median_end, label=round(median_end, 3), 
			line=-1, tick=F, cex.axis=.75, col.axis=grp_col_opaque[trt_ix]);

		points(trt_ends, col=grp_col_transp[trt_ix], cex=2, pch=19);
		points(trt_ends, col=grp_col_opaque[trt_ix], cex=.25, pch=19);
		
	}

	# Plot groups in single plot
	plot(0,0, type="n",
		main="Combined",
		xlab="", ylab="",
		xlim=end_ranges, ylim=end_ranges,
		xaxt="n", yaxt="n");
	abline(a=0, b=1, col="grey", lwd=2);
	axis(side=1, at=axis_ticks, label=axis_ticks);
	axis(side=2, at=axis_ticks, label=axis_ticks, las=2);
	title(xlab="Start", line=3);
	title(ylab="End", line=4);


	all_ends_and_groups=list();
	for(trt_ix in 1:num_trt){
                trt=trt_levels[trt_ix];
                trt_ind_ids=unique(offset_mat[trt==offset_mat[,"Group ID"], "Indiv ID"]);
                trt_ends=ends[trt_ind_ids,, drop=F];

		all_ends_and_groups[[paste(trt, ":start", sep="")]]=trt_ends[,"start"];
		all_ends_and_groups[[paste(trt, ":end", sep="")]]=trt_ends[,"end"];

		median_end=median(trt_ends[,"end"]);
		median_start=median(trt_ends[,"start"]);

		abline(h=median_end, col=grp_col_opaque[trt_ix]);
		abline(v=median_start, col=grp_col_opaque[trt_ix]);
		axis(side=3, at=median_start, label=round(median_start, 3), line=-1, 
			tick=F, cex.axis=.75, col.axis=grp_col_opaque[trt_ix]);
		axis(side=4, at=median_end, label=round(median_end, 3), line=-1, 
			tick=F, cex.axis=.75, col.axis=grp_col_opaque[trt_ix]);

		points(trt_ends, col=grp_col_transp[trt_ix], cex=2, pch=19);
		points(trt_ends, col=grp_col_opaque[trt_ix], cex=.25, pch=19);
	}

	#######################################################################
	# Compute pvalues between groups

	aeag_names=names(all_ends_and_groups);

	for(aeag_idx in aeag_names){
		cat(aeag_idx, ": ", median(all_ends_and_groups[[aeag_idx]]), "\n", sep="");	
	}

	pval_matrix=matrix(NA, nrow=num_trt*2, ncol=num_trt*2, dimnames=list(aeag_names, aeag_names));
	diff_matrix=matrix(NA, nrow=num_trt*2, ncol=num_trt*2, dimnames=list(aeag_names, aeag_names));
	for(aeag_idx_A in aeag_names){
		for(aeag_idx_B in aeag_names){
			res=wilcox.test(all_ends_and_groups[[aeag_idx_A]], all_ends_and_groups[[aeag_idx_B]]);
			pval_matrix[aeag_idx_A, aeag_idx_B]=res$p.value;
			diff_matrix[aeag_idx_A, aeag_idx_B]=(
				abs(diff(c(median(c(all_ends_and_groups[[aeag_idx_A]])), 
				median(all_ends_and_groups[[aeag_idx_B]])))));
		}
	}
	#print(pval_matrix);
	#print(diff_matrix);

	par(mfrow=c(1,1))
	paint_matrix(pval_matrix, 0, 1, log_col=T, title="End Point Difference: p-values", high_is_hot=F);
	paint_matrix(diff_matrix, log_col=F, title="End Point Differences: Effect Sizes", high_is_hot=T);

}

tail_statistic=function(x){
        sorted=sort(x, decreasing=TRUE);
        norm=sorted/sum(x);
        n=length(norm);
        tail=0;
        for(i in 1:n){
                tail=tail + norm[i]*((i-1)^2);
        }
        return(sqrt(tail));
}

###############################################################################

plot_text=function(strings){

        orig_par=par(no.readonly=T);
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

        par(orig_par);
}

###############################################################################
###############################################################################

offset_data=load_offset(OffsetFileName);
offset_mat=offset_data[["matrix"]];
print(offset_data);

###############################################################################

counts_mat=load_summary_file(InputFileName);
#print(counts_mat);

###############################################################################

offset_mat_samples=rownames(offset_mat);
counts_mat_samples=rownames(counts_mat);
shared=intersect(offset_mat_samples, counts_mat_samples);

cat("\n\n");
cat("Samples not represented in summary table file:\n");
print(setdiff(counts_mat_samples, shared));
cat("Samples not represented in offsets file:\n");
print(setdiff(offset_mat_samples, shared));
cat("\n\n");

offset_mat=offset_mat[shared,];
counts_mat=counts_mat[shared,];
indiv_ids=unique(offset_mat[,"Indiv ID"]);
num_indiv=length(indiv_ids);

###############################################################################

normalized_mat=normalize(counts_mat);
#print(normalized_mat);

if(DiversityType=="tail"){
	diversity_arr=apply(normalized_mat, 1, tail_statistic);
}else{
	diversity_arr=diversity(normalized_mat, DiversityType);
}
#print(diversity_arr);

###############################################################################

# simplify matrix
simplified_mat=simplify_matrix_categories(normalized_mat, top=3);
num_simp_cat=ncol(simplified_mat);


category_colors=get_colors(num_simp_cat);
ind_colors=get_colors(num_indiv);

col_assign=1:num_indiv;
names(col_assign)=indiv_ids;


###############################################################################

plot_text(c(paste("By Subject (", offset_data[["IndivID"]], "):  Diversity and Composition", sep="")));
plot_sample_distributions_by_individual(diversity_arr, DiversityType, 
	simplified_mat, offset_mat, col_assign, category_colors, ind_colors);

plot_text(c(paste("By Group (", offset_data[["GroupID"]], "):  Diversity", sep="")));
plot_sample_diversity_by_group(diversity_arr, DiversityType, simplified_mat, 
	offset_mat, col_assign, ind_colors, grp_name=offset_data[["GroupID"]]);

plot_text(c(paste("By Group (", offset_data[["GroupID"]], "):  Median Sample Diversity", sep="")));
plot_median_sample_diversity_by_group(diversity_arr, DiversityType, simplified_mat, 
	offset_mat, col_assign, ind_colors, grp_name=offset_data[["GroupID"]]);

plot_text(c(paste("By Group (", offset_data[["GroupID"]], "):  Composition", sep="")));
plot_sample_distributions_by_group(simplified_mat, offset_mat, category_colors,
	grp_name=offset_data[["GroupID"]]);

plot_text(c(paste("By Group (", offset_data[["GroupID"]], "):  Mean Composition", sep="")));
plot_mean_distributions_by_group(simplified_mat, offset_mat, category_colors,
	grp_name=offset_data[["GroupID"]]);

plot_text(c(paste("By Group (", offset_data[["GroupID"]], "):  Delta Scatter Plots", sep="")));
plot_change_scatter(diversity_arr, offset_mat,
	grp_name=offset_data[["GroupID"]]);

##############################################################################

cat("Done.\n")
dev.off();
warn=warnings();
if(length(warn)){
	print(warn);
}
q(status=0)
