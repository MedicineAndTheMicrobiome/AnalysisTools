#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');

source("~/git/AnalysisTools/Longitudinal/Longitudinal.r");

params=c(
	"summary_file", "s", 1, "character",

	"factor_file", "f", 1, "character",
	"offset_col", "t", 1, "character",
	"subject_id_col", "i", 1, "character",
	"output_root", "o", 1, "character",

	"diversity_type", "d", 2, "character",

	"group_col", "g", 2, "character",

        "dont_reset_offsets", "n", 2, "logical",
        "begin_offset", "b", 2, "numeric",
        "end_offset", "e", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

DEF_DIVERSITY="tail";

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <input summary_table.tsv file>\n",
	"\n",
	"	-f <factor file name>\n",
	"	-t <offset column name>\n",
	"	-i <subject identifier column name>\n", 
	"	-o <output file root name>\n",
	"\n",
	"	[-d <diversity, default=", DEF_DIVERSITY, ".]\n",
	"\n",
	"	[-g <group/cohort variable column name>]\n",
	"\n",
	"	[-n (do not reset earliest offsets to 0 to line up time points, default=reset offsets)]\n",
	"	[-b <begin offset, default=-Inf>]\n",
	"	[-e <end offset, default=Inf>]\n",
	"\n",
	"	Diversity types include:\n",
	"		shannon, tail, simpson, invsimpson\n",
	"\n");

if(
	!length(opt$summary_file) || 
	!length(opt$factor_file) || 
	!length(opt$offset_col) || 
	!length(opt$subject_id_col) || 
	!length(opt$output_root)
){
	cat(usage);
	q(status=-1);
}

# Required
SummaryFile=opt$summary_file;
FactorFile=opt$factor_file;
OutputRoot=opt$output_root;
SubjectIDCol=opt$subject_id_col;
TimeOffsetCol=opt$offset_col;

# Optional, i.e. with defaults
DiversityType=DEF_DIVERSITY;
GroupCol="";
ResetOffsets=T;
BeginOffset=-Inf;
EndOffset=Inf;

if(length(opt$diversity_type)){
	DiversityType=opt$diversity_type;
}
if(length(opt$group_col)){
	GroupCol=opt$group_col;
}
if(length(opt$dont_reset_offsets)){
	ResetOffsets=opt$dont_reset_offsets;
}
if(length(opt$begin_offset)){
	BeginOffset=opt$begin_offset;
}
if(length(opt$end_offset)){
	EndOffset=opt$end_offset;
}

###############################################################################

OutputRoot=paste(OutputRoot, ".", substr(DiversityType, 1, 4), sep="");

OutputPDF = paste(OutputRoot, ".div_ts.pdf", sep="");
cat("Output PDF file name: ", OutputPDF, "\n", sep="");
pdf(OutputPDF,width=8.5,height=11)

###############################################################################

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

load_factors=function(fname){
        factors=data.frame(read.table(fname,  sep="\t", header=TRUE, row.names=1,
                check.names=FALSE, comment.char=""));
        factor_names=colnames(factors);

        ignore_idx=grep("^IGNORE\\.", factor_names);

        if(length(ignore_idx)!=0){
                return(factors[-ignore_idx]);
        }else{
                return(factors);
        }
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
	offsets_rec, grp_to_sbj_info, col_assign, category_colors, ind_colors){

	cat("Running plot_sample_distributions_by_individual\n");

	offsets_mat=offsets_rec[["matrix"]];
	sorted_sids=offsets_rec[["SampleIDs"]];
	subject_ids=offsets_rec[["SubjectIDs"]];
	num_subjects=offsets_rec[["NumSubjects"]];

	def_par=par(no.readonly=T);
	par(mfrow=c(3,1));

	# Get range of offsets
	offset_ranges=range(offsets_rec[["Offsets"]]);
	cat("Offset Range:\n");
	print(offset_ranges);

	diversity_ranges=range(diversity_arr);
	cat("Diversity Range:\n");
	print(diversity_ranges);

	min_period=offsets_rec[["MinOffsetSep"]];
	cat("Minimum period: ", min_period, "\n");

	cat_names=colnames(normalized_mat);

	# Plot individual samples
	for(i in 1:num_subjects){

		# Grab members from group
		cat("Plotting: ", subject_ids[i], "\n");
		sbj_subset=which(offsets_mat[,"SubjectID"]==subject_ids[i]);

		grp_id=grp_to_sbj_info[["SbjToGrp"]][subject_ids[i]];
		#cat("Cohort ID:", cohort_id, "\n");

		num_members=length(sbj_subset);

		# Subset offsets, and sort by offset
		offset_info=offsets_rec[["OffsetsBySubject"]][[subject_ids[i]]];

		###############################################################

		# Subset diversity
		subset_samples=rownames(offset_info);
		subset_diversity=diversity_arr[subset_samples];
		#print(subset_diversity);

		title=paste(subject_ids[i], " [", grp_id, "]", sep="");

		# Plot Diversity lines
		palette(ind_colors);
		plot(offset_info[,"Offsets"], subset_diversity,
			main=paste(title, "Diversity"),
			xlab="Time", ylab=paste("Diversity (", div_type, ")", sep=""), type="l", 
			col=col_assign[subject_ids[i]], lwd=2.5,
			xlim=c(offset_ranges[1]-min_period/2, offset_ranges[2]+min_period/2), 
			ylim=c(diversity_ranges[1], diversity_ranges[2]+(diff(diversity_ranges)*.2))
		);

		char_height=par()$cxy[2]*.75;

		points(offset_info[c(1,1, num_members),"Offsets"], subset_diversity[c(1,1, num_members)], 
			type="p", pch=c(17, 1, 15), cex=c(1, 2, 1.25));
		points(offset_info[,"Offsets"], subset_diversity, type="b", pch=16, cex=.5, lwd=.1);

		# Mark Events
		#if(num_event_types>0){
		#	for(e_ix in 1:num_event_types){
		#		plot_event(offset_info[,"Offsets"], 
		#			subset_diversity+char_height*e_ix, events_mat[, e_ix]);
		#	}
		#}

		###############################################################

		# Plot abundances bar plots
		palette(category_colors);
		subset_norm=normalized_mat[subset_samples,];
		plot(offset_info[,"Offsets"], subset_diversity, main=paste(title, "Composition"),
			 xlab="Time", ylab="", type="n", col=i, lwd=2,
			 yaxt="n",
			 xlim=c(offset_ranges[1]-min_period/2, offset_ranges[2]+min_period/2), ylim=c(0,1+.25));

		abd_pos=c(0, .25, .33, .5, .66, .75, 1);
		axis(side=2, at=abd_pos, labels=sprintf("%2.2f", abd_pos), las=2, cex.axis=.75);
		abline(h=c(0, 1), lwd=.5, col="grey75");
		abline(h=c(.25, .5, .75), lwd=.4, col="grey85");
		abline(h=c(.33, .66), lwd=.4, col="grey95");

		char_height=par()$cxy[2]*.75;

		# Plot stacked abundances
		for(t in 1:num_members){
			plot_dist(offset_info[subset_samples[t],"Offsets"], y=0, 
				abundances=normalized_mat[subset_samples[t],], width=min_period*.95);
		}

		# Mark Events
		#if(num_event_types>0){
		#	for(e_ix in 1:num_event_types){
		#		plot_event(offset_info[subset_samples, "Offsets"], 
		#			rep(1.075+(e_ix-1)*char_height, num_members), events_mat[, e_ix]);
		#	}
		#}

		###############################################################

		# Plot color key/legend
		num_in_key=min(35, length(cat_names));
		plot(0, 0, main="",
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
	offsets_rec, grp_to_sbj_info, col_assign, ind_colors, grp_name=""){

	cat("Running plot_sample_diversity_by_group\n");

	offsets_mat=offsets_rec[["matrix"]];
	sorted_sids=offsets_rec[["SampleIDs"]];
	subject_ids=offsets_rec[["SubjectIDs"]];
	num_subjects=offsets_rec[["NumSubjects"]];

	# Get range of diversity
	diversity_ranges=range(diversity_arr);
	cat("Diversity Range:\n");
	print(diversity_ranges);

	offset_ranges=range(offsets_rec[["Offsets"]]);
	offset_span=diff(offset_ranges);
	x_plot_range=c(offset_ranges[1], offset_ranges[2]+offset_span*.15);

	# Set palette for individuals
	palette(ind_colors);

	num_groups=grp_to_sbj_info[["NumGroups"]];
	group_ids=grp_to_sbj_info[["Groups"]];

	# Set up plots per page
	def_par=par(no.readonly=T);
	par(mfrow=c(num_groups,1));

	for(g in 1:num_groups){

		cat("Plotting: ", group_ids[g], "\n");
		plot(0, 0, main=paste(grp_name, ": ", group_ids[g], sep=""),
			 xlab="Time", ylab=paste("Diversity (", div_type,")",sep=""), type="n", 
			 xlim=x_plot_range, ylim=diversity_ranges);

		# Get Unique Inidividuals in group
		indivs=grp_to_sbj_info[["GrpToSbj"]][[group_ids[g]]];
		num_indivs=length(indivs);
		cat("Number of Individuals: ", num_indivs, "\n");

		# Plot individual samples
		for(i in 1:num_indivs){

			# Grab from group
			cat("Plotting: ", indivs[i], "\n");

			offset_info=offsets_rec[["OffsetsBySubject"]][[indivs[i]]];
			num_timepts=nrow(offset_info);

			###############################################################

			# Subset diversity
			subset_samples=rownames(offset_info);
			subset_diversity=diversity_arr[subset_samples];

			# Plot Diversity
			points(offset_info[c(1,1, num_timepts),"Offsets"], subset_diversity[c(1,1, num_timepts)], 
				type="p", pch=c(17, 1, 15), cex=c(1, 2, 1.25), col=col_assign[indivs[i]]);
			points(offset_info[,"Offsets"], subset_diversity, type="l", 
				lwd=2.5, col=col_assign[indivs[i]]);
			points(offset_info[,"Offsets"], subset_diversity, type="l", lwd=.1, col="black");

			# Label on subject ID on right hand side
			text(offset_info[num_timepts, "Offsets"], subset_diversity[num_timepts],
				 indivs[i], pos=4);

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
	normalized_mat, offsets_rec, grp_to_sbj_info, col_assign, ind_colors, grp_name=""){

	cat("Running plot_median_sample_diversity_by_group\n");

	offsets_mat=offsets_rec[["matrix"]];
        sorted_sids=offsets_rec[["SampleIDs"]];
        subject_ids=offsets_rec[["SubjectIDs"]];
        num_subjects=offsets_rec[["NumSubjects"]];

	def_par=par(no.readonly=T);

	# Get Num Cohorts
	groups=grp_to_sbj_info[["Groups"]];
	num_groups=grp_to_sbj_info[["NumGroups"]];
	offset_pos=offsets_rec[["Offsets"]];
	offset_ranges=range(offset_pos);
	num_offset_positions=offsets_rec[["NumUniqOffsets"]];

	# Get range of diversity
	diversity_ranges=range(diversity_arr);
	cat("Diversity Range:\n");
	print(diversity_ranges);
	cat("\n\n");

	group_data=list();
	# Build data structure
	for(g in 1:num_groups){

		cat("Extracting: ", groups[g], "\n");

		# Get Unique Inidividuals
		ind_in_grp=grp_to_sbj_info[["GrpToSbj"]][[groups[g]]];
		num_indivs_in_grp=length(ind_in_grp);
		cat("Number of Subjects: ", num_indivs_in_grp, "\n");
		cat("Subjects:\n");
		print(ind_in_grp);
		cat("\n");

		group_id=as.character(groups[g]);
		group_data[[group_id]]=list();

		# Plot individual samples
		for(i in 1:num_indivs_in_grp){

			# Subset offsets, and sort by offset
			offset_info=offsets_rec[["OffsetsBySubject"]][[ind_in_grp[i]]];
			indiv_id=as.character(ind_in_grp[i]);
			group_data[[g]][[indiv_id]]=offset_info;
		}
	}

	#print(group_data);

	medians=list();
	loess=list();
	samps=list();
	
	# Compute Medians and compute loess per cohort
	for(g in 1:num_groups){

		group_id=as.character(groups[g]);
		indiv_ids=grp_to_sbj_info[["GrpToSbj"]][[group_id]];
		num_indv=length(indiv_ids);

		# per cohort medians across time
		med_mat=matrix(NA, nrow=num_offset_positions, ncol=3);
		colnames(med_mat)=c("Median", "LB_95CI", "UB_95CI");
		rownames(med_mat)=as.character(offset_pos);

		# per cohort div vs time
		div_offs_mat=matrix(NA, nrow=0, ncol=2);
		colnames(div_offs_mat)=c("Offsets", "Diversity");
		

		for(offset in offset_pos){
			div_arr=rep(NA,num_indv)
			
			for(i in 1:num_indv){

				idv=indiv_ids[i];
				ind_off=group_data[[group_id]][[idv]]		
				samp_names=rownames(ind_off);	
				targ_off_ix=which(ind_off[,"Offsets"]==offset);

				if(length(targ_off_ix)==1){
					div_arr[i]=diversity_arr[samp_names[targ_off_ix]];
					div_offs_mat=rbind(div_offs_mat, c(offset, div_arr[i]));
				}
			}	

			#cat("Offset: ", offset, "\n");

			bs_med=bootstrap_med(div_arr);
			qtl=quantile(bs_med, na.rm=T, probs=c(0.025, .5, 0.975));
			med_mat[as.character(offset), c("LB_95CI", "Median", "UB_95CI")]=qtl;
		}

		#print(div_offs_mat);
		samps[[group_id]]=div_offs_mat;

		# Fit loess	
		offset_tmp=div_offs_mat[,1];
		diversity_tmp=div_offs_mat[,2];

		loess_res=loess(diversity_tmp~offset_tmp);

		loess_fit_x=sort(unique(c(offset_pos, 
			seq(offset_ranges[1], offset_ranges[2], length.out=num_offset_positions*2))));
		loess_fit_y=tryCatch({
				predict(loess_res, loess_fit_x);
			}, error=function(e){
				return(NULL);
			});

		if(is.null(loess_fit_y)){
			# If loess not possible, just use original values.
			loess_fit_x=offset_tmp;
			loess_fit_y=diversity_tmp;
		}

		loess[[group_id]]=cbind(loess_fit_x, loess_fit_y);
		medians[[group_id]]=med_mat;

	}

	#print(medians);
	#print(loess);
	#print(samps);

	# Set up plots per page
	def_par=par(no.readonly=T);
	par(mfrow=c(3,1));

	colors_full  =rainbow(num_groups, start=0, end=4/6, s=1);
	colors_pastel=rainbow(num_groups, start=0, end=4/6, s=.3, alpha=.85);

	#----------------------------------------------------------------------
	# Plot Medians and PI's
	plot(0,0, type="n",
		xlim=c(offset_ranges[1], offset_ranges[2]),
		ylim=c(0, diversity_ranges[2]*1.2),
		xlab="Time",
		ylab="Median Diversity (w/ 95% CI)"
	);

	chwd=par()$cxy[1]/4;
	for(g in 1:num_groups){
		group_id=as.character(groups[g]);
		cur_med=medians[[group_id]]		

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
		#cat("[x,y]=[",x, ", ", y,"]\n");

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
	for(g in 1:num_groups){
		group_id=as.character(groups[g]);
		cur_samps=samps[[group_id]];

		x=cur_samps[,"Offsets"];
		y=cur_samps[,"Diversity"];
		points(x, y, type="p", cex=1, col=colors_pastel[g]);
	}
	
	# Draw loess line
	for(g in 1:num_groups){
		group_id=as.character(groups[g]);
		cur_loess=loess[[group_id]];

		x=cur_loess[,"loess_fit_x"];
		y=cur_loess[,"loess_fit_y"];
		points(x, y, type="l", col=colors_full[g], lwd=2);
	}

	#----------------------------------------------------------------------
	# Plot group legend
	plot(0,0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", bty="n", xaxt="n", yaxt="n");
	legend(0, 1, groups, fill=colors_full, bty="n", cex=2, title=grp_name);
	
	par(def_par);

}

###############################################################################

plot_sample_distributions_by_group=function(normalized_mat, offsets_rec, grp_to_sbj_info, cat_colors, grp_name=""){

	cat("Running plot_sample_distributions_by_group\n");

        offsets_mat=offsets_rec[["matrix"]];
        sorted_sids=offsets_rec[["SampleIDs"]];
        subject_ids=offsets_rec[["SubjectIDs"]];
        num_subjects=offsets_rec[["NumSubjects"]];

	groups=grp_to_sbj_info[["Groups"]];
	num_groups=grp_to_sbj_info[["NumGroups"]];
	offset_ranges=range(offsets_rec[["Offsets"]]);
	cat("Offset Range:\n");
	print(offset_ranges);

	min_period=offsets_rec[["MinOffsetSep"]];
	cat("Minimum period: ", min_period, "\n");

	palette(cat_colors);

	# Set up plots per page
	def_par=par(no.readonly=T);
	for(g in 1:num_groups){

		indivs=grp_to_sbj_info[["GrpToSbj"]][[groups[g]]];
		num_indivs=length(indivs);

		par(mfrow=c(num_indivs,1));
		par(oma=c(3,1,2,0));
		par(mar=c(0.1, 4.1, 0.1, 2.1)); # bot, left, top, right

		# Plot individual samples
		for(i in 1:num_indivs){

			# Grab from group
			cat("Plotting: ", as.character(indivs[i]), "\n");

			ind_offset=offsets_rec[["OffsetsBySubject"]][[indivs[i]]];
			num_timepts=nrow(ind_offset);
			cat("Num Time Pts: ", num_timepts, "\n");

			###############################################################
			if(i==num_indivs){
				xaxt_setting="s";	
			}else{
				xaxt_setting="n";
			}
		
			num_event_types=0;

			plot(0, 0, main="",
				xlab="", ylab=indivs[i], type="n", bty="n",
				xaxt=xaxt_setting, yaxt="n",
				xlim=c(offset_ranges[1]-min_period/2, offset_ranges[2]+min_period/2), 
				ylim=c(0, 1+.5*num_event_types));

			char_height=par()$cxy[2]*.5;

			subset_samples=rownames(ind_offset);

			# Draw guide lines
			abline(h=.5, col="grey", lwd=20);
			xaxp=par()$xaxp;
			xguides=seq(xaxp[1], xaxp[2], (xaxp[2]-xaxp[1])/xaxp[3]);
			#print(xguides);
			abline(v=xguides, col="grey", lwd=.5);

			# Draw stacked barplot
			for(t in 1:num_timepts){
				plot_dist(ind_offset[t,"Offsets"], y=0, 
					abundances=normalized_mat[subset_samples[t],], width=min_period*.95);

			}

			# Mark Events
			if(num_event_types>0){
				for(e_ix in 1:num_event_types){
					plot_event(ind_offset[subset_samples, "Offsets"], 
						rep(1.15+(e_ix-1)*char_height, num_timepts), events_mat[, e_ix]);
				}
			}
		}

		mtext(paste(grp_name, ": ", groups[g], sep=""), side=3, outer=T, font=2);
	}
	par(def_par);
}

###############################################################################

plot_mean_distributions_by_group=function(normalized_mat, offsets_rec, grp_to_sbj_info, cat_colors, grp_name=""){

	cat("Running plot_sample_distributions_by_group\n");

	num_categories=ncol(normalized_mat);

        offsets_mat=offsets_rec[["matrix"]];
        sorted_sids=offsets_rec[["SampleIDs"]];
        subject_ids=offsets_rec[["SubjectIDs"]];
        num_subjects=offsets_rec[["NumSubjects"]];

        groups=grp_to_sbj_info[["Groups"]];
        num_groups=grp_to_sbj_info[["NumGroups"]];
        offset_ranges=range(offsets_rec[["Offsets"]]);
	num_uniq_offsets=offsets_rec[["NumUniqOffsets"]];
	uniq_offsets=offsets_rec[["Offsets"]];
	
        cat("Offset Range:\n");
        print(offset_ranges);

        min_period=offsets_rec[["MinOffsetSep"]];
        cat("Minimum period: ", min_period, "\n");

	palette(cat_colors);

	# Set up plots per page
	def_par=par(no.readonly=T);
	par(oma=c(0,0,2,0));
	par(mfrow=c(num_groups+1, 1));

	for(g in 1:num_groups){

		# Save mean abundances in matrix
		avg_off_comp=matrix(0, nrow=num_uniq_offsets, ncol=num_categories);
		rownames(avg_off_comp)=as.numeric(uniq_offsets);
		colnames(avg_off_comp)=colnames(normalized_mat);

		# Save num samples per offset in array
		samp_size_at_offset=numeric(num_uniq_offsets);
		names(samp_size_at_offset)=uniq_offsets;

		subj_in_grp=grp_to_sbj_info[["GrpToSbj"]][[groups[g]]];

		# Calculate means and sample size per offset
		for(cur_off in uniq_offsets){

			samp_id_list=c();

			for(sbj in subj_in_grp){

				sbj_offset_info=offsets_rec[["OffsetsBySubject"]][[sbj]];
				off_ix=(cur_off==sbj_offset_info[,"Offsets"]);
				if(all(!off_ix)){
					next;
				}else{
					sbj_samp_id=rownames(sbj_offset_info[off_ix,]);
					samp_id_list=c(samp_id_list, sbj_samp_id);
				}
			}

			offset_key=as.character(cur_off);
			avg_off_comp[offset_key, ]=apply(normalized_mat[samp_id_list,,drop=F], 2, mean);
			samp_size_at_offset[offset_key]=length(samp_id_list);
			
		}

		#--------------------------------------------------------------
		# Plot times series for cohort
		plot(0, 0, main="",
			xlab="", 
			ylab=paste(grp_name, ": ", groups[g]), 
			type="n", bty="n",
			xaxt="s", yaxt="n",
			xlim=c(offset_ranges[1]-min_period/2, offset_ranges[2]+min_period/2), ylim=c(0, 1.1));

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
				abundances=avg_off_comp[offset_key,,drop=F], width=min_period*.95);

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

plot_change_scatter=function(diversity_arr, offsets_rec, grp_to_sbj_info, grp_name=""){

	cat("Plotting Change Scatter:\n");

        offsets_mat=offsets_rec[["matrix"]];
        sorted_sids=offsets_rec[["SampleIDs"]];
        subject_ids=offsets_rec[["SubjectIDs"]];
        num_subjects=offsets_rec[["NumSubjects"]];

        groups=grp_to_sbj_info[["Groups"]];
        num_groups=grp_to_sbj_info[["NumGroups"]];
        offset_ranges=range(offsets_rec[["Offsets"]]);
	num_uniq_offsets=offsets_rec[["NumUniqOffsets"]];
	uniq_offsets=offsets_rec[["Offsets"]];

	# Extract start/end diversity
	ends=matrix(NA, nrow=num_subjects, ncol=2, dimnames=list(subject_ids, c("start", "end")));
	for(sbj_ix in subject_ids){

		sbj_offset=offset_rec[["OffsetsBySubject"]][[sbj_ix]];
		num_offsets=nrow(sbj_offset);
		samp_ids=rownames(sbj_offset);

		ends[sbj_ix, "start"]=diversity_arr[samp_ids[1]];
		ends[sbj_ix, "end"]=diversity_arr[samp_ids[num_offsets]];
	}

	end_ranges=range(ends);	
	axis_ticks=round(seq(end_ranges[1], end_ranges[2], length.out=10), 3);

	grp_col_transp=get_colors(num_groups, alpha=.2);
	grp_col_opaque=get_colors(num_groups, alpha=1);


	par(mar=c(5,6,3,2));
	par(mfrow=c(num_groups+1, 1));

	for(grp_ix in 1:num_groups){

		cur_grp=groups[grp_ix];

		grp_subject_ids=grp_to_sbj_info[["GrpToSbj"]][[cur_grp]];
		grp_ends=ends[grp_subject_ids,, drop=F];

		plot(0,0, type="n",
			main=paste(grp_name, ": ", cur_grp, sep=""),
			xlab="", ylab="",
			xlim=end_ranges, ylim=end_ranges, 
			xaxt="n", yaxt="n");
		axis(side=1, at=axis_ticks, label=axis_ticks);
		axis(side=2, at=axis_ticks, label=axis_ticks, las=2);
		abline(a=0, b=1, col="grey", lwd=2);
		title(xlab="Start", line=3);
		title(ylab="End", line=4);

		median_end=median(grp_ends[,"end"]);
		median_start=median(grp_ends[,"start"]);

		abline(h=median_end, col=grp_col_opaque[grp_ix]);
		abline(v=median_start, col=grp_col_opaque[grp_ix]);

		axis(side=3, at=median_start, label=round(median_start, 3), 
			line=-1, tick=F, cex.axis=.75, col.axis=grp_col_opaque[grp_ix]);
		axis(side=4, at=median_end, label=round(median_end, 3), 
			line=-1, tick=F, cex.axis=.75, col.axis=grp_col_opaque[grp_ix]);

		points(grp_ends, col=grp_col_transp[grp_ix], cex=2, pch=19);
		points(grp_ends, col=grp_col_opaque[grp_ix], cex=.25, pch=19);
		
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
	for(grp_ix in 1:num_groups){

                grp=groups[grp_ix];
		grp_subject_ids=grp_to_sbj_info[["GrpToSbj"]][[grp]];

                grp_ends=ends[grp_subject_ids,, drop=F];

		all_ends_and_groups[[paste(grp, ":start", sep="")]]=grp_ends[,"start"];
		all_ends_and_groups[[paste(grp, ":end", sep="")]]=grp_ends[,"end"];

		median_end=median(grp_ends[,"end"]);
		median_start=median(grp_ends[,"start"]);

		abline(h=median_end, col=grp_col_opaque[grp_ix]);
		abline(v=median_start, col=grp_col_opaque[grp_ix]);
		axis(side=3, at=median_start, label=round(median_start, 3), line=-1, 
			tick=F, cex.axis=.75, col.axis=grp_col_opaque[grp_ix]);
		axis(side=4, at=median_end, label=round(median_end, 3), line=-1, 
			tick=F, cex.axis=.75, col.axis=grp_col_opaque[grp_ix]);

		points(grp_ends, col=grp_col_transp[grp_ix], cex=2, pch=19);
		points(grp_ends, col=grp_col_opaque[grp_ix], cex=.25, pch=19);
	}

	#######################################################################
	# Compute pvalues between groups

	aeag_names=names(all_ends_and_groups);

	for(aeag_idx in aeag_names){
		cat(aeag_idx, ": ", median(all_ends_and_groups[[aeag_idx]]), "\n", sep="");	
	}

	pval_matrix=matrix(NA, nrow=num_groups*2, ncol=num_groups*2, dimnames=list(aeag_names, aeag_names));
	diff_matrix=matrix(NA, nrow=num_groups*2, ncol=num_groups*2, dimnames=list(aeag_names, aeag_names));
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

factor_info=load_factors(FactorFile);

###############################################################################

counts_mat=load_summary_file(SummaryFile);
#print(counts_mat);

###############################################################################

factors_mat_samples=rownames(factor_info);
counts_mat_samples=rownames(counts_mat);
shared=intersect(factors_mat_samples, counts_mat_samples);

cat("\n\n");
cat("Samples not represented in factor file:\n");
print(setdiff(counts_mat_samples, shared));
cat("Samples not represented in summary table:\n");
print(setdiff(factors_mat_samples, shared));
cat("\n\n");

factor_info=factor_info[shared,,drop=F];
counts_mat=counts_mat[shared,,drop=F];

if(GroupCol==""){
        group_names=rep("All", nrow(factor_info));
        GroupCol="Group";
}else{
        group_names=factor_info[,GroupCol];
}
grp_to_sbj_info=create_GrpToSbj_map(factor_info[,SubjectIDCol], group_names);
unique_group_names=grp_to_sbj_info[["Groups"]];

# Extract offset info
offset_rec=extract_offset(factor_info, SubjectIDCol, TimeOffsetCol);
print(offset_rec);

print(grp_to_sbj_info);

indiv_ids=offset_rec[["SubjectIDs"]];
num_indiv=length(indiv_ids);

###############################################################################

plot_text(c(
	"Subject-Group Mapping:",
	"",
	capture.output(print(grp_to_sbj_info))
));

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

plot_text(c(paste("By Subject (", SubjectIDCol, "):  Diversity and Composition", sep="")));
plot_sample_distributions_by_individual(diversity_arr, DiversityType, 
	simplified_mat, offset_rec, grp_to_sbj_info, col_assign, category_colors, ind_colors);

plot_text(c(paste("By Group (", GroupCol, "):  Diversity", sep="")));
plot_sample_diversity_by_group(diversity_arr, DiversityType, simplified_mat, 
	offset_rec, grp_to_sbj_info, col_assign, ind_colors, grp_name=GroupCol);

plot_text(c(paste("By Group (", GroupCol, "):  Median Sample Diversity", sep="")));
plot_median_sample_diversity_by_group(diversity_arr, DiversityType, simplified_mat, 
	offset_rec, grp_to_sbj_info, col_assign, ind_colors, grp_name=GroupCol);

plot_text(c(paste("By Group (", GroupCol, "):  Composition as Stacked Bar Plots", sep="")));
plot_sample_distributions_by_group(simplified_mat, offset_rec, grp_to_sbj_info, category_colors,
	grp_name=GroupCol);

plot_text(c(paste("By Group (", GroupCol, "):  Mean Composition as Stacked Bar Plots", sep="")));
plot_mean_distributions_by_group(simplified_mat, offset_rec, grp_to_sbj_info, category_colors,
	grp_name=GroupCol);

plot_text(c(paste("By Group (", GroupCol, "):  Delta Scatter Plots", sep="")));
plot_change_scatter(diversity_arr, offset_rec, grp_to_sbj_info,
	grp_name=GroupCol);

##############################################################################

#print(diversity_arr);
#print(offset_rec);

diversity_mat=matrix(diversity_arr, nrow=length(diversity_arr), ncol=1);
colnames(diversity_mat)=c("diversity");
rownames(diversity_mat)=names(diversity_arr);

long_stats=calc_longitudinal_stats(offset_rec, diversity_mat);

#print(long_stats);

stat_names=names(long_stats);
for(stat_ix in stat_names){
        mat=long_stats[[stat_ix]];
        paint_matrix(mat, stat_ix);
}

#print(grp_to_sbj_info);

stat_table_grp_cmp=plot_pairwise_grp_comparisons(long_stats, grp_to_sbj_info, plots_pp=1);

#print(stat_table_grp_cmp);

output_stat_table_alternate_ordering(stat_table_grp_cmp, OutputRoot);

##############################################################################

cat("Done.\n")
dev.off();
warn=warnings();
if(length(warn)){
	print(warn);
}
q(status=0)
