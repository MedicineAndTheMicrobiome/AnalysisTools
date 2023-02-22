#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

options(useFancyQuotes=F);
options(digits=5)

params=c(
	"factors", "f", 1, "character",
	"targets", "t", 2, "character",
	"outputroot", "o", 1, "character",
	"maxcluster", "m", 2, "numeric",
	"subsample", "s", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <factors/metadata file name>\n",
	"	-t <target file name, for variables to compute clusters on>\n",
	"	-o <output filename root>\n",
	"	[-m <maximum number of clusters, default=max(log2(sample_size), sample_size/40)>]\n",
	"	[-s <subsample for testing, default=all samples>]\n",
	"\n",
	"This script will use the variables specified in the target file\n",
	"to assign subjects to clusters.  The cluster assignments will be\n",
	"output to a the a new factors/metadata file.\n",
	"\n",
	"Depending on your variables, you may want to first run them\n",
	"through some normality transformation or variable selection process\n",
	"ahead of time.\n",
	"\n");

if(
	!length(opt$factors) || 
	!length(opt$targets) || 
	!length(opt$outputroot)
){
	cat(usage);
	q(status=-1);
}

FactorsFname=opt$factors;
TargetsFname=opt$targets;
OutputFnameRoot=opt$outputroot;

MaxClusters=-1;
if(length(opt$maxcluster)){
	MaxClusters=opt$maxcluster;
}

Subsample=-1;
if(length(opt$subsample)){
	Subsample=opt$subsample;
}


param_text=capture.output({
	cat("\n");
	cat("Factor File Name: ", FactorsFname, "\n");
	cat("Targets File Name: ", TargetsFname, "\n");
	cat("Output File Name Root: ", OutputFnameRoot, "\n");
	cat("Maximum Cluster Cuts: ", MaxClusters, "\n");
	cat("Subsample:", Subsample, "\n");
});
print(param_text, quote=F);
cat("\n\n");

###############################################################################

load_factors=function(fname){
	factors=as.data.frame(read.delim(fname,  header=TRUE, row.names=1, check.names=FALSE, sep="\t"));
	return(factors);
}

load_list=function(fname){
	cat("Loading: ", fname, "\n");
	lst=read.delim(fname, header=F, check.names=F, comment.char="#", as.is=T);
	return(lst[,1]);	
}

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

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, deci_pts=4,
        label_zeros=T, counts=F, value.cex=1,
        plot_col_dendr=F,
        plot_row_dendr=F
){

        num_row=nrow(mat);
        num_col=ncol(mat);

        row_names=rownames(mat);
        col_names=colnames(mat);

        orig.par=par(no.readonly=T);

        cat("Num Rows: ", num_row, "\n");
        cat("Num Cols: ", num_col, "\n");

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
                plot_min=min(mat, na.rm=T);
        }
        if(is.na(plot_max)){
                plot_max=max(mat, na.rm=T);
        }

        if(plot_min>=-1 && plot_max<=1){
                fractions_only=T;
        }else{
                fractions_only=F;
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

        if(plot_col_dendr && plot_row_dendr){
                layoutmat=matrix(
                        c(
                        rep(c(rep(4, row_dend_width), rep(3, heatmap_width)), col_dend_height),
                        rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height)
                        ), byrow=T, ncol=row_dend_width+heatmap_width);

                col_dendr=get_dendrogram(mat, type="col");
                row_dendr=get_dendrogram(mat, type="row");

                mat=mat[row_dendr[["names"]], col_dendr[["names"]]];

        }else if(plot_col_dendr){
                layoutmat=matrix(
                        c(
                        rep(rep(2, heatmap_width), col_dend_height),
                        rep(rep(1, heatmap_width), heatmap_height)
                        ), byrow=T, ncol=heatmap_width);

                col_dendr=get_dendrogram(mat, type="col");
                mat=mat[, col_dendr[["names"]]];

        }else if(plot_row_dendr){
                layoutmat=matrix(
                        rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height),
                        byrow=T, ncol=row_dend_width+heatmap_width);

                row_dendr=get_dendrogram(mat, type="row");
                mat=mat[row_dendr[["names"]],];
        }else{

		if(heatmap_height < 200 &&  heatmap_width < 200){
			layoutmat=matrix(
				rep(1, heatmap_height*heatmap_width),
				byrow=T, ncol=heatmap_width);
		}else{
			larger_dim=max(heatmap_height, heatmap_width);

			hm_height_norm=heatmap_height/larger_dim;	
			hm_width_norm=heatmap_width/larger_dim;

			lo_height=100*hm_height_norm;
			lo_width=100*hm_width_norm;

			layoutmat=matrix(
				rep(1, lo_height*lo_width),
				byrow=T, ncol=lo_width);

		}
			
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

                        if(is.na(mat[y,x]) || mat[y,x]!=0 || label_zeros){
                                if(counts){
                                        text_lab=sprintf("%i", mat[y,x]);
                                }else{
                                        text_lab=sprintf(paste("%0.", deci_pts, "f", sep=""), mat[y,x]);
                                        if(fractions_only){
                                                if(!is.na(mat[y,x])){
                                                        if(mat[y,x]==-1 || mat[y,x]==1){
                                                                text_lab=as.integer(mat[y,x]);
                                                        }else{
                                                                text_lab=gsub("0\\.","\\.", text_lab);
                                                        }
                                                }
                                        }
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

plot_title_page=function(title, subtitle="", title_cex=3){

        orig.par=par(no.readonly=T);
        par(family="serif");
        par(mfrow=c(1,1));

        plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n",  xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                );

        # Title
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
# Main Program Starts Here!
##############################################################################

pdf(paste(OutputFnameRoot, ".cluster.pdf", sep=""), height=8.5, width=11);

plot_text(param_text);

# Load factors
cat("Loading Factors...\n");
loaded_factors=load_factors(FactorsFname);

if(Subsample!=-1){
	loaded_factors=loaded_factors[sample(nrow(loaded_factors), Subsample),];
}

loaded_factor_names=colnames(loaded_factors);
loaded_subject_names=rownames(loaded_factors);
num_subjects=length(loaded_subject_names);

cat("Loaded factors:\n");
print(loaded_factor_names);
cat("\n");

cat("Loaded subjects ids:\n");
print(head(loaded_subject_names, 10));
if(num_subjects>10){
	cat("...\n");
}
cat("Number of Subjects: ", num_subjects, "\n");
cat("\n");

# Load Targets
target_arr=load_list(TargetsFname);
num_targets=length(target_arr);
cat("Targets to include for distance measure:\n");
print(target_arr);

missing_targets=setdiff(target_arr, loaded_factor_names);
if(length(missing_targets)){
	cat("WARNING:  Missing targets:\n");
	print(missing_targets);
}else{
	cat("All targets found...\n");
}
cat("\n\n");

targeted_factors=loaded_factors[,target_arr];

is_numeric=apply(targeted_factors, 1, function(x){all(is.finite(x))});
if(any(!is_numeric)){
	cat("WARNING: Subjects removed because of non-finite values (e.g. NAs, NaN, etc)\n");
	print(targeted_factors[!is_numeric,]);
	cat("\n\n");
}
targeted_factors=targeted_factors[is_numeric,];
num_subjects=nrow(targeted_factors);

# Remove 0 information variables
targets_wstdev=apply(targeted_factors, 2, function(x){sd(x)>0});
targeted_factors=targeted_factors[,targets_wstdev,drop=F];
target_arr=colnames(targeted_factors);
num_targets=length(target_arr);

###############################################################################

cat("Standardizing...\n");

plot_title_page("Standardizing Targeted Variables", c(
	"Variables were first standardized (0 centered/mean, 1 stdev) so the distance",
	"measure that combines all the variables, will not be more greatly",
	"influenced by some variables that have nominally greater variances."
));

stats=matrix(0, nrow=num_targets, ncol=2);
colnames(stats)=c("Means", "StDev");
rownames(stats)=target_arr;

means=apply(targeted_factors, 2, mean);
stdev=apply(targeted_factors, 2, sd);
stats[,"Means"]=means;
stats[,"StDev"]=stdev;

print(stats);
plot_text(c(
	"Means and Standard Deviations before Standardization:",
	"",
	capture.output(print(stats))
));
cat("\n");

standardized_targets=matrix(NA, nrow=num_subjects, ncol=num_targets);
colnames(standardized_targets)=colnames(targeted_factors);
rownames(standardized_targets)=rownames(targeted_factors);
for(tar_ix in 1:num_targets){
	standardized_targets[,tar_ix]=(targeted_factors[,tar_ix]-means[tar_ix])/stdev[tar_ix];
}

#print(standardized_targets);

###############################################################################

#standardized_targets=standardized_targets[sample(200),];

calc_antioutlier_weight=function(std_val_mat){
	num_var=ncol(std_val_mat);
	num_sbj=nrow(std_val_mat);
	var_names=colnames(std_val_mat);

	weights=numeric(num_var);
	names(weights)=var_names;

	means=apply(std_val_mat, 2, mean);

	for(i in 1:num_var){
		prop_gt_mean=sum(std_val_mat[,i] > means[i])/num_sbj;
		spread=abs(1-2*prop_gt_mean);
		weight=1-spread;
		if(is.na(weight)){
			weight=0;
		}
		cat(var_names[i], " Prop: ", prop_gt_mean, " Spread: ", spread, " Weight: ", weight, "\n", sep="");
		weights[i]=weight;
	}
	weights=weights/(sum(weights));

	return(weights);
}

UseAntiOutlierWeighting=T;
if(UseAntiOutlierWeighting){

	plot_title_page("Anti-Outlier Weighting", c(
		"Anti-outlier weighting calculates a weighting for each variable so that",
		"variables appearing to represent rare observations (e.g. rare diseases)",
		"are weighted less than those that are more ubiquitous.",
		"",
		"To calculate the weighting:",
		"   1.) calculate the mean",
		"   2.) calculate the PGT: proportion of observations > mean",
		"   3.) calculate the spread = | 1 - 2*PGT |",
		"   4.) Weight = 1-spread",
		"",
		"The weights are then normalized across the variables to provide",
		"context for how much each variables is contributing to the distance.",
		"",
		"For a boolean variable, the anti-outlier weighting will give more",
		"weight to variables that are closer to 50%/50% distributed, and",
		"much less weight to those that are closer to 99%/1% distributed."
		));

	cat("Calculating anti-outlier weights...\n");
	anti_outlier_weights=calc_antioutlier_weight(standardized_targets);
	print(anti_outlier_weights);

	# Reorder variables by their weights
	weight_order=order(anti_outlier_weights);
	anti_outlier_weights=anti_outlier_weights[weight_order];
	standardized_targets=standardized_targets[,weight_order];
	targeted_factors=targeted_factors[,weight_order];

	plot_text(c(
		"Anti-outlier Weights (ordered):",
		"",
		capture.output(print(anti_outlier_weights))
	));

	decreasing_aow=sort(anti_outlier_weights, decreasing=T);
	par(mar=c(10,4,4,1));
	barplot(decreasing_aow, xlab="", ylab="Weighting", las=2,
		main="Anti-Outlier Weighting");


}


weighted_dist=function(variables, weight=numeric()){
	num_sbj=nrow(variables);
	num_var=ncol(variables);

	dist_mat=matrix(NA, ncol=num_sbj, nrow=num_sbj);
	colnames(dist_mat)=rownames(variables);
	rownames(dist_mat)=rownames(variables);

	if(length(weight)==0){
		weight=rep(1, num_var);
	}

	for(i in 1:num_sbj){
		for(j in 1:i){
			I=variables[i,];
			J=variables[j,];
			dist_mat[i,j]=sqrt(sum(weight*(I-J)^2));
			dist_mat[j,i]=dist_mat[i,j];
		}

	}
	return(dist_mat);

}


#distance_mat=weighted_dist(standardized_targets);
cat("Calculating intersubject distances...\n");
print(standardized_targets);
distance_mat=weighted_dist(standardized_targets, anti_outlier_weights);

###############################################################################

distance_subs=as.dist(distance_mat);
#print(distance_subs);

cat("Hierarchically clustering...\n");
hcl=hclust(distance_subs, method="ward.D2");
dend=as.dendrogram(hcl);

resize_labels=function(in_dend, in.lab.cex){

	dend_resize_labels=function(x){
		if(is.leaf(x)){
			leaf_attr=attributes(x);
			label=leaf_attr$label;
			#print(label);
			attr(x, "nodePar")=c(leaf_attr$nodePar, list(lab.cex=in.lab.cex, cex=in.lab.cex/2));
		}
		return(x);
	}

	out_dend=dendrapply(in_dend, dend_resize_labels);
	return(out_dend);	
}

color_labels=function(in_dend, colmap){

	dend_resize_labels=function(x){
		if(is.leaf(x)){
			leaf_attr=attributes(x);
			label=leaf_attr$label;
			#print(label);
			leafcol=colmap[label];
			attr(x, "nodePar")=c(leaf_attr$nodePar, list(lab.col=leafcol, col=leafcol));
		}
		return(x);
	}

	out_dend=dendrapply(in_dend, dend_resize_labels);
	return(out_dend);	
}

###############################################################################

find_height_at_k=function(hclust, k){
# Computes the height on the dendrogram for a particular k

        heights=hclust$height;
        num_heights=length(heights);
        num_clust=numeric(num_heights);
        for(i in 1:num_heights){
                num_clust[i]=length(unique(cutree(hclust, h=heights[i])));
        }
        height_idx=which(num_clust==k);
        midpoint=(heights[height_idx+1]+heights[height_idx])/2;
        return(midpoint);
}

get_clstrd_leaf_names=function(den){
# Recursively gets a list of the leaf names, from left to right from
#   the specified dendrogram

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

get_middle_of_groups=function(clustered_leaf_names, group_asgn){
# Finds middle of each group in the plot
        num_leaves=length(group_asgn);
        groups=sort(unique(group_asgn));
        num_groups=length(groups);

        reord_grps=numeric(num_leaves);
        names(reord_grps)=clustered_leaf_names;
        reord_grps[clustered_leaf_names]=group_asgn[clustered_leaf_names];

        mids=numeric(num_groups);
        names(mids)=1:num_groups;
        for(i in 1:num_groups){
                mids[i]=mean(range(which(reord_grps==i)));
        }
        return(mids);

}

###############################################################################

plot_annotated_dendrogram=function(in_dend, cluster_map, cut, cl_mids, h){

	num_samples=attr(dend, "members");

	# Calculate label size
	lab_size=min(1, 40/num_samples);
	cat("Label Size:", lab_size, "\n");

	tmp_dend=resize_labels(in_dend, in.lab.cex=lab_size);
	tmp_dend=color_labels(tmp_dend, colmap=cluster_map);

	plot(tmp_dend, main=paste("k = ", cut, sep=""));

	#text(cl_mids, rep(0, cut),  1:cut, pos=3, col=1:cut, font=2, cex=2);

	abline(h=h, col="blue", lty="dashed");

}

plot_variables_heatmap=function(factors, subj_order, cl_mids, num_colors=64){

	factors=factors[subj_order,];

	num_subjects=length(subj_order);
	num_variables=ncol(factors);
	variable_names=colnames(factors);

	plot(0,0, type="n", 
		xlim=c(0, num_subjects),
		ylim=c(0, num_variables),
		xaxt="n", yaxt="n", xlab="", ylab="", main="", bty="n");

	heat_col=rev(rainbow(num_colors, end=2/3));

	map_val_to_col=function(num_levels, values){
		# Remaps values to between 1 and num_levels
		minmax=range(values);
		span=minmax[2]-minmax[1];
		values=(values-minmax[1])/span;	
		values=ceiling(values*(num_levels-1))+1;
		return(values);
	}

	for(var_ix in 1:num_variables){

		mtext(variable_names[var_ix], at=var_ix-.5, side=4, las=2, line=-2, adj=0);

		values=factors[,var_ix];
		mapped_col=map_val_to_col(num_colors, values);

		for(sbj_ix in 1:num_subjects){

			rect(sbj_ix-1, var_ix-1, sbj_ix, var_ix,
				border=NA, 
				col=heat_col[mapped_col[sbj_ix]]);

		}
	}	

	cuts=length(cl_mids);
	mtext(1:cuts, side=3, line=-1, at=cl_mids, font=2, 
		cex=ifelse((1:cuts)<10, 2, 1.5),
		col=1:cuts);
	

	abline(v=c(0,num_subjects));
}

find_cluster_features=function(num_clusters, clst_map, factors_mat, pvalcutoff=0.1, bonf_adj=T){

	num_vars=ncol(factors_mat);
	sbj_ids=names(clst_map);

	features_rec=list();
	varnames=colnames(factors_mat);

	if(bonf_adj){
		pvalcutoff=pvalcutoff/(num_vars*num_clusters);
	}

	for(k in 1:num_clusters){
		members=sbj_ids[clst_map==k];
		non_members=sbj_ids[clst_map!=k];

		feat_tab=matrix(0, ncol=5, nrow=num_vars);		
		colnames(feat_tab)=c("isHigh", "p-value", "member_mean", "other_mean", "signif");
		rownames(feat_tab)=varnames;

		for(var in 1:num_vars){

			member_val=factors_mat[members,var];
			non_member_val=factors_mat[non_members,var];

			member_mean=mean(member_val);
			other_mean=mean(non_member_val);

			wc.res=wilcox.test(member_val, non_member_val);

			# 2 tail
			signif=ifelse(wc.res$p.value<(pvalcutoff/2), T, F);

			if(member_mean>other_mean){
				isHigh=T;
			}else{
				isHigh=F;
			}

			feat_tab[var,]=c(isHigh, wc.res$p.value, member_mean, other_mean, signif);
		}		

		order_ix=order(feat_tab[,"p-value"]);

		features_rec[[k]]=feat_tab[order_ix,];
	}

	print(features_rec);	
	#quit();
	return(features_rec);
	

}

plot_variables_heatmap_avgs=function(factors, grp_split_loc, cl_feat_rec, num_colors=64){

	num_variables=ncol(factors);
	num_subjects=nrow(factors);
	num_clusters=length(cl_feat_rec);

	# Use the order here, not the variable order in the cluster features rec
	variable_names=colnames(factors);

	plot(0,0, type="n", 
		xlim=c(0, num_subjects),
		ylim=c(0, num_variables),
		xaxt="n", yaxt="n", xlab="", ylab="", main="", bty="n");

	heat_col=rev(rainbow(num_colors, end=2/3));

	map_val_to_col=function(num_levels, ref_values, tar_val){
		# Remaps values to between 1 and num_levels
		minmax=range(ref_values);
		span=minmax[2]-minmax[1];
		tar_val=(tar_val-minmax[1])/span;	
		level=ceiling(tar_val*(num_levels-1))+1;
		return(level);
	}

	cl_loc=c(0, grp_split_loc, num_subjects);

	for(var_ix in 1:num_variables){

		cur_varname=variable_names[var_ix];
		mtext(cur_varname, at=var_ix-.5, side=4, las=2, line=-2, adj=0);

		for(cl_ix in 1:num_clusters){

			# Get the variable information for cluster
			cluster_feat_mat=cl_feat_rec[[cl_ix]];
			var_cl_mean=cluster_feat_mat[cur_varname, "member_mean"];
			var_cl_signf=cluster_feat_mat[cur_varname, "signif"];
			var_cl_ishigh=cluster_feat_mat[cur_varname, "isHigh"];

			# Map the mean to color
			mapped_col=map_val_to_col(num_colors, factors[, cur_varname], 
				var_cl_mean);

			# Get location of start/end for each cluster
			cl_start=cl_loc[cl_ix];
			cl_end=cl_loc[cl_ix+1];

			# Draw rectangle for the cluster/variable
			rect(cl_start, var_ix-1, cl_end, var_ix,
				border=NA, 
				col=heat_col[mapped_col]);

			# Draw annotation if significant
			if(var_cl_signf){
				ch=ifelse(var_cl_ishigh, 24, 25);
				points((cl_start+cl_end)/2, var_ix-.5, pch=ch, 
					bg=ifelse(var_cl_ishigh, heat_col[num_colors], 
					heat_col[1]));
			}

		}
	}	

	cuts=length(cl_mids);
	mtext(1:cuts, side=3, line=-1, at=cl_mids, font=2, 
		cex=ifelse((1:cuts)<10, 2, 1.5),
		col=1:cuts);
	

	abline(v=c(0,num_subjects));
}

plot_variables_lists=function(factors, grp_split_loc, cl_feat_rec, num_colors=64){

	num_variables=ncol(factors);
	num_subjects=nrow(factors);
	num_clusters=length(cl_feat_rec);

	plot(0,0, type="n", 
		xlim=c(0, num_subjects),
		ylim=c(0, num_variables+2),
		xaxt="n", yaxt="n", xlab="", ylab="", main="", bty="n");

	heat_col=rev(rainbow(num_colors, end=2/3));

	map_val_to_col=function(num_levels, ref_values, tar_val){
		# Remaps values to between 1 and num_levels
		minmax=range(ref_values);
		span=minmax[2]-minmax[1];
		tar_val=(tar_val-minmax[1])/span;	
		level=ceiling(tar_val*(num_levels-1))+1;
		return(level);
	}

	# Location of cluster widths
	cl_loc=c(0, grp_split_loc, num_subjects);
	# Location of colums widths
	column_pos=seq(0, num_subjects, length.out=num_clusters+1);
	# Location of column mids
	column_mids=(column_pos+(column_pos[2]-column_pos[1])/2)[1:num_clusters];
	# cluster sizes
	cluster_sizes=diff(cl_loc);


	for(cl_ix in 1:num_clusters){

		cat("\nCluster: ", cl_ix, "\n");
		# Get the variable information for cluster
		cluster_feat_mat=cl_feat_rec[[cl_ix]];
		print(cluster_feat_mat);

		los_ix=cluster_feat_mat[,"isHigh"]==0 & cluster_feat_mat[,"signif"]==1;
		his_ix=cluster_feat_mat[,"isHigh"]==1 & cluster_feat_mat[,"signif"]==1;

		lo_tab=cluster_feat_mat[los_ix,,drop=F];
		hi_tab=cluster_feat_mat[his_ix,,drop=F]

		lo_pval_ix=order(lo_tab[,"p-value"]);
		hi_pval_ix=order(hi_tab[,"p-value"]);

		lo_tab=lo_tab[lo_pval_ix,,drop=F];
		hi_tab=hi_tab[hi_pval_ix,,drop=F];

		#cat("\nLows:\n");
		#print(lo_tab);
		#cat("\nHighs:\n");
		#print(hi_tab);

		lo_varnames=rownames(lo_tab);
		hi_varnames=rownames(hi_tab);

		num_lo=nrow(lo_tab);
		num_hi=nrow(hi_tab);

		# Label cell with variable name
		if(num_lo>0){
			for(i in 1:num_lo){
				# Map the mean to color
				vname=lo_varnames[i];
				mapped_col=map_val_to_col(num_colors, 
					factors[, vname], lo_tab[vname, "member_mean"]);

				# Draw rectangle for the cluster/variable
				rect(column_pos[cl_ix], i, 
					column_pos[cl_ix+1], i+1,
					border=NA, 
					col=heat_col[mapped_col]);
				text(column_mids[cl_ix], i+.5, vname, adj=.5, font=2);
			}
		}

		if(num_hi>0){
			for(i in 1:num_hi){
				# Map the mean to color
				vname=hi_varnames[i];
				mapped_col=map_val_to_col(num_colors, 
					factors[, vname], hi_tab[vname, "member_mean"]);

				# Draw rectangle for the cluster/variable
				rect(column_pos[cl_ix], num_variables-i, 
					column_pos[cl_ix+1], num_variables-i+1,
					border=NA, 
					col=heat_col[mapped_col]);
				text(column_mids[cl_ix], num_variables-i-.5+1, vname, 
					adj=.5, font=2);
			}
		}

	}	

	# Label cluster numbers
	cuts=length(cl_mids);
	mtext(1:cuts, side=3, line=-1, at=cl_mids, font=2, 
		cex=ifelse((1:cuts)<10, 2, 1.5),
		col=1:cuts);

	# Label sample sizes
	if(num_clusters<5){
		mtext(paste("n=",cluster_sizes, sep=""), side=3, line=-2, at=cl_mids, font=3, cex=.7);
	}else if(num_clusters<10){
		mtext(cluster_sizes, side=3, line=-2, at=cl_mids, font=3, cex=.7);
	}else{
		mtext(cluster_sizes, side=3, line=-2, at=cl_mids, font=3, cex=.5);
	}


	# Draw lines and callouts
	for(i in 1:num_clusters){
		points(c(column_pos[i], column_pos[i]), c(0,num_variables), type="l", col="black");
		points(c(column_pos[i], cl_loc[i]), c(num_variables,num_variables+2), type="l", col="black");
		points(c(cl_loc[i], cl_loc[i]), c(num_variables+4,num_variables+2), type="l", col="black");
	}

	abline(v=c(0,num_subjects));
	
}

assign_cluster_names_from_features=function(cl_feat_rec, top_feat=3){
	# This function will look at the feature table and come up with
	# a 'short' name for the cluster

	num_clus=length(cl_feat_rec);
	cluster_name=character(num_clus);
	for(i in 1:num_clus){
		feat_tab=cl_feat_rec[[i]];
		
		ord_ix=order(feat_tab[,"p-value"]);
		num_signf=sum(feat_tab[,"signif"]);

		top_feat=min(top_feat, num_signf);
		top_tab=feat_tab[ord_ix[1:top_feat],,drop=F];

		feat_names=rownames(top_tab);
		feat_hilo=top_tab[,"isHigh"];

		cluster_name[i]=paste(
			"cl", sprintf("%02i", i), "_",
			paste(paste(ifelse(feat_hilo, "H", "L"), ".", feat_names, sep=""), collapse="_"),
			sep=""
		);
	}
	return(cluster_name);

}

append_export=function(fname, orig_mat, cl_mem, cl_names){

	num_cuts=length(cl_mem);
	subj_ids=rownames(orig_mat);
	orig_mat_colnames=colnames(orig_mat);

	# Appened columns to original factor file
	new_mat=orig_mat;
	found_cuts=c();
	for(k in 2:num_cuts){
	
		if(is.null(cl_mem[[k]])){
			next;
		}

		cur_mem=cl_mem[[k]];
		cur_names=cl_names[[k]];

		print(cur_mem);
		print(cur_names);

		new_mat=cbind(new_mat, cur_names[cur_mem[subj_ids]]);

		found_cuts=c(found_cuts, k);
	}

	# Rename column headers
	cut_names=paste(sprintf("K%02i", found_cuts));
	colnames(new_mat)=c(orig_mat_colnames, cut_names);

	# Write to file
	fh=file(fname, "w");
	cat(file=fh, "subject_id\t");
	close(fh);
	write.table(new_mat, file=fname, sep="\t", row.names=T, col.names=T, quote=F, append=T);

}

###############################################################################

plot_title_page("Dendrogram and Heatmaps", c(
	"The following plots are generated and grouped at various cuts (k).",
	"At the top of each page the hierarhical clustering is represented by a dendrogram at",
	"the top of the page.  A dotted blue line represents the height",
	"at which the tree was cut to achieve the k for that iteration.",
	"",
	"Below the dendrogram are heatmaps describing the nature of each cluster:",
	"  1.) Individual values",
	"",
	"  2.) Cluster values (average of individual values)",
	"If a cluster's values for a variable are significantly higher/lower then the other",
	"subjects not in the group, then a red/blue up/down triangle is drawn in the cell.",
	"",
	"  3.) Cluster descriptors",
	"These figures list the variables that were significantly greater/lower, ordered by",
	"p-value.  Note that the colors are assigned relative to all the values for a variable.",
	"This means a cluster may be high for variable, but it is still colored a shade",
	"of blue because overall the there are few observations of it.",
	"",
	"High values are colored red, and low values are colored blue.  Intermediate",
	"values are green."
));


if(MaxClusters<=2){
	max_cuts=ceiling(max(log2(num_subjects), num_subjects/40));
}else{
	max_cuts=MaxClusters;	
}
cat("Max Cuts: ", max_cuts, "\n");


cut_clusters=list();
clus_names=list();

palette_col=c("red", "green", "blue", "cyan", "magenta", "orange", "gray", "pink", "black", 
		"purple", "brown", "aquamarine");

palette(palette_col);

# Get the list of leaf names in dendrogram, from left to right
subj_ids_from_dendro_LtoR=get_clstrd_leaf_names(dend);
#print(subj_ids_from_dendro);

layout_mat=matrix(c(1,2,2,2), nrow=4);
layout(layout_mat);
par(oma=c(0,0,1,0));

left_label_margin=5;

for(k in 2:max_cuts){
	
	clmap=cutree(hcl, k=k);

	cl_mids=get_middle_of_groups(subj_ids_from_dendro_LtoR, clmap);
	#cat("Cluster Mids:\n");
	#print(cl_mids);

	# Change the membership IDs so they increase from left to right
	plot_order=order(cl_mids);
	clmap_tmp=numeric(num_subjects);
	for(gr_ix in 1:k){
		old_id=(clmap==plot_order[gr_ix]);
		clmap_tmp[old_id]=gr_ix;
	}
	names(clmap_tmp)=names(clmap);
	clmap=clmap_tmp;
	cl_mids=cl_mids[plot_order];
	cut_clusters[[k]]=clmap;

	# calculate cluster division points
	grp_sizes=table(clmap);
	print(grp_sizes);
	grp_split_loc=cumsum(grp_sizes)[1:(k-1)];

	# Find cut height
	h=find_height_at_k(hcl, k);

	# Compute cluster influencers
	cluster_features_rec=find_cluster_features(k, clmap, targeted_factors);

	#----------------------------------------------------------------------
	# By subject
	# Plot dendrogram with all the annotations
	par(mar=c(4,0,3,left_label_margin));
	plot_annotated_dendrogram(dend, clmap, k, cl_mids, h);
	abline(v=grp_split_loc, lwd=.5, col="grey");
	# Plot heatmap of variables
	par(mar=c(0,0,0,left_label_margin));
	plot_variables_heatmap(targeted_factors, subj_ids_from_dendro_LtoR, cl_mids=cl_mids);
	abline(v=grp_split_loc, lwd=1, col="black");

	#----------------------------------------------------------------------
	# By cluster
	# Plot heatmap with group averages and significances
	par(mar=c(4,0,3,left_label_margin));
	plot_annotated_dendrogram(dend, clmap, k, cl_mids, h);
	abline(v=grp_split_loc, lwd=.5, col="grey");
	# Plot heatmap of variables
	par(mar=c(0,0,0,left_label_margin));
	plot_variables_heatmap_avgs(targeted_factors, grp_split_loc, cluster_features_rec);
	abline(v=grp_split_loc, lwd=1, col="black");

	#----------------------------------------------------------------------
	# List cluster influencers
	# By cluster
	# Plot heatmap with group averages and significances
	par(mar=c(4,0,3,0));
	plot_annotated_dendrogram(dend, clmap, k, cl_mids, h);
	abline(v=grp_split_loc, lwd=.5, col="grey");
	# Plot heatmap of variables
	par(mar=c(0,0,0,0));
	plot_variables_lists(targeted_factors, grp_split_loc, cluster_features_rec);

	clus_names[[k]]=assign_cluster_names_from_features(cluster_features_rec);
			

}

cat("Exporting Clusters...\n");

append_export(paste(OutputFnameRoot,".loaded_factors.wClusters.tsv", sep=""),
	loaded_factors, cut_clusters, clus_names);

append_export(paste(OutputFnameRoot,".targeted_factors.wClusters.tsv", sep=""),
	targeted_factors, cut_clusters, clus_names);

###############################################################################
###############################################################################

# Calculate Calinkski-Harabasz statistic 

compute_pseudoF=function(dist_mat, memberships, resample_sequence){

        dist_mat=as.matrix(dist_mat);

        num_samples=length(memberships);
        num_clusters=length(unique(memberships));

        n=num_samples;
        k=num_clusters;

        sample_list=names(memberships);

        num_bs=nrow(resample_sequence);
        pseudo_f=numeric(num_bs);

	num_cells=num_samples*(num_samples-1)/2;

        for(bs_ix in 1:num_bs){
		cat(".");
                #SSB=0;
                #SSW=0;

		ssb_arr=numeric(num_cells);
		ssw_arr=numeric(num_cells);
		num_ssb=0;
		num_ssw=0;

                # Resample
                resamp_ix=resample_sequence[bs_ix,]; #sample(num_samples, replace=T);

                # Compute SSB and SSW
                for(i in 1:num_samples){
                        samp_i=sample_list[resamp_ix[i]];

                        for(j in 1:num_samples){

                                if(i<j){

                                        samp_j=sample_list[resamp_ix[j]];

                                        i_mem=memberships[samp_i];
                                        j_mem=memberships[samp_j];

                                        #cat(i_mem, "/", j_mem, "\n");
                                        if(i_mem==j_mem){
                                                #SSW=SSW+dist_mat[samp_i, samp_j];
                                                #cat("SSW: ", SSW, "\n");
						num_ssw=num_ssw+1;
						ssw_arr[num_ssw]=dist_mat[samp_i, samp_j];
                                        }else{
                                                #SSB=SSB+dist_mat[samp_i, samp_j];
                                                #cat("SSB: ", SSB, "\n");
						num_ssb=num_ssb+1;
						ssb_arr[num_ssb]=dist_mat[samp_i, samp_j];
                                        }
                                }
                        }
                }

                # Compute/Store this instance of bootstrap
                #pseudo_f[bs_ix]=SSB/(k-1)*(n-k)/SSW;
                #pseudo_f[bs_ix]=SSB/num_ssb*(num_ssw)/SSW;
                pseudo_f[bs_ix]=median(ssb_arr[1:num_ssb])/median(ssw_arr[1:num_ssw]);
        }
	cat("\n");

        return(pseudo_f);
}

#------------------------------------------------------------------------------

bootstrap_calculate_pseudoF_stats=function(distmat, membership_list, num_bs){

	num_samples=nrow(distmat);

	# Pre-randomize samples across all cluster sizes so finding max will be consistent
	resample_sequence=matrix(0, nrow=num_bs, ncol=num_samples);
	# Rows are each BS instance, Columns are the selected samples
	for(i in 1:num_bs){
		resample_sequence[i,]=sample(num_samples, replace=T);
	}
	#print(resample_sequence);

	max_clusters=length(membership_list);

	# Store the pseudoF stat for each cluster cut
	pseudoF_mat=matrix(NA, nrow=num_bs, ncol=max_clusters);

	# For each cluster cut, bootstrap calculate the pseudoF stat
	for(num_cl in 2:max_clusters){
		cat("Performing bootstrap pseudo-F on k=", num_cl, "\n");
		memberships=membership_list[[num_cl]];
		pseudoF_bs_res=compute_pseudoF(distmat, memberships, resample_sequence);
		pseudoF_mat[,num_cl]=pseudoF_bs_res;
	}

	# Remove infinite values, probably due to low sample size
	pseudoF_mat=apply(pseudoF_mat, 1:2, function(x){ if(is.infinite(x)){return(NA);}else{return(x);}});

	med_pseudoF=apply(pseudoF_mat, 2, function(x){quantile(x, .5, na.rm=T)});
	lb_pseudoF=apply(pseudoF_mat, 2, function(x){quantile(x, .025, na.rm=T)});
	ub_pseudoF=apply(pseudoF_mat, 2, function(x){quantile(x, .975, na.rm=T)});
	max_of_median_pseudoF=max(med_pseudoF);
	k_of_max=min(which(max_of_median_pseudoF==med_pseudoF));

	min_pseudoF=min(pseudoF_mat, na.rm=T);
	max_pseudoF=max(pseudoF_mat, na.rm=T);

	results=list();
	results[["medians"]]=med_pseudoF;
	results[["lb"]]=lb_pseudoF;
	results[["ub"]]=ub_pseudoF;
	results[["max_k"]]=k_of_max;
	results[["matrix"]]=pseudoF_mat;
	results[["min"]]=min_pseudoF;
	results[["max"]]=max_pseudoF;
	results[["num_bs"]]=num_bs;

	return(results);

}

#------------------------------------------------------------------------------

pf=bootstrap_calculate_pseudoF_stats(distance_mat, cut_clusters, num_bs=80);

print(pf);

# Plot cluster separation

par(mfrow=c(1,1));
par(mar=c(5,5,4,1));
plot(0, type="n", xlim=c(1, max_cuts+1), ylim=c(pf[["min"]], pf[["max"]]),
	xlab="k", ylab="Pseudo-F", main="Cluster Cut vs. Separation"
	);

scat=rnorm(pf[["num_bs"]], mean=0, sd=.10);

for(k in 2:max_cuts){
	cat("Working on plotting k: ", k, "\n");
	points(scat+k, pf[["matrix"]][,k], col="grey25");
	points(k, pf[["ub"]][k], pch="-", cex=2, col="red")
	points(k, pf[["lb"]][k], pch="-", cex=2, col="blue")
	points(k, pf[["medians"]][k], pch="-", cex=4, col="green")
}

points(2:max_cuts, pf[["medians"]][2:max_cuts], col="black", type="l", lwd=2.5);

###############################################################################

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
