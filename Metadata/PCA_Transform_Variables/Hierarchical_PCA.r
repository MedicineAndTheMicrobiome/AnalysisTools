#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

options(useFancyQuotes=F);
options(digits=5);
options(width=120);

params=c(
	"factors", "f", 1, "character",
	"outputroot", "o", 1, "character",
	"groupings", "g", 1, "character",
	"pc_contrib_cutoff", "p", 2, "numeric"
);

NO_CHANGE="orig";
NORM_PVAL_CUTOFF=0.20;
MIN_PC_PROP_CUTOFF=0.10;

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <factors/metadata file name>\n",
	"	-o <output filename root>\n",
	"	-g <variable groupings>\n",
	"\n",
	"	[-p <Min PC contribution cutoff, default=", MIN_PC_PROP_CUTOFF," >]\n",
	"\n",
	"This script will:\n",
	"	1.) Transform all variables to ensure they are normally distributed\n",
	"	2.) Perform PCA on each group of variables.  A variable can be in more than 1 group.\n",
	"	3.) Select the top PC's that contribute more than the min cutoff. e.g. >10%\n",
	"	4.) Annotate the top Selected PC's.\n",
	"	5.) Perform Overall PCA on Groups PCA\n",
	"	6.) Annotate the top Overall PCs with Group names\n",
	"	7.) Export Group and Overall PC Scores.\n",
	"\n");

if(
	!length(opt$factors) || 
	!length(opt$outputroot) || 
	!length(opt$groupings)
){
	cat(usage);
	q(status=-1);
}

FactorsFname=opt$factors;
OutputFnameRoot=opt$outputroot;
GroupingsFname=opt$groupings;

PCContribCutoff=MIN_PC_PROP_CUTOFF;

if(length(opt$pc_contrib_cutoff)){
	PCContribCutoff=opt$pc_contrib_cutoff;
}


param_text=capture.output({
	cat("\n");
	cat("Factor File Name: ", FactorsFname, "\n");
	cat("Output File Name Root: ", OutputFnameRoot, "\n");
	cat("Groupings File Name: ", GroupingsFname, "\n");
	cat("PC Min Coverage: ", PCContribCutoff, "\n");
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

test_and_apply_log_transform=function(mat_val, pval_cutoff=.2, plot_before_after=T){
	nrows=nrow(mat_val);
	ncols=ncol(mat_val);

	trans_mat=mat_val;
	orig_colnames=colnames(mat_val);
	new_colnames=character();

	if(plot_before_after){
		orig_par=par(no.readonly=T);
		par(mfrow=c(5,2));
	}

	delete_list=c();
	for(var in orig_colnames){
		values=mat_val[,var];


		transformed=F;

		num_unique_val=length(setdiff(unique(values), NA));
		values_nona=values[!is.na(values)];
		num_nona=length(values_nona);

		if(!is.numeric(values_nona)){
			cat("Error: Values not numeric for: ", var, "\n", sep="");
			print(values_nona);
		}

		if(num_nona<=3){
			cat("Not enough non NA values to measure normality.\n");
			new_colnames=c(new_colnames, var);
			next;
		}

		if(any(values_nona<0)){
			cat("Negative values.  Skipping transformation.\n");
			new_colnames=c(new_colnames, var);
			next;
		}

		cat("\n", var, ": Num Unique Values: ", num_unique_val, "\n");

		if(num_unique_val>1){

			test_res=shapiro.test(values);

			if(test_res$p.value<=pval_cutoff && num_unique_val>2){
				cat(" Not normal: ", test_res$p.value, "\n");
				if(any(values_nona==0)){
					log_values=log(values+1);
				}else{
					log_values=log(values);
				}

				test_log_res=shapiro.test(log_values);

				if(test_log_res$p.value < test_res$p.value){
					# Keep original
					cat("  No Improvement: ", test_log_res$p.value, "\n");
					new_colnames=c(new_colnames, paste("orig_", var, sep=""));
				}else{
					# Keep log transformed
					cat("  Transformation Effective: ", test_log_res$p.value, "\n");
					new_colnames=c(new_colnames, paste("log_", var, sep=""));
					trans_mat[, var]=log_values;
					transformed=T;		
				}
			}else{
				cat(" Normal enough. ", test_res$p.value, "\n");
				new_colnames=c(new_colnames, var);
			}

		}else{
			cat("  All values identical, removing...\n");
			new_name=paste("all_ident_", var, sep="");
			new_colnames=c(new_colnames, new_name);
			delete_list=c(delete_list, new_name);
		}

		if(plot_before_after){
			nclass=nclass.Sturges(values)*4;

			hist(values, main=var, breaks=nclass);
			title(main=sprintf("p-value: %4.4f", test_res$p.value), cex.main=.8, line=.5);
			
			if(transformed){
				hist(log_values, breaks=nclass, main=paste("log(", var,")", sep=""));
				title(main=sprintf("p-value: %4.4f", test_log_res$p.value), cex.main=.8, line=.5);
			}else{
				plot(0,0, xlab="", ylab="", main="", xaxt="n", yaxt="n", bty="n", type="n");
				title(main=sprintf("p-value: %4.4f", test_log_res$p.value), cex.main=.8, line=.5);
				text(0,0, "Transform not beneficial");
			}
		}

	}

	colnames(trans_mat)=new_colnames;

	mapping=new_colnames;
	names(mapping)=orig_colnames

	trans_mat=trans_mat[,setdiff(new_colnames, delete_list),drop=F];

	if(plot_before_after){
		par(orig_par);
	}

	trans_rec=list();
	trans_rec[["Transformed"]]=trans_mat;
	trans_rec[["NamesMap"]]=mapping;

	return(trans_rec);
}

remap_groupings=function(old_grp_rec, curtd_fac_rec){

	new_grp_rec=list();
	curtd_var=c();

	kept_var=colnames(curtd_fac_rec[["Transformed"]]);

	grps=old_grp_rec[["Groups"]];
	for(grp in grps){

		# Remove variables from group that are were removed from factor because they were degenerate
		translated=intersect(kept_var, curtd_fac_rec[["NamesMap"]][old_grp_rec[["GrpVarMap"]][[grp]]]);

		new_grp_rec[["GrpVarMap"]][[grp]]=translated;
		curtd_var=c(curtd_var, translated);
	}

	new_grp_rec[["Groups"]]=grps;
	new_grp_rec[["NumGroups"]]=old_grp_rec[["NumGroups"]];
	new_grp_rec[["Variables"]]=sort(unique(curtd_var));
	new_grp_rec[["NumUniqVar"]]=length(new_grp_rec[["Variables"]]);

	return(new_grp_rec);

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

compute_correlations=function(mat){
	num_col=ncol(mat);
	cor_mat=matrix(0, nrow=num_col, ncol=num_col);
	pval_mat=matrix(0, nrow=num_col, ncol=num_col);
	rownames(cor_mat)=colnames(mat);
	colnames(cor_mat)=colnames(mat);
	rownames(pval_mat)=colnames(mat);
	colnames(pval_mat)=colnames(mat);
	for(i in 1:num_col){
		for(j in 1:i){
			v1=mat[,i];
			v2=mat[,j];
			notna=!(is.na(v1) | is.na(v2));
			#cor_mat[i,j]=cor(v1[notna], v2[notna]);
			test=cor.test(v1[notna], v2[notna]);
			pval_mat[i,j]=test$p.value;
			pval_mat[j,i]=test$p.value;
			cor_mat[i,j]=test$estimate;
			cor_mat[j,i]=test$estimate;
		}
	}
	res=list();
	res[["val"]]=cor_mat;
	res[["pval"]]=pval_mat;;
	res[["dist"]]=as.dist(1-abs(cor_mat));
	return(res);
}

#-----------------------------------------------------------------------------#

calculate_grouped_correlations=function(cur_fact_rec, groupings_rec){

	groups=groupings_rec[["Groups"]];
	kept_factors=names(cur_fact_rec[["Transformed"]]);

	grp_cor_rec=list();

	for(cur_grp in groups){
		cat("Calculating Correl on: ", cur_grp, "\n");
		grp_mem=groupings_rec[["GrpVarMap"]][[cur_grp]];

		grp_var_subset=cur_fact_rec[["Transformed"]][grp_mem];
		if(ncol(grp_var_subset)>0){
			grp_cor_rec[[cur_grp]]=compute_correlations(grp_var_subset);
		}
	}

	return(grp_cor_rec);

}

correl_to_dendrogram=function(cor_rec, title){

	# Generate dendrogram
	# Barplot of PC contributions
	# Ordination of top 2 PCs
	# Annotated PCs selected PCs

	if(nrow(cor_rec$val) < 2){
		return();
	}

	par(mfrow=c(2,1));
	par(mar=c(20,4,1,2));
	hcl=hclust(cor_rec$dist, method="ward.D2");
	dend=as.dendrogram(hcl);

	plot(dend, main=title, ylab="1-|cor(x,y)|", ylim=c(0, 2));

}

correl_to_PC=function(cor_rec, values, title, contrib_cutoff){

	num_var=nrow(cor_rec$val);
	if(num_var < 2){
		return();
	}

	# Compute PCA
	eigen=eigen(cor_rec$val);

	# Compute variance contribution of each PC
	pca_propvar=eigen$values/sum(eigen$values);
	pca_propcumsum=cumsum(pca_propvar);

	keep_ix=pca_propvar>contrib_cutoff;
	num_pc_at_cutoff=sum(keep_ix);

	
	# Compute per sample scores
	scores=(scale(values, center=T, scale=T) %*% eigen$vectors);

	cat("Num PC's exceeding min cutoff: ", num_pc_at_cutoff, "\n");
	cat("Num Variables: ", num_var, "\n");
	var_names=colnames(values);
	pcname=character(num_pc_at_cutoff);
	for(pcix in 1:num_pc_at_cutoff){

		# Find correlation with original variables
		correls=numeric(num_var);
		for(valix in 1:num_var){
			correls[valix]=cor(scores[, pcix], values[, valix]);
		}

		# Find original variable most similar to PC
		mag_cor=abs(correls);
		max_cor=max(mag_cor);
		closest_var_to_pc=which(max_cor==mag_cor);

		# If PC is negative, multiple score by -1
		if(correls[closest_var_to_pc]<0){
			scores[, pcix]=-1*scores[, pcix];	
		}

		# Construct Name of PC
		#   As:  PC[ix]_[coverage].[varname]_[cor]
		pcname[pcix]=paste("PC",pcix,"_", round(pca_propvar[pcix]*100), ".",
			var_names[closest_var_to_pc],"_", round(max_cor*100), sep="");
		
	}

	kept_score=scores[,1:num_pc_at_cutoff, drop=F];
	colnames(kept_score)=pcname;

	########################################################################

	par(mfrow=c(2,1));
	par(mar=c(20,4,1,1));
	# Generate Dendrograms w/ PC
	value_and_pc=cbind(values, kept_score);
	pc_and_val=cor(value_and_pc);
	cor_as_dist=1-(abs(pc_and_val));
	hcl=hclust(as.dist(cor_as_dist), method="ward.D2");
        dend=as.dendrogram(hcl);
        plot(dend, main=title, ylab="1-|cor(x,y)|", ylim=c(0, 2));

	# Generate Barplots w/ PC Variance
	colors=rep("khaki", num_var);
	colors[1:num_pc_at_cutoff]="slategray2";
	barlabs=rep("", num_var);
	barlabs[1:num_pc_at_cutoff]=pcname;
	par(mar=c(20,4,0,1));
	mids=barplot(pca_propvar, ylim=c(0,1.1), names.arg=barlabs,  las=2,
		col=colors, ylab="Proportion of Variance");
	text(mids, pca_propvar, sprintf("%3.1f", pca_propvar*100), pos=3, cex=.9);
	abline(h=contrib_cutoff, col="red", lty=2);

	return(kept_score);
}

##############################################################################
# Main Program Starts Here!
##############################################################################

pdf(paste(OutputFnameRoot, ".pdf", sep=""), height=11, width=8.5);

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
groupings_rec=load_groupings(GroupingsFname);

# Reconcile Grouping Map with Factor Variables
num_grpmap_var=length(groupings_rec$Variables);
num_factor_var=length(loaded_factor_names);
shared_variables=sort(intersect(groupings_rec$Variables, loaded_factor_names));
num_shared_var=length(shared_variables);
unmapped_var=setdiff(loaded_factor_names, groupings_rec$Variables);

shared_text=capture.output({
	cat("Num Factor Variables: ", num_factor_var, "\n");
	cat("Num Unique Mapped Variables: ", num_grpmap_var, "\n");
	cat("Num Shared Variables: ", num_shared_var, "\n");
	cat("\n");
	cat("Unmapped Factor Variables: \n");
	print(unmapped_var);
});

plot_text(shared_text);
print(shared_text);

grouping_text=capture.output({
	cat("Num Groups: ", groupings_rec[["NumGroups"]], "\n");
	cat("Num (unique) Variables: ", groupings_rec[["NumUniqVar"]], "\n");
});
plot_text(grouping_text);
print(grouping_text);

shared_factors_mat=loaded_factors[,shared_variables, drop=F];
#print(shared_factors);

##############################################################################

cat("Testing and Apply Log Transforms...\n");
curated_fact_rec=test_and_apply_log_transform(shared_factors_mat, NORM_PVAL_CUTOFF);

remapped_groupings=remap_groupings(groupings_rec, curated_fact_rec);

group_cor_rec=calculate_grouped_correlations(curated_fact_rec, remapped_groupings);

for(grp_ix in names(group_cor_rec)){
	cat("------------------------------------------------------------\n");
	cat("Analyzing: ", grp_ix, "\n");

	grp_members=remapped_groupings[["GrpVarMap"]][[grp_ix]];
	subset_curated_fact=curated_fact_rec[["Transformed"]][,grp_members];
	correl_to_PC(group_cor_rec[[grp_ix]], subset_curated_fact, title=grp_ix, PCContribCutoff);
}

##############################################################################









quit();



##############################################################################

colnames(positive_scores)=pc_name;

correl_wpca=compute_correlations(cbind(positive_scores[,1:num_pc_at_cutoff], curated_pred_mat, curated_resp_mat));
hcl=hclust(correl_wpca$dist, method="ward.D2");
dend=as.dendrogram(hcl);

highlight_pcs=function(x){
	if(is.leaf(x)){
		leaf_attr=attributes(x);
		label=leaf_attr$label;
		print(label);
		if(any(label==pc_name)){
			color="red";
			font=2;
		}else{
			color="black";
			font=1;
		}
		attr(x, "nodePar")=c(leaf_attr$nodePar, list(lab.font=font, lab.col=color, cex=0));
	}
	return(x);
}

dend=dendrapply(dend, highlight_pcs);

par(mfrow=c(1,1));
par(mar=c(2,1,1,20));
plot(dend, horiz=T, main="Ward's Minimum Variance: dist(1-abs(cor)) with PCs");

##############################################################################

out_factors=loaded_factors;

##############################################################################

cat("Outputing New Factor File Values:\n");
fname=paste(OutputFnameRoot, ".pca.tsv", sep="");
fh=file(fname, "w");
cat(file=fh, "SampleID");
close(fh);
write.table(out_factors, file=fname, col.names=NA, append=T, quote=F, sep="\t");

##############################################################################

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
