#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

source("~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r");
source("~/git/AnalysisTools/Metadata/OutputFileLibrary/OutputFileLibrary.r");


options(useFancyQuotes=F);

params=c(
	"summary_file_pattern", "p", 2, "character",
	"summary_file_list_fn", "f", 2, "character",
	"output_fn_root", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

NUM_TOP_CAT=35;
ALPHA=.05;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	Summary File List:\n",
	"	[-p <Specify a pattern, e.g. \"GO.cellular_component/*names.summary_table.tsv\">]\n",
	"	[-f <Specify a list, e.g. summary_tables.lst>\n",
	"	-o <output filename root>\n",
	"\n",
	"This script will:\n",
	"   1.) read in a set/list of summary tables.\n",
	"   2.) perform CLR on each summary table independently.\n",
	"\n", sep="");

if(
	!length(opt$output_fn_root)
){
	cat(usage);
	q(status=-1);
}

OutputRoot=opt$output_fn_root;

STFilePattern="";
STFileList="";

if(length(opt$summary_file_pattern)){
	STFilePattern=opt$summary_file_pattern;
}

if(length(opt$summary_file_list)){
	STFileList=opt$summary_file_list;
}

cat("\n");
cat("Specified Summary File Pattern: ", STFilePattern, "\n");
cat("Specified Summary File List: ", STFileList, "\n");
cat("Output File: ", OutputRoot, "\n", sep="");
cat("\n");

##############################################################################

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

##############################################################################

find_minimal_unique_name=function(name_arr){
	# This function will take an array of names (paths and file names)
	# and trim redundant fields from the front and end, until the middle
	# fields make a unique string
	
	cat("In:\n");
	print(name_arr);

	if(length(name_arr)>1){

		toks=strsplit(name_arr, "[/\\.]");

		# Begining
		trim=T;
		while(trim){
			first=sapply(toks, function(x){x[1]});
			uniqs=unique(first);
			num_uniq=length(uniqs);
			if(num_uniq==1){
				toks=lapply(toks, function(x){x=x[2:length(x)]});
				trim=T;
			}else{
				trim=F;
			}
			#print(toks);
		}

		# Ending 
		trim=T;
		while(trim){
			last=sapply(toks, function(x){tail(x,1)});
			uniqs=unique(last);
			num_uniq=length(uniqs);
			if(num_uniq==1){
				toks=lapply(toks, function(x){x=x[1:(length(x)-1)]});
				trim=T;
			}else{
				trim=F;
			}
			#print(toks);
		}

		reconstr_names=unlist(sapply(toks, function(x){paste(x, collapse=".");}));

	}else{
		reconstr_names=name_arr;
	}

	cat("Out:\n");
	print(reconstr_names);

	return(reconstr_names);

}

##############################################################################

topN=function(norm_mat, max){
	mean_abund=apply(norm_mat, 2, mean);
	sort_ix=order(mean_abund, decreasing=T);
	sorted_norm_mat=norm_mat[,sort_ix,drop=F];
	sorted_mean_abund=mean_abund[sort_ix];
	sorted_cat_names=colnames(norm_mat)[sort_ix];

	# Take Remaining out for now
	remaining_ix=(sorted_cat_names=="Remaining");
	remaining_found=0;
	remaining_values=c();
	if(any(remaining_ix)){
		remaining_found=1;	
		remaining_values=sorted_norm_mat[,remaining_ix, drop=F];
		sorted_norm_mat=sorted_norm_mat[,!remaining_ix, drop=F];
	}
	
	# Extract top
	num_col=ncol(sorted_norm_mat);
	max=min(max, num_col);
	out_mat=sorted_norm_mat[,1:max, drop=F];

	rare_exists=(max<num_col);
	if(rare_exists){
		cat("Combining into Rare.\n");
		rare_mat=sorted_norm_mat[,(max+1):num_col, drop=F];
		rare_comb=apply(rare_mat, 1, sum);
		out_cname=colnames(out_mat);
		out_mat=cbind(out_mat, rare_comb);
		colnames(out_mat)=c(out_cname, "Rare");
	}

	if(remaining_found){
		cat("Re-including Remaining.\n");
		out_cname=colnames(out_mat);
		out_mat=cbind(out_mat, remaining_values);
		colnames(out_mat)=c(out_cname, "Remaining");
	}

	return(out_mat);
}

##############################################################################

set_zeros_globalmin=function(norm_mat){
	num_col=ncol(norm_mat);
	norm_vect=as.vector(norm_mat);

	z_ix=(norm_vect==0);
	nz_ix=!z_ix;
	
	min_nz=min(norm_vect[nz_ix]);
	zero_sub=min_nz/10;

	norm_vect[z_ix]=zero_sub;

	zsub_mat=matrix(norm_vect, ncol=num_col);
	colnames(zsub_mat)=colnames(norm_mat);
	rownames(zsub_mat)=rownames(norm_mat);

	return(zsub_mat);
}

set_zeros_catmin=function(norm_mat){
	num_col=ncol(norm_mat);

	for(i in 1:num_col){
		norm_vect=norm_mat[,i];
		z_ix=(norm_vect==0);
		nz_ix=!z_ix;
		min_nz=min(norm_vect[nz_ix]);
		zero_sub=min_nz/10;
		norm_vect[z_ix]=zero_sub;
		norm_mat[,i]=norm_vect;
	}

	return(norm_mat);
}

set_zeros_catdev=function(norm_mat, dev_sep=7){
	num_col=ncol(norm_mat);

	log_sd=apply(norm_mat, 2, function(x){
			nz_ix=(x!=0)
			nz_val=x[nz_ix];
			if(length(nz_val)<2){
				stdev=0;
			}else{
				stdev=sd(log(nz_val));
			}
			return(stdev);
		});

	mean_log_sd=mean(log_sd, na.rm=T);
	min_log_sd=min(log_sd, na.rm=T);
	max_log_sd=max(log_sd, na.rm=T);
	cat("Global Mean Log(sd): ", mean_log_sd, "\n");
	cat("Global Min Log(sd): ", min_log_sd, "\n");
	cat("Global Max Log(sd): ", max_log_sd, "\n");

	for(i in 1:num_col){
		norm_vect=norm_mat[,i];
		z_ix=(norm_vect==0);
		nz_ix=!z_ix;

		log_nz=log(norm_vect[nz_ix]);
		log_nz_mean=mean(log_nz);
		log_nz_sd=sd(log_nz);

		if(is.na(log_nz_sd)){
			log_nz_sd=mean_log_sd;
		}

		zero_sub=exp(log_nz_mean-dev_sep*log_nz_sd);
		norm_vect[z_ix]=zero_sub;
		norm_mat[,i]=norm_vect;
	}

	return(norm_mat);
}

lograt_trans=function(norm_mat){

	cat_names=colnames(norm_mat);
	num_col=ncol(norm_mat);
	remaining_ix=(cat_names=="Remaining");
	remaining_found=any(remaining_ix);


	if(remaining_found){
		# Do ALR, with Remaining as denominator
		cat("Performing ALR with Remaining...\n");

		# Split remaining from field
		norm_wo_rem=norm_mat[,!remaining_ix,drop=F];
		remainder=norm_mat[,remaining_ix];

		# Calc alr
		ratio_mat=norm_wo_rem / remainder;

	}else{
		# Do CLR, with geometric mean as denominator
		cat("Performing CLR with Geometric Mean...\n");

		geomean=function(x){
			prod(x)^(1/length(x));
		}

		gm=apply(norm_mat, 1, geomean);	
		
		# Calc clr
		ratio_mat=norm_mat / gm;
	}	

	lr_mat=log(ratio_mat);

	return(lr_mat);

}

##############################################################################

generate_overlayed_distributions=function(norm_mat){

	par(mfrow=c(4,1));
	num_cat=ncol(norm_mat);

	num_samp=nrow(norm_mat);
	lres=apply(norm_mat, 2, function(x){ 
		z_ix=(x==0);
		wgt=rep(1, num_samp);
		wgt[z_ix]=0;
		x[!z_ix]=log(x[!z_ix]);
		density(x, weights=wgt);		
	});

	plot(0, type="n", xlim=c(-14,4), ylim=c(0, 2000));
	for(i in 1:num_cat){
		points(lres[[i]], type="l", col=i);
	}
}

generate_histograms=function(norm_mat, lograt_mat){

	norm_cat=colnames(norm_mat);
	lograt_cat=colnames(lograt_mat);
	shared_cat=intersect(norm_cat, lograt_cat);
	
	par(mfrow=c(4,3));

	for(cat_ix in shared_cat){

		log_norm=log(norm_mat[,cat_ix]);
		lograt=lograt_mat[,cat_ix];

		hist(log_norm, xlab="Log(Norm)", main=cat_ix);
		hist(lograt, xlab="ALR or CLR", main=cat_ix);
		plot(log_norm, lograt, xlab="Log(Norm)", ylab="ALR or CLR", main=cat_ix);
	}

}

##############################################################################

make_fields_unique=function(lograt_lst){
	# Tries to make fields unique and associated with the the group
	# they go with.

	num_groups=length(lograt_lst);
	group_names=names(lograt_lst);

	out_list=list();

	for(i in 1:num_groups){

		tab=lograt_lst[[i]];

		new_tabcname=paste("G", sprintf("%02i", i), ".", colnames(tab), sep="");
		new_grpname=paste("G", sprintf("%02i", i), ".", group_names[i], sep="");
		
		colnames(tab)=new_tabcname;

		out_list[[new_grpname]]=tab;
	}	

	return(out_list);
}

export_groups=function(fn, lograt_lst){
	# Generates a 2 column file: <variable>\t<group>\n
	
	num_groups=length(lograt_lst);
	group_names=names(lograt_lst);
	
	cat("Exporting Groups [", num_groups, "]:\n");
	print(group_names);
	
	out_mat=matrix(NA, ncol=2, nrow=0);
	
	for(i in 1:num_groups){
		cur_grp=group_names[i]
		members=colnames(lograt_lst[[cur_grp]]);
		num_members=length(members);

		acc_mat=cbind(members, rep(cur_grp, num_members));
		out_mat=rbind(out_mat, acc_mat);
	}

	print(out_mat);

	write.table(out_mat, file=fn, quote=F, sep="\t", row.names=F, col.names=F);

}

export_joined_tables=function(fn, lograt_lst){
	# Generates full matrix output file by sample_id
	
	num_groups=length(lograt_lst);
	group_names=names(lograt_lst);
	
	cat("------------------------------------------------------------\n");
	cat("Exporting Tables [", num_groups, "]:\n");
	print(group_names);
	
	# Make sure all tables have the same sample IDs
	sample_ids=rownames(lograt_lst[[1]]);
	category_ids=c();
	num_samples=length(sample_ids);
	for(i in 1:num_groups){

		cur_tab=lograt_lst[[i]];

		cur_samp_ids=rownames(cur_tab);
		shared=intersect(sample_ids, cur_samp_ids)
		if(length(shared)!=num_samples){
			cat("Error:  Inconsistent Sample IDs.\n");
			cat("Group ID: ", group_names[i], "\n");
			
			if(length(sample_ids)>length(cur_samp_ids)){
				print(setdiff(sample_ids, cur_samp_ids));
			}else{
				print(setdiff(cur_samp_ids, sample_ids));
			}
			quit(status=-1);
		}else{
			category_ids=c(category_ids, colnames(cur_tab));
		}
	}

	num_categories=length(category_ids);

	cat("Combined Out Table Dimenension:\n");	
	cat("Rows: ", num_samples, "\n");
	cat("Cols: ", num_categories, "\n");

	out_matrix=matrix(NA, nrow=num_samples, ncol=num_categories);
	rownames(out_matrix)=sample_ids;
	colnames(out_matrix)=category_ids;

	# Populate
	for(i in 1:num_groups){
		cat("Accumulating Table: ", group_names[i], "\n");
		cur_tab=lograt_lst[[i]];
		cnames=colnames(cur_tab);
		cat("Columns:\n");
		print(cnames);
		out_matrix[sample_ids, cnames]=cur_tab[sample_ids, cnames];
		cat("\n");
	}

	# Move sample ID to first column
	out_matrix=cbind(sample_ids, out_matrix);
	colnames(out_matrix)=c("SampleID", category_ids);

	# Write table to file
	write.table(out_matrix, file=fn, quote=F, sep="\t", row.names=F, col.names=T);

	cat("------------------------------------------------------------\n");
}

generate_correlation_matrix=function(mat, name){
	cormat=cor(mat);
	paint_matrix(cormat, name, plot_min=-1, plot_max=1, deci_pts=2, show_leading_zero=F);
}

##############################################################################

pdf(paste(OutputRoot, ".lr_trans.pdf", sep=""), height=11, width=8.5);

st_files=c();
if(STFilePattern!=""){
	target_dir=dirname(STFilePattern);
	file_pat=basename(STFilePattern);

	cat("Target Dir: ", target_dir, "\n");
	cat("File Pattern: ", file_pat, "\n");

	st_files=list.files(path=target_dir, pattern=file_pat, full.names=T)
}

if(STFileList!=""){
	st_files=scan(STFileList, what=character());
}

cat("Files Found:\n");
print(st_files);

if(length(st_files)==0){
	cat("Error or Warning: No files found.\n");
}

targ_st_ids=find_minimal_unique_name(st_files);

filenames_list=list();
counts_list=list();
norm_list=list();
topN_list=list();
zsub_list=list();
lograt_list=list();

filenames_list=as.list(st_files)
names(filenames_list)=targ_st_ids;

num_targets=length(targ_st_ids);
index=1;

for(targ_ix in targ_st_ids){
	
	cat("-----------------------------------------------------------------\n");
	cat("Working on: ", targ_ix, " ", index, "/", num_targets, "\n");
	
	plot_title_page(paste(index, "/", num_targets, sep=""), gsub("\\.", "\n", targ_ix));

	counts_list[[targ_ix]]=load_summary_file(filenames_list[[targ_ix]]);
	norm_list[[targ_ix]]=normalize(counts_list[[targ_ix]]);
	topN_list[[targ_ix]]=topN(norm_list[[targ_ix]], max=25);
	zsub_list[[targ_ix]]=set_zeros_catdev(topN_list[[targ_ix]]);
	lograt_list[[targ_ix]]=lograt_trans(zsub_list[[targ_ix]]);

	#generate_overlayed_distributions(topN_list[[targ_ix]]);	
	generate_histograms(zsub_list[[targ_ix]], lograt_list[[targ_ix]]);

	generate_correlation_matrix(lograt_list[[targ_ix]], targ_ix);
	
	cat("-----------------------------------------------------------------\n");
	index=index+1;
	
}

# Append group and member names to make them unique
lograt_list=make_fields_unique(lograt_list);

# Export Groups
outgrp_fn=paste(OutputRoot, ".tr_trans.map", sep="");
export_groups(outgrp_fn, lograt_list);

# Export Combined LR 
outgrp_fn=paste(OutputRoot, ".tr_trans.tsv", sep="");
export_joined_tables(outgrp_fn, lograt_list);

##############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
