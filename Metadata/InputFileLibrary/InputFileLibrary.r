
# This is a collection of shared functions that make the paired analyses
# consistent in how the metadata and samples are reconciled.

source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');

print_se=function(x){
	sink(stderr());
	print(x);
	sink();
}

#------------------------------------------------------------------------------

specified=function(var){

	#cat("Checking: ", var, "\n");

	if(length(var)==1){
		if(all(is.na(var)) || all(var=="") || all(is.na(var))){
			return(F);
		}else{
			return(T);
		}
	}else if(length(var)>1){
		return(T);
	}else{
		return(F);
	}
}

#------------------------------------------------------------------------------

load_list=function(fname, memo=NULL){

	if(!specified(fname)){
		message("File name NULL, no list specified for: ", memo);
		return(NULL);
	}else{
		message("Loading List: ", fname, " (purpose: ", memo, ")", sep="");
 		val=scan(fname, what=character(), comment.char="#");
		num_vals=length(val);
		message("  Number of items loaded: ", num_vals, "\n");

		message();
		message("Load List: ", fname);
		message("  (", memo, ")");
		message("  Number of items loaded: ", num_vals);
		message();
		print_se(val);

		return(val);
	}
}


#------------------------------------------------------------------------------

remove_zero_count_categories_and_samples=function(mat, max_report=10){

	# Remove zero count samples
	cat("Checking for samples with no counts.\n");
	tot=apply(mat, 1, sum);
	nonzero=tot>0;
	if(!(all(nonzero))){

		num_found=sum(tot==0);

		cat("  WARNING: Zero count samples found: ", num_found, "\n");
		samp_names=rownames(mat);

		print(head(samp_names[!nonzero], max_report));
		if(num_found>max_report){
			cat("...\n");
		}

		cat("\n");
		mat=mat[nonzero,,drop=F];
	}else{
		cat("  OK.  None found.\n");
	}

	# Remove zero count categories 
	cat("Checking for categories with no counts.\n");
	tot=apply(mat, 2, sum);
	nonzero=tot>0;
	if(!(all(nonzero))){

		num_found=sum(tot==0);

		cat("  WARNING: Zero count categories found: ", num_found, "\n");
		cat_names=colnames(mat);
		
		print(head(cat_names[!nonzero], max_report));
		if(num_found>max_report){
			cat("...\n");
		}

		cat("\n");
		mat=mat[,nonzero,drop=F];
	}else{
		cat("  OK.  None found.\n");
	}

	return(mat);
}

#------------------------------------------------------------------------------

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

#------------------------------------------------------------------------------

reorder_by_decreasing_abundance=function(normalized_mat){

	mean_abund=apply(normalized_mat, 2, mean);
	ix=order(mean_abund, decreasing=TRUE);
	normalized_mat=normalized_mat[, ix, drop=F];
	return(normalized_mat);
}

#------------------------------------------------------------------------------

shorten_category_names=function(full_names, split_char){

        splits=strsplit(full_names, split_char);
        short_names=character();
        for(i in 1:length(full_names)){
                short_names[i]=tail(splits[[i]], 1);
                short_names[i]=gsub("_unclassified$", "_uncl", short_names[i]);
                short_names[i]=gsub("_group", "_grp", short_names[i]);
        }
	return(short_names);
}

#------------------------------------------------------------------------------

clean_category_names=function(cat_names){

        cat_names=gsub("-", "_", cat_names);
	cat_names=gsub("\\[", "", cat_names);
	cat_names=gsub("\\]", "", cat_names);
	cat_names=gsub("\\(", "", cat_names);
	cat_names=gsub("\\)", "", cat_names);

	return(cat_names);

}

#------------------------------------------------------------------------------

load_summary_file=function(fname){

	cat("Loading Summary Table: ", fname, "\n");

        inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, 
		comment.char="", quote="", row.names=1))

	# Extract out counts
        counts_mat=inmat[,2:(ncol(inmat))];
	loaded_counts_dim=dim(counts_mat);

	cat("  Loaded: Num Samples: ", loaded_counts_dim[1], "\n");
	cat("  Loaded: Num Categories: ", loaded_counts_dim[2], "\n");
	cat("\n");

	message("  Loaded: Num Samples: ", loaded_counts_dim[1]);
	message("  Loaded: Num Categories: ", loaded_counts_dim[2]);
	message("\n");

	# Remove alls zero categories and counts
	counts_mat=remove_zero_count_categories_and_samples(counts_mat);
	counts_dim=dim(counts_mat);

	cat("  Returned: Num Samples: ", counts_dim[1], "\n");
	cat("  Returned: Num Categories: ", counts_dim[2], "\n");

	message("  Returned: Num Samples: ", counts_dim[1]);
	message("  Returned: Num Categories: ", counts_dim[2], "\n");


        return(counts_mat);
}

#------------------------------------------------------------------------------

load_reference_levels_file=function(fname){

	cat("Loading Reference Releveling File: ", fname, "\n", sep="");
	if(!specified(fname)){
		cat("  Referencing Releveling File not specified.\n");
		return(NULL);
	}

        inmat=as.matrix(read.table(fname, sep="\t", header=F, check.names=FALSE, 
		comment.char="#", row.names=1))

        colnames(inmat)=c("ReferenceLevel");

	cat("  Reference Levels:\n");
        print(inmat);
        cat("\n");

        if(ncol(inmat)!=1){
		msg=capture.output({
			cat("  Error reading in reference level file: ", fname, "\n");
			cat("  Should be a 2 column tab-separated text file.\n");
			cat("  Format: <variable name>\\t<reference level name>\\n");
			});
		stop(msg);
        }
        return(inmat);
}

relevel_factors=function(factors, ref_lev_mat){

	cat("Releveling References...\n");
        num_factors_to_relevel=nrow(ref_lev_mat);
	cat("  Num factors to relevel: ", num_factors_to_relevel, "\n", sep="");
        relevel_names=rownames(ref_lev_mat);
        factor_names=colnames(factors);

	shared=setdiff(relevel_names, factor_names);
	if(length(shared)==num_factors_to_relevel){
		cat("  All targeted variables to relevel have been found.\n");
	}

        for(i in 1:num_factors_to_relevel){
                relevel_target=relevel_names[i];

                if(length(intersect(relevel_target, factor_names))){

                        target_level=ref_lev_mat[i, 1];

                        factor_val=factors[,relevel_target];

                        if(length(intersect(target_level, factor_val))){

                                factor_val=relevel(factor_val, target_level);
                                factors[,relevel_target]=factor_val;
				cat("  Level to set as reference (", relevel_target,
					 "), found and set for ", target_level, ".\n", sep="");

                        }else{

                                cat("  WARNING: Target level '", target_level,
                                        "' not found in '", relevel_target, "'!!!\n", sep="");

				top_levels=sort(table(factor_val), decreasing=T);
				cat("  Top Levels: \n");
				print(top_levels);
				most_prevalent_level=names(top_levels)[1];
				cat("  Most prevalent level: ", most_prevalent_level, "\n", sep="");

                                factor_val=relevel(factor_val, most_prevalent_level);
                                factors[,relevel_target]=factor_val;
				cat("  Level to set as reference (", most_prevalent_level,
					 "), chosen and set for ", target_level, ".\n", sep="");
				
                        }

                }else{
                        cat("  WARNING: Relevel Target Not Found: '", relevel_target, "'!!!\n", sep="");
                }
        }
        return(factors);
}

#------------------------------------------------------------------------------

load_factors_file=function(fname, prim_key_cname=NULL, relevel_fn=NULL){

	cat("Loading Factors File: ", fname, "\n", sep="");
	if(!specified(fname)){
		cat("  Factor File not specified.\n");
		return(NULL);
	}

	if(!specified(prim_key_cname)){
		cat("  Primary Key Column Name not specified.  ");
        	factors=data.frame(read.table(fname,  sep="\t", header=TRUE,
                	check.names=FALSE, comment.char="#", stringsAsFactors=T));
	}else{
		cat("  Primary Key Column Name: ", prim_key_cname, "\n", sep="");
        	factors=data.frame(read.table(fname,  sep="\t", header=TRUE, row.names=prim_key_cname,
                	check.names=FALSE, comment.char="#", stringsAsFactors=T));
	}

	if(!specified(relevel_fn)){
		cat("  Factor Releveling File Not Specified.\n", sep="");
	}else{
		cat("  Factor Releveling File Specified: ", relevel_fn, "\n", sep="");
	}


	cat("\n");
	cat("Primary Key Excerpt: \n");
	cat(head(rownames(factors), 10));
	cat("\n\n");

        factor_names=colnames(factors);
	fact_dim=dim(factors);

	cat("  Num Samples/Rows: ", fact_dim[1], "\n", sep="");
	cat("  Num Variables: ", fact_dim[2], "\n", sep="");

	if(specified(relevel_fn)){
		relevel_mat=load_reference_levels_file(relevel_fn);	
		cat("Releveling Matrix:\n");
		print(relevel_mat);
		factors=relevel_factors(factors, relevel_mat);
	}

	return(factors);

}


#------------------------------------------------------------------------------

load_mapping=function(fname, pA, pB, sbj_id_cname=""){
        # This function will return a 2 column mapping, as requested by the pA and pB
        # names.  The rownames of the mapping will be based on the colname of the sbj_id_cname

	cat("Loading Pair Mapping File: ", fname, "\n", sep="");
	cat("  Requested Columns: A: ", pA, ", B: ", pB, "\n", sep="");

	if(!specified(fname)){
		cat("  Mapping file not specifed.\n");
		return(NULL);
	}

        mapping=read.table(fname, sep="\t", header=T, comment.char="#", quote="", row.names=NULL);

	mapping_dim=dim(mapping);
	cat("  Mapping Rows (subject IDs): ", mapping_dim[1], "\n");
	cat("  Mapping Cols (sample IDs): ", mapping_dim[2], "\n");

	mapping_cnames=colnames(mapping);
	cat("\n");
	cat("  Mapping Column Names:\n");
	print(mapping_cnames);
	cat("\n");	

        if(!specified(sbj_id_cname)){
                subject_ids=mapping[,1];
        }else{
		cat("  Subject ID Column Name Specified.\n");
                subject_ids=mapping[,sbj_id_cname];
        }

        column_names=colnames(mapping);
        if(all(column_names!=pA)){
                stop("  Error: Could not find ", pA, " in header of map file.\n");
        }
        if(all(column_names!=pB)){
                stop("  Error: Could not find ", pB, " in header of map file.\n");
        }

        map=cbind(as.character(mapping[,pA]), as.character(mapping[,pB]));
        colnames(map)=c(pA, pB);
        rownames(map)=subject_ids;

	cat("  Loaded Pairing Map:\n");
	print(map);

        # Remove pairings with NAs
	cat("\n");
	cat("  Removing any incomplete mappings...\n");
        incomp=apply(map, 1, function(x){any(is.na(x))});
        map=map[!incomp,];
	num_incomplete=sum(incomp);
	cat("  Number of incomplete pairings removed: ", num_incomplete, "\n");

	cat("\n");
	cat("Complete Map:\n");
	print(map);

	map_dim=dim(map);
	cat("\n");
	cat("  Paired Rows (subject IDs): ", map_dim[1], "\n");
	cat("  Paired Cols (sample IDs): ", map_dim[2], "\n");

        return(map);
}

#------------------------------------------------------------------------------

load_sbj_smp_mapping=function(fname, sbj_cname, smp_cname){

	cat("Loading Subject ID to Sample ID Mapping File: ", fname, "\n", sep="");

	if(!specified(fname)){
		cat("  Subject to Sample Mapping File not specified.\n");
		return(NULL);
	}

        map=read.table(fname, sep="\t", header=T, comment.char="#", quote="", row.names=sbj_cname);

	cat("  Number of rows read in: ", nrow(map), "\n");

	map=mapping[, smp_cname, drop=F];

	incomp=is.na(map[,1]);
        map=map[!incomp,1,drop=F];
	num_incomplete=sum(incomp);

	cat("  Number of incomplete mappings removed: ", num_incomplete, "\n");
	cat("  Number of mappings returned: ", nrow(map), "\n");

        return(map);
}

#------------------------------------------------------------------------------

intersect_pairings_map_by_keep_list=function(
	pairs_map, 
	sample_id_keepers=-1,
	subject_id_keepers=-1){

	cat("Intersecting pairings map with keep lists.\n");
	num_smp_keepers=length(sample_id_keepers);
	num_sbj_keepers=length(subject_id_keepers);

        num_rows=nrow(pairs_map);
        if(num_rows==0){
		cat("  Pairs map is empty.\n");
		return(pairs_map);
	}

	delete_row_ix=c();

	# By sample id
	if(num_smp_keepers>0 && sample_id_keepers[1] != -1){
		cat("  Num Sample IDs available in Keep List: ", num_smp_keepers, "\n", sep="");
		message("  Num Sample IDs available in Keep List: ", num_smp_keepers, sep="");
		for(rix in 1:num_rows){
			if(!any(pairs_map[rix, 1]==sample_id_keepers) ||
				!any(pairs_map[rix, 2]==sample_id_keepers)){
					delete_row_ix=c(delete_row_ix, rix);
			}
		}
	}else{
		cat("No Sample ID Keepers Specified.  No filter applied.\n");
	}

	# By subject id
	if(num_sbj_keepers>0 && subject_id_keepers[1] != -1){
		cat("  Num Subject IDs available in Keep List: ", num_sbj_keepers, "\n", sep="");
		message("  Num Subject IDs available in Keep List: ", num_sbj_keepers, sep="");
		map_sbj_ids=rownames(pairs_map);
		for(rix in 1:num_rows){
			if(!any(map_sbj_ids[rix]==subject_id_keepers)){
				delete_row_ix=c(delete_row_ix, rix);
			}
		}
	}else{
		cat("No Subject ID Keepers Specified. No filter applied.\n");
	}

	num_rows_to_delete=length(delete_row_ix);
	cat("Number of Rows to delete: ", num_rows_to_delete, "\n", sep="");

	if(num_rows_to_delete>0){
		cat("Delete Rows: \n");
		print(delete_row_ix);
	}
	
	complete_row_ix=setdiff(1:num_rows, delete_row_ix);

        incomplete_pairs=pairs_map[delete_row_ix,,drop=F];
        complete_pairs=pairs_map[complete_row_ix,,drop=F];

	#cat("  Complete Pairs: \n");
	#print(complete_pairs);

	if(num_rows_to_delete>0){
		cat("  Incomplete Pairs: \n");
		print(incomplete_pairs);
	}

	cat("\n");
	cat("  Num Input Pairs: ", nrow(pairs_map), "\n");
	cat("  Num Complete Pairs: ", nrow(complete_pairs), "\n");
	cat("  Num Incomplete Pairs: ", nrow(incomplete_pairs), "\n");

	message("");
	message("  Num Input Pairs: ", nrow(pairs_map));
	message("  Num Complete Pairs: ", nrow(complete_pairs));
	message("  Num Incomplete Pairs: ", nrow(incomplete_pairs));

        return(complete_pairs);
}

#------------------------------------------------------------------------------

split_goodbad_pairings_map=function(pairs_map){

        num_rows=nrow(pairs_map);
        keepers=apply(pairs_map, 1, function(x){ all(!is.na(x))});

        good_pairs_map=pairs_map[keepers,,drop=F];
        bad_pairs_map=pairs_map[!keepers,,drop=F];
        num_good_collapsed_rows=nrow(good_pairs_map);
        cat("Collapsed ", num_rows, " pairs to ", num_good_collapsed_rows, " complete pairs.\n");

        res=list();
        res[["good_pairs"]]=good_pairs_map;
        res[["bad_pairs"]]=bad_pairs_map;

        return(res);

}

#------------------------------------------------------------------------------

extract_top_categories=function(ordered_normalized, top, additional_cat=c()){

	cat("Extracting top ", top, " categories.\n");

        num_samples=nrow(ordered_normalized);
        num_categories=ncol(ordered_normalized);

        cat("  Samples: ", num_samples, "\n");
        cat("  Categories: ", num_categories, "\n");

        num_top_to_extract=min(num_categories-1, top);

        cat("  Num Top Requested to Extract: ", top, "\n");
        message("  Num Top Requested to Extract: ", top);
        cat("  Columns Available for Extraction: ", num_top_to_extract, "\n");
        message("  Columns Available for Extraction: ", num_top_to_extract);

        # Extract top categories requested
        top_cat=ordered_normalized[,1:num_top_to_extract, drop=F];

        if(length(additional_cat)){
                cat("  Additional Categories to Include:\n");
                print(additional_cat);
        }else{
                cat("  No Additional Categories to Extract.\n");
        }

        # Extract additional categories
        # :: Make sure we can find the categories
        available_cat=colnames(ordered_normalized);
        missing_cat=setdiff(additional_cat, available_cat);
        if(length(missing_cat)){
		msg=capture.output({
			cat("  Error: Could not find categories: \n");
			print(missing_cat);
			});
		stop(msg);
        }

        # :: Remove categories we have already extracted in the top N
        already_extracted_cat=colnames(top_cat);
        extra_cat=setdiff(additional_cat, already_extracted_cat);

        num_extra_to_extract=length(extra_cat);
        cat("  Num Extra Categories to Extract: ", num_extra_to_extract, "\n");
        message("  Num Extra Categories to Extract: ", num_extra_to_extract);

        # Allocate/Prepare output matrix
        num_out_mat_cols=num_top_to_extract+num_extra_to_extract+1;
        out_mat=matrix(0, nrow=num_samples, ncol=num_out_mat_cols);
        rownames(out_mat)=rownames(ordered_normalized);
        colnames(out_mat)=c(already_extracted_cat, extra_cat, "Remaining");

        # Copy over top and additional categories, and compute remainding
        all_cat_names=c(already_extracted_cat, extra_cat);
        out_mat[,all_cat_names]=ordered_normalized[,all_cat_names];
        out_mat[,"Remaining"]=apply(out_mat, 1, function(x){1-sum(x)});

        return(out_mat);

}

extract_categories_for_counts=function(counts_matrix, cat_to_extr){

	# Extracts counts based on extract list, and computes
	# Remaining counts.

	cat_to_extr=setdiff(cat_to_extr, "Remaining");
	totals=apply(counts_matrix, 1, sum);
	extracted_mat=counts_matrix[,cat_to_extr,drop=F];
	extracted_totals=apply(extracted_mat, 1, sum);
	remaining_counts=totals-extracted_totals;
	extr_wrem=cbind(extracted_mat, remaining_counts);
	colnames(extr_wrem)=c(cat_to_extr, "Remaining");
	return(extr_wrem);
	
}


###############################################################################
###############################################################################

check_variables=function(covar, grp, req, avail){

	num_covar=length(covar);
	num_grp=length(grp);
	num_req=length(req);
	num_avail=length(avail);
	
	cat("Checking Variable Lists:\n");
	cat("  Num Covariates: ", num_covar, "\n");
	cat("  Num Group Var: ", num_grp, "\n");
	cat("  Num Required Var: ", num_req, "\n");
	cat("  Num Available Var (Factors): ", num_avail, "\n");

	# Checking for overlap between covariates and group lists
	overlap=intersect(covar, grp);
	if(length(overlap)==0){
		cat("  Good.  No overlap between Covariate and Group lists.\n");
	}else{
		msg=capture.output({
			cat("  Error:  Overlap detected between Covariates and Group Lists:\n");
			print(overlap);
			});
		stop(msg);
	}

	# Checking required are subset of covariates or group list
	combined=c(covar, grp);
	overlap=intersect(combined, req);
	if(length(overlap)==num_req){
		cat("  Good.  Required variables overlap with Covariates or Group List.\n");
	}else{
		msg=capture.output({
			cat("  Error: Required variables do not overlap with Covariates or Group List:\n");
			cat("  Covariates:\n");
			print(covar);
			cat("  Group:\n");
			print(grp);
			cat("  Required:\n");
			print(req);
			});
		stop(msg);
	}

	# Checking group and covariates are in available list
	overlap=intersect(combined, avail);
	if(length(overlap)==length(combined)){
		cat("  Good.  All variables (Covariates and Group) found/available.\n");
	}else{
		msg=capture.output({
			cat("  Error.  Some variables missing from available requested as Covariates/Group:\n");
			print(setdiff(combined, avail));
			});
		stop(msg);
	}
	
	return(1);
		
}

screen_variables=function(covar, grp, req, avail){

	num_covar=length(covar);
	num_grp=length(grp);
	num_req=length(req);
	num_avail=length(avail);

	cat("Screening Variable Lists:\n");
	cat("Inputs:\n");
	cat("  Num Covariates: ", num_covar, "\n");
	cat("  Num Group Var: ", num_grp, "\n");
	cat("  Num Required Var: ", num_req, "\n");
	cat("  Num Available Var (Factors): ", num_avail, "\n");
	cat("\n");

	fcovar=intersect(avail, covar);
	fgrp=intersect(avail, grp);
	freq=intersect(avail, req);

	num_fcovar=length(fcovar);
	num_fgrp=length(fgrp);
	num_freq=length(freq);

	if(num_fcovar!=num_covar){
		cat("Covariate list adjusted.  Removed:\n");
		print(setdiff(covar, fcovar));
	}
	if(num_fgrp!=num_grp){
		cat("Group Var list adjusted.  Removed:\n");
		print(setdiff(grp, fgrp));
	}
	if(num_freq!=num_req){
		cat("Required Var list adjusted.  Removed:\n");
		print(setdiff(req, freq));
	}

	cat("Outputs:\n");
	cat("  Num Covariates: ", num_fcovar, "\n");
	cat("  Num Group Var: ", num_fgrp, "\n");
	cat("  Num Required Var: ", num_freq, "\n");
	cat("  Num Available Var (Factors): ", num_avail, "\n");

	screened_var=list();
	screened_var[["covariates"]]=fcovar;
	screened_var[["group_var"]]=fgrp;
	screened_var[["required_var"]]=freq;

	return(screened_var);
		
}

# Reconcile Functions:

match_to_index=function(values_arr, targets_arr){
	# This function will make a boolean index based on matching
	# the values that are in the targets
	
	matches=sapply(values_arr, function(x){
		return(any(x==targets_arr));
		});

	return(matches);
}

recon_factor_to_sumtab=function(fact_mat, st_mat, fact_samp_id_cname){

	# Reconcile by Sample IDs

	st_samp_ids=rownames(st_mat);
	fact_samp_ids=as.character(fact_mat[, fact_samp_id_cname]);
	shared_samp_ids=intersect(st_samp_ids, fact_samp_ids);
	
	match_ix=match_to_index(fact_samp_ids, shared_samp_ids);	

	out_factor_mat=fact_mat[match_ix,,drop=F];
	
	out_summary_table_mat=st_mat[shared_samp_ids,,drop=F];	

	results=list();
	results[["factor_mat"]]=out_factor_mat;
	results[["summary_table"]]=out_summary_table_mat;

	return(results);

}



recon_factor_to_pairs_to_sumtab=function(fact_mat, pr_mat, st_mat, f_subj_id_cname){

	# Factor <-> Pairs (By Subject ID)
	fact_subj_ids=fact_mat[,f_subj_id_cname];
	pairs_subj_ids=rownames(pr_mat);
	shared_subj_ids=intersect(fact_subj_ids, pairs_subj_ids);

	#  keep shared subject IDs
	out_pr_mat=pr_mat[shared_subj_ids,,drop=F];
	fact_shared_subj_ix=match_to_index(fact_mat[,f_subj_id_cname], shared_subj_ids);
	out_fact_mat=fact_mat[fact_shared_subj_ix,,drop=F];

	# Pairs <-> Summary Table (By Sample ID)
	st_sample_ids=rownames(st_mat);
	out_pr_mat=intersect_pairings_map_by_keep_list(out_pr_mat, sample_id_keepers=st_sample_ids);
	shared_samp_ids=intersect(c(out_pr_mat[,1],out_pr_mat[,2]), st_sample_ids);
	out_st_mat=st_mat[shared_samp_ids,,drop=F];


	# Pairs <-> Factor (By Subject ID)
	pr_subject_ids=rownames(out_pr_mat);
	fact_shared_subj_ix=match_to_index(out_fact_mat[,f_subj_id_cname], pr_subject_ids);
	out_fact_mat=out_fact_mat[fact_shared_subj_ix,,drop=F];	 


	results=list();
	results[["factor_mat"]]=out_fact_mat;
	results[["pairs_mat"]]=out_pr_mat;
	results[["summary_table"]]=out_st_mat;

	return(results);

}



reconcile=function(param=list("summary_table_mat"=NULL, "factor_mat"=NULL, "pairs_mat"=NULL,
	"factor_sbj_cname"=NULL, "factor_smp_cname"=NULL)){	

	# Summary Table Mat and Factor Mat must be specified.
	# The Factor Mat must have a subject ID and sample ID column specified.
	# The Pairs Mat is Optional, but the first column must be the subject ID.

	# Keying:
	#  Summary Table is Keyed by Sample ID.
	#  Factor File is not keyed, it may have duplicated Subject IDs
	#  Pairs Mat is keyed by Subject ID.

	summary_table_mat=param[["summary_table_mat"]];
	factor_mat=param[["factor_mat"]];
	pairs_mat=param[["pairs_mat"]];

	# Get Smp/Sbj from the available data
	# Factor Smp/Sbj IDs
	factor_subj_id_cname=param[["factor_sbj_cname"]];
	factor_samp_id_cname=param[["factor_smp_cname"]];
	# Summary Table Smp IDs
	sumtab_sample_ids=rownames(summary_table_mat);

	message("----------------------------------------------------------");
	message("Files Specifed:");
	message(" Sum Tab File: ", specified(summary_table_mat));
	message(" Factor File: ", specified(factor_mat));
	message(" Pairs File: ", specified(pairs_mat));
	message("----------------------------------------------------------");
	message();


	if(!specified(pairs_mat)){
		message("Reconcile Summary Table/Factors by shared sample IDs.");
		reconned=recon_factor_to_sumtab(
			factor_mat, summary_table_mat, factor_samp_id_cname);

		summary_table_mat=reconned[["summary_table"]];
		factor_mat=reconned[["factor_mat"]];
	}else{
		message("Reconcile Factors <-> Pairs <-> SummaryTable")
		reconned=recon_factor_to_pairs_to_sumtab(
			factor_mat, pairs_mat, summary_table_mat, 
			factor_subj_id_cname);

		summary_table_mat=reconned[["summary_table"]];
		factor_mat=reconned[["factor_mat"]];
		pairs_mat=reconned[["pairs_mat"]];
	}

	results=list();
	results[["summary_table_mat"]]=summary_table_mat;
	results[["factor_mat"]]=factor_mat;
	results[["pairs_mat"]]=pairs_mat;
	results[["factor_sbj_cname"]]=param[["factor_sbj_cname"]]; 
	results[["factor_smp_cname"]]=param[["factor_smp_cname"]];

	# print_se(results);

	return(results);
}

load_and_reconcile_files=function(
	sumtab=list(fn=NULL, shorten_cat_names_char=NULL, return_top=NULL, 
		specific_cat_fn=NULL), 
	factors=list(fn=NULL, sbj_cname=NULL, smp_cname=NULL, ref_relvl_fn=NULL), 
	pairs=list(fn=NULL, a_cname=NULL, b_cname=NULL, subject_id_cname=NULL), 
	sbj_to_smp=list(fn=NULL, sbjid_cname=NULL, smpid_cname=NULL),
	covariates=list(fn=NULL), 
	grpvar=list(fn=NULL), 
	reqvar=list(fn=NULL)){

	message("Loading and Reconciling Files:\n");
	
	# 1.) Read in Summary Table, Factors, Sample Pairing
	# 2.) Keep/Subset target variables
	# 3.) Reconcile (remove subjects with missing samples)
	# 4.) Remove NAs (remove subjects/variables with NAs)
	# 5.) Reconcile (remove subjects with missing metadata)
	# 6.) Apply summary table parameters
	# 7.) Clean up category names

	report_list=list();

	#-----------------------------------------------------------------------------
	# 1.) Read in Summary Table, Factors, Sample Pairing

	message("\n***************************************************************");
	message("*** 1/6 Reading Summary Table, Factors, and Sample Pairing. ***");
	message("***************************************************************\n");

	message("Loading Summary Table: ", sumtab[["fn"]], "\n", sep="");
	report_list[["Summary Table"]]=capture.output({

		summary_table_mat=load_summary_file(fname=sumtab[["fn"]]);

		cat("Summary Table Dimensions: ", dim(summary_table_mat), "\n");
		message("Summary Table Dimensions: ", dim(summary_table_mat));
	});

	message("Loading Factors File: ", factors[["fn"]], "\n", sep="");
	report_list[["Factor File"]]=capture.output({

		factors_mat=load_factors_file(
			fname=factors[["fn"]],
			relevel_fn=factors[["ref_relvl_fn"]]
			);

		cat("Factor Table Dimensions: ", dim(factors_mat), "\n");
		message("   Factor Table Dimensions: ", dim(factors_mat));

	});


	message("Loading Mappings: ", pairs[["fn"]], " (pairs) and ", 
		sbj_to_smp[["fn"]], " (sbj_to_smp)\n", sep="");

	report_list[["Mappings"]]=capture.output({

		pairs_mat=load_mapping(
			fname=pairs[["fn"]], 
			pA=pairs[["a_cname"]], 
			pB=pairs[["b_cname"]], 
			sbj_id_cname=pairs[["subject_id_cname"]]
			);

		cat("Pairs Mapping Dimensions: ", dim(pairs_mat), "\n");
		message("  Pairs Mapping Dimensions: ", dim(pairs_mat));
		cat("\n");

		sbj_smp_mat=load_sbj_smp_mapping(
			fname=sbj_to_smp[["fn"]],
			sbj_cname=sbj_to_smp[["sbjid_cname"]],
			smp_cname=sbj_to_smp[["smpid_cname"]]
			);

		cat("Subject-Sample Mapping Dimensions: ", dim(sbj_smp_mat), "\n");
		message("  Subject-Sample Mapping Dimensions: ", dim(sbj_smp_mat));

	});

	#-----------------------------------------------------------------------------
	# 2.) Load and Keep/Subset target variables

	message("\n***************************************************************");
	message("*** 2/6 Loading lists.                                      ***");
	message("***************************************************************\n");

	message("\t", covariates[["fn"]], " (Covariates)\n", sep="");
	message("\t", grpvar[["fn"]], " (Group)\n", sep="");
	message("\t", reqvar[["fn"]], " (Required)\n", sep="");


	report_list[["Variable Lists"]]=capture.output({

		covariates_arr=load_list(fname=covariates[["fn"]], 
			memo="Covariates File");

		groupvar_arr=load_list(fname=grpvar[["fn"]], 
			memo="Group Variables File");

		requiredvar_arr=load_list(fname=reqvar[["fn"]], 
			memo="Required Variables File");

		cat("\n");
		
		if(length(covariates_arr)==0 && length(groupvar_arr)==0){
			cat("No covariates or group variables specified.\n");
			cat("  Continuing with all variables.\n");
			factors_subset_mat=factors_mat;
		}else{

			status=check_variables(
				covariates_arr, groupvar_arr, requiredvar_arr, 
				colnames(factors_mat));

			subset_col_arr=c();
			if(specified(factors[["sbj_cname"]])){
				subset_col_arr=c(subset_col_arr, factors[["sbj_cname"]], NULL);
			}

			if(specified(factors[["smp_cname"]])){
				subset_col_arr=c(subset_col_arr, factors[["smp_cname"]], NULL);
			}

			subset_col_arr=c(subset_col_arr, covariates_arr, groupvar_arr);

			cat("Subsetting requested variables from factors.\n");
			message("Subsetting requested variables from factors.");

			missing=setdiff(subset_col_arr, colnames(factors_mat));
			if(length(missing)>0){
				message("Missing variables from Factors:");
				print_se(missing);
			}

			factors_subset_mat=
				factors_mat[, subset_col_arr, drop=F];

		}

		factor_subset_dim=dim(factors_subset_mat);
		cat("  Factor Subset Rows/Subjects: ", factor_subset_dim[1], "\n");
		cat("  Factor Subset Cols/Variables: ", factor_subset_dim[2], "\n");

	});

	#-----------------------------------------------------------------------------
	# 3.) Reconcile

	message("\n***************************************************************");
	message("*** 3/6 Performing Pre-NA Removal Reconcile.                ***");
	message("***************************************************************\n");
	report_list[["Pre-NA Removal Reconcile"]]=capture.output({

		reconciled_files=reconcile(param=list(
			summary_table_mat=summary_table_mat, 
			factor_mat=factors_subset_mat, 
			pairs_mat=pairs_mat, 
			sbj_smp_mat=sbj_smp_mat,
			factor_sbj_cname=factors[["sbj_cname"]],
			factor_smp_cname=factors[["smp_cname"]]
			));

		recon_factors=reconciled_files[["factor_mat"]];

	});

	#-----------------------------------------------------------------------------
	# 4.) Remove NAs (remove subjects with NAs)

	message("\n***************************************************************");
	message("*** 4/6 Remove NAs from Factors (if necessary).            ***");
	message("***************************************************************\n");

	if(factors[["fn"]]!=""){
		message("Performing NA Removal.\n");
		report_list[["NA Removal"]]=capture.output({

			nona_factors_fn=paste(gsub("\\.tsv", "", factors[["fn"]]), ".noNAs", sep="");

			factors_wo_nas_res=remove_sample_or_factors_wNA_parallel(
				recon_factors,
				required=requiredvar_arr,
				outfile=nona_factors_fn,
				num_trials=64000, num_cores=64);

			factors_mat_nona=factors_wo_nas_res[["factors"]];

			# remove variables with no information
			variables_woNAs=colnames(factors_mat_nona);

			cat("\n");

			screened_var=screen_variables(
				covariates_arr,
				groupvar_arr,
				requiredvar_arr,
				variables_woNAs);

			covariates_arr=screened_var[["covariates"]];
			groupvar_arr=screened_var[["group_var"]];
			requiredvar_arr=screened_var[["required_var"]];
		});
	}else{
		report_list[["NA Removal"]]="Factor File not specified.  NA Removal Skipped.";
		factors_mat_nona=NULL;
	}

	#-----------------------------------------------------------------------------
	# 5.) Reconcile
	
	message("\n***************************************************************");
	message("*** 5/6 Reconciling after Factor NA remove (if necessary).  ***");
	message("***************************************************************\n");
	if(factors[["fn"]]!=""){
		message("Performing Post-NA Removal Reconcile.\n");
		report_list[["Post-NA Removal Reconcile"]]=capture.output({

			# Replace factor mat, but keep the other data/matrices the same
			reconciled_files[["factor_mat"]]=factors_mat_nona;
			reconciled_files=reconcile(reconciled_files);

		});
	}else{
		report_list[["Post-NA Removal Reconcile"]]=
			"Post-NA Removal Reconcile not performed.";
	}

	#-----------------------------------------------------------------------------
	# 6.) Apply Summary Table Parameters

	message("\n***************************************************************");
	message("*** 6/6 Apply Summary Table Parameters.                     ***");
	message("***************************************************************\n");

	message("Applying Summary Table Parameters.\n");
	report_list[["Summary Table Options"]]=capture.output({

		# Order categories by decreasing abundance
		cat("Reordering summary table counts by decreasing abundance.\n");
		counts_mat=reconciled_files[["summary_table_mat"]];

		normalized_mat=normalize(counts_mat);
		normalized_mat=reorder_by_decreasing_abundance(normalized_mat);
		counts_mat=counts_mat[,colnames(normalized_mat), drop=F];

		cat("\n");

		if(specified(sumtab[["return_top"]])){
			top_to_extract=as.numeric(sumtab[["return_top"]]);
			cat("Num Top Categories to extract: ", top_to_extract, "\n");
		}else{
			cat("Top Categories for Extraction not specified.\n");
			top_to_extract=NA;
		}

		spec_cat_arr=c();
		if(specified(sumtab[["specific_cat_fn"]])){
			cat("Specific Categories requested for extraction.\n");
			spec_cat_arr=load_list(
				fname=sumtab[["specific_cat_fn"]], 
				memo="Specific Categories File");
			cat("Num Targets: ", length(spec_cat_arr), "\n");
			cat("Targets: \n");
			print(spec_cat_arr);
		}else{
			cat("Specific categories not specified.\n");
		}

		cat("\n");

		if(is.na(top_to_extract) && length(spec_cat_arr)==0){
			cat("Top and/or Specific Summary Table categories not requested.\n");
			cat("  Returning all categories.\n");
		}else{
			if(length(spec_cat_arr)>0 && is.na(top_to_extract)){
				cat("Specific categories requested w/o accompanying top categories.\n");
				top_to_extract=0;
			}
	
			extracted_normalized=extract_top_categories(
				normalized_mat,
				top_to_extract,
				additional_cat=spec_cat_arr);

			normalized_mat=extracted_normalized;
			counts_mat=extract_categories_for_counts(counts_mat, colnames(normalized_mat));
		}


		norm_dim=dim(normalized_mat);
		cat("\n");
		cat("Output Dimensions:\n");
		cat("  Samples: ", norm_dim[1], "\n", sep="");
		cat("  Categories: ", norm_dim[2], "\n", sep="");
		cat("\n");

		#-----------------------------------------------------------------------------
		# Shorten categories
		if(specified(sumtab[["shorten_cat_names_char"]])){

			cat("Shortening category names specified. Delimiter: '", 
				sumtab[["shorten_cat_names_char"]], "'\n", sep="");

			colnames(normalized_mat)=
				shorten_category_names(colnames(normalized_mat), 
				sumtab[["shorten_cat_names_char"]]);

			colnames(counts_mat)=colnames(normalized_mat);
		}else{
			cat("Shortening category names not specified.\n");
		}

		# Clean category names a little
		cat("Cleaning Category names.\n");
		colnames(normalized_mat)=clean_category_names(colnames(normalized_mat));
		colnames(counts_mat)=colnames(normalized_mat);

	});

	#-----------------------------------------------------------------------------
	message("\n***************************************************************");
	message("*** Completing Load and Reconciliation, Returning Results.  ***");
	message("***************************************************************\n");
	
	results=list();
	results[["SummaryTable_counts"]]=counts_mat;
	results[["SummaryTable_normalized"]]=normalized_mat;
	results[["Factors"]]=factors_mat_nona;
	results[["PairsMap"]]=reconciled_files[["pairs_mat"]];
	results[["SubjectToSampleIDMap"]]=reconciled_files[["sbj_to_smp_id_mat"]];
	results[["Covariates"]]=covariates_arr;
	results[["GroupVariables"]]=groupvar_arr;
	results[["RequiredVariables"]]=requiredvar_arr;
	results[["Report"]]=report_list;

	return(results);
}

###############################################################################

write_file_report=function(report_list){
	# This will output a report to a pdf file.
	# To call use the return value of "load_and_reconcile_files"
	# E.g. write_file_report(loaded_files[["Report"]]);

	report_groups=names(report_list);

        par(family="Courier");
        par(oma=rep(.1,4));
        par(mar=rep(0,4));

	for(grp in report_groups){
	
		strings=report_list[[grp]];

		num_lines=length(strings);

		top=max(as.integer(num_lines), 52);

		par(mfcol=c(1,1));

		plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
			xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
			);

		text(0, top, grp, cex=1, font=2, pos=4);

		text_size=max(.01, min(.8, .8 - .003*(num_lines-52)));

		for(i in 1:num_lines){
			strings[i]=gsub("\t", "", strings[i]);
			text(0, top-i-1, strings[i], pos=4, cex=text_size);
		}

	}

}


