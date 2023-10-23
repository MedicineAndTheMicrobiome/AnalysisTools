
# This is a collection of shared functions that make the paired analyses
# consistent in how the metadata and samples are reconciled.

source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');

#------------------------------------------------------------------------------

load_list=function(filename){
        val=scan(filename, what=character(), comment.char="#");
        return(val);
}


#------------------------------------------------------------------------------

remove_zero_count_categories_and_samples=function(mat){

	# Remove zero count samples
	cat("Checking for samples with no counts.\n");
	tot=apply(mat, 1, sum);
	nonzero=tot>0;
	if(!(all(nonzero))){
		cat("WARNING: Zero count samples found:\n");
		samp_names=rownames(mat);
		print(samp_names[!nonzero]);
		cat("\n");
		mat=mat[nonzero,,drop=F];
	}

	# Remove zero count categories 
	cat("Checking for categories with no counts.\n");
	tot=apply(mat, 2, sum);
	nonzero=tot>0;
	if(!(all(nonzero))){
		cat("WARNING: Zero count categories found:\n");
		cat_names=colnames(mat);
		print(cat_names[!nonzero]);
		cat("\n");
		mat=mat[,nonzero,drop=F];
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
        }
	return(short_names);
}

#------------------------------------------------------------------------------

clean_category_names=function(catnames){

        cat_names=gsub("-", "_", cat_names);
	cat_names=gsub("\\[", "", cat_names);
	cat_names=gsub("\\]", "", cat_names);

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

	cat("Loaded: Num Samples: ", loaded_counts_dim[1], "\n");
	cat("Loaded: Num Categories: ", loaded_counts_dim[2], "\n");

	# Remove alls zero categories and counts
	counts_mat=remove_zero_count_categories_and_samples(counts_mat);
	counts_dim=dim(counts_mat);

	cat("Returned: Num Samples: ", counts_dim[1], "\n");
	cat("Returned: Num Categories: ", counts_dim[2], "\n");

        return(counts_mat);
}

#------------------------------------------------------------------------------

load_reference_levels_file=function(fname){
        inmat=as.matrix(read.table(fname, sep="\t", header=F, check.names=FALSE, comment.char="#", row.names=1))
        colnames(inmat)=c("ReferenceLevel");
        print(inmat);
        cat("\n");
        if(ncol(inmat)!=1){
                cat("Error reading in reference level file: ", fname, "\n");
                quit(status=-1);
        }
        return(inmat);
}

relevel_factors=function(factors, ref_lev_mat){

        num_factors_to_relevel=nrow(ref_lev_mat);
        relevel_names=rownames(ref_lev_mat);
        factor_names=colnames(factors);

        for(i in 1:num_factors_to_relevel){
                relevel_target=relevel_names[i];

                if(length(intersect(relevel_target, factor_names))){
                        target_level=ref_lev_mat[i, 1];
                        tmp=factors[,relevel_target];
                        if(length(intersect(target_level, tmp))){
                                tmp=relevel(tmp, target_level);
                                factors[,relevel_target]=tmp;
                        }else{
                                cat("WARNING: Target level '", target_level,
                                        "' not found in '", relevel_target, "'!!!\n", sep="");
                        }
                }else{
                        cat("WARNING: Relevel Target Not Found: '", relevel_target, "'!!!\n", sep="");
                }
        }
        return(factors);
}

#------------------------------------------------------------------------------

load_factors=function(fname, samp_id_colname=1, relevel_fn=NULL){

        factors=data.frame(read.table(fname,  sep="\t", header=TRUE, row.names=samp_id_colname,
                check.names=FALSE, comment.char="#", stringsAsFactors=T));
        factor_names=colnames(factors);

        ignore_idx=grep("^IGNORE\\.", factor_names);

        if(length(ignore_idx)!=0){
                return(factors[-ignore_idx]);
        }else{
                return(factors);
        }

	if(!is.null(relevel_fn)){
		relevel_mat=load_reference_levels_file(relevel_fn);	
		factors=relevel_factors(factors, relevel_mat);
	}

	return(factors);

}


#------------------------------------------------------------------------------

load_mapping=function(filename, pA, pB, sbj_id_cname=""){
        # This function will return a 2 column mapping, as requested by the pA and pB
        # names.  The rownames of the mapping will be based on the colname of the sbj_id_cname

        mapping=read.table(filename, sep="\t", header=T, comment.char="#", quote="", row.names=NULL);

        if(sbj_id_cname==""){
                subject_ids=mapping[,1];
        }else{
                subject_ids=mapping[,sbj_id_cname];
        }

        column_names=colnames(mapping);
        if(all(column_names!=pA)){
                cat("Error: Could not find ", pA, " in header of map file.\n");
                quit(status=-1);
        }
        if(all(column_names!=pB)){
                cat("Error: Could not find ", pB, " in header of map file.\n");
                quit(status=-1);
        }

        map=cbind(as.character(mapping[,pA]), as.character(mapping[,pB]));
        colnames(map)=c(pA, pB);
        rownames(map)=subject_ids;

        # Remove pairings with NAs
        incomp=apply(map, 1, function(x){any(is.na(x))});
        map=map[!incomp,];

        return(map);
}

#------------------------------------------------------------------------------

intersect_pairings_map_by_sample_id=function(pairs_map, keepers){

        missing=character();
        # Sets mappings to NA if they don't exist in the keepers array
        num_rows=nrow(pairs_map);
        if(num_rows>0){
                for(rix in 1:num_rows){
                        if(!any(pairs_map[rix, 1]==keepers) && !any(pairs_map[rix, 2]==keepers)){
                                missing=rbind(missing, pairs_map[rix, c(1,2)]);
                                pairs_map[rix, c(1,2)]=NA;
                        }
                }
        }

        results=list();
        results[["pairs"]]=pairs_map;
        results[["missing"]]=missing;
        return(results);
}

#------------------------------------------------------------------------------

intersect_pairings_map_by_subject_id=function(pairs_map, keepers){

        missing=character();

        # Sets mappings to NA if they don't exist in the keepers array
        num_rows=nrow(pairs_map);
	pairs_map_sbj_id=rownames(pairs_map);
        if(num_rows>0){
                for(rix in 1:num_rows){
			if(!any(pairs_map_sbj_id[rix]==keepers)){
                                missing=rbind(missing, pairs_map[rix, c(1,2)]);
                                pairs_map[rix, c(1,2)]=NA;
                        }
                }
        }

        results=list();
        results[["pairs"]]=pairs_map;
        results[["missing"]]=missing;
        return(results);
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

        num_samples=nrow(ordered_normalized);
        num_categories=ncol(ordered_normalized);

        cat("Samples: ", num_samples, "\n");
        cat("Categories: ", num_categories, "\n");

        num_top_to_extract=min(num_categories-1, top);

        cat("Top Requested to Extract: ", top, "\n");
        cat("Columns to Extract: ", num_top_to_extract, "\n");

        # Extract top categories requested
        top_cat=ordered_normalized[,1:num_top_to_extract, drop=F];

        if(length(additional_cat)){
                cat("Additional Categories to Include:\n");
                print(additional_cat);
        }else{
                cat("No Additional Categories to Extract.\n");
        }

        # Extract additional categories
        # :: Make sure we can find the categories
        available_cat=colnames(ordered_normalized);
        missing_cat=setdiff(additional_cat, available_cat);
        if(length(missing_cat)){
                cat("Error: Could not find categories: \n");
                print(missing_cat);
                quit(status=-1);
        }

        # :: Remove categories we have already extracted in the top N
        already_extracted_cat=colnames(top_cat);
        extra_cat=setdiff(additional_cat, already_extracted_cat);

        num_extra_to_extract=length(extra_cat);
        cat("Num Extra Categories to Extract: ", num_extra_to_extract, "\n");

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


###############################################################################
###############################################################################

check_variables=function(covar, grp, req, avail){

	num_covar=length(covar);
	num_grp=length(grp);
	num_req=length(req);
	num_avail=length(avail);
	
	cat("Checking Variable Lists:\n");
	cat("Num Covariates: ", num_covar, "\n");
	cat("Num Group Var: ", num_grp, "\n");
	cat("Num Required Var: ", num_req, "\n");
	cat("Num Available Var (Factors): ", num_avail, "\n");

	# Checking for overlap between covariates and group lists
	overlap=intersect(covar, grp);
	if(length(overlap)==0){
		cat("Good.  No overlap between Covariate and Group lists.\n");
	}else{
		cat("Error:  Overlap detected between Covariates and Group Lists:\n");
		print(overlap);
		quit(-1);
	}

	# Checking required are subset of covariates or group list
	combined=c(covar, grp);
	overlap=intersect(combined, req);
	if(length(overlap)==num_req){
		cat("Good.  Required variables overlap with Covariates or Group List:\n");
	}else{
		cat("Error: Required variables do not overlap with Covariates or Group List:\n");
		cat("Covariates:\n");
		print(covar);
		cat("Group:\n");
		print(grp);
		cat("Required:\n");
		print(req);
		quit(-1);
	}

	# Checking group and covariates are in available list
	overlap=intersect(combined, avail);
	if(length(overlap)==length(combined)){
		cat("Good.  All variables (Covariates and Group) found in available.\n");
	}else{
		cat("Error.  Some variables missing from available requested as Covariates/Group:\n");
		print(overlap);
		quit(-1);
	}
	
	return(1);
		
}


reconcile=function(param=list("summary_table_mat"=NULL, "factor_mat"=NULL, "pairs_mat"=NULL)){
	
	# Get sample IDs from summary table
	sumtab_sample_ids=rownames(summary_table_mat);

	# From the pairings, only keep pairs that have both samples available
	intersect_results=intersect_pairings_map_by_sample_id(pairs_mat, sumtab_sample_ids);

	# Get subject IDs from factors
	factor_subject_ids=rownames(factor_mat);

	# From the pairings, only keep pairs with subject IDs	
	pairs_mat=intersect_results[["good_pairs"]]
	intersect_results=intersect_pairings_map_by_subject_id(pairs_mat, factor_subject_ids);

	# Split out only the good pairs
	pairs_mat=intersect_results[["good_pairs"]]
	split_goodbad_pairings_map(pairs_mat);

	# Final subject_id list
	shared_sbj_ids=sort(rownames(pairs_mat));
	shared_smp_ids=sort(c(pairs_mat[,1], pairs_mat[,2]));

	results=list();
	results[["summary_table_mat"]]=summary_table_mat[shared_smp_ids,,drop=F];
	results[["factor_mat"]]=factor_mat[shared_sbj_ids,,drop=F];
	results[["pairs_mat"]]=pairs_mat[shared_sbj_ids,,drop=F];

	return(results);
}

load_and_reconcile_files=function(
	sumtab=list(fn=NULL, shorten_char=NULL, return_top=NULL, specific_cat_fn=NULL), 
	factors=list(fn=NULL, sbj_cname=NULL, ref_relvl_fn=NULL), 
	pairs=list(fn=NULL, a_cname=NULL, b_cname=NULL, subject_id_cname=NULL), 
	covariates=list(fn=NULL), 
	grpvar=list(fn=NULL), 
	reqvar=list(fn=NULL)){

	cat("Loading and Reconciling Files:\n");
	
	# 1.) Read in Summary Table, Factors, Sample Pairing
	# 2.) Keep/Subset target variables
	# 3.) Reconcile (remove subjects with missing samples)
	# 4.) Remove NAs (remove subjects/variables with NAs)
	# 5.) Reconcile (remove subjects with missing metadata)
	# 6.) Apply summary table parameters
	# 7.) Clean up category names

	#-----------------------------------------------------------------------------
	# 1.) Read in Summary Table, Factors, Sample Pairing

	summary_table_mat=load_summary_file(
		fn=sumtab[["fn"]], 
		shorten_cat_names_char=sumtab[["shorten_char"]],
		);

	load_list(sumtabl[["specific_cat_fn"]]);

	factors_mat=load_factors_file(
		fn=factors[["fn"]],
		subject_column_name=factors[["sbj_cname"]],
		ref_lvl_fn=factors[["ref_relvl_fn"]]
		);

	pairs_mat=load_mapping(
		filename=pairs[["fn"]], 
		pA=pairs[["a_cname"]], 
		pB=pairs[["b_cname"]], 
		sbj_id_cname=pairs[["subject_id_cname"]]
		);

	covariates_arr=load_list(filename=covariates[["fn"]]);

	groupvar_arr=load_list(filename=grpvar[["fn"]]);

	requiredvar_arr=load_list(filename=reqvar[["fn"]]);

	#-----------------------------------------------------------------------------
	# 2.) Keep/Subset target variables

	status=check_variables(covariates_arr, groupvar_arr, requiredvar_arr, colnames(factors_mat));
	
	cat("Subsetting requested variables from factors.\n");
	factors_subset_mat=factors_mat[, c(covariates_arr, groupvar_arr), drop=F];

	#-----------------------------------------------------------------------------
	# 3.) Reconcile

	reconcile(summary_table_mat, factors_mat, pairs_mat);

	#-----------------------------------------------------------------------------
	# 4.) Remove NAs (remove subjects with NAs)

	factors_wo_nas_res=remove_sample_or_factors_wNA_parallel(recon_factors,
        required=required_arr, num_trials=64000, num_cores=64, outfile=paste(OutputRoot, ".noNAs", sep=""));

	#-----------------------------------------------------------------------------
	# 5.) Reconcile

	reconcile(summary_table_mat, factors_mat);

	#-----------------------------------------------------------------------------
	# 6.) Apply Summary Table Parameters

	if(!is.null(return_top)){

		normalize_mat=normalize(counts_mat);
		ordered_normalized=reorder_by_decreasing_abundance(normalized_mat);

		extract_top_categories(
			ordered_normalized,
			sumtab[["return_top"]],
			additional_cat=sumtabl[["specific_cat_fn"]]);
	}


	#-----------------------------------------------------------------------------
	# 7.) Shorten categories
	if(!is.null(shorten_cat_names_char)){
		colnames(counts_mat)=
			shorten_category_names(colnames(counts_mat), shorten_cat_names_char);
	}

        # Clean category names a little
	colnames(counts_mat)=clean_category_names(colnames(counts_mat));

	#-----------------------------------------------------------------------------
	
	results=list();
	results[["SummaryTable"]];
	results[["Factors"]];
	results[["PairsMap"]];
	results[["Covariates"]];
	results[["GroupVariables"]];
	results[["RequiredVariables"]];
	return(results);
}


