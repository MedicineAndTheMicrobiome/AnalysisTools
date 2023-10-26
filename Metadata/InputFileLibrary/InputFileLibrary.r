
# This is a collection of shared functions that make the paired analyses
# consistent in how the metadata and samples are reconciled.

source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');

#------------------------------------------------------------------------------

specified=function(var){
	if(is.null(var) || is.na(var) || var==""){
		return(F);
	}else{
		return(T);
	}
}

#------------------------------------------------------------------------------

load_list=function(fname, memo=NULL){

	if(!specified(fname)){
		cat("File name NULL, no list specified for: ", memo, "\n");
		return(NULL);
	}else{
		cat("Loading List: ", fname, " (purpose: ", memo, ")\n", sep="");
 		val=scan(fname, what=character(), comment.char="#");
		return(val);
	}
}


#------------------------------------------------------------------------------

remove_zero_count_categories_and_samples=function(mat){

	# Remove zero count samples
	cat("Checking for samples with no counts.\n");
	tot=apply(mat, 1, sum);
	nonzero=tot>0;
	if(!(all(nonzero))){
		cat("  WARNING: Zero count samples found:\n");
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
		cat("  WARNING: Zero count categories found:\n");
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

	# Remove alls zero categories and counts
	counts_mat=remove_zero_count_categories_and_samples(counts_mat);
	counts_dim=dim(counts_mat);

	cat("  Returned: Num Samples: ", counts_dim[1], "\n");
	cat("  Returned: Num Categories: ", counts_dim[2], "\n");

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
                cat("  Error reading in reference level file: ", fname, "\n");
		cat("  Should be a 2 column tab-separated text file.\n");
		cat("  Format: <variable name>\\t<reference level name>\\n");
                quit(status=-1);
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

load_factors_file=function(fname, prim_key_cname=1, relevel_fn=NULL){

	cat("Loading Factors File: ", fname, "\n", sep="");
	cat("  Primary Key Column Name: ", prim_key_cname, "\n", sep="");
	cat("  Factor Releveling File Specified: ", relevel_fn, "\n", sep="");

	if(!specified(fname)){
		cat("  Factor File not specified.\n");
		return(NULL);
	}

	if(!specified(prim_key_cname)){
		cat("  Primary Key Column Name not specified.  Assuming 1st column is the key.\n");
		prim_key_cname=1;
	}

        factors=data.frame(read.table(fname,  sep="\t", header=TRUE, row.names=prim_key_cname,
                check.names=FALSE, comment.char="#", stringsAsFactors=T));

        factor_names=colnames(factors);
	fact_dim=dim(factors);

	cat("  Num Rows: ", fact_dim[1], "\n", sep="");
	cat("  Num Variables: ", fact_dim[2], "\n", sep="");

	if(specified(relevel_fn)){
		relevel_mat=load_reference_levels_file(relevel_fn);	
		factors=relevel_factors(factors, relevel_mat);
	}else{
		cat("  No Reference Releveling File specified.\n");
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
	cat("  Mapping Column Names:\n");
	print(mapping_cnames);
	

        if(!specified(sbj_id_cname)){
                subject_ids=mapping[,1];
        }else{
		cat("  Subject ID Column Name Specified.\n");
                subject_ids=mapping[,sbj_id_cname];
        }

        column_names=colnames(mapping);
        if(all(column_names!=pA)){
                cat("  Error: Could not find ", pA, " in header of map file.\n");
                quit(status=-1);
        }
        if(all(column_names!=pB)){
                cat("  Error: Could not find ", pB, " in header of map file.\n");
                quit(status=-1);
        }

        map=cbind(as.character(mapping[,pA]), as.character(mapping[,pB]));
        colnames(map)=c(pA, pB);
        rownames(map)=subject_ids;

	cat("  Loaded Pairing Map:\n");
	print(map);

        # Remove pairings with NAs
	cat("  Removing any incomplete mappings...\n");
        incomp=apply(map, 1, function(x){any(is.na(x))});
        map=map[!incomp,];
	num_incomplete=sum(incomp);
	cat("  Number of incomplete pairings removed: ", num_incomplete, "\n");

	cat("Complete Map:\n");
	print(map);

	map_dim=dim(map);
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
	sample_id_keepers=NULL,
	subject_id_keepers=NULL){

	cat("Intersecting pairings map with keep lists.\n");
	num_smp_keepers=length(sample_id_keepers);
	num_sbj_keepers=length(subject_id_keepers);

        num_rows=nrow(pairs_map);
        if(num_rows==0){
		cat("  Pairs map is empty.\n");
		return(pairs_map);
	}

	delete_row_ix=c();

	if(num_smp_keepers){
		for(rix in 1:num_rows){
			if(!any(pairs_map[rix, 1]==sample_id_keepers) ||
				!any(pairs_map[rix, 2]==sample_id_keepers)){
					delete_row_ix=c(delete_row_ix, rix);
			}
		}
	}else{
		cat("No Sample ID Keepers Specified.  No filter applied.\n");
	}

	
	if(num_sbj_keepers){
		cat("  Num Subjects to Keep: ", num_sbj_keepers, "\n", sep="");
		map_sbj_ids=rownames(pairs_map);
		for(rix in 1:num_rows){
			if(!any(map_sbj_ids[rix]==subject_id_keepers)){
				delete_row_ix=c(delete_row_ix, rix);
			}
		}
	}else{
		cat("No Subject ID Keepers Specified. No filter applied.\n");
	}

	cat("Delete Rows: \n");
	print(delete_row_ix);

	complete_row_ix=setdiff(1:num_rows, delete_row_ix);

        incomplete_pairs=pairs_map[delete_row_ix,,drop=F];
        complete_pairs=pairs_map[complete_row_ix,,drop=F];

	cat("  Complete Pairs: \n");
	print(complete_pairs);

	cat("  Incomplete Pairs: \n");
	print(incomplete_pairs);

	cat("\n");
	cat("  Num Input Pairs: ", nrow(pairs_map), "\n");
	cat("  Num Complete Pairs: ", nrow(complete_pairs), "\n");
	cat("  Num Incomplete Pairs: ", nrow(incomplete_pairs), "\n");

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
        cat("  Columns to Available for Extraction: ", num_top_to_extract, "\n");

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
                cat("  Error: Could not find categories: \n");
                print(missing_cat);
                quit(status=-1);
        }

        # :: Remove categories we have already extracted in the top N
        already_extracted_cat=colnames(top_cat);
        extra_cat=setdiff(additional_cat, already_extracted_cat);

        num_extra_to_extract=length(extra_cat);
        cat("  Num Extra Categories to Extract: ", num_extra_to_extract, "\n");

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
		cat("  Error:  Overlap detected between Covariates and Group Lists:\n");
		print(overlap);
		quit(-1);
	}

	# Checking required are subset of covariates or group list
	combined=c(covar, grp);
	overlap=intersect(combined, req);
	if(length(overlap)==num_req){
		cat("  Good.  Required variables overlap with Covariates or Group List.\n");
	}else{
		cat("  Error: Required variables do not overlap with Covariates or Group List:\n");
		cat("  Covariates:\n");
		print(covar);
		cat("  Group:\n");
		print(grp);
		cat("  Required:\n");
		print(req);
		quit(-1);
	}

	# Checking group and covariates are in available list
	overlap=intersect(combined, avail);
	if(length(overlap)==length(combined)){
		cat("  Good.  All variables (Covariates and Group) found/available.\n");
	}else{
		cat("  Error.  Some variables missing from available requested as Covariates/Group:\n");
		print(overlap);
		quit(-1);
	}
	
	return(1);
		
}


reconcile=function(param=list("summary_table_mat"=NULL, "factor_mat"=NULL, "pairs_mat"=NULL,
	"sbj_smp_mat"=NULL)){	

	# If sumtab, factors, and pairs are specified, then the keying must be
	#   sample_id -> sumtab
	#   subject_id -> factors
	#   subject_id -> pairs
	#
	# If sumtab and factors are specified, then
	#   they need to have the same IDs keyed, probably sample id
	
	summary_table_mat=param[["summary_table_mat"]];
	factor_mat=param[["factor_mat"]];
	pairs_mat=param[["pairs_mat"]];
	sbj_to_smp_mat=param[["sbj_smp_mat"]];

	# Get sample IDs from summary table
	sumtab_sample_ids=rownames(summary_table_mat);
	factor_subject_ids=rownames(factor_mat);

	if(specified(pairs_mat)){
	# If pairs mat specified, use it to bridge the factor to sample IDs

		# Determine which pairs are complete by sample ID
		pairs_mat=intersect_pairings_map_by_keep_list(
			pairs_mat, 
			sample_id_keepers=sumtab_sample_ids);

		# Determine which pairs are useable by subject ID
		pairs_mat=intersect_pairings_map_by_keep_list(
			pairs_mat, 
			subject_id_keepers=factor_subject_ids);

		# This are the subject/sample ids that are pairable
		pairable_sbj_ids=sort(rownames(pairs_mat));
		pairable_smp_ids=sort(c(pairs_mat[,1], pairs_mat[,2]));

		summary_table_mat=summary_table_mat[pairable_smp_ids,,drop=F];	
		factor_mat=factor_mat[pairable_sbj_ids,,drop=F];

	}else if(specified(sbj_to_samp_mat)){
	# If 1-to-1 mapping specified use it.

		invert_mapping=function(map){
			inv=matrix(rownames(map), ncol=1);	
			rownames(inv)=map[,1];
			return(inv);
		}

		map_sbj_ids=rownames(sbj_to_smp_mat);
		map_smp_ids=sbj_to_smp_mat[,1];

		# Remove mappings missing subjects from factors
		matched_sbj_ids=intersect(map_sbj_ids, factor_subject_ids);
		sbj_to_smp_mat=sbj_to_smp_mat[matched_sbj_ids,,drop=F];

		# Invert mappings
		smp_to_sbj_mat=invert_mapping(sbj_to_smp_mat);
		
		# Remove mappings missing samples from sumtab
		map_smp_ids=rownames(smp_to_sbj_mat);
		matched_smp_ids=intersect(map_smp_ids, sumtab_sample_ids);
		smp_to_sbj_mat=smp_to_sbj_mat[matched_smp_ids,,drop=F];

		matched_smp_ids=rownames(smp_to_sbj_mat);
		matched_sbj_ids=smp_to_sbj_mat[,1];

		# Invert mappings back
		sbj_to_smp_mat=invert_mapping(smp_to_sbj_mat);

		summary_table_mat=summary_table_mat[matched_smp_ids,,drop=F];	
		factor_mat=factor_mat[matched_sbj_ids,,drop=F];

	}else{
		# If only sumtab and factors specified, then assume both factors and 
		# summary table are keying by sample ID.

		shared_ids=sort(intersect(sumtab_sample_ids, factor_subject_ids));
		
		summary_table_mat=summary_table_mat[shared_ids,,drop=F];	
		factor_mat=factor_mat[shared_ids,,drop=F];

	}

	results=list();
	results[["summary_table_mat"]]=summary_table_mat;
	results[["factor_mat"]]=factor_mat;
	results[["pairs_mat"]]=pairs_mat;
	results[["sbj_smp_mat"]]=sbj_to_smp_mat;

	return(results);
}

load_and_reconcile_files=function(
	sumtab=list(fn=NULL, shorten_cat_names_char=NULL, return_top=NULL, specific_cat_fn=NULL), 
	factors=list(fn=NULL, sbj_cname=NULL, ref_relvl_fn=NULL), 
	pairs=list(fn=NULL, a_cname=NULL, b_cname=NULL, subject_id_cname=NULL), 
	sbj_to_smp=list(fn=NULL, sbjid_cname=NULL, smpid_cname=NULL),
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

	summary_table_mat=load_summary_file(fname=sumtab[["fn"]]);

	spec_cat_arr=load_list(fname=sumtab[["specific_cat_fn"]], memo="Specific Categories File");

	factors_mat=load_factors_file(
		fname=factors[["fn"]],
		prim_key_cname=factors[["sbj_cname"]],
		relevel_fn=factors[["ref_relvl_fn"]]
		);

	pairs_mat=load_mapping(
		fname=pairs[["fn"]], 
		pA=pairs[["a_cname"]], 
		pB=pairs[["b_cname"]], 
		sbj_id_cname=pairs[["subject_id_cname"]]
		);

	sbj_smp_mat=load_sbj_smp_mapping(
		fname=sbj_to_smp[["fn"]],
		sbj_cname=sbj_to_smp[["sbjid_cname"]],
		smp_cname=sbj_to_smp[["smpid_cname"]]
		);

	covariates_arr=load_list(fname=covariates[["fn"]], memo="Covariates File");

	groupvar_arr=load_list(fname=grpvar[["fn"]], memo="Group Variables File");

	requiredvar_arr=load_list(fname=reqvar[["fn"]], mem="Required Variables File");

	#-----------------------------------------------------------------------------
	# 2.) Keep/Subset target variables

	status=check_variables(covariates_arr, groupvar_arr, requiredvar_arr, colnames(factors_mat));
	
	cat("Subsetting requested variables from factors.\n");
	factors_subset_mat=factors_mat[, c(covariates_arr, groupvar_arr), drop=F];

	#-----------------------------------------------------------------------------
	# 3.) Reconcile

	reconciled_files=reconcile(param=list(
		summary_table_mat=summary_table_mat, 
		factor_mat=factors_subset_mat, 
		pairs_mat=pairs_mat, 
		sbj_smp_mat=sbj_smp_mat
		));

	recon_factors=reconciled_files[["factor_mat"]];

	#-----------------------------------------------------------------------------
	# 4.) Remove NAs (remove subjects with NAs)

	nona_factors_fn=paste(gsub("\\.tsv", "", factors[["fn"]]), ".noNAs", sep="");
	factors_wo_nas_res=remove_sample_or_factors_wNA_parallel(
		recon_factors,
        	required=requiredvar_arr,
		outfile=nona_factors_fn,
		num_trials=64000, num_cores=64);

	factors_mat_nona=factors_wo_nas_res[["factors"]];

	#-----------------------------------------------------------------------------
	# 5.) Reconcile
	
	# Replace factor mat, but keep the other data/matrices the same
	reconciled_files[["factor_mat"]]=factors_mat_nona;
	reconciled_files=reconcile(reconciled_files);

	#-----------------------------------------------------------------------------
	# 6.) Apply Summary Table Parameters

	# Order categories by decreasing abundance
	counts_mat=reconciled_files[["summary_table_mat"]];
	normalized_mat=normalize(counts_mat);
	normalized_mat=reorder_by_decreasing_abundance(normalized_mat);
	counts_mat=counts_mat[,colnames(normalized_mat), drop=F];

	if(specified(sumtab[["return_top"]])){

		top_cat_normalized=extract_top_categories(
			normalized_mat,
			sumtab[["return_top"]],
			additional_cat=spec_cat_arr);

		normalized_mat=top_cat_normalized;
		counts_mat=extract_categories_for_counts(counts_mat, colnames(normalized_mat));
	}

	#-----------------------------------------------------------------------------
	# 7.) Shorten categories
	if(specified(sumtab[["shorten_cat_names_char"]])){

		colnames(normalized_mat)=
			shorten_category_names(colnames(normalized_mat), 
			sumtab[["shorten_cat_names_char"]]);

		colnames(counts_mat)=colnames(normalized_mat);
	}

        # Clean category names a little
	colnames(normalized_mat)=clean_category_names(colnames(normalized_mat));
	colnames(counts_mat)=colnames(normalized_mat);

	#-----------------------------------------------------------------------------
	
	results=list();
	results[["SummaryTable_counts"]]=counts_mat;
	results[["SummaryTable_normalized"]]=normalized_mat;
	results[["Factors"]]=factors_mat_nona;
	results[["PairsMap"]]=reconciled_files[["pairs_mat"]];
	results[["SubjectToSampleIDMap"]]=reconciled_files[["sbj_to_smp_id_mat"]];
	results[["Covariates"]]=covariates_arr;
	results[["GroupVariables"]]=groupvar_arr;
	results[["RequiredVariables"]]=requiredvar_arr;

	return(results);
}


