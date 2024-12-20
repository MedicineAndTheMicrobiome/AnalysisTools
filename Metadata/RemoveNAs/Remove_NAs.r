

remove_sample_or_factors_wNA=function(factors, num_trials=500000, verbose=T){

	if(verbose){
		cat("Identifying Samples or Factors to remove to remove all NAs:\n");
	}

        num_col=ncol(factors);
        num_row=nrow(factors);

	if(verbose){
		cat("Categories before:\n");
		print(colnames(factors));

		cat("Num Factors Entering: ", num_col, "\n");
		cat("Num Samples Entering: ", num_row, "\n");
		cat("\n");
	}

	if(verbose){
		cat("Removing columns/factors with no variation/information...\n");
	}
	# Remove columns with no information
	factors=remove_no_variation_factors(factors, verbose);
	num_col=ncol(factors);
        num_row=nrow(factors);

        # Find rows and columns with NAs
        row_na_ix=which(apply(factors, 1, anyNA));
        col_na_ix=which(apply(factors, 2, anyNA));

        row_na_counts=apply(factors[row_na_ix,, drop=F], 1, function(x){sum(is.na(x))});
        col_na_counts=apply(factors[,col_na_ix, drop=F], 2, function(x){sum(is.na(x))});

        #row_na_counts=apply(factors[row_na_ix,, drop=F], 1, function(x){sum(is.na(x))/num_col});
        #col_na_counts=apply(factors[,col_na_ix], 2, function(x){sum(is.na(x))/num_row});
        combined_na_counts=c(row_na_counts, col_na_counts);

        num_row_na=length(row_na_ix);
        num_col_na=length(col_na_ix);

        if((num_row_na+num_col_na) == 0){
                return(factors);
        }

        # Assign an ID to rows and columns, so when we remove them, we won't have to
        # remap the indices
        na_ix=c(paste("r", row_na_ix, sep=""), paste("c", col_na_ix, sep=""));
        dim_ix=c(rep(1, num_row_na), rep(2, num_col_na));

	if(verbose){
		cat("NAs found in these rows/columns:\n");
		print(na_ix);
	}

        num_non_na=sum(!is.na(factors));
        #print(factors);
	if(verbose){
		cat("Maximum Number of non-NAs: ", num_non_na, "\n");
	}

        num_rowcol=sum(num_row_na, num_col_na);

        # Rename rows/columns
        renamed_factors=factors;
        colnames(renamed_factors)=paste("c", 1:num_col, sep="");
        rownames(renamed_factors)=paste("r", 1:num_row, sep="");

        # Keep track of the sequence of row/col removals that maximize the remaining data
        best_sequence_nonna=0;
        best_sequence_ix=numeric();

        max_no_improvement=0.05*num_trials;
	if(verbose){
		cat("Max allowable trials without improvement: ", max_no_improvement, "\n");
	}
	heartbeat_periodicity=ceiling(num_trials/10);

        # Repeatedly search for alternative sets that maximize remaining data
        last_improvement=0;
        for(i in 1:num_trials){

		if((i %% heartbeat_periodicity)==0){
			cat(".");
		}

		random_ix=sample(num_rowcol, replace=F, prob=combined_na_counts);
                # random_ix=sample(num_rowcol, replace=F);

                tmp_matrix=renamed_factors;
                cur_ix_sequence=numeric();

                # Step through a random sequence rows/columns to remove
                for(ix in random_ix){

                        rm_ix=na_ix[ix];
                        cur_ix_sequence=c(cur_ix_sequence, ix);

                        if(dim_ix[ix]==1){
                                names=rownames(tmp_matrix);
                                new_ix=which(names==rm_ix);
                                tmp_matrix=tmp_matrix[-new_ix,, drop=F];
                        }else{
                                names=colnames(tmp_matrix);
                                new_ix=which(names==rm_ix);
                                tmp_matrix=tmp_matrix[,-new_ix, drop=F];
                        }

                        # Count up number of NAs
                        num_na=sum(is.na(tmp_matrix));

                        # If there are no more NAs, stop
                        if(num_na==0){
                                break;
                        }
                }

                # Count up number of non NAs
                num_non_na=sum(!is.na(tmp_matrix));

                # If the number of non NAs is greater than before keep it
                if(num_non_na>best_sequence_nonna){
                        best_sequence_nonna=num_non_na;
                        best_sequence_ix=cur_ix_sequence;
                        last_improvement=0;
                        cat("Num Non NAs: ", best_sequence_nonna, "\n");
                }

                last_improvement=last_improvement+1;
                if(last_improvement>max_no_improvement){
			cat("No more significant improvements...\n");
                        break;
                }
        }
	cat("\n");

	if(best_sequence_nonna==0){
		# If no solution found, return empty matrix
		fact_subset=factors[c(),c(),drop=F];
	}else{
		best_na_ix=na_ix[best_sequence_ix];
		best_dim_ix=dim_ix[best_sequence_ix];

		# Strip r/c from name to recover original index
		rm_row=as.integer(gsub("r","", best_na_ix[best_dim_ix==1]));
		rm_col=as.integer(gsub("c","", best_na_ix[best_dim_ix==2]));

		# Remove rows/columns
		fact_subset=factors;
		if(length(rm_row)){
			fact_subset=fact_subset[-rm_row,,drop=F];
		}
		if(length(rm_col)){
			fact_subset=fact_subset[,-rm_col,drop=F];
		}

		fact_subset=remove_no_variation_factors(fact_subset, verbose);
	}

	return(fact_subset);
}

###############################################################################

summarize_as_text=function(orig_factors, after_req_var, after_na_removal){

	original_var=colnames(orig_factors);
	final_var=colnames(after_na_removal);
	removed_var=setdiff(original_var, final_var);

	text=c(
		"Original Factors Variables:",
		paste("  Num Samples: ", nrow(orig_factors), sep=""),
		paste("  Num Factors: ", ncol(orig_factors), sep=""),
		"After Requiring Variables:",
		paste("  Num Samples: ", nrow(after_req_var), sep=""),
		paste("  Num Factors: ", ncol(after_req_var), sep=""),
		"After NA Removal Variables:",
		paste("  Num Samples: ", nrow(after_na_removal), sep=""),
		paste("  Num Factors: ", ncol(after_na_removal), sep=""),
		"",
		"Kept Variables:",
		capture.output(print(final_var)),
		"",
		"Removed Factors:",
		capture.output(print(removed_var))
	);

	return(text);
}

remove_no_variation_factors=function(factors, verbose=T){
	no_info=c();
	cnames=colnames(factors);
	num_col=ncol(factors);

	if(ncol(factors)>0){
		for(i in 1:num_col){
			val=factors[,i];
			val=val[!is.na(val)];
			uniqs=unique(val);
			if(length(uniqs)==1){
				no_info=c(no_info, i);
				#if(verbose){
				cat(cnames[i], ": Has no useful information.  It is all:\n");
				print(uniqs);
				#}
			}
		}
		if(length(no_info)>0){
			factors=factors[,-no_info, drop=F];
		}
	}

	return(factors);
}

rem_missing_var_from_modelstring=function(model_string, kept_variables){

	model_string=gsub("\\s+","", model_string);
	# In case LHS was specified removed it
	formula_parts=strsplit(model_string, "~")[[1]];
	num_parts=length(formula_parts);

	if(num_parts==2){
		rhs=formula_parts[2];		
		lhs=formula_parts[1];
	}else{
		rhs=formula_parts[1];
		lhs="";
	}

	# Split formula into linear components
	lin_comp_arr=strsplit(rhs, "\\+")[[1]];

	# Split interaction terms into main effects
	kept_comp=character();
	for(lin_comp in lin_comp_arr){
		int_term_arr=unique(strsplit(lin_comp, "[:\\*]")[[1]]);
		num_terms=length(int_term_arr);
		num_intersect=length(unique(intersect(int_term_arr, kept_variables)));
		if(num_intersect==num_terms){
			kept_comp=c(kept_comp, lin_comp);	
		}
	}

	# Rebuild model string
	if(num_parts==2){
		# Attach response (LHS) to predictors (RHS)
		model_string=paste(lhs, " ~ ", paste(kept_comp, collapse=" + "), sep="");
	}else{
		# Only build RHS
		model_string=paste(kept_comp, collapse=" + ");
	}		
	return(model_string);

}

#rem_missing_var_from_modelstring("this ~ is + a + a*test+ b", c("is"))


get_var_from_modelstring=function(model_string){

	model_string=gsub("\\s+","", model_string);

	# In case LHS was specified removed it
	formula_parts=strsplit(model_string, "~")[[1]];
	if(length(formula_parts)==2){
		rhs=formula_parts[2];		
	}else{
		rhs=formula_parts[1];
	}

	# Split formula into linear components
	rhs=gsub("\\s+","", rhs);
	lin_comp_arr=strsplit(rhs, "\\+")[[1]];

	# Split interaction terms into main effects
	variables=character();
	for(lin_comp in lin_comp_arr){
		int_term_arr=unique(strsplit(lin_comp, "[:\\*]")[[1]]);
		int_term_arr=gsub("^\\s+","", int_term_arr);
		int_term_arr=gsub("\\s+$","", int_term_arr);
		variables=c(int_term_arr, variables);
	}

	return(unique(variables));

}

remove_samples_from_st=function(rm_na_result, summary_table){
	factor_samp_ids=rownames(rm_na_result$factors);
	summary_table_samp_ids=rownames(summary_table);
	shared_samp_ids=sort(intersect(factor_samp_ids, summary_table_samp_ids));
	return(summary_table[shared_samp_ids,,drop=F]);
}

# rem_missing_var_from_modelstring("apples+oranges+apples*oranges+grapes+grapes:oranges+grapes*oranges", c("apples", "oranges"));

###############################################################################


library(foreach);
library(doMC);


remove_samples_for_req_var=function(factors, req_var){
	# This function will remove all samples necessary to keep the required variables.
	fact_avail=colnames(factors);
	missing_var=setdiff(req_var, fact_avail);
	if(length(missing_var)>0){
		cat("\nError: Required variable(s) missing...\n");	
		print(missing_var);
		cat("\n");
		quit(status=-1);
	}
	req_var_mat=factors[, req_var, drop=F];
	na_mat=is.na(req_var_mat);
	rows_wnas=apply(na_mat, 1, any);
	return(factors[!rows_wnas,, drop=F]);
}

remove_sample_or_factors_wNA_parallel=function(factors, required_variables=NULL, num_trials=500000, num_cores=16, outfile=""){

	nsamples=nrow(factors);
	nfactors=ncol(factors);
	numNAs=sum(is.na(factors));
	nreqvar=length(required_variables);

	cat("Num Samples Entering: ", nsamples, "\n");
	cat("Num Factors Entering: ", nfactors, "\n");
	cat("Num NAs Entering:", numNAs, "\n");
	cat("Num Required Variables: ", nreqvar, "\n");
	cat("\n");

	orig_factors=factors;
	if(nreqvar>0){
		cat("Required Variables:\n");
		print(required_variables);

		factors=remove_samples_for_req_var(factors, required_variables);
	}

	factors=remove_no_variation_factors(factors, verbose);

	res=registerDoMC(num_cores);
	core_trials=ceiling(num_trials/num_cores);

	cat("Trials per core: ", core_trials, "\n");

	results=foreach(i = 1:num_cores) %dopar% {
		remove_sample_or_factors_wNA(factors=factors, num_trials=core_trials, verbose=F);
	}	

	sizes=numeric(num_cores);
	for(i in 1:num_cores){
		matdim=dim(results[[i]]);
		sizes[i]=matdim[1]*matdim[2];
	}

	cat("Non-NA Matrix Sizes: \n");
	print(sizes);

	which_max=min(which(max(sizes)==sizes));
	cat("Largest Matrix size: ", sizes[which_max], "\n");

	best_matrix=results[[which_max]];

	nsamples=nrow(best_matrix);
	nfactors=ncol(best_matrix);
	
	cat("Num Samples Leaving: ", nsamples, "\n");
	cat("Num Factors Leaving: ", nfactors, "\n");

	if(outfile!=""){
		outfile=gsub("\\.tsv$", "", outfile);
		write.table(
        		best_matrix,
			paste(outfile, ".noNAs.tsv", sep=""),
			quote=F, sep="\t", row.names=T, col.names=NA);
	}

	# Prep return values
	results=list();
	results$factors=best_matrix;
	results$summary_text=summarize_as_text(orig_factors, factors, best_matrix);

	return(results);
}



