
# This function computes the value for a subject (predicted y) based on
#   the subject's available x values that were not NAs (target_predictors).  
#   This is done by fitting a linear model based on the other 
#   subject's y's (responses) and x's (predictors), for which there is data.
#
# If there were fewer samples than predictors, then only the predictors
#   with the greatest correlation with responses will be used.
#
# If there were NAs in a predictor column, then they will be replaced with
#   the median value for that predictor.


impute_cell=function(target_predictors, responses, predictors, verbose=F){

	# target_predictors are the predictors/available data we have for the
	#    the subject's missing cell.
	
	# responses and predictors are the other subjects data for training

        num_samples=nrow(responses);
        avail_predictors=ncol(predictors);

        num_pred_to_use=min(num_samples-2, avail_predictors);

        cat("Num Samples Available for Imputation: ", num_samples, "\n");
        cat("Num Available Predictors: ", avail_predictors, "\n");
        cat("Num Predictors to Use: ", num_pred_to_use, "\n");

        # Compute correlation between response and predictors
        corr=apply(predictors, 2, function(x){
                non_na=!is.na(x);
                cor(x[non_na], responses[non_na,]);}
                );

        # Order predictors by decreasing correlation
        corr_order=order(abs(corr), decreasing=T);
        corr_ordered=(corr[corr_order]);
        #print(corr_ordered);
        top_pred=names(corr_ordered)[1:num_pred_to_use];


        # Replace NAs with median value from the same column/predictor
        predictors_nona=predictors;
        for(j in 1:ncol(predictors_nona)){
                cur_median=median(predictors[,j], na.rm=T);
                predictors_nona[is.na(predictors[,j]),j]=cur_median;
        }

        # Create model with the predictors with greatest correlaton with response first
        form_str=paste(colnames(responses), "~", paste(top_pred, collapse="+"));
        print(form_str);
        fit=lm(formula(form_str), data=as.data.frame(cbind(responses, predictors_nona)));


	coefficients=fit$coefficients;
	if(any(is.na(coefficients))){
		cat("\n");
		cat("Warning: NAs detected in estimated coefficients.\n");
		cat("  Refitting model with subset of predictors.\n");

		nonNA_predictors=setdiff(names(coefficients[!is.na(coefficients)]), "(Intercept)");
		cat("Non-NA Predictors:\n");
		print(nonNA_predictors);

		form_str=paste(colnames(responses), "~", paste(nonNA_predictors, collapse="+"));
		print(form_str);
		fit=lm(formula(form_str), data=as.data.frame(cbind(responses, predictors_nona)));
		
		top_pred=nonNA_predictors;

	}

        # Predict NA with the values that we have
        imputed_val=predict(fit, newdata=as.data.frame(target_predictors[,top_pred, drop=F]));

        #print(fit);
        if(verbose){
                print(summary(fit));
        }
        #eg=rbind(predictors, target_predictors);
        #print(eg);
        #x=predict(fit, new=eg);
        obs_resp_range=range(responses);
        cat("--------------------------------------------------------\n");
        cat("Name: ", colnames(responses), "\n");
        cat("Response Range: ", obs_resp_range[1], " - ", obs_resp_range[2], "\n");
        cat("Inputed Value: ", imputed_val, "\n");
        if(imputed_val<obs_resp_range[1] || imputed_val>obs_resp_range[2]){
                cat("WARNING: Inputed Value Outside Range of Observed Values\n");
        }
        cat("--------------------------------------------------------\n");
        return(imputed_val);
}

###############################################################################

# This is the main function which takes in a matrix with NAs dispersed it and
#   computes values to replace them, then returns the filled matrix.
# 
# Cells with NAs are identified, then the impute_cell function is called
#   

impute_matrix=function(mat_wna){

        num_rows=nrow(mat_wna);
        num_cols=ncol(mat_wna);

        cat("Original Matrix Dimensions: ", num_rows, " x ", num_cols,
                " matrix. (", num_rows*num_cols, ")\n", sep="");

        # remove rows with all NAs
        usable_rows=apply(mat_wna, 1, function(x){!all(is.na(x))});
        usable_cols=apply(mat_wna, 2, function(x){!all(is.na(x))});
        usable_mat_wna=mat_wna[usable_rows,usable_cols];

        usable_num_rows=nrow(usable_mat_wna);
        usable_num_cols=ncol(usable_mat_wna);

        cat("Removed rows/cols with all NAs: ", usable_num_rows, " x ", usable_num_cols,
                " matrix. (", usable_num_rows * usable_num_cols,")\n", sep="");

        #print(mat_wna);

        cat("Looking for NAs...\n");
        na_pos=numeric();

        if(1){
                for(rix in 1:usable_num_rows){
                        for(cix in 1:usable_num_cols){
                                if(is.na(usable_mat_wna[rix, cix])){
                                        na_pos=rbind(na_pos, c(rix, cix));
                                }
                        }
                }
        }else{
                # For validation
                for(rix in 1:usable_num_rows){
                        for(cix in 1:usable_num_cols){
                                na_pos=rbind(na_pos, c(rix, cix));
                        }
                }
        }

        if(length(na_pos)==0){
                num_nas_to_impute=0;
        }else{
                num_nas_to_impute=nrow(na_pos);
        }

        cat("Num NAs to try to impute:",  num_nas_to_impute, "\n");

        filled_matrix=usable_mat_wna;

        if(!is.null(num_nas_to_impute) && num_nas_to_impute>0){
                for(na_ix in 1:num_nas_to_impute){

                        target_row=na_pos[na_ix,1];
                        target_column=na_pos[na_ix,2];

                        cell=usable_mat_wna[target_row, target_column, drop=F];

                        if(!is.na(cell)){
                                cat("Error:  Trying to input cell not NA.\n");
                                quit();
                        }

                        cat("(", na_ix, "/", num_nas_to_impute, ") Imputing: ",
                                rownames(cell), " / ", colnames(cell), "\n");

                        # Identify data to impute with (i.e. excluding target row/col)
                        non_na_row=!is.na(usable_mat_wna[,target_column,drop=F]);
                        non_na_col=!is.na(usable_mat_wna[target_row,,drop=F]);

                        # target_predictors: the predictors x (from same subject)
                        #     to use to guess value (response y)  of missing data from subject
                        # responses: the y's of peers to build prediction model from
                        # predictors: the x's of peers to build prediction model from

                        imputed_val=impute_cell(
                                target_predictors=usable_mat_wna[target_row, non_na_col, drop=F],
                                responses=usable_mat_wna[non_na_row, target_column, drop=F],
                                predictors=usable_mat_wna[non_na_row, non_na_col, drop=F]
                                );

                        # Store imputed values, so we don't impute new values with previously
                        # imputed values.
                        filled_matrix[target_row, target_column]=imputed_val;

                }
        }

        repl_rows=rownames(filled_matrix);
        mat_wna[repl_rows,]=filled_matrix[repl_rows,];
        return(mat_wna);

}

###############################################################################
# Matrix_wNAs: contains a matrix with NAs in it
# variable_grouping_list: contains a R list, keyed with group names.  
#     Each list element is a vector of variables names in the Matrix_wNAs.

impute_by_groupings=function(matrix_wNAs, variable_grouping_list){

	#print(matrix_wNAs);
	#print(variable_grouping_list);

	num_var=ncol(matrix_wNAs);
	num_samples=nrow(matrix_wNAs);

	samp_names=rownames(matrix_wNAs);
	#cat("Sample Names:\n");
	#print(samp_names);
	#cat("\n");
	var_names=colnames(matrix_wNAs);
	#cat("Variable Names:\n");
	#print(var_names);
	#cat("\n");

	out_matrix=matrix(NA, nrow=num_samples, ncol=num_var);
	rownames(out_matrix)=samp_names;
	colnames(out_matrix)=var_names;

	group_names=names(variable_grouping_list);

	for(cur_grp in group_names){

		cat("Imputing Group: ", cur_grp, "\n");
		targ_var=variable_grouping_list[[cur_grp]];

		cat("Variables: \n");
		print(targ_var);
		
		targ_matrix_wNAs=matrix_wNAs[,targ_var,drop=F];

		noNA_matrix=impute_matrix(targ_matrix_wNAs);	

		# Not sure why I need to iterate to make assignment
		#   out_matrix[samp_names,targ_var]=noNA_matrix[samp_names,targ_var];
		for(var_ix in targ_var){
			out_matrix[samp_names, var_ix]=noNA_matrix[samp_names, var_ix];
		}

	}

	return(out_matrix);

}
