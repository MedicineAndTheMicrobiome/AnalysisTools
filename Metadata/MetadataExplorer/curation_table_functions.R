

VariableInfo.init=function(num_variables){

	var_info_colnames=c(
	
		# All Variables
		"VariableName",
		"IsVarnameValid",
		"NumberNAs",
		"PercentNAs",
		"NumberUnique",
		"IsDichotomous",
		
		# Numeric
		"IsNumeric",
		"IsBoolean",
		"PassedShapiroWilks",
		"ImprovedBySqrtTransform",
		"ImprovedByLogTransform",
		"ImprovedByLogitTransform",

		# Date
		"IsPotentialDate",
		"IsIsoDate",

		# Character
		"IsCategorical",
		"AreLevelNamesValid"
		
	);
	
	num_var_info_colnames=length(var_info_colnames);

	VariableInfo=as.data.frame(matrix(NA, nrow=num_variables, ncol=num_var_info_colnames));
	colnames(VariableInfo)=var_info_colnames;
	VariableInfo<<-VariableInfo;
	return;
}

VariableInfo.check_variable_name=function(vname){
	# Return T, if ok, F, if not ok.
	non_valid_chars=gsub("[a-zA-Z0-9\\.\\_]", "", vname);
	num_nonvalid=nchar(non_valid_chars);
	if(num_nonvalid>0){
		return(F);
	}else if(length(grep("^[0-9]", vname))>0){
		return(F);
	}else {
		return(T);
	}
}

VariableInfo.check_all_level_names=function(values){
	unique_levels=unique(values);
	num_uniq=length(unique_levels);
	any_invalid=0;
	for(i in 1:num_uniq){
		any_invalid=any_invalid+!VariableInfo.check_variable_name(unique_levels[i]);
	}
	return(any_invalid==0);
}

VariableInfo.check_date=function(values){
	
	num_values=length(values);
	slashes_found=length(grep("/", values))==num_values;
	commas_found=length(grep(",", values))==num_values;
	dashes_found=length(grep("-", values))==num_values;
	numbers_found=length(grep("[0-9]", values))==num_values;
	
	if(!numbers_found){
		return(-1); # Not a date
	}else if(dashes_found){
		if(all(!is.na(as.Date(values, "%Y-%m-%d")))){
			return(1); # Iso-8601
		}
		return(0); # Potential Date
	}else if(slashes_found || commas_found){
		return(0); # Potential Date
	}else{
		return(-1); # Not a date
	}
}

VariableInfo.build=function(metadata){
	
	num_variables=ncol(metadata);
	variable_names=colnames(metadata);
	num_samples=nrow(metadata);
	
	VariableInfo.init(num_variables);
	VariableInfo[,"VariableName"]=variable_names;
	VariableInfo[,"IsVarnameValid"]=sapply(VariableInfo[,"VariableName"], VariableInfo.check_variable_name);
	
	#print(metadata);
	
	for(i in 1:num_variables){
		
		cat("Examining ", variable_names[i], "\n");
		cur_val=metadata[,i];
		print(cur_val);
		
		# Check all variables
		na_ix=is.na(cur_val);
		non_na_val=cur_val[!na_ix];
		
		VariableInfo[i,"NumberNAs"]=sum(na_ix);
		VariableInfo[i,"PercentNAs"]=round((VariableInfo[i,"NumberNAs"]/num_samples*100), 2);
		VariableInfo[i,"NumberUnique"]=length(unique(non_na_val));
		VariableInfo[i,"IsDichotomous"]=VariableInfo[i,"NumberUnique"]==2;
		
		# Treat as number
		as_number=as.numeric(non_na_val);
		num_numerics=sum(!is.na(as_number));
		is_numeric=num_numerics>0;
		VariableInfo[i, "IsNumeric"]=is_numeric;
		
		# Treat as Boolean
		if(is_numeric && VariableInfo[i,"IsDichotomous"] && min(as_number)==0 && max(as_number)==1){
			# If 0's and 1's
			VariableInfo[i, "IsBoolean"]=T;
		}else if(is.logical(non_na_val)){
			# If T and F (logical)
			VariableInfo[i, "IsBoolean"]=T;
		}else{
			VariableInfo[i, "IsBoolean"]=F;
		}
		
		if(is_numeric){
		
			# Check as number
			raw_shpwlk_pvalue=shapiro.test(as_number)$p.value;
			VariableInfo[i, "PassedShapiroWilks"]=raw_shpwlk_pvalue>0.05;
			if(!VariableInfo[i, "PassedShapiroWilks"]){
				if(all(as_number>=0)){
				
					# Transform and Test
					#  sqrt
					sqrt_trans=sqrt(as_number);
					sqrt_shpwlk_pvalue=shapiro.test(sqrt_trans)$p.value;
					
					#  log
					log_trans=log(as_number+1);
					log_shpwlk_pvalue=shapiro.test(log_trans)$p.value;
					
					#  logit
					val_range=range(as_number);
					if(val_range[1]>0 && val_range[2]<1){
						logit_trans=log(as_number/(1-as_number));
						logit_shpwlk_pvalue=shapiro.test(logit_trans)$p.value;
					}else{
						logit_shpwlk_pvalue=0;
					}
					
					# Look for least not normal transform
					pvals=c(raw_shpwlk_pvalue, sqrt_shpwlk_pvalue, log_shpwlk_pvalue, logit_shpwlk_pvalue);
					print(pvals);
					max_pval=max(pvals);
					best_transform=min(which(pvals==max_pval));
					
					VariableInfo[i,"ImprovedBySqrtTransform"]=F;
					VariableInfo[i,"ImprovedByLogTransform"]=F;
					VariableInfo[i,"ImprovedByLogitTransform"]=F;
					
					if(best_transform==2){
						VariableInfo[i,"ImprovedBySqrtTransform"]=T;
					}else if(best_transform==3){
						VariableInfo[i,"ImprovedByLogTransform"]=T;
					}else if(best_transform==4){
						VariableInfo[i,"ImprovedByLogitTransform"]=T;
					}
				}
			}
			
		}else{
		
			# Check as date
			date_type=VariableInfo.check_date(non_na_val);
			if(date_type==1){
				VariableInfo[i,"IsIsoDate"]=T;
				VariableInfo[i,"IsPotentialDate"]=T;
			}else if(date_type==0){
				VariableInfo[i,"IsIsoDate"]=F;
				VariableInfo[i,"IsPotentialDate"]=T;
			}else{
				VariableInfo[i,"IsIsoDate"]=F;
				VariableInfo[i,"IsPotentialDate"]=F;
			}
			
			if(date_type==-1){
			
				# Check as categorical
				VariableInfo[i,"IsCategorical"]=T;
				
				all_clear=VariableInfo.check_all_level_names(non_na_val);
				if(all_clear){
					VariableInfo[i,"AreLevelNamesValid"]=T;
				}else{
					VariableInfo[i,"AreLevelNamesValid"]=F;
				}
				
			}else{
				VariableInfo[i,"IsCategorical"]=F;
			}

		}
	}
	
	VariableInfo<<-VariableInfo;
	
}

#------------------------------------------------------------------------------
WarningsTable.warnings_str=list();

WarningsTable.warnings_str[["INV_VAR_NAME"]]="Invalid variable name";
WarningsTable.warnings_str[["NA_GT50"]]="Significant proportion of NAs: >50%";
WarningsTable.warnings_str[["NA_GT25"]]="Noteworthy proportion of NAs: >25%";
WarningsTable.warnings_str[["NA_GT10"]]="Non-negligible proportion of NAs: >10%";
WarningsTable.warnings_str[["DICH_NOTBOOL"]]="Dichotomous but not Boolean";
WarningsTable.warnings_str[["ALL_IDENT"]]="All values identical";
WarningsTable.warnings_str[["REC_SQRT_TRANS"]]="Consider sqrt() transformation";
WarningsTable.warnings_str[["REC_LOG_TRANS"]]="Consider log() transformation";
WarningsTable.warnings_str[["REC_LOGIT_TRANS"]]="Consider logit() transformation";
WarningsTable.warnings_str[["NON_ISODATE"]]="Non-ISO 8601 Date";
WarningsTable.warnings_str[["INV_LVL_NAMES"]]="Invalid level names";
WarningsTable.warnings_str[["EXCESS_CATS"]]="Excessive number of categories: >7";
WarningsTable.warnings_str[["UNDERREP_CATS"]]="Under-represented categories";

WarningsTable.FindWarnings=function(ind_var_info){
	# Identify and assign warnings based on the variable information
	
	num_var_info_col=ncol(ind_var_info);
	
	# Get variable name first or else unlisting will convert everything to character()
	var_name=unlist(ind_var_info["VariableName"]);
	var_info=unlist(ind_var_info[2:num_var_info_col]);
	
	warnings_list_matrix=matrix(NA, ncol=3, nrow=0);
	colnames(warnings_list_matrix)=c("VariableName", "WarningCode", "WarningMessage");
	
	cat("Analyzing: ", var_name, "\n");
	
	if(!var_info["IsVarnameValid"]){
		warnings_list_matrix=rbind(warnings_list_matrix, c(var_name, "INV_VAR_NAME"));
	}
	
	if(var_info["PercentNAs"]>50){
		warnings_list_matrix=rbind(warnings_list_matrix, c(var_name, "NA_GT50"));
	}else if(var_info["PercentNAs"]>25){
		warnings_list_matrix=rbind(warnings_list_matrix, c(var_name, "NA_GT25"));
	}else if(var_info["PercentNAs"]>10){
		warnings_list_matrix=rbind(warnings_list_matrix, c(var_name, "NA_GT10"));
	}
	
	if(var_info["IsDichotomous"] & !var_info["IsBoolean"]){
		warnings_list_matrix=rbind(warnings_list_matrix, c(var_name, "DICH_NOTBOOL"));
	}
	
	if(var_info["NumberUnique"]==1){
		warnings_list_matrix=rbind(warnings_list_matrix, c(var_name, "ALL_IDENT"));
	}
	
	if(!is.na(var_info["ImprovedBySqrtTransform"]) & var_info["ImprovedBySqrtTransform"]){
		warnings_list_matrix=rbind(warnings_list_matrix, c(var_name, "REC_SQRT_TRANS"));
	} else if(!is.na(var_info["ImprovedByLogTransform"]) & var_info["ImprovedByLogTransform"]){
		warnings_list_matrix=rbind(warnings_list_matrix, c(var_name, "REC_LOG_TRANS"));
	} else if(!is.na(var_info["ImprovedByLogitTransform"]) & var_info["ImprovedByLogitTransform"]){
		warnings_list_matrix=rbind(warnings_list_matrix, c(var_name, "REC_LOGIT_TRANS"));
	}
	
	if(!is.na(var_info["IsPotentialDate"]) & var_info["IsPotentialDate"]){
		warnings_list_matrix=rbind(warnings_list_matrix, c(var_name, "NON_ISODATE"));
	}
	
	if(!is.na(var_info["AreLevelNamesValid"]) & !var_info["AreLevelNamesValid"]){
		warnings_list_matrix=rbind(warnings_list_matrix, c(var_name, "INV_LVL_NAMES"));
	}
	
	return(warnings_list_matrix);
}

#------------------------------------------------------------------------------

WarningsTable.GenerateFromVariableInfo=function(var_info_table){
	# Based on the variable information, identify warnings

	warnings_list_matrix=matrix(NA, ncol=3, nrow=0);
	num_variables=nrow(var_info_table);
	cat("Num Variables: ", num_variables, "\n");
	
	# Find the warning for each variable
	for(i in 1:num_variables){	
		warnings_list_matrix=rbind(
			warnings_list_matrix,
			WarningsTable.FindWarnings(var_info_table[i,])
		);
	}
	
	# Attach the full message to the warning table
	num_warnings=nrow(warnings_list_matrix);
	for(i in 1:num_warnings){
		warnings_list_matrix[i,"WarningMessage"]=WarningsTable.warnings_str[[warnings_list_matrix[i,"WarningCode"]]];
	}
	
	return(warnings_list_matrix);

}

###############################################################################

if(!exists("integration")){
	setwd("D:\\work_git\\AnalysisTools\\Metadata\\MetadataExplorer");
	test_mat=as.data.frame(read.table("TestExample.tsv", as.is=T, check.names=F, header=T, sep="\t"));
	
	VariableInfo.build(test_mat);
	print(VariableInfo);

	WarningsTable=WarningsTable.GenerateFromVariableInfo(VariableInfo);
	print(WarningsTable);

}

