
test_and_apply_normalizing_transforms=function(mat_val, pval_cutoff=.2, plot_before_after=T){

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

		log_transformed=F;
		sqrt_transformed=F;
		negrecip_transformed=F;

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
			test_log_res=NULL;
			test_sqrt_res=NULL;
			test_negrecip_res=NULL;

			if(test_res$p.value<=pval_cutoff && num_unique_val>2){
				cat(" Not normal: ", test_res$p.value, "\n");

				# Apply transform
				if(any(values_nona==0)){
					log_values=log(values+1);
				}else{
					log_values=log(values);
				}
				sqrt_values=sqrt(values);
				negrecip_values=-1/values;

				# Test transforms for normality
				test_log_res=shapiro.test(log_values);
				test_sqrt_res=shapiro.test(sqrt_values);
				test_negrecip_res=shapiro.test(negrecip_values);

				# Look for transformation least not-normal
				test_pvals=c(
					test_res$p.value, 
					test_log_res$p.value,
					test_sqrt_res$p.value,
					test_negrecip_res$p.value
					);

				max_pvalue=max(test_pvals, rm.na=T);
				best_trans=which(test_pvals==max_pvalue);
			
				# Store results of best transform
				if(best_trans==1){
					# Keep original
					cat("  No Improvement: ", test_log_res$p.value, "\n");
					new_colnames=c(new_colnames, paste("orig_", var, sep=""));
				}else{

					cat("     Log p-val : ", test_log_res$p.value, "\n");
					cat("    Sqrt p-val : ", test_sqrt_res$p.value, "\n");
					cat("negrecip p-val : ", test_negrecip_res$p.value, "\n");
					
					if(best_trans==2){
						# Keep log transformed
						cat("  Log Transformation Effective: ", 
							test_log_res$p.value, "\n");
						new_colnames=c(new_colnames, paste("log_", var, sep=""));
						trans_mat[, var]=log_values;
						log_transformed=T;		
					}else if(best_trans==3){
						# Keep sqrt transformed
						cat("  Sqrt Transformation Effective: ", 
							test_sqrt_res$p.value, "\n");
						new_colnames=c(new_colnames, paste("sqrt_", var, sep=""));
						trans_mat[, var]=sqrt_values;
						sqrt_transformed=T;		
					}else if(best_trans==4){
						# Keep negrecip transformed
						cat("  NegRecip Transformation Effective: ", 
							test_sqrt_res$p.value, "\n");
						new_colnames=c(new_colnames, paste("negrecip_", var, sep=""));
						trans_mat[, var]=negrecip_values;
						negrecip_transformed=T;
					}
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
			
			if(log_transformed){
				hist(log_values, breaks=nclass, main=paste("log(", var,")", sep=""));
				title(main=sprintf("p-value: %4.4f", test_log_res$p.value), cex.main=.8, line=.5);
			}else if(sqrt_transformed){
				hist(sqrt_values, breaks=nclass, main=paste("sqrt(", var,")", sep=""));
				title(main=sprintf("p-value: %4.4f", test_sqrt_res$p.value), cex.main=.8, line=.5);
			}else if(negrecip_transformed){
				hist(negrecip_values, breaks=nclass, main=paste("-1/(", var,")", sep=""));
				title(main=sprintf("p-value: %4.4f", test_negrecip_res$p.value), cex.main=.8, line=.5);
			}else{
				plot(0,0, xlab="", ylab="", main="", xaxt="n", yaxt="n", bty="n", type="n");

				if(test_res$p.value>pval_cutoff){
					text(0,0, "Transform not necessary");
				}else{
					text(0,0, "Transform not beneficial");
				}
			}
		}

	}

	colnames(trans_mat)=new_colnames;

	trans_mat=trans_mat[,setdiff(new_colnames, delete_list),drop=F];

	if(plot_before_after){
		par(orig_par);
	}

}
