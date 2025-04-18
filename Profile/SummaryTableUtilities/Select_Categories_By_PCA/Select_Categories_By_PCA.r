#!/usr/bin/env Rscript

library('getopt');

source("~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r");
source("~/git/AnalysisTools/Metadata/OutputFileLibrary/OutputFileLibrary.r");

###############################################################################

params=c(
	"summary_table", "s", 1, "character",
	"min_cuml_cov", "c", 2, "numeric",
	"min_indv_cov", "v", 2, "numeric",
	"output_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-s <input summary_table.tsv>\n",
	"	[-c <min cumulative coverage of variance, e.g. 0.90>]\n",
	"	[-v <min individual coverage of variance, e.g. 0.05>]\n",
	"	-o <output file root>\n",
	"\n",	
	"This script will read in the summary table, look for the Remaining\n",
	"column, perform the ALR if found, or compute the CLR. \n",
	"\n",
	"For each of the transformed columns, the correlation matrix will be\n",
	"calculated and PCA will be performed.\n",
	"\n",
	"The minimum number of PCs necessary to exceed the specified\n",
	"minimum coverage of variance, will be identified, and the closest\n",
	"variable will be selected.\n",
	"\n",
	"Selected variables will be exported.  The following files will be generated:\n",
	"	<output_file_root>.top_pc_prox.summary_table.tsv\n",
	"	<output_file_root>.top_pc_prox.report.pdf\n",
	"\n",
	"\n");

if(!length(opt$summary_table)){
	cat(usage);
	q(status=-1);
}

MinCumlCov=NULL;
if(length(opt$min_cuml_cov)){
	MinCumlCov=opt$min_cuml_cov;
}

MinIndvCov=NULL;
if(length(opt$min_indv_cov)){
	MinIndvCov=opt$min_indv_cov;
}


options(width=200);

###############################################################################

InputSummaryTable=opt$summary_table;
OutputRoot=opt$output_root;
cnt_adjmt=0.01;

cat("\n");
cat("Input Summary Table: ", InputSummaryTable, "\n");
cat("Minimum PC Cumulative Coverage: ", MinCumlCov, "\n");
cat("Minimum PC Individual Coverage: ", MinIndvCov, "\n");
cat("Output Root: ", OutputRoot, "\n");
cat("\n");

pdf(paste(OutputRoot, ".top_pc_prox.report.pdf", sep=""), height=8.5, width=11);

###############################################################################

calc_non_remaining_proportions=function(norm_mat){
	
	categories=colnames(norm_mat);
	rem_ix=(categories=="Remaining");
	if(any(rem_ix)){
		categories=setdiff(categories, "Remaining");
	}
	norm_mat=norm_mat[,categories,drop=F];
	non_rem_prop=apply(norm_mat, 1, sum);
	return(non_rem_prop);
	
}

###############################################################################

orig_counts=load_summary_file(InputSummaryTable);

sample_depths=apply(orig_counts, 1, sum);
cat("Normalizing...\n");
normalized=normalize(orig_counts);

non_rem_prop=calc_non_remaining_proportions(normalized);
original_mean_non_remaining_prop=mean(non_rem_prop);
cat("Mean Non-Remaining: ", original_mean_non_remaining_prop, "\n");

num_samples=nrow(normalized);
num_categories=ncol(normalized);

###############################################################################

ratio_trans=function(norm_mat){

        cat_names=colnames(norm_mat);
        num_col=ncol(norm_mat);
        remaining_ix=(cat_names=="Remaining");
        remaining_found=any(remaining_ix);


        if(remaining_found){
                # Do ALR, with Remaining as denominator
                cat("Performing Additive Ratio with Remaining...\n");

                # Split remaining from field
                norm_wo_rem=norm_mat[,!remaining_ix,drop=F];
                remainder=norm_mat[,remaining_ix];

                # Calc alr
                ratio_mat=norm_wo_rem / remainder;

        }else{
                # Do CLR, with geometric mean as denominator
                cat("Performing Centered Ratio with Geometric Mean...\n");

                geomean=function(x){
                        prod(x)^(1/length(x));
                }

                gm=apply(norm_mat, 1, geomean);

                # Calc clr
                ratio_mat=norm_mat / gm;
        }

        #lr_mat=log(ratio_mat);

        return(ratio_mat);

}

###############################################################################

plot_text(c(
	script_name,
	"",
	paste("Input Summary Table: ", InputSummaryTable),
	paste("Minimum PC Cuml Coverage: ", MinCumlCov),
	paste("Minimum PC Indv Coverage: ", MinIndvCov),
	paste("Output Root: ", OutputRoot)
	));

cat("Computing Ratios...\n");
ratio_tab=ratio_trans(normalized);
num_avail_var=ncol(ratio_tab);
#print(ratio_tab);

cat("Computing *Spearman* correlation coefficients on Ratios...\n");
correl_mat=cor(ratio_tab, method="spearman");

paint_matrix(correl_mat, "Spearman Correlation Matrix", plot_min=-1,
	plot_max=1, deci_pts=2, label_zeros=F, show_leading_zero=F)

#print(correl_mat);
cat("Calculating Eigen Values...\n");
eigen_rec=eigen(correl_mat);

pca_propvar=eigen_rec$values/sum(eigen_rec$values);
pca_propcumsum=cumsum(pca_propvar);

if(!is.null(MinCumlCov)){
	cat("Using Minimum Cumulative Coverage (", MinCumlCov, ")\n");
	num_pc_at_cutoff=sum(pca_propcumsum<MinCumlCov)+1;
}else if(!is.null(MinIndvCov)){
	cat("Using Minimum Individual Coverage (", MinIndvCov, ")\n");
	num_pc_at_cutoff=sum(pca_propvar>MinIndvCov);
}else{
	cat("Coverage Minimum not specified.\n");
	exit();
}

cumlsum_at_cutoff=pca_propcumsum[num_pc_at_cutoff];
cat("Num Selected PCs: ", num_pc_at_cutoff, "\n");
cat("Cuml Variance Selected: ", cumlsum_at_cutoff, "\n");


ratio_tab_as_rank=apply(ratio_tab, 2, rank);
ratio_tab_rank_scaled=scale(ratio_tab_as_rank, center=T, scale=T);

# Calculate the ordination positions
scores=(ratio_tab_rank_scaled %*% eigen_rec$vectors);

print(scores);
if(ncol(scores)==1){
	scores=cbind(scores,0);
}

par(mar=c(5,5,7,1));
plot(scores[,1], scores[,2], 
	xlab=paste("PC01:", round(pca_propvar[1]*100, 2), "%"),
	ylab=paste("PC02:", round(pca_propvar[2]*100, 2), "%"),
	main="Top 2 PCs", asp=1
	);


#cat("PC Proportion of Variance (PDF):\n");
#print(pca_propvar);
#cat("PC Cumulative Variance (CDF):\n");
#print(pca_propcumsum);

pc_table=cbind(pca_propvar, pca_propcumsum)
colnames(pc_table)=c("Prop", "Cuml");

pc_table_rounded=apply(pc_table, 1:2, function(x){round(x,5);});
rownames(pc_table_rounded)=paste("PC", sprintf("%02i", 1:nrow(pc_table_rounded)), sep="");
print(pc_table_rounded);
plot_text(c(
	"PCA Variance Allocation:",
	"",
	capture.output(print(pc_table_rounded, quote=F))
	));

#------------------------------------------------------------------------------

par(mar=c(6,6,4,1));
mids=barplot(pc_table[,"Prop"],
	yaxt="n", xaxt="n",
	ylim=c(0,1.1), col="green",
	ylab="Proportion of Variance");
points(mids, pc_table[,"Cuml"], type="b", col="blue");

abline(v=mids[num_pc_at_cutoff], col="red", lwd=3, lty="dashed");
abline(h=pc_table[num_pc_at_cutoff,"Cuml"], col="red", lwd=3, lty="dashed");
points(mids[num_pc_at_cutoff], pc_table[num_pc_at_cutoff,"Cuml"], col="red", pch=1, cex=3);

abline(h=c(MinIndvCov, MinCumlCov), col="black", lwd=1,lty="dotted");

yticks=sort(c(MinIndvCov, MinCumlCov, 0,.25, .33, .5, .66, .75, .9, .95, 1));
axis(side=2, at=yticks, labels=yticks, las=2);
axis(side=1, at=mids, rownames(pc_table_rounded), las=2);

#------------------------------------------------------------------------------

cat("\n\n");
cat("Num PCs Kept: ", num_pc_at_cutoff, "\n");
cat("\n\n");

pc_to_proxy_table=matrix(character(), nrow=num_pc_at_cutoff, ncol=2);
rownames(pc_to_proxy_table)=paste("PC", sprintf("%02i", 1:num_pc_at_cutoff), sep="");
colnames(pc_to_proxy_table)=c("ClosestProxy", "Correl");

varnames=colnames(ratio_tab);
for(pc_ix in 1:num_pc_at_cutoff){

	cat("\n\n");
	cat("Looking for proxy for PC: ", pc_ix, "\n");
	cur_pc_val=scores[,pc_ix];

	correl_w_ratios=as.vector(cor(cur_pc_val, ratio_tab, method="spearman"));
	mag_correl=abs(correl_w_ratios);	

	correl_mat=cbind(correl_w_ratios, mag_correl);
	rownames(correl_mat)=colnames(ratio_tab);
	print(correl_mat);
	
	max_cor_mag=max(mag_correl);
	cat("Max Magnitude of Cor\n");
	print(max_cor_mag);

	max_cor_ix=min(which(max_cor_mag==mag_correl));
	cat("Most Cor Variable:\n");
	print(max_cor_ix);
	best_proxy=varnames[max_cor_ix];

	pc_to_proxy_table[pc_ix, ]=c(
		best_proxy, 
		sprintf("%7.4f", correl_w_ratios[max_cor_ix])
		);


	#----------------------------------------------------------------------

	plot_text(c(
		paste("PC:", pc_ix, sep=""),
		"",
		capture.output(print(correl_mat, quote=F)),
		"",
		paste("Greatest Magnitude of Correlation: ", max_cor_mag),
		paste("(Index: ", max_cor_ix, ")"),
		paste("Best Proxy: ", best_proxy)
		));

}

cat("\nSelected Variables\n");
print(pc_to_proxy_table);

non_redund_variables=unique(pc_to_proxy_table[,"ClosestProxy"]);
num_nrvar=length(non_redund_variables);
cat("\n");
cat("Non-Redundant Selected Variables: \n");
print(non_redund_variables);

###############################################################################

selected_var_norm_mat=normalized[,non_redund_variables,drop=F];
selected_non_rem_prop=mean(selected_var_norm_mat);
excluded_prop=original_mean_non_remaining_prop-selected_non_rem_prop;
cat("\n");
cat("Mean Original Non-Remaining Abundance: ", original_mean_non_remaining_prop, "\n");
cat("Mean Selected Non-Remaining Abundance: ", selected_non_rem_prop, "\n");
cat("  Excluded: ", excluded_prop, "\n");


pc_table_formatted=apply(pc_table[1:num_pc_at_cutoff,,drop=F], 1:2, function(x){sprintf("%7.5f", x)});
plot_text(c(
	"Selected Proxies for each PC",
	"",
	capture.output(print(cbind(pc_table_formatted, pc_to_proxy_table), quote=F)),
	"",
	"",
	"Non-Redundant Selected Variables:",
	"",
	paste(1:num_nrvar, ". ", non_redund_variables, sep=""),
	"",
	"",
	paste("Num Available Variables: ", num_avail_var, sep=""),
	paste("   ", num_nrvar, " of ", num_avail_var, " kept.", sep=""),
	paste("   ", signif(num_nrvar/num_avail_var, 5), " proportion retained.", sep=""),
	"",
	"",
	paste("Original Non-Remaining Abundance (Mean): ", signif(original_mean_non_remaining_prop,5)),
	paste("Selected Non-Remaining Abundance (Mean): ", signif(selected_non_rem_prop,5)),
	paste("   Excluded Abundance (Mean): ", signif(excluded_prop,5)),
	paste("   Proportion Retained: ", 
		signif(selected_non_rem_prop/original_mean_non_remaining_prop,5))
	));

#out_fn=paste(OutputRoot, ".pca_select.pdf", sep="");
#pdf(out_fn, height=8.5, width=11);

cat("\n");
cat("Exporting selected counts (PC Proxies)...\n");
selected_counts=orig_counts[,non_redund_variables,drop=F];
selected_depths=apply(selected_counts,1,sum);
recalc_remaining=sample_depths-selected_depths;
output_counts=cbind(selected_counts, recalc_remaining);
colnames(output_counts)=c(non_redund_variables, "Remaining");

outfn=paste(OutputRoot, ".top_pc_prox.summary_table.tsv", sep="");
write_summary_file(output_counts, outfn);

###############################################################################
# Generate before rank abundance plot

mean_norm=apply(normalized, 2, mean);
par(mar=c(15,3,2,10));
plot_rank_abundance(mean_norm, num_disp_max=50, title="Original Distribution");

# Generate after rank abundance plot
norm_selected=normalize(output_counts);
mean_norm_selected=apply(norm_selected, 2, mean);
plot_rank_abundance(mean_norm_selected, num_disp_max=50, title="Selected Distribution");

plot_rank_abundance(mean_norm, num_disp_max=50, 
	title="Original Distribution", exclude_Remaining=T);
plot_rank_abundance(mean_norm_selected, num_disp_max=50, 
	title="Selected Distribution", exclude_Remaining=T);

###############################################################################
# Output stats

cat("\n");
cat("Exporting selection statistics...\n");

select_stats=c(
	"NumSelected", "NumAvail", "PropVarRetained",
	"NumBasePCs", "BasePCCumlVar", 
	"SelectedNonRemAbund", "OrigNonRemAbund", "PropAbundRetained"
	);

num_select_stats=length(select_stats);
out_stats=matrix(0, nrow=1, ncol=num_select_stats);
colnames(out_stats)=select_stats;
rownames(out_stats)=tail(strsplit(InputSummaryTable, "/")[[1]], 1);

out_stats[1,]=c(
	num_nrvar, num_avail_var, signif(num_nrvar/num_avail_var, 5),
	num_pc_at_cutoff, signif(cumlsum_at_cutoff, 5),
	signif(selected_non_rem_prop, 5), signif(original_mean_non_remaining_prop,5),
	signif(selected_non_rem_prop/original_mean_non_remaining_prop,5)
	);

out_mat=cbind(rownames(out_stats), out_stats);
colnames(out_mat)=c("InputFile", colnames(out_stats));

outfn=paste(OutputRoot, ".top_pc_prox.select_stats.tsv", sep="");
write.table(out_mat, outfn, sep="\t", row.names=F, quote=F);


###############################################################################

cat("Done.\n")
print(warnings());

q(status=0)
