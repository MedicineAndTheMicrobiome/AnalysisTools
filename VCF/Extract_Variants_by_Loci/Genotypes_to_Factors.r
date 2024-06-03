#!/usr/bin/env Rscript

###############################################################################

library(getopt);
library(parallel);

params=c(
	"genotype_dir_name", "d", 1, "character",
	"output_fn_root", "o", 1, "character",
	"qual_cutoff", "q", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

DEF_QUAL_CUTOFF=30;

options(width=129);

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-d <directory containing all the *.genotype.tsv files>\n",
	"	-o <output file root>\n",
	"\n",
	"This script will read in a list of genotype files and generate a factor\n",	
	"file, which can be used in subsequent analyses depending on a\n",
	"subject_id x features (rows x columns) matrix.\n",
	"\n",
	"	1.) Read in all the genotype files.\n",
	"	2.) For each locus (chrom, pos), accumulate allele frequencies\n",
	"		and confirm that the allele of the reference is identical.\n",
	"	  5.) Plot the variant frequencies vs. quality\n", 
	"	  6.) For each alternate allele do:\n",
	"	    7.) Calculate the hardy-weinberg equilibrium stats.\n",
	"	    8.) Create a 'column': <chrom>_<pos>_<ref_al>_<var_al>,\n",
	"		with values: 0,.5,1, depending hom/het/hom of each subject\n",
	"\n",
	"\n",	
	"An example of a genotype file is: \n",
	"\n",
	"Chromosome  Pos        Ref Var VarFreq Qual    MeanRegionQual  NumRegionVariants\n",
	"chr2        143551211  C   G   1       701.77  471.77          6\n",
	"chr1        234217153          0       NaN     NaN             0\n",
	"chr3        15461174           0       340.88  340.88          9\n",
	"chr19       38266675   CA  C   0.5     75.73   272.16          5\n",
	"chr1        219121228          0       776.76  776.76          3\n",
	"...\n",
	"\n",
	"\n");

if(
	!length(opt$genotype_dir_name) || 
	!length(opt$output_fn_root)
){
	cat(usage);
	q(status=-1);
}

GenotypeDir=opt$genotype_dir_name;
OutputFNameRoot=opt$output_fn_root;


QualCutoff=DEF_QUAL_CUTOFF;
if(length(opt$qual_cutoff)){
	QualCutoff=opt$qual_cutoff;
}

cat("\n");
cat("Input Genotypes Dir: ", GenotypeDir, "\n");
cat("Output Filename Root: ", OutputFNameRoot, "\n");
cat("Quality Cutoff: ", QualCutoff, "\n");
cat("\n");

##############################################################################
# Load Genotypes Files
##############################################################################

load_genotype_file=function(filepath){
	cat("Loading: ", filepath, "\n");
	geno_tab=read.table(filepath, header=T, sep="\t", as.is=T, 
		stringsAsFactors=F, colClasses="character", na.strings=NULL);

	# For some reason, if all the column values are "", the values are read
	# in as NA, instead of "".

	if(all(is.na(geno_tab[,"Ref"]))){
		geno_tab[,"Ref"]="";
	}

	if(all(is.na(geno_tab[,"Var"]))){
		geno_tab[,"Var"]="";
	}

	rownames(geno_tab)=paste(geno_tab[,"Chromosome"], "_", geno_tab[,"Pos"], sep="");
	return(geno_tab);
}

#-----------------------------------------------------------------------------

cat("Loading Genotype Files from: ", GenotypeDir, "\n");
genotype_fn_arr=list.files(GenotypeDir, pattern="*.genotype.tsv");

num_genotype_files=length(genotype_fn_arr);
#print(genotype_fn_arr);
cat("Found ", num_genotype_files, " genotype.tsv files in : ", GenotypeDir, "\n");
cat("\n");

genotypes_list=list();
for(genotype_path in genotype_fn_arr){
	path_comp=strsplit(genotype_path, "/")[[1]];	
	genotype_fn=path_comp[length(path_comp)];
	genotype_fn=gsub("\\.genotype\\.tsv$", "", genotype_fn);
	genotypes_list[[genotype_fn]]=
		load_genotype_file(paste(GenotypeDir, "/", genotype_path, sep=""));

	
}

num_variants_list=lapply(genotypes_list, nrow)
num_variants_arr=unlist(num_variants_list);
minmax_num_var=range(num_variants_arr);
if(minmax_num_var[1]!=minmax_num_var[2]){
	cat("\n");
	cat("Error:  The number of variants is not the same across all the\n");
	cat("genotype files.  There should be a placeholder for all reference\n");
	cat("alleles.\n");
	cat("\n");

	cat("The following genotype files (subjects) have fewer than ", 
		minmax_num_var[2], " loci.\n", sep="");
	print(num_variants_arr[num_variants_arr<minmax_num_var[2]]);
	
	quit(status=-1);
}
num_variants=minmax_num_var[2];

cat("Number of rows/loci to consolidate: ", num_variants, "\n");
cat("\n");

##############################################################################

calculate_chisq_HWE=function(var_freq, epsilon=0.001){
	# Given list of 0, .5, and 1's, will calculate
	# Hardy Weinberg Distribution and chi^2 test, and use epsilon
	# for the continuity correction (0's for Exp will give Inf)

	num_subj=length(var_freq);
	
	if(num_subj==0){
		stats=c(1, "All Reference", rep(0,6));
	}else{
		tot_var_freq=sum(var_freq);
		var_prop=tot_var_freq/num_subj;
		ref_prop=1-var_prop;

		# p = prop(Reference)
		# q = prop(Variant)

		p_sqrd=ref_prop^2;
		two_pq=2*var_prop*ref_prop;
		q_sqrd=var_prop^2;
		expected=(c(p_sqrd, two_pq, q_sqrd)*num_subj)+epsilon;

		obs_p_sqrd=sum(var_freq==0);
		obs_two_pq=sum(var_freq==0.5);
		obs_q_sqrd=sum(var_freq==1);
		observed=c(obs_p_sqrd, obs_two_pq, obs_q_sqrd)+epsilon;

		cat("Expected:\n");
		print(expected);
		cat("Observed:\n");
		print(observed);

		chisq_val=sum(((expected-observed)^2)/expected);
		pval=1-pchisq(chisq_val, df=1);

		msg="";
		if(pval<0.01){
			if(observed[1]>expected[1]){
				msg="Excess Homo Reference";
			}else if(observed[2]>expected[2]){
				msg="Excess Heterozygous";
			}else if(observed[3]>expected[3]){
				msg="Excess Homo Variants";
			}
		}else{
			msg="Passed HWE Test";
		}

		stats=c(round(pval, 3), msg, round(observed,2), round(expected,2));
	}

	names(stats)=c("pval", "msg",
			"obs_psqrd", "obs_2pq", "obs_qsqrd",
			"exp_psqrd", "exp_2pq", "exp_qsqrd"
			);

	return(stats);
}

#------------------------------------------------------------------------------

create_variant_focused_table=function(lim, allele){
	# Given a genotype table with multiple alleles, it will extract a
	# new table, that is focused on the specified var/allele.  The
	# resultant table is still with respect to the reference, but
	# but any alleles not matching the specified allele will be set
	# to zero.

	num_rows=nrow(lim);
	for(i in 1:num_rows){
		if(lim[i, "VarFreq"]=="0"){
			next;
		}else{
			variants=strsplit(lim[i, "Var"], ",")[[1]];
			freqs=strsplit(lim[i, "VarFreq"], ",")[[1]];

			lim[i,"Var"]="";
			lim[i,"VarFreq"]="0";

			# Copy over variant into if it matches our specified allele
			num_variants=length(variants);
			for(vix in 1:num_variants){
				if(variants[vix]==allele){
					lim[i, "Var"]=variants[vix];
					lim[i, "VarFreq"]=freqs[vix];
				}
			}
		}
	}

	return(lim);
}

#------------------------------------------------------------------------------

split_table_multiallele=function(lim){
	# This function takes a genotype table and splits it into multiple tables
	# one for each allele type.  

	# Look for multi-allele variants
	variants=c();
	num_rows=nrow(lim);
	for(i in 1:num_rows){
		variants=c(variants, strsplit(lim[i,"Var"], ",")[[1]]);
	}

	variants=setdiff(unique(variants), c(NA, ""));
	num_variants=length(variants);

	cat("Variants found at locus:\n");
	print(variants);
	
	if(length(variants)==0){
		variants="";
	}

	table_list=list();
	if(num_variants<=1){
		# No variants or single allele variant
		table_list[[variants]]=lim;
	}else{
		cat("Multiple Variant Alleles Found.\n");
		for(var in variants){
			cat("  Extracting out: ", var, "\n");
			table_list[[var]]=create_variant_focused_table(lim, var);
		}
	}
	return(table_list);
}

#------------------------------------------------------------------------------

genotype_names=names(genotypes_list);
loci_ids=rownames(genotypes_list[[1]]);

cat("Identified Loci IDs:\n");
print(loci_ids);

loci_list=list();
loci_info=list();

for(cur_loci_id in loci_ids){
	cat("\n\n");
	cat("Working on: ", cur_loci_id, "\n");
	
	locus_tab_header= c("Ref", "Var", "VarFreq", "Qual", "MeanRegionQual", "NumRegionVariants");
	locus_info_matrix=matrix(NA, nrow=num_genotype_files, ncol=length(locus_tab_header));
	colnames(locus_info_matrix)=locus_tab_header;
	rownames(locus_info_matrix)=genotype_names;

	for(cur_genotype_name in genotype_names){
		gt=unlist(genotypes_list[[cur_genotype_name]][cur_loci_id, locus_tab_header]);
		locus_info_matrix[cur_genotype_name, locus_tab_header]= gt[locus_tab_header];
	}	

	median_regional_variants=median(as.numeric(locus_info_matrix[,"NumRegionVariants"]),na.rm=T);
	median_variant_qual=median(as.numeric(locus_info_matrix[,"Qual"]), na.rm=T);

	print(locus_info_matrix);
	cat("Median Num Regional Variants: ", median_regional_variants, "\n");
	cat("Median Variant Quality: ", median_variant_qual, "\n");
	
	# Confirm all Ref alleles are the same:
	consistent_ref=T;
	ref_all=setdiff(unique(locus_info_matrix[,"Ref"]), c("", NA));
	if(length(ref_all)<=1){
		cat("Unique Reference Allele: ", ref_all, "\n");
	}else{
		cat("Error.  Multiple Reference Alleles...\n");
		consistent_ref=F;
		print(ref_all);

		info=c(paste(ref_all, collapse=","), "",
			NA, "Inconsistent Reference",
			c(NA, NA, NA), c(NA, NA, NA),
			consistent_ref,
			median_variant_qual,median_regional_variants);

		names(info)=c("pval", "msg",
			"obs_psqrd", "obs_2pq", "obs_qsqrd",
			"exp_psqrd", "exp_2pq", "exp_qsqrd"
			);
			
		
		loci_info[[cur_loci_id]]=info;
		next;
	}

	# Split table if more than 1 variant allele
	loc_inf_mat_list=split_table_multiallele(locus_info_matrix);
	#print(loc_inf_mat_list);

	for(var_allele in names(loc_inf_mat_list)){
		# If there are not variants, then we don't have any alleles.
		# We don't even know what the reference is.

		loc_var_names=paste(
			cur_loci_id, "_", ref_all, "_", var_allele, sep="");
		cat("Extracting: ", loc_var_names, "\n");

		cur_loc_inf_mat=loc_inf_mat_list[[var_allele]];

		var_freq=as.numeric(cur_loc_inf_mat[,"VarFreq"]);
		names(var_freq)=rownames(cur_loc_inf_mat);

		#cat("Variant Frequencies:\n");
		#print(var_freq);

		hwe_stats=calculate_chisq_HWE(var_freq);
		cat("Hardy Weinberg Equilibrium Statistics:\n");
		print(hwe_stats);

		cat("P-value: ", hwe_stats["pval"], "\n");
		if(hwe_stats["pval"]<0.1){
			cat("Locus fails Hardy-Weinberg Equilibrium Chi-squared Test.\n");
		}

		if(length(var_freq)){
			loci_list[[loc_var_names]]=var_freq;
		}

		if(length(ref_all)==0){
			ref_all="";
		}

		info=c(ref_all, var_allele, hwe_stats, consistent_ref, median_variant_qual, 
			median_regional_variants);
		info_names=c("ref", "var", names(hwe_stats), "cons_ref", 
			"med_var_qual", "med_num_reg_var");

		names(info)=info_names;
		loci_info[[loc_var_names]]=info;
	}
}

##############################################################################
# Write loci as factor file

#print(loci_list);

loci_names=names(loci_list);
num_variants=length(loci_list);

subject_ids=names(genotypes_list);
num_subjects=length(subject_ids);

#print(subject_ids);

subj_genotype_tab=matrix(NA, nrow=num_subjects, ncol=num_variants);
rownames(subj_genotype_tab)=subject_ids;
colnames(subj_genotype_tab)=loci_names;

for(i in 1:num_variants){
	subj_genotype_tab[subject_ids,i]=loci_list[[i]][subject_ids];
}

##############################################################################
# Write loci information 

#print(loci_info);
num_variants=length(loci_info);
variant_names=names(loci_info);
num_info_fields=length(loci_info[[1]]);
info_field_names=names(loci_info[[1]]);
loci_info_tab=matrix(NA, nrow=num_variants, ncol=num_info_fields);
rownames(loci_info_tab)=variant_names;
colnames(loci_info_tab)=info_field_names;

for(i in 1:num_variants){
	loci_info_tab[i,]=loci_info[[i]];
}

# Calculate HWE
exp_psqrd=as.numeric(loci_info_tab[,"exp_psqrd"]);
exp_2pq=as.numeric(loci_info_tab[,"exp_2pq"]);
exp_qsqrd=as.numeric(loci_info_tab[,"exp_qsqrd"]);
exp_variants=exp_2pq+exp_qsqrd;

any_var_exp_prob=round(exp_variants/(exp_variants+exp_psqrd), 4);
names(any_var_exp_prob)=variant_names;

loci_info_tab=cbind(variant_names, loci_info_tab, any_var_exp_prob);
colnames(loci_info_tab)=c("locus_id", info_field_names, "any_var_exp_prob");
rownames(loci_info_tab)=variant_names;

write.table(loci_info_tab, file=paste(OutputFNameRoot, ".loci_report.tsv", sep=""),
	quote=F, sep="\t", row.names=F, col.names=T);

##############################################################################
##############################################################################
# Estimate the missed variant detection probability for each subject

subj_ids=rownames(subj_genotype_tab);
exported_variant_names=colnames(subj_genotype_tab);

cat("Exported Variants, Any Prob:\n");
exported_any_var_prob=any_var_exp_prob[exported_variant_names];
print(exported_any_var_prob);

estimate_pr_detect=function(obs_det, exp_det){
	# Given observed detection for a subject (obs_det) and
	# the expected detection rate for the cohort, estimate
	# the probability of detection.  In other words,
	# Estimate Pr( Ind Det[i] | Obs Det[i], Exp Pr(Det))
	# Note that the expected pr(detect)= prob(hets) + prob(homo variant)
	# Similarly, the obs detection includes both hets and homo variant.

	loglik=function(pr_detect, detect_loci, detect_freq){
		pr_obs=pr_detect*detect_freq;
		ll=sum(detect_loci*log(pr_obs) + (1-detect_loci)*log(1-pr_obs));
		cat("Param: Pr(Det)=", pr_detect, " LogLik=", ll, "\n");
		return(ll);
	}

	optim_nll=function(param){
		return(-loglik(param, obs_det, exp_det));
	}

	pr_detect=0.9;
	opt_res=optim(pr_detect, fn=optim_nll, method="L-BFGS-B", lower=1e-300, upper=1)
	#print(opt_res);	
	return(opt_res);
}

subj_detect_mat=apply(subj_genotype_tab, c(1,2), 
	function(x){as.numeric((as.numeric(x)>0))}
	);

detect_rates=list();
for(subj in subj_ids){
	cat("Working on: ", subj, "\n");
	cur_detects=subj_detect_mat[subj,];	
	
	#cat("Observed:\n");
	#print(cur_detects);
	#cat("Cohort Probabilities:\n");
	#print(exported_any_var_prob);

	pr_det=estimate_pr_detect(obs_det=cur_detects, exp_det=exported_any_var_prob);
	detect_rates[[subj]]=pr_det$par;

	cat("\n");
}
detect_rate_arr=unlist(detect_rates);

pdf(paste(OutputFNameRoot, ".detect_rates.hist.pdf", sep=""), height=11, width=8.6);
hist_res=hist(detect_rate_arr, breaks=20);

#------------------------------------------------------------------------------
# Output Estimated Pr(Detect) for each VCF/Subject

vcfname=names(detect_rate_arr);
sort_by_vcfname_ix=order(vcfname);
sort_by_detect_rate_ix=order(detect_rate_arr);

outmat_by_vcfname=cbind(vcfname[sort_by_vcfname_ix], 
	sprintf("%3.4f", detect_rate_arr[sort_by_vcfname_ix]));

outmat_by_detrat=cbind(vcfname[sort_by_detect_rate_ix], 
	sprintf("%3.4f", detect_rate_arr[sort_by_detect_rate_ix]));

colnames(outmat_by_vcfname)=c("VCFName", "DetectRate");
colnames(outmat_by_detrat)=c("VCFName", "DetectRate");

write.table(outmat_by_vcfname, file=paste(OutputFNameRoot, ".detect_rates.vcfname_order.tsv", sep=""),
	quote=F, sep="\t", row.names=F, col.names=T);

write.table(outmat_by_detrat, file=paste(OutputFNameRoot, ".detect_rates.detrat_order.tsv", sep=""),
	quote=F, sep="\t", row.names=F, col.names=T);

#------------------------------------------------------------------------------
# Output distribution of Estimated Pr(Detect) 

num_subj=sum(hist_res$counts);
stats=cbind(hist_res$mids, hist_res$counts, 
	cumsum(hist_res$counts),
	sprintf("%3.4f", (hist_res$counts/num_subj)), 
	sprintf("%3.4f", cumsum(hist_res$counts/num_subj)));
colnames(stats)=c("PrDet_BinMids", "Num_VCFs", "CumulSum", "Prop_VCFs", "CDF_VCFs");
write.table(stats, file=paste(OutputFNameRoot, ".detect_rates_distr_tab.tsv", sep=""),
	quote=F, sep="\t", row.names=F, col.names=T);

#------------------------------------------------------------------------------
# Export snp table with subjects filtered out

for(cutoff in c(0.00, 0.60, 0.75, .90)){
	
	cat("Exporting at cutoff: ", cutoff, "\n");

	if(cutoff>0.00){
		cutoff_str=sprintf("%.2f", cutoff);
	}else{
		cutoff_str="all";
	}
	detect_rate_arr_gte=detect_rate_arr[detect_rate_arr>=cutoff];
	detect_rate_arr_excluded=detect_rate_arr[detect_rate_arr<cutoff];
	kept_vcfs=names(detect_rate_arr_gte);

	# Write snp tsv
	out_tab=cbind(
			kept_vcfs,
			detect_rate_arr_gte[kept_vcfs],
			subj_genotype_tab[kept_vcfs,,drop=F]
		);
	colnames(out_tab)=c("VCF_ID", "DetectRate", colnames(subj_genotype_tab));
	
	write.table(
		out_tab,
		file=paste(OutputFNameRoot, ".", cutoff_str, ".snps.tsv", sep=""), quote=F,
		sep="\t", row.names=F, col.names=T);

	# Write exclusion list
	order_ix=order(detect_rate_arr_excluded);
	detect_rate_arr_excluded=detect_rate_arr_excluded[order_ix];
	excl_tab=cbind(names(detect_rate_arr_excluded), detect_rate_arr_excluded);
	colnames(excl_tab)=c("vcf_name", "detection_rate");
	write.table(excl_tab, 
		file=paste(OutputFNameRoot, ".", cutoff_str, ".exclusion.lst", sep=""), quote=F,
		sep="\t", row.names=F, col.names=T);

}

##############################################################################

cat("\nDone.\n");

##############################################################################

print(warnings());
q(status=0);
