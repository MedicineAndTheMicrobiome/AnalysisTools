#!/usr/bin/env Rscript

library('getopt');

source("~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r");

params=c(
	"summary_table", "s", 1, "character",
	"fba_xml_file", "x", 1, "character",
	"override_file", "r", 1, "character",
	"output_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-s <input summary_table.tsv>\n",
	"	-x <FBA model list, Genus\\tXML>\n",
	"	-r <Override FBA model list>\n",
	"	-o <output file root>\n",
	"\n",	
	"This script will load the summary table, then identify the top\n",
	"taxa that need to be included so that every sample is covered\n",
	"by a minimum threshold, or has at least a minimum abundance.\n",
	"\n",
	"The 'Override FBA model list' is the same format as the -x option\n",
	"The purpose is to remap some taxa to specific xml files in case\n",
	"the taxa (genus) doesn't have a mapping in the model list.\n",
	"\n",
	"The output are files at various cutoffs:\n",
	"	<taxa name> \\t <xml file> \\n\n",
	"\n",
	"\n");

if(
	!length(opt$summary_table) ||
	!length(opt$fba_xml_file) ||
	!length(opt$override_file) ||
	!length(opt$output_root)
	){
	cat(usage);
	q(status=-1);
}

options(width=200);

###############################################################################

InputSummaryTable=opt$summary_table;
FBAModelList=opt$fba_xml_file;
OverrideFBAModelList=opt$override_file;
OutputRoot=opt$output_root;

cat("\n");
cat("Input Summary Table: ", InputSummaryTable, "\n");
cat("FBA Model List: ", FBAModelList, "\n");
cat("Override List File: ", OverrideFBAModelList, "\n");
cat("Output Root: ", OutputRoot, "\n");
cat("\n");

pdf(paste(OutputRoot, ".coverage.fba_model_select.pdf", sep=""), height=11, width=8.5);

###############################################################################

remove_categories=function(targets, mat){

	num_targets=length(targets);

	all_cat=colnames(mat);
	found=intersect(all_cat, targets);
	cat("Found:\n");
	print(found);
	clean_cat=setdiff(all_cat, found);
	clean_mat=mat[,clean_cat,drop=F];

	return(clean_mat);
}

remove_zero_count_samples=function(mat){
	sums=apply(mat, 1, sum);
	nz_ix=sums>0;
	num_zero=sum(sums==0);
	cat("Num Zero Count Samples: ", num_zero, "\n");
	return(mat[nz_ix,,drop=F]);
}

load_mappings=function(existing_mapping, map_fn){

	dat=read.table(map_fn, header=F, stringsAsFactors=F);
	map=as.list(dat[,2]);
	names(map)=dat[,1];

	if(!is.null(existing_mapping)){
		map=modifyList(existing_mapping, map);
	}

	return(map);
}

lookup=function(keys, map){
	out=character(length(keys));
	names(out)=keys;
	for(k in keys){
		if(is.null(map[[k]])){
			out[k]=NA;	
		}else{
			out[k]=map[[k]];	
		}
	}
	return(out);

}

###############################################################################

cat("Loading Mapping:\n");
model_mapping=load_mappings(NULL, FBAModelList);
#print(model_mapping);

cat("Loading Overrides:\n");
model_mapping=load_mappings(model_mapping, OverrideFBAModelList);

###############################################################################

# Load Summary Table
orig_counts=load_summary_file(InputSummaryTable);
num_samples=nrow(orig_counts);
num_categories=ncol(orig_counts);

cat("Num Categories: ", num_categories, "\n", sep="");
cat("Num Samples: ", num_samples, "\n", sep="");

#------------------------------------------------------------------------------
# Filter categories and samples

orig_catnames=colnames(orig_counts);

# Remove virsues
virus_idx=grep("virus$", orig_catnames);
virus_names=orig_catnames[virus_idx];

cat("Viruses found:\n");
print(virus_names);

# Remove other abundant non-relevent or unwanted taxa
other_taxa=c("Homo", "Benincasa",  "Borrelia");

clean_counts=remove_categories(c(virus_names, other_taxa), orig_counts);
clean_counts=remove_zero_count_samples(clean_counts);

num_samples=nrow(clean_counts);
num_categories=ncol(clean_counts);
cat("Cleaned Num Categories: ", num_categories, "\n", sep="");
cat("Cleaned Num Samples: ", num_samples, "\n", sep="");

sample_depths=apply(clean_counts, 1, sum);
par(mfrow=c(2, 1));
hist(sample_depths, main="Sample Depths");
hist(log(sample_depths), main="Log(Sample Depths)");

#------------------------------------------------------------------------------
# Normalize and sort

cat("Normalizing...\n");
#print(normalize);
normalized=normalize(clean_counts);

mean_normalized=apply(normalized, 2, mean);
#print(mean_normalized);
order_ix=order(mean_normalized, decreasing=T);

ordered_means=mean_normalized[order_ix];
cat("\n");
cat("Top Categories:\n");
head(ordered_means);
cat("\n");
	 
###############################################################################
# Generate prop/model files

samp_ids=rownames(normalized);

coverage_cutoffs=c(0.75, 0.80, 0.85, 0.90, 0.95);
abd_cutoffs=c(0.01, 0.02, 0.025, 0.03, 0.04, 0.05);

num_cov_cutoffs=length(coverage_cutoffs);
num_abd_cutoffs=length(abd_cutoffs);

num_cov_matrix=matrix(NA, nrow=num_samples, ncol=num_cov_cutoffs,
	dimnames=list(samp_ids, coverage_cutoffs));

num_abd_matrix=matrix(NA, nrow=num_samples, ncol=num_abd_cutoffs,
	dimnames=list(samp_ids, abd_cutoffs));

cumul_names_by_cov=list();
cumul_names_by_abd=list();
for(cix in 1:num_cov_cutoffs){
	cumul_names_by_cov[[cix]]="";
}
names(cumul_names_by_cov)=coverage_cutoffs;

for(cix in 1:num_abd_cutoffs){
	cumul_names_by_abd[[cix]]="";
}
names(cumul_names_by_abd)=abd_cutoffs;

for(smp_nm in samp_ids){

	cat("\n\n\n");
	cat("Working on: ", smp_nm, "\n");

	# Sorting sample profile
	samp_prof=normalized[smp_nm,,drop=F];
	names(samp_prof)=colnames(normalized);
	order_ix=order(samp_prof, decreasing=T);
	samp_prof=samp_prof[order_ix];
	taxa_names=names(samp_prof);
	cumul=cumsum(samp_prof);

	# Count taxa need to achieve coverage
	num_taxa_cov=numeric(num_cov_cutoffs);
	for(cix in 1:num_cov_cutoffs){
		num_taxa_cov[cix]=sum(cumul<=coverage_cutoffs[cix]);
		cumul_names_by_cov[[cix]]=
				unique(c(cumul_names_by_cov[[cix]], taxa_names[1:num_taxa_cov[cix]]));
	}

	cat("\n");
	cat("Num Taxa by coverage:\n");
	names(num_taxa_cov)=coverage_cutoffs;	
	print(num_taxa_cov);

	# Count taxa by abundance
	#print(samp_prof);	
	num_taxa_abv_cutoff=numeric(num_abd_cutoffs);
	for(cix in 1:num_abd_cutoffs){
		num_taxa_abv_cutoff[cix]=sum(samp_prof>=abd_cutoffs[cix]);
		cumul_names_by_abd[[cix]]=
				unique(c(cumul_names_by_abd[[cix]], taxa_names[1:num_taxa_abv_cutoff[cix]]));
	}

	cat("\n");
	cat("Num Taxa above abundance:\n");
	names(num_taxa_abv_cutoff)=abd_cutoffs;
	print(num_taxa_abv_cutoff);

	cat("\n");

	num_cov_matrix[smp_nm,]=num_taxa_cov;
	num_abd_matrix[smp_nm,]=num_taxa_abv_cutoff;

	cat("\n\n\n");
}

cumul_cov_counts=(unlist(lapply(cumul_names_by_cov, length)));
cumul_abd_counts=(unlist(lapply(cumul_names_by_abd, length)));

#------------------------------------------------------------------------------

plot_cutoffs=function(mat, cmct, title){
	x_val=as.numeric(colnames(mat));

	x_max=max(x_val);
	x_min=min(x_val);
	x_range=x_max-x_min;

	y_max=max(mat);
	
	num_rows=nrow(mat);

	plot(NA, type="n", xlim=c(x_min-x_range/10, x_max+x_range/10), ylim=c(0, y_max),
		xlab="Cutoff", ylab="Num Taxa Acquired",
		main=title
		);

	title("[Cumulative Distinct Taxa Count]", col.main="blue", font.main=2, line=1, cex.main=.5);
	title("[Median Sample Taxa Count]", col.main="red", font.main=2, line=0.5, cex.main=.5);

	x_spacing=min(diff(x_val));

	num_col=ncol(mat);
	for(i in 1:num_col){
		jitter=rnorm(num_rows, 0, x_spacing/10);
		points(rep(x_val[i], num_rows)+jitter, mat[,i], col="grey");
		med=median(mat[,i]);
		points(x_val[i], med, pch="-", cex=3, col="red");
		text(x_val[i], x_max, med, font=2, col="red");
		text(x_val[i], y_max, cmct[i], font=2, col="blue");
	}
}

plot_cutoffs(num_cov_matrix, cumul_cov_counts, "Cumulative Abundance of Admitted");
plot_cutoffs(num_abd_matrix, cumul_abd_counts, "Minimum Abundance for Admittance");

#------------------------------------------------------------------------------

export_targeted=function(list_by_cutoff, output_root, name_to_xml_map){
	
	num_cutoffs=length(list_by_cutoff);
	cutoffs_arr=names(list_by_cutoff);

	for(i in 1:num_cutoffs){

		cutoff=cutoffs_arr[i];

		cat("\n\n");
		cat("Exporting: ", cutoff, "\n");

		target_list=list_by_cutoff[[i]];
		target_list=setdiff(target_list, "");

		xmls=unlist(name_to_xml_map)[target_list];
		mat=cbind(target_list, xmls);

		fname=paste(output_root, ".", cutoff, ".map", sep="");
		cat("Filename: ", fname, "\n");

		print(mat);
		colnames(mat)=c("Name", "Model");
		write.table(mat, file=fname, quote=F, sep="\t", row.names=F, col.names=T);
	}

}

#print(model_mapping);
export_targeted(cumul_names_by_cov, paste(OutputRoot, ".mincov", sep=""), model_mapping);
export_targeted(cumul_names_by_abd, paste(OutputRoot, ".minabd", sep=""), model_mapping);

###############################################################################

cat("Done.\n")
print(warnings());

q(status=0)
