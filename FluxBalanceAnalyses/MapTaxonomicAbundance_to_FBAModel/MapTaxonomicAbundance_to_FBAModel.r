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
	"10, 15, 20, 25, and 30 taxa that can be mapped to an available\n",
	"FBA model.\n",
	"\n",
	"The 'Override FBA model list' is the same format as the -x option\n",
	"The purpose is to remap some taxa to specific xml files in case\n",
	"the taxa (genus) doesn't have a mapping in the model list.\n",
	"\n",
	"The output contains a row for each sample, and a column for each\n",
	"model that was selected.\n",
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

OutputRoot=opt$output_root;

cat("\n");
cat("Input Summary Table: ", InputSummaryTable, "\n");
cat("FBA Model List: ", FBAModelList, "\n");
cat("Override List File: ", OverrideFBAModelList, "\n");
cat("Output Root: ", OutputRoot, "\n");
cat("\n");

pdf(paste(OutputRoot, ".fba_model_select.pdf", sep=""), height=11, width=8.5);

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

clean_counts=remove_categories(c("Homo"), orig_counts);
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

# Calculating pres/abs
pres_abs=apply(clean_counts, 2, function(s){sum(s>0)});
pres_abs=pres_abs[order_ix]

#------------------------------------------------------------------------------
# Generate summary stats and join with available XML files

cumulative_matrix=cbind(1:num_categories, ordered_means, cumsum(ordered_means), pres_abs);
colnames(cumulative_matrix)=c("Index", "Abundance", "Cumulative", "NumPresent");
#print(cumulative_matrix);

Model=lookup(rownames(cumulative_matrix), model_mapping);
cumulative_matrix=cbind(cumulative_matrix, Model);

# Report Average Statistics
outfn=paste(OutputRoot, ".joined.info.table.tsv", sep="");
outmat=cbind(rownames(cumulative_matrix), cumulative_matrix);
colnames(outmat)=c("Taxa", colnames(cumulative_matrix));
print(outmat[1:50,]);
write.table(
	x=outmat,
	file=outfn,
	quote=F, sep="\t", row.names=F, col.names=T);
	 
###############################################################################
# Generate prop/model files

topn=c(10, 15, 20, 25, 30, 35);

# Remove genera without models
found_models=cumulative_matrix[!is.na(cumulative_matrix[,"Model"]),];
print(found_models);

for(cutoff in topn){

	top_models=found_models[1:cutoff, "Model"];
	print(top_models);

	kept=normalized[,names(top_models),drop=F];
	remaining=apply(kept, 1, function(x){ 1-sum(x)});

	outmat=cbind(rownames(kept), kept, remaining);
	colnames(outmat)=c("SampleID", top_models, "Remaining");

	outfn=paste(OutputRoot, ".top", cutoff, ".unnorm_model_prop_tab.tsv", sep="");
	write.table(
		x=outmat,
		file=outfn,
		quote=F, sep="\t", row.names=F, col.names=T);

}

###############################################################################

cat("Done.\n")
print(warnings());

q(status=0)
