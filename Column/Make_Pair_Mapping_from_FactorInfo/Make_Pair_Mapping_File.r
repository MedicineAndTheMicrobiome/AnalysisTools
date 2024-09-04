#!/usr/bin/env Rscript

###############################################################################

library('getopt');

options(useFancyQuotes=F);

params=c(
	"factors", "f", 1, "character",
	"subject_id_colname", "s", 1, "character",
	"sample_type_colname","t", 1, "character",
	"sample_id_colname", "r", 1, "character",
	"output", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <factors file, contains covariates and factors>\n",
	"	-s <subject id column name>\n",
	"	-t <sample type column name>\n",
	"	[-r <sample id column name (row name)>]\n",
	"	[-o <output map file name>]\n",
	"\n",
	"This script will read in a factor file, then based on the subject\n",
	"id, look for the corresponding sample with a different sample type.\n",
	"\n",
	"A 2D matrix is created with the subject ID as rows, and the\n",	
	"sample type as columns.  If there are any missing pairings\n",
	"the sample ID missing will have an NA.\n",
	"\n", sep="");

if(
	!length(opt$factors) || 
	!length(opt$subject_id_colname) || 
	!length(opt$sample_type_colname)
){
	cat(usage);
	q(status=-1);
}

FactorFilename=opt$factors;
SubjectIDColname=opt$subject_id_colname;
SampleTypeColname=opt$sample_type_colname;
SampleIDColname=opt$sample_id_colname;

if(length(opt$output)){
	OutputRoot=opt$output;
}else{
	OutputRoot=gsub(".tsv$", "", FactorFilename);
}

cat("           Factors File: ", FactorFilename, "\n", sep="");
cat("  Sample ID Column Name: ", SampleIDColname, "\n", sep="");
cat(" Subject ID Column Name: ", SubjectIDColname, "\n", sep="");
cat("Sample Type Column Name: ", SampleTypeColname, "\n", sep="");
cat("            Output File: ", OutputRoot, "\n", sep="");
cat("\n");

options(width=100);
cat("Text Line Width: ", options()$width, "\n", sep="");

##############################################################################
##############################################################################

load_factors=function(fname, rowname){
	factors=data.frame(read.table(fname,  sep="\t", header=TRUE, 
		row.names=rowname, as.is=T, check.names=FALSE, comment.char=""));
	factor_names=colnames(factors);

	ignore_idx=grep("^IGNORE\\.", factor_names);

	if(length(ignore_idx)!=0){
		return(factors[-ignore_idx]);
	}else{
		return(factors);
	}
}

##############################################################################

# Load factors
cat("Loading Factors...\n");
factors=load_factors(FactorFilename, SampleIDColname);
factor_names=colnames(factors);
num_factors=ncol(factors);
factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);

###############################################################################

# Extract info we need for mapping

cat("Num Samples to pair: ", nrow(factors), "\n");

if(length(intersect(SubjectIDColname, factor_names))==0){
	cat("\nCould not find Subject ID column: ", SubjectIDColname, "\n");
	cat("Available Columns: \n");
	print(factor_names);	
	quit();
}

subj_ids=factors[,SubjectIDColname];
nona_ix=!is.na(subj_ids);

# Remove rows without subject ids
factors=factors[nona_ix,,drop=F];
subj_ids=factors[,SubjectIDColname];

uniq_subj_ids=unique(subj_ids);
num_subj_ids=length(uniq_subj_ids);

sample_types=as.character(factors[,SampleTypeColname]);
uniq_sample_types=sort(unique(sample_types));
num_samp_types=length(uniq_sample_types);

# Make sure sample types are not numeric
if(all(!is.na(as.numeric(as.character(uniq_sample_types))))){
	initial_char=substr(SampleTypeColname,1,1);
	cat("Numbers found as sample type.  Appending ", initial_char, " in front.\n");
	asn_sample_types=as.numeric(as.character(sample_types));
	fmt_str=paste("%0", log10(asn_sample_types+1)+1, "i", sep="");
	sample_types=paste(initial_char, sprintf(fmt_str, asn_sample_types), sep="");
	uniq_sample_types=sort(unique(sample_types));
	num_samp_types=length(uniq_sample_types);
}

# Make sure sample IDs are not numeric
if(all(!is.na(as.numeric(as.character(uniq_subj_ids))))){
	initial_char=substr(SubjectIDColname,1,1);
	cat("Numbers found as Subject ID.  Appending ", initial_char, " in front.\n");
	fmt_str=paste("%0", log10(as.numeric(as.character(subj_ids))+1)+1, "i", sep="");
	subj_ids=paste(initial_char, sprintf(fmt_str, subj_ids), sep="");
	uniq_subj_ids=sort(unique(subj_ids));
	num_subj_ids=length(uniq_subj_ids);
}


cat("\nNum Sample Types (", SampleTypeColname, "):", num_samp_types, "\n");
print(head(uniq_sample_types));
cat("\n");

cat("\nNum Subject IDs (", SubjectIDColname, "): ",  num_subj_ids, "\n");
print(head(uniq_subj_ids)); 
cat("\n");

mapping=matrix(NA, nrow=num_subj_ids, ncol=num_samp_types);
rownames(mapping)=uniq_subj_ids;
colnames(mapping)=uniq_sample_types;

samp_ids=rownames(factors);
for(i in 1:num_factor_samples){
	if(!is.na(sample_types[i]) && !is.na(subj_ids[i])){
		mapping[subj_ids[i], sample_types[i]]=samp_ids[i];

	}
}

###############################################################################

mapping_column_names=c(SubjectIDColname, colnames(mapping));
mapping=cbind(as.character(uniq_subj_ids), mapping);
colnames(mapping)=mapping_column_names;

print(mapping);
write.table(mapping, paste(OutputRoot, ".", SampleTypeColname, ".map.tsv", sep=""), quote=F, sep="\t", row.names=F);

###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
