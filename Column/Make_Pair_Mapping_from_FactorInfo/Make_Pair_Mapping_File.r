#!/usr/bin/env Rscript

###############################################################################

library('getopt');

options(useFancyQuotes=F);

params=c(
	"factors", "f", 1, "character",
	"subject_id_colname", "s", 1, "character",
	"sample_type_colname","t", 1, "character",
	"output", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-f <factors file, contains covariates and factors>\n",
	"	-s <subject id column name>\n",
	"	-t <sample type column name>\n",
	"	[-o <output map file name>]\n",
	"\n",
	"This script will read in a factor file, then based on the subject\n",
	"id, look for the corresponding sample with a different sample type.\n",
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

if(length(opt$output)){
	OutputRoot=opt$output;
}else{
	OutputRoot=gsub(".tsv$", "", FactorFilename);
}

cat("           Factors File: ", FactorFilename, "\n", sep="");
cat(" Subject ID Column Name: ", SubjectIDColname, "\n", sep="");
cat("Sample Type Column Name: ", SampleTypeColname, "\n", sep="");
cat("            Output File: ", OutputRoot, "\n", sep="");
cat("\n");

options(width=100);
cat("Text Line Width: ", options()$width, "\n", sep="");

##############################################################################
##############################################################################

load_factors=function(fname){
	factors=data.frame(read.table(fname,  sep="\t", header=TRUE, 
		row.names=1, check.names=FALSE, comment.char=""));
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
factors=load_factors(FactorFilename);
factor_names=colnames(factors);
num_factors=ncol(factors);
factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);

###############################################################################

# Extract info we need for mapping

cat("Num Samples to pair: ", nrow(factors), "\n");

subj_ids=factors[,SubjectIDColname];
uniq_subj_ids=unique(subj_ids);
num_subj_ids=length(uniq_subj_ids);

sample_types=factors[,SampleTypeColname];
uniq_sample_types=sort(unique(sample_types));
num_samp_types=length(uniq_sample_types);

cat("\nNum Sample Types (", SampleTypeColname, "):", num_samp_types, "\n");
print(head(uniq_sample_types));

cat("\nNum Subject IDs (", SubjectIDColname, "): ",  num_subj_ids, "\n");
print(head(uniq_subj_ids)); 

if(num_samp_types!=2){
	cat("Cannot build mapping file when there are not exactly 2 sample types.\n");	
	quit(-1);
}

sample_ids=rownames(factors);

# Split 

type1_samples_ix=(sample_types==uniq_sample_types[1]);
type2_samples_ix=(sample_types==uniq_sample_types[2]);

num_type1=sum(type1_samples_ix);
num_type2=sum(type2_samples_ix);

cat("\n");
cat("Num samples from ", SampleTypeColname, "(", uniq_sample_types[1], ")", num_type1, "\n");
cat("Num samples from ", SampleTypeColname, "(", uniq_sample_types[2], ")", num_type2, "\n");
cat("\n");

###############################################################################

unmatched_sample_ids=character();

sample_map=c();

for(si in uniq_subj_ids){
	#cat("\n");
	#cat("Working on: ", si, "\n");
	cur_ic=(si==subj_ids);

	#print(sample_ids[cur_ic]);
	#print(sample_types[cur_ic]);

	num_samps_from_subj=sum(cur_ic);
	if(num_samps_from_subj==1){
		cat("Unmatched: ", si, "\n");
		unmatched_sample_ids=c(unmatched_sample_ids, sample_ids[cur_ic]);
	}else if(num_samps_from_subj>2){
		cat("Error:  Ambiguous mapping for: ", si, "\n");
		print(sample_ids[cur_ic]);	
		quit(-1);
	}else{
		cur_samp_ids=sample_ids[cur_ic];
		cur_samp_types=sample_types[cur_ic];
		s1=which(cur_samp_types[1]==uniq_sample_types);
		s2=which(cur_samp_types[2]==uniq_sample_types);
		sample_map=rbind(sample_map, c(cur_samp_ids[s1], cur_samp_ids[s2]));
	}

}

cat("Unmatched Samples (", length(unmatched_sample_ids), "):\n");
print(unmatched_sample_ids);

cat("\n");
colnames(sample_map)=uniq_sample_types;
cat("Matched Samples (", nrow(sample_map), "):\n");
print(sample_map);

###############################################################################

#fh=open(OutputRoot, "w");
write.table(sample_map, paste(OutputRoot, ".", SampleTypeColname, ".map.tsv", sep=""), quote=F, sep="\t", row.names=F);

#print(mapping_info);

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
