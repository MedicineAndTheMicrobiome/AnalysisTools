#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

source('~/git/AnalysisTools/Metadata/InputFileLibrary/InputFileLibrary.r');
source('~/git/AnalysisTools/Metadata/OutputFileLibrary/OutputFileLibrary.r');

options(useFancyQuotes=F);

params=c(
	"summary_file", "s", 1, "character",
	"num_top_cat", "p", 2, "character",
	"shorten_category_names", "x", 2, "character",
	"pairings_map_fn", "P", 1, "character",
	"samp_a_cname", "A", 1, "character",
	"samp_b_cname", "B", 1, "character",
	"subject_id_cname", "S", 1, "character",
	"empir_count_fn", "n", 1, "character",
	"ec_sampid_cname", "i", 1, "character",
	"ec_counts_cname", "j", 1, "character",
	"output_root_fn", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

NUM_TOP_CAT=35;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file .tsv table>\n",
	"	[-p <num top taxa to include, default=", NUM_TOP_CAT, " >]\n",
	"	[-x <shorten category names, with separator in double quotes>]\n",
	"\n",
	"	-P <pairings map, pairing Resp and Pred sample IDs. Must have header/column names>\n",
	"	-A <Sample Type A (Pos Fraction) Colname>\n",
	"	-B <Sample Type B (Neg Fraction) Colname>\n",
	"	-S <Subject ID Colname>\n",
	"\n", 
	"	-n <empirical (qPCR) count, (sample_id , counts)>\n",
	"	-i <Sample ID Colname>\n",
	"	-j <Sample Counts Colname>\n",
	"\n",
	"	[-o <output tsv name>]\n",
	"\n",
	"This script will compute the composition for each sample, then based on the qPCR\n",
	"counts, will estimate the probability a taxon will be found in Sample A or \n",
	"Sample B.  The assumption is that a sample was passed through a filter/sorter.\n",
	"The filter has an effect, but it's not reliable, so are trying to find the probability\n",
	"of a being in A or B.  The summary file will tell us the composition in A or B,\n",
	"but we need to know the qPCR counts (the actual relative recovered items in A vs. B)\n",
	"\n", sep="");

if(
	!length(opt$summary_file) || 
	!length(opt$shorten_category_names) || 
	!length(opt$pairings_map_fn) || 
	!length(opt$samp_a_cname) || 
	!length(opt$samp_b_cname) || 
	!length(opt$subject_id_cname) || 
	!length(opt$empir_count_fn) ||
	!length(opt$ec_sampid_cname) ||
	!length(opt$ec_counts_cname)
){
	cat(usage);
	q(status=-1);
}


SummaryFile=opt$summary_file;
NumTopCat=opt$num_top_cat;
PairingsFile=opt$pairings;
SampAColname=opt$samp_a_cname;
SampBColname=opt$samp_b_cname;
SubjIDColname=opt$subject_id_cname;
EmpirCountFile=opt$empir_count_fn;
EmpirCountSampIDColname=opt$ec_sampid_cname;
EmpirCountCountsColname=opt$ec_counts_cname;


if(!length(opt$output_root_fn)){
	OutputRoot=gsub(".summary_table.xls", "", opt$summary_file);
	OutputRoot=gsub(".summary_table.tsv", "", OutputRoot);
}else{
	OutputRoot=opt$output_root_fn;
}

if(!length(opt$num_resp_var)){
	NumTopALR=NUM_TOP_CAT;
}else{
	NumTopALR=opt$num_top_cat;
}

if(!length(opt$shorten_category_names)){
	ShortenCategoryNames="";
}else{
	ShortenCategoryNames=opt$shorten_category_names;
}


cat("\n");
cat("                      Summary File: ", SummaryFile, "\n", sep="");
cat("                Num Top Categories: ", NumTopCat, "\n", sep="");
cat("\n");
cat("                     Pairings File: ", PairingsFile, "\n", sep="");
cat("                    Samp A Colname: ", SampAColname, "\n", sep="");
cat("                    Samp B Colname: ", SampBColname, "\n", sep="");
cat("                Subject ID Colname: ", SubjIDColname, "\n", sep="");
cat("\n");
cat("             Empirical Counts File: ", EmpirCountFile, "\n", sep="");
cat("Empirical Counts Sample ID Colname: ", EmpirCountSampIDColname, "\n", sep="");
cat("    Empirical Counts Couts Colname: ", EmpirCountCountsColname, "\n", sep="");
cat("\n");
cat("                       Output File: ", OutputRoot, "\n", sep="");
cat("            Shorten Category Names: '", ShortenCategoryNames, "'\n", sep="");
cat("\n");

if(ShortenCategoryNames==TRUE){
	cat("Error:  You need to specify a delimitor to split the category names.\n");
	cat("        i.e., this -x option is not a flag, it requires a parameter.\n");
	quit(status=-1);
}

options(width=100);
cat("Text Line Width: ", options()$width, "\n", sep="");

##############################################################################
##############################################################################

# Count Adjusted Logit

pdf(paste(OutputRoot, ".cal.pdf", sep=""), height=11, width=9.5);

input_files=load_and_reconcile_files(
                sumtab=list(fn=SummaryFile, 
			shorten_cat_names_char=ShortenCategoryNames,
                        return_top=NumTopCat,
			specific_cat_fn=""
			),
                factors=list(fn=""),
                pairs=list(fn=PairingsFile, a_cname=SampAColname, b_cname=SampBColname),
                covariates=list(fn=""),
                grpvar=list(fn=""),
                reqvar=list(fn="")
        );

#print(names(input_files));
normalized=input_files[["SummaryTable_normalized"]];
pairs_map=input_files[["PairsMap"]];
sample_ids=rownames(normalized);

print(input_files[["Report"]]);
write_file_report(input_files[["Report"]]);

###############################################################################
# Load empirical counts

emp_cnts_in=read.table(EmpirCountFile, header=T, row.names=EmpirCountSampIDColname, sep="\t");
emp_cnts_wna=emp_cnts_in[,EmpirCountCountsColname,drop=F];
emp_cnts=emp_cnts_wna[!is.na(emp_cnts_wna),,drop=F];
emp_samp_ids=rownames(emp_cnts);

cat("Example Counts:\n");
print(head(emp_cnts));
cat("...\n\n");

###############################################################################

num_subject_ids=nrow(pairs_map);
subject_ids=rownames(pairs_map);
cat("Num Pairs/Subjects in Map: ", num_subject_ids, "\n");

categories=colnames(normalized);
num_cat=length(categories);
cat("Number of Categories: ", num_cat, "\n", sep="");
cat("Categories:\n");
print(categories);

###############################################################################
# Apply emp counts to each sample ID

count_adj_mat=normalized;
epsi=.5;
for(smp_ix in sample_ids){
	cat("Sample ID: ", smp_ix, "\n");
	cnt=emp_cnts[smp_ix,];
	composition=normalized[smp_ix,];
	print(cnt);
	print(composition);
	count_adj_mat[smp_ix,]=(cnt*composition)+epsi;
}

###############################################################################
# Compute ratios per subject id

cat("\n");
cat("Calculating Logits:\n");

logit_mat=matrix(0, nrow=num_subject_ids, ncol=ncol(normalized));
colnames(logit_mat)=colnames(normalized);
rownames(logit_mat)=subject_ids;

print(logit_mat);

for(sbj_ix in subject_ids){
	
	cat("Subject ID: ", sbj_ix, "\n");
	a_smp_id=pairs_map[sbj_ix,SampAColname];
	b_smp_id=pairs_map[sbj_ix,SampBColname];

	sampA_cnts=count_adj_mat[a_smp_id,];
	sampB_cnts=count_adj_mat[b_smp_id,];

	sums=sampA_cnts+sampB_cnts;

	sampA_norm=sampA_cnts/sums;
	sampB_norm=sampB_cnts/sums;

	logit_mat[sbj_ix,]=log(sampA_norm/sampB_norm);
}

print(logit_mat);

###############################################################################

par(mar=c(4,4,4,1));
par(mfrow=c(4,3));
for(cat_ix in categories){
	cat("Generating Histogram for: ", cat_ix, "\n");
	logit=logit_mat[,cat_ix];
	logit=logit[!is.na(logit)];
	mav=max(abs(range(logit)));
	hist(logit, main=cat_ix, breaks=seq(-mav, mav, length.out=10));
	axis(side=1, at=c(mav, -mav), labels=c(SampAColname, SampBColname), 
		line=1, cex=.5, tick=F);
	abline(v=0, col="blue", lty="dashed");

}


###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
