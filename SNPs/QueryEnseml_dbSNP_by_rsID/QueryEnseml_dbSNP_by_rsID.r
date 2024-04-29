#!/usr/bin/env Rscript

###############################################################################

library(getopt);

# Need to install biomaRt:
# dzdo yum install libxml2-devel
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("biomaRt")

options(useFancyQuotes=F);
options(width=120);

params=c(
	"input_filename", "i", 1, "character",
	"target_colname", "t", 1, "character",
	"output_filename", "o", 1, "character",
	"mart_dataset", "d", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

DEF_DATASET="hsapiens_snp";

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input file with multiple tables>\n",
	"	-t <column name of rs SNP IDs to look up>\n",
	"	-o <output file with chrom and coordinates appended to input>\n",
	"	[-d <dataset, default=hsapiens_snp>]\n",
	"\n",
	"This script will read in a list of rs IDs\n",
	"and then extract the chromosome and coordinates\n",
	"of the SNP from the database.\n",
	"\n",
	"The output will consist of the downloaded data appended to the input table.\n", 
	"\n",
	"The input can have multiple columns, in case you want to presever some comments.\n",
	"\n",
	"The dataset database may also be mmusculus_snp, rnorvegicus_snp, and others.\n",
	"\n");

if(
	!length(opt$input_filename) || 
	!length(opt$target_colname) ||
	!length(opt$output_filename)
){
	cat(usage);
	q(status=-1);
}

Dataset=DEF_DATASET;
if(length(opt$mart_dataset)){
	Dataset=opt$mart_dataset;
}

InputFile=opt$input_filename;
OutputFile=opt$output_filename;
TargetColname=opt$target_colname;

cat("\n");
cat("Input File: ", InputFile, "\n");
cat("Target Column (w/ rs IDs): ", TargetColname, "\n");
cat("Output File: ", OutputFile, "\n");
cat("Dataset: ", Dataset, "\n");
cat("\n");

##############################################################################
# Loading targeted SNPs

cat("Loading input file: ", InputFile, "\n");
intab=read.delim(InputFile, header=T, check.names=F, row.names=NULL,
	as.is=T, comment.char="#", quote="", sep="\t");

incolnames=colnames(intab);
num_targets=nrow(intab);
cat("Num Rows Read: ", num_targets, "\n");
cat("\n");

cat("Column names in input file:\n");
print(incolnames);
cat("\n");

# These are the IDs that will be queried
rs_ids=intab[,TargetColname];
cat("IDs found in ", TargetColname, " Column\n", sep="");
print(rs_ids);
cat("\n\n");

##############################################################################
# Querying SNPs

query_attrib_columns=c('refsnp_id', 'chr_name', 'chrom_start','chrom_end', 'allele');

cat("Starting up biomaRt library... [", date(), "]\n");
cat("  (may take 10s)\n");
library(biomaRt);
cat("ok. [", date(), "]\n\n");

cat("Starting up SNP mart with dataset = ", Dataset, " [", date(), "]\n");
cat("  (may take 10s)\n");
snp_mart=useEnsembl(biomart="snp", dataset=Dataset);
cat("ok. [", date(), "]\n\n");

dataset_list=listDatasets(snp_mart);
dataset_info=dataset_list[
	dataset_list[,"dataset"]==Dataset, c("dataset", "description", "version")];

cat("Running query... [", date(), "]\n");
res_tab=getBM(attributes = query_attrib_columns,
      	filters = c('snp_filter'), 
	values=rs_ids, mart=snp_mart);
cat("ok. [", date(), "]\n\n");

cat("\n\n");
cat("Returned from query:\n");
print(res_tab);

##############################################################################
# Cleaning up table

num_res_rows=nrow(res_tab);
keep=rep(T, num_res_rows);
for(i in 1:num_res_rows){
	chr_name=as.character(res_tab[i,"chr_name"]);
	if(length(grep("_CTG\\d+$", chr_name))){
		keep[i]=F;
	}else if(length(grep("_PATCH$", chr_name))){
		keep[i]=F;
	}
}

clean_res_tab=res_tab[keep,,drop=F];
num_clean=nrow(clean_res_tab);

res_tab_refsnp_ids=clean_res_tab[,"refsnp_id"];
tryCatch({
	rownames(clean_res_tab)=res_tab_refsnp_ids;
	}, error=function(msg){
		# This should catch duplicated or insufficiently cleaned up results
		cat("Error: \n");
		print(msg);
		cat("\n");
		cat("Cleaned Result Table:\n");
		print(clean_res_tab);
		quit(status=-1);
	});
clean_res_tab=clean_res_tab[rs_ids,,drop=F];

cat("\n");
cat("Cleaned up and reordered table:\n");
print(clean_res_tab);

cat("\n");
# This will catch missing SNPs
cat("Number of Targets: ", num_targets, "\n");
cat("Number of Cleaned Up IDs: ", num_clean, "\n");
if(num_clean!=num_targets){
	cat("******************************************************************\n");
	cat("Error: Number of Targets don't match returned in cleaned up table.\n");
	cat("Please check input file or nature of returned values from the query.\n");
	cat("******************************************************************\n");
	cat("\n");
	cat("Missing SNP IDs: \n");
	print(setdiff(rs_ids, res_tab_refsnp_ids));
	cat("******************************************************************\n");
	quit(status=-1);
}

##############################################################################
# Split up ref/variant

ref_allele=character(num_clean);
var_allele=character(num_clean);

for(i in 1:num_clean){
	alleles=clean_res_tab[i,"allele"];
	allele_arr=strsplit(alleles, "/")[[1]];
	ref_allele[i]=allele_arr[1];
	var_allele[i]=paste(allele_arr[2:length(allele_arr)], collapse="/");
}

clean_res_tab=cbind(clean_res_tab, ref_allele, var_allele);
print(clean_res_tab);

##############################################################################
# Append to input file

cat("\n");
cat("Appending queried information to input table...\n");
#print(intab[,TargetColname]);
rownames(intab)=intab[,TargetColname];
outtab=cbind(
	clean_res_tab[rs_ids, 
		c("refsnp_id", "chr_name", 
		"chrom_start", "chrom_end", "ref_allele", "var_allele"), drop=F],
	intab[rs_ids,,drop=F]
)

#print(outtab);
cat("Writing output table...\n");
write.table(x=outtab, file=OutputFile, row.names=F, col.names=T, sep="\t", quote=F);
cat("Writing version as comment to end of table....\n");
fh=file(OutputFile, "a");
cat(file=fh, "# ", paste(dataset_info, collapse="/"), "\n", sep="");
close(fh);


cat("\nDone.\n");

##############################################################################

print(warnings());
q(status=0);
