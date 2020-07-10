#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.tsv>\n",
	"	[-o <output summary table root name>]\n",
	"\n",	
	"This script will read in a summary file and \n",
	"clean up the categories names so they are R friendly\n",
	"for variable names.\n",
	"\n",
	"These will be converted:\n",
	"	- to _\n",
	"	[ to ''\n",
	"	] to ''\n",
	"	( to ''\n",
	"	) to ''\n",
	"\n",
	"If there is an 'Incertae_Sedis', replace with <taxa above>_IncSed\n",
	"\n",
	"For sake of brevity, these components will also be abbreviated:\n",
	"	Unknown_	to Unk\n",
	"	_unclassified	to _uncl\n",
	"	_group		to _grp\n",
	"\n",
	"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$output_file)){
	outputroot=gsub("\\.summary_table\\.xls", "", opt$input_file);
	outputroot=gsub("\\.summary_table\\.tsv", "", opt$input_file);
	OutputFileName = paste(outputroot, ".tx_cln.summary_table.tsv", sep="");
}else{
	outputroot=gsub("\\.summary_table\\.xls", "", opt$output_file);
	outputroot=gsub("\\.summary_table\\.tsv", "", opt$output_file);
	OutputFileName=opt$output_file;
}

###############################################################################

InputFileName=opt$input_file;

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFileName, "\n");       
cat("\n");

###############################################################################
###############################################################################

load_summary_table=function(summary_table_fn){
        # Load data
        cat("Loading Matrix (", summary_table_fn, ") ...\n", sep="");
        inmat=as.matrix(read.table(summary_table_fn, sep="\t", header=TRUE, check.names=FALSE, row.names=1, quote=""))

        #cat("\nOriginal Matrix:\n")
        #print(inmat);

        # Grab columns we need into a vector, ignore totals, we won't trust it.
        counts_mat=inmat[,2:(ncol(inmat))];
        #cat("\nCounts Matrix:\n");
        #print(counts_mat);

        num_samples=nrow(counts_mat);
        num_categories=ncol(counts_mat);
        sample_names=rownames(counts_mat);

        cat("\n");
        cat("Num Samples: ", num_samples, "\n");
        cat("Num Categories: ", num_categories, "\n");
        cat("\n");
        return(counts_mat);
}

###############################################################################

# Load counts matrix
counts_mat=load_summary_table(InputFileName);
num_categories=ncol(counts_mat);
num_samples=nrow(counts_mat);

# Get the category names
category_names=colnames(counts_mat);

adjusted_cat_names=character(num_categories);

num_cat_adj=0;
for(i in 1:num_categories){
	cur_cat=category_names[i];
	
	res=cur_cat;

	res=gsub("-", "_", res); 
	res=gsub("\\[", "", res); 
	res=gsub("\\]", "", res); 
	res=gsub("\\(", "", res); 
	res=gsub("\\)", "", res); 

	res=gsub(",", "_", res); 
	res=gsub("'", "p", res); 
	res=gsub("\\+", "p", res); 
	res=gsub("/", "_", res); 
	res=gsub(":", "_", res); 
	res=gsub(" ", "_", res); 
	res=gsub("_+", "_", res);

	if(length(grep("^[0-9]", res))){
		res=paste("n",res, sep="");
	}

	res=gsub("_unclassified", "_uncl", res); 
	res=gsub("_group", "_grp", res); 

	inssed_found=length(grep("Incertae_Sedis", res))>0;
	unknown_found=length(grep("Unknown_", res))>0;
	uncultured_found=length(grep("uncultured", res))>0;

	if(inssed_found || unknown_found || uncultured_found){ 

		splits=strsplit(res, ";")[[1]];
		num_splits=length(splits);
		for(s in 1:num_splits){
			
			tok=splits[s];

			if(tok=="Unknown_Phylum"){
				tok=paste(splits[s-1], "_UnkPhylum", sep="");
			}else if(tok=="Unknown_Class"){
				tok=paste(splits[s-1], "_UnkClass", sep="");
			}else if(tok=="Unknown_Order"){
				tok=paste(splits[s-1], "_UnkOrder", sep="");
			}else if(tok=="Unknown_Family"){
				tok=paste(splits[s-1], "_UnkFamily", sep="");
			}else if(tok=="Unknown_Genus"){
				tok=paste(splits[s-1], "_UnkGenus", sep="");

			}else if(tok=="Phylum_Incertae_Sedis"){
				tok=paste(splits[s-1], "_PhylumIncSed", sep="");
			}else if(tok=="Class_Incertae_Sedis"){
				tok=paste(splits[s-1], "_ClassIncSed", sep="");
			}else if(tok=="Order_Incertae_Sedis"){
				tok=paste(splits[s-1], "_OrderIncSed", sep="");
			}else if(tok=="Family_Incertae_Sedis"){
				tok=paste(splits[s-1], "_FamilyIncSed", sep="");
			}else if(tok=="Genus_Incertae_Sedis"){
				tok=paste(splits[s-1], "_GenusIncSed", sep="");

			}else if(tok=="uncultured"){
				tok=paste(splits[s-1], "_Uncltrd", sep="");

			}else if(tok=="Incertae_Sedis"){
				tok=paste(splits[s-1], "_IncSed", sep="");
			}

			splits[s]=tok;
		}

		res=paste(splits, collapse=";");
	}

	if(res!=cur_cat){
		cat(cur_cat, " ->\n");
		if(inssed_found){
			cat("\t** Incertae Sedis Adjusted **\n");
		}
		cat("\t", res, "\n\n", sep="");	
		num_cat_adj=num_cat_adj+1;
	}

	adjusted_cat_names[i]=res;
}

cat("Number of Categories Adjusted: ", num_cat_adj, "\n", sep="");

# Rename the counts_mat
colnames(counts_mat)=adjusted_cat_names;
outmat=counts_mat;

###############################################################################
# Output
cat("\nWriting New Matrix...\n");
fc=file(paste(OutputFileName), "w");

write(paste("sample_id", "total", paste(colnames(outmat), collapse="\t"), sep="\t"), file=fc);
sample_names=rownames(counts_mat);
for(samp_idx in 1:num_samples){
	total=sum(outmat[samp_idx,]);
	outline=paste(sample_names[samp_idx],total,paste(outmat[samp_idx,], collapse="\t"), sep="\t");
	write(outline, file=fc);
}
close(fc);	

###############################################################################

cat("Done.\n")
warns=warnings();
if(!is.null(warns)){
	print(warnings());
}
q(status=0)
