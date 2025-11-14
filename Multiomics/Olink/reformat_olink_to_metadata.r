#!/usr/bin/env Rscript

###############################################################################

library('getopt');

options(useFancyQuotes=F);

params=c(
	"input_olink_fn", "i", 1, "character",
	"output_fn_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input olink file>\n",
	"	-o <output filename root>\n",
	"\n",
	"This script will read in a olink file and extract out\n",
	"SampleID (rows), Assay (protein names, columns), and \n",
	"PCNormalizedNPX (cell values).\n",
	"\n", sep="");

if(
	!length(opt$input_olink_fn) || 
	!length(opt$output_fn_root)
){
	cat(usage);
	q(status=-1);
}

OlinkInputFile=opt$input_olink_fn;
OutputFileRoot=opt$output_fn;

cat("           Input Olink File: ", OlinkInputFile, "\n", sep="");
cat("           Output File Root: ", OutputFileRoot, "\n", sep="");
cat("\n");

options(width=100);
cat("Text Line Width: ", options()$width, "\n", sep="");

##############################################################################

extract_to_matrix=function(raw_dat, rown, coln, valn){
	sample_ids=unique(raw_dat[,rown]);
	categories=unique(raw_dat[,coln]);

	num_samples=length(sample_ids);
	num_categories=length(categories);

	cat("Num Samples: ", num_samples, "\n");
	cat("Num Categories: ", num_categories, "\n");

	temp_mat=matrix(NA, nrow=num_samples, ncol=num_categories);
	colnames(temp_mat)=categories;
	rownames(temp_mat)=sample_ids;

	num_raw_dat_lines=nrow(raw_dat);
	for(i in 1:num_raw_dat_lines){
		samp=raw_dat[i,rown];
		cat=raw_dat[i,coln];
		val=raw_dat[i,valn];
		temp_mat[samp, cat]=val;
	}

	return(temp_mat);
}

##############################################################################
# Load Raw Data 

ROWNAME="SampleID";
COLNAME="Assay";
VALNAME="PCNormalizedNPX";

cat("Loading Olink File:\n");

olinkin=read.table(OlinkInputFile, header=T, sep="\t");

coln=colnames(olinkin);
num_rows=nrow(olinkin);

cat("Number of Rows: ", num_rows, "\n");
cat("Header Names:\n");
print(coln);

# Move data to 2D matrix
mat=extract_to_matrix(raw_dat=olinkin, rown=ROWNAME, coln=COLNAME, valn=VALNAME);
#print(mat);

# Calculate stats
cat("Calc Stats:\n");
cat_avg=apply(mat, 2, function(x){mean(x,na.rm=T)});
cat_sd=apply(mat, 2, function(x){sd(x,na.rm=T)});
cat_wc=apply(mat, 2, function(x){x=x[!is.na(x)]; wilcox.test(x)$p.val});

# Sort by decreasing average
cat("Reorder Stats:\n");
order_ix=order(cat_avg, decreasing=T);
cat_avg=cat_avg[order_ix];
cat_sd=cat_sd[order_ix];
cat_wc=cat_wc[order_ix];
mat_order=mat[,order_ix,drop=F];

# Find NAs
num_nas=apply(mat_order, 1, function(x){sum(is.na(x))});
num_categories=ncol(mat_order);
cat("Num Categories: ", num_categories, "\n");
prop_na=num_nas/num_categories;
names(prop_na)=rownames(mat_order);

cat("Prop NA: ", prop_na, "\n");
print(prop_na);

# Extra samples with no NAs
kept_samples=(prop_na==0);
mat_order_nonas=mat_order[kept_samples,,drop=F];

##############################################################################
# Export Histograms

pdf(paste(OutputFileRoot, ".pdf", sep=""), height=11, width=8.5);

par(mfrow=c(4,3));
category_names=colnames(mat_order_nonas);
range=range(mat_order_nonas);
hist_bins=seq(range[1], range[2], length.out=20);
for(i in 1:num_categories){
	hist(mat_order_nonas[,i], 
		main=paste(i, ": ", category_names[i], sep=""), 
		xlab=VALNAME,
		breaks=hist_bins
		);
}

###############################################################################
# Export samples x categories

outtab=cbind(rownames(mat_order_nonas), mat_order_nonas);
cnames=c("SampleID", colnames(mat_order));
colnames(outtab)=cnames;
write.table(
	x=outtab,
	file=paste(OutputFileRoot, ".meta.tsv", sep=""),
	quote=F, sep="\t", row.names=F, col.names=T);

###############################################################################
# Export stats on each category

stats_tab=cbind(colnames(mat_order), cat_avg, cat_sd, cat_wc);
colnames(stats_tab)=c("Category", "Mean", "Stdev", "WilcoxPval");
write.table(
	x=stats_tab,
	file=paste(OutputFileRoot, ".stats.tsv", sep=""),
	quote=F, sep="\t", row.names=F, col.names=T);

###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
