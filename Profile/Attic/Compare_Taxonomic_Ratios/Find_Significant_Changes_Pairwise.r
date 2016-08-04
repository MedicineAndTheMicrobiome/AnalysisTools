#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"association_file", "a", 1, "character",
	"output_file", "o", 2, "character",
	"p_value", "p", 2, "numeric",
	"no_correction", "d", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input summary table.xls>\n",
	"	-a <input association file>\n",
	"	[-p <family wise acceptable error rate, default p-value = 0.05>]\n",
	"	[-o <output file name>]\n",
	"	[-d <turn off correction for multiple testing>]\n",
	"\n",	
	"This script will go through all the pairs in your association\n",
	"file and see if there is any statistically significant changes\n",
	"between the pairs among all associations.\n",
	"\n",
	"\n");

if(!length(opt$input_file)||!(length(opt$association_file))){
	cat(usage);
	q(status=-1);
}

###############################################################################

InputFileName=opt$input_file;
OutputFilename=opt$output_file;
AssociationFile=opt$association_file;
Pvalue=opt$p_value;
NoCorrection=TRUE;

if(length(opt$output_file)==0){
	OutputFilenameCSV=paste(gsub("\\.summary_table\\.xls$", "", InputFileName), ".paired_ratios_diff.csv", sep="");
	OutputFilenamePDF=paste(gsub("\\.summary_table\\.xls$", "", InputFileName), ".paired_ratios_diff.pdf", sep="");
}

if(length(opt$p_value)==0){
	Pvalue=0.05;
}

if(length(opt$no_correction)==0){
	NoCorrection=FALSE;
}

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFilenameCSV, "\n");
cat("Association File Name: ", AssociationFile, "\n");
cat("Family-wise P-value: ", Pvalue, "\n");
cat("Correction is turned off: ", NoCorrection, "\n");

###############################################################################
###############################################################################

load_summary_file_table=function(input_filename){
        # Load data
        cat("Loading: ", input_filename, "\n");
        inmat=as.matrix(read.table(input_filename, sep="\t", header=TRUE, check.names=FALSE, comment.char="*", row.names=1))
        #cat("Original Matrix:\n")
        #print(inmat);

        # Grab columns we need into a vector, ignore totals, we won't trust it.
        counts_mat=inmat[,2:(ncol(inmat))];
        #print(counts_mat);

        # Summary what we've loaded
        num_samples=nrow(counts_mat);
        num_categories=ncol(counts_mat);
        cat("Num Samples: ", num_samples, "\n");
        cat("Num Categories: ", num_categories, "\n");

        return(counts_mat);
}

#-----------------------------------------------------------------------------

normalize=function(inmat){
        samp_sums=apply(inmat, 1, sum);
        norm_mat=matrix(0, nrow=nrow(inmat), ncol=ncol(inmat));
        for(i in 1:length(samp_sums)){
                norm_mat[i,]=inmat[i,]/samp_sums[i];
        }

        rownames(norm_mat)=rownames(inmat);
        colnames(norm_mat)=colnames(inmat);

        #cat("Normalized:\n");
        #print(norm_mat);
        return(norm_mat);

}

#-----------------------------------------------------------------------------

load_associations=function(input_filename){
        cat("Loading association file: ", input_filename, "\n", sep="");

        inmat=as.matrix(read.table(input_filename, sep=",", header=TRUE, check.names=FALSE, row.names=1))

        associations=list();
        associations$sets=inmat;
        associations$num_groups=nrow(inmat);
        associations$num_types=ncol(inmat);
        associations$type_names=colnames(inmat);
        associations$group_ids=rownames(inmat);

        return(associations);
}

#-----------------------------------------------------------------------------

associate_pairs=function(num_groups){
	nrows=((num_groups-1)*num_groups)/2;
	pairs=matrix(0,nrow=nrows,ncol=2);
	k=1;
	for(i in 1:(num_groups-1)){
		for(j in (i+1):num_groups){
			pairs[k,1]=i;
			pairs[k,2]=j;
			k=k+1;
		}
	}
	return(pairs);
}

#-----------------------------------------------------------------------------

combined_histogram=function(data1, data2, name1, name2, title, comments){
	bothhist=hist(c(data1,data2),plot=FALSE)
	ahist=hist(data1,breaks=bothhist$breaks,plot=FALSE)
	bhist=hist(data2,breaks=bothhist$breaks,plot=FALSE)
	bothmat=t(matrix(c(ahist$intensities, bhist$intensities),length(bothhist$breaks)-1,2))
	dimnames(bothmat)=list(NULL, bothhist$breaks[1:(length(bothhist$breaks)-1)])
	barplot(bothmat,beside=TRUE,col=c("blue","red"), ylim=c(0,max(bothmat)*1.1),space=c(0,.3), main=title)
	median1=sprintf("%3.4e",median(data1));
	median2=sprintf("%3.4e",median(data2));
	mtext( paste(name1, ": ", median1, sep=""), col="blue", side=3, line=-1, cex=.7);
	mtext( paste(name2, ": ", median2, sep=""), col="red", side=3, line=-2, cex=.7);
	
	for(i in 1:length(comments)){
		mtext(comments[[i]], col="black", side=3, line=-3-i+1, cex=.7);
	}
}


###############################################################################

# Load associations
associations=load_associations(AssociationFile);
cat("Num Groups (N): ", associations$num_groups, "\n");
#print(associations);

# Load data
in_matrix=load_summary_file_table(InputFileName);
#cat("Original Matrix:\n")
#print(in_matrix);
category_names=colnames(in_matrix);
sample_names=rownames(in_matrix);
num_samples=nrow(in_matrix);
num_categories=ncol(in_matrix);

# Normalize
normalized=normalize(in_matrix);
#cat("Normalized Matrix:\n");
#print(normalized);

###############################################################################

# Compute association pairs
association_pairs=associate_pairs(associations$num_types);
#print(association_pairs);

num_comparisons=nrow(association_pairs)*length(category_names);
cat("Num comparisons to be made: ", num_comparisons, "\n");
corrected_alpha=Pvalue/num_comparisons;
cat("Corrected alpha = ", corrected_alpha, "\n");

kept_association=numeric();
kept_category_id=character();
kept_pvalues=numeric();
kept_counter=1;

# Go through every pair of associations
for(association in 1:nrow(association_pairs)){

	type1=association_pairs[association,1];
	type2=association_pairs[association,2];

	cat("Looking at relationship between: ", associations$type_names[type1], " and ", associations$type_names[type2], "\n", sep="");

	type1_samples=associations$sets[,type1];
	type2_samples=associations$sets[,type2];
	#print(type1_samples);
	#print(type2_samples);

	# Go through every category, ie. taxonomic classification
	for(category_id in category_names){
		#cat("\tAnalysing ", category_id, "\n", sep="");
		type1_values=normalized[type1_samples, category_id];
		type2_values=normalized[type2_samples, category_id];
		#print(type1_values);
		#print(type2_values);

		# Compute wilcoxon rank sum test
		wilcox_result=wilcox.test(type1_values, type2_values, paired=TRUE);
		#print(wilcox_result);
		pval=wilcox_result$p.value;

		# Test p-value against corrected
		if(!is.nan(pval)){
			if(pval<=corrected_alpha || (NoCorrection && pval<=Pvalue)){
				type1_median=median(type1_values);
				type2_median=median(type2_values);
				cat("\t", category_id, ":", pval, "\n", sep="");
				cat("\t\t", type1_median, " vs ", type2_median, "\n");

				# Store results of interest
				kept_association[kept_counter]=association;
				kept_category_id[kept_counter]=category_id;
				kept_pvalues[kept_counter]=pval;
				kept_counter=kept_counter+1;			
			}
		}

	}
}

#-----------------------------------------------------------------------------

num_kept=kept_counter-1;

fc=file(OutputFilenameCSV, "w");
pdf(OutputFilenamePDF, height=8.5, width=11);

unique_assoc=unique(kept_association);

for(assoc in unique_assoc){

	assoc_indcs=which(kept_association==assoc);

	# Get associations, eg. Lesion vs Normal
	type1=association_pairs[assoc,1];
	type2=association_pairs[assoc,2];
	cat(file=fc, associations$type_names[type1], " vs ", associations$type_names[type2], "\n", sep="");

	par(mfrow=c(1,1));
	plot(0,0, xlab="", ylab="", xaxt="n", yaxt="n",type="n", bty="n", oma=c(0,0,0,0), mar=c(0,0,0,0));
	text(0,0, associations$type_names[type1], pos=3, offset=1.3, cex=4, font=2);
	text(0,0, " vs ", cex=3, font=1);
	text(0,0, associations$type_names[type2], pos=1, offset=1.3, cex=4, font=2);
	par(mfrow=c(2,3));

	# Only look at data within associations
	categories=kept_category_id[assoc_indcs];
	pvalues=kept_pvalues[assoc_indcs];

	# Sort the categories by increasing pvalues
	sort_result=sort(pvalues, decreasing=FALSE, index.return=TRUE);

	# For each top hit
	for(sortedix in sort_result$ix){
		type1_samples=associations$sets[,type1];
		type2_samples=associations$sets[,type2];
		type1_values=normalized[type1_samples, categories[sortedix]];
		type2_values=normalized[type2_samples, categories[sortedix]];
		type1_median=median(type1_values);
		type2_median=median(type2_values);		

		combined_histogram(
			type1_values,type2_values, 
			associations$type_names[type1], associations$type_names[type2],
			title=tail(strsplit(categories[sortedix], " ")[[1]],1),
			comment=list(
				paste("Uncorrected p-value: ", sprintf("%3.4e",pvalues[sortedix])),
				paste("  Corrected p-value: ", sprintf("%3.4e",pvalues[sortedix]*num_comparisons))
			)
		);

		cat(file=fc, "\t", categories[sortedix], ", ", pvalues[sortedix], ", ", type1_median, ", ", type2_median, "\n", sep="");
	}
	
}


###############################################################################

cat("Done.\n")
#print(warnings());

q(status=0)
