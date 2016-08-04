#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file_a", "a", 1, "character",
	"input_file_b", "b", 1, "character",
	"category_list", "t", 2, "character",
	"output_file", "o", 2, "character",
	"pvalue_cutoff", "p", 2, "numeric",
	"use_difference", "d", 2, "logical",
	"plot_histograms", "h", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-a <input summary table.xls>\n",
	"	-b <input summary table.xls>\n",
	"	[-t <list of taxa to compare between>]\n",
	"	[-o <output summary table file name>]\n",
	"	[-p <p-value cutoff, default=.05 >]\n",
	"	[-d (Flag to use difference (A-B) instead of (A/B)]\n",
	"	[-h (Flag to plot histograms]\n",
	"\n",	
	"This program will use the Wilcoxon rank sum test to compute the p-values\n",
	"between the proportions of the two sets of samples, based on all the\n",
	"taxa/categories specified in the list file.\n",
	"\n",
	"If the list file is not provide, then all pairwise comparisons will be made.\n",
	"\n");

if(!length(opt$input_file_a) || !length(opt$input_file_b)){
	cat(usage);
	q(status=-1);
}

UseDifferenceAsComparator=FALSE;
if(length(opt$use_difference)){
	UseDifferenceAsComparator=TRUE;
	comp_ext="diff";
	cat("Using difference as a comparator.\n");
}else{
	comp_ext="ratio";
	cat("Using ratios as a comparator.\n");
}

if(!length(opt$output_file)){
	root_a=gsub("\\.summary_table\\.xls$", "", opt$input_file_a);
	root_b=gsub("\\.summary_table\\.xls$", "", opt$input_file_b);
	root_a=tail(strsplit(root_a, "/")[[1]],1);
	root_b=tail(strsplit(root_b, "/")[[1]],1);
	OutputFileName = paste(root_a, "_vs_", root_b, ".", comp_ext, ".stats" , sep="");
}else{
	OutputFileName=opt$output_file;
}

###############################################################################

InputFileNameA=opt$input_file_a;
InputFileNameB=opt$input_file_b;
CategoryList=opt$category_list;
DoAll=FALSE;
if(length(CategoryList)==0){
	DoAll=TRUE;
}
	
PvalueCutoff=opt$pvalue_cutoff;
if(length(PvalueCutoff)==0){
	PvalueCutoff=.05;
}

PlotHistograms=FALSE;
if(length(opt$plot_histograms)){
	PlotHistograms=TRUE;
}

cat("\n")
cat("Input File Name A: ", InputFileNameA, "\n");
cat("Input File Name B: ", InputFileNameB, "\n");
cat("Output File Name: ", OutputFileName, "\n");       

if(!DoAll){
	cat("Category List: ", CategoryList, "\n");
}

cat("p-value cutoff (before correction for multiple testing): ", PvalueCutoff, "\n");
if(PvalueCutoff<0 || PvalueCutoff>1){
	cat("Invalid p-value cutoff specified.  Must be between 0 and 1.\n");
	q(status=-1);	
}

cat("\n");

###############################################################################
###############################################################################

load_data=function(input_filename){
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

load_list=function(input_filename){
	list=scan(input_filename, what="character", sep="\n");
	return(list);
}

#-----------------------------------------------------------------------------

compute_comparisons=function(matrix, category_list, diff){
	num_categories_to_analyze=length(category_list);

	num_samples=nrow(matrix);
	num_categories=ncol(matrix);
	results=array(0, dim=c(num_categories_to_analyze, num_categories_to_analyze, num_samples));

	stderr_fc=stderr();
	for(s in 1:num_samples){
		cat(".", file=stderr_fc);
		for(i in 1:num_categories_to_analyze){
			for(j in 1:num_categories_to_analyze){
				#cat("Computing: ", category_list[i], " / ", category_list[j], "\n");

				if(diff==TRUE){
					results[i,j,s]=matrix[s,category_list[i]]-matrix[s,category_list[j]];
				}else{
					results[i,j,s]=matrix[s,category_list[i]]/matrix[s,category_list[j]];
				}
			}
		}
	}

	cat("\n", file=stderr_fc);
	return(results);
}

#-----------------------------------------------------------------------------

extract_vector=function(array, cat1, cat2){
	return(array[cat1, cat2,]);
}

#-----------------------------------------------------------------------------

ci=function(x, ci){
	n=length(x);
	bound_idx=floor(n*(1-ci)/2);
	x_sorted=sort(x);
	conf_int=list();

	if(bound_idx>0){
		conf_int$ub=x_sorted[n-bound_idx];
		conf_int$lb=x_sorted[bound_idx+1];
	}else{
		conf_int$ub=NA;
		conf_int$lb=NA;
	}

	return(conf_int);
}
 
significant_pvalues=function(pvalues, alpha){

	# Get dimensions
	ncol=ncol(pvalues);
	nrow=nrow(pvalues);
	if(ncol!=nrow){
		cat("Error: matrix is not square.\n");
		q(status=-1);
	}

	# Adjust alpha for multiple comparisions
	num_tests=((ncol-1)*nrow)/2
	alpha_corrected=alpha/num_tests;
	cat("Num tests: ", num_tests, "\n");
	cat("Corrected alpha=", alpha_corrected, "\n");

	# Find signficant members
	significant_members=which(pvalues<=alpha_corrected);
	
	# Find coordinates of members
	x=((significant_members-1) %% nrow)+1;
	y=((significant_members-1) %/% nrow)+1;

	# Store, so we can return it
	val=pvalues[significant_members];
	significant=cbind(x,y,val);
	signif_rec=list();
	signif_rec$significant=significant;
	signif_rec$alpha_cor=alpha_corrected;

	return(signif_rec);	
}

remove_redundant=function(coordinates){
	nrows=nrow(coordinates);
	ncols=ncol(coordinates);

	if(ncols<2){
		cat("Error: input needs to consist of at least 2 columns");
		q(status=-1);
	}

	save=numeric();
	if(nrows>1){
		seen=list();
		for(i in 1:nrows){
			key=sprintf("%i,%i",coordinates[i,1], coordinates[i,2]);
			key_alt=sprintf("%i,%i",coordinates[i,2], coordinates[i,1]);
			if(is.null(seen[[key]]) && is.null(seen[[key_alt]])){
				seen[[key]]=1;
				save=rbind(save, coordinates[i,]);
			}
		}
	}

	return(save);
}


print_mat=function(matrix, name, categories, fc){
	nrow=nrow(matrix);
	ncol=ncol(matrix);

	cat(name, ":", sep="", file=fc);

	for(j in 1:ncol){
		cat(",[",j,": ", categories[j], "]",sep="",file=fc);
	}
	cat("\n", file=fc);
	
	for(i in 1:nrow){
		cat("[",i,": ", categories[i], "]", ",", sep="", file=fc);
		for(j in 1:ncol){
			cat(matrix[i,j], ",", sep="", file=fc);
		}
		cat("\n", file=fc);
	}
	cat("\n", file=fc);
}

abbreviate_categories=function(category_list, allowed_char){
	num_cats=length(category_list);
	abbreviated=character(0);
	
	for(i in 1:num_cats){
		last=tail(strsplit(category_list[i], " ")[[1]],1);
		if(nchar(last)<=allowed_char){
			abbreviated[i]=last;
		}else{
			no_vowels=paste(substr(last,1,1), gsub("[aeiouy]", "", substr(last,2,nchar(last))), sep="");
			abbreviated[i]=substr(no_vowels, 1, allowed_char);
		}
	}
	return(abbreviated);
}

###############################################################################

inmatA=load_data(InputFileNameA);
inmatB=load_data(InputFileNameB);

if(!DoAll){
	category_list=sort(load_list(CategoryList));
}else{
	category_list=intersect(colnames(inmatA), colnames(inmatB));
}

abbrev_categories=abbreviate_categories(category_list, 8);
print(abbrev_categories);
num_categories_to_analyze=length(category_list);
cat("Number of categories to compare between: ", num_categories_to_analyze, "\n", sep="");

normA=normalize(inmatA);
normB=normalize(inmatB);

cat("Computing comparisons for A:\n");
compA=compute_comparisons(normA, category_list,UseDifferenceAsComparator);
cat("Computing comparisons for B:\n");
compB=compute_comparisons(normB, category_list,UseDifferenceAsComparator);

#print(compA);
#print(compB);

mean_arrayA=matrix(1, nrow=num_categories_to_analyze, ncol=num_categories_to_analyze);
med_arrayA=matrix(1, nrow=num_categories_to_analyze, ncol=num_categories_to_analyze);
sd_arrayA=matrix(1, nrow=num_categories_to_analyze, ncol=num_categories_to_analyze);
lb_arrayA=matrix(1, nrow=num_categories_to_analyze, ncol=num_categories_to_analyze);
ub_arrayA=matrix(1, nrow=num_categories_to_analyze, ncol=num_categories_to_analyze);
nA=nrow(inmatA);

mean_arrayB=matrix(1, nrow=num_categories_to_analyze, ncol=num_categories_to_analyze);
med_arrayB=matrix(1, nrow=num_categories_to_analyze, ncol=num_categories_to_analyze);
sd_arrayB=matrix(1, nrow=num_categories_to_analyze, ncol=num_categories_to_analyze);
lb_arrayB=matrix(1, nrow=num_categories_to_analyze, ncol=num_categories_to_analyze);
ub_arrayB=matrix(1, nrow=num_categories_to_analyze, ncol=num_categories_to_analyze);
nB=nrow(inmatB);

pvalues=matrix(1, nrow=num_categories_to_analyze, ncol=num_categories_to_analyze);

if(PlotHistograms){
	pdf(paste(OutputFileName, ".pdf", sep=""), height=8.5, width=11);
}

for(i in 1:num_categories_to_analyze){
	cat(abbrev_categories[i], "\n");
	for(j in 1:num_categories_to_analyze){
		
		# Do not process self/self comparisons
		if(i==j){
			next;
		}

		title=paste(category_list[i], " / ", category_list[j], sep="");
		abbr_title=paste(abbrev_categories[i], " / ", abbrev_categories[j], sep="");
		#cat(abbr_title, "\n", sep="");

		# Grab vector of comparisons for each set of samples
		vectorsA=extract_vector(compA, i, j);
		vectorsB=extract_vector(compB, i, j);

		# Plot histogram
		if(PlotHistograms){	

			if(!UseDifferenceAsComparator){
				plot_vectorsA=log(vectorsA);
				plot_vectorsB=log(vectorsB);
				xlabel_str="Log(Ratio)";
			}else{
				plot_vectorsA=vectorsA;
				plot_vectorsB=vectorsB;
				xlabel_str="Difference";
			}

			ab=c(plot_vectorsA, plot_vectorsB);
			if(all(!is.finite(plot_vectorsA)) || all(!is.finite(plot_vectorsB))){
				next;
			}
			ab_hist=hist(ab, plot=F);
			#print(ab_hist);
			a_hist=hist(plot_vectorsA, breaks=ab_hist$breaks, plot=F);
			b_hist=hist(plot_vectorsB, breaks=ab_hist$breaks, plot=F);
			both=rbind(a_hist$counts, b_hist$counts);

			barplot(both, beside=TRUE, space=c(0,.1), 
				col=c("blue", "green"),  names.arg=ab_hist$mids, xlab=xlabel_str, ylab="Counts");
			mtext(paste("[", i, "] ", category_list[i]),  line=3);
			mtext("vs", line=2);
			mtext(paste("[", j, "] ", category_list[j]),  line=1);

			mtext(InputFileNameA, col="blue", side=4, line=-1, cex=.8);
			mtext(InputFileNameB, col="green", side=4, line=0, cex=.8);
		}

		# Compute stats on A
		mean_arrayA[i,j]=mean(vectorsA);
		med_arrayA[i,j]=median(vectorsA);
		sd_arrayA[i,j]=sd(vectorsA);
		ci95A=ci(vectorsA,.95);
		lb_arrayA[i,j]=ci95A$lb;
		ub_arrayA[i,j]=ci95A$ub;

		# Compute stats on B
		mean_arrayB[i,j]=mean(vectorsB);
		med_arrayB[i,j]=median(vectorsB);
		sd_arrayB[i,j]=sd(vectorsB);
		ci95B=ci(vectorsB,.95);
		lb_arrayB[i,j]=ci95B$lb;
		ub_arrayB[i,j]=ci95B$ub;

		# Compute Wilcoxon rank sum test, and keep the p-value
		#cat("[", i, "] vs [", j, "]\n");
		#print(vectorsA);
		#print(vectorsB);
		wilcox_result=wilcox.test(x=vectorsA, y=vectorsB, pair=FALSE);
		pvalues[i,j]=wilcox_result$p.value;
		#print(pvalues[i,j]);

	}
}

significant_comparisons=significant_pvalues(pvalues, PvalueCutoff);
print(significant_comparisons);
nonred_significant_comparisons=remove_redundant(significant_comparisons$significant);

if(0){
	cat("A:", InputFileNameA, "\n");
	cat("Mean\n");
	print(mean_arrayA);
	cat("Median\n");
	print(med_arrayA);
	cat("St. Dev\n");
	print(sd_arrayA);
	cat("95% Lowerbound\n");
	print(lb_arrayA);
	cat("95% Upperbound\n");
	print(ub_arrayA);
	cat("\n");

	cat("B:", InputFileNameB, "\n");
	cat("Mean\n");
	print(mean_arrayB);
	cat("Median\n");
	print(med_arrayB);
	cat("St. Dev\n");
	print(sd_arrayB);
	cat("95% Lowerbound\n");
	print(lb_arrayB);
	cat("95% Upperbound\n");
	print(ub_arrayB);
	cat("\n");

	cat("P-values:\n");
	print(pvalues);

	print(significant_comparisons);
	print(nonred_significant_comparisons);
}

###############################################################################
# Output
fc=file(paste(OutputFileName, ".csv", sep=""), "w");

# Category lists
cat(file=fc, "Categories: \n");
for(i in 1:num_categories_to_analyze){
	cat(file=fc, "[", i, "] (", abbrev_categories[i], ") ", category_list[i], "\n", sep=""); 
}

cat(file=fc, "\n\n\n");
cat(InputFileNameA, "\n", sep="", file=fc);
cat("N=", nA, "\n", sep="", file=fc);
print_mat(mean_arrayA, "Mean", abbrev_categories, fc);
print_mat(med_arrayA, "Median", abbrev_categories, fc);
print_mat(sd_arrayA, "St. Dev", abbrev_categories, fc);
print_mat(lb_arrayA, "LB 95%", abbrev_categories, fc);
print_mat(ub_arrayA, "UB 95%", abbrev_categories, fc);

cat(file=fc, "\n\n\n");
cat(InputFileNameB, "\n", sep="", file=fc);
cat("N=", nB, "\n", sep="", file=fc);
print_mat(mean_arrayB, "Mean", abbrev_categories, fc);
print_mat(med_arrayB, "Median", abbrev_categories, fc);
print_mat(sd_arrayB, "St. Dev", abbrev_categories, fc);
print_mat(lb_arrayB, "LB 95%", abbrev_categories, fc);
print_mat(ub_arrayB, "UB 95%", abbrev_categories, fc);

cat(file=fc, "\n\n\n");
cat(file=fc, "\n");
cat("Wilcoxon Rank Sum Test:\n", file=fc);
print_mat(pvalues, "p-values", abbrev_categories, fc);

cat(file=fc, "\n\n\n");
cat(file=fc, "Corrected alpha = ", significant_comparisons$alpha_cor, "\n");
num_nonred_sig_comp=nrow(nonred_significant_comparisons);
if(length(num_nonred_sig_comp)==0){
	num_nonred_sig_comp=0;
}
cat("Num statistically significant comparisions: ", num_nonred_sig_comp, "\n");
cat(file=fc, "Num Non-Redundant Significant Comparisions = ", num_nonred_sig_comp, "\n");
cat(file=fc, "Idx,CatA,CatB,NameA,NameB,Uncor p-val,MedA,MedB\n");
if(num_nonred_sig_comp>0){
	for(i in 1:num_nonred_sig_comp){
		x_idx=nonred_significant_comparisons[i,1];
		y_idx=nonred_significant_comparisons[i,2];

		#vectorsA=extract_vector(compA, x_idx, y_idx);
		#vectorsB=extract_vector(compB, x_idx, y_idx);
		
		#cat("A:\n");
		#print(vectorsA);
		#cat("B:\n");
		#print(vectorsB);

		cat(file=fc, 
			paste("<", i,">", sep=""),
			x_idx,
			y_idx,
			abbrev_categories[x_idx],
			abbrev_categories[y_idx],
			nonred_significant_comparisons[i,3],
			med_arrayA[x_idx,y_idx],
			med_arrayB[x_idx,y_idx],
			sep=","
		);
		cat(file=fc, "\n");
	}
}

close(fc);

###############################################################################

writeLines("Done.\n")
#print(warnings());

q(status=0)
