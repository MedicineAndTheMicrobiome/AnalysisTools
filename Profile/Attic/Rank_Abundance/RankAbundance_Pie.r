#!/usr/local/bin/Rscript

###############################################################################

progname <- commandArgs(FALSE)[4]
args <- commandArgs(TRUE)

arg_count=1

if(is.na(args[arg_count])){

	script_name <- unlist(strsplit(progname,"="))[2]
	usage <- paste (
		"\nUsage:\n\t", script_name, "\n\t\t<Input summary_table.xls FileName>\n\n",
		"Generates a pie chart\n",
		"\n")

	writeLines(usage)
	writeLines("Input FileName not defined.\n")
	quit(status=0)
}

###############################################################################
# Main program loop

while(!(is.na(args[arg_count]))){
	InputFileName=args[arg_count]

	OutputFileRoot=InputFileName;
	OutputFileRoot=gsub("\\.summary_table\\.tsv$", "", OutputFileRoot);
	OutputFileRoot=gsub("\\.summary_table\\.xls$", "", OutputFileRoot);

	PiePDF= paste(OutputFileRoot, ".pie.pdf", sep="")
	PieXLS= paste(OutputFileRoot, ".pie.xls", sep="")

	cat("\n")
	cat("Input File Name: ", InputFileName, "\n")
	cat(" Pie Output PDF: ", PiePDF, "\n")
	cat(" Pie Output XLS: ", PieXLS, "\n")

	###############################################################################
	###############################################################################

	# Load data
	mat<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, row.names=1, comment.char=""));

	print(mat);

	print(dim(mat));

	if(is.null(ncol(mat))){
		mat=t(mat);
	}
	
	count_mat=mat[,2:ncol(mat)];

	if(is.null(ncol(count_mat))){
		count_mat=t(count_mat);
	}

	cat("Rows: ", nrow(count_mat), "\n");
	cat("Cols: ", ncol(count_mat), "\n");

	categories=as.vector(colnames(count_mat));
	num_categories=length(categories);

	short_names=character(num_categories);
	#for(i in 1:num_categories){
	#	taxonomy=unlist(strsplit(categories[i], " "));
	#	short_names[i]=taxonomy[length(taxonomy)];
	#}
	short_names=categories;

	#print(short_names);

	NumSamples=nrow(count_mat);
	NumCategories=ncol(count_mat);

	cat("\n");
	cat("Num Samples: ", NumSamples, "\n");
	cat("Num Categories: ", NumCategories, "\n");
	cat("\n");

	pdf(PiePDF, width=11,height=8.5)

	category_totals=numeric(NumCategories);
	for(cat_idx in 1:NumCategories){
		category_totals[cat_idx]=sum(count_mat[,cat_idx]);
	}

	all_total=sum(category_totals);
	
	totals_sort=sort(category_totals, index.return=TRUE, decreasing=TRUE);
	sorted_totals=totals_sort$x;

	normalized=sorted_totals/all_total;
	#print(normalized);
	lt1p=normalized[normalized >= 0 & normalized < 0.01];
	lt2p=normalized[normalized >= .01 & normalized < 0.02];
	lt3p=normalized[normalized >= .02 & normalized < 0.03];

	lt1p_sum=sum(lt1p);
	lt2p_sum=sum(lt2p);
	lt3p_sum=sum(lt3p);

	confined=normalized[normalized >=.03];
	confined=c(confined,lt3p_sum,lt2p_sum,lt1p_sum);
	
	sorted_shortnames=short_names[totals_sort$ix];
	confined_names=sorted_shortnames[normalized >=.03];
	confined_names=c(confined_names, "2-3%", "1-2%", "<1%");

	labels=character(length(confined_names));
	for(i in 1:length(confined_names)){
		labels[i]=sprintf("%s (%1.1f%%)", confined_names[i], confined[i]*100);
		print(labels[i]);
	}
	
	cdf=numeric(length(confined_names));
	cdf[1]=0;
	for(i in 1:(length(confined_names)-1)){
		cdf[i+1]=cdf[i]+confined[i];
	}	
	allcolors=rainbow(100, start=0, end=.65);
	colors=allcolors[as.integer(cdf*100+1)];

	pie(confined, labels, col=colors, init.angle=180, clockwise=TRUE,
		main=InputFileName);

	dev.off();

	###############################################################################
	# Output table

	fc=file(PieXLS, "w");
	outline=paste("Category Name", "Percentage", sep=",");
	write(outline, file=fc);

	for(i in 1:length(confined)){
		outline=paste(confined_names[i], sprintf("%.2f", confined[i]*100), sep=",");
		write(outline, file=fc);
	}

	###############################################################################

	arg_count=arg_count+1;
}

writeLines("Done.\n")

q(status=0)
