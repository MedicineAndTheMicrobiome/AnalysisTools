#!/usr/local/bin/Rscript

###############################################################################

progname <- commandArgs(FALSE)[4]
args <- commandArgs(TRUE)

arg_count=1

if(is.na(args[arg_count])){

	script_name <- unlist(strsplit(progname,"="))[2]
	usage <- paste (
		"\nUsage:\n\t", script_name, "\n\t\t<Input summary_table.xls FileName>\n\n",
		"  [ -r (to remove unknowns)]\n",
		"\n");

	writeLines(usage)
	writeLines("Input FileName not defined.\n")
	quit(status=0)
}

###############################################################################
# Main program loop

REMOVE_UNKNOWNS=FALSE;
options=logical(0);
for(i in 1:length(args)){
       options[i]=FALSE;
       if(args[i]=="-r"){
               cat ("Removing Unknowns options set.\n");
               REMOVE_UNKNOWN=TRUE;
               options[i]=TRUE;
       }
}
args=args[!options];

while(!(is.na(args[arg_count]))){
	InputFileName=args[arg_count]

	OutputFileRoot=InputFileName;
	OutputFileRoot=gsub("\\.summary_table\\.tsv$", "", OutputFileRoot);
	OutputFileRoot=gsub("\\.summary_table\\.xls$", "", OutputFileRoot);

	if(REMOVE_UNKNOWNS){
		RAPlot= paste(OutputFileRoot, ".rank_abundancy.boxplot.no_unknowns.pdf", sep="")
		RAXLS= paste(OutputFileRoot, ".rank_abundancy.boxplot.no_unknowns.xls", sep="")
	}else{
		RAPlot= paste(OutputFileRoot, ".rank_abundancy.boxplot.pdf", sep="")
		RAXLS= paste(OutputFileRoot, ".rank_abundancy.boxplot.xls", sep="")
	}

	cat("\n");
	cat("             Input File Name: ", InputFileName, "\n")
	cat("Rank Abundancy Box Plots PDF: ", RAPlot, "\n")
	cat("\n");

	###############################################################################
	###############################################################################

	# Load data
	inmat=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, row.names=1, check.names=FALSE, comment.char="", quote=NULL))
	#cat("Original Matrix:\n")
	#print(inmat);

	# Get input dimensions
	num_orig_samples=nrow(inmat);
	num_orig_categories=ncol(inmat)-2;
	cat("Num samples:", num_orig_samples, "\n");
	cat("Num categories:", num_orig_categories, "\n");

	# Get sample/category names
	#orig_sample_names=inmat[,1];
	orig_sample_names=row.names(inmat)[1:(nrow(inmat))];
	
	#print(orig_sample_names);
	orig_category_names=colnames(inmat)[3:(num_orig_categories+2)];
	#print(category_names);

	# Only keep the taxon with the finest classification
	shortened_category_names=rep("",num_orig_categories);
	for(catidx in 1:num_orig_categories){
		#print(category_names[catidx]);
		components=unlist(strsplit(orig_category_names[catidx], " "));
		#shortened_category_names[catidx]=tail(components,1);
		shortened_category_names[catidx]=orig_category_names[catidx];
	}
	#print(shortened_category_names);
	
	# Convert counts into integers
	orig_counts=matrix(nrow=num_orig_samples, ncol=num_orig_categories);
	orig_counts[1:num_orig_samples,1:num_orig_categories]=as.numeric(inmat[1:num_orig_samples,3:(num_orig_categories+2)]);
#	print(orig_counts);

	# Compute sample totals
	orig_totals=rep(0,num_orig_samples);
	for(sampidx in 1:num_orig_samples){
		orig_totals[sampidx]=sum(orig_counts[sampidx,]);
	}
	#print(orig_totals);

	# Normalized counts by sample totals.
	orig_normalized=matrix(nrow=num_orig_samples, ncol=num_orig_categories);
	for(sampidx in 1:num_orig_samples){
		orig_normalized[sampidx,]=orig_counts[sampidx,]/orig_totals[sampidx];
	}
#	print(orig_normalized);
	#print(orig_counts);

	# Remove any samples with 0 counts
	nonzero_samples=orig_totals>0;
	nonzero_samples_normalized=orig_normalized[nonzero_samples,];
	nonzero_samples_names=orig_sample_names[nonzero_samples];
	num_nonzero_samples=nrow(nonzero_samples_normalized);
	#print(nonzero_samples);
	#print(nonzero_samples_names);
	#print(orig_sample_names);

	# Compute median percentages and sort
	median=rep(0,num_orig_categories);
	for(catidx in 1:num_orig_categories){
		median[catidx]=median(nonzero_samples_normalized[,catidx]);
	}	
	median_sort=sort(median,decreasing=TRUE,index.return=TRUE);
	sorted_median=median_sort$x;
	#print(sorted_median);
	#print(sum(sorted_median));

	# Re-order the matrix by median normalized counts 
	nonzero_sample_normalized_sorted=nonzero_samples_normalized[,median_sort$ix];
	category_names_sorted=shortened_category_names[median_sort$ix];
	#print(nonzero_sample_normalized_sorted);
	#print(category_names_sorted);

	# Apply a cutoff so that infrequently represented categories can be ignored
	cutoff=0.000;
	selected_categories_normalized=nonzero_sample_normalized_sorted[,sorted_median>cutoff];
	selected_sorted_median=sorted_median[sorted_median>cutoff];
	num_selected_categories=ncol(selected_categories_normalized);
	#print(num_selected_categories);
	selected_category_names=category_names_sorted[sorted_median>cutoff];

	# If the median abundance for all taxa is 0, then just take the top 10.
	if(num_selected_categories==0){
		selected_categories_normalized=nonzero_sample_normalized_sorted[,1:10];
		selected_sorted_median=sorted_median[1:10];
		num_selected_categories=ncol(selected_categories_normalized);
		selected_category_names=category_names_sorted[1:10];
		cat("\n\tWARNING: The median for all categories is 0.  Reporting top 10 taxa.\n\n");
	}

	if(REMOVE_UNKNOWNS){
		keep=logical(0);
		for(i in 1:num_selected_categories){
			split_result=strsplit(selected_category_names[i], ":");
			keep[i]=(length(split_result[[1]])==1);
		}
		selected_sorted_median=selected_sorted_median[keep];
		selected_categories_normalized=selected_categories_normalized[,keep];
		selected_category_names=selected_category_names[keep];
		num_selected_categories=sum(keep);
	}

	#print(selected_category_names);
	cat("Num selected categories: ", num_selected_categories, "\n");

	# Compute colors
	total_median=sum(sorted_median);
	normalized_median=sorted_median/total_median;
	cdf=rep(0,num_selected_categories);
	#print(num_selected_categories);
	for(i in 1:num_selected_categories){
		#cat("i =", i, "\n");
		cdf[i+1]=cdf[i]+normalized_median[i];
	}
	colors=rainbow(100,start=0,end=.65)[as.integer(cdf*100+1)];
	#print(cdf);


	# Generate overall box plot
	pdf(RAPlot, height=8.5, width=11);
	label_scale=40/num_selected_categories;
	if(label_scale>1){label_scale=1;}
	par(mar=c(10,4,3,2)); # bottom, left, top, right)

	shift=.5;
	boxplot(selected_categories_normalized[,1], at=1-shift, show.names=TRUE, 
		boxwex=.8, las=2, cex.axis=label_scale, cex=.7, col=colors[1],
		names=selected_category_names[1], xlim=c(0, num_selected_categories), ylim=c(0,1),
	    	main="Rank Abundancy"
		);
	for(i in 2:num_selected_categories){
		boxplot(selected_categories_normalized[,i], add=TRUE, at=i-shift, show.names=TRUE,
		boxwex=.8, las=2, cex.axis=label_scale, cex=.7,  col=colors[i],
		names=selected_category_names[i],
		yaxt="n"
		);
	}
	#boxplot(data.frame(selectedn_sort_categories), 
	#	    boxwex=.8, las=2, cex.axis=label_scale, cex=.7, 
	#	    main=paste("Rank Abundancy ( proportion >", cutoff, ")"));
	mtext(sprintf("Median Proportion > %i", cutoff), side=3, at=-.5, adj=0);


	######################################################################################
	# Plot top 10
	top_subset=10;
	if(num_selected_categories<top_subset){
		top_subset=num_selected_categories;
	}

	label_scale=40/top_subset;
	if(label_scale>1){label_scale=1;}
	par(mar=c(10,4,3,2)); # bottom, left, top, right)

	shift=.5;
	boxplot(selected_categories_normalized[,1], at=1-shift, show.names=TRUE, 
		boxwex=.8, las=2, cex.axis=label_scale, cex=.7, col=colors[1],
		names=selected_category_names[1], xlim=c(0, top_subset), ylim=c(0,1),
	    	main=paste("Rank Abundancy")
		);
	for(i in 2:top_subset){
		boxplot(selected_categories_normalized[,i], add=TRUE, at=i-shift, show.names=TRUE,
		boxwex=.8, las=2, cex.axis=label_scale, cex=.7,  col=colors[i],
		names=selected_category_names[i],
		yaxt="n");
	}
	mtext(sprintf("Top %i of Median Proportion>%i", top_subset, cutoff), side=3, at=-.5, adj=0);
	
	med=numeric(0);
	q1=numeric(0);
	q3=numeric(0);
	for(i in 1:top_subset){
		smry=summary(selected_categories_normalized[,i]);
		med=smry["Median"];	
		q1=smry["1st Qu."];	
		q3=smry["3rd Qu."];	

		medst=sprintf("%.2f",smry["Median"]);	
		q1st=sprintf("%.2f",smry["1st Qu."]);	
		q3st=sprintf("%.2f",smry["3rd Qu."]);	

		text(i-.10,med,labels=medst, cex=.9, srt=0, font=2);
		text(i-.24,q3+.015,labels=q3st, cex=.7, srt=0, font=3);
		text(i-.24,q1-.015,labels=q1st, cex=.7, srt=0, font=3);
	}
	


	######################################################################################

	cat("Plotting sample distributions...\n");
	sort_samples=sort(nonzero_samples_names, index.return=TRUE);
#	print(nonzero_samples_names);
	for(i in sort_samples$ix){

		output_names=rep("",num_selected_categories);
		for(catidx in 1:num_selected_categories){
			output_names[catidx]=paste(selected_category_names[catidx], 
				" [", sprintf("%0.4f",selected_categories_normalized[i,catidx]), "]", sep="");
		}

		barcenters=barplot(selected_categories_normalized[i,], main=nonzero_samples_names[i], col=colors,
		    ylim=c(0,1), las=2, cex.axis=label_scale, cex=.7,names=output_names);
		points(barcenters[,1],selected_sorted_median, pch="+", col="black");
	}

	dev.off();

	###############################################################################

	cat("Generating spreadsheet...\n");
	fc=file(RAXLS, "w");
	outline=paste("Name", "Minimum", "1st_Quartile", "Median", "Mean", "3rd_Quartile", "Maximum", sep=",");
	write(outline, file=fc);
	for(i in 1:num_selected_categories){
		smry=summary(selected_categories_normalized[,i]);
		outline=paste(
			selected_category_names[i],
			smry[["Min."]],
			smry[["1st Qu."]],
			smry[["Median"]],
			smry[["Mean"]],
			smry[["3rd Qu."]],
			smry[["Max."]],
			sep=",");
		write(outline, file=fc);
		
	}
	

	###############################################################################

	arg_count=arg_count+1;
}

writeLines("Done.\n")

q(status=0)
