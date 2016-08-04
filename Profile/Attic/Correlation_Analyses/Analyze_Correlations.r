#!/usr/local/bin/Rscript

###############################################################################

progname <- commandArgs(FALSE)[4]
args <- commandArgs(TRUE)

arg_count=1

if(is.na(args[arg_count])){

	script_name <- unlist(strsplit(progname,"="))[2]
	usage <- paste (
		"\nUsage:\n\t", script_name, "\n\t\t<Index file>\n\n",
		"Format is:\n",
		"SampleID\\tx-Axis\\tIndices...\\n",
		"\n")

	writeLines(usage)
	writeLines("Input FileName not defined.\n")
	quit(status=0)
}

###############################################################################
# Main program loop

while(!(is.na(args[arg_count]))){
	InputFileName=args[arg_count]

	CorrelPlot= paste(InputFileName, ".correlation.pdf", sep="")

	cat("\n")
	cat(" Input File Name: ", InputFileName, "\n")
	cat("Correlation Plot: ", CorrelPlot, "\n")

	###############################################################################
	###############################################################################

	# Load data
	in_tab<-read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, row.names=1)
	#cat("Original Matrix:\n")
	#print(in_tab);

	# Get input dimensions
	num_orig_samples=nrow(in_tab);
	num_orig_categories=ncol(in_tab);
	cat("Num samples:", num_orig_samples, "\n");
	cat("Num categories:", num_orig_categories, "\n");

	# Determine which column should be on the x-axis
	XAXIS_COL=1
	xaxis_name= colnames(in_tab)[XAXIS_COL];
	cat("Using ", xaxis_name, " on x-axis\n", sep="");
	mat=as.matrix(in_tab);
	x=mat[,XAXIS_COL];
	
	# Start the plot
	pdf(CorrelPlot, height=8.5, width=8.5);
	for(plot_idx in 1:num_orig_categories){

		if(plot_idx==XAXIS_COL){
			next;
		}

		plot_name=colnames(in_tab)[plot_idx];
		cat("Working on: ", plot_name, "\n", sep="");

		y=mat[,plot_idx];

		symsize=100/num_orig_samples;
		if(symsize>1){
			symsize=1;
		}
		plot(x, y, xlab=xaxis_name, ylab=plot_name, main=paste(xaxis_name, "vs", plot_name, sep=" "), cex=symsize);
		lin_mod=lm(y~x);
		
		#print(lin_mod);
		abline(lin_mod, col="grey", lty=2);
		cor_coef=cor(x,y);
		cat("Correl: ", cor_coef, "\n");
		
		yintercept=lin_mod$coefficients[1];		
		slope=lin_mod$coefficients[2];		

		cat("slope: ", slope, "\n");
		cat("y-intercept:", yintercept, "\n");

		model_summary=summary(lin_mod);
		rsquared=model_summary$r.squared;

		mtext(paste("Correlation Coefficient:  ", sprintf("%.2f", cor_coef)), side=3, line=0, cex=.9);
		mtext(paste("R^2 (unadjusted):  ", sprintf("%.2f", rsquared)), side=3, line=-1, cex=.6);
		mtext(paste("slope:  ", sprintf("%.2f", slope)), side=3, line=-1.5, cex=.6);
		mtext(paste("y-intercept:  ", sprintf("%.2f", yintercept)), side=3, line=-2, cex=.6);

		cat("\n");
	}

	dev.off();

	###############################################################################

	arg_count=arg_count+1;
}

writeLines("Done.\n")

q(status=0)
