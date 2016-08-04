#!/usr/bin/env Rscript

###############################################################################

progname <- commandArgs(FALSE)[4]
args <- commandArgs(TRUE)

arg_count=1

if(is.na(args[arg_count])){

	script_name <- unlist(strsplit(progname,"="))[2]
	usage <- paste (
		"\nUsage:\n\t", script_name, "\n\t\t<Input FileName>\n\n",
		"Reads in summary_table.xls and generates a histogram of the sample distributions across RDP hits.",
		"\n")

	writeLines(usage)
	writeLines("Input FileName not defined.\n")
	quit(status=0)
}

###############################################################################
# Main program loop

while(!(is.na(args[arg_count]))){
	InputFileName=args[arg_count]
	FileNameSplit<-unlist(strsplit(InputFileName,"/"))
	FileName<-sub(".summary_table.xls","",FileNameSplit[length(FileNameSplit)])
	OutputFileName=paste(FileNameSplit[1:(length(FileNameSplit)-1)],collapse="/")	
	OutputFileName=paste(OutputFileName,FileName,sep="/")
	HistogramPlotPDF = paste(OutputFileName,".hist_plot.pdf", sep="")

	cat("\n")
	cat("             Input File Name: ", InputFileName, "\n")
	cat("    	  Histogram Plot PDF: ", HistogramPlotPDF, "\n")

	# Example input:List of total species per sample

	# 4443
	# 553254
	# 1334
	# 85676
	# 568413
	# 235
	# 42444
	# 54


	###############################################################################
	###############################################################################

	# Load data
	counts_table=read.table(file=InputFileName, sep="\t", header=TRUE);
	A=counts_table[,2];
	#cat("Original Matrix:\n")
	#print(A)


	###############################################################################
	###############################################################################
	# Draw Histogram Plot
	pdf(HistogramPlotPDF,width=8.5,height=11)

	# hist
	hist(A, br=15, main=FileName, xlab="Counts per Sample", ylab="Frequency of Sample");

	dev.off();

	###############################################################################
	##############################################################################

	arg_count=arg_count+1;
}

writeLines("Done.\n")

q(status=0)
