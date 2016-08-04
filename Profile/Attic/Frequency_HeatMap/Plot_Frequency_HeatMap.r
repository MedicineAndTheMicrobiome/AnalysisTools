#!/usr/local/bin/Rscript

###############################################################################

progname <- commandArgs(FALSE)[4]
args <- commandArgs(TRUE)

arg_count=1

if(is.na(args[arg_count])){

	script_name <- unlist(strsplit(progname,"="))[2];
	usage <- paste (
		"\nUsage:\n\t", script_name, "\n\t\t<Input FileName>\n\n",
		"Reads in a table and plots a frequency heat map. No reordering is done.\n",
		"\n")

	writeLines(usage)
	writeLines("Input FileName not defined.\n")
	quit(status=0)
}

###############################################################################
# Main program loop

while(!(is.na(args[arg_count]))){

	InputFileName=args[arg_count];
	FrequencyHeatMapPDF= paste(InputFileName, ".freq_heatmap.pdf", sep="");

	cat("\n");
	cat("            Input File Name: ", InputFileName, "\n");
	cat("Frequency Heat Map PDF File: ", FrequencyHeatMapPDF, "\n");
	cat("\n");

	###############################################################################
	###############################################################################

	# Load data
	mat<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, row.names=1, check.names=FALSE));
	#print(mat);
	ncols=ncol(mat);
	nrows=nrow(mat);

	# Scale labels (so that they don't overlap)
	cscale=42/ncols;
	rscale=42/nrows;
	if(cscale>1){
		cscale=1;
	}
	if(rscale>1){
		rscale=1;
	}

	###############################################################################
	###############################################################################
	# Draw Hierarchical Clustering Plot
	pdf(FrequencyHeatMapPDF,width=8.5,height=11)

	# Pad labels with space (so it looks better)
	cnames=character();
	rnames=character();
	for(i in 1:ncols){
		cnames[i]=paste(colnames(mat)[i], "  ", sep="");
	}
	for(i in 1:nrows){
		rnames[i]=paste(rownames(mat)[i], "  ", sep="");
	}

	# Estimate margins based on max label length
	left_mar=max(nchar(cnames))/1.5;
	bottom_mar=max(nchar(rnames))/1.25;
	par(mar=c(bottom_mar,left_mar,.5,.5)); #b,l,t,r
	
	# Pick colors from blue to red, only.
	colors=(rainbow(2^16, start=0, end=0.65));
	#colors=rev(rainbow(2^16, start=0, end=0.65));

	# Layout
	#laymat=matrix(c(rep(1,18),rep(2,2)),1,20);
	HM_WID=18;
	laymat=matrix(rep(1,HM_WID^2),nrow=HM_WID, ncol=HM_WID); # heatmap area
	laymat=rbind(rep(2, HM_WID), laymat); #top margin
	#laymat=rbind(rep(2, HM_WID), laymat); #top margin
	laymat=cbind(laymat, c(5,rep(3, HM_WID)) ); # Left margin
	#laymat=cbind(laymat, c(5,5,rep(3, HM_WID)) ); # Left margin
	laymat=cbind(laymat, 4 ); # scale
	laymat=cbind(laymat, 4 ); # scale
	#laymat=cbind(laymat, 4 ); # scale
	layout(laymat);

	nanfree_mat=mat[!is.na(mat)];
	matmax=max(nanfree_mat);
	matmin=min(nanfree_mat);

	# Plot matrix
	image(1:nrow(mat),1:ncol(mat), mat,
		xaxt="n", yaxt="n",
		xlab="", ylab="",
		col=colors
	);	

	# Plot values
	max_width=max(nchar(sprintf("%.4g",mat)));
	for(i in 1:nrow(mat)){
		for(j in 1:ncol(mat)){
			text(i,j,labels=sprintf("%.4g",mat[i,j]), cex=(2/max_width)*(cscale+rscale), srt=45);
		}
	}

	# Plot the labels
	mtext(cnames, at=1:ncols, side=2, las=2, cex=cscale);
	mtext(rnames, at=1:nrows, side=1, las=2, cex=rscale);

	print("Main Done.\n");

	# Plot top row margin
	row_mean=numeric(0);
	row_sd=numeric(0);
	isnum=numeric(0);
	for(r in 1:nrow(mat)){
		isnum=mat[r,];
		row_mean[r]=mean(isnum[!is.na(isnum)]);
		row_sd[r]=sd(isnum[!is.na(isnum)]);
	}
	par(mar=c(0,left_mar,2,.5)); #b,l,t,r
	image(1:nrow(mat), 1,matrix(row_mean,nrow=nrow(mat),ncol=1), zlim=c(matmin, matmax),xaxt="n", yaxt="n", xlab="", ylab="",c=colors);
	for(i in 1:nrow(mat)){
		#text(i,1,labels=sprintf("%.4g\n+/-%.4g",row_mean[i],row_sd[i]), cex=(1.7/max_width)*(cscale+rscale), srt=45, font=3);
		text(i,1,labels=sprintf("%.4g",row_mean[i]), cex=(1.7/max_width)*(cscale+rscale), srt=45, font=3);
	}
	mtext("mean   ", at=1, side=2, las=2, cex=rscale*.7, font=3);
	print("Row margins done.");

	# Title
	mtext(InputFileName, padj=-1, side=3, cex=rscale*.7);

	# Plot right col margin
	col_mean=numeric(0);
	col_sd=numeric(0);
	isnum=numeric(0);
	for(c in 1:ncol(mat)){
		isnum=mat[,c];
		col_mean[c]=mean(isnum[!is.na(isnum)]);
		col_sd[c]=sd(isnum[!is.na(isnum)]);
	}
	par(mar=c(bottom_mar,0,.5,2)); #b,l,t,r
	image(1,1:ncol(mat),matrix(col_mean,nrow=1,ncol=ncol(mat)), zlim=c(matmin, matmax), xaxt="n", yaxt="n", xlab="", ylab="", c=colors)
	for(i in 1:ncol(mat)){
		#text(1,i,labels=sprintf("%.4g\n%.4g",col_mean[i],col_sd[i]), cex=(1.7/max_width)*(cscale+rscale), srt=45, font=3);
		text(1,i,labels=sprintf("%.4g",col_mean[i]), cex=(1.7/max_width)*(cscale+rscale), srt=45, font=3);
	}
	mtext("mean   ", at=1, side=1, las=2, cex=cscale*.7, font=3);
	print("Col margins done.");

	
	# Plot Color Key
	scale=matrix(1:100,1,100);

	par(mar=c(bottom_mar+10,2,1+10,3)); #b,l,t,r
	image(1:nrow(scale),1:ncol(scale),scale,c=colors,
		xaxt="n", yaxt="n", xlab="", ylab="");
	cat("Min = ", matmin, " Max = ", matmax, "\n");
	num_labels=5;
	num_divisions=num_labels-1;
	incr=(matmax-matmin)/(num_divisions);
	label=character(0);
	pos=numeric(0);
	for(i in 1:num_labels){
		label[i]=sprintf("%1.1f ", ((i-1)*incr)+matmin); 
		pos[i]=((i-1)*100/(num_labels-1))+1;
		if(pos[i]>100){pos[i]=100};
		cat(i, pos[i], label[i], "\n", sep=" ");
	}
	mtext(label, at=pos, las=2, side=2, cex=.6);
	
	matmean=mean(nanfree_mat);
	meanpos=100*(matmean-matmin)/(matmax-matmin);
	matmed=median(nanfree_mat);
	medpos=100*(matmed-matmin)/(matmax-matmin);
	cat("Mean = ", matmean, " Median = ", matmed, "\n");
	maxpos=100*(matmax-matmin)/(matmax-matmin);
	mtext(c(" min"," mean"," median"," max"), at=c(1,meanpos,medpos,maxpos), las=2, side=4, cex=.5);
	print("Legend Done.\n");

	arg_count=arg_count+1;
}

dev.off();

writeLines("Done.\n")

q(status=0)
