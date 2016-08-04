#!/usr/local/bin/Rscript

###############################################################################

progname <- commandArgs(FALSE)[4]
args <- commandArgs(TRUE)

arg_count=1

if(is.na(args[arg_count])){

	script_name <- unlist(strsplit(progname,"="))[2]
	usage <- paste (
		"\nUsage:\n\t", script_name, "\n\t\t<.stat file from computing differences>\n\n",
		"This script generates a graph similar to rank abundancy, except instead of abundancy\n",
		"the difference between two samples is generated.  The input is the .stat from from\n",
		"computing the venn between two samples.",
		"\n")

	writeLines(usage)
	writeLines("Input FileName not defined.\n")
	quit(status=0)
}

###############################################################################
# Main program loop

while(!(is.na(args[arg_count]))){
	InputFileName=args[arg_count]
	RDPlot= paste(InputFileName, ".rank_difference.pdf", sep="")

	cat("\n")
	cat("           Input File Name: ", InputFileName, "\n")
	cat(" Rank Difference Plots PDF: ", RDPlot, "\n")

	###############################################################################
	###############################################################################

	# Load data
	A<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="*"))
	#cat("Original Matrix:\n")
	#print(A)

	pdf(RDPlot, width=11, height=8.5);

	for(sorttype in 1:4){
	
		# Grab columns we need into a vector
		categories=as.character(A[,1]);
		num_cat=length(categories);
		cat("Num categories: ", num_cat, "\n");

		diff=as.double(A[,7]);
		zstat=as.double(A[,8]);
		pval=as.double(A[,9]);
		pval[is.na(pval)]=1;
	
		a_counts=as.integer(A[,2]);
		b_counts=as.integer(A[,4]);

		a_counts[is.na(a_counts)]=0;
		b_counts[is.na(b_counts)]=0;

		a_prob=as.double(A[,3]);
		b_prob=as.double(A[,5]);
		ab_ratio=as.double(A[,6]);
		ab_ratio[a_prob==0 & b_prob==0]=0;


		cat("Working on ");
		if(sorttype==1){
			comptype="A-B";
			sorted=sort(diff, index=TRUE,decreasing=F);
			comparator=diff[sorted$ix];
		}else if(sorttype==2){
			comptype="A/B";
			ab_ratio[a_prob==0]=0;
			sorted=sort(ab_ratio, index=TRUE,decreasing=T);
			comparator=ab_ratio[sorted$ix];
		}else if(sorttype==3){
			comptype="B/A";
			ba_ratio=1/ab_ratio;
			ba_ratio[b_prob==0]=0;
			sorted=sort(ba_ratio, index=TRUE,decreasing=T);
			comparator=ba_ratio[sorted$ix];
		}else if(sorttype==4){
			comptype="1 - p-value";
			comparator[is.na(comparator)]=1;
			sorted=sort(pval, index=TRUE,decreasing=F);
			comparator=1-pval[sorted$ix];
		}
		cat(comptype,"\n");
		#print(comparator);

		categories=categories[sorted$ix];
#		cat("Categories" , categories,"\n");
		zstat=zstat[sorted$ix];
		pval=pval[sorted$ix];
		a_prob=a_prob[sorted$ix];
		b_prob=b_prob[sorted$ix];
		a_counts=a_counts[sorted$ix];
		b_counts=b_counts[sorted$ix];
		a_total=sum(a_counts);
		b_total=sum(b_counts);

#		cat("Categories" , categories,"\n");
		#print(categories);
		#print(zstat);
		#print(pval);

		# Find max difference so we know how much we need to plot
		range=c(min(comparator),max(comparator));
		max_not_inf=range[2];
		if(is.infinite(range[2])){
			max_not_inf=max(comparator[is.finite(comparator)]);
			range[2]=max_not_inf;
			infinity_height=1.2;
		}else{
			infinity_height=1;
		}
		padding_height=1.1
		range[2]=range[2]*infinity_height*padding_height;
		if(range[1]>0){
			range[1]=0;
		}

		# Shorten the name of the category
		shortnames=rep("",num_cat);
		for(i in 1:num_cat){
			splitnames=strsplit(as.character(categories[i])," ");
			shortnames[i]=splitnames[[1]][length(splitnames[[1]])];
		}

		#print(shortnames);
	#	print(categories);		
		# Append differences symbol
		#for(i in 1:num_cat){
		#	if(diff[i]>0){
		#		shortnames[i]=paste(shortnames[i], "+", sep=" ");
		#	}else if(diff[i]<0){
		#		shortnames[i]=paste(shortnames[i], "-", sep=" ");
		#	}else{
		#		shortnames[i]=paste(shortnames[i], "=", sep=" ");
		#	}
		#}

		# Compute colors based on p-value
		colors=rep("",num_cat);
		for(i in 1:num_cat){
			#print(pval[i]);
			if(is.na(pval[i])){
				colors[i]="darkgreen";
			}else{
				if(pval[i]<0.001){
					colors[i]="darkred";	
				}else if (pval[i]<0.01){
					colors[i]="darkorange";
				}else if (pval[i]<0.05){
					colors[i]="gold";
				}else{
					colors[i]="darkgreen";
				}
			}
		}

		
		# Generate plot
		resize=120/length(categories);

		layout(matrix(c(1,2,3,3,3,3,3,3,3)));
		axis_incr=.2;
		axis_ticks=seq(0,1,axis_incr);

		#(bottom, left, top, right);
		mar_right=4
		par(oma=c(0,1,4,4));

		par(mar=c(1,4,0,mar_right));
		pos=barplot(a_prob, ylim=c(0,1), las=2, ylab="Sample A");
		text(pos,a_prob,sprintf("%.3f",a_prob), srt=90, pos=3, offset=1.5);
		mtext(InputFileName, side=3, line=1, col="black");

		a_label=as.integer(axis_ticks*a_total);
		axis(side=4, at=axis_ticks, a_label,las=2);
		#print(a_prob);

		par(mar=c(1,4,0,mar_right));
		pos=barplot(b_prob, ylim=c(0,1), las=2, ylab="Sample B");
		text(pos,b_prob,sprintf("%3.3f",b_prob), srt=90, pos=3, offset=1.5);
		b_label=as.integer(axis_ticks*b_total);
		axis(side=4, at=axis_ticks, b_label,las=2);
		#print(b_prob);

		par(mar=c(15,4,0,mar_right));

		infinite_elements=is.infinite(comparator);
		comparator[infinite_elements]=max_not_inf*infinity_height;
		labels=sprintf("%3.3f",comparator);
		labels[infinite_elements]="Inf";

print(categories);
#print(range);
		pos=barplot(comparator, ylim=range, col=colors, border=NA, las=2, ylab=comptype,
			names.arg=shortnames,  cex.names=resize, las=2);
#		pos=barplot(comparator, ylim=range, col=colors, border=NA, las=2, ylab=comptype,
#			names.arg=sprintf("%s",categories),  cex.names=resize, las=2);
		text(pos[comparator>=0],comparator[comparator>=0],labels[comparator>=0], srt=90, pos=3, offset=1.5);
		if(length(pos[comparator<0])>0){
			text(pos[comparator<0],comparator[comparator<0],labels[comparator<0], srt=90, pos=1, offset=1.7);
		}
		#print(comparator);

		#mtext("Red: p < 0.001", side=1, line=-6, col="darkred");
		#mtext("Orange: p < 0.01", side=1, line=-5, col="darkorange");
		#mtext("Yellow: p < 0.05", side=1, line=-4, col="gold");
		#mtext("Green: p >= 0.05", side=1, line=-3, col="darkgreen");


		#barplot(rep(0,length(shortnames)), names.arg=shortnames,  cex.names=resize, las=2, mar=c(oma_b,0,oma_b,0));

	}

	dev.off();

	###############################################################################

	arg_count=arg_count+1;
}

writeLines("Done.\n")
print(warnings());

q(status=0)
