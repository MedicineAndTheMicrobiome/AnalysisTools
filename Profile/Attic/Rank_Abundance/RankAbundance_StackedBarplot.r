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

REMOVE_UNKNOWN=FALSE;
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

	if(REMOVE_UNKNOWN){
		RAPlot= paste(OutputFileRoot, ".rank_abundancy.stacked_presence.no_weak_class.pdf", sep="")
	}else{
		RAPlot= paste(OutputFileRoot, ".rank_abundancy.stacked_presence.pdf", sep="")
	}

	PresenceTableCSV= paste(InputFileName, ".presence", sep="");

	cat("\n")
	cat("              Input File Name: ", InputFileName, "\n")
	cat(" Rank Abundancy Box Plots PDF: ", RAPlot, "\n")
	cat("           Presence Table CSV: ", PresenceTableCSV, "\n")

	###############################################################################
	###############################################################################

	# Load data
	inmat<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="", quote=NULL))
	#cat("Original Matrix:\n")
	#print(inmat);

	# Get input dimensions
	num_orig_samples=nrow(inmat);
	num_orig_categories=ncol(inmat)-2;
	cat("Num samples:", num_orig_samples, "\n");
	cat("Num categories:", num_orig_categories, "\n");

	# Get sample/category names
	orig_sample_names=inmat[,1];
	#print(sample_names);
	orig_category_names=colnames(inmat)[3:(num_orig_categories+2)];
	#print(category_names);

	# Only keep the taxon with the finest classification
	shortened_category_names=rep("",num_orig_categories);
	for(catidx in 1:num_orig_categories){
		#print(category_names[catidx]);
		components=unlist(strsplit(orig_category_names[catidx], " "));
	#	shortened_category_names[catidx]=tail(components,1);
		shortened_category_names[catidx]=orig_category_names[catidx];
	}
	#print(shortened_category_names);
	
	# Convert counts/strings into integers
	orig_counts=matrix(nrow=num_orig_samples, ncol=num_orig_categories);
	orig_counts[1:num_orig_samples,1:num_orig_categories]=as.numeric(inmat[1:num_orig_samples,3:(num_orig_categories+2)]);
	#print(orig_counts);

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
	#cat("Normalized:\n");
	#print(orig_normalized);

	# Remove any samples with 0 counts
	nonzero_samples=orig_totals>0;
	nonzero_samples_normalized=orig_normalized[nonzero_samples,];
	nonzero_samples_names=orig_sample_names[nonzero_samples];
	num_nonzero_samples=nrow(nonzero_samples_normalized);
	#print(nonzero_samples_normalized);
	#print(nonzero_samples_names);

	pdf(RAPlot, height=11, width=8.5);
	plots_generated=F;

	min_cutoffs=c(0, .01, .02, .03, .04, .05);
	for(min_cutoff in min_cutoffs){
		# Set up increment thresholds for presence/absence
		INCR=.1;
		thresholds=seq(0,1,INCR);
		thresholds[1]=min_cutoff;
		num_thresholds=length(thresholds);
		num_bins=num_thresholds-1;
		cat("Thresholds:\n");
		print(thresholds);

		labels=sprintf("<%.2f",thresholds[2:num_thresholds]);
		labels[1]=sprintf("%.2f-%.2f", min_cutoff, thresholds[2]);

		# Count up presence based on thresholds, save in matrix
		mat=numeric(0);
		for(t in 1:num_bins){
			gt=nonzero_samples_normalized>thresholds[t] & nonzero_samples_normalized<=thresholds[t+1];
			cat(thresholds[t], " - ", (thresholds[t+1]), "\n");
			thres_counts=numeric(0);
			for(cat in 1:num_orig_categories){
				thres_counts[cat]=sum(gt[,cat]);
			}
			mat=rbind(mat,thres_counts);
		}
		if(length(dim(mat))!=2){
			mat=t(mat);
		}

		# Sum up counts for each category
		totals=rep(0,num_orig_categories);
		for(cat in 1:num_orig_categories){
			totals[cat]=sum(mat[,cat]);	
		}
		max_count=max(totals);

		# Sort by abundance
		sorted=sort(totals, decreasing=TRUE, index.return=TRUE);
		num_nonzero=sum(sorted$x>0);
		if(num_nonzero>0){
			
			# Compute colors based on number of bins	
			colors=rev(rainbow(num_bins, start=0, end=0.65));

			# Layout space for scale and plot
			vect=c(1,2,2,2,2,2,2,2,2,2,2,2,2);
			ly=matrix(vect,nrow=length(vect),ncol=1);
			layout(ly);

			# Compute percentage tick positions
			percs=c(0,5,10,25,33,50,66,75,90,95,100);
			perc_pos=percs*num_nonzero_samples/100;

			# Sort the matrix
			sorted_mat=mat[,sorted$ix];
			nonzero_mat=sorted_mat[,1:num_nonzero, drop=F];
			sorted_cat_names=shortened_category_names[sorted$ix];
			nonzero_cat_names=sorted_cat_names[1:num_nonzero];

			if(REMOVE_UNKNOWN){
				#print(nonzero_cat_names);

				remove_array=numeric(0);
				for(i in 1:length(nonzero_cat_names)){
					split_result=strsplit(nonzero_cat_names[i], ":");
					remove_array[i]=length(split_result[[1]])>1;
				}
				nonzero_mat=nonzero_mat[,!remove_array, drop=F];
				nonzero_cat_names=nonzero_cat_names[!remove_array];
				num_nonzero=ncol(nonzero_mat);
				if(!length(num_nonzero) || num_nonzero==0){
					next;	
				}
				#print(nonzero_cat_names);
			}

			#-------------------------------------------------------------------
			# Plot all
			# Plot scale
			cat("Plotting stacked bar plot for all categories at threhold: ", min_cutoff, ".\n");
			plots_generated=T;

			par(mar=c(0,20,3,20)); # bottom, left, top, right)
			scale_mat=matrix(c(1:num_bins), ncol=1, nrow=num_bins);
			image(1:nrow(scale_mat), 1:ncol(scale_mat), scale_mat,col=colors, 
				xaxt="n", yaxt="n", xlab="", ylab="", main=InputFileName);
			text(1:num_bins, y=rep(1,num_bins),labels=labels);

			# Estimate label/margin sizes
			label_scale=70/num_nonzero;
			if(label_scale>1){label_scale=1;}
			par(mar=c(15*label_scale, 4.5, .2, 4)); # bottom, left, top, right)

			# Build empty plot, so we can draw reference lines behind the actual plot
			barplot(rep(0,num_nonzero), yaxt="n", xaxt="n", ylim=c(0,num_nonzero_samples));
			abline(h=perc_pos, lty=2, col="grey");

			# Plot actual rank abundance graph
			barplot(nonzero_mat, xlab="", ylab="Absolute Sample Presence Count",
				names=nonzero_cat_names, las=2, cex.names=label_scale, cex.axis=1, 
				col=colors, ylim=c(0,num_nonzero_samples), border=NA,
				add=TRUE);

			# Draw percentage axes
			axis(4, at=perc_pos, label=percs, las=2);
			mtext("Percentage of Samples", side=4, line=2.5);
			mtext(sprintf("All, Min Cutoff = %.2f", min_cutoff), side=3, at=-.5, adj=0);
		
			#-------------------------------------------------------------------
			# Plot top 10

			top_to_plot=10;
			if(num_nonzero<top_to_plot){
				top_to_plot=num_nonzero;
			}

			cat("Top to plot: ", top_to_plot, "\n");
			if(top_to_plot>1){

				# Plot scale
				par(mar=c(0,20,3,20)); # bottom, left, top, right)
				scale_mat=matrix(c(1:num_bins), ncol=1, nrow=num_bins);
				image(1:nrow(scale_mat), 1:ncol(scale_mat), scale_mat,col=colors, xaxt="n", yaxt="n", 
					xlab="", ylab="", main=InputFileName);
				labels=sprintf("<%.2f",thresholds[2:num_thresholds]);
				labels[1]=sprintf("%.2f-%.2f", min_cutoff, thresholds[2]);
				text(1:num_bins, y=rep(1,num_bins),labels=labels);

				# Estimate label/margin sizes
				label_scale=70/top_to_plot;
				if(label_scale>1.5){label_scale=1.5;}
				par(mar=c(10*label_scale,4.5,.2,4)); # bottom, left, top, right)

				cat("Plotting stacked bar plot for top: ", top_to_plot, " categories.\n");

				# Build empty plot, so we can draw reference lines behind the actual plot
				barplot(rep(0,top_to_plot), yaxt="n", xaxt="n", ylim=c(0,num_nonzero_samples));
				abline(h=perc_pos, lty=2, col="grey");

			
				# Plot actual rank abundance graph
				barplot(nonzero_mat[,1:top_to_plot, drop=F], xlab="", ylab="Absolute Sample Presence Count",
					names=nonzero_cat_names[1:top_to_plot], las=2, cex.names=label_scale, cex.axis=1, 
					col=colors, ylim=c(0,num_nonzero_samples), 
					add=TRUE);

				# Draw percentage axes
				axis(4, at=perc_pos, label=percs, las=2);
				mtext("Percentage of Samples", side=4, line=2.5);
				mtext(sprintf("Top %i, Min Cutoff = %.2f", top_to_plot, min_cutoff), side=3, at=-.5, adj=0);
			}else{
				cat("Not plotting for top: ", top_to_plot, "\n");
			}
		}

		if(!plots_generated){
			par(mar=c(1,1,1,1));
			plot(0,0, type="n");
			text(0,0, "No plots generated because no abundances exceeded cutoffs.\n");
		}

		###########################################################################
		# Output table

		fc=file(sprintf("%s.%02i.csv", PresenceTableCSV, min_cutoff*100), "w");

		outline=paste("Category", "Percent_Presence",  "Total_Count", paste(labels, collapse=","), sep=",");
		write(outline, file=fc);
		for(i in 1:length(sorted_cat_names)){
			perc=100*sorted$x[i]/num_nonzero_samples;
			counts_str=paste(sorted_mat[,i], collapse=",");
			org_total=sum(sorted_mat[,i]);
			outline=paste(sorted_cat_names[i], sprintf("%.2f", perc), org_total, counts_str, sep=",");
			write(outline, file=fc);
		}


		close(fc);


	}

	dev.off();

	###############################################################################

	arg_count=arg_count+1;
}

writeLines("Done.\n")

q(status=0)
