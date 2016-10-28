#!/usr/bin/env Rscript

###############################################################################

progname = commandArgs(FALSE)[4]
args = commandArgs(TRUE)

arg_count=1

if(is.na(args[arg_count])){

	script_name <- unlist(strsplit(progname,"="))[2]
	usage <- paste (
		"\nUsage:\n\t", script_name, "\n\t\t<Input summary_table.tsv FileName>\n\n",
		"Generates a pie chart\n",
		"\n")

	writeLines(usage)
	writeLines("Input summary table FileName not defined.\n")
	quit(status=0)
}

###############################################################################

get_colors=function(num_col, alpha=1){
        colors=hsv(seq(0,1,length.out=num_col+1), c(1,.5), c(1,.75,.5), alpha=alpha);
        color_mat_dim=ceiling(sqrt(num_col));
        color_pad=rep("grey", color_mat_dim^2);
        color_pad[1:num_col]=colors[1:num_col];
        color_mat=matrix(color_pad, nrow=color_mat_dim, ncol=color_mat_dim);
        colors=as.vector(t(color_mat));
        colors=colors[colors!="grey"];
}

###############################################################################
# Main program loop

while(!(is.na(args[arg_count]))){
	InputFileName=args[arg_count]

	OutputFileRoot=InputFileName;
	OutputFileRoot=gsub("\\.summary_table\\.tsv$", "", OutputFileRoot);
	OutputFileRoot=gsub("\\.summary_table\\.xls$", "", OutputFileRoot);

	PiePDF= paste(OutputFileRoot, ".pie.pdf", sep="");
	PieXLS= paste(OutputFileRoot, ".pie.tsv", sep="");

	cat("\n");
	cat("Input File Name: ", InputFileName, "\n");
	cat(" Pie Output pdf: ", PiePDF, "\n");
	cat(" Pie Output tsv: ", PieXLS, "\n");
	cat("\n");

	###############################################################################
	###############################################################################

	# Load data
	mat=as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE, row.names=1, comment.char=""));
	#print(mat);

	count_mat=mat[,2:ncol(mat), drop=F];

	num_samples=nrow(count_mat);
	num_categories=ncol(count_mat);
	cat("Num Samples:    ", num_samples, "\n");
	cat("Num Categories: ", num_categories, "\n");
	categories=colnames(count_mat);

	short_names=character(num_categories);
	#for(i in 1:num_categories){
	#	taxonomy=unlist(strsplit(categories[i], " "));
	#	short_names[i]=taxonomy[length(taxonomy)];
	#}
	short_names=categories;

	###############################################################################

	pdf(PiePDF, width=20, height=8.5)

	# Normalized before combining
	samp_tot=apply(count_mat, 1, sum);
	norm_mat=count_mat;
	for(samp_ix in 1:num_samples){
		norm_mat[samp_ix,]=count_mat[samp_ix,]/samp_tot[samp_ix];
	}

	# Combine all samples
	all_samp_norm=apply(norm_mat, 2, mean);
	print(all_samp_norm);
	order_ix=order(all_samp_norm, decreasing=T);

	# Reorder values and categories names
	all_samp_norm=all_samp_norm[order_ix];
	short_names=short_names[order_ix];

	# Combine small abundances to reduce clutter
	lt1p=all_samp_norm[all_samp_norm >= 0 & all_samp_norm < 0.01];
	lt2p=all_samp_norm[all_samp_norm >= .01 & all_samp_norm < 0.02];
	lt3p=all_samp_norm[all_samp_norm >= .02 & all_samp_norm < 0.03];

	lt1p_sum=sum(lt1p);
	lt2p_sum=sum(lt2p);
	lt3p_sum=sum(lt3p);

	combined_val=all_samp_norm[all_samp_norm >=.03];
	combined_val=c(combined_val, lt3p_sum, lt2p_sum, lt1p_sum);

	# Sort by abundance
	combined_names=short_names[all_samp_norm >=.03];
	combined_names=c(combined_names, "2-3%", "1-2%", "<1%");
	num_slices=length(combined_names);

	# Append percentages to labels
	labels=character(num_slices);
	for(i in 1:num_slices){
		labels[i]=sprintf("%s (%1.1f%%)", combined_names[i], combined_val[i]*100);
		cat("Slice Label [", i, "] : ", labels[i], "\n");
	}
	
	colors=get_colors(num_slices);

	pie(combined_val, labels, col=colors, init.angle=180, clockwise=TRUE,
		main=OutputFileRoot);

	dev.off();

	###############################################################################
	# Output table

	fc=file(PieXLS, "w");
	outline=paste("Category Name", "Percentage", sep=",");
	write(outline, file=fc);

	for(i in 1:num_slices){
		outline=paste(combined_names[i], sprintf("%.2f", combined_val[i]*100), sep=",");
		write(outline, file=fc);
	}

	###############################################################################

	arg_count=arg_count+1;
}

writeLines("Done.\n")

q(status=0)
