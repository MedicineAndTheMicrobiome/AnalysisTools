#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
        "input_file", "i", 1, "character",
        "output_file", "o", 2, "character",
	"shorten_name_char", "s", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
        "\nUsage:\n\n", script_name,
        "\n",
        "	-i <input summary table.xls>\n",
	"	[-o <output file name root>]\n",
	"	[-s \"<name splitter character, e.g. ;>\"]\n",
        "\n",
	"This script will read in a summary table and generate\n",
	"rank abundance plots, diversity indices, relative abundance, and\n",
	"descriptive stats.\n",
	"\n",
	"The -s shorten names flag will split the input category names\n",
	"by spaces, and chose the one furthest to the right.\n",
        "\n",
        "\n", sep="");

if(!length(opt$input_file)){
        cat(usage);
        q(status=-1);
}

InputFileName=opt$input_file;

if(!length(opt$output_file)){
        OutputFileName=gsub("\\.summary_table\\.xls$", "", opt$input_file);
        OutputFileName=gsub("\\.summary_table\\.tsv$", "", OutputFileName);
}else{
        OutputFileName=opt$output_file;
}

ShortenNameChar="";
if(length(opt$shorten_name_char)){
	ShortenNameChar=opt$shorten_name_char;
	cat("Shortening names by splitting by: ", ShortenNameChar, "\n");
	OutputFileName=paste(OutputFileName, ".shrt", sep="");
}

###############################################################################

tail_stat=function(x){
	# Takes pmf or categorical data
	sorted=sort(x, decreasing=TRUE);
	norm=sorted/sum(x);
	n=length(norm);
	tail=0;
	for(i in 1:n){
		tail=tail + norm[i]*((i-1)^2);
	}
	return(sqrt(tail));
}

normalize=function(x){
	sums=apply(x, 1, sum);
	norm=matrix(0, ncol=ncol(x), nrow=nrow(x));
	colnames(norm)=colnames(x);
	rownames(norm)=rownames(x);
	for(i in 1:nrow(x)){
		norm[i,]=counts[i,]/sums[i];
	}	
	return(norm);
}

###############################################################################

# Rank abundance plot
RA_All_Plot= paste(OutputFileName, ".rnk_abn.all.pdf", sep="")
RA_Top40_Plot= paste(OutputFileName, ".rnk_abn.top40.pdf", sep="")
RA_AllNr_Plot= paste(OutputFileName, ".rnk_abn.all.nr.pdf", sep="")
RA_Top100Nr_Plot= paste(OutputFileName, ".rnk_abn.top100.nr.pdf", sep="")
RA_Top75Nr_Plot= paste(OutputFileName, ".rnk_abn.top75.nr.pdf", sep="")
RA_Top40Nr_Plot= paste(OutputFileName, ".rnk_abn.top40.nr.pdf", sep="")
RA_Top25Nr_Plot= paste(OutputFileName, ".rnk_abn.top25.nr.pdf", sep="")
RA_Top15Nr_Plot= paste(OutputFileName, ".rnk_abn.top15.nr.pdf", sep="")
RA_Top10Nr_Plot= paste(OutputFileName, ".rnk_abn.top10.nr.pdf", sep="")

# Indices
DIFile= paste(OutputFileName, ".indices.csv", sep="")
CIFile= paste(OutputFileName, ".confidence.csv", sep="");
CI1LFile= paste(OutputFileName, ".confidence.1line.csv", sep="");
IDPlot= paste(OutputFileName, ".index_distributions.pdf", sep="")

cat("\n");
cat("             Input File Name: ", InputFileName, "\n");
cat("           Diversity Indices: ", DIFile, "\n");
cat("        Confidence Intervals: ", CIFile, "\n");
cat("Confidence Intervals (1 row): ", CI1LFile, "\n");
cat("         Index Distributions: ", IDPlot, "\n");
cat("\n");
cat("Shortening Category Char: ", ShortenNameChar, "\n");
cat("\n");

###############################################################################
###############################################################################

# Load data
inmat=as.matrix(read.delim(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="", 
	quote="", row.names=1))

#cat("Original Matrix:\n")
#print(inmat)

counts=inmat[,2:ncol(inmat)];

# Remove zero count samples
tot=apply(counts, 1, sum);
non_zero_samp=tot>0;
if(!all(non_zero_samp)){
	samp_names=rownames(counts);
	removed_samp_names=samp_names[!non_zero_samp];
	cat("WARNING: Zero Counts Samples Found:\n");
	print(removed_samp_names);
	cat("\n");	
	counts=counts[non_zero_samp,, drop=F];	
}


norm=normalize(counts);

num_samples=nrow(counts);
num_categories=ncol(counts);
cat("Num Samples: ", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");
cat("\n");

category_names=colnames(counts);
sample_names=rownames(counts);

# Shorten display name if requested
display_name=character(num_categories);
if(ShortenNameChar!=""){
	for(i in 1:length(category_names)){
		display_name[i]=utils::tail(unlist(strsplit(category_names[i], ShortenNameChar)),1);
		display_name[i]=gsub("_unclassified", "", display_name[i]);
	}
}else{
	display_name=category_names;
}
names(display_name)=category_names;

norm_mod_names=norm;
colnames(norm_mod_names)=display_name;


# List for diversity indices
num_distrib=num_samples + 2; # 1 for all reads, and 1 for equal weighted
Entropy=rep(0, num_distrib);
Evenness=rep(0, num_distrib);
Simpsons=rep(0, num_distrib);
SimpsonsRecip=rep(0, num_distrib);
Tail=rep(0, num_distrib);
DiscTaxa=rep(0, num_distrib);

# Per sample category/distribution manipulations
samp_norm_list=list();
samp_norm_sorted_list=list();
samp_count_list=list();

# Place all distributions we want to analyze in sample list
ALL_TOTCTS_NAME="< All Total Counts >";
ALL_NORM_NAME="< All Normalized >";
totcts_idx=num_samples+1;
allnorm_idx=num_samples+2;
non_sample_idx=c(totcts_idx, allnorm_idx);
sample_names=c(sample_names, ALL_TOTCTS_NAME, ALL_NORM_NAME);

norm_vect=function(x){
	return(x/sum(x));
}

samp_norm_list[[ALL_TOTCTS_NAME]]=norm_vect(apply(counts, 2, sum));
samp_norm_list[[ALL_NORM_NAME]]  =norm_vect(apply(norm, 2, sum));

# Convert count matrix into list of vectors
for(i in 1:num_samples){
	cur_name=sample_names[i];
	samp_norm_list[[cur_name]]=norm[cur_name,];
}

# Compute stats on each of the distributions
for(i in 1:num_distrib){

	cur_name=sample_names[i]
	cat("Working on: ", cur_name, "\n");

	cur_norm=samp_norm_list[[cur_name]];
	nz_norm=cur_norm[cur_norm>0]

	Entropy[i]=-sum(nz_norm*log(nz_norm));
	Evenness[i]=Entropy[i]/log(length(nz_norm));
	Simpsons[i]=1-sum(cur_norm^2);
	SimpsonsRecip[i]=1/sum(cur_norm^2);
	Tail[i]=tail_stat(nz_norm);
	DiscTaxa[i]=length(nz_norm);
	samp_norm_sorted_list[[cur_name]]=sort(cur_norm, decreasing=TRUE, method="shell");

	if(1){
		cat("Sample Name = ", sample_names[i], "\n");
		cat("Shannon's Diversity = ", Entropy[i], "\n");
		cat("Simpson's = ",Simpsons[i], "\n");
		cat("Evenness = ", Evenness[i], "\n");
		cat("Simpson's Recip= ",SimpsonsRecip[i], "\n");
		cat("Tail = ", Tail[i], "\n");
		cat("DiscTaxa = ", DiscTaxa[i], "\n");

		cat("\n");
	}
}

div_indices=list();
div_indices[["Entropy"]]=Entropy;
div_indices[["Evenness"]]=Evenness;
div_indices[["Simpsons"]]=Simpsons;
div_indices[["SimpsonsRecip"]]=SimpsonsRecip;
div_indices[["Tail"]]=Tail;
div_indices[["DiscTaxa"]]=DiscTaxa;

# Assign color increment according to CDF
all_norm_sorted_dist=samp_norm_sorted_list[[sample_names[allnorm_idx]]];
cum_sum_prob=c(0, cumsum(all_norm_sorted_dist));
cum_sum_prob=cum_sum_prob[1:(length(cum_sum_prob)-1)];
names(cum_sum_prob)=names(all_norm_sorted_dist);
allcolors=rainbow(100, start=0, end=.65);
indices=as.integer(cum_sum_prob*100)+1;
indices=sapply(indices, function(x){min(100, x);});
color_map=allcolors[indices];
names(color_map)=names(all_norm_sorted_dist);
#print(color_map);

###############################################################################

plot_rank_abundance=function( 
		output_filename, top_to_plot=-1,
		normalize_sorted_abundances, indices,
		color_map, display_name, ignore_remaining=F){

	cat("\n");
	cat("Generating Rank Abundance Plot in: ", output_filename, "\n");
	pdf(output_filename, width=11, height=9.5);

	num_categories=length(display_name);
	sample_names=names(normalize_sorted_abundances);
	num_samples=length(sample_names);

	if(top_to_plot==-1){
		top_to_plot=num_categories;
	}else{
		top_to_plot=min(c(top_to_plot, num_categories));
	}

	resize=30/top_to_plot;
	if(resize > 1.0){
		resize=1.0;
	}
	
	par(mar=c(25*resize+1,4,4,15));

	cat("Num categories to plot: ", top_to_plot, "\n");
	cat("Num samples to plot: ", num_samples, "\n"); 

	sorted_by_sample_name=sort(sample_names, method="shell", index.return=TRUE);

	for(i in sorted_by_sample_name$ix){

		cur_samp_name=sample_names[i];
		cur_values=normalize_sorted_abundances[[cur_samp_name]];

		if(ignore_remaining){
			# Remove "Unknown" and "Remaining"
			categories=names(cur_values);
			uc_cat=toupper(categories);
			ignore_ix=which(uc_cat=="UNKNOWN" | uc_cat=="REMAINING" | uc_cat=="REMAINDER");
			if(length(ignore_ix)>0){
				cur_values=cur_values[-ignore_ix];
			}
			top_to_plot=min(length(cur_values), top_to_plot);
		}

		category_names=names(cur_values);
		cur_disp_names=display_name[category_names];

		cat("   Plotting: ", cur_samp_name, "\n");

		# Generate labels that include abundance and category name
		category_str=character();
		abund_str=character();
		for(j in 1:length(cur_values)){
			if(cur_values[j]>0){
				abund_str[j]=sprintf("%3.3f",cur_values[j]);
			}else{
				abund_str[j]="----";
			}
			
			category_str[j]=paste("             ", cur_disp_names[j], sep="");
		}

		#cat("Categories: \n");
		#print(category_str);
		#cat("Abundances: \n");
		#print(abund_str);

		# plot bar plot without any labeling
		assigned_color=color_map[category_names];
		ylimit=max(cur_values[1:top_to_plot])*1.2;

		barmids=barplot(
			cur_values[1:top_to_plot],
			ylim=c(0,ylimit), 
			names.arg="",
			main=sample_names[i], 
			cex.names=resize, las=2, col=assigned_color[1:top_to_plot]
		);

		# plot labels
		graph_bottom=0;
		spacing=(barmids[2]-barmids[1]);
		for(j in 1:top_to_plot){
		       text(barmids[j]-spacing*.65, graph_bottom-ylimit*.035, abund_str[j], 
				srt=-45, xpd=T, cex=resize, pos=4);

		       text(barmids[j]-spacing*.65, graph_bottom-ylimit*.035, category_str[j], 
				srt=-45, xpd=T, cex=resize, pos=4, font=2);
		}


		# Plot diversity annotation
		cumul_rep=sum(cur_values[1:top_to_plot]);
		mtext(paste("Top ", top_to_plot, " Categories", sep=""), side=3, line=0);
		mtext(paste("Representing ", sprintf("%3.2f",cumul_rep*100), "% of all data", sep=""), 
			side=3, line=-1);

		mtext(paste("Shannon:  ", sprintf("%4.4f", indices[["Entropy"]][i])), 
			side=3, line=-3, cex=.8);
		mtext(paste("Simpson:  ", sprintf("%4.4f", indices[["Simpsons"]][i])), 
			side=3, line=-4, cex=.8);
		mtext(paste("Evenness: ", sprintf("%4.4f", indices[["Evenness"]][i])), 
			side=3, line=-5, cex=.8);
		mtext(paste("Simpson's Recip: ", sprintf("%4.4f", indices[["SimpsonsRecip"]][i])), 
			side=3, line=-6, cex=.8);
		mtext(paste("Tail: ", sprintf("%4.4f", indices[["Tail"]][i])), 
			side=3, line=-7, cex=.8);
		mtext(paste("DiscTaxa: ", sprintf("%i", indices[["DiscTaxa"]][i])), 
			side=3, line=-8, cex=.8);

	}

	dev.off();

}

# Plot remaining/remainder/unknown
plot_rank_abundance(RA_All_Plot, -1, samp_norm_sorted_list, div_indices, color_map, display_name, F);
plot_rank_abundance(RA_Top40_Plot, 40, samp_norm_sorted_list, div_indices, color_map, display_name, F);

# Remove remaining/remainder/unknown
plot_rank_abundance(RA_AllNr_Plot, -1, samp_norm_sorted_list, div_indices, 
	color_map, display_name, T);
plot_rank_abundance(RA_Top100Nr_Plot, 100, samp_norm_sorted_list, div_indices, 
	color_map, display_name, T);
plot_rank_abundance(RA_Top75Nr_Plot, 75, samp_norm_sorted_list, div_indices, 
	color_map, display_name, T);
plot_rank_abundance(RA_Top40Nr_Plot, 40, samp_norm_sorted_list, div_indices, 
	color_map, display_name, T);
plot_rank_abundance(RA_Top25Nr_Plot, 25, samp_norm_sorted_list, div_indices, 
	color_map, display_name, T);
plot_rank_abundance(RA_Top15Nr_Plot, 15, samp_norm_sorted_list, div_indices, 
	color_map, display_name, T);
plot_rank_abundance(RA_Top10Nr_Plot, 10, samp_norm_sorted_list, div_indices, 
	color_map, display_name, T);

###############################################################################

# Output indices

fc=file(DIFile, "w")

outline=paste("Sample ID", "Shannon", "Simpson", "Evenness", 
	"SimpsonsRecip", "Tail", "DiscTaxa", sep=",");

write(outline, file=fc);

sorted_by_entropy=sort(Entropy, decreasing=TRUE, method="shell", index.return=TRUE);

for(i in sorted_by_entropy$ix){
	outline=paste(sample_names[i], 
		Entropy[i], Simpsons[i], Evenness[i], SimpsonsRecip[i], Tail[i], DiscTaxa[i], sep=",");
	write(outline, file=fc);
}

close(fc);

###############################################################################

# Output Confidence Intervals and Standard Deviations
ci = function(x, alpha){
	n=length(x);

	med=median(x);

	ba=1-(n-2)/n;
	#cat("Best alpha = ", ba, "\n");
	if(ba <= (alpha+.0000000000000001)){
		sorted=sort(x);
		lb=sorted[floor(n*(alpha/2))+1];
		ub=sorted[ceiling(n*(1-(alpha/2)))];
		return(c(med,lb,ub));
	}else{
		return(c(med,NA,NA))
	}
}

# Compute statistics, remove the < ALL > entry first
#print(SampleNames[-(NumSamples+1)]);
entropy_samponly=Entropy[-non_sample_idx];
simpsons_samponly=Simpsons[-non_sample_idx];
evenness_samponly=Evenness[-non_sample_idx];
simpsonsrecip_samponly=SimpsonsRecip[-non_sample_idx];
tail_samponly=Tail[-non_sample_idx];
disctaxa_samponly=DiscTaxa[-non_sample_idx];

entropy_ci=ci(entropy_samponly, .05);
simpsons_ci=ci(simpsons_samponly, .05);
evenness_ci=ci(evenness_samponly, .05);
simpsonsrecip_ci=ci(simpsonsrecip_samponly, .05);
tail_ci=ci(tail_samponly, .05);
disctaxa_ci=ci(disctaxa_samponly, .05);

# Plot distributions
pdf(IDPlot, width=11, height=8.5);
par(mfrow=c(2,3), oma=c(0,0,1,0));
hist(entropy_samponly, main="Entropy");
hist(simpsons_samponly, main="Simpson's");
title=utils::tail(strsplit(OutputFileName,"/")[[1]], n=1);
mtext(title, side=3, line=3, cex=1.3, font=2);
hist(evenness_samponly, main="Evenness");
hist(simpsonsrecip_samponly, main="Simpson's Reciprocal");
hist(tail_samponly, main="Tail");
hist(disctaxa_samponly, main="Discovered Taxa");

# Write out statistics
fc=file(CIFile,"w");
write(paste("Statistic","Median","LB_95","UB_95", sep=","), file=fc);
write(paste("Entropy", paste(entropy_ci, collapse=","), sep=","), file=fc);
write(paste("Simpsons", paste(simpsons_ci, collapse=","), sep=","), file=fc);
write(paste("Evenness", paste(evenness_ci, collapse=","), sep=","), file=fc);
write(paste("SimpsonsRecip", paste(simpsonsrecip_ci, collapse=","), sep=","), file=fc);
write(paste("Tail", paste(tail_ci, collapse=","), sep=","), file=fc);
write(paste("DiscTaxa", paste(disctaxa_ci, collapse=","), sep=","), file=fc);
write(paste("N",num_samples, sep=","), file=fc);
close(fc);

fc=file(CI1LFile,"w");
file_path_comp=strsplit(InputFileName, "/")[[1]];
filename=file_path_comp[length(file_path_comp)];
cat(filename, file=fc);
cat(paste(",Entropy", paste(entropy_ci, collapse=","), sep=","), file=fc);
cat(paste(",Simpsons", paste(simpsons_ci, collapse=","), sep=","), file=fc);
cat(paste(",Evenness", paste(evenness_ci, collapse=","), sep=","), file=fc);
cat(paste(",SimpsonsRecip", paste(simpsonsrecip_ci, collapse=","), sep=","), file=fc);
cat(paste(",Tail", paste(tail_ci, collapse=","), sep=","), file=fc);
cat(paste(",DiscTaxa", paste(disctaxa_ci, collapse=","), sep=","), file=fc);
cat(paste(",N",num_samples, sep=","), file=fc);
cat("\n", file=fc);
close(fc);

###############################################################################
###############################################################################

mean_abund=apply(norm_mod_names, 2, mean);
order_ix=order(mean_abund, decreasing=T);

mean_abund=mean_abund[order_ix];
norm_mod_names=norm_mod_names[,order_ix];

outmatrix=matrix(NA, nrow=ncol(norm_mod_names), ncol=8);
rownames(outmatrix)=colnames(norm_mod_names);
colnames(outmatrix)=c("Mean", "StdDev", "StdErr", "Median", "CI95_LB", "CI95_UB", "Min", "Max");
N=nrow(norm_mod_names);

outmatrix[,"Mean"]=mean_abund;
outmatrix[,"StdDev"]=apply(norm_mod_names, 2, sd);
outmatrix[,"StdErr"]=outmatrix[,"StdDev"]/sqrt(N);
outmatrix[,"Median"]=apply(norm_mod_names, 2, median);
outmatrix[,"CI95_UB"]=apply(norm_mod_names, 2, function(x){quantile(x, .975)});
outmatrix[,"CI95_LB"]=apply(norm_mod_names, 2, function(x){quantile(x, .025)});
outmatrix[,"Min"]=apply(norm_mod_names, 2, min);
outmatrix[,"Max"]=apply(norm_mod_names, 2, max);

print(outmatrix);

pdf(paste(OutputFileName, ".combined_barplots.pdf", sep=""), height=8.5, width=11);

barplot_categories=function(center_values, top=10, lb, ub, data=NULL, title){

	top_center_val=center_values[1:top];

	top_lb=lb[1:top];
	top_ub=ub[1:top];

	if(!is.null(data)){
		datamax=max(data);
	}else{
		datamax=0;
	}

	par(mar=c(10,5,5,10));
	mids=barplot(top_center_val,
		names.arg="",
		ylim=c(0, max(ub, datamax)*1.1),
		main=paste("Top ", top, ": ", title, sep=""),
		xlab="", ylab="Proportion",
		fill="grey75"
	);
	bar_width=mids[2]-mids[1];
	ebw=bar_width/4;

	# Draw scatters
	if(!is.null(data)){
		num_samp=nrow(data);
		x_scatter=rnorm(num_samp, 0, ebw/4);
		for(i in 1:top){
			points(mids[i]+x_scatter, data[,i], col="grey50", cex=.5);
		}
	}

	# Draw bounds
	for(i in 1:top){
		points(c(mids[i]-ebw, mids[i]+ebw), c(top_ub[i], top_ub[i]), col="red", lwd=2, type="l"); 
		points(c(mids[i]-ebw, mids[i]+ebw), c(top_lb[i], top_lb[i]), col="blue", lwd=2, type="l"); 
		points(c(mids[i], mids[i]), c(top_lb[i], top_ub[i]), col="grey5", lwd=.5, type="l"); 
	}

	# Label 
        plot_range=par()$usr;
        plot_height=plot_range[4];
        label_size=min(c(1,.7*bar_width/par()$cxy[1]));
	cat_names=names(top_center_val);
        text(mids-par()$cxy[1]/2, rep(-par()$cxy[2]/2, top), cat_names, srt=-45, xpd=T, pos=4, cex=label_size);

}

for(num_disp_cat in c(10, 15, 20, 25, 30)){

	barplot_categories(outmatrix[,"Mean"], top=num_disp_cat, 
		lb=outmatrix[,"Mean"]-outmatrix[,"StdDev"],
		ub=outmatrix[,"Mean"]+outmatrix[,"StdDev"],
		title="Means and Standard Deviations");

	barplot_categories(outmatrix[,"Mean"], top=num_disp_cat, 
		lb=outmatrix[,"Mean"]-outmatrix[,"StdDev"],
		ub=outmatrix[,"Mean"]+outmatrix[,"StdDev"],
		data=norm_mod_names,
		title="Means and Standard Deviations");



	barplot_categories(outmatrix[,"Mean"], top=num_disp_cat, 
		lb=outmatrix[,"Mean"]-outmatrix[,"StdErr"],
		ub=outmatrix[,"Mean"]+outmatrix[,"StdErr"],
		title="Means and Standard Errors");

	barplot_categories(outmatrix[,"Mean"], top=num_disp_cat, 
		lb=outmatrix[,"Mean"]-outmatrix[,"StdErr"],
		ub=outmatrix[,"Mean"]+outmatrix[,"StdErr"],
		data=norm_mod_names,
		title="Means and Standard Errors");



	barplot_categories(outmatrix[,"Median"], top=num_disp_cat, 
		lb=outmatrix[,"CI95_LB"],
		ub=outmatrix[,"CI95_UB"],
		title="Medians and 95% CI");

	barplot_categories(outmatrix[,"Median"], top=num_disp_cat, 
		lb=outmatrix[,"CI95_LB"],
		ub=outmatrix[,"CI95_UB"],
		data=norm_mod_names,
		title="Medians and 95% CI");



	barplot_categories(outmatrix[,"Median"], top=num_disp_cat, 
		lb=outmatrix[,"Min"],
		ub=outmatrix[,"Max"],
		title="Medians and Min/Max");

	barplot_categories(outmatrix[,"Median"], top=num_disp_cat, 
		lb=outmatrix[,"Min"],
		ub=outmatrix[,"Max"],
		data=norm_mod_names,
		title="Medians and Min/Max");

}

dev.off();

# Write stats to file
fname=paste(OutputFileName, ".combined_statistcs.tsv", sep="");
cat(file=fname, "Category\t");
write.table(outmatrix, file=fname, append=T, quote=F, sep="\t");

###############################################################################

cat("Done.\n")
warnings();

q(status=0)
