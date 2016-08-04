#!/usr/local/bin/Rscript

###############################################################################

progname <- commandArgs(FALSE)[4]
args <- commandArgs(TRUE)

arg_count=1

if(is.na(args[arg_count])){

	script_name <- unlist(strsplit(progname,"="))[2]
	usage <- paste (
		"\nUsage:\n\t", script_name, "\n\t\t<Input summary_table.xls FileName>\n\n",
		"\n",
		" Generates a stacked barplot, where each bar is a sample and the colored components\n",
		" within each barplot are the categories.\n\n"
	);

	writeLines(usage)
	quit(status=0)
}

###############################################################################
# Main program loop

InputFileName=args[arg_count]

PlotPDF= paste(InputFileName, ".dom_prof.pdf", sep="")
cat("\n")
cat("                        Input File Name: ", InputFileName, "\n")
cat(" Dominance Profile PDF Output File Name: ", PlotPDF, "\n")

###############################################################################
###############################################################################

# Load data
inmat<-as.matrix(read.table(InputFileName, sep="\t", header=TRUE, check.names=FALSE))
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
	shortened_category_names[catidx]=tail(components,1);
}
#print(shortened_category_names);

# Convert counts/strings into integers
orig_counts=matrix(nrow=num_orig_samples, ncol=num_orig_categories);
orig_counts[1:num_orig_samples,1:num_orig_categories]=as.numeric(inmat[1:num_orig_samples,3:(num_orig_categories+2)]);
#cat("Original counts\n");
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

mean_dominance=apply(orig_normalized,2,mean);
#cat("Mean dominance:\n");
print(mean_dominance);
mean_sort_dec=sort(mean_dominance,decreasing=TRUE,index.return=TRUE);
mean_sort_inc=sort(mean_dominance,decreasing=FALSE,index.return=TRUE);

# Sort the categories by the dominance the "average" sample
normalized_dec=orig_normalized[,mean_sort_dec$ix];
category_names_dec=shortened_category_names[mean_sort_dec$ix];

#assign colors
mean_sort_wo_last_element=(mean_sort_dec$x[-num_orig_categories]);
renormal_mean=mean_sort_wo_last_element/sum(mean_sort_wo_last_element);
csum=cumsum(renormal_mean);
allcolors=(rainbow(10000, start=0, end=.65)); 		
color_idx=as.integer(csum*10000);
#print(color_idx);
colors_dec=c(allcolors[1],allcolors[color_idx]);

#colors_dec=sample(rainbow(num_orig_categories,start=0, end=1), num_orig_categories);

colors_inc=rev(colors_dec);


# Start plot
pdf(PlotPDF, height=8.5, width=11);

normalized_wAverage=rbind(orig_normalized, mean_dominance);
sample_names=c(orig_sample_names, "All");

# Sort so that most dominant are on top
normalized_wAverage_inc=matrix(0, nrow=nrow(normalized_wAverage), ncol=ncol(normalized_wAverage));
for(sample_idx in 1:(num_orig_samples+1)){
	normalized_wAverage_inc[sample_idx,]= normalized_wAverage[sample_idx,mean_sort_inc$ix];
}

#print(normalized_wAverage_inc);

###############################################################################
# Plot stack barplot
#(bottom, left, top, right);
par(oma=c(12,3,2,10));
par(mar=c(1,0,2,0));
resize=30/num_orig_samples;
if(resize>1){
	resize=1;
}

# Transpose for graphing
transposed=t(normalized_wAverage_inc);

barplot(transposed, names.arg=sample_names, axes = FALSE,
	las=2, cex.names=resize, col=colors_inc, main=InputFileName, ylim=c(0,1)
);

limit_cat=33;
display_cat=category_names_dec; 
display_col=colors_dec;
#print(display_cat);
if(num_orig_categories>=limit_cat){
	display_cat[limit_cat]="...";
	display_col[limit_cat]="black";
}
#print(display_cat);

#print(colors_dec);
mtext(display_cat[1:limit_cat], side=4, at=seq(.9,by=-.025,length.out=limit_cat), col=display_col[1:limit_cat], outer=TRUE, las=1);

perc=seq(0,1,.1);
axis(2, at=perc, label=sprintf("%i",as.integer(100*perc)), las=2);
 
###############################################################################


normalized_wAverage=t(normalized_wAverage);
#print(normalized_wAverage);

neg_log_norm=-log10(normalized_wAverage);
#print(neg_log_norm);

# Remove infinities
for(i in 1:nrow(neg_log_norm)){
	for(j in 1:ncol(neg_log_norm)){
		if(is.infinite(neg_log_norm[i,j])){
			neg_log_norm[i,j]=0;
		}
	}
}

#print(neg_log_norm);

# quantize/integerize 
ceilinged=ceiling(neg_log_norm);
#print(ceilinged);

max_lognorm=max(ceilinged);
cat("Maximum ceiling(-log(X)): ", max_lognorm, "\n");

num_samples_wAvg=ncol(neg_log_norm);
logcounts=matrix(0, nrow=max_lognorm, ncol=num_samples_wAvg);
colnames(logcounts)=sample_names;

for(samp_idx in 1:num_samples_wAvg){
	hist=as.data.frame(table(ceilinged[,samp_idx]));
	logval=as.integer(as.vector(hist$Var1));
	freq=as.integer(as.vector(hist$Freq));
	if(logval[1]==0){
		logval=logval[-1];
		freq=freq[-1];
	}
	for(i in 1:length(logval)){
		logcounts[logval[i],samp_idx]=freq[i];
	}
}

#cat("Log counts:\n");
#print(logcounts);
num_cat=sum(logcounts[,num_samples_wAvg]);
totals=apply(logcounts, 2, sum);

cat("Max height (sum of all category counts): ", num_cat, "\n");

logcols=rainbow(max_lognorm, end=.65);


#(bottom, left, top, right);
par(oma=c(12,3,2,7));
par(mar=c(1,4,2,5));
layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,2), nrow=1));

perc_mark=c(0,.1,.25,.333,.50,.666,.75,.9,1);
perc_lab=sprintf("%2i%%", perc_mark*100);
perc_pos=num_cat*perc_mark;

barpos=barplot(rep(0, ncol(logcounts)), ylim=c(0, num_cat*1.1), axes=F, axisnames=F);

axis(side=4, at=perc_pos, labels=perc_lab, las=2);
abline(h=perc_pos, col="grey70", lty=2);
barpos=barplot(logcounts, names.arg=sample_names,
	las=2, cex.names=resize, col=logcols, main=InputFileName, ylim=c(0, num_cat*1.1),
	ylab="Unique Taxa Counts", add=T
);

text(barpos, totals, totals, adj=c(.5,-1), cex=resize);

par(mar=c(1,0,6,1));
legbarpos=barplot(rep(1, max_lognorm), col=logcols, horiz=T, axes=F);

rightmost=max(barpos);
cutoffs=1:max_lognorm;
cutoffs=sprintf(">1.0e-%i", cutoffs);
print(cutoffs);
text(rep(0,max_lognorm) , legbarpos, cutoffs, adj=c(-.15,0));


dev.off();

###############################################################################

writeLines("Done.\n")

q(status=0)
