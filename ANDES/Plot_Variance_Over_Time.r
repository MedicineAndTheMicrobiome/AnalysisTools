#!/usr/bin/env Rscript

###############################################################################
#                                                                             # 
#       Copyright (c) 2009 J. Craig Venter Institute.                         #     
#       All rights reserved.                                                  #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                             #
###############################################################################

###############################################################################

library('getopt');

params=c(
		"input_file", "i", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage <- paste (
		"\nUsage:\n", script_name, "\n",
		"	-i <input std deviation over time file>\n",
		"\n",
		"\n");

if(!length(opt$input_file)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_file;
OutputFleName=paste(InputFileName, ".plot.pdf", sep="");

cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name: ", OutputFleName, "\n");

library(MASS)

###############################################################################

plot_matrix=function(matrix, caption=NULL){

	# Transpose matrix and reverse order of columns so 0,0 is at top left instead of bottom left
	img_matrix=t(matrix);
	img_matrix=img_matrix[,ncol(img_matrix):1]
	image(img_matrix, xaxt="n", yaxt="n",
		xlab=caption,
		col=hsv(1,1:100/100)
		#col=rev(grey(0:100/100))
		#col=rev(rainbow(500, start=0, end=0.65, alpha=.4))
	);
	
	rownames=rownames(matrix);
	colnames=colnames(matrix);
	nrow=nrow(matrix);
	ncol=ncol(matrix);

	# Plot labels
	img_x_pos=((1:ncol)-1)/(ncol-1)
	img_y_pos=((1:nrow)-1)/(nrow-1)
	axis(side=3, at=img_x_pos, labels=colnames, tick=FALSE);
	axis(side=2, at=img_y_pos, labels=rev(rownames), las=2, tick=FALSE);

	# Go through and label cells with values
	for(i in 1:nrow){
		for(j in 1:ncol){
			xpos=(i-1)/(nrow-1);
			ypos=(j-1)/(ncol-1);
			text(ypos, 1-xpos, sprintf("%3.4f",matrix[i,j]));	
		}
	}

	
}

###############################################################################

compute_mean_reversion_intervals=function(ci_perc, x){
	bounds=numeric(2);
	num_series=length(x);
	names(x)=c();
	alpha=1-ci_perc;

	cat("---------------------------------------------------------------\n");
	cat("Computing ", ci_perc *100, "% Confidence Intervals.\n", sep="");
	cat("Mean of series: ", mean(x), "\n");
	cat("St. Dev of series: ", sd(x), "\n");
	cat("\n");

	# Regress the relationship between change and value
	diffs=diff(x);
	prev=x[1:(num_series-1)];
	regr=lm(diffs~prev);
	regr_summary=summary(regr);

	b=regr_summary$coefficients[1,1];
	m=regr_summary$coefficients[2,1];
	sigma=regr_summary$sigma;

	reversion_speed=-m;
	long_run_mean=b/reversion_speed;

	cat("Reversion Speed: ", reversion_speed, "\n");
	cat("Long Run Mean: ", long_run_mean, "\n");
	cat("Volatility: ", sigma, "\n");
	cat("R-squared: ", regr_summary$r.squared, "\n");
	cat("Spearman's Cor: ", cor(prev, diffs, method="spearman"), "\n");
	cat("Pearson's Cor: ", cor(prev, diffs, method="pearson"), "\n");

	plot(prev,diffs, xlab="Nominal Value", ylab="Change", main="Regression Results");
	abline(regr, col="blue", lty=2);
	abline(h=0, col="grey");
	abline(v=long_run_mean, col="grey");

	# Lower bound
	last_val=tail(x,1);
	bounds[1]=reversion_speed*(long_run_mean-last_val)+sigma*qnorm(alpha/2, mean=0, sd=1)+last_val; 
	# Upper bound
	bounds[2]=reversion_speed*(long_run_mean-last_val)+sigma*qnorm(1-alpha/2, mean=0, sd=1)+last_val;

	return(bounds);
}

###############################################################################

# Load data
in_data<-as.matrix(read.table(InputFileName, header=TRUE, check.names=FALSE));

# Convert text to numeric
data=apply(in_data[,2:ncol(in_data)], 2, as.numeric);
rownames(data)=in_data[,1];

sample_names=rownames(data);
region_names=colnames(data);

print(data);

pdf(OutputFleName, height=8.5, width=11);

max_var=max(data)
num_periods=nrow(data);
num_samples=ncol(data);

# Compute basic stats
means=apply(data, 2, mean);
stdev=apply(data, 2, sd);

# Plot times series
plot(0,0, ylim=c(0, max_var*1.05), xlim=c(0, num_periods+1), type="n", xaxt="n", xlab="Time", ylab="Variation");
points(rep(0,num_samples), means, pch=15, col=1:num_samples);
axis(1, at=0, labels="Means", font=2);

num_regions=ncol(data);
for(i in 1:num_regions){
	lines(data[,i], col=i);
}
for(i in 1:num_periods){
	abline(v=i, col="grey", lty=2);
}
legend(0, max_var, region_names, fill=1:num_regions);
axis(1, at=1:num_periods, labels=sample_names);

# Output raw table table/matrix
plot_matrix(data, caption="Variance over Time");

# Output correlation table/matrix
cor=cor(data);
plot_matrix(cor, caption="Correlation Coefficient");

# Compute autocorrelation
#par(mfrow=c(3,2));
#for(i in 1:num_regions){
#	acf(data[,i], main=region_names[i]);
#}

# Mean reversion analysis
par(mfrow=c(2,1));
for(i in 1:num_regions){
	bounds=compute_mean_reversion_intervals(.3, data[,i]);

	data_wPred=c(data[,i], (bounds[1]+bounds[2])/2);
	plot(data_wPred, main=region_names[i], ylab="St. Dev", type="b", xaxt="n");
	points(num_periods+1, bounds[1], pch=3, col="blue");
	points(num_periods+1, bounds[2], pch=3, col="blue");
	segments(num_periods+1, bounds[1], num_periods+1, bounds[2], col="blue")

	abline(h=means[i], lty=2, col="red")
	abline(h=means[i]+stdev[i], lty=3, col="red");
	abline(h=means[i]-stdev[i], lty=3, col="red");
	axis(1, at=1:(num_periods+1), labels=c(sample_names, "Next"));

}


######################################################################################################

dev.off()
print(warnings());
cat("Done.\n");
