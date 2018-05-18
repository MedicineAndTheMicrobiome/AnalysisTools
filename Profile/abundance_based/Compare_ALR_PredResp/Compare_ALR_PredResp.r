#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);

params=c(
	"pred_file", "x", 1, "character",
	"resp_file", "y", 1, "character",
	"output_file", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];
script_path=paste(head(strsplit(script_name, "/")[[1]], -1), collapse="/");

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-x <as_pred p-value matrix input file>\n",
	"	-y <as_resp p-value matrix input file>\n",
	"	-o <output file name root>\n",
	"\n",
	"This script will read in two matrices:\n",
	" 1.) Predictor p-values\n",
	" 2.) Response p-values\n",
	"\n",
	"The matrices should have the format:\n",
	"\n",
	"Rows: categories (taxa)\n",
	"Cols: factors (variables)\n",
	"\n",
	"The numbers factors and categories should be the same\n",
	"(as well as the covariates) so that the values can\n",
	"be comparable.\n",
	"\n", sep="");

if(!length(opt$pred_file) || !length(opt$resp_file) || !length(opt$output_file)){
	cat(usage);
	q(status=-1);
}


PredFile=opt$pred_file;
RespFile=opt$resp_file;
OutputRoot=opt$output_file;

input_summary=capture.output({
	cat("\n");
	cat("Input As Predictor File: ", PredFile, "\n");
	cat("Input As Response File: ", RespFile, "\n");
	cat("Output File Root: ", OutputRoot, "\n", sep="");
	cat("\n");
});
cat(input_summary, sep="\n");

##############################################################################

plot_text=function(strings){
	par(family="Courier");
	par(oma=rep(.5,4));
	par(mar=rep(0,4));

	num_lines=length(strings);
	
	top=max(as.integer(num_lines), 52);

	plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
		xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
		);

	text_size=max(.01, min(.8, .8 - .003*(num_lines-52)));
	#print(text_size);

	for(i in 1:num_lines){
		#cat(strings[i], "\n", sep="");
		strings[i]=gsub("\t", "", strings[i]);
		text(0, top-i, strings[i], pos=4, cex=text_size); 
	}
}


##############################################################################
# Open PDF output

pdf(paste(OutputRoot,".pred_vs_resp.pdf", sep=""), height=8.5, width=8.5);

plot_text(input_summary);

#############################################################################	

as_pred=read.table(PredFile);
as_resp=read.table(RespFile);

num_pred_fact=ncol(as_pred);
num_pred_cat=nrow(as_pred);
pred_fact_names=colnames(as_pred);
pred_cat_names=rownames(as_pred);

num_resp_fact=ncol(as_resp);
num_resp_cat=nrow(as_resp);
resp_fact_names=colnames(as_resp);
resp_cat_names=rownames(as_resp);

shrd_fact_names=intersect(pred_fact_names, resp_fact_names);
shrd_cat_names=intersect(pred_cat_names, resp_cat_names);

var_summary=capture.output({
	cat("As Predictor:\n");
	cat("Categories:\n");
	print(pred_cat_names);
	cat("Factors:\n");
	print(pred_fact_names);
	cat("\n");

	cat("As Response:\n");
	cat("Categories:\n");
	print(resp_cat_names);
	cat("Factors:\n");
	print(resp_fact_names);
	cat("\n");
});

shrd_var_summary=capture.output({
	cat("Shared:\n");
	cat("Categories:\n");
	print(shrd_cat_names);
	cat("Factors:\n");
	print(shrd_fact_names);
});


cat(var_summary, sep="\n");
cat(shrd_var_summary, sep="\n");
plot_text(var_summary);
plot_text(shrd_var_summary);

#############################################################################	

plot_resp_pred_scatter=function(factor_name, categories, as_pred_mat, as_resp_mat){


	par(family="");
	par(mar=c(5,5,8,3));
	x_val=-log10(as_pred_mat[categories, factor_name]);
	y_val=-log10(as_resp_mat[categories, factor_name]);

	# Specify reference significance plots
	siglines_val=c(.1, .05, .01, .001);
	siglines_pos=-log10(siglines_val);

	xrange=range(c(x_val, siglines_pos));
	yrange=range(c(y_val, siglines_pos));

	plot(x_val, y_val, xlab="As Predictor", ylab="As Response",
		main=factor_name,
		xlim=c(xrange[1]-.25, xrange[2]+1), 
		ylim=c(yrange[1]-.25, yrange[2]+1)
	);

	# Above untransformed p-values
	axis(side=3, at=siglines_pos, labels=siglines_val, cex.axis=.7);
	axis(side=4, at=siglines_pos, labels=siglines_val, cex.axis=.7);

	abline(a=0, b=1, lwd=2, col="#8888FF");
	abline(h=siglines_pos, col="grey");
	abline(v=siglines_pos, col="grey");

	# Grow label size if it is more significant
	labelsize=sqrt(x_val^2 + y_val^2)/2;
	num_cat=length(categories);
	for(i in 1:num_cat){
		labelsize[i]=max(labelsize[i], .4);
		labelsize[i]=min(labelsize[i], 2);
	}

	# Label categories/taxa
	text(x_val, y_val, categories, pos=3, cex=labelsize);

}


#############################################################################	

for(cur_fact in shrd_fact_names){
	plot_resp_pred_scatter(cur_fact, shrd_cat_names, as_pred, as_resp);
}


#############################################################################	
# Close PDF output

dev.off();

#############################################################################

cat("Done.\n");
print(warnings());
q(status=0);
