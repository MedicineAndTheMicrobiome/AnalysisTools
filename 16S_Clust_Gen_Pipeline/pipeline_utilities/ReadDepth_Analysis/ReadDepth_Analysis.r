#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_depth_file", "i", 1, "character",
	"output_fname_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input depth file>\n",
	"	-o <output root>\n",
	"\n",	
	"\n",
	"This script will read in a file that contains\n",
	"sample ids, Read_Depth, Sample_Type, and generate\n",
	"bar plots and pairwise comparisons of their medians.\n",
	"\n",
	"Generates: \n",
	"  <output root>.ctrl_depth_analysis.pdf\n",
	"\n");

if(!length(opt$input_depth_file) && !length(opt$output_fname_root)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_depth_file;
OutputFileNameRoot=opt$output_fname_root;

###############################################################################

cat("\n")
cat("Input File Name: ", InputFileName, "\n");
cat("Output File Name Root: ", OutputFileNameRoot, "\n");       

pdf(paste(OutputFileNameRoot, ".ctrl_depth_analysis.pdf", sep=""), height=8.5, width=11);

###############################################################################

plot_text=function(strings){
        par(family="Courier");
        par(oma=rep(.1,4));
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

###############################################################################

plot_text(c(
	paste("Input File: ", InputFileName, sep="")
));	

###############################################################################
# Load data
cat("Loading Report: ", InputFileName, "\n");
mat=read.table(InputFileName, sep="\t", header=TRUE, row.names=1);

cat("Data:\n");
print(mat);

sample_names=rownames(mat);
sample_types=mat[,"Sample_Type"];
sample_depth=mat[,"Read_Depth"];
num_samples=nrow(mat);
max_depth=max(sample_depth);


cat("Sample Types:\n");
print(sample_types);
sample_colors=as.numeric(sample_types);
unique_sample_types=as.character(levels(sample_types));
num_sample_types=length(unique_sample_types);

#------------------------------------------------------------------------------
# Split samples and calc stats

by_sample_type=list();
stat_mat=matrix(0, nrow=num_sample_types, ncol=6);
rownames(stat_mat)=unique_sample_types;
colnames(stat_mat)=c("median", "N", "min", "max", "lb95", "ub95");

for( smp_typ in unique_sample_types){
	smp_ix=(sample_types==smp_typ);
	by_sample_type[[smp_typ]]=mat[smp_ix,];

	depths=by_sample_type[[smp_typ]][,"Read_Depth"];
	stat_mat[smp_typ, "median"]=median(depths);
	stat_mat[smp_typ, "N"]=length(depths);
	stat_mat[smp_typ, "min"]=min(depths);
	stat_mat[smp_typ, "max"]=max(depths);
	if(stat_mat[smp_typ, "N"]>=40){
		stat_mat[smp_typ, "lb95"]=quantile(depths, .0250);
		stat_mat[smp_typ, "ub95"]=median(depths, .9750);
	}
}
print(stat_mat);

#------------------------------------------------------------------------------
# Calculate wilcoxon pvalues
pvals=matrix(NA, nrow=num_sample_types, ncol=num_sample_types);
rownames(pvals)=unique_sample_types;
colnames(pvals)=unique_sample_types;
for(smp_typ_a in unique_sample_types){
	for(smp_typ_b in unique_sample_types){
		if(smp_typ_a==smp_typ_b){next;}
		depth_a=by_sample_type[[smp_typ_a]][,"Read_Depth"];
		depth_b=by_sample_type[[smp_typ_b]][,"Read_Depth"];
		res=wilcox.test(depth_a, depth_b);
		pvals[smp_typ_a, smp_typ_b]=res$p.value;
	}
}

###############################################################################

layout_mat=matrix(c(1,1,1,2), nrow=1, ncol=4);

#------------------------------------------------------------------------------
# Scatter plot
layout(layout_mat);
par(mar=c(2,4,4,1));
plot(0,0, type="n", xlim=c(-4,4), ylim=c(0, max_depth*1.1), 
	xaxt="n", xlab="", ylab="Read Depth", main="All Samples: Points Only");
jitter=rnorm(num_samples, 0, sqrt(1.5));
points(jitter, sample_depth, col=sample_colors);
par(mar=c(0,0,2,0));
plot(0,0, type="n", xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", bty="n",
	xlab="", ylab="", main="");
legend(x=0,y=1, legend=unique_sample_types, fill=1:num_sample_types, bty="n");

#------------------------------------------------------------------------------
# Labeled scatter plot
par(mfrow=c(1,1));
par(mar=c(2,4,4,1));
# Scatter plot
plot(0,0, type="n", xlim=c(-4,4), ylim=c(0, max_depth*1.2), 
	xaxt="n", xlab="", ylab="Read Depth", main="All Samples: Labeled");
points(jitter, sample_depth, col=sample_colors);
text(jitter, sample_depth, labels=sample_names, col=sample_colors, cex=.7, pos=3);


#------------------------------------------------------------------------------
# bar plots
par(mar=c(20,5,5,1));
par(lwd=2);
barmids=barplot(stat_mat[unique_sample_types, "median"], col="white", 
	border=1:num_sample_types,
	ylab="Read Depth", ylim=c(0, max_depth*1.1), las=2,
	main="Median Read Depth by Sample Type");
spacing=diff(barmids)[1];
cat("Bar spacings: ", spacing, "\n");

#------------------------------------------------------------------------------
# plot points
for(i in 1:num_sample_types){
	smp_type=unique_sample_types[i];
	depths=by_sample_type[[smp_type]][,"Read_Depth"];
	num_pts=length(depths);
	jitter=rnorm(num_pts, 0, spacing/8);
	print(barmids[i]+jitter);
	points(barmids[i]+jitter, depths, col=i);
}

#------------------------------------------------------------------------------
# Output stat tables
options(width=2000);
plot_text(c(
	"Sample Type Stats:",
	"",
	capture.output(print(stat_mat)),
	"",
	"",
	"Wilcoxon Rank Sum Test Pairwise P-values:",
	"",
	capture.output(print(pvals))
));

###############################################################################

dev.off();

q(status=0);
