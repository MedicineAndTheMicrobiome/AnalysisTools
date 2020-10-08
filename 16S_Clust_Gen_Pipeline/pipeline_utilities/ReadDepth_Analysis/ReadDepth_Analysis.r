#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
        "input_depth_file", "m", 1, "character",
        "input_groups_file", "g", 1, "character",
        "output_fname_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-m <input depth/metadata.tsv file>\n",
	"	-g <input groups file>\n",
	"	-o <output root>\n",
	"\n",	
	"\n",
	"This script will read in a file that contains\n",
	"sample ids, Read_Depth, Sample_Type, and generate\n",
	"bar plots and pairwise comparisons of their medians.\n",
	"Then it will read in the groups file, which contains\n",
	"the input read depth and analyze the attrition.\n",
	"\n",
	"Generates: \n",
	"  <output root>.ctrl_depth_analysis.pdf\n",
	"\n");

if(
        !length(opt$input_depth_file) &&
        !length(opt$input_groups_file) &&
        !length(opt$output_fname_root)
){
        cat(usage);
        q(status=-1);
}

InputDepthFileName=opt$input_depth_file;
InputGroupsFileName=opt$input_groups_file;
OutputFileNameRoot=opt$output_fname_root;

###############################################################################

cat("\n")
cat("Input Depth File Name: ", InputDepthFileName, "\n");
cat("Input Groups File Name: ", InputGroupsFileName, "\n");
cat("Output File Name Root: ", OutputFileNameRoot, "\n");

pdf(paste(OutputFileNameRoot, ".ctrl_depth_analysis.pdf", sep=""), height=8.5, width=14);

###############################################################################

project_ids_map=list(
        "0000"="MoBio.Powersoil.DNA.Ext.Neg",
        "0001"="PCR.Negative",
        "0002"="PCR.Positive.Zymo",
        "0003"="PCR.Positive.Ecoli.DH5a",
        "0004"="Stool.Test.Extraction",
        "0005"="Qiagen.FastStool.DNA.Ext.Neg",
        "0006"="MoBio.Ultraclean.DNA.Ext.Neg",
        "0007"="PCR.Positive.Other",
        "0008"="Vehicle.Control.Saline",
        "0009"="Vehicle.Control.Swab",
        "0010"="MoBio.PowerMicrobiome.DNA.Ext.Neg",
        "0011"="Zymo.IntctCells.Pos.Ctrl.For.Extr",
        "0012"="Ecoli.IntctCells.Pos.Ctrl.For.Extr",
        "0013"="Pos.Ctrl.For.Extr.Other",
        "0015"="Neg.Ctrl.Genomic.DNA.by.Invstgtr"
);

cat("Mapped Project IDs:\n");
print(project_ids_map);

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
	"Read Depth and Sample Attrition Analysis:",
	"",
	"",
        "Input Depth File: ", 
	paste("  ", InputDepthFileName, sep=""),
	"",
        "Input Groups File: ",
        paste("  ", InputGroupsFileName, sep=""),
	"",
        "Output File Root: ",
        paste("  ", OutputFileNameRoot, sep=""),
	"",
	"",
	paste("Date:", date())
));

###############################################################################

#------------------------------------------------------------------------------
# Load output read depths
cat("Loading Metadata/Depth Report: ", InputDepthFileName, "\n");
mat=read.table(InputDepthFileName, sep="\t", header=TRUE, row.names=1);

cat("\n");
cat("Output Read Depth:\n");
print(mat);

cat("\n");

#------------------------------------------------------------------------------
# Load input read depths from group file
input_read_count_fn=paste(OutputFileNameRoot, ".input_read_counts", sep="");
cmd=paste("cut -f 2 ", InputGroupsFileName,
        " | sort | uniq -c | gawk '{print $2 \"\\t\" $1}' > ", input_read_count_fn, sep="");

cat("\n");
print(cmd);
system(cmd);
cat("\n");

input_reads_counts=read.table(input_read_count_fn, sep="\t", header=F, row.names=1);

read_info_mat=cbind(input_reads_counts, 0, NA, NA);
colnames(read_info_mat)=c("InputReadDepth", "OutputReadDepth", "Difference", "PropLoss");
all_sample_ids=rownames(read_info_mat);

matching_sample_ids=intersect(rownames(mat), all_sample_ids);
read_info_mat[matching_sample_ids, "OutputReadDepth"]=mat[matching_sample_ids,"Read_Depth"];
read_info_mat[, "Difference"]=read_info_mat[,"InputReadDepth"]-read_info_mat[,"OutputReadDepth"];
read_info_mat[, "PropLoss"]=read_info_mat[, "Difference"]/read_info_mat[,"InputReadDepth"];

# Append parsed out project ID and descriptions
samp_ids=rownames(input_reads_counts);
proj_ids=c();
proj_dsc=c();
for(smp in samp_ids){
	pid=strsplit(smp, "\\.")[[1]][1];
	proj_ids=c(proj_ids, pid);
	desc=project_ids_map[[pid]];
	if(is.null(desc)){
		desc=paste("P_", pid, sep="");
	}
	proj_dsc=c(proj_dsc, desc);
}
read_info_mat[,"ProjID"]=proj_ids;
read_info_mat[,"ProjDesc"]=proj_dsc;

read_info_mat=as.data.frame(read_info_mat);
print(read_info_mat);

#------------------------------------------------------------------------------
# Extract overall sample info

num_samples=nrow(read_info_mat);
max_depth=max(read_info_mat[,"InputReadDepth"]);
sample_names=rownames(read_info_mat);
sample_types=read_info_mat[,"ProjDesc"];
unique_sample_types=unique(sample_types);
num_sample_types=length(unique_sample_types);
colors=1:num_sample_types;
names(colors)=unique_sample_types;

sample_colors=numeric();
for(i in 1:num_samples){
	sample_colors[i]=colors[sample_types[i]];
}

###############################################################################

#------------------------------------------------------------------------------
# Split samples and calc stats

by_sample_type=list();
stat_mat=matrix(NA, nrow=num_sample_types, ncol=6);
rownames(stat_mat)=unique_sample_types;
colnames(stat_mat)=c("median", "N", "min", "max", "lb95", "ub95");

for( smp_typ in unique_sample_types){
	smp_ix=(sample_types==smp_typ);
	by_sample_type[[smp_typ]]=read_info_mat[smp_ix,];

	depths=by_sample_type[[smp_typ]][,"OutputReadDepth"];

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
		depth_a=by_sample_type[[smp_typ_a]][,"OutputReadDepth"];
		depth_b=by_sample_type[[smp_typ_b]][,"OutputReadDepth"];
		res=wilcox.test(depth_a, depth_b);
		pvals[smp_typ_a, smp_typ_b]=res$p.value;
	}
}

###############################################################################
###############################################################################

layout_mat=matrix(c(1,1,1,2), nrow=1, ncol=4);

#------------------------------------------------------------------------------
# Scatter plot
layout(layout_mat);
par(mar=c(2,4,4,1));
plot(0,0, type="n", xlim=c(-4,4), ylim=c(0, max_depth*1.1), 
	xaxt="n", xlab="", ylab="Read Depth", main="All Samples Post-Mothur: Points Only");

for(smp_typ in unique_sample_types){
	medval=stat_mat[smp_typ, "median"];
	abline(h=medval, col=colors[smp_typ], lwd=.8, lty=2);
	axis(side=4, at=medval, labels=medval, las=2, col.axis=colors[smp_typ], cex=.8);
}

jitter=rnorm(num_samples, 0, sqrt(1.5));
points(jitter, read_info_mat[,"OutputReadDepth"], col=sample_colors);
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
	xaxt="n", xlab="", ylab="Read Depth", main="All Samples Post-Mothur: Labeled");
points(jitter, read_info_mat[,"OutputReadDepth"], col=sample_colors);
text(jitter, read_info_mat[,"OutputReadDepth"], labels=sample_names, col=sample_colors, cex=.7, pos=3);


#------------------------------------------------------------------------------
# Bar plots
par(mar=c(20,8,5,1));
par(lwd=2);
barmids=barplot(stat_mat[unique_sample_types, "median"], col="white", 
	border=1:num_sample_types,
	ylim=c(0, max_depth*1.1), las=2,
	main="Median Read Depth Post-Mothur by Sample Type");
title(ylab="Read Depth", line=6);
spacing=diff(barmids)[1];
cat("Bar spacings: ", spacing, "\n");

#------------------------------------------------------------------------------
# Plot points
for(i in 1:num_sample_types){
	smp_type=unique_sample_types[i];
	depths=by_sample_type[[smp_type]][,"OutputReadDepth"];
	num_pts=length(depths);
	jitter=rnorm(num_pts, 0, spacing/8);
	print(barmids[i]+jitter);
	points(barmids[i]+jitter, depths, col=i);
	axis(side=1, at=barmids[i], tick=F, line=-1.2, cex.axis=.6, 
		labels=sprintf("Med Depth=%g", stat_mat[smp_type, "median"]), col.axis=colors[i]);
}

format=function(mat){
	outmat=apply(mat, 1:2, function(x){
			if(!is.na(x)){
				sprintf("%6.5f",x);
			}else{
				return("-");
			}
		}
	);
	return(outmat);
}

signif=function(mat){
	outmat=apply(mat, 1:2, function(x){
			if(!is.na(x)){
				if(x<=.001){return("****")}
				if(x<=.01){return("***")}
				if(x<=.05){return("**")}
				if(x<=.1){return("*")}
				return(".");
			}else{
				return(".");
			}
		}
	);
	return(outmat);
}

#------------------------------------------------------------------------------
# Output stat tables
options(width=2000);
plot_text(c(
	"Post-Mothur Sample Type Stats:",
	"",
	capture.output(print(stat_mat)),
	"",
	"",
	"Wilcoxon Rank Sum Test Pairwise P-values:",
	"",
	capture.output(print(format(pvals), quote=F)),
	"",
	"",
	"Significantly Different (p-value<0.10):",
	"",
	capture.output(print(signif(pvals), quote=F))
));

#------------------------------------------------------------------------------

max_rows=40;
num_pages=ceiling(num_samples/max_rows);
for(i in 1:num_pages){
	start=((i-1)*max_rows)+1;
	stop=min(start+(max_rows-1), num_samples);
	cat("Range: ", start, "-", stop, "\n");
	plot_text(c(
		"Read Attrition Table:",
		paste("Sample Range: ", start, " - ", stop, sep=""),
		"",
		capture.output(print(read_info_mat[start:stop,,drop=F]))
	));

}

#------------------------------------------------------------------------------
# Plot lines between In/Out

all_depths=c(read_info_mat[,"OutputReadDepth"], read_info_mat[,"InputReadDepth"]);
min_depth=min(all_depths);
max_depth=max(all_depths);
depth_range=max_depth-min_depth;

par(mar=c(3,5,5,5));
plot(0, xlim=c(0-.05,1+.05), ylim=c(min_depth-.1*depth_range, max_depth+.1*depth_range),
	main="Pre/Post Mothur Read Depth", xlab="", ylab="", xaxt="n", type="n");
mtext("Pre-Mothur Read Depth", side=2, line=3);
mtext("Post-Mothur Read Depth", side=4, line=3);
for(i in 1:num_samples){
	points(
		c(0,1),
		c(read_info_mat[i,"InputReadDepth"], read_info_mat[i,"OutputReadDepth"]), 
		col=sample_colors[i], type="b", lwd=.5);
}

#------------------------------------------------------------------------------
# Plot scatter In vs Out
par(mar=c(5,5,5,2));
plot(0, type="n", 
	xlim=c(min_depth-.1*depth_range, max_depth+.1*depth_range),
	ylim=c(min_depth-.1*depth_range, max_depth+.1*depth_range),
	xlab="Pre-Mothur Read Depth",
	ylab="Post-Mothur Read Depth",
	main="Pre/Post Mothur Read Scatter Plot"
);
abline(a=0, b=1, lwd=.5, col="grey", lty=2);
points(read_info_mat[,"InputReadDepth"], read_info_mat[,"OutputReadDepth"], 
	col=sample_colors, type="p");

#------------------------------------------------------------------------------
# histogram prop loss by projdesc
med_prop_loss=numeric(num_sample_types);
names(med_prop_loss)=unique_sample_types;

# Calculate medians
for(smp_type in unique_sample_types){
	typ_ix=read_info_mat[,"ProjDesc"]==smp_type;
	med=median(read_info_mat[typ_ix, "PropLoss"]);
	med_prop_loss[smp_type]=med;
}

par(mar=c(20,8,6,1));
mids=barplot(med_prop_loss, ylim=c(-.05, max(read_info_mat[, "PropLoss"])*1.05), 
	ylab="Proportion Loss",
	col="white", border=colors, las=2, main="Percent Read Loss by Sample Type (Less is Better)");
bar_sep=diff(mids)[1];
i=1;
for(smp_type in unique_sample_types){
	typ_ix=read_info_mat[,"ProjDesc"]==smp_type;
	vals=read_info_mat[typ_ix, "PropLoss"];
	num_vals=length(vals);
	jitter=rnorm(num_vals, 0, bar_sep/7);
	points(mids[i]+jitter, vals, col=colors[i]);
	text(mids[i], 0, labels=sprintf("Med Loss=%5.4f", med_prop_loss[smp_type]), pos=1, col=colors[i]);
	i=i+1;
}

#------------------------------------------------------------------------------
# scatter diff vs prop loss

par(mar=c(5,5,5,2));
plot(0, type="n", 
	main="Pre-Mothur Read Depth vs. Attrition",
	xlab="Pre-Mothur Read Depth", ylab="Proportion of Read Attrition",
	xlim=c(0, max(read_info_mat[, "InputReadDepth"])*1.05),
	ylim=c(0, max(read_info_mat[, "PropLoss"])*1.05),
);
abline(h=median(read_info_mat[, "PropLoss"]), lwd=.5, col="grey", lty=2);
points(read_info_mat[,"InputReadDepth"], read_info_mat[, "PropLoss"],
	col=sample_colors
);

#------------------------------------------------------------------------------
# List samples that were "lost"

zero_output_depth_ix=read_info_mat[,"OutputReadDepth"]==0;
sample_names=rownames(read_info_mat);
lost_samples=sample_names[zero_output_depth_ix];
if(length(lost_samples)==0){
	lost_samples="<none>";
}
plot_text(c(
	"Samples with zero Post-Mothur Reads:",
	"",
	lost_samples
));
	
###############################################################################

dev.off();

q(status=0);
