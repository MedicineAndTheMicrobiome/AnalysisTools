#!/usr/bin/env Rscript

###############################################################################
#                                                                             #
#       Copyright (c) 2014 J. Craig Venter Institute.                         #
#                                                                             #
#	Copyright (c) 2017 The Center for Medicine and the Microbiome         #
#                          University of Pittsburgh                           #
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

COMBO_TOP=100;
IND_TOP=50;
NUM_BS=2000;
REMAINING="Remaining";

library('getopt');

params=c(
	"profA", "a", 1, "character",
	"profB", "b", 1, "character",
	"nameA", "A", 2, "character",
	"nameB", "B", 2, "character",
	"output_root", "o", 1, "character",
	"resize_off", "r", 2, "logical",
	"num_bs", "s", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-a <input a summary_table.tsv>\n",
	"	-b <input b summary_table.tsv>\n",
	"	[-A <Preferred name for A, default=A>\n",
	"	[-B <Preferred name for B, default=B>\n",
	"	[-r (turn off autoresize)]\n",
	"	-o <output filename root>\n",
	"	[-s <number of bootstraps, default=", NUM_BS, ">]\n",
	"\n",	
	"This script will generate a probability based Venn diagram\n",
	"for each taxa in the files.  Probabilities are computed with\n",
	"bootstrapping both the depth and samples between the two input\n",
	"profiles.\n",
	"\n",
	"Three graphs will be generated:\n",
	"   Skewered Venn:\n",
	"	Each category will be broken down into the probability\n",
	"	of it being resampled upon rerunning the experiment.\n",
	"		A will be on the left off center.\n",
	"		B will be on the right off center.\n",
	"		A&B will be in the center.\n",
	"		!(A&B) will be on the right side.\n",
	"\n",
	"   Compressed Venn:\n",
	"	Similar to the skewered venn excepted there is no spacing\n",
	"	between each of the 4 set categories.\n",
	"\n",
	"   Most Likely Outcome Venn Square:\n",
	"	Based on the most likely of the 4 set categories, a\n",
	"	count is generated and the portions of the unit square\n",
	"	are sectioned off.\n",
	"\n",
	"\n",
	"Three tables will also be generated:\n",
	"   1.) The values underlying the first two bar venn graphs (.shared_prob.tsv)\n",
	"   2.) The values underlying the most likely outcome venn square (.likely_counts.tsv)\n",
	"   3.) A list of taxa for each of the 4 venn categories (.category_lists.tsv)\n",
	"   4.) A set of summary tables with only the categories in group A that are exclusive\n",
	"       to A at various likelihood cutoffs (###.summary_table.tsv)\n",
	"\n");

if(!length(opt$profA) || !length(opt$profB) || !length(opt$output_root)){
	cat(usage);
	q(status=-1);
}

###############################################################################

ProfAFname=opt$profA;
ProfBFname=opt$profB;
OutputRoot=opt$output_root;

nameA=opt$nameA;
nameB=opt$nameB;

if(length(nameA)==0){
	nameA=tail(strsplit(ProfAFname, "/")[[1]], 1);
}
if(length(nameB)==0){
	nameB=tail(strsplit(ProfBFname, "/")[[1]], 1);
}

if(length(opt$resize_off)==0){
	Resize=T;
}else{
	Resize=F;
}

if(length(opt$num_bs)==0){
	NumBS=NUM_BS;
}else{
	NumBS=opt$num_bs;
}

cat("\n")
cat("Prof A: ", ProfAFname, "\n");
cat("Prof B: ", ProfBFname, "\n");
cat("Output Filename Root: ", OutputRoot, "\n");
cat("Number of Bootstraps: ", NumBS, "\n");
cat("\n");


###############################################################################
###############################################################################

load_idmap=function(MapFname){
	
	cat("Loading ID Mapping: ", MapFname, "\n");

	mat=as.matrix(read.delim(MapFname, header=F, sep="\t", na.strings=""));
	nrows=nrow(mat);

	cat("Num rows: ", nrows, "\n");

	mat[,1]=gsub("^ +", "", mat[,1]);
	mat[,1]=gsub(" +$", "", mat[,1]);

	mapping=list(nrows);

	for(i in 1:nrows){
		id=mat[i,1];
		name=mat[i,2];
		mapping[[id]]=name;
	}
	#print(head(mapping));
	cat("Done loading ID Map.\n");
	return(mapping);
}

rename=function(names, mapping){
	num_names=length(names);
	mapnames=names(mapping);
	
	new_names=character(num_names);
	for(i in 1:num_names){
		idx=which(mapnames==names[i]);
		#cat("idx: ", idx, "\n");
		if(length(idx)!=0){
			new_names[i]=mapping[[names[i]]];
		}else{
			new_names[i]=names[i];
		}
	}
	#print(new_names);
	return(new_names);

}

###############################################################################
###############################################################################

load_summary_file=function(InputFileName){

	cat("Loading: ", InputFileName, "\n");

	# Load data
	inmat=as.matrix(read.delim(InputFileName, sep="\t", header=TRUE, check.names=FALSE,
		na.strings="", comment.char="", quote="", row.names=1))

	#cat("Original Matrix:\n")
	#print(inmat);

	# Grab columns we need into a vector, ignore totals, we won't trust it.
	counts_mat=inmat[,2:(ncol(inmat)), drop=FALSE];
	#print(counts_mat);

	rownames(counts_mat)=rownames(inmat);
	colnames(counts_mat)=colnames(inmat)[2:(ncol(inmat))];
	return(counts_mat);
}

###############################################################################

normalize=function(counts_mat){
	sample_totals=apply(counts_mat, 1, sum);
	normalized=matrix(0, nrow=nrow(counts_mat), ncol=ncol(counts_mat));
	for(i in 1:nrow(counts_mat)){
		normalized[i,]=counts_mat[i,]/sample_totals[i];
	}
	colnames(normalized)=colnames(counts_mat);
	rownames(normalized)=rownames(counts_mat);
	return(normalized);	
}

###############################################################################

resample_st=function(norm_mat, counts){

	num_cat=ncol(norm_mat);
	num_samp=nrow(norm_mat);

	new_st=matrix(NA, nrow=num_samp, ncol=num_cat);
	samp_ix=sample(1:num_samp, num_samp, replace=T);

	for(i in 1:num_samp){
		samp_prob=norm_mat[samp_ix[i], ];
		samp_counts=counts[samp_ix[i]];
		new_st[i,]=as.vector(rmultinom(1, size=samp_counts, prob=samp_prob));
	}
	
	colnames(new_st)=colnames(norm_mat);
	rownames(new_st)=paste("bs_", samp_ix, sep="");
	return(new_st);
}

###############################################################################

compute_shared_matrix=function(a_counts_mat, b_counts_mat){

	num_a_samples=nrow(a_counts_mat);
	num_b_samples=nrow(b_counts_mat);

	if(!setequal(colnames(a_counts_mat), colnames(b_counts_mat))){
		cat("Error: categories in A and B, do not match.\n");
	}

	# Sort by names
	sorted_names=sort(colnames(a_counts_mat));
	a_counts_mat=a_counts_mat[,sorted_names];
	b_counts_mat=b_counts_mat[,sorted_names];
	num_categories=ncol(a_counts_mat);

	# Normalize
	prof_a_norm_mat=normalize(prof_a_count_mat);
	prof_b_norm_mat=normalize(prof_b_count_mat);

	#print(prof_a_norm_mat);
	#print(prof_b_norm_mat);

	# Counts per sample
	prof_a_counts=apply(prof_a_count_mat, 1, sum);
	prof_b_counts=apply(prof_b_count_mat, 1, sum);

	#print(prof_a_counts);
	#print(prof_b_counts);

	num_bootstraps=NumBS;

	prop_matrix=matrix(0, nrow=num_categories, ncol=4);
	colnames(prop_matrix)=c("A", "B", "AB", "Neither");
	rownames(prop_matrix)=sorted_names;

	for(i in 1:num_bootstraps){

		a_resample=resample_st(prof_a_norm_mat, prof_a_counts);
		b_resample=resample_st(prof_b_norm_mat, prof_b_counts);

		#print(a_resample);
		#print(b_resample);
		
		a_presence=apply(a_resample, 2, function(x){sum(x)>0});		
		b_presence=apply(b_resample, 2, function(x){sum(x)>0});		

		a_only  = a_presence & !b_presence;
		b_only  = !a_presence & b_presence;
		a_and_b = a_presence & b_presence;
		neither = !a_presence & !b_presence;

		prop_matrix[,"A"]      =prop_matrix[,"A"]       +a_only;
		prop_matrix[,"B"]      =prop_matrix[,"B"]       +b_only;
		prop_matrix[,"AB"]     =prop_matrix[,"AB"]      +a_and_b;
		prop_matrix[,"Neither"]=prop_matrix[,"Neither"] +neither;
	
	}	
	
	prop_matrix=prop_matrix/num_bootstraps;
	return(prop_matrix);

}

###############################################################################

plot_shared_matrix_skewer=function(shared_matrix, aname="A", bname="B", exp_fact=1){
	num_cat=nrow(shared_matrix);

	par(mar=c(0,17,1,0));

	graph_margin=0.25
	right_pos=1.25;
	if(num_cat>50){
		num_rows=num_cat;
	}else{
		num_rows=50;
	}
	plot(x=c(-1-graph_margin, right_pos+graph_margin), y=c(1-.5, num_rows+.5), type="n",
		xlab="", ylab="", bty="n", yaxt="n", xaxt="n"
	);

	abline(v=0, col="grey");

	for(i in 1:num_cat){
		points(c(-1-graph_margin, right_pos),
			c(i,i), 
			type="l", col="grey", lty="dashed", lwd=.75);
	}

	txt_size=35/(num_rows)*exp_fact;
	txt_size=min(1, txt_size);
	bar_width=exp_fact*750/num_rows;

	for(i in 1:num_cat){
		
		# Shared
		ab=shared_matrix[i,"AB"];
		points(c(-ab/2,ab/2), c(i,i), col="purple", lend=1, type="l", lwd=bar_width);
		if(ab>0){
			text(0, i, sprintf("%3.3f", ab), col="white", cex=txt_size);
		}

		# AOnly
		a=shared_matrix[i, "A"];
		points(c(-ab/2-a,-ab/2), c(i,i), col="red", lend=1, type="l", lwd=bar_width);
		if(a>0){
			text(-ab/2-a, i, sprintf("%3.3f", a), col="black", pos=2, cex=txt_size, offset=.25);
		}

		# BOnly
		b=shared_matrix[i, "B"];
		points(c(ab/2,ab/2+b), c(i,i), col="blue", lend=1, type="l", lwd=bar_width);
		if(b>0){
			text(ab/2+b, i, sprintf("%3.3f", b), col="black", pos=4, cex=txt_size, offset=.25);
		}

		# Neither
		n=shared_matrix[i, "Neither"];
		points(c(right_pos,right_pos-n), c(i,i), col="black", lend=1, type="l", lwd=bar_width);
		if(n>0){
			text(right_pos, i, sprintf("%3.3f", n), col="black", pos=4, cex=txt_size, offset=.25);
		}

	}

	# Plot category name
	axis(side=2, at=1:num_cat, labels=rownames(shared_matrix), las=2, cex.axis=txt_size*1.05);

	legend(1+graph_margin, num_rows,
		legend=c(aname, bname, "Both", "Neither"),
		fill=c("red", "blue", "purple", "black"),
		xjust=1,
		yjust=0,
		cex=.6,
		bty="n"
	);

}

###############################################################################

plot_shared_matrix_compressed=function(shared_matrix, aname="A", bname="B", exp_fact=1){
	num_cat=nrow(shared_matrix);

	par(mar=c(0,17,0,0));

	graph_margin=0.05
	right_pos=1.25;

	if(num_cat>50){
		num_rows=num_cat;
	}else{
		num_rows=50;
	}

	plot(x=c(0-graph_margin, 1+graph_margin), y=c(1-.5, num_rows+.5), type="n",
		xlab="", ylab="", bty="n", yaxt="n", xaxt="n"
	);

	for(i in 1:num_cat){
		points(c(-1-graph_margin, right_pos),
			c(i,i), 
			type="l", col="grey", lty="dashed", lwd=.75);
	}

	txt_size=35/(num_rows)*exp_fact;
	txt_size=min(1, txt_size);

	bar_width=exp_fact*750/num_rows;

	for(i in 1:num_cat){

		# AOnly
		a=shared_matrix[i, "A"];
		begin=0;
		end=a;
		points(c(begin,end), c(i,i), col="red", lend=1, type="l", lwd=bar_width);
		if(a>.15){
			text((end-begin)/2, i, sprintf("%3.3f", a), col="green", cex=txt_size, offset=.25);
		}
		
		# Shared
		ab=shared_matrix[i,"AB"];
		begin=end;
		end=begin+ab;
		points(c(begin,end), c(i,i), col="purple", lend=1, type="l", lwd=bar_width);
		if(ab>.15){
			text(begin+(end-begin)/2, i, sprintf("%3.3f", ab), col="yellow", cex=txt_size);
		}


		# BOnly
		b=shared_matrix[i, "B"];
		begin=end;
		end=begin+b;
		points(c(begin,end), c(i,i), col="blue", lend=1, type="l", lwd=bar_width);
		if(b>.15){
			text(begin+(end-begin)/2, i, sprintf("%3.3f", b), col="orange", cex=txt_size, offset=.25);
		}

		# Neither
		n=shared_matrix[i, "Neither"];
		begin=end;
		end=begin+n;
		points(c(begin,end), c(i,i), col="black", lend=1, type="l", lwd=bar_width);
		if(n>.15){
			text(begin+(end-begin)/2, i, sprintf("%3.3f", n), col="white", cex=txt_size, offset=.25);
		}

	}

	# Plot category name
	axis(side=2, at=1:num_cat, labels=rownames(shared_matrix), las=2, cex.axis=txt_size*1.05);

	legend(1+graph_margin, num_rows,
		legend=c(aname, bname, "Both", "Neither"),
		fill=c("red", "blue", "purple", "black"),
		xjust=1,
		yjust=0,
		cex=.6,
		bty="n"
	);
}

###############################################################################

identify_most_likely_category=function(shared_matrix){
	num_categories=nrow(shared_matrix);

	most_like=function(probs_vect){
		max_prob=max(probs_vect);
		max_ix=which(max_prob==probs_vect);
		if(length(max_ix)==1){
			return(max_ix);
		}else{
			return(sample(max_ix, 1));
		}
	}

	most_likely=apply(shared_matrix, 1, most_like);
	most_likely=c("A", "B", "AB", "Neither")[most_likely];
	names(most_likely)=rownames(shared_matrix);
	return(most_likely);

}

###############################################################################

plot_square_venn=function(venn_cat, aname="A", bname="B", title="Square Venn"){

	margin=.25;
	par(mar=c(7,1,7,1));
	plot(0,0, xlim=c(0-margin, 1+margin), ylim=c(0-margin, 1+margin), type="n",
		xlab="", ylab="", bty="n", yaxt="n", xaxt="n"
	);
	rect(0,0, 1,1, lwd=3);

	sum=sum(venn_cat);
	norm_cat=venn_cat/sum;
	
	n=norm_cat["Neither"];
	rect(0,0, 1, n, col="black");

	a=norm_cat["A"]/(1-n);
	ab=norm_cat["AB"]/(1-n);
	b=norm_cat["B"]/(1-n);

	rect(0, n, a, 1, col="red");
	rect(a, n, a+ab, 1, col="purple");
	rect(a+ab, n, 1, 1, col="blue");

	legend(0, -.05, paste(
		c(aname, bname, "Both", "Unrecovered"), " (",
		c(venn_cat["A"], venn_cat["B"], venn_cat["AB"], venn_cat["Neither"]),
		")", sep=""
		),
		fill=c("red", "blue", "purple", "black")
	);

	text(.5, 1.2, title, cex=2, font=2);
}

###############################################################################


prof_a_count_mat=load_summary_file(ProfAFname);
prof_b_count_mat=load_summary_file(ProfBFname);

#if(length(MapFname)>0){
#	cat("Performing remapping of IDs...\n")
#	idmap=load_idmap(MapFname);
#	colnames(prof_a_count_mat)=rename(colnames(prof_a_count_mat), idmap);
#	colnames(prof_b_count_mat)=rename(colnames(prof_b_count_mat), idmap);
#	cat("Ok.\n");
#}

cat("\n");
num_a_samples=nrow(prof_a_count_mat);
num_a_categories=ncol(prof_a_count_mat);
cat("A: Num Samples:", num_a_samples, "\n");
cat("A: Num Categories: ", num_a_categories, "\n");
cat("\n");

num_b_samples=nrow(prof_b_count_mat);
num_b_categories=ncol(prof_b_count_mat);
cat("B: Num Samples:", num_b_samples, "\n");
cat("B: Num Categories: ", num_b_categories, "\n");
cat("\n");

in_A=apply(prof_a_count_mat, 2, function(x){sum(x)>0});
in_B=apply(prof_b_count_mat, 2, function(x){sum(x)>0});

# Make sure both profiles have the same columns in the same order
names=colnames(prof_a_count_mat);
names_sorted=sort(names);
prof_a_count_mat=prof_a_count_mat[,names_sorted];
prof_b_count_mat=prof_b_count_mat[,names_sorted];

# Get idea of absolute shared and unique counts
abs_A=(in_A & !in_B);
abs_B=(!in_A & in_B);
abs_AB=(in_A & in_B);
abs_Neither=(!in_A & !in_B);

cat("Observed in A:\n");
print(names[abs_A]);
cat("\n");
cat("Observed in B:\n");
print(names[abs_B]);
cat("\n");
cat("Observed in AB:\n");
print(names[abs_AB]);
cat("\n");
cat("Observed in Neither:\n");
print(names[abs_Neither]);
cat("\n\n");

###############################################################################

shared_matrix=compute_shared_matrix(prof_a_count_mat, prof_b_count_mat);

venn_cats_names=c("A", "B", "AB", "Neither");

###############################################################################

#rnd=sprintf(".%03i", sample(1000,1));
rnd="";

if(Resize){
	exp_fact=num_a_categories/100;
}else{
	exp_fact=1;
}
cat("Expansion factor: ", exp_fact, "\n");
pdf(paste(OutputRoot, rnd, ".pdf", sep=""), height=11*exp_fact, width=8.5);

# Print name of output file on top of each page
par(oma=c(0,0,.3,0));
setHook("plot.new", function(){mtext(OutputRoot, outer=T, line=-.5, cex=.5);});

###############################################################################
# Sort taxa by probabilities and sor

ix=order(shared_matrix[,"Neither"], decreasing=F);
shared_matrix=shared_matrix[ix,];
ix=order(shared_matrix[,"B"], decreasing=T);
shared_matrix=shared_matrix[ix,];
ix=order(shared_matrix[,"A"], decreasing=F);
shared_matrix=shared_matrix[ix,];
ix=order(shared_matrix[,"AB"], decreasing=F);
shared_matrix=shared_matrix[ix,];


plot_shared_matrix_skewer(shared_matrix, nameA, nameB, exp_fact);
plot_shared_matrix_compressed(shared_matrix, nameA, nameB, exp_fact);

###############################################################################
# Identify most likely venn category for each taxa

most_likely=identify_most_likely_category(shared_matrix);
venn_cat_counts=table(most_likely);

missing_venn_cat=setdiff(c("A", "B", "AB", "Neither"), names(venn_cat_counts));
venn_cat_counts[missing_venn_cat]=0;

print(venn_cat_counts);
plot_square_venn(venn_cat_counts, nameA, nameB, "Most Likely Outcomes Upon Resampling");

###############################################################################
# Output values gone into the skewer and compressed plot

fh=file(paste(OutputRoot, ".shared_prob.tsv", sep=""), "w");

hdr=paste("Category", nameA, nameB, "Both", "Neither", "MostLikely", sep="\t");
cat(file=fh, hdr, "\n", sep="");

cat_names=rownames(shared_matrix);
likely_names=venn_cats_names;

for(i in nrow(shared_matrix):1){

	cat(file=fh, cat_names[i], paste(shared_matrix[i, venn_cats_names], collapse="\t"), 
		most_likely[i], sep="\t");
	cat(file=fh, "\n");
}

close(fh);

###############################################################################
# Output values gone into the Venn square plot

fh=file(paste(OutputRoot, ".likely_counts.tsv", sep=""), "w");

cat(file=fh, "# Name", nameA, nameB, "Both", "Neither", sep="\t");
cat(file=fh, "\n");
cat(file=fh, OutputRoot, venn_cat_counts[venn_cats_names], sep="\t");
cat(file=fh, "\n");

close(fh);

###############################################################################
# Output Lists of taxa in each category

Alist=cat_names[most_likely=="A"];
Blist=cat_names[most_likely=="B"];
Bothlist=cat_names[most_likely=="AB"];
Neitherlist=cat_names[most_likely=="Neither"];

longest_list=max(venn_cat_counts);

fh=file(paste(OutputRoot, ".category_lists.tsv", sep=""), "w");

cat(file=fh, nameA, nameB, "Both", "Neither", sep="\t");
cat(file=fh, "\n\n");

for(i in 1:longest_list){

	c1=ifelse(i<=venn_cat_counts["A"], Alist[i], "");
	c2=ifelse(i<=venn_cat_counts["B"], Blist[i], "");
	c3=ifelse(i<=venn_cat_counts["AB"], Bothlist[i], "");
	c4=ifelse(i<=venn_cat_counts["Neither"], Neitherlist[i], "");

	cat(file=fh, paste(c(c1, c2, c3, c4), collapse="\t"));
	cat(file=fh, "\n");
} 

close(fh);

###############################################################################
# Generate summary tables at various likelihoods for "A"

write_summary_file=function(out_mat, fname){
        fc=file(fname, "w");
        cat(file=fc, paste("sample_id\ttotal", paste(colnames(out_mat), collapse="\t"), sep="\t"));
        cat(file=fc, "\n");
        sample_names=rownames(out_mat);
        num_samples=nrow(out_mat);
        for(samp_idx in 1:num_samples){
                total=sum(out_mat[samp_idx,]);
                outline=paste(sample_names[samp_idx], total,
                        paste(out_mat[samp_idx,], collapse="\t"), sep="\t");
                cat(file=fc, outline);
                cat(file=fc, "\n");
        }
        close(fc);
}

cat_names=rownames(shared_matrix);
for(lik in seq(.50,1.00,.10)){

	cat("Extracting categories with likelihoods greater than ", lik, "\n", sep="");

	lik_str=sprintf("%03.0f", lik*100);
	stname=paste(OutputRoot, ".", lik_str, ".summary_table.tsv", sep="");
	
	keep_categories=cat_names[shared_matrix[ , "A"]>=lik];

	num_kept_cat=length(keep_categories);
	cat("Number of categories found: ", num_kept_cat, "\n", sep=""); 
	keep_counts=prof_a_count_mat[, keep_categories, drop=F];

	write_summary_file(keep_counts, stname);
}



###############################################################################

dev.off();
cat("Done.\n")
if(!is.null(warnings())){
	print(warnings());
}
q(status=0)
