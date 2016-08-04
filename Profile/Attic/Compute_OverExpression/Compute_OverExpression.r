#!/usr/bin/env Rscript

###############################################################################
#                                                                             #
#       Copyright (c) 2014 J. Craig Venter Institute.                         #
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
	"genomic", "G", 1, "character",
	"transcriptomic", "T", 1, "character",
	"output_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-G <input genomic summary_table.csv>\n",
	"	-T <input transcriptomic summary_table.csv>\n",
	"	-o <output filename root>\n",
	"\n",	
	"This script will take in two summary tables, one for transcriptomic data\n",
	"and one for metagenomic data, and generate a heat map for which categories\n",
	"are over expressed.  The two files should have the same number of columns\n",
	"and matching sample IDs with the following format.\n",
	"\n",
	"For example:\n",
	"	<sample ID>.<label1>\n",
	"	<sample ID>.<label2>\n",
	"	<sample ID>.<label3>\n",
	"\n",
	"Where the <label> will be used to label the column in the output graphics.\n",
	"\n",
	"\n");

if(!length(opt$genomic) || !length(opt$transcriptomic) || !length(opt$output_root)){
	cat(usage);
	q(status=-1);
}

###############################################################################

GenomicFname=opt$genomic;
TranscriptomicFname=opt$transcriptomic;
OutputRoot=opt$output_root;

cat("\n")
cat("             Genomic: ", GenomicFname, "\n");
cat("      Transcriptomic: ", TranscriptomicFname, "\n");
cat("Output Filename Root: ", OutputRoot, "\n");

pdf(paste(OutputRoot, ".pdf", sep=""), height=8.5, width=11);

###############################################################################
###############################################################################

load_summary_file=function(InputFileName){
	# Load data
	inmat=as.matrix(read.delim(InputFileName, sep="\t", header=TRUE, check.names=FALSE,
		 comment.char="", quote="", row.names=1))

	#cat("Original Matrix:\n")
	#print(inmat);

	if(ncol(inmat)==1){
		cat("\nError: There are no categories in this file:\n\n");
		cat("\t", InputFileName, "\n\n\n", sep="");
		quit(status=-1);
	}

	if(ncol(inmat)==2){
		cat("\nError: There is only one category in this file:\n\n");
		cat("\t", InputFileName, "\n\n\n", sep="");
		quit(status=-1);
	}

	# Grab columns we need into a vector, ignore totals, we won't trust it.
	counts_mat=inmat[,2:(ncol(inmat)), drop=F];
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

compute_logratio=function(trans_norm_mat, genom_norm_mat){
	ncol=ncol(trans_norm_mat);
	nrow=nrow(trans_norm_mat);
	lr_mat=log(trans_norm_mat/genom_norm_mat);
	colnames(lr_mat)=colnames(trans_norm_mat);
	rownames(lr_mat)=rownames(trans_norm_mat);
	return(lr_mat);
}

###############################################################################

pad_missing=function(inmat, complete_cat){
	inmat_cat=colnames(inmat);
	num_complete=length(complete_cat);

	nsamples=nrow(inmat);
	outmat=matrix(-1, nrow=nsamples, ncol=num_complete);

	for(i in 1:num_complete){

		matching=which(complete_cat[i]==inmat_cat);
		if(length(matching)>1){
			cat("\nWARNING: Duplicate categories: ", paste(matching, collapse=", "), "\n", sep="");
			cat("\t", complete_cat[i], "\n");
			print(inmat[,matching]);
			cat("\nKeeping first column\n");
			matching=min(matching);
		}

		if(any(matching)){
			outmat[,i]=inmat[,matching];
		}else{
			outmat[,i]=rep(0, nsamples);
		}
	}

	colnames(outmat)=complete_cat;
	rownames(outmat)=rownames(inmat);

	#cat("Inmat\n");
	#print(inmat)
	#cat("Outmat\n");
	#print(outmat);

	return(outmat);
		
}

###############################################################################

trim_common=function(names){
	num_names=length(names);

	splits=strsplit(names, ":");

	# Determine number of splits
	split_counts=numeric();
	for(i in 1:num_names){
		split_counts[i]=length(splits[[i]]);
	}
	max_splits=max(split_counts);

	# Determine how many splits are identical
	rem_ix=0;
	for(i in 1:max_splits){
		ref=splits[[1]][i];
		all_eq=T;
		for(j in 2:num_names){
			if(splits[[j]][i]!=ref){
				all_eq=F;
				break;
			}
		}
		if(all_eq){
			rem_ix=rem_ix+1;
		}else{
			break;
		}
	}

	cat("Number of redundant category name components: ", rem_ix, "\n");

	if(rem_ix>0){
		common_str=paste(splits[[1]][1:rem_ix], sep=":");
	}else{
		common_str="";
	}
	cat("Common prefix: ", common_str, "\n");	

	results=list();
	# Return common string and cleaned up array of names
	results[["common"]]=common_str;
	results[["different"]]=gsub( paste(common_str, ":", sep=""), "", names);
	return(results);
}

###############################################################################

trans_count_mat=load_summary_file(TranscriptomicFname);
genom_count_mat=load_summary_file(GenomicFname);

unique_categories=unique(c(colnames(trans_count_mat), colnames(genom_count_mat)));

trans_count_mat=pad_missing(trans_count_mat, unique_categories);
genom_count_mat=pad_missing(genom_count_mat, unique_categories);

if(nrow(trans_count_mat)!=nrow(genom_count_mat)){
	cat("Error, number of samples (rows) do not match.\n");
	quit(status=-1);
}

if(0){
	cat("\n");
	cat("Transcriptomic Counts:\n");
	print(trans_count_mat);
	cat("\n");
	cat("Genomic Counts:\n");
	print(genom_count_mat);
	cat("\n");
}

num_samples=nrow(trans_count_mat);
num_categories=ncol(trans_count_mat);
category_names=colnames(trans_count_mat);
sample_names=rownames(trans_count_mat);

cat("Num Samples:", num_samples, "\n");
cat("Num Categories: ", num_categories, "\n");
cat("\n");
cat("Sample Names: \n");
print(sample_names);

# Extract labels
	
trim_results=trim_common(category_names);
common_prefix=trim_results$common;
category_names=trim_results$different;
colnames(trans_count_mat)=category_names;
colnames(genom_count_mat)=category_names;

cat("New category names: \n");
print(category_names);


###############################################################################

# Normalize
trans_norm_mat=normalize(trans_count_mat);
genom_norm_mat=normalize(genom_count_mat);

if(0){
	cat("\n");
	cat("Transcriptomic Normalized:\n");
	print(trans_norm_mat);
	cat("\n");
	cat("Genomic Normalized:\n");
	print(genom_norm_mat);
	cat("\n");
}

###############################################################################

plot_lr_mat=function(lr_mat, title=""){
	nrows=nrow(lr_mat);
	ncols=ncol(lr_mat);

	par(mar=c(20,.5,3,15));

	plot(0,0, ylim=c(.5, nrows+.5), xlim=c(.5, ncols+.5), type="n",
		xaxt="n", yaxt="n", 
		bty="n",
		xlab="", ylab="",
		main=title
	);

	vect=as.vector(lr_mat);
	vect=vect[is.finite(vect)];
	max_diff=max(abs(vect));
	cat("Max ratio: ", max_diff, "\n");

	num_levels=10;
	neg_colors=rainbow(num_levels, start=1/6, end=2/6);
	pos_colors=rev(rainbow(num_levels, start=0/6, end=1/6));
	thresholds=seq(0, max_diff, length.out=num_levels);
	#print(thresholds);

	colmat=matrix(0, nrow=nrows, ncol=ncols);
	for(i in 1:nrows){
		for(j in 1:ncols){
			val=lr_mat[i,j];
			if(is.finite(val)){
				col_idx=max(which(abs(val)>thresholds));
				if(val>0){	
					colmat[i,j]=pos_colors[col_idx];
				}else if(val<0){
					colmat[i,j]=neg_colors[col_idx];
				}else{
					colmat[i,j]=pos_color[0];
				}
			}else{
				if(is.nan(val)){
					colmat[i,j]="grey";
				}else{
					if(val==Inf){
						colmat[i,j]=pos_colors[num_levels];
					}
					if(val==-Inf){
						colmat[i,j]=neg_colors[num_levels];
					}
				}
			}
			
		}
	}

	maxdim=max(nrows, ncols);
	for(y in 1:nrows){
		for(x in 1:ncols){
			points(x, nrows-y+1, pch=15, cex=70/maxdim, col=colmat[y, x]); 
			if(is.infinite(abs(lr_mat[y, x]))){
				points(x, nrows-y+1, pch=15, cex=.6*70/maxdim, col="white"); 
			}
		}

	}

	# Label right axises
	sam_names=rownames(lr_mat);
	for(y in 1:nrows){
		text(ncols+.5, y, sam_names[nrows-y+1], srt=-45, xpd=T, cex=1.5, pos=4, font=2);
	}

	# Label bottom axises
	cat_names=colnames(lr_mat);
	labelsize=min(c(1.5, 20/ncols));
	for(x in 1:ncols){
		text(x, .5, cat_names[x], srt=-45, xpd=T, cex=labelsize, pos=4, font=3, offset=c(0,-.5));
	}


}

###############################################################################

# compute logs
lr_mat=compute_logratio(trans_norm_mat, genom_norm_mat);

lr_vect=as.vector(lr_mat);
finite_lr=is.finite(lr_vect);

min_nzero=min(lr_vect[finite_lr]);
min_nzero=min(c(-.1, min_nzero));

max_nzero=max(lr_vect[finite_lr]);
max_nzero=max(c(.1, max_nzero));

# Assign scores to -Inf, Inf, and NA, so categories can be sorted
score_matrix=matrix(0, nrow=num_samples, ncol=num_categories);
for(j in 1:num_categories){
	for(i in 1:num_samples){
		score=lr_mat[i, j];
		if(is.na(score)){
			score=0;
		}else{
			if(score==-Inf){
				score=min_nzero*1.5;
			}else if(score==Inf){
				score=max_nzero*1.5;
			}
		}
		score_matrix[i, j]=score;	
	}	
}
print(score_matrix)
col_score=apply(score_matrix, 2, sum);

cat("\nAssigned sort scores:\n");
print(col_score)

lr_mat=lr_mat[,order(col_score, decreasing=T)];


# Plot all
plot_lr_mat(lr_mat, title=common_prefix);

# Plot top and bottom 
TOP_TO_PLOT=20;
if(num_categories>20){
	plot_lr_mat(lr_mat[, 1:TOP_TO_PLOT], 
		title=paste(common_prefix, ": Top ", TOP_TO_PLOT, sep=""));
	plot_lr_mat(lr_mat[, (num_categories-TOP_TO_PLOT+1):num_categories], 
		title=paste(common_prefix, ": Bottom ", TOP_TO_PLOT , sep=""));
}

quit(status=0);



###############################################################################

dev.off();
cat("Done.\n")
print(warnings());

q(status=0)
