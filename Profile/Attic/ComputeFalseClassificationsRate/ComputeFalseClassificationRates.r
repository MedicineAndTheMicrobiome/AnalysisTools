#!/usr/bin/env Rscript

###############################################################################

library('getopt');
library('MASS');

params=c(
	"input_matrix", "i", 1, "character",
	"regular_expression", "r", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage= paste(
	"\nUsage:\n\t", script_name, "\n\t\t<Input summary_table.xls FileName>\n\n",
	"	-i <input summary_table.xls file>\n",
	"	-r <regular expression, eg. \"z.\">\n",
	"\n",
	"\n");

if(!length(opt$input_matrix) || !length(opt$regular_expression)){
	cat(usage);
	q(status=-1);
}

InputFileName=opt$input_matrix;
RegularExp=opt$regular_expression;
OutputFileName=InputFileName;


###############################################################################
# Load counts from file

cat("Working on ", InputFileName, "\n", sep="");
cat("Regex: ", RegularExp, "\n", sep="");

# Load summary_table.xls
mat=as.matrix(read.table(InputFileName, sep=",", header=TRUE, check.names=FALSE, row.names=1))

if(sum(mat<0)!=0){
	cat("There are negative distances in here!\n");	
	q(status=-1);
}

#print(mat);
nrow=nrow(mat);
ncol=ncol(mat);
if(nrow != ncol){
	cat("Error, matrix is not symmetric.\n");
	q(status=-1);
}
num_total_samples=nrow;
cat("Num samples in matrix:", num_total_samples, "\n");

# Get names
sample_names=rownames(mat);
print(sample_names);

###############################################################################

find_regex_groups=function(names, RegularExp){
	num_names=length(names);
	groups=c();
	for(i in 1:num_names){
		x=regexpr(RegularExp, names[i]);	
		pos=x[1];
		mlen=attr(x,"match.length");
		match_str=substr(names[i], pos, pos+mlen-1);
		#cat("Found: ", match_str, " in ", names[i], "\n");
		groups=c(groups, match_str);
	}
	return(groups);
}

compute_percent_misclassified=function(A, B){
	lenA=length(A);
	lenB=length(B);
	if(lenA!=lenB){
		cat("Error, length of classification vectors are not the same length.\n");
		cat(lenA, "!=", lenB, "\n");
		q(status=-1);
	}else{
		len=lenA;
	}

	mismatch=sum(abs(A-B));
	perc_mism=mismatch/len;
	if(perc_mism>.5){
		perc_mism=1-perc_mism;
	}

	return(perc_mism);
}

###############################################################################

groups=find_regex_groups(sample_names, RegularExp);
#print(groups);
unique_groups=unique(groups);
num_groups=length(unique_groups);
cat("Unique Groups: ", unique_groups, "\n");
groups_table=(table(groups));

group_names=(names(groups_table));

assigned_colors=rep(0, num_total_samples);

pdf(paste(OutputFileName, ".pdf", sep=""), height=8.5, width=14);
outfile=file(paste(OutputFileName, ".miscl_rate", sep=""), "w");
misclass_names_file=file(paste(OutputFileName, ".miscl_names", sep=""), "w");

for(grp_idx1 in 1:(num_groups-1)){
	for(grp_idx2 in (grp_idx1+1):num_groups){

		cur_group1=unique_groups[grp_idx1];
		cur_group2=unique_groups[grp_idx2];

		title=paste("Group ", grp_idx1, " (", cur_group1, ") vs Group ", grp_idx2, " (", cur_group2, ")", sep="");
		cat("Working on: ", title, "\n", sep="");

		group1_indices=(groups==cur_group1);
		group2_indices=(groups==cur_group2);

		group1_size=sum(group1_indices);
		group2_size=sum(group2_indices);

		#print(group1_size);
		#print(group2_size);
		if(group1_size<1 || group2_size<1){
			cat("Cannot analyze pair when either group has no members.\n");
			next;
		}
		if((group1_size + group2_size) < 3){
			cat("Cannot analyze groups if there are less than 3 members.\n");
			next;
		}

		combined_groups=(group1_indices | group2_indices);
		combined_groups_names=sample_names[combined_groups];

		# Assign colors 
		assigned_colors[group1_indices]="red";
		assigned_colors[group2_indices]="green";
		used_colors=assigned_colors[combined_groups];
		#print(sample_names[combined_groups]);
		#print(used_colors);

		# Assign samples to groups empirically
		group_ids=rep(NULL, num_total_samples);
		group_ids[group1_indices]=1;
		group_ids[group2_indices]=2;
		empirical_groups=group_ids[combined_groups];
		#print(empirical_groups);

		subset_dist=as.dist(mat[combined_groups, combined_groups]);
		hc=hclust(subset_dist);

		# Assign samples to groups computationally
		cuts=cutree(hc, k=2);
		computational_groups=as.vector(cuts);
		#print(as.vector(cuts));
		
		mism_rate=compute_percent_misclassified(empirical_groups, computational_groups);
		cat("Misclassification Rate: ", mism_rate, "\n", sep="");
		cat(file=outfile, paste(cur_group1, cur_group2, sprintf("%4.4f",mism_rate), group1_size, group2_size, sep=","), "\n");

		mismatch=(empirical_groups!=computational_groups);

		cat(file=misclass_names_file, title, "\n");
		cat(file=misclass_names_file, combined_groups_names[mismatch], sep="\n");
		cat(file=misclass_names_file, "\n");

		# Plot Dendrogram
		dendro_scale=34/num_total_samples;
		plot(hc,cex=dendro_scale, main=title);

		# Compute MDS plot
		subset_dist[subset_dist==0]=1e-323;
		mds=isoMDS(subset_dist);

		# Plot Points
		plot(mds$points, col=used_colors, main=title);

		# Plot Names
		plot(mds$points, col=used_colors, main=title, type="n");
		text(mds$points, labels=sample_names[combined_groups], col=used_colors, cex=.6);

	}
}




cat("Done.\n");
q(status=0);
