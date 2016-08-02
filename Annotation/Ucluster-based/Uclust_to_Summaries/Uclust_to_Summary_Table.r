#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_file_list", "i", 1, "character",
	"output_file_root", "o", 2, "character",
	"truncate_category_name", "t", 2, "logical"
	
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input uc clust file list>\n",
	"	-o <output file name root>\n",
	"	[-t (truncate category name to ID)]\n",
	"\n");

if(	!length(opt$input_file_list) || 
	!length(opt$output_file)
){
	cat(usage);
	q(status=-1);
}

InputFileList=opt$input_file_list;
OutputFileNameRoot=opt$output_file;

TruncateCategoryName=ifelse(length(opt$truncate_category_name)>0, T, F);

###############################################################################

cat("\n")
cat("Input File List: ", InputFileList, "\n");
cat("Output File Name Root: ", OutputFileNameRoot, "\n");       
cat("Truncate Category Names? ", TruncateCategoryName, "\n");
cat("\n");

###############################################################################

read_in_hits=function(fname){

	cat("Reading in: ", fname, ", ");
	data=read.delim(fname, header=F, sep="\t", quote=NULL, comment.char="#", as.is=T);
	
	# Only keep hits
	hits_idx=(data[,1]=="H");
	
	# Num no hits
	num_no_hits=sum(data[,1]=="N");

	cat("Num records acquired: ", nrow(data), "\n");

	# Only keep columns of interest
	keep_col=c(3, 4, 5, 9, 10);
	data_wHits=data[hits_idx, keep_col];
	colnames(data_wHits)=c("Length", "PctID", "Strand", "QryID", "Database");
	
	# Hit Record
	hit_rec=list();
	hit_rec[["num_no_hits"]]=num_no_hits;
	hit_rec[["hit_matrix"]]=data_wHits;
	return(hit_rec);

}

#------------------------------------------------------------------------------

read_in_file_list=function(fname){
	data=read.delim(fname, header=F, sep="\t", as.is=T);
	if(ncol(data)==1){
		assumed_name=character();
		for(i in 1:nrow(data)){
			assumed_name[i]=tail(strsplit(data[i,1], "/")[[1]],1);
			assumed_name[i]=gsub("\\.uc$", "", assumed_name[i]);
		}
		data=cbind(data,assumed_name);
	}
	colnames(data)=c("FileName", "GroupName");
	return(data);
}

#------------------------------------------------------------------------------

combine_summary_file_by_group=function(counts_matrix, group_matrix){
	cat("Combining counts by group...\n");
	print(counts_matrix);
	print(group_matrix);

	unique_groups=sort(unique(group_matrix[,"GroupName"]));
	num_uniq_grps=length(unique_groups);

	grp_cmb_matrix=matrix(NA, nrow=num_uniq_grps, ncol=ncol(counts_matrix));
	colnames(grp_cmb_matrix)=colnames(counts_matrix);
	rownames(grp_cmb_matrix)=unique_groups;

	for(i in 1:num_uniq_grps){
		grp_name=unique_groups[i];
		members=group_matrix[,"GroupName"]==grp_name;

		grp_cmb_matrix[grp_name,]=apply(counts_matrix[members,,drop=F], 2, sum);
	}
	
	return(grp_cmb_matrix);
}

###############################################################################

file_list=read_in_file_list(InputFileList);
#print(file_list);

num_files=nrow(file_list);
cat("Num files to read: ", num_files, "\n");

hit_list=list();
hits_tables=list();

unique_database_hits=character();
for(i in 1:num_files){
	hit_list[[i]]=read_in_hits(file_list[i, 1]);
	hits=(hit_list[[i]][["hit_matrix"]]["Database"])[,1];

	hits_tables[[i]]=table(hits);
	unique_database_hits=unique(c(unique_database_hits, hits));
}

cat("\n");
cat("Database Hits Found: \n");
print(unique_database_hits);

# Combine tables
combined_count_table=matrix(0, ncol=length(unique_database_hits)+1, nrow=num_files);
colnames(combined_count_table)=c(unique_database_hits, "Remaining");
rownames(combined_count_table)=file_list[,1];

cat("Making summary table out of counts.\n");

for(i in 1:num_files){
	table_names=names(hits_tables[[i]]);

	combined_count_table[i, table_names]=
		combined_count_table[i, table_names]+hits_tables[[i]][table_names];

	combined_count_table[i, "Remaining"]=
		hit_list[[i]][["num_no_hits"]];
}

#pdf(paste(OutputFileNameRoot, ".cluster_stats.pdf", sep=""), height=11, width=8.5);

#par(mfrow=c(3,2));

#num_unique_database_hits=length(unique_database_hits);
#for(dbix in 1:num_unique_database_hits){
	# Scatter plot length vs pct ID

#	for(sampix in 1:num_files){
#		cat("Plotting across samples for: ", unique_database_hits[i], "\n");
#		hits=(hit_list[[i]][["hit_matrix"]]["Database"])[,1];
		#hits=(hit_list[[i]][["hit_matrix"]]["Length"])[,1];
		#hits=(hit_list[[i]][["hit_matrix"]]["PctID"])[,1];
		#print(hits);
#	}
	
#}
#quit();

# Truncate names
if(TruncateCategoryName){
	cat("Truncating Names...\n");
	old=colnames(combined_count_table);
	new=character();
	for(i in 1:length(old)){
		new[i]=strsplit(old[i], " ")[[1]][1];
	}
	if(length(unique(new))!=length(old)){
		cat("Error: Can not truncate category names to ID, because the IDs are not unique.\n");
		quit(-1);
	}
	colnames(combined_count_table)=new;
}

# Collapse counts by group
grouped_count_table=combine_summary_file_by_group(combined_count_table, file_list);

###############################################################################

write_summary_file=function(out_mat, fname){
	cat("Writing summary table: ", fname, "\n");
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

###############################################################################

# Output Individual
write_summary_file(combined_count_table, paste(OutputFileNameRoot, ".by_file.summary_table.tsv", sep=""));
write_summary_file(grouped_count_table, paste(OutputFileNameRoot, ".by_group.summary_table.tsv", sep=""));

###############################################################################

cat("Done.\n")
warns=warnings();
if(length(warns)>0){
	print(warnings());
}

q(status=0)
