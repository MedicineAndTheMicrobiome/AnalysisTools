#!/usr/bin/env Rscript

###############################################################################

library(getopt);

options(useFancyQuotes=F);
options(width=120);

params=c(
	"input", "i", 1, "character",
	"parent_id_fn", "p", 1, "character",
	"output_dir", "d", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input OBO Table/Mapping .tsv file>\n",
	"	-p <input target parent GO IDs filename>\n",
	"	-o <output directory>\n",
	"\n",
	"This script will:\n",
	"  1.) Read in the parent GO IDs from the target parent file.\n",
	"  2.) For each GO parent ID, look up it's immediate children. \n",
	"  2.) For each immediate child, all it's decendants will be identified.\n",
	"\n",
	"A set of files will be created.\n",
	"  1.) GO<GO ID>.chld.descrp.map\n",
	"         Format: <GO ID>\\t<textual description\\n\n",
	"  2.) GO<GO ID>.chld.descnd.map\n",
	"         Format: <GO ID>\\t<descendent1>\\n\n",
	"                 <GO ID>\\t<descendent2>\\n\n",
	"                 <GO ID>\\t<descendent...>\\n\n",
	"                 <GO ID>\\t<descendentn>\\n\n",
	"\n",
	"	 A single child of the specified parent ID will have multiple descendants.\n",
	"\n",
	"All of these files will be placed into the specified output directory.\n",
	"\n");

if(
	!length(opt$input) |
	!length(opt$parent_id_fn) |
	!length(opt$output_dir)
){
	cat(usage);
	q(status=-1);
}

InputOBOMap=opt$input;
ParentIDFilename=opt$parent_id_fn;
OutputDir=opt$output_dir;

##############################################################################


cat("\n");
cat("Targeted Parent ID Filename: ", ParentIDFilename, "\n");
cat("OBO Map Name: ", InputOBOMap, "\n");
cat("Output Dir: ", OutputDir, "\n");
cat("\n");

if(!dir.exists(OutputDir)){
	cat("Creating: ", OutputDir, "\n");
	dir.create(OutputDir);
}else{
	cat(OutputDir, " already exists.  Using it.\n");
}

##############################################################################

load_list=function(fn){
	dat=read.table(fn, header=F, sep="\t", quote=NULL, comment.char="#");
	return(dat[,1]);
}

##############################################################################

load_ids_to_names=function(obo_tab){

	cat("Loading ID to Name Map...\n");
	start_time=Sys.time();
	id_to_name_map=list();
	num_rows=nrow(obo_tab);
	for(i in 1:num_rows){
		id_to_name_map[[obo_tab[i, "id"]]]=obo_tab[i, "name"];
	}
	end_time=Sys.time();
	cat("ok.\n");
	exec_time=end_time-start_time;
	cat("Time: ", exec_time, " secs\n", sep="");
	return(id_to_name_map);
}


##############################################################################

lookup_name=function(id, id_name_map){

	#cat("Looking up: \n");
	#print(id);
	num_ids=length(id);
	out_arr=character(num_ids);
	for(i in 1:num_ids){
		out_arr[i]=id_name_map[[id[i]]];
	}
	return(out_arr);
}

##############################################################################

load_obo_to_tree=function(obo_tab){
	
	cat("Loading OBO Table into Tree...\n");
	start_time=Sys.time();

	num_ids=nrow(obo_tab);
	parent_to_child_tree=vector("list", length=num_ids);
	names(parent_to_child_tree)=obo_tab[,"id"];

	parent_presplit=strsplit(obo_tab[,"is_a"], ";");

	for(i in 1:num_ids){
		child_id=obo_tab[i,"id"];
		parent_ids=parent_presplit[[i]];
		for(pid in parent_ids){
			parent_to_child_tree[[pid]]=c(
				parent_to_child_tree[[pid]],
				child_id);
		}
	}

	end_time=Sys.time();
	cat("ok.\n");
	exec_time=end_time-start_time;
	cat("Time: ", exec_time, " secs\n", sep="")
	return(parent_to_child_tree);

}

##############################################################################

find_all_descendants=function(parent_id, parent_to_child_tree){

	# For a particular parent_id, all the descendants will be returned.
	# The parent is not included.

	children=parent_to_child_tree[[parent_id]];	

	if(!length(children)){
		# Base case, when there is not children (leaf node)
		return(c());
	}else{
		all_descendants=children;
		for(child in children){
			child_descendants=find_all_descendants(child, parent_to_child_tree);
			all_descendants=c(all_descendants, child_descendants);
		}
		return(all_descendants);
	}
}

##############################################################################

find_immediate_children=function(parent_id, parent_to_child_tree){
	return(parent_to_child_tree[[parent_id]]);
}

##############################################################################

write_list_to_mapfile=function(fname, parent_id, child_list, idnm_map){

	fh=file(fname, "w");

	cat(file=fh, parent_id, "\t", parent_id, "\t",
		idnm_map[[parent_id]], "\n", sep="");

	immed_child_ids=names(child_list);
	
	for(immed_child_idx in immed_child_ids){
		cat(file=fh, immed_child_idx, "\t", immed_child_idx, "\t", 
			idnm_map[[immed_child_idx]], "\n", sep="");

		desc_arr=child_list[[immed_child_idx]];
		
		for(descend_ix in desc_arr){
			cat(file=fh, immed_child_idx, "\t", descend_ix, "\t",
			idnm_map[[descend_ix]], "\n", sep="");
		}
	}

	close(fh);
}

##############################################################################

load_obo_map=function(obo_map_fn){

	cat("Loading OBO Map File...\n");

	start_time=Sys.time();

	mat=read.table(InputOBOMap, quote="", comment.char="", sep="\t", 
		row.names=NULL, header=T);

	end_time=Sys.time();
	exec_time=end_time-start_time;

	cat("Time: ", exec_time, " secs\n", sep="");
	cat("Num Rows: ", nrow(mat), " Loaded.\n");
	cat("OBO Map Column Names:\n");
	print(colnames(mat));
	cat("\n");

	return(mat);
}

##############################################################################

write_descriptions=function(fname, grps, lookup_map){

	cat("Writing Descriptions: ", fname, "\n", sep="");
	grp_names=names(grps);
	num_names=length(grp_names);

	colnm=c("GO_ID", "Name", "[Num_Members]");
	name_matrix=matrix("", nrow=num_names, ncol=length(colnm));
	colnames(name_matrix)=colnm;
	
	for(i in 1:num_names){
		num_members=length(grps[[grp_names[i]]]);
		name_matrix[i,]=c(
			grp_names[i], 
			lookup_name(grp_names[i], lookup_map), 
			paste("[", num_members, "]", sep=""));
	}

	write.table(name_matrix, fname, row.names=F, col.names=T, quote=F);
	cat("ok.\n");
	
}

##############################################################################

write_parent_child_map=function(fname, grps, lookup_map){

	cat("Writing Groupings Map: ", fname, "\n", sep="");

	grp_names=names(grps);
	num_names=length(grp_names);

	fh=file(fname, "w");
	cat(file=fh, paste(c("GroupID", "Member", "Description"), collapse="\t"), "\n", sep="");
	close(fh);

	for(nm in grp_names){

		ids=grps[[nm]];
		num_ids=length(ids);
		out_mat=matrix("", nrow=num_ids, ncol=3);
		out_mat[,1]=rep(nm, num_ids);
		out_mat[,2]=ids;
		out_mat[,3]=lookup_name(ids, lookup_map);
	
		cat("\tWriting descendants for: ", nm, "\n");
		write.table(out_mat, fname, row.names=F, col.names=F, quote=F, append=T);
	}
	cat("ok.\n");

}

##############################################################################

extract_descendants=function(pid, parent_to_child_tree, lookup_map){


}

##############################################################################
# Load inputs

target_parent_ids=load_list(ParentIDFilename);
num_target_parent_ids=length(target_parent_ids);

cat("Target Parent IDs [", num_target_parent_ids, "]:\n");
print(target_parent_ids);
cat("\n");

obo_mat=load_obo_map(InputOBOMap);
id_to_name_map=load_ids_to_names(obo_mat);
go_tree=load_obo_to_tree(obo_mat);

for(pidx in target_parent_ids){

	cat("-------------------------------------------------------------\n");

	parent_id_name=lookup_name(pidx, id_to_name_map);
	cat("\n");
	cat("Name of Parent ID: ", pidx, ":\n", sep="");
	cat("\t\"", parent_id_name, "\"\n", sep="");
	cat("\n");

	# Lookup immediate children
	immed_children_arr=find_immediate_children(pidx, go_tree);
	num_immed_children=length(immed_children_arr);

	cat("Num Immediate Children: ", num_immed_children, "\n");
	print(immed_children_arr);
	cat("\n");

	# Foreach immediate child look up descendants.
	# Include parents for a placeholder
	groupings=list();
	groupings[[pidx]]=pidx;

	for(imm_child in immed_children_arr){

		cat("Looking up Immediate Child: ", imm_child, "\n");	
		start_time=Sys.time();

		# Include immediate child for placeholder
		descendants=find_all_descendants(imm_child, go_tree);
		end_time=Sys.time();
		search_time=end_time-start_time;

		num_desc=length(descendants);
		cat("  Num of Descendants: ", num_desc, "\n", sep="");
		cat("  Search time: ", search_time, " secs.\n", sep="");

		groupings[[imm_child]]=c(imm_child, descendants);

	}

	##############################################################################

	cat("\n");

	cln_pid=gsub(":", "_", pidx);

	# Write descriptions
	child_names_fn=paste(OutputDir, "/", cln_pid, ".immediate_descriptions.tsv", sep="");
	write_descriptions(child_names_fn, groupings, id_to_name_map);

	# Write mappings
	imm_child_to_desc_map_fn=paste(OutputDir, "/", cln_pid, ".groupings.map", sep="");
	write_parent_child_map(imm_child_to_desc_map_fn, groupings, id_to_name_map);

	##############################################################################

	cat("-------------------------------------------------------------\n");
}


cat("\nDone.\n");

print(warnings());
q(status=0);
