
###############################################################################

load_distance_matrix=function(distmat_name){
# Loads a distance matrix
	cat("Loading Distance Matrix: ", distmat_name, "\n");
        distance_matrix=as.matrix(read.table(distmat_name, header=TRUE, check.names=FALSE));
        dims=dim(distance_matrix);
        cat("Distance Matrix: Rows: ", dims[1], " x Cols: ", dims[2], "\n");
        return(distance_matrix);
}

###############################################################################

load_map=function(filename, numeric=FALSE){
# Loads a two column matrix so that the first column can be translated into the second column
        map=list();

	if(filename!=""){
		cat("Loading Mapping File: ", filename, "\n");
		members=as.matrix(read.table(filename, sep="\t", header=FALSE));
		num_members=nrow(members);
		cat("Num members in map: ", num_members, "\n");

		if(numeric){
			for(i in 1:num_members){
				map[[members[i,1]]]=as.numeric(members[i,2]);
			}
		}else{
			for(i in 1:num_members){
				map[[members[i,1]]]=members[i,2];
			}
		}
	}

        return(map);
}

###############################################################################

rename_distance_matrix=function(distmat, new_name_map=list()){
# This function will rename the rows/column of the distance matrix based
# on the new_name_map.  Note that the names in the distance matrix must be
# unique.  So, if there are any degenerate oldname to newname maps, the
# new name, will be appended with the old accession.

	# Get map information
	oldnames=names(new_name_map);             # Essently the hash keys
	num_map_entries=length(new_name_map);
	cat("Num map entries: ", num_map_entries, "\n");

	if(num_map_entries>0){

		# Generate vector of new names
		distmat_names=colnames(distmat);
		num_dm_names=length(distmat_names);
		new_names=character(num_dm_names);

		for(i in 1:num_dm_names){
			dm_name=distmat_names[i];		# Eg: ABW23353.1
			new_name=new_name_map[[dm_name]];	# Eg: A/BRISBANE/10/2007
			#cat("Renaming: ", dm_name, " -> ", new_name, "\n");

			if(!is.null(new_name)){
				new_names[i]=paste(new_name, ":", dm_name,  sep="");
				# So it will look like A/BRISBANE/10/2007:ABW23353.1
			}else{
				new_names[i]=dm_name;
				# Do not change any names
			}
		}

		# Save names back into distance matrix
		rownames(distmat)=new_names;
		colnames(distmat)=new_names;
	}

	# Return updated distance matrix
	return(distmat);
}

###############################################################################

rename_weighting_map=function(weighting_map,  new_name_map=list()){

	new_weighting_map=list();
	wt_map_length=length(weighting_map);
	keys=names(weighting_map);

	for(i in 1:wt_map_length){

		old_key=keys[i];

		new_name=new_name_map[[old_key]];
		if(!is.null(new_name)){
			new_key=paste(new_name, ":", old_key, sep="");
		}

		new_weighting_map[[new_key]]=weighting_map[[old_key]];
	}
	
	return(new_weighting_map);
}

###############################################################################

cat("Sequence-Distancing Libraries Loaded...\n");

if(0){
	distmat=load_distance_matrix("short_dist");
	print(distmat)

	map=list();
	map[["ACS71642.1"]]="Apples";
	map[["ACI26318.1"]]="Apples";
	map[["EPI125925"]]="Bananas";

	newdistmat=rename_distance_matrix(distmat, map);
	print(newdistmat)
}

