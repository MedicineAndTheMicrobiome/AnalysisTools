
###############################################################################

load_antigenic_distances=function(distance_file){

	# Load the antigenic distances
	cat("Loading Antigenic Distances: ", distance_file, "\n");
	distance_list=as.matrix(read.table(distance_file, header=F, check.names=F));
	dims=dim(distance_list);

	# Double check the number of columns make sense
	if(!(dims[2]==4 || dims[2]==3)){
		cat("Error: Unexpected columns count.\n"); 
	}

	# If the 4 column doesn't contain standard deviations, assume they are 0's.
	if(dims[2]==3){
		stdevs=rep(0, dims[1]);
	}else{
		stdevs=distance_list[,4];
	}

	# Label columns
	pairwise_dist=list();
	pairwise_dist[["AntigenA"]]=distance_list[,1];
	pairwise_dist[["AntigenB"]]=distance_list[,2];
	pairwise_dist[["Distance"]]=as.numeric(distance_list[,3]);
	pairwise_dist[["StdDev"]]  =as.numeric(stdevs);

	#print(pairwise_dist);

	return(pairwise_dist);
}

###############################################################################

load_list=function(list_file){
	cat("Loading List: ", list_file, "\n");
        list=as.matrix(read.table(list_file, header=F, check.names=F));
	return(as.vector(list));
}

###############################################################################

get_list_of_all_antigens=function(pairwise_dists){
	unique_antigens=unique(c(pairwise_dists[["AntigenA"]], pairwise_dists[["AntigenB"]]));
	cat("Num antigens: ", length(unique_antigens), "\n");
	return(unique_antigens);
}

###############################################################################

extract_distances_by_antigen=function(antigen, pairwise_dists){

	cat("Extracting by: ", antigen, "\n");

	#print(pairwise_dists);
	#for(i in 1:length(pairwise_dists$Distance)){
	#	cat(pairwise_dists$AntigenA[i], " ", pairwise_dists$AntigenB[i], " ", pairwise_dists$Distance[i], "\n");
	#}

	A_idx=which(pairwise_dists[["AntigenA"]]==antigen);
	B_idx=which(pairwise_dists[["AntigenB"]]==antigen);
	combined_idx=c(A_idx, B_idx);

	# If To matches, return From.  If From matches, return To
	Aantigens=pairwise_dists[["AntigenB"]][A_idx];
	Bantigens=pairwise_dists[["AntigenA"]][B_idx];

	# Build structure, with antigens of interest
	oneway_dists=list();
	oneway_dists[["FromAntigen"]]=antigen;
	oneway_dists[["ToAntigen"]]=c(Aantigens, Bantigens, antigen);
	oneway_dists[["Distance"]]=c(pairwise_dists[["Distance"]][combined_idx], 0);
	oneway_dists[["StdDev"]]=c(pairwise_dists[["StdDev"]][combined_idx], 0);

	#for(i in 1:length(oneway_dists$Distance)){
	#	cat(oneway_dists$ToAntigen[i], " ", oneway_dists$Distance[i], "\n");
	#}

	return(oneway_dists);
}

###############################################################################

perturb_distances=function(oneway_dists){

	# Resample with replacement from distances that are available
	num_distances=length(oneway_dists[["Distance"]]);
	resampled_dist_idx=sample(1:num_distances, num_distances, replace=T);
	
	# Perturb the distances based on the specified standard deviation
	stdevs=oneway_dists[["StdDev"]];
	distances=oneway_dists[["Distance"]];
	perturbed_distances=numeric(num_distances);
	for(i in 1:num_distances){
		idx=resampled_dist_idx[i];
		perturbed_distances[i] = rnorm(1, distances[idx], stdevs[idx]);
	}

	# Rebuild structure with perturbed antigens and distances
	#   We need to attach the reference antigen back to the list because we always know
	#   the distance to itself is 0.
	perturbed_oneway_dist=list();
	perturbed_oneway_dist[["FromAntigen"]]=oneway_dists[["FromAntigen"]];
	perturbed_oneway_dist[["ToAntigen"]]=c(oneway_dists[["ToAntigen"]][resampled_dist_idx], oneway_dists[["FromAntigen"]]);
	perturbed_oneway_dist[["Distance"]]=c(perturbed_distances, 0);
	perturbed_oneway_dist[["StdDev"]]=-1;	# Means perturbation already applied

	# Return new distances
	return(perturbed_oneway_dist);
}

###############################################################################

get_antigens_within_cutoff=function(cutoff, oneway_dists){

	num_distances=length(oneway_dists[["Distance"]]);
	distances=oneway_dists[["Distance"]];
	toAntigenNames=oneway_dists[["ToAntigen"]];

	# Collect all the distances by antigen name
	grouped_distances=list();
	for(i in 1:num_distances){
		toAntigenName=toAntigenNames[i];
		dist=distances[i];
		if(length(grouped_distances[[toAntigenName]])>0){
			grouped_distances[[toAntigenName]]=c(grouped_distances[[toAntigenName]], dist);
		}else{
			grouped_distances[[toAntigenName]]=dist;
		}
	}

	num_collapsed_dist=length(grouped_distances);
	collapsed_names=names(grouped_distances);

	inout_map=list();
	for(i in 1:num_collapsed_dist){

		# Average distances together if there is more than one
		ag_name=collapsed_names[i];
		mean_dist=mean(grouped_distances[[ag_name]]);

		# Determine IN or OUT
		if(mean_dist<=cutoff){
			inout_map[[ag_name]]="IN";
		}else{
			inout_map[[ag_name]]="OUT";
		}

	}

	return(inout_map);
}

###############################################################################

cat("Antigenic Distancing Library Loaded...\n");

if(0){
	a=c("apple", "banana", "orange", "strawberry");
	b=c("dairy", "elephant", "frenchfry", "goat");
	d=c(10 ,30, 40, 50);
	s=c(0, 20, 5, 1);

	sample_distances=data.frame(a,b,d,s);

	pd=perturb_distances(sample_distances);
	print(pd);
}

