
cat("Loading Wilcoxon Cluster Difference...\n");

###############################################################################

wilcoxon_cluster_difference=function(a_idc, b_idc, distmat){
	# a_idc: indices of members in cluster A
	# b_idc: indices of members in cluster B
	# distmat: distance matrix. 

	num_a=length(a_idc);
	num_b=length(b_idc);
	num_total_samples=num_a+num_b;
	samples_to_analyze=c(a_idc, b_idc);

	if(class(distmat)=="dist"){
		distmat=as.matrix(distmat);
	}

	if(num_a==0){
		cat("Error: Num sequences in A is 0\n");
		return(NULL);
	}
	if(num_b==0){
		cat("Error: Num sequences in B is 0\n");
		return(NULL);
	}

	# Store the average distance between each sample and the cluster its in, and the cluster that its not in.
	intra_distances=numeric(num_total_samples);
	inter_distances=numeric(num_total_samples);

	count=1;  # samples_to_analyze is not sequential so we can't use that as an index
	for(i in samples_to_analyze){

		if(any(i==a_idc)){
			intra=distmat[i, a_idc];
			inter=distmat[i, b_idc];
		}else if(any(i==b_idc)){
			intra=distmat[i, b_idc];
			inter=distmat[i, a_idc];
		}else{
			cat("Error: sample index not in either cluster.\n");
		}

		#cat("Num Intra: ", length(intra), "\n");
		#cat("Num Inter: ", length(inter), "\n");
		
		intra_distances[count]=sum(intra)/(length(intra)-1); # get rid of 0 from self vs self.
		inter_distances[count]=mean(inter);
		count=count+1;

	}

	#cat("Distances:\n");
	#print(intra_distances);
	#cat("vs\n");
	#print(inter_distances);

	# Run wilcoxon rank sum test (Mann-Whitney)
	wilcox_result=wilcox.test(intra_distances, inter_distances, paired=FALSE, alternative="two.sided");
	
	# Return a list
	results=list();
	results$wilcoxon=wilcox_result;
	results$intra_distances=intra_distances;
	results$inter_distances=inter_distances;

	return(results);	
}

###############################################################################

plot_mds=function(a_idc, b_idc, distmat, title=""){

	# If distmat is half matrix, then make a full matrix out of it
	if(class(distmat)=="dist"){
		distmat=as.matrix(distmat);
	}

	# Disalllow 0's off the diagonal, self vs self
	num_samples=ncol(distmat);
	for(i in 1:num_samples){
		for(j in 1:num_samples){
			if(distmat[i,j]==0 && i!=j){
				distmat[i,j]=1e-323;
			}
		}
	}

	# Color a and b and blacken others
	colors=rep("black", num_samples);
	colors[a_idc]="green";
	colors[b_idc]="blue";

	# Size a and b, and downsize others
	text_scale=42/num_samples;
	if(text_scale>1){
		text_scale=1;
	}
	sizes=rep(.8, num_samples);
	sizes[c(a_idc,b_idc)]=1.2;
	sizes=sizes*text_scale;

	# Get sample names
	sample_names=colnames(distmat);	

	# Compute and plot MDS
	iso=isoMDS(distmat);
	plot(iso$points[,1], iso$points[,2], col=colors, cex=sizes, xlab="", ylab="", main=title);
	text(iso$points[,1], iso$points[,2], labels=sample_names, col=colors, cex=sizes, pos=3);

}

###############################################################################

plot_combo_histogram=function(a,b, notes=c()){
# Plots two histograms in the same graph.

	# Estimate range of combine histogram first
	combine=c(a,b);
	comb_hist_info=hist(combine, plot=FALSE);

	# Use bin sizes of combined histogram to bin individual 
	ahist=hist(a, breaks=comb_hist_info$breaks, plot=FALSE);
	bhist=hist(b, breaks=comb_hist_info$breaks, plot=FALSE);

	# Do a side-by-side bar plot plot each histogram separately
	barplot(rbind(ahist$counts, bhist$counts), 
		beside=TRUE, names.arg=sprintf("%3.2f",comb_hist_info$mids),
		main="Intra (Dark) vs Inter (Light) Distances",
		ylab="Count",
		xlab="Mean Distance"
	);

	# Add notes to top of plot
	l=0
	for(str in notes){
		mtext(str, side=3, line=l);
		l=l-1
	}
}

###############################################################################

cat("Loaded.\n");

if(0){

	library(MASS);

	distances1=c(
		0,1,1,1,2,3,2,3,
		1,0,1,1,1,2,2,2,
		1,1,0,1,2,3,2,3,
		1,1,1,0,1,2,1,2,
		2,1,2,1,0,1,1,1,
		3,2,3,2,1,0,1,1,
		2,2,2,2,1,1,0,1,
		3,2,3,2,1,1,1,0
	);

	distances2=c(
		0,2,2,3,1,3,4,4,
		2,0,3,3,2,1,4,4,
		2,3,0,3,2,3,1,4,
		3,3,3,0,3,2,2,1,
		1,2,2,3,0,2,3,3,
		3,1,3,2,2,0,3,2,
		4,4,1,2,3,3,0,2,
		4,4,4,1,3,2,2,0
	);


	matrix=matrix(distances2,nrow=8);
	names=1:8;
	rownames(matrix)=names;
	colnames(matrix)=names;
	dist=as.dist(matrix, diag=TRUE, upper=TRUE);
	print(dist);
	
	a=c(3,4,8)
	b=c(1,2,5,6);


	result=wilcoxon_cluster_difference(a,b, dist);
	print(result);

	par(mfrow=c(2,1));
	plot_mds(a,b, dist);
	plot_combo_histogram(result$intra_distances, result$inter_distances, notes=c("test", "test2"));

}

