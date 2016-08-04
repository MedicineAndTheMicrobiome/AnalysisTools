#!/usr/bin/env Rscript

cat("Loading Wilcoxon ANOVA.\n");

###############################################################################

wilcoxanova=function(dist, grpA_idx, grpB_idx){
	
	# Get number of samples in each group for sanity check
	nA=length(grpA_idx);
	nB=length(grpB_idx);
	cat("GrpA Size: ", nA, "\n");
	cat("GrpB Size: ", nB, "\n");

	#print(grpA_idx);
	#print(grpB_idx);

	# Convert dist structure (half matrix) to full matrix
	distmat=as.matrix(dist);
	
	# Extract the INTER distances
	inter_dist=distmat[grpA_idx, grpB_idx];
	#print(inter_dist);
	#print(as.vector(inter_dist));

	# Extract the INTRA distances
	intra_dist_A=distmat[grpA_idx, grpA_idx];
	intra_dist_B=distmat[grpB_idx, grpB_idx];

	# Just grab the bottom triangle, excluding self-self
	intra_dist_A_vect=as.vector(as.dist(intra_dist_A));
	intra_dist_B_vect=as.vector(as.dist(intra_dist_B));
	intra_dist=c(intra_dist_A_vect, intra_dist_B_vect);

	# Compare median intra to intra distances.
	# Test if INTER is greater than INTRA.
	# hist(c(inter_dist, intra_dist))
	wctest=wilcox.test(inter_dist, intra_dist, alternative="greater");

	# Prepare return value
	res=list();
	res$inter_med=median(inter_dist);
	res$intra_med=median(intra_dist);
	res$pval=wctest$p.value;
	res$grpA_med=median(intra_dist_A_vect);
	res$grpB_med=median(intra_dist_B_vect);
	return(res);
}

###############################################################################

# Set to 1 to test
if(0){

	#set.seed(10);
	n1=100;
	n2=100;
	cl1x=rnorm(n1, 0, 1);
	cl1y=rnorm(n1, 0, 1);
	cl2x=rnorm(n2, 1, 1);
	cl2y=rnorm(n2, 1, 1);

	x_all=c(cl1x, cl2x);
	y_all=c(cl1y, cl2y);

	mat_all=cbind(x_all, y_all);
	rownames(mat_all)=c(sprintf("A%i", 1:n1), sprintf("B%i", 1:n2));
	#print(mat_all);
	
	dist_mat=dist(mat_all);
	#print(dist_mat);
	plot(x_all, y_all, col=c(rep(1, n1), rep(2, n2)));

	res=wilcoxanova(dist_mat, 1:n1, (n1+1):(n1+n2));
	
	print(res);
	
}

###############################################################################

add_centroids_to_distmat=function(dist, grpA_idx, grpB_idx){
	
	distmat=as.matrix(dist);
	
	#print(distmat);
	cA=apply(distmat[grpA_idx,],2,mean);
	cB=apply(distmat[grpB_idx,],2,mean);
	#print(cA);
	#print(cB);

	cvc_dist=mean(distmat[grpA_idx,grpB_idx]);
	newdistmat=rbind(distmat,cA,cB);
	newdistmat=cbind(newdistmat, c(cA,0,cvc_dist), c(cB,cvc_dist,0));
	colnames(newdistmat)=rownames(newdistmat);

	#print(newdistmat);
	return(as.dist(newdistmat))
}

if(0){
	set.seed(10);
	pts=matrix(runif(15,0,10), nrow=5);
	distmat=dist(pts);

	print(distmat);
	wcent=add_centroids_to_distmat(distmat, 1:2, 3:5);
	print(wcent);
}

###############################################################################

level_rotate=function(pt1, pt2, points){

	x1=pt1[1];
	y1=pt1[2];
	x2=pt2[1];
	y2=pt2[2];
	# Compute angle between pt1 and pt2
	rad_ang=atan((y2-y1)/(x2-x1));
	print(rad_ang*(360/(2*pi)))
	rot_radang=rad_ang;
	if(x1<x2){
		rot_radang=rot_radang+pi;
	}	
	cat("Counter Rotation Angle: ", rot_radang*360/(2*pi), "\n")
	rot_mat=matrix(c(cos(rot_radang), -sin(rot_radang), sin(rot_radang), cos(rot_radang)), nrow=2, byrow=T);
	print(rot_mat);
	rotated=points %*% rot_mat;

	return(rotated);
	
}

if(1){

	pdf("Rplots.pdf", height=8.5, width=17);
	par(mfrow=c(1,2));
	points=matrix(runif(20,-1,1)*7, ncol=2);
	plot(points[,1], points[,2], col=1:20, ylim=c(-10,10), xlim=c(-10,10), type="l");
	rotated=level_rotate(c(0,0), c(1,0), points);
	plot(rotated[,1], rotated[,2], col=1:20, ylim=c(-10,10), xlim=c(-10,10), type="l");
	

}

###############################################################################
