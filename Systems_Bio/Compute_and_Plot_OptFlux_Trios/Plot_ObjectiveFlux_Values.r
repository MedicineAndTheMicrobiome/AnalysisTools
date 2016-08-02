#!/usr/bin/env Rscript

###############################################################################

library('getopt');
library('plot3D');

options(useFancyQuotes=F);

params=c(
		"input", "i", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

PTS_PD=15;

usage = paste (
		"\nUsage:\n", script_name, "\n",
		"	-i <Input coordinates/values file>\n",
		"\n",
		"This script will read in the coordinates/values file from\n",
		"runing the compute intensive Compute_ObjectiveFlux_Trios.r\n",
		"script and generate graphics.\n",
		"\n");

if(!length(opt$input)){
	cat(usage);
	q(status=-1);
}

InputFile=opt$input;
OutputFile=gsub("\\.tsv", "", InputFile);

###############################################################################

plot_flux_ranges=function(coord, values){

	coord_name=colnames(coord);	

	# Find max slices
	max_val=max(values);
	max_ix=min(which(max_val==values));
	coords_at_max=coord[max_ix,];
	cat("Coords at max: ", paste(coords_at_max, collapse=","), "\n");
	cat("Val at max:  ", max_val, "\n");
	cat("\n");

	# Assign colors
	color_res=20;
	colors=rev(rainbow(color_res, start=0, end=4/6));	
	val_range=range(c(values,0));
	cat("Value Ranges: ", val_range[1], " - ", val_range[2], "\n");
	val_span=val_range[2]-val_range[1];
	scaled_values=(values-val_range[1])/val_span;
	color_ix= as.integer((color_res-1)*scaled_values)+1;
	assigned_color=colors[color_ix];

	almost_equal=function(x, y){
		return(abs(x-y)<1e-6);
	}

	# Extract slice
	yz_slice_ix=which(almost_equal(coords_at_max[1], coord[,1]));
	xz_slice_ix=which(almost_equal(coords_at_max[2], coord[,2]));
	yx_slice_ix=which(almost_equal(coords_at_max[3], coord[,3]));

	par(oma=c(0,0,3,0));	

	# Plot slices
	par(mfrow=c(2,2));
	#par(mar=c(4,4,2,1));
	lay_mat=matrix(c(
		1,1,1, 2,2,2,
		1,1,1, 2,2,2,
		1,1,1, 2,2,2,
	
		3,3,3, 5,5,5,
		3,3,3, 4,4,4,
		3,3,3, 6,6,6
	), byrow=T, ncol=6);
	layout(lay_mat);

	trans_red=rgb(.7,.7,.7, .7);
	trans_green=rgb(.7,.7,.7, .7);	
	trans_blue=rgb(.7,.7,.7,.7);

	pt_size=2;
	pt_shp=16;

	plot(coord[yz_slice_ix, 2], coord[yz_slice_ix, 3], xlab=coord_name[2], ylab=coord_name[3], type="n");
	rect(-1000, -1000, 0, 1000, col=trans_blue, border=NA);
	rect(-1000, -1000, 1000, 0, col=trans_green, border=NA);
	abline(h=0, lty="dotted");
	abline(v=0, lty="dotted");
	points(coord[yz_slice_ix, 2], coord[yz_slice_ix, 3], col=assigned_color[yz_slice_ix], pch=pt_shp, cex=pt_size);
	abline(h=0, lty="dotted", col="grey");
	abline(v=0, lty="dotted", col="grey");

	plot(coord[xz_slice_ix, 1], coord[xz_slice_ix, 3], xlab=coord_name[1], ylab=coord_name[3], type="n");
	rect(-1000, -1000, 0, 1000, col=trans_red, border=NA);
	rect(-1000, -1000, 1000, 0, col=trans_green, border=NA);
	abline(h=0, lty="dotted");
	abline(v=0, lty="dotted");
	points(coord[xz_slice_ix, 1], coord[xz_slice_ix, 3], col=assigned_color[xz_slice_ix], pch=pt_shp, cex=pt_size);
	abline(h=0, lty="dotted", col="grey");
	abline(v=0, lty="dotted", col="grey");

	plot(coord[yx_slice_ix, 2], coord[yx_slice_ix, 1], xlab=coord_name[2], ylab=coord_name[1], type="n");
	rect(-1000, -1000, 0, 1000, col=trans_blue, border=NA);
	rect(-1000, -1000, 1000, 0, col=trans_green, border=NA);
	abline(h=0, lty="dotted");
	abline(v=0, lty="dotted");
	points(coord[yx_slice_ix, 2], coord[yx_slice_ix, 1], col=assigned_color[yx_slice_ix], pch=pt_shp, cex=pt_size);
	abline(h=0, lty="dotted", col="grey");
	abline(v=0, lty="dotted", col="grey");


	#plot(coord[xz_slice_ix, 1], coord[xz_slice_ix, 3], col=assigned_color[xz_slice_ix],
	# 	pch=15,cex=3, xlab=coord_name[1], ylab=coord_name[3]);
	#plot(coord[yx_slice_ix, 2], coord[yx_slice_ix, 1], col=assigned_color[yx_slice_ix],
	# 	pch=15,cex=3, xlab=coord_name[2], ylab=coord_name[1]);

	barplot(rep(1,color_res), col=colors, main="Flux Levels at Optimal Solution",
		names.arg=sprintf("%6.5f", seq(0, max_val, length.out=color_res)),
		las=2, yaxt="n"
	);

	plot(0,0, type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n");
	text(0,0, sprintf("Maximum Flux: %6.5f", max_val), cex=2.5);

	mtext(side=3, outer=T, OutputFile, cex=2, font=2);
	
	for(mode in c("all", "opt_slice")){
		
		# Identify non-zero coordinates
		nz=(values>1e-6);

		#phi=rotates around x axis (away/towards), theta=rotates around z axis (away/towards)
		par(mfrow=c(4,3));
		par(mar=c(1,1,1,1));

		if(mode=="all"){
			draw_ix=nz;
		}else if(mode=="opt_slice"){
			draw_ix=intersect(which(nz), c(yz_slice_ix, xz_slice_ix, yx_slice_ix));
		}

		for(thetav in seq(0,360, length.out=12)){
			scatter3D(coord[draw_ix,1], coord[draw_ix,2], coord[draw_ix,3], ticktype="detailed", 
				pch=21, col="black", bg=assigned_color[draw_ix],
				xlab=coord_name[1], ylab=coord_name[2], zlab=coord_name[3],
				colvar=color_ix[draw_ix], theta=thetav, phi=45, colkey=list(plot=F));
		}
		mtext(side=3, outer=T, OutputFile, cex=2, font=2);

		for(phiv in seq(0,360, length.out=12)){
			scatter3D(coord[draw_ix,1], coord[draw_ix,2], coord[draw_ix,3], ticktype="detailed", 
				pch=21, col="black", bg=assigned_color[draw_ix],
				xlab=coord_name[1], ylab=coord_name[2], zlab=coord_name[3],
				colvar=color_ix[draw_ix], theta=45, phi=phiv, colkey=list(plot=F));
		}
		mtext(side=3, outer=T, OutputFile, cex=2, font=2);

	}

}

data=as.matrix(read.table(InputFile, header=T, sep="\t"));
#print(data);

pdf(paste(OutputFile, ".pdf", sep=""), height=11, width=8.5);

cat("\n");
cat("Example Coordinates: \n");
print(head(data[,1:3]));
cat("\n");
cat("Example Fluxes: \n");
print(head(data[,4]));
cat("\n");
plot_flux_ranges(data[,1:3], data[,4]);







###############################################################################

print(warnings());
cat("Done.\n");


























