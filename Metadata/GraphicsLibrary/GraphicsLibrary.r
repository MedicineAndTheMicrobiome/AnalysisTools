
draw_chevron_line=function(start, end, cex=1, col="black", arrow_dens=1){
	# This function draws a line from the start to end with >>>>>
	# The start and end arguments are vectors, e.g.  c(x, y)
	# Since the directionality of the line is not implied by the
	#    end of the line with an arrow, the end coordinates do
	#    not need to have precise coordinates.
	# You can draw between the centers of objects, and if the objects
	#    are opaque, the object edges will appear to be connect.

	cat("Draw chevron line from: (", 
		start[1], ", ", start[2], ") to (", 
		end[1], ", ", end[2], ")\n", sep="");

        vec=end-start;
        distance=sqrt(sum(vec^2));
	cat("Dist: ", distance, "\n");

        num_arrows=40*distance;
        arrow_spacing=seq(0,1, length.out=num_arrows);
        arrow_pos=t(sapply(arrow_spacing, function(i){ return(start + i*vec);}));

        for(i in 1:num_arrows){
                text(arrow_pos[i,1], arrow_pos[i,2], ">",
                        srt=atan2(vec[2], vec[1])*180/pi,
                        cex=cex*1.25,
                        col=col,
                        adj=0.5);
        }

}

