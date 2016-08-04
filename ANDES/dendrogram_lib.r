
###############################################################################

infer_leaf_coverage=function(group_map, dend, reference){
# This is the DASH algorithm.
#
# Given a dendrogram constructed based on a sequence-based distance matrix
#  mark the leaves (antigens) that are specified in the group map according to
#  whether they are "in" or "out".  Identify whether the edges coming
#  out of the leaves are "in" or "out".  If there is disagreement, then mark the
#  edge as unknown.  Then from top down, if there is consensus between the
#  two children of each node, then mark all the children as "in" or "out".

	color_map=list();
	color_map$IN="green";
	color_map$OUT="red";
	color_map$UNKNOWN="grey";
	color_map$CONFLICT="black";
	par(mfrow=c(1,1));
	par(mar=c(16,4,4,2));
	
	# Mark leaf nodes with group identifier
	#cat("Marking leaf nodes with group identifier.\n");
	dend=dendrapply(dend, mark_leaf_nodes_with_group_id, group_map, reference);
	#cat("done.\n");

	#dend=dendrapply(dend, color_dendrogram_by_group_id, color_map);
	#plot(dend);

	# Propagate edges with leaf group id
	#cat("Propagating group ids up edges.\n");
	dend=dendrapply(dend, propagate_up_edges_with_leaf_group_id);
	#cat("done.\n");
	
	#dend=dendrapply(dend, color_dendrogram_by_group_id, color_map);
	#plot(dend);

	# Color nodes with homogenous parents
	#cat("Determing group under homogenous parents.\n");
	dend=mark_leafs_with_agreeing_edges(dend);
	#cat("done.\n");

	#dend=dendrapply(dend, color_dendrogram_by_group_id, color_map);
	#plot(dend, main=reference);

	# Get group id of all leaves
	inferred_leaf_groups=get_inferred_leaf_mapping(dend);
	
	return(inferred_leaf_groups);
	
}

#-------------------------------------------------------------------------------

color_dendrogram_by_group_id=function(n, colormap){
# If sets the end/leaf based on the group_id 
# TESTED: OK

	group_id=attr(n, "group_id");
	att=attributes(n);

	if(length(group_id)!=0){
		color=colormap[[group_id]];

		if(is.leaf(n)){
			#cat("Leaf:\n");
			attr(n, "nodePar")$lab.col=color;
			attr(n, "edgePar")$col=color;

			# Mark the antigens we have HI data for
			explicit_group_id=attr(n, "explicit_group_id");
			if(length(explicit_group_id)){
				if(explicit_group_id){
					attr(n, "nodePar")$pch=15;
					attr(n, "nodePar")$cex=.5;
				}else{
					attr(n, "nodePar")$cex=0;
				}
			}

			# Mark the reference antigen
			reference=attr(n, "reference");
			if(length(reference)){
				if(reference){
					attr(n, "nodePar")$pch="*";
					attr(n, "nodePar")$cex=2;
				}
			}
		}else{
			#cat("Edge:\n");
			attr(n, "edgePar")$col=color;
		}
	}
	
	return(n);
}


#-------------------------------------------------------------------------------

mark_leaf_nodes_with_group_id=function(n, group_map, reference){
# Sets all the leaves to the group_id in the group_map
# Also, sets the explicit flag, so we can later differentiate whether the group was inferred.
# TESTED: OK

	if(is.leaf(n)){
		leaf_attr=attributes(n);
		leaf_name=strsplit(leaf_attr$label, ":")[[1]][1];	# Remove the accession

		if(length(group_map[[leaf_name]])){
			group_id=group_map[[leaf_name]];
			explicit=TRUE;
		}else{
			group_id="UNKNOWN";
			explicit=FALSE;
		}

		if(leaf_name==reference){
			reference=TRUE;
		}else{
			reference=FALSE;
		}

		att=attributes(n);
		attributes(n)=c(att, list(group_id=group_id, explicit_group_id=explicit, reference=reference));
	}
	return(n);
}

#-------------------------------------------------------------------------------

propagate_up_edges_with_leaf_group_id=function(n){
	consensus_group_id=get_consensus_group_id(n);
	attr(n, "group_id") = consensus_group_id;
	return(n);
}

get_consensus_group_id=function(n){
# For a given node, it will determine if the group id's match.
# If they don't and one of them is unknown, then then the known
# one dominates
# TESTED: OK

        if(is.leaf(n)){
		group_id=attr(n, "group_id");
                if(length(group_id)==0){
                        return("UNKNOWN");
                }else{
                        return(group_id);
                }
        }else{
                num_nodes=length(n);
		groups=character(num_nodes);
                for(i in 1:num_nodes){
                        groups[i]=get_consensus_group_id(n[[i]]);
                }
                unique_group=unique(setdiff(groups, "UNKNOWN"));
                if(length(unique_group)==0){
                        return("UNKNOWN");
                }else if(length(unique_group)==1){
                        return(unique_group);
                }else{
                        return("CONFLICT"); # Because it's conflicting
                }
        }
}

#-------------------------------------------------------------------------------

assign_all_descendents_to_group=function(den, group_id){
	attr(den, "group_id")=group_id;
	return(den);
}

mark_leafs_with_agreeing_edges=function(n){
	if(is.leaf(n)){
		# Do nothing, end of recursion
		return(n);
	}else{
		#num_children=length(n);

		# If a edge, if both children are the same color, set their children to the same color.
		left_group_id=attr(n[[1]], "group_id");
		right_group_id=attr(n[[2]], "group_id");
	
		if((left_group_id==right_group_id) &&
			(left_group_id!="UNKNOWN") &&
			(left_group_id!="CONFLICT")){

			# If left and right edges agree, then assign every child below to that group.
			n[[1]]=dendrapply(n[[1]], assign_all_descendents_to_group, left_group_id);
			n[[2]]=dendrapply(n[[2]], assign_all_descendents_to_group, left_group_id);
		}else{
			
			# If left and right edges don't agree, then go down to children
			n[[1]]=mark_leafs_with_agreeing_edges(n[[1]]);
			n[[2]]=mark_leafs_with_agreeing_edges(n[[2]]);
		}

		return(n);
	}
}

#-------------------------------------------------------------------------------

get_inferred_leaf_mapping=function(n){
# Extracts out a mapping from name to group id
# TESTED: OK

	# Define recursive function locally because the mapping variable needs to be
	# accessible.
	save_name_to_group_mapping=function(den, mapping){
		if(is.leaf(den)){
			label=attr(den,"label");
			mapping[[label]]=attr(den, "group_id");
		}else{
			num_children=length(den)
			for(i in 1:num_children){
				mapping=save_name_to_group_mapping(den[[i]], mapping);
			}
		}
		return(mapping);
	}

	mapping=list();
	mapping=save_name_to_group_mapping(n, mapping);

	return(mapping);
}

#-------------------------------------------------------------------------------

accumulate_coverage=function(covered_sequences_list, coverage_counts){

	num_covered_sequences=length(covered_sequences_list);
	seq_names=names(covered_sequences_list);

	# Initialize the coverage counts
	if(length(coverage_counts)==0){
		coverage_counts=list();
		for(i in 1:num_covered_sequences){
			cur_name=seq_names[i]
			coverage_counts[[cur_name]][["IN"]]=0;
			coverage_counts[[cur_name]][["OUT"]]=0;
			coverage_counts[[cur_name]][["UNKNOWN"]]=0;
		}
	}

	# Increment the coverage counts
	for(i in 1:num_covered_sequences){
		cur_name=seq_names[i];
		if(covered_sequences_list[i]=="IN"){
			coverage_counts[[cur_name]][["IN"]]=
				coverage_counts[[cur_name]][["IN"]]+1;
		}else if(covered_sequences_list[i]=="OUT"){
			coverage_counts[[cur_name]][["OUT"]]=
				coverage_counts[[cur_name]][["OUT"]]+1;
		}else if(covered_sequences_list[i]=="UNKNOWN"){
			coverage_counts[[cur_name]][["UNKNOWN"]]=
				coverage_counts[[cur_name]][["UNKNOWN"]]+1;
		}
	}

	return(coverage_counts);
}

######################################################################################################

set_node_attributes=function(dend, cex="", pch="", col="", lab.fon="", lab.cex=""){
# Applies the the attributes to the node for the entire dendrogram unconditionally.
# If the entries are empty, the original value is left the same.
        if(is.leaf(dend)){
                att=attr(dend, "nodePar");

		if(cex!=""){
			att$cex=cex;
		}
		if(pch!=""){
			att$pch=pch;
		}
		if(col!=""){
			att$col=col;
		}
		if(lab.fon!=""){
			att$lab.fon=lab.fon;
		}
		if(lab.cex!=""){
			att$lab.cex=lab.cex;
		}

                attr(den, "nodePar")=att;
        }else{
		if(col!=""){
			attributes(dend)$edgePar[["col"]]=col;
		}
	}
        return(dend);

}

######################################################################################################

print_dendrogram=function(dend){
# Outputs all the attributes for the specified dendrogram
# For debugging.

	pden=function(dend){
		cat("----------------------------------------------------\n");
		print(attributes(dend));
	}
	dendrapply(dend, pden);
}

######################################################################################################

compute_percent_coverage=function(weighting_map, covered_sequences_list){

	num_sequences=length(covered_sequences_list);
	antigen=names(covered_sequences_list);

	if(length(weighting_map)>0){
		apply_wt=T;	
	}else{
		apply_wt=F;
	}

	n_in=0;
	n_out=0;
	n_unk=0;

	for(i in 1:num_sequences){
		type=covered_sequences_list[i];

		if(apply_wt){
			wt=weighting_map[[antigen[i]]];	
			if(length(wt)==0){
				wt=1;
			}
		}else{
			wt=1;
		}

		if(type=="IN"){
			n_in=n_in + wt;
		}else if(type=="OUT"){
			n_out=n_out + wt;
		}else if(type=="UNKNOWN"){
			n_unk=n_unk + wt;
		}else{
			cat("Error: unknown coverage type: ", type, ". Expeced IN/OUT/UNKNOWN.\n");
		}
	}

	n_total=n_in + n_out + n_unk;	

	coverage=list();
	coverage$num_in=n_in;
	coverage$num_out=n_out;
	coverage$num_unk=n_unk;
	coverage$total=n_total;
	coverage$perc_in=n_in/n_total;
	coverage$perc_out=n_out/n_total;
	coverage$perc_unk=n_unk/n_total;

	return(coverage);

}

######################################################################################################

get_CI = function(x, alpha){
        n=length(x);
        med=median(x);
        ba=1-(n-2)/n;
        #cat("Best alpha = ", ba, "\n");
        if(ba <= (alpha+.0000000000000001)){
                sorted=sort(x);
                lb=sorted[floor(n*(alpha/2))+1];
                ub=sorted[ceiling(n*(1-(alpha/2)))];
                return(c(med,lb,ub));
        }else{
                return(c(med,NA,NA))
        }
}

######################################################################################################

plot_summary_dendrogram=function(dend, coverage, title="", cov_ci, cov_obs, weighting_map, group_map, label_scale=-1){

	par(mar=c(15, 3, 3, .5));

	# Get attributes in dendrogram, from left to right
	get_attributes=function(dend, att){
		if(is.leaf(dend)){
			return(attr(dend, att));
		}else{
			val1=get_attributes(dend[[1]], att);
			val2=get_attributes(dend[[2]], att);
			return(c(val1,val2));
		}
	}

	# Remove labels
	hide_labels=function(dend){
		if(is.leaf(dend)){
			attr(dend, "nodePar")$cex=0;
			attr(dend, "nodePar")$lab.col=0;
		}
		return(dend);
	}
	

	labels=get_attributes(dend, "label");
	group_id=get_attributes(dend, "group_id");
	explicit=get_attributes(dend, "explicit");
	reference=get_attributes(dend, "reference");

	num_samples=length(labels);
	#print(labels);
	#print(group_id);
	#print(explicit);
	#print(reference);


	# Normalize coverage
	percs=list();
	for(i in 1:num_samples){
		cur_label=labels[i];
		n_in=coverage[[cur_label]][["IN"]];
		n_out=coverage[[cur_label]][["OUT"]];
		n_unk=coverage[[cur_label]][["UNKNOWN"]];
		n_tot= n_in + n_out + n_unk;

		if(length(n_tot)==0){
			percs[[cur_label]][["p_in"]]=0;
			percs[[cur_label]][["p_out"]]=0;
			#percs[[cur_label]][["p_unk"]]=1;
		}else{
			percs[[cur_label]][["p_in"]]=n_in/n_tot;
			percs[[cur_label]][["p_out"]]=n_out/n_tot;
			#percs[[cur_label]][["p_unk"]]=n_unk/n_tot;
		}
	}


	# Set ups plot
	dend_height=attr(dend, "height");	
	bar_height=dend_height*.15;
	spec_ind_height=dend_height*.05;
	padding=dend_height*.05;

	# Hide labels
	dend=dendrapply(dend, hide_labels);

	# Plot dendrogram and make room for barplots
	plot(dend, ylim=c(-(bar_height+spec_ind_height+padding), dend_height), main=title);
	
	# Plot each barplot
	bar_width=.80/2; # Half the thickness of the bar width you want
	line_width=bar_width/300;
	cat("Line Width: ", line_width, "\n");

	if(label_scale==-1){
		label_scale=50/num_samples;
		label_scale=min(c(1, label_scale));
	}
	cat("Label Scale: ", label_scale, "\n");
	for(i in 1:num_samples){

		cur_label=labels[i];		

		# unknown
		top=0;
		bottom=-bar_height;
		left=i-bar_width;
		right=i+bar_width;
		rect(left, bottom, right, top, col="grey", border=NA);

		# Ins
		top=-bar_height+(percs[[cur_label]][["p_in"]])*bar_height;
		bottom=-bar_height;
		left=i-bar_width;
		right=i+bar_width;
		rect(left, bottom, right, top, col="green", border=NA);

		# Outs 
		top=0
		bottom=-(percs[[cur_label]][["p_out"]])*bar_height;
		left=i-bar_width;
		right=i+bar_width;
		rect(left, bottom, right, top, col="red", border=NA);

		# Border
		top=0;
		bottom=-bar_height;
		left=i-bar_width;
		right=i+bar_width;
		rect(left, bottom, right, top, border="black", lwd=line_width);

		# Weights
		weight=weighting_map[[cur_label]];
		if(is.null(weight)){
			cat("\nError: ", cur_label, " not found in weighting map.\n\n");
		}
		text(i, -bar_height*.95, sprintf("%2.1f", weight), pos=1, cex=label_scale*.9, srt=90);
		
		# Special Indicator
		if(explicit[i]){

			if(group_id[i]=="IN"){
				spc_id_col="darkgreen";
			}else if(group_id[i]=="OUT"){
				spc_id_col="darkred";
			}else if(group_id[i]=="UNKNOWN"){
				spc_id_col="grey";
			}

			top=-bar_height-padding;
			bottom=top-spec_ind_height;
			left=i-bar_width;
			right=i+bar_width;
			
			if(reference[i]){
				marker_height=top-bottom;
				ref_top=-bar_height-padding +marker_height/10;
				ref_bottom=top-spec_ind_height -marker_height/10;
				rect(left, ref_bottom, right, ref_top, col="green", border="black", lwd=line_width)
			}

			# Mark the observed in/out
			rect(left, bottom, right, top, col=spc_id_col, border="black", lwd=line_width)

			if(reference[i]){
				points(i, (top+bottom)/2, pch="*", cex=2*label_scale, col="green");
			}

		}
	}

	# Plot labels on bottom of graph

	#axis(side=1, at=1:num_samples, labels=labels, las=2, cex.axis=label_scale);



	# Assign group a color from 1 to number of groups
	if(length(group_map)){
		groups=character();

		# Get the groups referenced in the labels
		for(i in 1:num_samples){
			groups[i]=group_map[[labels[i]]];
		}

		uniq_groups=sort(unique(groups));
		grp_color_map=list();

		# Assign group a color/index
		for(i in 1:length(uniq_groups)){
			grp_color_map[[ uniq_groups[i] ]]=i;
		}
		
		# Map label to group to color
		color=numeric();
		for(i in 1:length(labels)){
			color[i]=grp_color_map[[ group_map[[  labels[i] ]]  ]];
		}


		colors=as.integer(as.factor(groups));	
		legend(0,dend_height*.95, fill=1:length(uniq_groups), legend=uniq_groups, bg="white");
	}else{
		colors=rep(1, num_samples);
	}

	# Draw each label separately so we can assign a color
	for(i in 1:num_samples){
		axis(side=1, at=i, labels=labels[i], las=2, cex.axis=label_scale, col=colors[i], line=-.5, col.axis=colors[i]);
	}

	# Coverage Legend
	obs_str=sprintf(
		"Observed:\nIn: %3.1f%%  Out: %3.1f%%  Unknown: %3.1f%%",
			 100*cov_obs$perc_in, 100*cov_obs$perc_out, 100*cov_obs$perc_unk);
	
	ci_str=sprintf(
		"Bootstrapped Medians (95%% CI):\nIn: %3.1f%% (%3.1f, %3.1f)\nOut: %3.1f%% (%3.1f, %3.1f)\nUnknown: %3.1f%% (%3.1f, %3.1f)",
		100*cov_ci$perc_in[1], 100*cov_ci$perc_in[2], 100*cov_ci$perc_in[3],
		100*cov_ci$perc_out[1], 100*cov_ci$perc_out[2], 100*cov_ci$perc_out[3],
		100*cov_ci$perc_unk[1], 100*cov_ci$perc_unk[2], 100*cov_ci$perc_unk[3]);

	text(.85*num_samples, dend_height*.85, paste(obs_str, "\n\n", ci_str, sep=""), cex=.65);
	

}

######################################################################################################

plot_coverage_statistics=function(coverage_info, distances_from_candidate, inout_cutoff, ag_dist_range, antigen_name){

	#par(mfrow=c(2,3));
	layout_mat=matrix(c(
		1,1,2,2,3,3,
		4,4,4,5,5,5
		),
		byrow=T, ncol=6);
	layout(layout_mat);
	par(mar=c(4.5,4,4,2));
	par(oma=c(0,0,3,0));

	#-----------------------------------------------------------------------------

	# Compute the max counts so we can have the same limits on the 3 histograms
	bins=seq(0,1, length.out=11);
	inh=hist(coverage_info$perc_in, breaks=bins, plot=F);
	outh=hist(coverage_info$perc_out, breaks=bins, plot=F); 
	unkh=hist(coverage_info$perc_unk, breaks=bins, plot=F);
	max_counts=max(c(inh$counts, outh$counts, unkh$counts));

	# Plot the 3 histograms
	hist(coverage_info$perc_in, xlim=c(0,1), ylim=c(0, max_counts),  main="In", xlab="Coverage", col="green", breaks=bins);
	hist(coverage_info$perc_out, xlim=c(0,1), ylim=c(0, max_counts), main="Out", xlab="Coverage", col="red", breaks=bins);
	hist(coverage_info$perc_unk, xlim=c(0,1), ylim=c(0, max_counts), main="Unknown", xlab="Coverage", col="grey", breaks=bins);

	#-----------------------------------------------------------------------------

	# Compute the medians, so we have a single value for the pie chart
	inPerc=100*median(coverage_info$perc_in);
	outPerc=100*median(coverage_info$perc_out);
	unkPerc=100*median(coverage_info$perc_unk);

	# Generate the pie chart
	pie(c(inPerc, outPerc, unkPerc),
		labels=c(
			sprintf(c("In (%3.1f%%)", "Out (%3.1f%%)", "Unknown (%3.1f%%)"),
				c(inPerc, outPerc, unkPerc))
		), 
		col=c("green", "red", "grey"),
		main="Median Proportions"
	);

	#-----------------------------------------------------------------------------

	# Plot histogram of distances
	hist_break=seq(floor(ag_dist_range[1]), ceiling(ag_dist_range[2]), .25);
	hist_info=hist(distances_from_candidate$Distance, breaks=hist_break, plot=F);
	
	# Set up blank plot
	max_count=max(hist_info$counts);
	barpos=barplot(hist_info$counts, space=0, ylim=c(0, max_count*1.2), ylab="Frequency", xlab="Antigenic Distance");

	# Build translation from ag distance to plot location
	slope=(barpos[2]-barpos[1])/(hist_info$mids[2]-hist_info$mids[1]);
	intercept=barpos[1]-slope*hist_info$mids[1];
	conv=function(ag_dist){
		return(slope*ag_dist+intercept);
	}

	# Color in/out
	colors=rep("darkgreen", length(hist_info$mids));
	colors[hist_info$mids>inout_cutoff]="darkred";

	# Label cutoff
	abline(v=conv(inout_cutoff), col="grey", lty="dashed");
	text(conv(inout_cutoff), max_count*1.2, "Antigenic Distance Cutoff", srt=-90, adj=c(0,-.75), cex=.75);

	# Label median 
	median_ag=median(distances_from_candidate$Distance);
	abline(v=conv(median_ag), col="blue", lty="dashed");
	text(conv(median_ag), max_count*1.2, sprintf("Median = %2.2f", median_ag), srt=-90, adj=c(0,-.75), cex=.75);

	# Plot histogram
	barplot(hist_info$counts, space=0, add=T, col=colors);
	axis_labels=floor(ag_dist_range[1]):ceiling(ag_dist_range[2]);
	axis(side=1, conv(axis_labels), labels=sprintf("%2.1f", axis_labels));

	# Label page
	mtext(antigen_name, side=3, cex=1.2, font=2, outer=TRUE);
}

######################################################################################################

plot_all_relative_values=function(candidate_cov_list){

	par(mar=c(0,4,2,0));
	
	# Extract out information from intervals_list
	num_samples=length(candidate_cov_list);
	labels=names(candidate_cov_list);
	medians=numeric(num_samples);
	lb=numeric(num_samples);
	ub=numeric(num_samples);
	for(i in 1:num_samples){
		medians[i]=candidate_cov_list[[i]]$perc_in[1];		
		lb[i]=candidate_cov_list[[i]]$perc_in[2];		
		ub[i]=candidate_cov_list[[i]]$perc_in[3];		
	}

	# Sort the data from largest to smallest
	ord=order(medians, decreasing=T);
	medians=medians[ord];
	labels=labels[ord];
	labels=sprintf("%s (%3.1f%%)", labels, medians*100);

	# Plot axis
	padding=.5; # Padding between top and bottom of plot, in case we need spaces for labels
	x_pos=rep(0,num_samples);
	plot(x_pos, medians, ylab="Median Coverage", main="Sequence Coverage",
		ylim=c(-padding, 1+padding), xlim=c(-padding, 10), type="n", xlab="", xaxt="n", yaxt="n", bty="n");
	axis(side=2, at=seq(0,1,.1), labels=seq(0,100,10), las=2);
	
	# Plot tick marks
	for(i in 1:num_samples){
		segments(x_pos-.1, medians[i], x_pos+.1, medians[i], lwd=1, col="black");
	}

	# Produce spacings for labels so they don't overlap
	TXT_SPC=.03; INC=.001;
	rchange=T; fchange=T;
	text_pos=medians; 

	if(num_samples>1){
		while(rchange || fchange){

			# Adjust labels upwards
			for(i in 1:(num_samples-1)){
				diff=abs(text_pos[i]-text_pos[i+1]);
			
				fchange=F;
				if(diff<TXT_SPC){
					# Compute increment
					if((i-1)==0){
						inc=INC;
					}else{
						inc=min(abs(text_pos[i]-text_pos[i-1])/2, INC*1.05);
					}

					text_pos[i]=text_pos[i]+inc;
					fchange=T;
					break;
				}
			}

			# Adjust labels downwards
			for(i in seq(num_samples,2, -1)){
				diff=abs(text_pos[i-1]-text_pos[i]);

				rchange=F;
				if(diff<TXT_SPC){
					# Compute increment
					if((i+1)>num_samples){
						inc=INC;
					}else{
						inc=min(abs(text_pos[i]-text_pos[i+1])/2, INC/1.05);
					}

					text_pos[i]=text_pos[i]-inc;
					rchange=T;
					break;
				}
			}
		}
	}

	
	# Draw labels
	label_offset=.5;
	text(x_pos+label_offset, text_pos, labels, pos=4);
	for(i in 1:num_samples){
		segments(x_pos+label_offset-.1, text_pos[i], x_pos+label_offset, text_pos[i], lwd=1, col="black");
	}

	# Draw connectors
	for(i in 1:num_samples){
		segments(x_pos+.1, medians[i], x_pos+label_offset-.1, text_pos[i], lwd=1, col="black");
	}
	
}

######################################################################################################

output_coverage_statistics=function(fname, candidate_cov_list){
        num_samples=length(candidate_cov_list);
	
	fh=file(fname, "w");

	headers=c("# Antigen", "ObsIn", "ObsOut", "ObsUnk", 
		"NumSeqs", "DistIn", "DistOut",		
		"InMedian", "InLB", "InUB", "OutMedian", "OutLB", "OutUB", "UnkMedian", "UnkLB", "UnkUB");
	cat(file=fh, paste(headers, collapse=","));
	cat(file=fh, "\n");
	for(i in 1:num_samples){
		name=candidate_cov_list[[i]]$antigen;
	
		observed_cov=candidate_cov_list[[i]]$observed_coverage;	
		obs=c(observed_cov$perc_in, observed_cov$perc_out, observed_cov$perc_unk);
		distinout=c(candidate_cov_list[[i]]$inout$dist_in, candidate_cov_list[[i]]$inout$dist_out);

		num_seq=observed_cov$total;
		num_dist=candidate_cov_list[[i]]$num_distances;

		pin=candidate_cov_list[[i]]$perc_in;		
		pout=candidate_cov_list[[i]]$perc_out;		
		punk=candidate_cov_list[[i]]$perc_unk;		

		info=c(name, 
			obs, 
			num_seq,distinout, 
			pin,pout,punk);
		cat(file=fh, paste(info, collapse=","));
		cat(file=fh, "\n")
	}

	close(fh);

}



######################################################################################################

cat("Dendrogram/DASH Libraries Loaded...\n");

if(0){
	x=list();
	x[["Apples"]]=c(.10, .15, .20);
	x[["Bananas"]]=c(.15, .12, .20);
	x[["Cucumbers"]]=c(.90,.89 , .92);
	x[["Dinner"]]=c(.85, .80, .90);
	x[["Elephants"]]=c(.45, .40, .60);
	x[["Fortress"]]=c(.50, .3, .610);

	v=rbeta(20, 1, 2);
	vlab=1:20;
	x=list();
	for(i in 1:20){
		x[[sprintf("%i", vlab[i])]]=v[i];
	}

	plot_all_relative_values(x);
}

if(0){
	test=matrix(c(
		1,1,
		1,3,
		2,2,
		4,1,
		6,6,
		2,8,
		3,9,
		3,7,
		4,8
		),
		byrow=T, ncol=2);
	rownames(test)=c("A","B","C","D","E","F","G","This is a test", "I");
	#plot(test);

	hcl=hclust(dist(test));

	dend=as.dendrogram(hcl);
	#plot(dend);

	group_map=list();
	group_map[["A"]]="IN";
	group_map[["B"]]="IN";
	group_map[["D"]]="OUT";
	group_map[["E"]]="OUT";
	group_map[["F"]]="IN";
	group_map[["I"]]="OUT";
	
	coverage=list();

	coverage$A$IN=10;
	coverage$A$OUT=0;
	coverage$A$UNKNOWN=0;

	coverage$B$IN=0;
	coverage$B$OUT=0;
	coverage$B$UNKNOWN=10;

	coverage$D$IN=0;
	coverage$D$OUT=10;
	coverage$D$UNKNOWN=0;

	coverage$E$IN=50;
	coverage$E$OUT=50;
	coverage$E$UNKNOWN=10;

	coverage$F$IN=40;
	coverage$F$OUT=40;
	coverage$F$UNKNOWN=20;

	coverage$I$IN=10;
	coverage$I$OUT=10;
	coverage$I$UNKNOWN=10;


	plot_summary_dendrogram(dend, coverage);
	#infer_leaf_coverage(group_map, dend);
}

