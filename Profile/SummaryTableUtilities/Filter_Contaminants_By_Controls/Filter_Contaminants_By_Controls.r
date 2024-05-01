#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
	"input_summary_table", "i", 1, "character",
	"match_contam", "m", 2, "character",
	"avg_contam", "a", 2, "character",
	"mixture_contam", "x", 2, "character",
	"paired_contam", "p", 2, "character",
	"plevel", "l", 2, "numeric",
	"output_fname_root", "o", 2, "character",
	"short_name_delim", "d", 2, "character",
	"num_bs", "b", 2, "numeric",
	"counts", "c", 2, "numeric",
	"sim_anneal", "s", 2, "logical",
	"quantile", "q", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

options(width=200);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

DEFAULT_PLEVEL=2;
DEFAULT_DELIM=";";
DEFAULT_COUNTS=400;
DEFAULT_NUM_BS=8000;
DEFAULT_SIM_ANNEAL=F;
DEFAULT_QUANTILE=0.025;

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input samples summary table.tsv>\n",
	"\n",
	"	One of:\n",
	"	[-m <contamiment/negative control sample IDs list, to find closest match>]\n",
	"	[-a <contamiment/negative control sample IDs list, to pool/average together>]\n",
	"	[-x <contamiment/negative control sample IDs list, for mixture model.>]\n",
	"	[-p <list of experimental/contaminant paired sample IDs>]\n",
	"\n",
	"	[-l <p-norm used in the objective function, default=", DEFAULT_PLEVEL, ">]\n",
	"	[-o <output filename root>]\n",
	"	[-d <name shortening delimitor, default=", DEFAULT_DELIM, ">]\n",
	"\n",
	"	Bootstrapping parameters:\n",
	"	[-c <counts, estimated number of DNA molecules sampled, default=", DEFAULT_COUNTS, ">]\n",
	"	[-b <number of bootstraps, default=", DEFAULT_NUM_BS, ">]\n",
	"	[-s <simulated annealing, default=", DEFAULT_SIM_ANNEAL, ">]\n",
	"	[-q <quantile of best fit to keep, default=", DEFAULT_QUANTILE, ">]\n",
	"\n",	
	"This script will use the contaminant profiles/negative controls\n",
	"to adjust the taxa proportions from the experimental samples.\n",
	"The goal is essentially to try to find the distribution of the\n",
	"contaminant sample inside the experimental sample and then substract\n",
	"the contaminant distribution from the experimental distribution.\n",
	"\n",
	"The output is a summary table consisting only of the filtered\n",
	"experimental samples, and some diagnostic pdf file.\n",
	"The counts in the summary table is be adjusted to reflected\n",
	"the reduced proportion of reads per sample.\n",
	"\n",
	"The algorithm works by using a mixture model consisting of\n",
	"an unknown ratio of experimental to contaminant.\n",
	"\n",
	"	PureExperimental+ c*Contaminant = ObservedExperimental\n",
	"\n",
	"The goal is to try to calculate the coefficient c.\n",
	"We will calculate c by minimizing the following objective function:\n",
	"\n",
	"	Sum (abs(Exp[i]-c*Cont[i]))^p\n",
	"	 where i is only the categories in the Cont\n",
	"\n",
	"	With the constraint that c is in [0,Inf].\n",
	"\n",
	"	Alternatively, the objective function can be reformulated as:\n",
	"	Sum (abs((x*Exp[i]-(1-x)*Cont[i]))^p\n",
	"	so that x is between [0,1], which makes the search easier and\n",
	"	easier to plot.\n",
	"\n",
	"	p is the p-norm degree:\n",
	"	When p=2, the fit is a least squares.\n",
	"	When p=1, the fit is absolute value.\n",
	"\n",
	"	The smaller the p, e.g. .5, the less the effect of outliers,\n",
	"	but the fit can become less robust.\n",
	"	For contaminant removal, outliers may be actual levels of\n",
	"	non-contaminant, so we want the fit to be spread out across\n",
	"	all categories, and not dominated by a few categories\n",
	"	that fit poorly.\n",
	"	\n",
	"	You can think of it this way: smaller p levels are like electoral votes, \n",
	"	whereas larger p levels are popular votes.\n",
	"\n",
	"\n",
	"There are two options to apply the contaminant filter:\n",
	"	1.) Averaging (-a): Take an average of contaminants and apply a single\n",
	"		profile as the contaminent to each experimental sample.\n",
	"	2.) Paired (-p): If each sample is paired with a negative control,\n",
	"		then apply a single negative control with a single\n",
	"		experimental sample, e.g. BAL / BAL Control.\n",
	"\n",
	"The format for the -a negative control file:\n",
	"	<contam1 sample id>\\n\n",
	"	<contam2 sample id>\\n\n",
	"	<contam3 sample id>\\n\n",
	"	...\n",
	"	<contamn sample id>\\n\n",
	"\n",
	"The format for the -p paired exp/cont file:\n",
	"	<exper1 sample id>\\t<contam1 sample id>\\n\n",	
	"	<exper2 sample id>\\t<contam2 sample id>\\n\n",	
	"	<exper3 sample id>\\t<contam3 sample id>\\n\n",	
	"	...\n",
	"	<expern sample id>\\t<contamn sample id>\\n\n",	
	"\n",
	"The bootstrapping comes with several parameters:\n",
	"	* Counts (-c) specifies the number of sequences that you think the PCR\n",
	"	reaction actually amplified.  The larger this number the more\n",
	"	robust you think the distribution of the control to be and\n",
	"	less variation will be introduced into the perturbed control\n",
	"	distribution upon bootstrapping.\n",
	"	* Num Boostraps (-b) is the number of variations of distributions to try.\n",
	"	* Simulated Annealing (-s) is a tweak to the resampling algorithm.\n",
	"	Instead of assuming the counts is always the same, the counts become\n",
	"	the uppower bound, so that we start from 1 and iterate up to the Counts.\n",
	"	Thus the more bootstrap you allow, the more repeats at the same count\n",
	"	you will have.  \n",
	"	* Quantile (-q) is the best matching distribution you want to keep.\n",
	"	Smaller values will take the better fit.  Don't use -q 0.0, because\n",
	"	the best may be a outlier of perturbed distributions.\n",	
	"\n", sep="");

if(!length(opt$input_summary_table)){
	cat(usage);
	q(status=-1);
}

###############################################################################

InputSummaryTable=opt$input_summary_table;
AvgContamFName=opt$avg_contam;
PairedContamFName=opt$paired_contam;
MatchContamFName=opt$match_contam;
MixtureContamFName=opt$mixture_contam;

PLevel=opt$plevel;
OutputFileRoot=opt$output_fname_root;
ShortNameDelim=opt$short_name_delim;
NumBS=opt$num_bs;
Counts=opt$counts;
SimAnneal=opt$sim_anneal;
Quantile=opt$quantile;


if(!length(OutputFileRoot)){
	OutputFileRoot=gsub("\\.summary_table\\.xls", "", InputSummaryTable);
	OutputFileRoot=gsub("\\.summary_table\\.tsv", "", OutputFileRoot);
}

if(!length(PLevel)){
	PLevel=DEFAULT_PLEVEL;
}

if(!length(ShortNameDelim)){
	ShortNameDelim=DEFAULT_DELIM;
}

if(!length(NumBS)){
	NumBS=DEFAULT_NUM_BS;
}

if(!length(Counts)){
	Counts=DEFAULT_COUNTS;
}

if(!length(SimAnneal)){
	SimAnneal=DEFAULT_SIM_ANNEAL;
}

if(!length(Quantile)){
	Quantile=DEFAULT_QUANTILE;
}

if(SimAnneal){
	Quantile=0.0;
	cat("Quantile should be 0.0 for Simulated Annealing, because the resampling size keeps changing.\n");
	cat("Just taking the best resample...\n");
}

percentile_str=sprintf("%3.2f%%", 100*(1-Quantile));

cat("\n")
cat("Input Summary File Table: ", InputSummaryTable, "\n", sep="");
cat("Output Filename Root: ", OutputFileRoot, "\n", sep="");
cat("P-norm level: ", PLevel, "\n", sep="");
cat("Short Name Delim: ", ShortNameDelim, "\n", sep="");
cat("Num Bootstraps: ", NumBS, "\n", sep="");
cat("Counts: ", Counts, "\n", sep="");
cat("Perform 'simulated annealing':", SimAnneal, "\n", sep="");
cat("Quantile (1-alpha), alpha: ", Quantile, "\n", sep="");
cat("\n");

doPaired=NULL;
doMatched=NULL;
doAveraged=NULL;
doMixture=NULL;

method="";

if(length(PairedContamFName)){
	cat("Performing paired contaminant removal: ", PairedContamFName, "\n");
	doPaired=TRUE;
	doMatched=FALSE;
	doAveraged=FALSE;
	doMixture=FALSE;
	method="paired";
}else if(length(AvgContamFName)){
	cat("Performing averaged contaminant removal: ", AvgContamFName, "\n");
	doPaired=FALSE;
	doMatched=FALSE;
	doAveraged=TRUE;
	doMixture=FALSE;
	method="averaged";
}else if(length(MatchContamFName)){
	cat("Performing matched contaminant removal: ", MatchContamFName, "\n");
	doPaired=FALSE;
	doMatched=TRUE;
	doAveraged=FALSE;
	doMixture=FALSE;
	method="matched";
}else if(length(MixtureContamFName)){
	cat("Performing mixture contaminant removal: ", MixtureContamFName, "\n");
	doPaired=FALSE;
	doMatched=FALSE;
	doAveraged=FALSE;
	doMixture=TRUE;
	method="mixture";
}

OutputFileRoot=paste(OutputFileRoot, ".", method, sep="");


###############################################################################
###############################################################################

load_summary_table=function(summary_table_fn){
	# Load data
	cat("Loading Matrix (", summary_table_fn, ") ...\n", sep="");
	inmat=as.matrix(read.table(summary_table_fn, sep="\t", header=TRUE, 
		check.names=FALSE, row.names=1))

	#cat("\nOriginal Matrix:\n")
	#print(inmat);

	# Grab columns we need into a vector, ignore totals, we won't trust it.
	counts_mat=inmat[,2:(ncol(inmat))];
	#cat("\nCounts Matrix:\n");
	#print(counts_mat);

	num_samples=nrow(counts_mat);
	num_categories=ncol(counts_mat);
	sample_names=rownames(counts_mat);

	cat("\n");
	cat("Num Samples: ", num_samples, "\n");
	cat("Num Categories: ", num_categories, "\n");
	cat("\n");
	return(counts_mat);
}

###############################################################################

load_ids=function(list_fn){
	cat("Loading List (", list_fn, ") ...\n", sep="");
        list=as.vector(read.table(list_fn)[,1]);
	#print(list);
	return(list);
}

###############################################################################

load_pairs=function(pairs_fn){
	cat("Loading Pairs (", pairs_fn, ") ...\n", sep="");
	pairs=read.delim(pairs_fn, sep="\t", row.names=1, header=F, as.is=T);
	src_names=rownames(pairs);
	dst_names=pairs[,1];
	names(dst_names)=src_names;
	return(dst_names);
}

###############################################################################

normalize=function(counts){
	# Sum sample totals
	sample_totals=numeric();
	num_samples=nrow(counts);
	num_categories=ncol(counts);

	for(i in 1:num_samples){
		sample_totals[i]=sum(counts_mat[i,]);
	}
	#print(sample_totals);

	# normalize, to compute probabilities
	normalized=matrix(0, nrow=num_samples, ncol=num_categories);
	for(i in 1:num_samples){
		normalized[i,]=counts_mat[i,]/sample_totals[i];
	}

	# Preserve the names
	colnames(normalized)=colnames(counts);
	rownames(normalized)=rownames(counts);

	return(normalized);
}

###############################################################################

plot_text=function(strings, max_lines_pp=Inf){

        orig.par=par(no.readonly=T);

        par(mfrow=c(1,1));
        par(family="Courier");
        par(oma=rep(.5,4));
        par(mar=rep(0,4));

        num_lines=length(strings);
        num_pages=max(1, ceiling(num_lines/max_lines_pp));

        cat("Num Pages for ", num_lines, " lines: ", num_pages, "\n", sep="");

        lines_pp=min(num_lines, max_lines_pp);
        for(p in 1:num_pages){

                top=max(as.integer(lines_pp), 52);

                plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
                        xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                        );

                text_size=max(.01, min(.7, .7 - .003*(lines_pp-52)));
                #print(text_size);

                start=(p-1)*lines_pp+1;
                end=start+lines_pp-1;
                end=min(end, num_lines);
                line=1;
                for(i in start:end){
                        #cat(strings[i], "\n", sep="");
                        strings[i]=gsub("\t", "", strings[i]);
                        text(0, top-line, strings[i], pos=4, cex=text_size);
                        line=line+1;
                }

        }

        par(orig.par);
}

###############################################################################

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, deci_pts=4,
        label_zeros=T, counts=F, value.cex=1,
        plot_col_dendr=F,
        plot_row_dendr=F
){

        num_row=nrow(mat);
        num_col=ncol(mat);

        row_names=rownames(mat);
        col_names=colnames(mat);

        orig.par=par(no.readonly=T);

        cat("Num Rows: ", num_row, "\n");
        cat("Num Cols: ", num_col, "\n");

        # Flips the rows, so becuase origin is bottom left
        mat=mat[rev(1:num_row),, drop=F];

        # Generate a column scheme
        num_colors=50;
        color_arr=rainbow(num_colors, start=0, end=4/6);
        if(high_is_hot){
                color_arr=rev(color_arr);
        }

        # Provide a means to map values to an (color) index
        remap=function(in_val, in_range, out_range){
                in_prop=(in_val-in_range[1])/(in_range[2]-in_range[1])
                out_val=in_prop*(out_range[2]-out_range[1])+out_range[1];
                return(out_val);
        }

        # If range is not specified, find it based on the data
        if(is.na(plot_min)){
                plot_min=min(mat, na.rm=T);
        }
        if(is.na(plot_max)){
                plot_max=max(mat, na.rm=T);
        }

        if(plot_min>=-1 && plot_max<=1){
                fractions_only=T;
        }else{
                fractions_only=F;
        }
        cat("Plot min/max: ", plot_min, "/", plot_max, "\n");

        # Get Label lengths
        row_max_nchar=max(nchar(row_names));
        col_max_nchar=max(nchar(col_names));
        cat("Max Row Names Length: ", row_max_nchar, "\n");
        cat("Max Col Names Length: ", col_max_nchar, "\n");

        ##################################################################################################

        get_dendrogram=function(in_mat, type){
                if(type=="row"){
                        dendist=dist(in_mat);
                }else{
                        dendist=dist(t(in_mat));
                }

                get_clstrd_leaf_names=function(den){
                # Get a list of the leaf names, from left to right
                        den_info=attributes(den);
                        if(!is.null(den_info$leaf) && den_info$leaf==T){
                                return(den_info$label);
                        }else{
                                lf_names=character();
                                for(i in 1:2){
                                        lf_names=c(lf_names, get_clstrd_leaf_names(den[[i]]));
                                }
                                return(lf_names);
                        }
                }


                hcl=hclust(dendist, method="ward.D2");
                dend=list();
                dend[["tree"]]=as.dendrogram(hcl);
                dend[["names"]]=get_clstrd_leaf_names(dend[["tree"]]);
                return(dend);
        }


        ##################################################################################################
        # Comput Layouts
        col_dend_height=ceiling(num_row*.1);
        row_dend_width=ceiling(num_col*.2);

        heatmap_height=num_row;
        heatmap_width=num_col;

        if(plot_col_dendr && plot_row_dendr){
                layoutmat=matrix(
                        c(
                        rep(c(rep(4, row_dend_width), rep(3, heatmap_width)), col_dend_height),
                        rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height)
                        ), byrow=T, ncol=row_dend_width+heatmap_width);

                col_dendr=get_dendrogram(mat, type="col");
                row_dendr=get_dendrogram(mat, type="row");

                mat=mat[row_dendr[["names"]], col_dendr[["names"]]];

        }else if(plot_col_dendr){
                layoutmat=matrix(
                        c(
                        rep(rep(2, heatmap_width), col_dend_height),
                        rep(rep(1, heatmap_width), heatmap_height)
                        ), byrow=T, ncol=heatmap_width);

                col_dendr=get_dendrogram(mat, type="col");
                mat=mat[, col_dendr[["names"]]];

        }else if(plot_row_dendr){
                layoutmat=matrix(
                        rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height),
                        byrow=T, ncol=row_dend_width+heatmap_width);

                row_dendr=get_dendrogram(mat, type="row");
                mat=mat[row_dendr[["names"]],];
        }else{
                layoutmat=matrix(
                        rep(1, heatmap_height*heatmap_width),
                        byrow=T, ncol=heatmap_width);
        }

        #print(layoutmat);
        layout(layoutmat);

        ##################################################################################################

        par(oma=c(col_max_nchar*.60, 0, 3, row_max_nchar*.60));
        par(mar=c(0,0,0,0));
        plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), xaxt="n", 
		yaxt="n", bty="n", xlab="", ylab="");
        mtext(title, side=3, line=0, outer=T, font=2);

        # x-axis
        axis(side=1, at=seq(.5, num_col-.5, 1), labels=colnames(mat), las=2, line=-1.75);
        axis(side=4, at=seq(.5, num_row-.5, 1), labels=rownames(mat), las=2, line=-1.75);

        if(log_col){
                plot_min=log10(plot_min+.0125);
                plot_max=log10(plot_max+.0125);
        }

        for(x in 1:num_col){
                for(y in 1:num_row){

                        if(log_col){
                                col_val=log10(mat[y,x]+.0125);
                        }else{
                                col_val=mat[y,x];
                        }

                        remap_val=remap(col_val, c(plot_min, plot_max), c(1, num_colors));
                        col_ix=ceiling(remap_val);

                        rect(x-1, y-1, (x-1)+1, (y-1)+1, border=NA, col=color_arr[col_ix]);

                        if(is.na(mat[y,x]) || mat[y,x]!=0 || label_zeros){
                                if(counts){
                                        text_lab=sprintf("%i", mat[y,x]);
                                }else{
                                        text_lab=sprintf(paste("%0.", deci_pts, "f", sep=""), mat[y,x]);
                                        if(fractions_only){
                                                if(!is.na(mat[y,x])){
                                                        if(mat[y,x]==-1 || mat[y,x]==1){
                                                                text_lab=as.integer(mat[y,x]);
                                                        }else{
                                                                text_lab=gsub("0\\.","\\.", text_lab);
                                                        }
                                                }
                                        }
                                }
                                text(x-.5, y-.5, text_lab, srt=atan(num_col/num_row)/pi*180, 
					cex=value.cex, font=2);
                        }
                }
        }

        ##################################################################################################

        par(mar=c(0, 0, 0, 0));

        if(plot_row_dendr && plot_col_dendr){
                rdh=attributes(row_dendr[["tree"]])$height;
                cdh=attributes(col_dendr[["tree"]])$height;
                plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", 
			bty="n", xlim=c(rdh, 0));
                plot(col_dendr[["tree"]], leaflab="none",xaxt="n", yaxt="n", bty="n", ylim=c(0, cdh));
                plot(0,0, type="n", bty="n", xaxt="n", yaxt="n");
                #text(0,0, "Placeholder");
        }else if(plot_row_dendr){
                rdh=attributes(row_dendr[["tree"]])$height;
                plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", 
			bty="n", xlim=c(rdh, 0));
                #text(0,0, "Row Dendrogram");
        }else if(plot_col_dendr){
                cdh=attributes(col_dendr[["tree"]])$height;
                plot(col_dendr[["tree"]], leaflab="none", xaxt="n", yaxt="n", bty="n", ylim=c(0, cdh));
                #text(0,0, "Col Dendrogram");
        }

        par(orig.par);

}

###############################################################################

plot_RAcurve=function(ordered_composition, title="", top=0, max=1, overlay_dist=NULL){
	cat("Plotting: ", title, "\n", sep="");
	names=names(ordered_composition);
	num_cat=length(names);
	if(top<=0){
		top=num_cat;
	}
	mids=barplot(ordered_composition[1:top], names.arg="", las=2, 
		cex.names=.75, cex.axis=.75,
		main=paste("\n\n",title, sep=""), ylim=c(0,max));

	if(!is.null(overlay_dist)){
		points(mids, overlay_dist[1:top], pch="o", col="blue");	
	}

	cw=par()$cxy;
	labsize=min(c(1, 1/cw[1]));
	text(mids-cw[1]/2, rep(-cw[2]/2, top), names[1:top], srt=-45, xpd=T, pos=4, cex=labsize);
}

plot_RAstrip=function(profile_mat, title="", top=10){

	mean_arr=apply(profile_mat, 2, mean);
	order_ix=order(mean_arr, decreasing=T);

	top_prof_mat=profile_mat[,order_ix[1:top]];
	num_samp=nrow(top_prof_mat);

	plot(0, type="n", xlim=c(0, top+1), ylim=c(0,1), xaxt="n", bty="n", xlab="",
		main=title);

	jitter=rnorm(num_samp, 0, .1);

	for(i in 1:top){
		points(i+jitter, top_prof_mat[,i], cex=.5);
		points(i, median(top_prof_mat[,i]), cex=1.5, col="green", pch="-");
		points(i, mean(top_prof_mat[,i]), cex=1.5, col="blue", pch="-");
	}

	cw=par()$cxy;
	labsize=min(c(1, 1/cw[1]));
	text((1:top)-cw[1]/2, rep(-cw[2]/2, top), colnames(top_prof_mat), srt=-45, xpd=T, pos=4, cex=labsize);
}

plot_RAprop=function(profile_mat, title="", top=10){

	mean_arr=apply(profile_mat, 2, mean);
	order_ix=order(mean_arr, decreasing=T);

	top_prof_mat=profile_mat[,order_ix[1:top]];
	num_samp=nrow(top_prof_mat);

	num_any_removed=apply(top_prof_mat, 2, function(x){sum(x>0)});

	par(mar=c(10, 5, 5, 5));
	mids=barplot(num_any_removed, names.arg="", ylim=c(0, num_samp), main=title,
		xlab="", ylab="Number of Samples", las=2
		);

	tick_pos=c(0, .25, .33, .50, .66, .75, 1);
	axis(side=4, at=tick_pos*num_samp, labels=tick_pos, las=2);

	cw=par()$cxy;
	labsize=min(c(1, 1/cw[1]));
	text(mids-cw[1]/2, rep(-cw[2]/2, top), colnames(top_prof_mat), srt=-45, xpd=T, pos=4, cex=labsize);
}

plot_RAbox=function(ordered_composition, title="", top=0, max=1){
	#cat("Plotting: ", title, "\n", sep="");
	names=colnames(ordered_composition);
	num_cat=length(names);
	if(top<=0){
		top=num_cat;
	}
	subsamp=sample(1:nrow(ordered_composition), min(80, nrow(ordered_composition)), replace=F);
	boxplot(ordered_composition[subsamp,1:top], names=rep("", top), las=2, 
		cex.names=.75, cex.axis=.75,
		bty="n",
		main=paste("\n\n",title, sep=""), ylim=c(0,max));

	cw=par()$cxy;
	labsize=min(c(1, 1/cw[1]));
	text((1:top)-cw[1]/2, rep(-cw[2], top), names[1:top], srt=-45, xpd=T, pos=4, cex=labsize);
}

###############################################################################

fit_contaminant_mixture_model=function(contam_comp, exper_comp, plevel){

	num_cat=length(contam_comp);
	cont_idx=contam_comp>0;

	# Extract the categories that are non-zero in contaminant
	expr_sub=exper_comp[cont_idx];
	cont_sub=contam_comp[cont_idx];

	# Define the objective function
	objective=function(p){
		diff=(abs((1-p)*expr_sub-p*cont_sub))^plevel;
		score=sum(diff);
		return(score);
	}

	# Perform the optimization
	optim_res=optim(par=0, objective, method="Brent", lower=0, upper=1);
	contam_prop=optim_res$par;
	#cat("Contaminant Proportion: ", contam_prop, "\n", sep="");
	contam_scale=contam_prop/(1-contam_prop);
	#cat("Contaminant Scale: ", contam_scale, "\n", sep="");

	# Remove contaminant from experimentals
	adjusted_comp=exper_comp - contam_scale*contam_comp;

	# Remove/zero out negative composition
	negative_cat_ix=adjusted_comp<0;
	adjusted_comp[negative_cat_ix]=0;

	# Compute normalized filtered composition (so distribution adds up to 1 again)
	normalized_filt_comp=adjusted_comp/sum(adjusted_comp);
	names(normalized_filt_comp)=names(exper_comp);

	# Compute normalized contam composition, after filtred exp has been re-normalized
	removed_contam_prof=exper_comp-normalized_filt_comp;
	removed_contam_prof[removed_contam_prof<0]=0;
	normalized_contam_comp=removed_contam_prof/sum(removed_contam_prof);

	# Estimate Proportion removed (after zero-ing out negatives)
	proportion_removed=sum(exper_comp-adjusted_comp);

	# Output results in structure/list
	fit=list();
	fit$multiplier=contam_scale;
	fit$proportion=contam_prop;
	fit$removed=proportion_removed;
	fit$cleaned=normalized_filt_comp;
	fit$removed_contam_profile=normalized_contam_comp;
	return(fit);
}

###############################################################################

counts_mat=load_summary_table(InputSummaryTable);
counts_total=apply(counts_mat, 1, sum);
summary_table_sample_ids=rownames(counts_mat);
num_categories=ncol(counts_mat);

# Clean names
long_names=colnames(counts_mat);
split_names=strsplit(long_names, ShortNameDelim);
num_names=length(long_names);
short_names=character();
for(i in 1:num_names){
	short_names[i]=tail(split_names[[i]],1);
}

colnames(counts_mat)=short_names;
normalized_mat=normalize(counts_mat);
original_names=long_names;

# Load IDs
avg_cont_dist=NULL;
experm_samples=character();
contam_samples=character();

if(doPaired){

	pairs=load_pairs(PairedContamFName);
	experm_samples=intersect(names(pairs), summary_table_sample_ids);
	contam_samples=intersect(pairs, summary_table_sample_ids);

	# Remove pairs not complete in summary table
	valid_pairs_list=list();
	for(exp in experm_samples){
		if(any(pairs[exp]==contam_samples)){
			valid_pairs_list[[exp]]=pairs[[exp]];
		}else{
			cat("Missing contaminant id for experimental id: ", exp, "\n");
		}
	}
	pairs=unlist(valid_pairs_list);

	experm_samples=names(pairs);
	contam_samples=pairs;
	names(contam_samples)=NULL;

}else if(doAveraged){

	ids_array=load_ids(AvgContamFName);
	overlapping_ids=intersect(ids_array, summary_table_sample_ids);

	contam_mat=normalized_mat[overlapping_ids,, drop=F];
	avg_cont_dist=apply(normalized_mat, 2, mean);

	experm_samples=setdiff(summary_table_sample_ids, ids_array);
	contam_samples=overlapping_ids;

}else if(doMatched){

	ids_array=load_ids(MatchContamFName);
	overlapping_ids=intersect(ids_array, summary_table_sample_ids);

	contam_mat=normalized_mat[overlapping_ids,, drop=F];

	experm_samples=setdiff(summary_table_sample_ids, ids_array);
	contam_samples=overlapping_ids;

}else if(doMixture){

	ids_array=load_ids(MixtureContamFName);
	overlapping_ids=intersect(ids_array, summary_table_sample_ids);

	contam_mat=normalized_mat[overlapping_ids,, drop=F];

	experm_samples=setdiff(summary_table_sample_ids, ids_array);
	contam_samples=overlapping_ids;

	mixture_component_matrix=matrix(NA, nrow=length(experm_samples), ncol=length(contam_samples));
	colnames(mixture_component_matrix)=contam_samples;
	rownames(mixture_component_matrix)=experm_samples;
}

cat("Experimental Samples:\n");
print(experm_samples);

cat("Contaminant Samples:\n");
print(contam_samples);

# Reorder the normalized matrix so that we'll see both exp and ctl 
# categories in the rank abundance plot
mean_exp_abund=apply(normalized_mat[experm_samples,,drop=F], 2, mean);
mean_con_abund=apply(normalized_mat[contam_samples,,drop=F], 2, mean);
mean_both_abund=(mean_exp_abund+mean_con_abund)/2;
mean_order=order(mean_both_abund, decreasing=T);

# Need to keep original names, because we will eventually export a summary table with full names
normalized_mat=normalized_mat[,mean_order, drop=F];
original_names=original_names[mean_order];

contam_mat=normalized_mat[contam_samples,,drop=F];

if(doAveraged){
	avg_cont_dist=avg_cont_dist[mean_order];
}


###############################################################################

perturb_dist_classical=function(distr, counts, num_bs){

	cat("Classical Perturbation:\n");

	# Generate random distributions based input dist
	pert=t(rmultinom(num_bs, size=counts, prob=distr));	

	# Normalize
	for(i in 1:num_bs){
		pert[i,]=pert[i,]/counts;	
	}
	return(pert);
}

perturb_dist_sim_anneal=function(distr, counts, num_bs){

	cat("Simulated Annealing Perturbation:\n");

	num_bs=max(counts, num_bs);
	reps=ceiling(num_bs/counts);

	pert_mat=matrix(0, nrow=reps*counts, ncol=length(distr));
	for(cts in 1:counts){
		pert_reps=t(rmultinom(reps, size=cts, prob=distr));	
		pert_reps=pert_reps/cts;

		pert_mat[
			(((cts-1)*reps)
			:
			(cts*reps-1))+1,]=pert_reps;
	}

	return(pert_mat);
}

perturb_dist=function(distr, counts, num_bs, sim_anneal=F){
	if(sim_anneal){
		pert=perturb_dist_sim_anneal(distr, counts, num_bs);
	}else{
		pert=perturb_dist_classical(distr, counts, num_bs);
	}
	return(pert);
}

bootstrp_fit=function(exper_dist, pert_ctl_dist_matrix, plevel){

	# Number of BS to perform depends on number of rows/pertrubations of ctrl
	num_bs=nrow(pert_ctl_dist_matrix);

	# Statistic of removal
	results_mat=matrix(0, nrow=num_bs, ncol=3);
	colnames(results_mat)=c("multiplier", "proportion", "removed");
	
	# Distribution after removing contaminants
	cleaned_dist_mat=matrix(-1, nrow=num_bs, ncol=length(exper_dist));
	colnames(cleaned_dist_mat)=names(exper_dist);

	# Loop for fitting mixture model
	for(bs_ix in 1:num_bs){
		fit=fit_contaminant_mixture_model(pert_ctl_dist_matrix[bs_ix,], exper_dist, plevel);

		results_mat[bs_ix, "multiplier"]=fit$multiplier;
		results_mat[bs_ix, "proportion"]=fit$proportion;
		results_mat[bs_ix, "removed"]=fit$removed;
		cleaned_dist_mat[bs_ix, ]=fit$cleaned;
	}
	
	# Package in list before returning
	results=list();
	results$cleaned=cleaned_dist_mat;
	results$stats=results_mat;

	return(results);
}

###############################################################################

find_closest=function(query, references, p=1){

	num_ref=nrow(references);
	distances=numeric(num_ref);
	ref_names=rownames(references);
	names(distances)=ref_names;

	for(i in 1:num_ref){

		cur_ref=references[i,];
		non_zero_ref=(cur_ref>0);

		dist=sum((abs(query[non_zero_ref]-cur_ref[non_zero_ref]))^p)^(1/p);
		distances[i]=dist;
	}
	
	#cat("Distances:\n");
	#print(distances);
	min_dist=min(distances);
	closest_ref=min(which(distances==min_dist));
	return(ref_names[closest_ref]);
	
}

find_mixture=function(query, references, p=1){

	num_ref=nrow(references);
	ref_names=rownames(references);

	num_cat=ncol(references);

	logistic_transform=function(real, s=250){
		# converts from from -inf to info to 0 to 1
		# so we don't have to deal with bounds in optimizer
		prop=1/(1+exp(-real/s));
		return(prop);
	}

	exponential_transform=function(prop, s=250){
		real=-s*(log((1/prop)-1));
		return(real);
	}

	make_mixture=function(weights, prof_mat){
		weighted_mat=prof_mat;
		for(i in 1:nrow(prof_mat)){
			weighted_mat[i,]=weighted_mat[i,]*weights[i];
		}
		sum_prof=apply(weighted_mat, 2, sum);
		norm_prof=sum_prof/sum(sum_prof);
		return(norm_prof);
	}

	obj_contam_dist_fun=function(param){
		weights=logistic_transform(param);
		#print(weights);
		mix_prof=make_mixture(weights, references);
		non_zero_ref=(mix_prof>0);
		dist=sum((abs(query[non_zero_ref]-mix_prof[non_zero_ref]))^p)^(1/p);
		return(dist);
	}

	# Instead of using bounded BFGS, with upper/lower bounds between 0-1, use 
	# 1/(1+exp(-x)) to allow x to be unbounded while leaving the transformed value bounded
	init_param=rep(exponential_transform(1/num_ref), num_ref);
	opt_res=optim(init_param, obj_contam_dist_fun, method="Nelder-Mead");

	opt_weights=logistic_transform(opt_res$par);
	opt_weights_norm=opt_weights/sum(opt_weights);
	names(opt_weights_norm)=ref_names;
	greatest_weight=max(opt_weights_norm);
	key_ref=ref_names[min(which(greatest_weight==opt_weights_norm))];

	results_rec=list();
	results_rec[["categorical_profile"]]=make_mixture(opt_weights_norm, references);
	results_rec[["mixture_name"]]=
		paste("Mix_", key_ref, "_", sprintf("%2.0f", greatest_weight*100), sep="");
	results_rec[["component_profile"]]=opt_weights_norm;

	print(results_rec);

	return(results_rec);

}


###############################################################################

pdf(paste(OutputFileRoot, ".dist_contam.pdf", sep=""), height=14, width=8.5);

top_cat_to_plot=45;

###############################################################################

par(mfrow=c(8,2));
mar=par()$mar;
mar=c(8.1, 3.1, 2.1, 1);
#mar=c(6.1, 3.1, 3.1, 1);
par(mar=mar);

layout_mat=matrix(c(
1,1,1,1,
2,2,3,3,
4,4,5,5,
6,6,7,7,
8,8,9,9
), byrow=T, ncol=4);

layout(layout_mat);

num_exp_cat=ncol(normalized_mat);
num_exp_samples=length(experm_samples);


cleaned_obs_matrix=matrix(-1, nrow=num_exp_samples, ncol=num_exp_cat, 
	dimnames=list(experm_samples, original_names));
cleaned_bs_matrix=matrix(-1, nrow=num_exp_samples, ncol=num_exp_cat, 
	dimnames=list(experm_samples, original_names));
removed_contam_matrix=matrix(0, nrow=num_exp_samples, ncol=num_exp_cat, 
	dimnames=list(experm_samples, colnames(normalized_mat)));
selected_control_matrix=matrix(0, nrow=num_exp_samples, ncol=num_exp_cat,
	dimnames=list(experm_samples, colnames(normalized_mat)));

obs_prop_removed=numeric(num_exp_samples);
bs_prop_removed=numeric(num_exp_samples);
names(obs_prop_removed)=experm_samples;
names(bs_prop_removed)=experm_samples;

#num_exp_samples=length(exp_ids);
fits=list();

num_samp=length(experm_samples);
counter=1;

selected_controls=character();

for(exp_samp_id in experm_samples){

	cat("\nWorking on: ", exp_samp_id, " (", counter, "/", num_samp, ")\n", sep="");
	exp_dist=normalized_mat[exp_samp_id,];

	if(doPaired){
		ctl_name=pairs[exp_samp_id];
		ctl_dist=normalized_mat[ctl_name,];
	}else if(doAveraged){
		ctl_name="average control";
		ctl_dist=avg_cont_dist;
	}else if(doMatched){
		closest_name=find_closest(exp_dist, contam_mat);
		cat("Closest control was: ", closest_name, "\n");
		ctl_name=paste("closest control: ", closest_name, sep="");
		ctl_dist=contam_mat[closest_name,];
		selected_controls=c(selected_controls, closest_name);
	}else if(doMixture){
		mixture_distr_rec=find_mixture(exp_dist, contam_mat);
		ctl_dist=mixture_distr_rec[["categorical_profile"]];
		ctl_name=mixture_distr_rec[["mixture_name"]];
		mixture_component_matrix[exp_samp_id,]=mixture_distr_rec[["component_profile"]];
	}

	# Get Num Taxa:
	num_exp_cat=sum(exp_dist>0);
	num_ctl_cat=sum(ctl_dist>0);
	# Get Min Abund
	min_exp_abd=min(exp_dist[exp_dist>0]);
	min_ctl_abd=min(ctl_dist[ctl_dist>0]);


	# Observed fit
	obs_fit=fit_contaminant_mixture_model(ctl_dist, exp_dist, PLevel);
	selected_control_matrix[exp_samp_id,]=ctl_dist;

	# Observed prof of removed contaminants
	removed_contam_arr=obs_fit$removed_contam_profile;
	removed_contam_arr[is.na(removed_contam_arr)]=0;
	removed_contam_matrix[exp_samp_id,]=removed_contam_arr;

	# Bootstrap fit
	cat("Perturbing...\n");
	pert_ctrl=perturb_dist(ctl_dist, Counts, NumBS, sim_anneal=SimAnneal);
	cat("Num Perturbations: ", nrow(pert_ctrl), "\n", sep="");
	cat("Fitting...\n");
	fits=bootstrp_fit(exp_dist, pert_ctrl, PLevel);

	#print(quantile(fits$stats[,"removed"]));
	perc95_ix=min(which(fits$stats[,"removed"]==quantile(fits$stats[,"removed"], 
		1-Quantile, type=1)));

	# Save cleaned to matrix for export
	cleaned_obs_matrix[exp_samp_id,]=obs_fit$cleaned;
	cleaned_bs_matrix[exp_samp_id,]=fits$cleaned[perc95_ix,];
	obs_prop_removed[exp_samp_id]=obs_fit$removed;
	bs_prop_removed[exp_samp_id]=fits$stats[perc95_ix, "removed"];

	# Get the max abundance expect across all fits
	max_abund=max(exp_dist, ctl_dist, obs_fit$cleaned, pert_ctrl[perc95_ix,], 
		fits$cleaned[perc95_ix,]);
	max_disp_y=max_abund*1.1;
	
	# 1.) Plot obs remove
	plot_RAcurve(exp_dist, title=paste("Obs. Experimental.:", exp_samp_id), 
		top=top_cat_to_plot, max=max_disp_y);
	mtext(paste("Num Categories:", num_exp_cat), line=-1.75, outer=F, cex=.5);
	mtext(paste("Min Abundance:", min_exp_abd), line=-2.5, outer=F, cex=.5);

	# 2.) Plot obs ctrl
	plot_RAcurve(ctl_dist, title=paste("Obs. Control:", ctl_name), top=top_cat_to_plot, 
		max=max_disp_y, overlay_dist=exp_dist);
	mtext(paste("Num Categories:", num_ctl_cat), line=-1.75, outer=F, cex=.5);
	mtext(paste("Min Abundance:", min_ctl_abd), line=-2.5, outer=F, cex=.5);

	# 3.) Plot straight obs filter
	plot_RAcurve(obs_fit$cleaned, title=paste("Obs. Cleaned:"), top=top_cat_to_plot, 
		max=max_disp_y, overlay_dist=exp_dist);
	mtext(paste("Proportion Removed:", round(obs_fit$removed, 3)), line=-1.75, outer=F, cex=.5);
	mtext(paste("Multiplier:", round(obs_fit$multiplier, 3)), line=-2.5, outer=F, cex=.5);

	# 4.) Plot perturbation instance at 95% best 
	plot_RAcurve(pert_ctrl[perc95_ix,], 
		title=paste(percentile_str, " Most Removed Pert. Control Instance", sep=""), 
		top=top_cat_to_plot, max=max_disp_y, overlay_dist=exp_dist);
	# 5.) Plot filtered instance at 95% best
	plot_RAcurve(fits$cleaned[perc95_ix,], 
		title=paste(percentile_str, " Most Removed Cleaned", sep=""), 
		top=top_cat_to_plot, max=max_disp_y, overlay_dist=exp_dist);
	mtext(paste("Proportion Removed:", round(fits$stats[perc95_ix, "removed"], 3)), 
		line=-1.75, outer=F, cex=.5);
	mtext(paste("Multiplier:", round(fits$stats[perc95_ix, "multiplier"], 3)), 
		line=-2.5, outer=F, cex=.5);

	# 6.) Plot range of pertubation
	plot_RAbox(pert_ctrl, title="Perturbed Control", top=top_cat_to_plot, max=max_disp_y);

	# 7.) Plot range of filtered
	plot_RAbox(fits$cleaned, title="Range of Cleaned", top=top_cat_to_plot, max=max_disp_y);

	# 8.) Plot histogram of percent removed
	hist(fits$stat[,"removed"], main="Bootstrapped Proportions Removed", 
		xlab="Bootstrapped Proportions Removed", 
		breaks=seq(0,1,.025), xlim=c(0,1));
	abline(v=fits$stat[perc95_ix,"removed"], col="blue");
	
	# 9.) Plot histogram of multiplier
	hist(fits$stat[,"proportion"], main="Bootstrapped Mixture: (c/(1+c))", xlab="Multipliers",
		breaks=seq(0,1,.025), xlim=c(0,1));
	abline(v=fits$stat[perc95_ix,"proportion"], col="blue");

	counter=counter+1;
}

###############################################################################

write_summary_file=function(out_mat, fname){
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
# Output summary file table

# Adjust counts
obs_saved=1-obs_prop_removed;
bs_saved=1-bs_prop_removed;

# exp totals
exp_tot=apply(counts_mat[experm_samples,],1, sum);

obs_saved_totals=numeric(num_exp_samples);
bs_saved_totals=numeric(num_exp_samples);

names(obs_saved_totals)=experm_samples;
names(bs_saved_totals)=experm_samples;

# Adjust counts based on proportion of non-contaminant
obs_saved_totals[experm_samples]=exp_tot[experm_samples]*obs_saved[experm_samples];
bs_saved_totals[experm_samples]=exp_tot[experm_samples]*bs_saved[experm_samples];

# Compute individual category counts based on non-contam * category abundance
cleaned_obs_counts=matrix(-1, ncol=num_categories, nrow=num_exp_samples, 
	dimnames=list(experm_samples,original_names));
cleaned_bs_counts=matrix(-1, ncol=num_categories, nrow=num_exp_samples, 
	dimnames=list(experm_samples,original_names));

for(samp in experm_samples){
	cleaned_obs_counts[samp,]=floor(cleaned_obs_matrix[samp,]*obs_saved_totals[samp]);
	cleaned_bs_counts[samp,]=floor(cleaned_bs_matrix[samp,]*bs_saved_totals[samp]);
}

write_summary_file(cleaned_obs_counts, paste(OutputFileRoot, ".obsrv_cln.summary_table.tsv", sep=""));
write_summary_file(cleaned_bs_counts, paste(OutputFileRoot,  ".btstp_cln.summary_table.tsv", sep=""));

###############################################################################

# Keep track of proportion and number of reads virtually removed
obs_removed_totals=numeric(num_exp_samples);
bs_removed_totals=numeric(num_exp_samples);

names(obs_removed_totals)=experm_samples;
names(bs_removed_totals)=experm_samples;

obs_removed_totals[experm_samples]=exp_tot[experm_samples]*obs_prop_removed[experm_samples];
bs_removed_totals[experm_samples]=exp_tot[experm_samples]*bs_prop_removed[experm_samples];

obs_removed_total_median=median(obs_removed_totals);
obs_removed_prop_median=median(obs_prop_removed);
bs_removed_total_median=median(bs_removed_totals);
bs_removed_prop_median=median(bs_prop_removed);

# Plot histogram of what was removed

par(mfrow=c(2,2));

hist(obs_removed_totals, breaks=nclass.Sturges(obs_removed_totals)*2,
	 main="Observed Reads Removed (Counts)", xlab="Number of Reads");
abline(v=obs_removed_total_median, col="blue");
text(obs_removed_total_median, 0, adj=c(-.4,-.1), srt=90, 
	sprintf("Median=%g", obs_removed_total_median), col="blue", font=2);

hist(obs_prop_removed*100, breaks=nclass.Sturges(obs_prop_removed)*2, 
	main="Observed Reads Removed (Percentage)", xlab="Percent of Reads");
abline(v=obs_removed_prop_median*100, col="blue");
text(obs_removed_prop_median*100, 0, adj=c(-.4,-.1), srt=90,
	sprintf("Median=%3.2f", 100* obs_removed_prop_median), col="blue", font=2);

hist(log10(obs_removed_totals), breaks=nclass.Sturges(log10(obs_removed_totals))*2,
	 main="Observed Reads Removed Log10(Counts)", xlab="Log10(Number of Reads)");
abline(v=log10(obs_removed_total_median), col="blue");
text(log10(obs_removed_total_median), 0, adj=c(-.4,-.1), srt=90, 
	sprintf("Median=%g", log10(obs_removed_total_median)), col="blue", font=2);

plot(log10(obs_removed_totals), obs_prop_removed*100, 
	main="Log10 Counts Removed vs. Percentage Removed", xlab="Log10(Counts)", ylab="Percent Removed");



hist(bs_removed_totals, breaks=nclass.Sturges(bs_removed_totals)*2, 
	main="Bootstrapped Reads Removed (Counts)", xlab="Number of Reads");
abline(v=bs_removed_total_median, col="blue");
text(bs_removed_total_median, 0, adj=c(-.4,-.1), srt=90,
	sprintf("Median=%g", bs_removed_total_median), col="blue", font=2);

hist(bs_prop_removed*100, breaks=nclass.Sturges(bs_prop_removed)*2, 
	main="Bootstrapped Reads Removed (Percentage)", xlab="Percentage of Reads");
abline(v=bs_removed_prop_median*100, col="blue");
text(bs_removed_prop_median*100, 0, adj=c(-.4,-.1), srt=90,
	sprintf("Median=%3.2f", 100* bs_removed_prop_median), col="blue", font=2);

hist(log10(bs_removed_totals), breaks=nclass.Sturges(log10(bs_removed_totals))*2, 
	main="Bootstrapped Reads Removed Log10(Counts)", xlab="Log10(Number of Reads)");
abline(v=log10(bs_removed_total_median), col="blue");
text(log10(bs_removed_total_median), 0, adj=c(-.4,-.1), srt=90,
	sprintf("Median=%g", log10(bs_removed_total_median)), col="blue", font=2);

plot(log10(bs_removed_totals), bs_prop_removed*100, 
	main="Log10 Counts Removed vs. Percentage Removed", xlab="Log10(Counts)", ylab="Percent Removed");

###############################################################################

par(mfrow=c(4,1));

# Plot the average of what the contamination profile looked like for each sample
mean_removed_arr=apply(removed_contam_matrix, 2, function(x){mean(x, na.rm=T);});
mean_removed_ordered_ix=order(mean_removed_arr, decreasing=T);

mean_removed_ordered_arr=mean_removed_arr[mean_removed_ordered_ix];
num_cat_wmean_gt0=sum(mean_removed_ordered_arr>0);

plot_RAcurve(mean_removed_ordered_arr, 
	title="Mean Contaminant Profile",
	top=min(top_cat_to_plot, num_cat_wmean_gt0));

# Plot the median of what the contamination profile looked like for each sample
median_removed_arr=apply(removed_contam_matrix, 2, function(x){median(x, na.rm=T);});
median_removed_ordered_ix=order(median_removed_arr, decreasing=T);

median_removed_ordered_arr=median_removed_arr[median_removed_ordered_ix];
plot_RAcurve(median_removed_ordered_arr, 
	title="Median Contaminant Profile",
	top=min(top_cat_to_plot, num_cat_wmean_gt0));

# Plot scatter
plot_RAstrip(removed_contam_matrix, title="Contaminant Strip Plot (Proportion of Contaminants)", 
	top=min(top_cat_to_plot, num_cat_wmean_gt0));

# Plot proportion any removed
plot_RAprop(removed_contam_matrix, title="Proportion of Samples (With Any Amount of Category Removed)", 
	top=min(top_cat_to_plot, num_cat_wmean_gt0));

###############################################################################
# Plot the original, filtered, controls used, and removed in the order of the original abundance

top_to_plot=50;

mean_experm_samples=apply(normalized_mat[experm_samples,], 2, mean);
mean_selected_control=apply(selected_control_matrix, 2, mean);
mean_cleaned_obs=apply(cleaned_obs_matrix, 2, mean);
# removed_contam_matrix, calculated above

order_by_cleanobs_ix=order(mean_cleaned_obs, decreasing=T);

max_in_top=max(c(
	mean_experm_samples, mean_selected_control,
	mean_cleaned_obs, mean_removed_arr));

ylimit=max_in_top*1.10;

# Plot the average input/original profile
plot_RAcurve(mean_experm_samples[order_by_cleanobs_ix], 
	title="Mean of Original Experimental Samples", 
	top=top_to_plot, max=ylimit); 

# Plot the average control that was selected
plot_RAcurve(mean_selected_control[order_by_cleanobs_ix], 
	title="Mean of Selected Controls", 
	top=top_to_plot, max=ylimit); 

# Plot the average output/filtered profile
plot_RAcurve(mean_removed_arr[order_by_cleanobs_ix], 
	title="Mean of Removed Contaminants (Obs)", 
	top=top_to_plot, max=ylimit); 

# Plot the removed
names(mean_cleaned_obs)=names(mean_removed_arr);
plot_RAcurve(mean_cleaned_obs[order_by_cleanobs_ix], 
	title="Mean of Filtered Experimental Samples (Obs)", 
	top=top_to_plot, max=ylimit); 


###############################################################################

plot_category_specific_contam_stats=function(prop_removed_mat, sample_counts, num_to_plot){

	samp_ids=names(sample_counts);
	categories=colnames(prop_removed_mat);

	par(mfrow=c(5,2));
	par(mar=c(4.1,4.1,5,1));

	hist_breaks=seq(0,1,length.out=20+1);

	for(i in 1:num_to_plot){

		plot(
			log10(sample_counts+1),
			prop_removed_mat[samp_ids,i],
			main=paste(categories[i], ":\nSample Depth vs. Proportion of Contaminants", 
				sep=""),
			ylim=c(0,1),
			xlab="Log10(Sample Depth)",
			ylab="Proportion of Contaminants"
		);

		hist(prop_removed_mat[samp_ids,i], breaks=hist_breaks,
			main=paste(categories[i], ":\nFrequency of Proportion of Contaminants", sep=""),
			xlab="Proportions of Contaminants", ylab="Num Samples");
	}

}

plot_category_specific_contam_stats(removed_contam_matrix, exp_tot, top_cat_to_plot);

###############################################################################
# Calculate correlation across taxa identified as contamination

contam_correl_mat=cor(removed_contam_matrix[,1:top_cat_to_plot]);
contam_correl_mat[is.na(contam_correl_mat)]=0;
paint_matrix(contam_correl_mat, title="Correlation Among Contamination Categories",
	plot_min=-1, plot_max=1, deci_pts=2, value.cex=.6, label_zeros=F);

paint_matrix(contam_correl_mat, title="Correlation Among Contamination Categories",
	plot_min=-1, plot_max=1, deci_pts=2, value.cex=.6, label_zeros=F, plot_row_dendr=T);



###############################################################################

# Generate a histogram of the frequency by which controls have been matched
if(doMatched){
	par(mfrow=c(1,1));
	par(mar=c(4,10,4,1));
	matched_profile=table(selected_controls);
	print(matched_profile);
	sort_ix=order(matched_profile, decreasing=T);
	matched_profile=matched_profile[sort_ix];
	barplot(matched_profile, main="Frequency of Control Matching", 
		horiz=T, las=1, cex.names=.5);

}

if(doMixture){
	formatted_mat=apply(mixture_component_matrix, c(1,2), function(x){sprintf("%3.3f", x)});
	print(formatted_mat, quote=F, width=200);

	paint_matrix(mixture_component_matrix, title="Mixture Contributions",
		plot_min=0, plot_max=1, high_is_hot=T, deci_pts=2, label_zeros=F);

}

###############################################################################


# Write parameters
fc=file(paste(OutputFileRoot, ".cleaned.stats.tsv", sep=""), "w");

cat(file=fc, "Input Summary File Table: ", InputSummaryTable, "\n", sep="");
cat(file=fc, "Output Filename Root: ", OutputFileRoot, "\n", sep="");
cat(file=fc, "P-norm level: ", PLevel, "\n", sep="");
cat(file=fc, "Num Bootstraps: ", NumBS, "\n", sep="");
cat(file=fc, "Counts: ", Counts, "\n", sep="");
cat(file=fc, "Perform 'simulated annealing':", SimAnneal, "\n", sep="");
cat(file=fc, "Quantile (1-alpha), alpha: ", Quantile, "\n", sep="");

close(fc);

###############################################################################

cat("Done.\n\n")
if(!is.null(warnings())){
	print(warnings());
}

q(status=0)
