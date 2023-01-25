#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');

options(useFancyQuotes=F);
options(digits=5);
options(width=300);

params=c(
	"tmr_results", "r", 1, "character",
	"output_root", "o", 1, "character",
	"pval_cutoff", "p", 2, "character"
);

PVAL_CUTOFF=0.1;

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

script_path=paste(head(strsplit(script_name, "/")[[1]], -1), collapse="/");

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-r <Results from TMR analysis>\n",	
	"	-o <Output Root>\n",
	"	[-p <p-value cutoff, default=", PVAL_CUTOFF, ">]\n",
	"\n",
	"This script will read in the results from a TMR analysis\n",
	"that has been run with the repeated measures identifier\n",
	"and generate longitudinal plots.\n",
	"\n",
	"\n", sep="");

if(
	!length(opt$tmr_results) || 
	!length(opt$output_root) 
){
	cat(usage);
	q(status=-1);
}

TMR_ResultsFname=opt$tmr_results;
OutputFnameRoot=opt$output_root;
PValCutoff=PVAL_CUTOFF;

if(length(opt$pval_cutoff)){
	PValCutoff=opt$pval_cutoff;
}

OutputFnameRoot=paste(OutputFnameRoot, ".p", sprintf("%.4f", as.numeric(PValCutoff)), sep="");

param_text=capture.output({
	cat("\n");
	cat("TMR Results Input File: ", TMR_ResultsFname, "\n");
	cat("Output File Name Root: ", OutputFnameRoot, "\n");
	cat("P-value Cutoff: ", PValCutoff, "\n");
	cat("\n");
});

cat(paste(param_text, collapse="\n"), "\n");

###############################################################################

plot_title_page=function(title, subtitle=""){

        orig.par=par(no.readonly=T);
        par(family="serif");
        par(mfrow=c(1,1));

        plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n",  xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                );

        # Title
        title_cex=3;
        title_line=1;
        text(0.5, title_line, title, cex=title_cex, font=2, adj=c(.5,1));

        # Subtitle
        num_subt_lines=length(subtitle);
        cxy=par()$cxy;
        for(i in 1:num_subt_lines){
                text(.5, title_line -title_cex*cxy[2] -i*cxy[2], subtitle[i], adj=.5);
        }

        par(orig.par);
}


###############################################################################

load_tmr_results=function(fname, pval_cutoff=0.01){

	cat("Loading: ", fname, " with cutoff at: ", pval_cutoff, "\n", sep="");
	longit_tmr_table=as.data.frame(read.delim(fname,  header=TRUE, row.names=NULL, check.names=FALSE, sep="\t"));
	#print(longit_tmr_table);

	model_types=sort(unique(longit_tmr_table[,"model_type"]));

	cat("Model Types:\n");
	print(model_types);

	signf_header=c("model_type", "model_name", "predictor", "response", "min_pval");
	significant_hits=matrix(character(), nrow=0, ncol=length(signf_header));
	colnames(significant_hits)=signf_header;

	longit_tmr_rec=list();
	
	for(mt in model_types){
		cat("Extracting Model Type: ", mt, "\n");
		mt_ix=longit_tmr_table[,"model_type"]==mt;
		mt_tab=longit_tmr_table[mt_ix,];
		#print(mt_tab);	
		longit_tmr_rec[[mt]]=list();

		model_names=sort(unique(mt_tab[,"model_name"]));
		for(mn in model_names){
			cat("  Extracting Model Name: ", mn, "\n");
			mn_ix=mt_tab[,"model_name"]==mn;
			mn_tab=mt_tab[mn_ix,];
			#print(mn_tab);
			longit_tmr_rec[[mt]][[mn]]=list();

			predictors=sort(unique(mn_tab[, "predictor"]));
			for(prd in predictors){
				cat("   Extracting Predictors:", prd, "\n");
				prd_ix=mn_tab[, "predictor"]==prd;
				prd_tab=mn_tab[prd_ix,];
				#print(tab);
				longit_tmr_rec[[mt]][[mn]][[prd]]=list();

				responses=sort(unique(prd_tab[, "response"]));
				for(rsp in responses){
					cat("    Extracting Responses:", rsp, "\n");
					rsp_ix=prd_tab[, "response"]==rsp;
					rsp_tab=prd_tab[rsp_ix,];

					offset=as.numeric(rsp_tab[, "RepeatedMeasureAnalysisID"]);
					order_ix=order(offset, decreasing=F);
					#print(rsp_tab[order_ix,]);
					longit_tmr_rec[[mt]][[mn]][[prd]][[rsp]]=rsp_tab[order_ix,];

					min_pvals=min(rsp_tab[,"pval"]);
					if(min_pvals<=pval_cutoff){
						significant_hits=rbind(significant_hits, 
							c(mt, mn, prd, rsp, min_pvals));
					}
				}
			}
		}
	}

	#print(longit_tmr_rec);
	#print(significant_hits);
	results=list();
	results[["tables"]]=longit_tmr_rec;
	results[["signif_index"]]=significant_hits;
	results[["pvalue_cutoff"]]=pval_cutoff;

	return(results);
}

#-----------------------------------------------------------------------------#


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

		text_size=max(.01, min(.8, .8 - .003*(lines_pp-52)));
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

#-----------------------------------------------------------------------------#

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, deci_pts=4,
        label_zeros=T, counts=F, value.cex=2,
        plot_col_dendr=F,
        plot_row_dendr=F
){

        num_row=nrow(mat);
        num_col=ncol(mat);

        row_names=rownames(mat);
        col_names=colnames(mat);

        orig.par=par(no.readonly=T);

        cat("Painting Matrix: ", title, "\n");
        cat("Num Rows: ", num_row, "\n");
        cat("Num Cols: ", num_col, "\n");


        if(num_row==0 || num_col==0){
                plot(0, type="n", xlim=c(-1,1), ylim=c(-1,1), xaxt="n", yaxt="n", bty="n", xlab="", ylab="",
                        main=title);
                text(0,0, "No data to plot...");
                return();
        }

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
                plot_min=min(mat);
        }
        if(is.na(plot_max)){
                plot_max=max(mat);
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

        if(num_row==1){
                plot_row_dendr=F;
        }
        if(num_col==1){
                plot_col_dendr=F;
        }

        # Don't plot dendrogram if there are any NAs in the matrix
        if(any(is.na(mat))){
                plot_col_dendr=F;
                plot_row_dendr=F;
        }

        if(plot_col_dendr && plot_row_dendr){
                layoutmat=matrix(
                        c(
                        rep(c(rep(4, row_dend_width), rep(3, heatmap_width)), col_dend_height),
                        rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height)
                        ), byrow=T, ncol=row_dend_width+heatmap_width);

                col_dendr=get_dendrogram(mat, type="col");
                row_dendr=get_dendrogram(mat, type="row");

                mat=mat[row_dendr[["names"]], col_dendr[["names"]], drop=F];

        }else if(plot_col_dendr){
                layoutmat=matrix(
                        c(
                        rep(rep(2, heatmap_width), col_dend_height),
                        rep(rep(1, heatmap_width), heatmap_height)
                        ), byrow=T, ncol=heatmap_width);

                col_dendr=get_dendrogram(mat, type="col");
                mat=mat[, col_dendr[["names"]], drop=F];

        }else if(plot_row_dendr){
                layoutmat=matrix(
                        rep(c(rep(2, row_dend_width), rep(1, heatmap_width)), heatmap_height),
                        byrow=T, ncol=row_dend_width+heatmap_width);

                row_dendr=get_dendrogram(mat, type="row");
                mat=mat[row_dendr[["names"]],,drop=F];
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
        plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), xaxt="n", yaxt="n", bty="n", xlab="", ylab="");
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

                        if(mat[y,x]!=0 || label_zeros){
                                if(counts){
                                        text_lab=sprintf("%i", mat[y,x]);
                                }else{
                                        text_lab=sprintf(paste("%0.", deci_pts, "f", sep=""), mat[y,x]);
                                }
                                text(x-.5, y-.5, text_lab, srt=atan(num_col/num_row)/pi*180, cex=value.cex, font=2);
                        }
                }
        }

        ##################################################################################################

        par(mar=c(0, 0, 0, 0));

        if(plot_row_dendr && plot_col_dendr){
                rdh=attributes(row_dendr[["tree"]])$height;
                cdh=attributes(col_dendr[["tree"]])$height;
                plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", bty="n", xlim=c(rdh, 0));
                plot(col_dendr[["tree"]], leaflab="none",xaxt="n", yaxt="n", bty="n", ylim=c(0, cdh));
                plot(0,0, type="n", bty="n", xaxt="n", yaxt="n");
                #text(0,0, "Placeholder");
        }else if(plot_row_dendr){
                rdh=attributes(row_dendr[["tree"]])$height;
                plot(row_dendr[["tree"]], leaflab="none", horiz=T, xaxt="n", yaxt="n", bty="n", xlim=c(rdh, 0));
                #text(0,0, "Row Dendrogram");
        }else if(plot_col_dendr){
                cdh=attributes(col_dendr[["tree"]])$height;
                plot(col_dendr[["tree"]], leaflab="none", xaxt="n", yaxt="n", bty="n", ylim=c(0, cdh));
                #text(0,0, "Col Dendrogram");
        }

        par(orig.par);

}

##############################################################################
# Main Program Starts Here!
##############################################################################

pdf(paste(OutputFnameRoot, ".tmr_longitudinal.pdf", sep=""), height=8.5, width=11);
plot_text(param_text);

# Load TMR Results 
cat("Loading TMR Results", TMR_ResultsFname, "\n");
loaded_results=load_tmr_results(TMR_ResultsFname, PValCutoff);

##############################################################################

tally_by_model_name=function(tally_in, tab, pval_cutoff){

	num_offsets=nrow(tab);
	offsets=as.character(tab[,"RepeatedMeasureAnalysisID"]);

	empty_rec=c(0,0,0);
	names(empty_rec)=c("Total", "Negative", "Positive");

	if(length(tally_in)==0){
		for(i in 1:num_offsets){
			tally_in[[offsets[i]]]=empty_rec;
		}
	}

	for(i in 1:num_offsets){
		if(tab[i,"pval"]<pval_cutoff){
			if(tab[i,"coef"]<0){
				current=c(1,1,0);
			}else{
				current=c(1,0,1);
			}
			tally_in[[offsets[i]]]=tally_in[[offsets[i]]]+current;
		}
	}

	
	return(tally_in);
}

plot_tally=function(tally_in, title){
	
	offsets_char=names(tally_in);
	offsets=as.numeric(offsets);
	num_offsets=length(offsets_char);

	tally_matrix=matrix(NA, ncol=3, nrow=num_offsets);
	colnames(tally_matrix)=c("Total", "Negative", "Positive");
	for(i in 1:num_offsets){
		tally_matrix[i,]=tally_in[[offsets_char[i]]];	
	}
	#print(tally_matrix);
	max_count=max(tally_matrix[,"Total"]);

	# Stacked barplot
	plot(offsets, tally_matrix[,"Total"], main=title, 
		ylim=c(0, max_count*1.1),
		type="h", lend=1,
		col="green",
		xlab="offset", ylab="Number of Associations",
		lwd=10, las=1);
	title(line=.5, main="Significant Associations", 
		font.main=3, cex.main=.8, col.main="black");
	points(offsets, tally_matrix[,"Total"], type="p", col="black", pch="-");
	points(offsets, tally_matrix[,"Negative"], type="h", lend=1, col="red", lwd=10);

	# Connnected dots
	plot(offsets, tally_matrix[,"Total"], main=title,
		ylim=c(0, max_count*1.1),
		type="b", xlab="offset", ylab="Number of Associations",
		lwd=1, las=1);
	title(line=.5, main="Significant Associations", 
		font.main=3, cex.main=.8, col.main="black");

}

##############################################################################

plot_title_page("Longitudinal Accumulation of TMR Results", c(
	"The results from this accumulation are organized hierarchically.",
	"There could be (if the experiment supports a complete TMR model) up to 3 main sections.", 
	"1.) Treatments/Covariates-to-Measured",
	"2.) Measured-to-Measured",
	"3.) Measured-to-Response", 
	"",
	"Within each model, the compartment-to-compartment associations are plotted.",
	"",
	"The left-hand-side plot illustrate coefficients and the right-hand-side plot illustrate",
	"p-values for the same intercompartment associations, over time.",
	"",
	"The coefficients plot are annotated with a horizontal dashed line at 0.",
	"The p-values plot are annotated with multiple dashed lines at various significance levels.",
	"The p-values plot are negative log10-scaled with with the -log10 values labeled on the",
	"left y-axis, and the original untransformed p-values on the right.",
	"In both the coefficient and p-value plots, when a timepoint yields a significant association",
	"the point is highlighted with a glyph.  A red or green glyph is used to signify whether the",
	"association is negative or positive, respectively.",
	"",
	"The name of the predictor [P] variable and response [R] variable are annotated above each plot.",
	"",
	"At the end of each compartment's plot, there are two summary plots.  The first summary plot",
	"is a barplot with the number of significant assocations over time.  The number of positive and",
	"negative associations at each time point are also differentiated by green and red coloring,",
	"respectively.  The second summary plot is a simple line plot without any colors."
));

signif_index=loaded_results[["signif_index"]];
tables=loaded_results[["tables"]];

num_signf=nrow(signif_index);

pval_cutoffs=c(1, .1, .05, .01, .001);
neglog10_cutoffs=-log10(pval_cutoffs);

overall_tally=list();

model_types=unique(signif_index[,"model_type"]);
for(mt in model_types){
	plot_title_page(paste("Model Type:\n", mt, sep=""));	
	par(mfrow=c(4,2));
	par(las=2);

	mt_ix=signif_index[,"model_type"]==mt;

	table_by_model_type=signif_index[mt_ix,];

	model_names=sort(unique(table_by_model_type[,"model_name"]));

	overall_tally[[mt]]=list();

	max_signf_assoc=0;

	for(mn in model_names){
		mn_ix=table_by_model_type[,"model_name"]==mn;
		table_by_model_name=table_by_model_type[mn_ix,,drop=F];

		# Section divider by model name
		par(mar=rep(0,4));
		plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n",  xaxt="n", yaxt="n",
			xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
		);
		text(.5, .5, mn, pos=1, cex=2, font=2);
		plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n",  xaxt="n", yaxt="n",
			xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
		);
		text(.5, .5, mt, pos=1, cex=2, font=2, col="grey");

		#print(table_by_model_name);

		by_model_name_tally=list();
		par(mar=c(3,4,4,4));
		for(ix in 1:nrow(table_by_model_name)){
			signf_keys=table_by_model_name[ix,];
			
			t=signf_keys["model_type"];
			n=signf_keys["model_name"];
			p=signf_keys["predictor"];
			r=signf_keys["response"];
			minitab=tables[[t]][[n]][[p]][[r]];

			#print(minitab);
			by_model_name_tally=tally_by_model_name(by_model_name_tally, minitab, PValCutoff);

			offsets=minitab[,"RepeatedMeasureAnalysisID"];
			coef=minitab[,"coef"];
			pval=minitab[,"pval"];

			model_typename=paste("type: ", t, "    name: ", n, sep="");
			compare_name=paste( "[P] ", p, " -> [R] ", r, sep="");

			coef_max_mag=max(abs(range(coef)));

			# Make significant glyphs look different
			above_cutoff=pval<PValCutoff;
			pos_assoc=coef>0;

			num_offsets=length(offsets);
			pch=rep(1, num_offsets);
			pch[above_cutoff]=22;
			bg=rep(1, num_offsets);
			bg[above_cutoff & pos_assoc]="green";
			bg[above_cutoff & !pos_assoc]="red";
			cex=rep(1, num_offsets);
			cex[above_cutoff]=1.7;
				

			# Plot coefficients
			plot(offsets, coef, ylim=c(-coef_max_mag*1.1, coef_max_mag*1.1),
				pch=pch, bg=bg, cex=cex, 
				las=1,
				type="b",
				ylab="Coefficient", xlab="Time", main="");
			abline(h=0, col="blue", lty="dashed");

			title(line=2.5, main=model_typename, cex.main=.9, font.main=1, col.main="grey25");
			title(line=.5, main="Associations / Coefficients", 
				font.main=3, cex.main=.8, col.main="black");
			title(line=1.5, main=compare_name, cex.main=1.1, font.main=2, col.main="black");

			
			# Plot p-values
			cutoffs_range=range(c(-log10(min(pval)), neglog10_cutoffs));	
			plot(offsets, -log10(pval), 
				pch=pch, bg=bg, cex=cex, 
				las=1,
				ylim=cutoffs_range,
				type="b",
				ylab="-Log10(p-value)", xlab="Time", main="", cex.axis=1);
			axis(side=4, at=neglog10_cutoffs, label=pval_cutoffs, las=2, cex.axis=.8);
			abline(h=setdiff(neglog10_cutoffs, 0), col="blue", lty="dashed");

			title(line=2.5, main=model_typename, cex.main=.9, font.main=1, col.main="grey25");
			title(line=.5, main="Significances / P-values", 
				font.main=3, cex.main=.8, col.main="black");
			title(line=1.5, main=compare_name, cex.main=1.1, font.main=2, col.main="black");
		}

		#cat("Model Name: ", mn, "\n");
		#print(by_model_name_tally);
		# Plot tallies
		plot_tally(by_model_name_tally, mn);
		overall_tally[[mt]][[mn]]=by_model_name_tally;

	}

}

plot_title_page("Overall Summaries of Significance Counts", c(
	"These plots recapitulation the barplots that were previously provided at the end of each",
	"compartment summary.  These plots have a common scale and are grouped together so they can be",
	"more readily compared.  Plots are grouped by model type.",
	"",
	"These plots provide a quick means to identify the timepoints when some compartments may be",
	"more predictive than others."
));

plot_overall_tally=function(tally_rec){

	converted=list();
	max_total=0;
	max_offset=0;

	model_types=names(tally_rec);
	for(mt in model_types){
		cat("Model Type: ", mt, "\n");
		model_names=names(tally_rec[[mt]]);

		converted[[mt]]=list();

		for(mn in model_names){

			acc=matrix(NA, nrow=0, ncol=4);
			colnames(acc)=c("Offset", "Total", "Negative", "Positive");

			cat("Model Name: ", mn, "\n");
			offsets=names(tally_rec[[mt]][[mn]]);
			for(offset in offsets){

				print(tally_rec[[mt]][[mn]][[offset]]);
				acc=rbind(acc, c(as.numeric(offset), tally_rec[[mt]][[mn]][[offset]]));
			}

			converted[[mt]][[mn]]=acc;
			max_total=max(max_total, acc[,"Total"]);
			max_offset=max(max_offset, acc[,"Offset"]);

		}
	}		

	#print(converted);
	#print(max_total);

	for(mt in model_types){
		
		
		model_names=names(tally_rec[[mt]]);
		num_models=length(model_names);
		cat("Num Models in ", mt, ": ", num_models, "\n");

		par(mfrow=c(5, 1));
		par(mar=c(2.5, 4, 2,.5));

		# Place holder model type
		plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n",  xaxt="n", yaxt="n",
                        xlab="", ylab="", bty="n", mar=c(0,0,0,0));
                text(.5, .5, mt, pos=1, cex=2, font=2);	

		for(mn in model_names){

			signf_tab=converted[[mt]][[mn]]

			plot(0, type="n", xlim=c(0, max_offset), ylim=c(0, max_total*1.1),
				ylab="Counts", xlab="Offset",
				main="", las=1);
			title(main=mn, line=.5);

			# lines connecting totals
			points(signf_tab[,"Offset"], signf_tab[,"Total"], 
				type="l", col="blue", lty="dotted", lwd=2, cex=2);

			# Positive Prop
			points(signf_tab[,"Offset"], signf_tab[,"Total"], 
				type="h", col="green", lwd=20, lend=1);
			# Negative Prop
			points(signf_tab[,"Offset"], signf_tab[,"Negative"], 
				type="h", col="red", lwd=20, lend=1);
			# Glyph at totals
			points(signf_tab[,"Offset"], signf_tab[,"Total"], 
				type="p", pch="-", col="black", cex=2);

		}
		

	}

}


plot_overall_tally(overall_tally);

##############################################################################

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
