#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);

params=c(
	"summary_file", "s", 1, "character",
	"factors", "f", 1, "character",
	"outputroot", "o", 2, "character",
	"trim_names", "t", 2, "logical"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table>\n",
	"	-f <factors>\n",
	"	[-t (trim long names flag)]\n",
	"	[-o <output filename root>]\n",
	"\n",
	"This script will read in the summary file table and the factors\n",
	"file and generate a page for each factor column in the summary table.\n",
	"In each page, for each factor, a compositional barplot will be\n",
	"generated for each factor level.\n",
	"\n",
	"For each factor level, all the samples at that level will be averaged\n",
	"together.\n",
	"\n");

if(!length(opt$summary_file) || !length(opt$factors)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputRoot=gsub(".summary_table.xls", "", opt$summary_file);
	OutputRoot=gsub(".summary_table.tsv", "", OutputRoot);
}else{
	OutputRoot=opt$outputroot;
}

TrimLongNames=length(opt$trim_names)>0;

SummaryFile=opt$summary_file;
FactorsFile=opt$factors;

cat("\n");
cat("Summary File : ", SummaryFile, "\n", sep="");
cat("Factors File: ", FactorsFile, "\n", sep="");
cat("Output File: ", OutputRoot, "\n", sep="");
cat("Trim Long Names: ", TrimLongNames, "\n", sep="");
cat("\n");

##############################################################################

load_factors=function(fname){
	cat("Loading Factors: ", fname, "\n");
	factors=data.frame(read.table(fname,  header=TRUE, row.names=1, check.names=FALSE, comment.char="", quote=""));
	factor_names=colnames(factors);

	ignore_idx=grep("^IGNORE\\.", factor_names);

	if(length(ignore_idx)!=0){
		return(factors[-ignore_idx]);
	}else{
		return(factors);
	}
}

#------------------------------------------------------------------------------

load_summary_file=function(fname){
	cat("Loading Summary Table: ", fname, "\n");
	inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, comment.char="", row.names=1))
	counts_mat=inmat[,2:(ncol(inmat))];
	return(counts_mat);
}

trim_long_names=function(names){
	new_names=sapply(names, function(x){tail(strsplit(x, " ")[[1]],1)});
	return(new_names);
}

#------------------------------------------------------------------------------

normalize=function(counts){
	totals=apply(counts, 1, sum);
	num_samples=nrow(counts);
	normalized=matrix(0, nrow=nrow(counts), ncol=ncol(counts));

	for(i in 1:num_samples){
		normalized[i,]=counts[i,]/totals[i];
	}
	
	colnames(normalized)=colnames(counts);
	rownames(normalized)=rownames(counts);	
	return(normalized);
}

#------------------------------------------------------------------------------

plot_text=function(strings){
	par(mfrow=c(1,1));
        par(family="Courier");
        par(oma=rep(.5,4));
        par(mar=rep(0,4));

        num_lines=length(strings);

        top=max(as.integer(num_lines), 52);

        plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                );

        norm=sorted/sum(x);
        n=length(norm);
        tail=0;
        for(i in 1:n){
                tail=tail + norm[i]*((i-1)^2);
        }
        return(sqrt(tail));
}

#------------------------------------------------------------------------------

generate_colors=function(num_col){

	if(num_col<=64){inc=c(4,4,4);}
	if(num_col<=48){inc=c(4,4,3);}
	if(num_col<=36){inc=c(4,3,3);}
	if(num_col<=27){inc=c(3,3,3);}
	if(num_col<=18){inc=c(3,3,2);}
	if(num_col<=12){inc=c(3,2,2);}
	if(num_col<=8) {inc=c(2,2,2);}
	
	r_inc=seq(0,1,length.out=inc[1]);
	g_inc=seq(0,1,length.out=inc[2]);
	b_inc=seq(0,1,length.out=inc[3]);

	num_colors=prod(inc);
	colors=character(num_colors);

	i=1;
	for(r in r_inc){
		for(g in g_inc){
			for(b in b_inc){
				colors[i]=rgb(r,b,g);
				i=i+1;
			}
		}
	}

	cat("Generated colors: ", num_colors, "\n");
	cat("Targeted colors: ", num_col, "\n");

	colors=(colors[seq(1, num_colors, num_colors/num_col)]);

	return(sample(colors, replace=F));

}

plot_compositions=function(sample_means, num_samples, title="", top=10){

	num_fact_level=ncol(sample_means);
	num_categories=nrow(sample_means);

	factor_levels=colnames(sample_means);
	category_names=rownames(sample_means);

	# Get top 10 for each factor level
	cumulative_top=character();
	for(fact_level in factor_levels){
		order_ix=order(sample_means[,fact_level, drop=F], decreasing=T);
		sorted_top=sample_means[order_ix, fact_level, drop=F];
		cumulative_top=c(cumulative_top, rownames(sorted_top)[1:top]);
	}
	cumulative_top=unique(cumulative_top);
	cat("Cumulative Top ", top, ":\n");
	print(cumulative_top);

	# Get top to plot
	top_matrix=sample_means[cumulative_top, , drop=F];

	# Compute and insert Average column
	avg=apply(top_matrix, 1, mean);
	top_matrix=cbind(top_matrix, avg);
	colnames(top_matrix)=c(factor_levels, "Average");
	factor_levels=colnames(top_matrix);
	
	# Sort by average
	avg_ord=order(top_matrix[,"Average"], decreasing=T);
	top_matrix=top_matrix[avg_ord,];
	nsamp_names=names(num_samples);
	num_samples=c(num_samples, sum(num_samples));
	names(num_samples)=c(nsamp_names, "Average");

	# Compute and insert Remaining row
	remaining=apply(top_matrix, 2, function(x){1-sum(x)});
	top_categories=rownames(top_matrix);
	top_matrix=rbind(top_matrix, remaining);
	rownames(top_matrix)=c(top_categories, "Remaining");
	top_categories=rownames(top_matrix);

	print(top_matrix);

	colors=generate_colors(nrow(top_matrix));

	plot_sbs=function(tm, num_samples, bar_width=1, bar_spacing=.3, category_label_spacing=4, title=""){

		num_bars=ncol(tm);
		num_cat=nrow(tm);
		factor_levels=colnames(tm);

		pos_of_avg=which(factor_levels=="Average");

		cat("Num bars to plot: ", num_bars, "\n");
		cat("Num categories to plot: ", num_cat, "\n");

		def_par=par(no.readonly=T);

		margin=def_par$mar;
		margin[1]=12;
		margin[3]=1;
		par(mar=margin);

		# Compute placement of bars
		bar_offset=numeric();
		bar_offset[1]=0;	
		for(i in 2:num_bars){
			bar_offset[i]=bar_offset[i-1]+bar_width+bar_spacing;
		}	
		bar_x_end=bar_offset[num_bars]+bar_width+bar_spacing;

		# Label categories on the average bar plot
		label_y_pos=top_matrix[,pos_of_avg];
		num_top_cat=nrow(top_matrix);
		component_top=0;
		for(i in 1:num_top_cat){
			label_y_pos[i]=top_matrix[num_top_cat-i+1, pos_of_avg]/2+component_top;
			component_top=component_top+top_matrix[num_top_cat-i+1, pos_of_avg];
		}
		label_y_pos=rev(label_y_pos);

		#--------------------------------------------------------------
		# Adjust position of labels so there is good spacing

		adj_pos=rev(label_y_pos);
		min_spacing=.04;
		num_values=length(label_y_pos);

	        adj=min_spacing/100;
		cat("Adjustment Increment: ", adj, "\n");

		for(adj_ix in 1:100000){
			adjusted=F
			for(forw in 1:(num_values-1)){
				if(abs(adj_pos[forw]-adj_pos[forw+1])<min_spacing){
					adj_pos[forw]=adj_pos[forw]-adj;
					adjusted=T;
					break;
				}
			}
			for(rev in (num_values:2)){
				if(abs(adj_pos[rev]-adj_pos[rev-1])<min_spacing){
					adj_pos[rev]=adj_pos[rev]+adj;
					adjusted=T;
					break;
				}
			}
			if(!adjusted){
				cat("Done adjusting at iteration: ", adj, "\n");
				break;
			}
		}
		adj_pos=rev(adj_pos);
		label_range=range(c(adj_pos, 0, 1));

		#--------------------------------------------------------------

		# Initialize plot
		plot(0, type="n", xlim=c(0, bar_x_end+category_label_spacing), 
			ylim=label_range, xaxt="n", xlab="", ylab="Proportion", main=title);

		# Draw individual bar plot bars
		for(bar_ix in 1:num_bars){
			y_offset=0;
			for(cat_ix in 1:num_cat){
				rect_bottom=y_offset;
				rect_top=rect_bottom+tm[num_cat-cat_ix+1, bar_ix];

				rect(bar_offset[bar_ix], rect_bottom,
					bar_offset[bar_ix]+bar_width, rect_top,
					col=colors[num_cat-cat_ix+1], 
					border="black",
					lwd=.1
					#border=colors[num_cat-cat_ix+1]
				);
				
				y_offset=rect_top;	
			}
		}

		# Draw outline
		for(bar_ix in 1:num_bars){
			border_thck=1.25;
			if(bar_ix==pos_of_avg){ border_thck=2.5; }
			rect(bar_offset[bar_ix], 0, bar_offset[bar_ix]+bar_width, 1, lwd=border_thck);
		}

		# Label factor levels
		font_arr=rep(1, num_bars);
		font_arr[pos_of_avg]=2;
		label_pos=bar_offset+bar_width/2;
		label_val=paste(factor_levels, " (", num_samples[factor_levels], ")", sep="");
		bars_ix=1:num_bars;
		# Label regular labels
		axis(side=1, at=label_pos[bars_ix[-pos_of_avg]], labels=label_val[bars_ix[-pos_of_avg]], las=2);
		# Label average in bold
		axis(side=1, at=label_pos[pos_of_avg], labels=label_val[pos_of_avg], las=2, font.axis=2, cex=1.1);

		# Draw tick marks from center of bar
		for(i in 1:num_top_cat){
			points(c(bar_x_end-bar_spacing, bar_x_end-bar_spacing/2), c(label_y_pos[i], label_y_pos[i]), type="l");
		}

		# Draw line from center of bar to position of adjusted labels
		for(i in 1:num_top_cat){
			points(c(bar_x_end-bar_spacing/2, bar_x_end), c(label_y_pos[i], adj_pos[i]), type="l");
		}
		for(i in 1:num_top_cat){
			points(c(bar_x_end, bar_x_end+bar_spacing/3), c(adj_pos[i], adj_pos[i]), type="l");
		}

		# Label the categories at their adjusted positions
		text(bar_x_end, adj_pos, top_categories, pos=4, cex=.9);
			
			
		par(def_par);

	}

	plot_sbs(top_matrix, num_samples, title=title);
	

	#quit();

}

##############################################################################

# Load matrix
counts=load_summary_file(SummaryFile);
num_taxa=ncol(counts);
num_samples=nrow(counts);
if(TrimLongNames){
	colnames(counts)=trim_long_names(colnames(counts));
}
category_names=colnames(counts);

# Normalize
normalized=normalize(counts);
#print(normalized);

##############################################################################

# Load factors
factors=load_factors(FactorsFile);
factor_names=colnames(factors);
num_factors=ncol(factors);
factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);

cat("\n");
cat(num_factors, " Factor(s) Loaded:\n", sep="");
cat("\n");

##############################################################################

# Reconcile factors with samples
factor_sample_ids=rownames(factors);
counts_sample_ids=rownames(counts);

#print(factor_sample_id);
#print(counts_sample_id);

shared_sample_ids=intersect(factor_sample_ids, counts_sample_ids);
num_shared_sample_ids=length(shared_sample_ids);
num_factor_sample_ids=length(factor_sample_ids);
num_counts_sample_ids=length(counts_sample_ids);

cat("Num counts sample IDs: ", num_counts_sample_ids, "\n");
cat("Num factor sample IDs: ", num_factor_sample_ids, "\n");
cat("Num shared sample IDs: ", num_shared_sample_ids, "\n");
cat("\n");

cat("Samples missing from count information:\n");
print(setdiff(factor_sample_ids, counts_sample_ids));
cat("\n");
cat("Samples missing from factor information:\n");
print(setdiff(counts_sample_ids, factor_sample_ids));
cat("\n");
cat("Total samples shared: ", num_shared_sample_ids, "\n");

shared_sample_ids=sort(shared_sample_ids);

# Reorder data by sample id
normalized=normalized[shared_sample_ids,];
num_samples=nrow(normalized);
factors=factors[shared_sample_ids,, drop=F];

##############################################################################

#print(normalized);

#print(factors);


pdf(paste(OutputRoot, ".sbs_plot.pdf", sep=""), width=11, height=8.5);

factor_names=colnames(factors);
factor_samples=rownames(factors);

for(cur_factor in factor_names){
	cat("Working on: ", cur_factor, "\n");

	# Identify levels for current factor
	found_levels=as.character(sort(unique(factors[,cur_factor])));
	found_levels=setdiff(found_levels, NA);
	num_found_levels=length(found_levels);
	cat("Found Levels: \n");
	print(found_levels);
	cat("Number of levels found: ", num_found_levels, "\n");

	# Allocate space to store mean, CI and N
	mean_matrix=matrix(0, nrow=num_taxa, ncol=num_found_levels);
	CI95_matrix=matrix(0, nrow=num_taxa, ncol=num_found_levels);
	num_samples_arry=numeric(num_found_levels);

	# Name the cells so we won't used indices
	rownames(mean_matrix)=colnames(normalized);
	rownames(CI95_matrix)=colnames(normalized);
	colnames(mean_matrix)=found_levels;
	colnames(CI95_matrix)=found_levels;
	names(num_samples_arry)=found_levels;

	for(cur_level in found_levels){
		cat("Cur level: ", cur_level, "\n");

		# Pull out samples for each level that aren't NA
		cur_lvl_samp=factor_samples[factors[,cur_factor]==cur_level];
		cur_lvl_samp=cur_lvl_samp[!is.na(cur_lvl_samp)];
		cur_level_norm_abund=normalized[cur_lvl_samp,, drop=F];

		# Compute mean for each level
		categ_mean_acr_samples=apply(cur_level_norm_abund, 2, mean);
		mean_matrix[, cur_level]=categ_mean_acr_samples;
		num_samples_arry[cur_level]=nrow(cur_level_norm_abund);
	}

	plot_compositions(mean_matrix, num_samples_arry, title=cur_factor);

}

##############################################################################
	
dev.off();

##############################################################################

cat("Done.\n");
print(warnings());
q(status=0);
