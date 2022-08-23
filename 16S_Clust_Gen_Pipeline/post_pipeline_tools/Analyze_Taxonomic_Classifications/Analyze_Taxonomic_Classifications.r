#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library('getopt');
library(stringr);

options(useFancyQuotes=F);
options(digits=5)

params=c(
	"taxonomy", "t", 1, "character",
	"counts", "c", 1, "character",
	"outputroot", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-t <taxonomy classification file>\n",
	"	-c <counts table (to adjust stats for representatives)>\n",
	"	-o <output root filename>\n",
	"\n",
	"This script will generate statistics on the RDP classification\n",
	"confidences.\n",
	"\n");

if(
	!length(opt$taxonomy) || 
#	!length(opt$counts) ||
	!length(opt$outputroot)
){
	cat(usage);
	q(status=-1);
}

TaxonomyFile=opt$taxonomy;
CountsFile=opt$counts;
OutputFnameRoot=opt$outputroot;

cat("\n");
cat("Taxonomy Filename: ", TaxonomyFile , "\n", sep="");
cat("Counts Filename: ", CountsFile, "\n", sep="");
cat("Output Filename Root: ", OutputFnameRoot, "\n", sep="");
cat("\n");

##############################################################################

load_taxonomy_file=function(fname){
	# Read ID, Classification
	data=as.data.frame(read.delim(fname,  header=TRUE, row.names=NULL, check.names=FALSE, sep="\t"));
	return(data);
}

##############################################################################

load_count_table=function(fname){

	# Read in first few lines, then table
	fh=file(fname, "r");
	comments=readLines(fh,n=1);			# Skip comments
	sampleoffsets_example=readLines(fh,n=1);	# Skip examples

	# Read in order/offset of each sample id being represented
	sampleoffset_line=readLines(fh,n=1);
	sampleoffsets=strsplit(sampleoffset_line, "\\t")[[1]];
	sample_ids=sampleoffsets[3:length(sampleoffsets)];	# First 2 fields are headers

	# Accumulate to return
	rep_to_sample_map=list();
	
	cat("Loading Rep to Sample Mapping...\n");
	line=readLines(fh,n=1);
	while(length(line)>0){
	
		# Split the lines into representative/num_reps/represented,counts,etc.
		lineinfo=strsplit(line, "\\t")[[1]];

		representative=lineinfo[1];
		num_represented=lineinfo[2];
		source_ids=lineinfo[3:length(lineinfo)];

		#cat("\nRep: ", representative, "\n");
		#print(source_ids);
		
		translated=c();
		for(i in 1:length(source_ids)){
			src_info=strsplit(source_ids[i], ",")[[1]];
			repres=as.numeric(src_info[1]);
			counts=src_info[2];
			translated=c(translated, paste(sample_ids[repres], ",", counts, sep=""));
		}

		rep_to_sample_map[[representative]]=translated;
		line=readLines(fh,n=1);
	}

	close(fh);

	# Flatten to get sample-to-represented map
	repres=names(rep_to_sample_map);
	sample_id_to_repmap=list();
	for(rep_ix in repres){
		reps=rep_to_sample_map[[rep_ix]];

		for(sample_info in reps){
			sample_info=strsplit(sample_info, ",")[[1]];
			samp_id=sample_info[1];
			rep_num=sample_info[2];

			sample_id_to_repmap[[samp_id]]=
				c(sample_id_to_repmap[[samp_id]], rep(rep_ix, rep_num));
		}

	}

	# Sort the list items
	sample_id_to_repmap=sample_id_to_repmap[order(names(sample_id_to_repmap))];

	return(sample_id_to_repmap);
}

##############################################################################

plot_text=function(strings, max_lines=75){

	plot_page=function(strings){
		orig_par=par(no.readonly=T);

		par(mfrow=c(1,1));
		par(family="Courier");
		par(oma=rep(.5,4));
		par(mar=rep(0,4));

		num_lines=length(strings);

		top=max_lines;

		plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
			xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
			);
		for(i in 1:num_lines){
			#cat(strings[i], "\n", sep="");
			text(0, top-i, strings[i], pos=4, cex=.8);
		}

		par(orig_par);
	}

	num_lines=length(strings);
	num_pages=ceiling(num_lines / max_lines);
	#cat("Num Pages: ", num_pages, "\n");
	for(page_ix in 1:num_pages){
		start=(page_ix-1)*max_lines+1;
		end=start+max_lines-1;
		end=min(end, num_lines);
		##print(c(start,end));
		plot_page(strings[start:end]);
	}
}

##############################################################################

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, counts=F){

	orig_par=par(no.readonly=T);
	par(mfrow=c(1,1));

        num_row=nrow(mat);
        num_col=ncol(mat);

        cat("Num Rows: ", num_row, "\n");
        cat("Num Cols: ", num_col, "\n");

        mat=mat[rev(1:num_row),, drop=F];

        num_colors=50;
        color_arr=rainbow(num_colors, start=0, end=4/6);
        if(high_is_hot){
                color_arr=rev(color_arr);
        }

        remap=function(in_val, in_range, out_range){
                in_prop=(in_val-in_range[1])/(in_range[2]-in_range[1])
                out_val=in_prop*(out_range[2]-out_range[1])+out_range[1];
                return(out_val);
        }

        if(is.na(plot_min)){
                plot_min=min(mat);
        }
        if(is.na(plot_max)){
                plot_max=max(mat);
        }
        cat("Plot min/max: ", plot_min, "/", plot_max, "\n");
	par(mar=c(10,10,1,1));
        plot(0, type="n", xlim=c(0,num_col), ylim=c(0,num_row), 
		xaxt="n", yaxt="n", bty="n", xlab="", ylab="", main=title);

        # x-axis
        axis(side=1, at=seq(.5, num_col-.5, 1), labels=colnames(mat), las=2, cex.axis=.7);
        axis(side=2, at=seq(.5, num_row-.5, 1), labels=rownames(mat), las=2, cex.axis=.7);

        if(log_col){
                plot_min=log10(plot_min+.0125);
                plot_max=log10(plot_max+.0125);
        }

	text_size=min(1, 17.5/num_row);
	cat("Text Size: ", text_size, "\n");

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

                        if(counts){
                                text_lab=sprintf("%i", mat[y,x]);
                        }else{
                                text_lab=mat[y,x];
                        }
                        text(x-.5, y-.5, text_lab, srt=45, cex=text_size, font=2);
                }
        }

	par(orig_par);

}

##############################################################################

##############################################################################
# Main Program Starts Here!
##############################################################################

pdf(paste(OutputFnameRoot, ".taxa_class_analysis.pdf", sep=""), height=11, width=8.5);

param_msg=capture.output({
	cat("Taxonomy Filename:\n");
	cat("  ", TaxonomyFile, "\n");
	cat("\n");
	cat("Counts Filename:\n");
	cat("  ", CountsFile, "\n");
	cat("\n");
	cat("Output Filename Root: ", OutputFnameRoot, "\n");
});

# Load Counts
cat("Loading Counts File...\n");
sample_to_rep_map=load_count_table(CountsFile);

# Load Taxonomy classifications
cat("Loading Taxonomy File...\n");
taxa_info=load_taxonomy_file(TaxonomyFile);
dimen=dim(taxa_info);
cat("Row: ", dimen[1], "\n");
cat("Col: ", dimen[2], "\n");

repr_ids=as.character(taxa_info[,1]);
repr_class=as.character(taxa_info[,2]);

##############################################################################

num_reps=nrow(taxa_info);
num_sample_ids=length(sample_to_rep_map);
sample_ids=names(sample_to_rep_map);

cat("Number of Unique Sample IDs: ", num_sample_ids, "\n");

plot_text(c(
	param_msg,
	"",
	"Taxonomy Classifications:",
	paste("Num Representatives: ", dimen[1], sep=""),
	"",
	"Count Table:",
	paste("Num Samples: ", num_sample_ids, sep="")
)); 

##############################################################################

# Parse out genus and confidence for taxa class
taxa_splits=strsplit(repr_class, ";");
num_repre=nrow(taxa_info);

rep_table=as.data.frame(matrix(0, nrow=num_repre, ncol=2));
rownames(rep_table)=repr_ids;
colnames(rep_table)=c("genera", "confidence");

for(i in 1:num_repre){
	tok=tail(taxa_splits[[i]], 1);
	captures=str_match(tok, "^(.+)\\((\\d+)\\)$");
	match=captures[,1];
	rep_table[i,"genera"]=captures[,2];
	rep_table[i,"confidence"]=as.numeric(captures[,3]);
}

# Combine count table with classifications through representatives
taxa_conf_list_by_sample=list();
taxa_conf_list_overall=list();
j=1;
for(sid in sample_ids){

	cat("\n\nWorking on expanding sample: ", sid, "(", j, "/", num_sample_ids, ")\n");

	reps=sample_to_rep_map[[sid]];
	expanded_reads_table=rep_table[reps,];

	num_reads_in_sample=nrow(expanded_reads_table);
	taxa_hash=list();
	for(i in 1:num_reads_in_sample){

		cur_taxa=expanded_reads_table[i,"genera"];
		cur_conf=expanded_reads_table[i,"confidence"];

		# Keep sample specific table
		taxa_hash[[cur_taxa]]=
			c(
				taxa_hash[[cur_taxa]],
				cur_conf
			);

		# Add to overall table
		taxa_conf_list_overall[[cur_taxa]]=
			c(
				taxa_conf_list_overall[[cur_taxa]],
				cur_conf
			);

	}

	taxa_conf_list_by_sample[[sid]]=taxa_hash;

	j=j+1;
}

#print(taxa_conf_list_by_sample);
#print(taxa_conf_list_overall);

conf_list_to_stats=function(conf_list, subsamp_size=300){
	num_taxa=length(conf_list);
	taxa_names=names(conf_list);

	median_conf=numeric(num_taxa);
	lb95_conf=numeric(num_taxa);
	ub95_conf=numeric(num_taxa);
	num_reads=numeric(num_taxa);
	total_reads=0;

	names(median_conf)=taxa_names;
	names(lb95_conf)=taxa_names;
	names(ub95_conf)=taxa_names;
	names(num_reads)=taxa_names;

	for(i in 1:num_taxa){
		qnts=quantile(conf_list[[i]], c(0.025, 0.5, 0.975), na.rm=T);
		lb95_conf[i]=qnts[1];
		median_conf[i]=qnts[2];
		ub95_conf[i]=qnts[3];
		num_reads[i]=length(conf_list[[i]]);
		total_reads=total_reads+num_reads[i];
	}

	sort_ix=order(num_reads, decreasing=T);
	median_conf=median_conf[sort_ix];
	num_reads=num_reads[sort_ix];
	normalized=num_reads/sum(num_reads);
	num_top_taxa=sum(cumsum(normalized)<.95);
	top_taxa_names=names(normalized)[1:num_top_taxa];

	remaining=1-sum(normalized[1:num_top_taxa]);
	top_w_remaining=c(normalized[1:num_top_taxa], remaining);
	top_w_remaining_names=c(names(normalized[1:num_top_taxa]), "remaining");
	names(top_w_remaining)=top_w_remaining_names;

	# Calculate for remaining taxa
	remaining_taxa=setdiff(taxa_names, top_taxa_names);
	rem_conf=c();
	for(rt in remaining_taxa){
		rem_conf=c(rem_conf, conf_list[[rt]]);
	}
	qnts=quantile(rem_conf, c(0.025, 0.5, 0.975), na.rm=T);
	top_med_conf_w_rem=c(median_conf[top_taxa_names], qnts[2]);
	top_lb95_conf_w_rem=c(lb95_conf[top_taxa_names], qnts[1]);
	top_ub95_conf_w_rem=c(ub95_conf[top_taxa_names], qnts[3]);
	names(top_med_conf_w_rem)=top_w_remaining_names;
	names(top_lb95_conf_w_rem)=top_w_remaining_names;
	names(top_ub95_conf_w_rem)=top_w_remaining_names;

	# Subsample (for later plotting)
	subsamp=list();
	for(i in 1:num_taxa){
		cur_taxa_name=taxa_names[i];
		num_reads_for_taxa=length(conf_list[[i]]);
		target_sss=min(num_reads_for_taxa, subsamp_size);
                subsamp[[cur_taxa_name]]=sample(conf_list[[i]], target_sss);
        }
	# Subsample for remaining
	num_reads_for_taxa=length(rem_conf);
	target_sss=min(num_reads_for_taxa, subsamp_size);
	subsamp[["remaining"]]=sample(rem_conf, target_sss);
	
	# Set up return list/values
	results=list();
	# All
	results[["MedianConfidence"]]=median_conf;	
	results[["LB95Confidence"]]=lb95_conf;	
	results[["UB95Confidence"]]=ub95_conf;	
	results[["NumReads"]]=num_reads;
	results[["Normalized"]]=normalized;
	results[["Subsampled"]]=subsamp;
	results[["TotalReads"]]=total_reads;
	# Top + Remaining
	results[["NumTop"]]=num_top_taxa;
	results[["TopTaxaNames"]]=top_w_remaining_names;
	results[["TopWRemaining"]]=top_w_remaining;
	results[["TopConfwRem_Med"]]=top_med_conf_w_rem;
	results[["TopConfwRem_LB95"]]=top_lb95_conf_w_rem;
	results[["TopConfwRem_UB95"]]=top_ub95_conf_w_rem;

	return(results);
}

overall_stats=conf_list_to_stats(taxa_conf_list_overall);
#print(overall_stats);

##############################################################################

plot_stats=function(stats_record, title){

	cat("Generating plots for: ", title, "\n");
	layout_mat=matrix(c(
		1,1,
		2,2,
		2,2), byrow=T, nrow=3, ncol=2);
	layout(layout_mat);

	# Generate Rank Abundance Plot
	par(mar=c(2,5,2,1));
	cat("Plotting Abundances...\n");
	barplot(stats_record[["TopWRemaining"]], 
		ylim=c(0, max(stats_record[["TopWRemaining"]], na.rm=T)*1.1),
		#names.arg=c(names(stats_record[["TopWRemaining"]])), las=2,
		names.arg="",
		ylab="Abundance",
		main=paste("Top Abundances: ", title, sep=""));
	title(main=paste("Total Reads: ", stats_record[["TotalReads"]]), line=-.5, cex.main=1);

	# Plot confidences
	par(mar=c(15,5,2,1));
	cat("Plotting Confidences...\n");
	toptaxa=stats_record[["TopTaxaNames"]];
	print(toptaxa);
	print(stats_record[["TopConfwRem_UB95"]]);
	mids=barplot(stats_record[["TopConfwRem_Med"]][toptaxa],
		ylim=c(0, max(stats_record[["TopConfwRem_UB95"]][toptaxa], na.rm=T)*1.1),
		names.arg=toptaxa, las=2,
		ylab="Confidence",
		main=paste("Top Class. Confid.: ", title, sep=""));

	# Scatter draw samples
	i=1;
	for(tx in toptaxa){
		samp=stats_record[["Subsampled"]][[tx]];
		scat_x=rnorm(length(samp), 0, (mids[2]-mids[1])/10);
		scat_y=rnorm(length(samp), 0, .25);
		points(mids[i]+scat_x, samp+scat_y, cex=.7, lwd=.5);
		i=i+1;
	}

	# Mark intervals
	points(mids, stats_record[["TopConfwRem_LB95"]][toptaxa], col="Green", pch="-", cex=2);
	points(mids, stats_record[["TopConfwRem_UB95"]][toptaxa], col="Red", pch="-", cex=2);

}

plot_stats(overall_stats, "Overall Reads");

##############################################################################

samp_ids=names(taxa_conf_list_by_sample);
num_sample_ids=length(samp_ids);

# Accumulate medians into single matrix
cat("Accumulating Medians into single matrix...\n");
overall_top_taxa=setdiff(overall_stats[["TopTaxaNames"]], "remaining");
cat("Overall Top Taxa:\n");
print(overall_top_taxa);
num_overall_top_taxa=overall_stats[["NumTop"]];
sample_medians_mat=matrix(NA, nrow=num_sample_ids, ncol=num_overall_top_taxa);
rownames(sample_medians_mat)=samp_ids;
colnames(sample_medians_mat)=overall_top_taxa;

# Calculate stats on individual samples, plot them and accumulate medians
i=1;
for(smp_ix in samp_ids){
	cat("Calculating Stats for: ", smp_ix, "(", i, "/", num_sample_ids, ")\n");
	smp_stats=conf_list_to_stats(taxa_conf_list_by_sample[[smp_ix]]);
	plot_stats(smp_stats, smp_ix);

	for(tx_id in overall_top_taxa){
		sample_medians_mat[smp_ix, tx_id]=smp_stats[["MedianConfidence"]][tx_id];
	}
	i=i+1;
}

cat("Sampled Medians Matrix:\n");
print(sample_medians_mat);

##############################################################################

# Generate heatmap for across samples
cat("Generating heatmap for confidences across samples...\n");
max_samples_per_page=30;
num_pages=ceiling(num_sample_ids/max_samples_per_page);
for(i in 1:num_pages){
	start=(i-1)*max_samples_per_page;
	end=start+max_samples_per_page;

	end=min(end, num_sample_ids);

	paint_matrix(title=paste("Median Confidences by Sample: (", start+1, " - ", end, ")", sep=""),
		sample_medians_mat[(start+1):end,], plot_min=50, plot_max=100);
}

##############################################################################
# Generate barplot of medians with median sample confidences overlaid.
# summarize medians
cat("Generating barplots of medians of medians...\n");
median_of_medians=apply(sample_medians_mat, 2, function(x){median(x, na.rm=T);});

par(mfrow=c(1,1));
mids=barplot(median_of_medians, names.arg=colnames(sample_medians_mat), las=2,
	ylim=c(0,110), ylab="Median Confidence",
	main="Median Confidences Across Samples"
);
scatter=rnorm(num_sample_ids, 0, (mids[2]-mids[1])/10);
for(tix in 1:(num_overall_top_taxa)){
	points(rep(mids[tix], num_sample_ids)+scatter, sample_medians_mat[,tix], col="blue");
}

##############################################################################

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
