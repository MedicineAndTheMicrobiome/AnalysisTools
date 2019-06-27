#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
library(car);


source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');

options(useFancyQuotes=F);


params=c(
	"summary_file", "s", 1, "character",
	"summary_file2", "S", 2, "character",

	"factors", "f", 1, "character",
	"factor_samp_id_name", "F", 1, "character",
	"model_var", "M", 1, "character",
	"required", "q", 2, "character",

	"pairings", "p", 1, "character",
	"B_minuend", "B", 1, "character",
	"A_subtrahend", "A", 1, "character",

	"outputroot", "o", 2, "character",

	"reference_levels", "c", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

NUM_TOP_RESP_CAT=35;

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-s <summary file table>\n",
	"	[-S <second summary file table, in case pairings were in different files.>]\n",
	"\n",
	"	-f <factors file, contains covariates and factors>\n",
	"	-F <column name of sample ids in factor file>\n",
	"	-M <list of covariate X's names to include in the model from the factor file>\n",
	"	[-q <required list of variables to include after NA removal>]\n",
	"\n",
	"	-p <pairings map, pairing sample IDs from two groups. Must have header/column names>\n",
	"	-B <Sample Group B (minuend: B-A=diff), (column name in pairings file)\n",
	"	-A <Sample Group A (subtrahend: B-A=diff), (column name in pairings file)\n",
	"\n",
	"	[-o <output filename root>]\n",
	"\n",
	"	[-c <reference levels file for Y's in factor file>]\n",
	"\n",
	"This script will fit the following model for diversity:\n",
	"\n",
	"	1.) (Diversity_B_minuend[i] - Diversity_A_subtrahend[i]) = covariates\n",
	"\n",
	"The intercept of the model represents the difference between B and A.\n",
	"\n", sep="");

if(
	!length(opt$summary_file) || 
	!length(opt$factors) || 
	!length(opt$model_var) || 
	!length(opt$A_subtrahend) || 
	!length(opt$B_minuend) || 
	!length(opt$pairings) ||
	!length(opt$factor_samp_id_name)
){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputRoot=gsub(".summary_table.xls", "", opt$summary_file);
	OutputRoot=gsub(".summary_table.tsv", "", OutputRoot);
}else{
	OutputRoot=opt$outputroot;
}

if(!length(opt$summary_file2)){
	SecondSummaryTable="";
}else{
	SecondSummaryTable=opt$summary_file2;
}

if(!length(opt$reference_levels)){
        ReferenceLevelsFile="";
}else{
        ReferenceLevelsFile=opt$reference_levels;
}

if(length(opt$required)){
	RequiredFile=opt$required;
}else{
	RequiredFile="";
}

SummaryFile=opt$summary_file;
FactorsFile=opt$factors;
ModelVarFile=opt$model_var;
PairingsFile=opt$pairings;
A_subtrahend=opt$A_subtrahend;
B_minuend=opt$B_minuend;

FactorSampleIDName=opt$factor_samp_id_name;

OutputRoot=paste(OutputRoot, ".a_", A_subtrahend, ".b_", B_minuend, sep="");

cat("\n");
cat("         Summary File: ", SummaryFile, "\n", sep="");
if(SecondSummaryTable!=""){
	cat("     2nd Summary File: ", SecondSummaryTable, "\n", sep="");
}
cat("         Factors File: ", FactorsFile, "\n", sep="");
cat("Factor Sample ID Name: ", FactorSampleIDName, "\n", sep="");
cat(" Model Variables File: ", ModelVarFile, "\n", sep="");
cat("        Pairings File: ", PairingsFile, "\n", sep="");
cat("            A Minuend: ", A_subtrahend, "\n", sep="");
cat("         B Subtrahend: ", B_minuend, "\n", sep="");
cat("          Output File: ", OutputRoot, "\n", sep="");
cat("\n");
cat("Reference Levels File: ", ReferenceLevelsFile, "\n", sep="");
cat("\n");

options(width=100);
cat("Text Line Width: ", options()$width, "\n", sep="");

##############################################################################
##############################################################################

load_factors=function(fname){
	factors=data.frame(read.table(fname,  sep="\t", header=TRUE, row.names=1, 
		check.names=FALSE, comment.char=""));
	factor_names=colnames(factors);

	ignore_idx=grep("^IGNORE\\.", factor_names);

	if(length(ignore_idx)!=0){
		return(factors[-ignore_idx]);
	}else{
		return(factors);
	}
}

load_summary_file=function(fname){
	inmat=as.matrix(read.table(fname, sep="\t", header=TRUE, check.names=FALSE, 
	comment.char="", quote="", row.names=1))

	counts_mat=inmat[,2:(ncol(inmat))];

	# Clean category names a little
	cat_names=colnames(counts_mat);
	cat_names=gsub("-", "_", cat_names);
	colnames(counts_mat)=cat_names;
	
	cat("Num Categories in Summary Table: ", ncol(counts_mat), "\n", sep="");
	return(counts_mat);
}

load_reference_levels_file=function(fname){
        inmat=as.matrix(read.table(fname, sep="\t", header=F, check.names=FALSE, comment.char="#", row.names=1))
        colnames(inmat)=c("ReferenceLevel");
        print(inmat);
        cat("\n");
        if(ncol(inmat)!=1){
                cat("Error reading in reference level file: ", fname, "\n");
                quit(status=-1);
        }
        return(inmat);
}

relevel_factors=function(factors, ref_lev_mat){
        num_factors_to_relevel=nrow(ref_lev_mat);
        relevel_names=rownames(ref_lev_mat);
	factor_names=colnames(factors);
        for(i in 1:num_factors_to_relevel){
		
		target_relev_name=relevel_names[i];
		if(any(target_relev_name==factor_names)){
			tmp=factors[,target_relev_name];
			#print(tmp);
			tmp=relevel(tmp, ref_lev_mat[i, 1]);
			#print(tmp);
			factors[,target_relev_name]=tmp;
		}else{
			cat("Note: ", target_relev_name, " not in model.  Ignoring reference releveling.\n\n", sep="");
		}
        }
        return(factors);
}

merge_summary_tables=function(st1, st2){
	st1_cat_names=colnames(st1);
	st2_cat_names=colnames(st2);
	st1_samp=rownames(st1);
	st2_samp=rownames(st2);

	samp_names=sort(c(st1_samp, st2_samp));
	num_samp=length(samp_names);

	cat_names=sort(unique(c(st1_cat_names, st2_cat_names)));
	num_cat=length(cat_names);
	
	# Allocate
	merged_st=matrix(0, nrow=num_samp, ncol=num_cat);
	colnames(merged_st)=cat_names;
	rownames(merged_st)=samp_names;

	# Copy over
	merged_st[st1_samp, st1_cat_names]=st1[st1_samp, st1_cat_names];
	merged_st[st2_samp, st2_cat_names]=st2[st2_samp, st2_cat_names];
	return(merged_st);
}

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

load_list=function(filename){
	val=scan(filename, what=character(), comment.char="#");
	return(val);
}

load_mapping=function(filename, src, dst){
	mapping=read.table(filename, sep="\t", header=T, comment.char="#", quote="", row.names=NULL);
	column_names=colnames(mapping);
	if(all(column_names!=src)){
		cat("Error: Could not find ", src, " in header of map file.\n");
		quit(status=-1);
	}
	if(all(column_names!=dst)){
		cat("Error: Could not find ", dst, " in header of map file.\n");
		quit(status=-1);
	}

	map=cbind(as.character(mapping[,src]), as.character(mapping[,dst]));
	colnames(map)=c(src, dst);

	# Remove pairings with NAs
	incomp=apply(map, 1, function(x){any(is.na(x))});
	map=map[!incomp,];

	return(map);
}

intersect_pairings_map=function(pairs_map, keepers){

	missing=character();
	# Sets mappings to NA if they don't exist in the keepers array
	num_rows=nrow(pairs_map);
	for(cix in 1:2){
		for(rix in 1:num_rows){
			if(!any(pairs_map[rix, cix]==keepers)){
				missing=c(missing, pairs_map[rix, cix]);
				pairs_map[rix, cix]=NA;
			}
		}
	}
	results=list();
	results[["pairs"]]=pairs_map;
	results[["missing"]]=missing;
	return(results);
}

split_goodbad_pairings_map=function(pairs_map){
	
	num_rows=nrow(pairs_map);
	keepers=apply(pairs_map, 1, function(x){ all(!is.na(x))});

	good_pairs_map=pairs_map[keepers,,drop=F];
	bad_pairs_map=pairs_map[!keepers,,drop=F];
	num_good_collapsed_rows=nrow(good_pairs_map);
	cat("Collapsed ", num_rows, " pairs to ", num_good_collapsed_rows, " complete pairs.\n");

	res=list();
	res[["good_pairs"]]=good_pairs_map;
	res[["bad_pairs"]]=bad_pairs_map;
	return(res);
	
}

mask_matrix=function(val_mat, mask_mat, mask_thres, mask_val){
	masked_matrix=val_mat;
	masked_matrix[mask_mat>mask_thres]=mask_val;
	return(masked_matrix);
}

plot_text=function(strings){
	par(family="Courier");
	par(oma=rep(.1,4));
	par(mar=rep(0,4));

	num_lines=length(strings);
	
	top=max(as.integer(num_lines), 52);

	plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
		xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
		);

	text_size=max(.01, min(.8, .8 - .003*(num_lines-52)));
	#print(text_size);

	for(i in 1:num_lines){
		#cat(strings[i], "\n", sep="");
		strings[i]=gsub("\t", "", strings[i]);
		text(0, top-i, strings[i], pos=4, cex=text_size); 
	}
}

paint_matrix=function(mat, title="", plot_min=NA, plot_max=NA, log_col=F, high_is_hot=T, deci_pts=4, 
	label_zeros=T, counts=F, value.cex=1, 
	plot_col_dendr=F,
	plot_row_dendr=F
){

	cat("Working on: ", title, "\n");

        num_row=nrow(mat);
        num_col=ncol(mat);

	if(num_row==0 || num_col==0){
		cat("Nothing to plot.\n");
		return();
	}

	any_nas=any(is.na(mat));

	if(num_row==1 || any_nas){
		plot_row_dendr=F;
	}
	if(num_col==1 || any_nas){
		plot_col_dendr=F;
	}

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
		mat=mat[row_dendr[["names"]],, drop=F];
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
        axis(side=1, at=seq(.5, num_col-.5, 1), labels=colnames(mat), las=2, line=-1.75, cex.axis=value.cex);
        axis(side=4, at=seq(.5, num_row-.5, 1), labels=rownames(mat), las=2, line=-1.75, cex.axis=value.cex);

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

plot_histograms=function(table){
	num_cols=ncol(table);	
	orig.par=par(no.readonly=T);

	par(mfrow=c(5,3));
	par(mar=c(2,2,2,2));
	par(oma=c(2,2,2,2));
	colname=colnames(table);
	for(i in 1:num_cols){
		vals=table[,i];
		if(mode(vals)!="numeric" || is.factor(vals)){
			vals=as.factor(vals);
			barplot(prop.table(table(vals)), main=colname[i], col="white");
		}else{
			hist(vals, main=colname[i], xlab="values", ylab="");
		}
	}

	par(orig.par);
}

add_sign_col=function(coeff){
	cnames=colnames(coeff);
	pval_ix=which(cnames=="Pr(>|t|)");
	pval=coeff[,pval_ix];

	sig_char=function(val){
		if(val == 0){ return("***");}
		if(val <= .001){ return("** ");}
		if(val <= .01){ return("*  ");}
		if(val <= .05){ return(".  ");}
		return(" ");
	}

	sig_arr=sapply(pval, sig_char);
	fdr=round(p.adjust(pval, method="fdr"), 4);
	fdr_sig_char=sapply(fdr, sig_char);
	#out_mat=cbind(coeff, fdr);
	out_mat=cbind(coeff, sig_arr, fdr, fdr_sig_char);
	colnames(out_mat)=c(cnames, "Signf", "FDR", "Signf");
	
	return(formatC(out_mat, format="f", digits=5,));
	
}

mask_matrix=function(val_mat, mask_mat, mask_thres, mask_val){
        masked_matrix=val_mat;
        masked_matrix[mask_mat>mask_thres]=mask_val;
        return(masked_matrix);
}


tail_statistic=function(x){
        sorted=sort(x, decreasing=TRUE);
        norm=sorted/sum(x);
        n=length(norm);
        tail=0;
        for(i in 1:n){
                tail=tail + norm[i]*((i-1)^2);
        }
        return(sqrt(tail));
}

##############################################################################
##############################################################################

# Open main output file
pdf(paste(OutputRoot, ".paird_diff_div.pdf", sep=""), height=11, width=9.5);

# Load summary file table counts 
cat("\n");
cat("Loading summary table...\n");
counts1=load_summary_file(SummaryFile);

counts2=NULL;
if(SecondSummaryTable!=""){
	cat("Loading second summary table...\n");
	counts2=load_summary_file(SecondSummaryTable);	
	cat("Merging second summary table..\n");
	counts=merge_summary_tables(counts1, counts2);
	cat("Merged Summary Table Samples: ", nrow(counts), "\n", sep="");
	cat("Merged Summary Table Categories: ", ncol(counts), "\n", sep="");
}else{
	counts=counts1;
}

# Remove zero count samples
tot=apply(counts, 1, sum);
nonzero=tot>0;
if(!(all(nonzero))){
	cat("WARNING: Zero count samples found:\n");
	samp_names=rownames(counts);
	print(samp_names[!nonzero]);
	cat("\n");
	counts=counts[nonzero,,drop=F];
}

num_st_categories=ncol(counts);
num_st_samples=nrow(counts);

# Load resp/pred sample mappings
all_pairings_map=load_mapping(PairingsFile, A_subtrahend, B_minuend);
st_samples=rownames(counts);
num_pairings_loaded=nrow(all_pairings_map);
cat("\n");
cat("Pairing entries loaded: ", num_pairings_loaded, "\n");
print(all_pairings_map);

cat("Intersecting with samples in summary table:\n");
intersect_res=intersect_pairings_map(all_pairings_map, st_samples);
pairs=intersect_res[["pairs"]];
split_res=split_goodbad_pairings_map(pairs);
good_pairs_map=split_res$good_pairs;
bad_pairs_map=split_res$bad_pairs;

cat("Available pairs:\n");
print(good_pairs_map);
num_complete_pairings=nrow(good_pairs_map);
num_incomplete_pairings=nrow(bad_pairs_map);
cat("  Number of complete pairings: ", num_complete_pairings, "\n");
cat("Number of incomplete pairings: ", num_incomplete_pairings, "\n");

loaded_sample_info=c(
	"Summary Table Info:",
	paste(" 1st Summary Table Name: ", SummaryFile, sep=""),
	paste("    Samples: ", nrow(counts1), " Categories: ", ncol(counts1), sep=""),
	paste(" 2nd Summary Table Name: ", SecondSummaryTable, sep=""),
	paste("    Samples: ", nrow(counts2), " Categories: ", ncol(counts2), sep=""),
	"",
	paste("  Total Number Samples Loaded: ", num_st_samples, sep=""),
	paste("  Total Number Categories Loaded: ", num_st_categories, sep=""),
	"",
	"Sample Pairing Info:",
	paste("  Mapping Name: ", PairingsFile, sep=""),
	paste("  Number of Possible Pairings Loaded: ", num_pairings_loaded, sep=""),
	"",
	paste("Number of Complete/Matched Pairings: ", num_complete_pairings, sep=""),
	paste("Number of InComplete/UnMatched Pairings: ", num_incomplete_pairings, sep="")
);

if(length(intersect_res[["missing"]])){
	missing_info=c(
		"",
		"Missing:",
		capture.output(print(intersect_res[["missing"]]))
	);
}else{
	missing_info=c();
}

incomplete_pairing_info=c(
	"Incomplete Pairings: ",
	capture.output(print(bad_pairs_map)),
	"",
	paste("  Num complete pairings: ", num_complete_pairings, sep=""),
	paste("Num incomplete pairings: ", num_incomplete_pairings, sep=""),
	"",
	"(Double NA entries mean that the samples are missing from both groups.)",
	"(Known incomplete pairings (i.e. NAs in map file) are not included.)",
	missing_info
);

plot_text(loaded_sample_info);
plot_text(incomplete_pairing_info);

##############################################################################

# Normalize
cat("Normalizing counts...\n");
normalized=normalize(counts);

##############################################################################

# Load factors
cat("Loading Factors...\n");
factors=load_factors(FactorsFile);
factor_names=colnames(factors);
num_factors=ncol(factors);
factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);

num_loaded_factors=num_factors;
num_loaded_factor_samp=num_factor_samples;

cat("\n");
cat(num_factors, " Factor(s) Loaded:\n", sep="");
print(factor_names);
cat("\n");

# Load predictorss to include in model
if(ModelVarFile!=""){
	model_var_arr=load_list(ModelVarFile);
	cat("Model Variables:\n");
	print(model_var_arr);
	cat("\n");
}else if(ModelString!=""){
	cat("Error: Model String option features not implemented yet.\n");
	quit(status=-1);
}

# Load variables to require after NA removal
required_arr=NULL;
if(""!=RequiredFile){
	required_arr=load_list(RequiredFile);
	cat("Required Variables:\n");
	print(required_arr);
	cat("\n");
	missing_var=setdiff(required_arr, factor_names);
	if(length(missing_var)>0){
		cat("Error: Missing required variables from factor file:\n");
		print(missing_var);
	}
}else{
	cat("No Required Variables specified...\n");
}

plot_text(c(
	"Variables Targeted:",
	"",
	"Predictors:",
	capture.output(print(model_var_arr)),
	"",
	"Required Variables:",
	 capture.output(print(required_arr))
));

# Confirm we can find all the factors
missing_fact=setdiff(model_var_arr, factor_names);
if(length(missing_fact)>0){
	cat("Error: Factors in model, missing Factor File.\n");
	print(missing_fact);
	quit(status=-1);
}else{
	cat("All model variables found in factor file...\n");
}
factors=factors[,model_var_arr, drop=F];
factor_names=colnames(factors);
num_factors=ncol(factors);

# Relevel factor levels
if(ReferenceLevelsFile!=""){
        ref_lev_mat=load_reference_levels_file(ReferenceLevelsFile)
        factors=relevel_factors(kept_factors, ref_lev_mat);
}else{
        cat("No Reference Levels File specified.\n");
}

##############################################################################
# Reconcile factors with samples

cat("\nReconciling samples between factor file and paired mapping...\n");
factor_sample_ids=rownames(factors);
counts_sample_ids=rownames(counts);

shared_sample_ids=sort(intersect(factor_sample_ids, counts_sample_ids));

num_shared_sample_ids=length(shared_sample_ids);
num_factor_sample_ids=length(factor_sample_ids);
num_counts_sample_ids=length(counts_sample_ids);

cat("Num counts (paired responses) sample IDs: ", num_counts_sample_ids, "\n");
cat("Num factor sample IDs: ", num_factor_sample_ids, "\n");
cat("Num shared sample IDs: ", num_shared_sample_ids, "\n");
cat("\n");

cat("Samples missing from factor information:\n");
samples_missing_factor_info=(setdiff(counts_sample_ids, factor_sample_ids));
num_samp_missing_fact_info=length(samples_missing_factor_info);
cat("\n");
cat("Total samples shared: ", num_shared_sample_ids, "\n");

# Remove samples not in summary table 
cat("Adjusting pairings map based on factor/summary table reconciliation...\n");
intersect_res=intersect_pairings_map(good_pairs_map, shared_sample_ids);
pairs=intersect_res[["pairs"]];
split_res=split_goodbad_pairings_map(pairs);
good_pairs_map=split_res$good_pairs;
bad_pairs_map=split_res$bad_pairs;
paired_samples=as.vector(good_pairs_map);
num_shared_sample_ids=length(paired_samples);

# Reorder data by sample id
normalized=normalized[paired_samples,];
num_samples=nrow(normalized);
recon_factors=factors[paired_samples,,drop=F];

factor_file_info=c(
	paste("Factor File Name: ", FactorsFile, sep=""),
	"",
	paste("Num Loaded Factors/Variables: ", num_loaded_factors, sep=""),
	paste("Num Samples in Factor File: ", num_loaded_factor_samp, sep=""),
	"",
	paste("Num Samples Shared between Factors and Pairable Samples: ", num_shared_sample_ids, sep=""),
	"",
	paste("Num Samples Missing Factor Information: ", num_samp_missing_fact_info, sep=""),
	"",
	"Samples missing info: ",
	capture.output(print(samples_missing_factor_info))
);

plot_text(factor_file_info);

##############################################################################
# Remove samples with NAs

cat("Identifying samples/factors to keep with NAs...\n");
num_samples_recon=nrow(recon_factors);
num_factors_recon=ncol(recon_factors);
num_samples_before_na_removal=num_samples_recon;
num_factors_before_na_removal=num_factors_recon;

factors_wo_nas_res=remove_sample_or_factors_wNA_parallel(recon_factors, 
	required=required_arr, num_trials=64000, num_cores=64, outfile=paste(OutputRoot, ".noNAs", sep=""));

factors_wo_nas=factors_wo_nas_res$factors;
factor_names_wo_nas=colnames(factors_wo_nas);
factor_sample_ids_wo_nas=rownames(factors_wo_nas);
model_var_arr=intersect(model_var_arr, factor_names_wo_nas);

# Subset pairing map based on factor sample IDs
cat("Adjusting pairings map based on post-NA removal samples...\n");
intersect_res=intersect_pairings_map(good_pairs_map, factor_sample_ids_wo_nas);
pairs=intersect_res[["pairs"]]
split_res=split_goodbad_pairings_map(pairs);
good_pairs_map=split_res$good_pairs;
bad_pairs_map=split_res$bad_pairs;
paired_samples=as.vector(good_pairs_map);
num_shared_sample_ids=length(paired_samples);

# Subset the normalized counts based on pairing map
normalized=normalized[paired_samples,, drop=F];
num_samples_wo_nas=nrow(factors_wo_nas);
num_factors_wo_nas=ncol(factors_wo_nas);

#cat("Num Samples w/o NAs: ", num_samples_wo_nas, "\n");
#cat("Num Factors w/o NAs: ", num_factors_wo_nas, "\n");
#cat("\n");

##############################################################################

plot_text(c(
	paste("Num (Reconciled) Samples before NA removal: ", num_samples_before_na_removal, sep=""),
	paste("Num Factors before NA removal: ", num_factors_before_na_removal, sep=""),
	"",
	"Acceptable Variables after NA Removal:",
	"",
	capture.output(print(factor_names_wo_nas)),
	"",
	paste("Num Samples w/o NAs: ", num_samples_wo_nas, sep=""),
	paste("Num Factors w/o NAs: ", num_factors_wo_nas, sep="")
));

##############################################################################

# Order/Split the data....
# Order the pairings map by response sample IDs
A_sample_ids=good_pairs_map[,A_subtrahend];
B_sample_ids=good_pairs_map[,B_minuend];

factors=factors_wo_nas[A_sample_ids,,drop=F];

##############################################################################
# Compute diversity indices
cat("Computing diversity indices:\n");

div_names=c("Tail", "Shannon", "Simpson", "Evenness", "SimpsonsRecip");
num_div_idx=length(div_names);

div_mat=matrix(0, nrow=num_samples, ncol=num_div_idx);
colnames(div_mat)=div_names;
rownames(div_mat)=rownames(normalized);

cat("Computing diversity indices across samples.\n");
for(i in 1:num_samples){
        curNorm=normalized[i,];
        zeroFreeNorm=curNorm[curNorm>0]
        div_mat[i,"Tail"]=tail_statistic(zeroFreeNorm);
        div_mat[i,"Shannon"]=-sum(zeroFreeNorm*log(zeroFreeNorm));
        div_mat[i,"Simpson"]=1-sum(curNorm^2);
        div_mat[i,"Evenness"]=div_mat[i,"Shannon"]/log(length(zeroFreeNorm));
        div_mat[i,"SimpsonsRecip"]=1/sum(curNorm^2);
}

# Set evenness to 0 if there is only 1 category.
evenness=div_mat[,"Evenness"];
div_mat[is.na(evenness),"Evenness"]=0;

A_div_val=div_mat[A_sample_ids,,drop=F];
B_div_val=div_mat[B_sample_ids,,drop=F];

##############################################################################

num_model_pred=length(model_var_arr);

cov_model_string=paste(model_var_arr, collapse="+");
mmat=model.matrix(as.formula(paste("~", cov_model_string, "-1")), data=as.data.frame(factors));
cov_coeff_names=c("(Intercept)", colnames(mmat));
num_cov_coeff_names=length(cov_coeff_names);

cat("Anticipated Coefficient Names:\n");
print(cov_coeff_names);
cat("\n");

wilcoxon_pval_mat=matrix(NA, nrow=num_div_idx, ncol=1,
        dimnames=list(div_names, "Wilcoxon Paired"));

covariates_coef_mat  =matrix(NA, nrow=num_div_idx, ncol=num_cov_coeff_names, 
	dimnames=list(div_names, cov_coeff_names));
covariates_pval_mat  =matrix(NA, nrow=num_div_idx, ncol=num_cov_coeff_names, 
	dimnames=list(div_names, cov_coeff_names));

rsqrd_mat	=matrix(NA, nrow=num_div_idx, ncol=2, 
	dimnames=list(div_names[1:num_div_idx], c("R^2", "Adj-R^2")));

model_pval_mat	=matrix(NA, nrow=num_div_idx, ncol=1, 
	dimnames=list(div_names[1:num_div_idx], c("F-statistic P-value")));

###############################################################################

plot_ab_comparisons=function(a, b, aname, bname, pval, title){
# Plot histogram of differences and paired samples
	layout_mat=matrix(c(
		1,2,
		3,3,
		4,5
	), ncol=2, byrow=T);
	layout(layout_mat);
	par(oma=c(1,1,3,1));
	par(mar=c(4,4,3,1));

	ranges=range(c(a,b));
	ab=c(a,b);
	ab_hrec=hist(ab, plot=F);

	# Plot individual A/B
	ahrec=hist(a, breaks=ab_hrec$breaks, plot=F);
	bhrec=hist(b, breaks=ab_hrec$breaks, plot=F);
	max_count=max(c(ahrec$counts, bhrec$counts));
	ahrec=hist(a, breaks=ab_hrec$breaks, main=aname, ylim=c(0, max_count), xlab="ALR");
	bhrec=hist(b, breaks=ab_hrec$breaks, main=bname, ylim=c(0, max_count), xlab="ALR");

	# Plot difference
	dif=b-a;
	magn=max(abs(dif));
	diff_plot_range=c(-magn*1.2, magn*1.2);
	print(diff_plot_range);
	hist(dif, xlim=diff_plot_range,
		xlab="",
		main=paste("diff = ", bname, "-", aname, "\n", 
			"Wilcoxon Paired P-value:", round(pval, 4), sep=""), breaks=30);
	abline(v=0, col="blue", lty=2);

	# Plot scatter
	plot(a,b, xlim=ranges, ylim=ranges, xlab=aname, ylab=bname, cex=1.2);
	abline(a=0, b=1, col="blue", lty=2);

	# Plot barplot
	amean=mean(a);
	bmean=mean(b);

	bs_mean=function(x){
		if(length(x)<40){
			return(rep(mean(x),2));
		}
		num_bs=160;
		means=numeric(num_bs);
		for(i in 1:num_bs){
			means[i]=mean(sample(x, replace=T));
		}
		return(quantile(means, c(.025, .975)));
	}
	
	a95ci=bs_mean(a);
	b95ci=bs_mean(b);

	val_max=max(c(a, b));
	mids=barplot(c(amean, bmean), ylim=c(0, val_max*1.2), col="white");
	title(ylab="Bootstrapped 95% CI around Mean", line=2, cex.lab=.8, col.lab="blue");
	axis(1, at=mids, labels=c(aname, bname));
	barsep=diff(mids);
	
	draw95ci=function(ci95, x, width){
		if(ci95[1]==ci95[2]){return};
		points(c(x,x), ci95, type="l", col="blue", lwd=.7);
		points(c(x-width/2,x+width/2), rep(ci95[1],2), type="l", col="blue", lwd=.8);
		points(c(x-width/2,x+width/2), rep(ci95[2],2), type="l", col="blue", lwd=.8);
	}

	drawsamples=function(val, x, width){
		points(x+rnorm(length(val),0,width/5), val, cex=.5, col="grey");
	}

	drawsigbar=function(xs, height, pval){
		if(pval<=.1){
			if(pval<=.001){
				sigch="***";
			}else if(pval<=.01){
				sigch="**";
			}else if(pval<=.05){
				sigch="*";
			}else if(pval<=.1){
				sigch="";
			}
	
			points(xs, rep(height,2), type="l", col="black", lwd=.9);
			points(rep(xs[1],2), c(height,height*.97), type="l", col="black", lwd=.9);
			points(rep(xs[2],2), c(height,height*.97), type="l", col="black", lwd=.9);
			text(mean(xs), height, sigch, adj=c(.5,-1.2));
	
		}
	}

	drawsamples(a, mids[1], barsep/2);
	drawsamples(b, mids[2], barsep/2);

	draw95ci(a95ci, mids[1], barsep/2);
	draw95ci(b95ci, mids[2], barsep/2);

	#label mean and CI
	text(mids[1], a95ci[2], sprintf("mean = %5.3f", amean), col="blue", adj=c(.5,-2.4));
	text(mids[1], a95ci[2], sprintf("95%% CI = (%5.3f, %5.3f)", a95ci[1], a95ci[2]), 
		col="blue", adj=c(.5,-.8), cex=.8);

	text(mids[2], b95ci[2], sprintf("mean = %5.3f", bmean), col="blue", adj=c(.5,-2.4));
	text(mids[2], b95ci[2], sprintf("95%% CI = (%5.3f, %4.3f)", b95ci[1], b95ci[2]), 
		col="blue", adj=c(.5,-.8), cex=.8);

	drawsigbar(mids, val_max*1.1, pval);
	
	mtext(title, side=3, line=0, outer=T, font=2, cex=2);
}

###############################################################################

for(div_ix in 1:num_div_idx){

	A_div=A_div_val[,div_ix,drop=F];
	B_div=B_div_val[,div_ix,drop=F];

	cur_div_name=div_names[div_ix];
	
	cat("Fitting: ", div_ix, ".) ", cur_div_name, "\n");

	div_dif=(B_div-A_div);

	model_pred=cbind(factors, div_dif);

	model_str=paste("div_dif ~ ", paste(model_var_arr, collapse=" + "), sep="");

	lm_fit=lm(as.formula(model_str), data=model_pred);
	sum_fit=summary(lm_fit);

	print(sum_fit);
	ab_wilcox_res=wilcox.test(A_div, B_div, paired=T);
	wilcoxon_pval_mat[cur_div_name, 1]=ab_wilcox_res$p.value;

	plot_ab_comparisons(A_div, B_div, A_subtrahend, B_minuend, 
		pval=ab_wilcox_res$p.value,
		title=paste(div_ix, ".) ", cur_div_name, sep=""));


	par(mfrow=c(1,1));

	# ANOVA on full model
	anova_res=anova(lm_fit);

	plot_text(c(
		paste(div_ix, ".) ", cur_div_name, ":", sep=""),
		"",
		capture.output(print(sum_fit))
		)
	);

	plot_text(c(
		paste(div_ix, ".) ", cur_div_name, ":", sep=""),
                "",
		capture.output(print(anova_res))
		)
	);

	mmps(lm_fit);

	#model_coef_names=setdiff(rownames(sum_fit$coefficients), "(Intercept)");
	# The intercept will be zero if there is no difference, so we should report this.
	model_coef_names=rownames(sum_fit$coefficients);

	print(model_coef_names);

	# Save the covariate result stats
	coef_names=intersect(model_coef_names, cov_coeff_names);
	covariates_coef_mat[cur_div_name, coef_names]=sum_fit$coefficients[coef_names,"Estimate"];
	covariates_pval_mat[cur_div_name, coef_names]=sum_fit$coefficients[coef_names,"Pr(>|t|)"];

	rsqrd_mat[cur_div_name, "R^2"]=sum_fit$r.squared;
	rsqrd_mat[cur_div_name, "Adj-R^2"]=sum_fit$adj.r.squared;

	model_pval_mat[cur_div_name, "F-statistic P-value"]=
		1-pf(sum_fit$fstatistic[1], sum_fit$fstatistic[2], sum_fit$fstatistic[3]);

	cat("\n\n*************************************************\n\n");

}

###############################################################################

all.nas=apply(covariates_coef_mat, 2, function(x){all(is.na(x))});
covariates_coef_mat=covariates_coef_mat[,!all.nas,drop=F];

all.nas=apply(covariates_pval_mat, 2, function(x){all(is.na(x))});
covariates_pval_mat=covariates_pval_mat[,!all.nas,drop=F];

print(covariates_coef_mat);
print(covariates_pval_mat);
print(rsqrd_mat);
print(model_pval_mat);

par(oma=c(2,1,5,2));

rename_intercept=function(mat, old_n, new_n){
	cname=colnames(mat);
	old_ix=which(cname==old_n);
	cname[old_ix]=new_n;
	colnames(mat)=cname;
	return(mat);
}

covariates_coef_mat=rename_intercept(covariates_coef_mat, "(Intercept)", "\"Difference\"");
covariates_pval_mat=rename_intercept(covariates_pval_mat, "(Intercept)", "\"Difference\"");

print(covariates_coef_mat);

paint_matrix(wilcoxon_pval_mat, plot_min=0, plot_max=1,
	title="Wilcoxon Difference in ALR: P-values",
	high_is_hot=F, deci_pts=2, value.cex=1); 
mtext(text="(No controlling for covariates)", side=3, cex=.8, font=3, line=2, outer=T);

paint_matrix(covariates_coef_mat, 
        title="Regression Model Coefficient Values",
        high_is_hot=T, deci_pts=2, value.cex=1);

paint_matrix(covariates_pval_mat, plot_min=0, plot_max=1,
        title="Regression Model Coeff P-Values", 
        high_is_hot=F, deci_pts=2, value.cex=1);
mtext(text="(H0: Coefficients equal zero, H1: Non-zero Coefficients)", side=3, cex=.9, font=3, line=2, outer=T);

signf_coef_mat=mask_matrix(covariates_coef_mat, covariates_pval_mat, .1, 0);
paint_matrix(signf_coef_mat,
        title="Regression Model Significant Coefficients", 
        high_is_hot=T, deci_pts=2, value.cex=1, label_zeros=F);
mtext(text="(P-Values < 0.10 Shown)", side=3, cex=.9, font=3, line=2, outer=T);

paint_matrix(rsqrd_mat, plot_min=0, plot_max=1,
        title="Regression R^2's", 
        high_is_hot=T, deci_pts=2, value.cex=1);

paint_matrix(model_pval_mat, plot_min=0, plot_max=1,
        title="Regression Model Fit P-values", 
        high_is_hot=F, deci_pts=2, value.cex=1);
mtext(text="(H0: Predictors have no contribution to model fit)", side=3, cex=.8, font=3, line=2, outer=T);



###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
