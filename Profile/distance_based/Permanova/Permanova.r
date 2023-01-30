#!/usr/bin/env Rscript

###############################################################################

library(MASS);
library(vegan);
library('getopt');
options(useFancyQuotes=F);
options(width=120);

params=c(
	"distmat", "d", 1, "character",
	"factors", "f", 1, "character",
	"model_formula", "m", 2, "character",
	"model_variables_file", "M", 2, "character",
	"required_var", "q", 2, "character",
	"blocking", "b", 2, "character",
	"outputroot", "o", 2, "character",
	"xrange", "x", 2, "character",
	"yrange", "y", 2, "character",
	"testing", "T", 2, "logical",
	"strip_samples_nas", "s", 2, "logical",
	"tag_name", "t", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

script_path=paste(head(strsplit(script_name, "/")[[1]], -1), collapse="/");
source(paste(script_path, "/../../../Metadata/RemoveNAs/Remove_NAs.r", sep=""));

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-d <distance matrix>\n",
	"	-f <factors>\n",
	"\n",
	"	[-o <output filename root>]\n",
	"	[-m \"model formula string\"]\n",
	"	[-M <model variables filename>]\n",
	"\n",
	"	[-q <required variables list>]\n",
	"\n",
	"	[-b <factor to use as blocking variable>]\n",
	"	[--xrange=<MDS Dim1 Range, eg. -2,2>]\n",
	"	[--yrange=<MDS Dim2 Rnage, eg. -2,2>]\n",
	"	[-T (testing flag)]\n",
	"	[-t <tag name>]\n",
	"\n",
	"	[-s (Flag to strip samples with NAs, default=F)]\n",
	"\n",
	"This script will run Permutational Analysis of Variance (PERMANOVA)\n",
	"on your specified distance matrix, with the factors that are available.\n",
	"\n",
	"The distance matrix should have the sample IDs matching the first\n",
	"column of the factors. ",
	"\n",
	"The -x and -y command specify the range of how to generate the MDS Plot.\n",
	"\n",
	"Each pair of MDS plots (labeled and with centroids marked), are reoriented\n",
	"so that the centroid of first factor level samples is left of the second factor level.\n",
	"\n",
	"If you specify the -b option, it will assume that you want to perform a single\n",
	"level nesting, where you are essentially blocking by the groups specified in the\n",
	"variable in the blocking variable.\n",
	"\n",
	"If using the -t flag is for testing, so don't use it for production runs.\n",
	"\n",
	"If -s flag is set, samples are removed if any factors have NAs.\n",
	"By default, samples/factors are removed automatically to maximize the\n",	
	"the number of non-NA values.\n",
	"\n");

if(!length(opt$distmat) || !length(opt$factors)){
	cat(usage);
	q(status=-1);
}

if(!length(opt$outputroot)){
	OutputFnameRoot=gsub(".distmat", "", opt$distmat);
}else{
	OutputFnameRoot=opt$outputroot;
}
OutputFnameRoot=paste(OutputFnameRoot, ".perm", sep="");

if(!length(opt$model_formula)){
	ModelFormula="";
}else{
	ModelFormula=opt$model_formula;
}

if(length(opt$model_variables_file)){
        ModelVariablesFile=opt$model_variables_file;
}else{
        ModelVariablesFile="";
}

DistmatFname=opt$distmat;
FactorsFname=opt$factors;

Xrange=numeric(2);
if(!length(opt$xrange)){
	Xrange=c(-Inf, Inf);
}else{
	Xrange=as.numeric(strsplit(opt$xrange, ",")[[1]]);
}

Yrange=numeric(2);
if(!length(opt$yrange)){
	Yrange=c(-Inf, Inf);
}else{
	Yrange=as.numeric(strsplit(opt$yrange, ",")[[1]]);
}

Blocking="";
if(length(opt$blocking)){
	Blocking=opt$blocking;
	cat("Blocking variable: ", Blocking, "\n");
}

Testing=F;
if(length(opt$testing)){
	Testing=T;
	cat("**************************************************************\n");
	cat("*  Testing Flag Set...                                       *\n");
	cat("**************************************************************\n");
	rand=sprintf(".%04i", sample(1000,1));
}else{
	rand="";
}

StripSamplesWithNAs=F;
if(length(opt$strip_samples_nas)){
	StripSamplesWithNAs=T;
}

RequiredFile="";
if(length(opt$required_var)){
        RequiredFile=opt$required_var;
}

if(length(opt$tag_name)){
        TagName=opt$tag_name;
        cat("Setting TagName Hook: ", TagName, "\n");
        setHook("plot.new",
                function(){
                        #cat("Hook called.\n");
                        if(par()$page==T){
                                oma_orig=par()$oma;
                                exp_oma=oma_orig;
                                exp_oma[1]=max(exp_oma[1], 1);
                                par(oma=exp_oma);
                                mtext(paste("[", TagName, "]", sep=""), side=1, line=exp_oma[1]-1,
                                        outer=T, col="steelblue4", font=2, cex=.8, adj=.97);
                                par(oma=oma_orig);
                        }
                }, "append");

}else{
        TagName="";
}

###############################################################################

cat("\n");
cat("Distance Matrix Filename: ", DistmatFname, "\n", sep="");
cat("Factors Filename: ", FactorsFname, "\n", sep="");
cat("Output Filename Root: ", OutputFnameRoot, "\n", sep="");
cat("\n");
cat("X Range for MDS plot: ", Xrange[1], ", ", Xrange[2], "\n", sep="");
cat("Y Range for MDS plot: ", Yrange[1], ", ", Yrange[2], "\n", sep="");
cat("\n");

if(ModelFormula!=""){
	cat("Model Formula specified: ", ModelFormula, "\n\n");
}

cat("NAs in Metadata Policy:\n");
if(StripSamplesWithNAs){
	cat("Stripping out samples with NAs.\n");
}else{
	cat("Maximizing non-NAs.\n");
}

###############################################################################

load_distance_matrix=function(fname){
	distmat=as.matrix(read.delim(fname, sep=" ",  header=TRUE, row.names=1, 
		check.names=FALSE, comment.char="", quote=""));
	#print(distmat);
	mat_dim=dim(distmat);
	cat("Read in distance matrix: \n");
	cat("  Rows: ", mat_dim[1], "\n");
	cat("  Cols: ", mat_dim[2], "\n");

	# Remove NAs
	diag(distmat)=NA;
	non_na_rows=apply(distmat, 1, function(x){!all(is.na(x))});
	non_na_cols=apply(distmat, 2, function(x){!all(is.na(x))});
	diag(distmat)=0;
	distmat=distmat[non_na_rows, non_na_cols];

	if(mat_dim[1]!=mat_dim[2]){
		cat("Error: Distance Matrix is not squared.\n");
		print(colnames(distmat));
		print(rownames(distmat));
		q(status=-1);
	}
	return(distmat);
}

##############################################################################

load_factors=function(fname){

	cat("Loading factor file: ", fname, "\n");
	factors=data.frame(read.table(fname,  header=TRUE, row.names=1, 
		check.names=FALSE, sep="\t", quote="", comment.char=""));

	#print(factors);

	dimen=dim(factors);
	cat("Rows Loaded: ", dimen[1], "\n");
	cat("Cols Loaded: ", dimen[2], "\n");

	return(factors);
}

##############################################################################

orient_points_by_centroid=function(x, y, fact_col){
	num_points=length(fact_col);
	num_groups=length(unique(fact_col));
	x_centroids=numeric(num_groups);
	y_centroids=numeric(num_groups);

	# Compute centroids for all groups
	for(g in 1:num_groups){
		#cat("\n\nComputing centroid for group: ", g, "\n");
		members=which(fact_col==g)
		#print(members);
		x_centroids[g]=mean(x[members]);
		y_centroids[g]=mean(y[members]);
	}

	# Compute angle between first two factor levels
	arc=atan2(y_centroids[2]-y_centroids[1], x_centroids[2]-x_centroids[1]);
	
	rot=function(x, y, arc){
		rotated=list();
		rotated$x = x*cos(arc)-y*sin(arc);
		rotated$y = x*sin(arc)+y*cos(arc);
		return(rotated);
	}

	rotx=numeric(num_points);
	roty=numeric(num_points);
	rotcentx=numeric(num_groups);
	rotcenty=numeric(num_groups);

	for(i in 1:num_points){
		rot_res=rot(x[i], y[i], -arc);	
		rotx[i]=rot_res$x;
		roty[i]=rot_res$y;
	}
	for(i in 1:num_groups){
		rot_res=rot(x_centroids[i], y_centroids[i], -arc);	
		rotcentx[i]=rot_res$x;
		rotcenty[i]=rot_res$y;
	}

	result=list();	
	result$x=rotx;
	result$y=roty;
	result$x_centroids=rotcentx;
	result$y_centroids=rotcenty;
	return(result);
}

##############################################################################

compute_dispersion=function(variation, groups, group_names){

	cur_levels=sort(unique(as.character(groups)));
	num_levels=length(cur_levels);

	if(num_levels == 0){
		cur_levels=as.character(sort(unique(as.numeric(samp_levels[,1]))));
		num_levels=length(cur_levels);
		cat("Warning: Treating ordinal/continuous levels as categories for dispersion analysis.\n");
	}

	pval_matrix=matrix(0, nrow=num_levels, ncol=num_levels);		
	colnames(pval_matrix)=cur_levels;
	rownames(pval_matrix)=cur_levels;

	points=list();	

	for(li1 in 1:num_levels){

		l1_variation=variation[groups==cur_levels[li1]];
		points[[cur_levels[li1]]]=l1_variation;

		for(li2 in 1:num_levels){

			if(li1>li2){
				pval_matrix[li1, li2]=pval_matrix[li2, li1];	
			}else if(li1==li2){
				pval_matrix[li1, li2]=1;
			}else{

				#cat(cur_levels[li1], " vs. ", cur_levels[li2], "\n");
				l2_variation=variation[groups==cur_levels[li2]];

				#print(l1_variation);
				#print(l2_variation);
				result=wilcox.test(l1_variation, l2_variation);
				pval_matrix[li1, li2]=result$p.value;
			
			}
		}
	}

	#print(pval_matrix);
	colnames(pval_matrix)=group_names;
	rownames(pval_matrix)=group_names;
	names(points)=group_names;
	
	factor_dispersions=list();
	factor_dispersions[["pvals"]]=pval_matrix;
	factor_dispersions[["points"]]=points;

	return(factor_dispersions);
}

##############################################################################

plot_pval_heatmap=function(mat, title=""){

	orig_par=par(no.readonly=T);

        #par(family="Courier");
        par(mar=c(15.1, 14.1, 1.5, 1.5));

        # Generate colors from red to blue
        colors=(rainbow(2^16, start=0, end=0.65));

        # Remember that rows and columsn are reversed in the image
        image(1:nrow(mat),1:ncol(mat), mat,
                xaxt="n", yaxt="n",
                xlab="", ylab="",
                col=colors
        );

        # Pad strings
        cnames=paste(colnames(mat), " ", sep="");
        rnames=paste(rownames(mat), " ", sep="");

        # Get longest length of each column or row label
        cname_max_len=max(nchar(cnames));
        rname_max_len=max(nchar(rnames));

        # Get the number of rows and columns
        ncols=ncol(mat);
        nrows=nrow(mat);

	base_sf=10;
        cscale=min(c(base_sf/cname_max_len, base_sf/ncols));
        rscale=min(c(base_sf/rname_max_len, base_sf/nrows));

        max_width=max(nchar(sprintf("%.2f",mat)));
        cell_cex=(3.5/max_width)*sqrt(min(c(cscale, rscale))^2);

        for(i in 1:nrow(mat)){
                for(j in 1:ncol(mat)){
                        str=sprintf("%.2f",mat[i,j]);
                        str=gsub("0\\.",".", str);
                        text(i,j,labels=str, cex=cell_cex, srt=45);
                }
        }

        # Plot the labels
        mtext(cnames, at=1:ncols, side=2, las=2, cex=cscale);
        mtext(rnames, at=1:nrows, side=1, las=2, cex=rscale);

        # Plot the title
        mtext(title, line=0, at=nrows*.5, side=3, font=2);

	par(orig_par);
}

##############################################################################

plot_text=function(strings){
	orig_par=par(no.readonly=T);

        par(mfrow=c(1,1));
        par(family="Courier");
        par(oma=rep(.5,4));
        par(mar=rep(0,4));

        num_lines=length(strings);

        top=max(as.integer(num_lines), 40);

        plot(0,0, xlim=c(0,top), ylim=c(0,top), type="n",  xaxt="n", yaxt="n",
                xlab="", ylab="", bty="n", oma=c(1,1,1,1), mar=c(0,0,0,0)
                );
        for(i in 1:num_lines){
                #cat(strings[i], "\n", sep="");
                text(0, top-i, strings[i], pos=4, cex=.8);
        }

	par(orig_par);
}

##############################################################################

remove_samples_wNA=function(factors){
	
	cat("Identifying Samples to remove because factors have NAs.\n");
	isnas=is.na(factors);
		isnas=is.na(factors);
	samples_wNAs=apply(isnas, 1, any);
	return(factors[!samples_wNAs,,drop=F]);
}

##############################################################################

load_list=function(filename){
        val=scan(filename, what=character(), comment.char="#");
        return(val);
}

##############################################################################

sig_char=function(val){
        if(!is.null(val) && !is.nan(val) && !is.na(val)){
                if(val <= .0001){ return("***");}
                if(val <= .001 ){ return("** ");}
                if(val <= .01  ){ return("*  ");}
                if(val <= .05  ){ return(":  ");}
                if(val <= .1   ){ return(".  ");}
        }
        return(" ");
}

##############################################################################

subset_model_string=function(model_string, avail_factors){
	lin_var_str=gsub(" ", "", model_string);
	lin_comp=strsplit(lin_var_str, "\\+")[[1]];
	num_components=length(lin_comp);
	keep_comp=c();
	for(i in 1:num_components){
		vars=strsplit(lin_comp[i], "[\\*\\:]")[[1]];
		shared=intersect(vars, avail_factors);
		if(setequal(shared, vars)){
			keep_comp=c(keep_comp, lin_comp[i]);
		}
	}
	new_model=paste(keep_comp, collapse="+");
	return(new_model);
}

#subset_model_string("This+ is + a:test + yes + it:is+This:test:yes", c("This", "test", "yes"));

##############################################################################

pdf(paste(OutputFnameRoot, rand, ".pdf", sep=""), height=5.5, width=11);

##############################################################################

# Load distance matrix
distmat=load_distance_matrix(DistmatFname);
num_distmat_samples=ncol(distmat);
distmat_sample_names=colnames(distmat);
#print(distmat):

# Load factors
factors=load_factors(FactorsFname);
factor_sample_names=rownames(factors);
num_factor_samples=length(factor_sample_names);
factor_names=colnames(factors);
num_factors=ncol(factors);
cat(num_factors, " Factor(s) Loaded:\n", sep="");
print(factor_names);
num_factor_samples=length(factor_sample_names);
cat(num_factor_samples, " Samples in factor file.\n", sep="");
cat("\n");

###############################################################################

# Confirm/Reconcile that the samples in the matrix and factors file match
cat("DistMat Samples:\n");
print(distmat_sample_names);
cat("\n");
cat("Factor Samples:\n");
print(factor_sample_names);

common_sample_names=intersect(distmat_sample_names, factor_sample_names);
num_common_samples=length(common_sample_names);
if(num_common_samples < num_distmat_samples || num_common_samples < num_factor_samples){
	cat("\n");
	cat("*** Warning: The number of samples in factors file does not match those in your distance matrix. ***\n");
	cat("Taking intersection (common) sample IDs between both.\n");
	cat("Please confirm this is what you want.\n");
	cat("\tNum Distmat Samples: ", num_distmat_samples, "\n");
	cat("\tNum Factor  Samples: ", num_factor_samples, "\n");
	cat("\tNum Common  Samples: ", num_common_samples, "\n");	
	cat("\n");
}


# Set the working distance matrix to the same order
distmat=distmat[common_sample_names, common_sample_names];
factors=factors[common_sample_names, , drop=F];

###############################################################################

if(ModelVariablesFile!=""){
        model_variables_file_list=load_list(ModelVariablesFile);
        ModelFormula=paste(model_variables_file_list, collapse=" + ");
}

if(ModelFormula!=""){
	# Based on factors in model string, identity which factors are used
	model_vars_str=ModelFormula;
	model_vars_str=gsub(" ", "", model_vars_str);
	model_vars_str=gsub("[\\+\\:\\*]", " ", model_vars_str);
	model_var=unique(strsplit(model_vars_str, " ")[[1]]);

	avail_factors=colnames(factors);
	if(!setequal(model_var, intersect(model_var, avail_factors))){
		cat("ERROR: Could not find model variables in factor file.\n\n");
		cat("Missing Model Variables:\n");
		print(setdiff(model_var, avail_factors));
		cat("\nFactor File Variables:\n");
		print(avail_factors);
		cat("\n");
		quit(status=-1);
	}

	factors=factors[,model_var, drop=F];
	num_factors=ncol(factors);
}else{
	model_var=factor_names;
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

if(ncol(factors)==1){
	cat("Only one variable selected.\n");
	required_arr=c(colnames(factors));
	cat("Making ", required_arr, " required.\n");
}

# Decide what to do with NAs.
if(StripSamplesWithNAs){
	factors=remove_samples_wNA(factors);
}else{

	noNA_result=remove_sample_or_factors_wNA_parallel(factors, required=required_arr,
		num_trials=500000, num_cores=64, outfile=OutputFnameRoot);
	factors=noNA_result$factors;
	plot_text(noNA_result$summary_text);
}
factor_names=colnames(factors);
num_factors=ncol(factors);
factor_sample_names=rownames(factors);

if(ModelFormula!=""){
	ModelFormula=subset_model_string(ModelFormula, factor_names);
	cat("Adjusted Model Formula: ", ModelFormula, "\n");
}

# Reconcile samples between distance matrix and factor file again
common_sample_names=intersect(distmat_sample_names, factor_sample_names);
distmat=distmat[common_sample_names, common_sample_names];

num_samples=ncol(distmat);
sample_names=colnames(distmat);
cat("Num Samples used: ", num_samples, "\n\n");

for(i in 1:num_factors){
	categories=sort(unique(unique(factors[,i])));
	cat("'", factor_names[i], "' has ", length(categories), " categories.\n", sep="");
	cat("\t", paste(head(categories), collapse=", "), sep="");
	if(length(categories)>10){
		cat(" ...");
	}
	cat("\n");
}

##############################################################################

# Perform Permanova
dist=(as.dist(distmat));

if(ModelFormula==""){
	model_string=paste("dist ~", paste(factor_names, collapse=" + "));
}else{
	model_string=paste("dist ~", ModelFormula);
}

num_linear_components=length(strsplit(model_string, "\\+")[[1]]);

cat("\nFitting this model: ", model_string, "\n");

cat("\n--------------------------------------------------------------------------\n");
cat("Before invoking Adonis:\n");

if(Blocking!=""){
	stratify=factors[[Blocking]];
}else{
	stratify=NULL;
}

perm_factor=max(10, num_linear_components);


# Old version using type I sum of squares
old_res=adonis(as.formula(model_string), data=as.data.frame(factors), strata=stratify, permutations=perm_factor*1000);
cat("Old Adonis Results:\n");
print(names(old_res));
print(old_res);


res2=adonis2(as.formula(model_string), data=as.data.frame(factors), strata=stratify, permutations=perm_factor*1000, by="margin");
cat("New Adonis Results:\n");
print(names(res2));
print((res2));


# Merge new AOV table into old result to maintain other statistics we need
merge_new_with_old_adonis=function(old, new){

	numrows_old=nrow(old[["aov.tab"]]);
	numrows_new=nrow(new);

	cat("Number of Rows in old: ", numrows_old, " vs. new: ", numrows_new, "\n");

	old_row_names=rownames(old[["aov.tab"]]);
	new_row_names=rownames(new);
	new_row_names[(new_row_names=="Residual")]="Residuals";
	rownames(new)=new_row_names;

	print(old_row_names);
	print(new_row_names);	

	out_res=old_res;

	if(numrows_old==numrows_new && all(old_row_names==new_row_names)){
		cat("Moving AOV table from adonis2 into old adonis result.\n");

		aov.tab=old[["aov.tab"]];
		aov.tab[old_row_names, "Df"]=new[old_row_names, "Df"];
		aov.tab[old_row_names, "SumsOfSqs"]=new[old_row_names, "SumOfSqs"];
		aov.tab[old_row_names, "MeanSqs"]=new[old_row_names, "SumOfSqs"]/new[old_row_names, "Df"];
		aov.tab[old_row_names, "F.Model"]=new[old_row_names, "F"];
		aov.tab[old_row_names, "R2"]=new[old_row_names, "R2"];
		aov.tab[old_row_names, "Pr(>F)"]=new[old_row_names, "Pr(>F)"];
		attr(aov.tab, "heading")="";
		out_res[["aov.tab"]]=aov.tab;
		
		#print(out_res[["call"]]);
		out_res[["call"]]="Original Adonis results with Adonis2 margin (type II) Sum of Squares AOV";
		return(out_res);		
	}else{
		cat("Number of rows do not match between old and new Adonis results.\n");
		cat("Going forward with old Adonis results only...\n");
		return(NULL);
	}

}

res=merge_new_with_old_adonis(old_res, res2);

if(is.null(res)){

	adonis_out=capture.output(print(old_res));
	adonis2_out=capture.output(print(res2));

	plot_text(c(
		"WARNING: The adonis2 output could not be merged into prior format.",
		"The handling of degenerate variables could explain the discrepancy.",
		"You may need to remove/curate the model to remove factor coefficients that could not be estimated.",
		"",
		"-----------------------------------------------------------------------------------------------",
		"",
		"Original Adonis:",
		"",
		adonis_out,
		"",
		"-----------------------------------------------------------------------------------------------",
		"",
		"Adonis2:",
		"",
		adonis2_out
	));

	res=old_res;
}

print(res);

cat("--------------------------------------------------------------------------\n");
cat("\n\n");

if(0){
	cat("\ncall:\n");
	print(res$call);
	cat("\ncoefficients:\n");
	print(res$coefficients);
	cat("\ncoef.sites:\n");
	print(res$coef.sites);
	#cat("\nx.perms:\n");
	#print(res$f.perms);
	cat("\nmodel.matrix:\n");
	print(res$model.matrix);
	cat("\nterms:\n");
	print(res$terms);
}

# Not sure if these are the residuals or what the scale are.  But they may be a good proxy.

cat("\nFactors:\n");
print(factors);

cat("\nModel Matrix: \n");
cat("The factors/levels are coded into this model matrix.\n\n");
print(res$model.matrix);
#dim(res$model.matrix);

cat("\nCoefficients: \n");
cat("For each sample a linear model is fit for all the distances to that sample of interest.\n");
cat("It is 'multivariate' in the number of samples (n), not the number of categories.\n");
cat("The number of categories is lost anyway, because this is just a (n x n) distance matrix.\n\n");

print(res$coef.sites);
#dim(res$coef.sites);

cat("\nUnexplained distance between samples:\n");
cat("Multiplying the model matrix by the estimated coefficients, would equal 0 for each sample\n");
cat("if all the components of inter sample distances for that sample of interest could be accounted for.\n");
cat("We don't really have residuals for each sample, but these values capture their essence.\n\n");

#cat("Coefficients:\n");
#print(res$coef.sites);

#cat("Model Matrix:\n");
#print(t(res$model.matrix));

cat("Removing coefficients that were not estimable...\n");
non_na_coeff=apply(res$coef.sites, 1, function(x){all(!is.na(x))});

non_na_coef.sites=res$coef.sites[non_na_coeff,,drop=F];

cat("Calculating residuals...\n");
fit=non_na_coef.sites * t(res$model.matrix);
variation=apply(fit, 2, function(x){ abs(sum(x))});
names(variation)=rownames(res$model.matrix);
variation_order=order(variation, decreasing=T);

print(variation[variation_order]);
cat("\n");

##############################################################################

##############################################################################
# Output summary input and results

# Output factor summaries

used_factors=intersect(model_var, factor_names);

factor_summary=capture.output(summary(factors[used_factors]));
out_text=c(
	"PERMANOVA analysis for:",
	OutputFnameRoot, 
	"", "", 
	paste("Num Samples Used: ", num_samples, sep=""), 
	"", "", 
	"Factor Level Summary:",
	"",
	factor_summary, 
	""
);
plot_text(out_text);

# Output model and ANOVA table
anova_lines=capture.output(print(res$aov.tab));
out_text=c(
	"Model: ",
	paste("    ", model_string),
	"",
	"Stratified Resampling (Blocking): ",
	paste("    ", ifelse(Blocking!="", Blocking, "No blocking performed.")),
	"", "",
	anova_lines
);
plot_text(out_text);

##############################################################################
# Plot SS barplots

plot_sumsqr_barplot=function(adonis_res){

	#print(res);
	#print(names(res));
	#print(res[["aov.tab"]]);
	#print(rownames(res[["aov.tab"]]));
	#print(colnames(res[["aov.tab"]]));

	orig_par=par(no.readonly=T);

	ss_val=adonis_res[["aov.tab"]][,"R2"];
	pval=adonis_res[["aov.tab"]][,"Pr(>F)"];
	pval_lt10=(pval<=.10 & !is.na(pval));

	ss_varname=rownames(adonis_res[["aov.tab"]]);
	
	signf_varnames=ss_varname[pval_lt10];
	cat("Significant variable names:\n");
	print(signf_varnames);

	res_ix=which(ss_varname=="Residuals");
	ss_varname[res_ix]="\"Unexplained\"";

	ss_arr_len=length(ss_val)-1; # exclude Total

	ss_val=ss_val[1:ss_arr_len];
	ss_varname=ss_varname[1:ss_arr_len];
	ss_val_order=order(ss_val, decreasing=T);
	ss_val=ss_val[ss_val_order];
	ss_varname=ss_varname[ss_val_order];


	print(ss_val);
	print(ss_varname);

	par(mar=c(10,5,5,10));
	barcol=rep("grey", ss_arr_len);
	textcol=rep("grey33", ss_arr_len);

	names(barcol)=ss_varname;
	names(textcol)=ss_varname;

	barcol[signf_varnames]="blue";
	barcol["\"Unexplained\""]="red";

	textcol[signf_varnames]="black";
	textcol["\"Unexplained\""]="darkred";

	mids=barplot(ss_val, main="R^2 By Factor", xlab="", las=2, col=barcol,
		ylim=c(0,1.1), ylab="Proportion of Sum of Squares (SS)");

	text(mids, ss_val, sprintf("%2.3f", ss_val), pos=3, cex=.7, col=textcol);

        bar_width=mids[2]-mids[1];
        plot_range=par()$usr;
        label_size=min(c(1,.7*bar_width/par()$cxy[1]));
        text(
		mids-par()$cxy[1]/2, 
		rep(-par()$cxy[2]/2, ss_arr_len), 
		ss_varname, srt=-45, xpd=T, pos=4, cex=label_size, col=textcol
	);

	par(orig_par);

	return;
}

plot_sumsqr_barplot(res);

##############################################################################
# Pre-Compute PCA and MDS

# Compute PCA
cor_mat=cor(distmat);
eigen_out=eigen(cor_mat);

PC_contributions=eigen_out$values/sum(eigen_out$values);
pc1=eigen_out$vectors[,1];
pc2=eigen_out$vectors[,2];

# Compute MDS
for(i in 1:length(distmat)){
        if(distmat[i]==0){
                distmat[i]=1e-323;
        }
}
#mds=isoMDS(distmat);
#mds1=mds$points[,1];
#mds2=mds$points[,2];

clmds=cmdscale(distmat);
mds1=clmds[,1];
mds2=clmds[,2];

##############################################################################
# Define page layout for the analyses

layout_mat=matrix(c(
	1,1,1,2,2,2,3,
	1,1,1,2,2,2,3,
	1,1,1,2,2,2,3
), byrow=T, ncol=7);

variations_layout_mat=matrix(c(
	1,1,1,1,1,1,2,
	1,1,1,1,1,1,2,
	1,1,1,1,1,1,2
), byrow=T, ncol=7);

variation_comparison_layout_mat=matrix(c(
	1,1,1,2,2,2,	
	1,1,1,2,2,2,	
	1,1,1,2,2,2
), byrow=T, ncol=6);


##############################################################################
# Plot with points colored by group

simple_colors=c(
	"blue", "red", "green", "orange", "violet", "pink", "deepskyblue", "black");
num_simple_colors=length(simple_colors);

##############################################################################

factor_sample_names=rownames(factors);
num_factor_sample_names=length(factor_sample_names);

anova_terms=rownames(res$aov.tab);
non_variable=which(anova_terms=="Residuals" | anova_terms=="Total");
anova_terms=anova_terms[-non_variable];
num_anova_terms=length(anova_terms);

mm_var_names=colnames(res$model.matrix);
print(mm_var_names);

bin_continuous_values=function(values, num_bins=10){
	minv=min(values);	
	maxv=max(values);
	range=maxv-minv;
	# Map values between 0 and 1
	prop=(values-minv)/range;
	# Scale value up to bin, and round, to quantize
	closest=round(prop*num_bins,0);
	log10range=log10(range);
	trunc=signif(closest/num_bins*range+minv, 5)
	# Remap values to original range and location
	return(trunc);
}

par(oma=c(0,0,4,0));

for(fact_id in 1:num_anova_terms){

	term_name=anova_terms[fact_id];

	cat("Working on: ", term_name, "\n");

	# Get ANOVA information
	df=res$aov.tab[term_name,"Df"];
	SS=res$aov.tab[term_name,"SumsOfSqs"];
	MS=res$aov.tab[term_name,"MeanSqs"];
	F=res$aov.tab[term_name,"F.Model"];
	R2=res$aov.tab[term_name,"R2"];
	pval=res$aov.tab[term_name,"Pr(>F)"];
	eta_sqrd=SS/res$aov.tab["Total","SumsOfSqs"]; # eta^2: .1 small, .25 medium, .4 large

	signf_char="";
	if(pval<.001){
		signf_char=" ***";
	}else if(pval<0.01){
		signf_char=" **";
	}else if(pval<0.05){
		signf_char=" *";
	}else if(pval<0.10){
		signf_char=" .";
	}

	cur_factor=term_name;
	if(length(grep(":", term_name))){

		# Break down main effects from term name, F1:F2
		main_eff=strsplit(term_name, ":")[[1]];

		# Build unique identifer for F1 x F2 combination
		samp_to_crossfact=apply(factors[,main_eff], 1, function(x){
			str=paste(as.character(x), collapse="_x_");
			return(str);
		});

		# Determine number of unique combinations
		unique_crosses=sort(unique(samp_to_crossfact));
		num_unique_crosses=length(unique_crosses);

		cur_factor=samp_to_crossfact;
		factor_levels=unique_crosses;	
		num_levels=num_unique_crosses;	

	}else{
		# Get factor information
		cur_factor=factors[,term_name];
		is_factor=is.factor(cur_factor);
		is_ordered=is.ordered(cur_factor);
		names(cur_factor)=rownames(factors);

		if(!is_factor){
			num_unique=length(unique(cur_factor));
			cat("Number of unique 'levels':", num_unique, "\n");
			if(num_unique>2){

				if(is.numeric(cur_factor)){
					cur_factor=bin_continuous_values(cur_factor, num_bins=10);
				}

				factor_levels=sort(unique(cur_factor));
			}
			factor_levels=sort(unique(cur_factor));
		}else{
			factor_levels=unique(cur_factor);
		}
		num_levels=length(factor_levels);

	}

	cat("Factor Levels for ", term_name, "\n");
	cat("Is Factor? ", is_factor, "\n");
	cat("Is Ordered? ", is_ordered, "\n");
	print(factor_levels);
	cat("\n");

	# allocate/assign colors to palette
	if(num_levels>num_simple_colors || !is_factor){
		palette(rev(rainbow(num_levels, start=0, end=4/6)));
	}else{
		palette(simple_colors);
	}

	# Map factor levels to colors
	fact_col=rep(-1,num_samples);
	names(fact_col)=sample_names;
	for(i in 1:num_samples){
		fact_col[sample_names[i]]=which(factor_levels==cur_factor[sample_names[i]]);
	}	

	# Set up layout
	layout(layout_mat);

	XPAD=0.15;
	YPAD=0.05

	# Plot PCA
	xrange=range(pc1); xspan=abs(diff(xrange));
	yrange=range(pc2); yspan=abs(diff(yrange));
	plot(pc1, pc2, type="n", 
		xlab=sprintf("PC 1 (%3.1f%%)", PC_contributions[1]*100), 
		ylab=sprintf("PC 2 (%3.1f%%)", PC_contributions[2]*100),
		main=sprintf("PCA: (%3.1f%%)", (PC_contributions[1]+PC_contributions[2])*100),
		xlim=c(xrange[1]-XPAD*xspan, xrange[2]+XPAD*xspan),
		ylim=c(yrange[1]-YPAD*yspan, yrange[2]+YPAD*yspan)
	);
	text(pc1, pc2, labels=sample_names, cex=.7, col=fact_col[sample_names]);

	# Plot MDS
	xrange=range(mds1); xspan=abs(diff(xrange));
	yrange=range(mds2); yspan=abs(diff(yrange));
	plot(mds1, mds2, type="n",
		xlab="Dim 1",
		ylab="Dim 2",
		main="MDS",
		xlim=c(xrange[1]-XPAD*xspan, xrange[2]+XPAD*xspan),
		ylim=c(yrange[1]-YPAD*yspan, yrange[2]+YPAD*yspan)
	);
	text(mds1, mds2, labels=sample_names, cex=.7, col=fact_col[sample_names]);

	# Plot Legend
	mar=par()$mar;
	par(mar=c(0,0,0,0));
	plot(0,0, type="n", xlim=c(0,10), ylim=c(0,10), ylab="", xlab="", xaxt="n", yaxt="n", bty="n");
	legend(0,9, legend=factor_levels, fill=1:num_levels, bty="n", title=term_name);
	text(5,2, sprintf("df = %i\nSS = %5.4f\nMS = %5.4f\nF = %5.4f\nR^2 = %5.4f\np-value = %5.4f%s\neta^2 = %5.4f", 
		df, SS, MS, F, R2, pval, signf_char, eta_sqrd));
	par(mar=mar);

	# Label this page with the input file name
	mtext(OutputFnameRoot, side=3, outer=T, line=1.5, cex=1.2, font=2);

	###############################################################################
	# Plot Reoriented MDS with centroids

	# Reorient points so that second factor level is on the the right of the first factor level
	reoriented=orient_points_by_centroid(mds1, mds2, fact_col);
	mds1_reori=reoriented$x;
	mds2_reori=reoriented$y;
	mds1_centoid=reoriented$x_centroids;
	mds2_centoid=reoriented$y_centroids;

	if(is.finite(Xrange[1])){
		xrange=Xrange;
	}else{
		xrange=range(mds1_reori); xspan=abs(diff(xrange));
	}

	if(is.finite(Yrange[1])){
		yrange=Yrange;
	}else{
		yrange=range(mds2_reori); yspan=abs(diff(yrange));
	}

	# Plot oriented with labels
	plot(mds1_reori, mds2_reori, type="n",
		xlab="Dim 1",
		ylab="Dim 2",
		main="Rotated MDS: Samples Labeled",
		xlim=xrange,
		ylim=yrange
	);
	text(mds1_reori, mds2_reori, labels=sample_names, cex=.7, col=fact_col[sample_names]);

	# Plot reoriented with glyphs
	plot(mds1_reori, mds2_reori, type="n",
		xlab="Dim 1",
		ylab="Dim 2",
		main="Rotated MDS: Centroids Labeled",
		xlim=xrange,
		ylim=yrange
	);

	if(num_samples>100){
		pt_size=.5;
	}else if(num_samples>50){
		pt_size=.75;
	}else{
		pt_size=1;
	}
	points(mds1_reori, mds2_reori, cex=pt_size, col=fact_col[sample_names]);

	# bull eye
	points(mds1_centoid, mds2_centoid, cex=2, col=1:num_levels, pch=19);
	points(mds1_centoid, mds2_centoid, cex=2, col="black", pch=21);

	# Only plot category labels if the labels are factors, i.e. no continuous values
	if(is_factor || is_ordered){
		text(mds1_centoid, mds2_centoid, labels=factor_levels, cex=1.2, font=2, pos=1);
	}


	# Plot Legend
	mar=par()$mar;
	par(mar=c(0,0,0,0));
	plot(0,0, type="n", xlim=c(0,10), ylim=c(0,10), ylab="", xlab="", xaxt="n", yaxt="n", bty="n");
	legend(0,9, legend=factor_levels, fill=1:num_levels, bty="n", title=term_name);
	text(5,2, sprintf("df = %i\nSS = %5.4f\nMS = %5.4f\nF = %5.4f\nR^2 = %5.4f\np-value = %5.4f%s\neta^2 = %5.4f", 
		df, SS, MS, F, R2, pval, signf_char, eta_sqrd));
	par(mar=mar);

	###############################################################################

	# Plot "residuals"
	layout(variations_layout_mat);
	par(mar=c(10, 4.1, 2.1, 2.1));
	names_by_variation=names(variation)[variation_order];
	barplot(variation[names_by_variation], names=names_by_variation, 
		col=fact_col[names_by_variation], las=2, cex.names=.6,
		main="Unexplained (Residuals) Distances");

	# Plot Legend
	mar=par()$mar;
	par(mar=c(0,0,0,0));
	plot(0,0, type="n", xlim=c(0,10), ylim=c(0,10), ylab="", xlab="", xaxt="n", yaxt="n", bty="n");
	legend(0,9, legend=factor_levels, fill=1:num_levels, bty="n", title=term_name);
	par(mar=mar);

	###############################################################################

	# Plot dispersion analyses
	factor_dispersion=compute_dispersion(variation, fact_col, factor_levels);
	#print(factor_dispersion);	
	layout(variation_comparison_layout_mat);
	ymax=max(variation);
		
	par(mar=c(10, 4.1, 4.1, 2.1));
	boxplot(factor_dispersion$points, col=1:num_levels, 
		main=paste("Dispersion Ranges by Factor Level:\n", term_name, "\n", sep=""),
		ylim=c(0, ymax),
		ylab="Unexplained (Residual) Distances",
		xaxt="n",
		cex=0
	);
	stripchart(factor_dispersion$points, vertical=T, method="jitter", add=T, pch=1, col="grey40");
	abline(h=0, col="grey");
	
	# Label levels under boxplot
	num_levels=length(factor_dispersion$points);
	level_names=names(factor_dispersion$points);
	for(i in 1:num_levels){
                text(i, -ymax*.1, level_names[i], pos=4,
                        srt=-45, xpd=T, cex=min(c(1, (35/num_levels)), pos=4));
        }

	# Plot the heat map
	plot_pval_heatmap(factor_dispersion$pvals, "Differences in Dispersion, P-values");

	###############################################################################

}

##############################################################################
# Output Log

sink(paste(OutputFnameRoot, rand, ".log.txt", sep=""));

cat("\n");
cat("Permanova Run: ", date(), "\n");
cat("\n");
cat("Distance Matrix Filename: ", DistmatFname, "\n", sep="");
cat("Factors Filename: ", FactorsFname, "\n", sep="");
cat("Output Filename Root: ", OutputFnameRoot, "\n", sep="");
cat("\n");
cat("Number of Samples: ", num_samples, "\n", sep="");
cat("\n");

num_factors=ncol(factors);
num_samples=nrow(factors);

factor_names=colnames(factors);

cat("Factor Names:\n");
cat("\t", paste(factor_names, collapse=", "), "\n", sep="");

cat("\n");
for(i in 1:num_factors){

	unique_levels=unique(factors[,i]);
	num_levels=length(unique_levels);
	cat(num_levels, " levels in \"", factor_names[i], "\"\n", sep="");
	#cat("\t", paste(unique_levels, collapse=", "), "\n", sep="");
	ftable=table(factors[,i]);
	table_hdr=names(ftable);
	for(j in 1:length(table_hdr)){
		cat("\t", table_hdr[j], ": ", ftable[j], "\n", sep="");
	}
	cat("\n");
}

print(res);

cat("\n");
sink();

##############################################################################
# Output pvalue table

fh=file(paste(OutputFnameRoot, rand, ".pval.tsv", sep=""), "w");

# Get factor names and pvalues from ANOVA table
pval=res$aov.tab[["Pr(>F)"]];
fact=rownames(res$aov.tab);

num_fact=length(fact)-2;	# 2: Residuals and Total

# Output factor names
cat(file=fh, "# Dataset");
for(i in 1:num_fact){
	cat(file=fh, "\t", fact[i], sep="");
}
cat(file=fh, "\n");

# Output p-values
cat(file=fh, OutputFnameRoot);
for(i in 1:num_fact){
	cat(file=fh, "\t", pval[i], sep="");
}
cat(file=fh, "\n");

close(fh);

##############################################################################
# Output MANOVA files

print(res$aov.tab);
num_variables=nrow(res$aov.tab)-2;
print(num_variables);

if(TagName==""){
	TagName="Variables";
}

outmat=matrix("", nrow=num_variables, ncol=3);
colnames(outmat)=c(TagName, "Pr(>F)", "Signf");
varnames=rownames(res$aov.tab);
pvals=res$aov.tab[,"Pr(>F)"];

outmat[,TagName]=varnames[1:num_variables];
outmat[,"Pr(>F)"]=sprintf("%4.4f", pvals[1:num_variables]);
outmat[,"Signf"]=sapply(pvals[1:num_variables], sig_char);

print(outmat);

write.table(outmat, file=paste(OutputFnameRoot, ".anova.summary.tsv", sep=""),
	sep="\t", quote=FALSE, col.names=T, row.names=FALSE);

##############################################################################

if(Testing){
	cat("**************************************************************\n");
	cat("*  Reminder: Testing Flag was Set...                         *\n");
	cat("**************************************************************\n");
}

cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
