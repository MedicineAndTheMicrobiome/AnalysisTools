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
	"strip_samples_nas", "s", 2, "logical",
	"tag_name", "t", 2, "character",
	"bootstrap_override", "S", 2, "numeric"
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
	"	[-t <tag name>]\n",
	"\n",
	"	[-s (Flag to strip samples with NAs, default=F)]\n",
	"	[-S <bootstrap override>]\n",
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
#OutputFnameRoot=paste(OutputFnameRoot, ".perm", sep="");

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

Blocking="";
if(length(opt$blocking)){
	Blocking=opt$blocking;
	cat("Blocking variable: ", Blocking, "\n");
}

StripSamplesWithNAs=F;
if(length(opt$strip_samples_nas)){
	StripSamplesWithNAs=T;
}

RequiredFile="";
if(length(opt$required_var)){
        RequiredFile=opt$required_var;
}

BootstrapOverride=NULL;
if(length(opt$bootstrap_override)){
	BootstrapOverride=opt$bootstrap_override;
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
	factors=read.delim(fname,  header=TRUE, row.names=1, 
		stringsAsFactors=T,
		check.names=FALSE, sep="\t", quote="", comment.char="");
	# Returns data from of character strings and numbers

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

	# Compute angle between first and last factor levels
	arc=atan2(y_centroids[num_groups]-y_centroids[1], x_centroids[num_groups]-x_centroids[1]);
	
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

compute_dispersion=function(residuals, groups, group_names){

	samp_ids=names(groups);
	residuals=residuals[samp_ids];

	#cat("----------------------------------------------------\n");
	#print(residuals);
	#print(groups);
	#print(group_names);
	#cat("----------------------------------------------------\n");


	num_levels=length(group_names);	

	# Create matrix to store p-values
	pval_matrix=matrix(0, nrow=num_levels, ncol=num_levels);		
	colnames(pval_matrix)=group_names;
	rownames(pval_matrix)=group_names;

	# Store for export
	points=list();	

	# Compute pval for difference between residuals
	for(li1 in 1:num_levels){

		l1_res=residuals[groups==li1];

		# Store the grouping now, even we don't need it for this function.
		points[[group_names[li1]]]=l1_res;

		for(li2 in 1:num_levels){

			if(li1>li2){
				pval_matrix[li1, li2]=pval_matrix[li2, li1];	
			}else if(li1==li2){
				pval_matrix[li1, li2]=1;
			}else{

				#cat(cur_levels[li1], " vs. ", cur_levels[li2], "\n");
				l2_res=residuals[groups==li2];

				result=wilcox.test(l1_res, l2_res);
				pval_matrix[li1, li2]=result$p.value;
			
			}
		}
	}

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

pdf(paste(OutputFnameRoot, ".perm.pdf", sep=""), height=5.5, width=11);

plot_text(c(
	script_name,
	"",
	"",
	"Distance Matrix Filename: ",
	paste(" ", DistmatFname),
	"",
	"Factors Filename: ", 
	paste(" ", FactorsFname),
	"",
	"Output Filename Root: ",
	paste(" ", OutputFnameRoot)
));

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
# Construct and Compute PERMANOVA

dist=(as.dist(distmat));

if(ModelFormula==""){
	model_string=paste("dist ~", paste(factor_names, collapse=" + "));
}else{
	model_string=paste("dist ~", ModelFormula);
}

num_linear_components=length(strsplit(model_string, "\\+")[[1]]);

cat("\nFitting this model: ", model_string, "\n");

cat("\n--------------------------------------------------------------------------\n");

if(Blocking!=""){
	stratify=factors[[Blocking]];
}else{
	stratify=NULL;
}

# When we have fewer degrees of freedom available, run more bootstraps
min_perms=1000;
num_permutations_to_run=min_perms*max(1, 80/(num_samples-num_linear_components-1));

cat("num_perm = ", 
	min_perms, " * max(1, 80/(", num_samples, "-", num_linear_components, "-1))\n", sep="");
cat("Num Permutations to run: ", num_permutations_to_run, "\n");


if(!is.null(BootstrapOverride)){
	cat("Boostrap Overrided to: ", BootstrapOverride, "\n");
	num_permutations_to_run=BootstrapOverride;
}

#------------------------------------------------------------------------------

cat("-----------------------------------------------------------------\n");
cat("Running adonis2()...\n");
cat("-----------------------------------------------------------------\n");
adonis2_res=adonis2(as.formula(model_string), data=as.data.frame(factors), strata=stratify, 
	permutations=num_permutations_to_run, by="margin");
print(names(adonis2_res));
cat("Completed...\n");

#------------------------------------------------------------------------------

cat("-----------------------------------------------------------------\n");
cat("Running db-RDA...\n");
cat("-----------------------------------------------------------------\n");
dbrda_res=dbrda(as.formula(model_string), data=as.data.frame(factors), strata=stratify,
	permutations=num_permutations_to_run, by="margin");
print(names(dbrda_res));
cat("db-RDA completed.\n");

# These are the Residual Dissimilarities in the Original Distance Scale
estimate_residuals=function(dbrda_res_in){

	resid=residuals(dbrda_res_in, "response");
	resid_mat=as.matrix(resid);

	persamp_rms_resid_dist=apply(resid_mat, 1, function(x){
			ss=sum(x^2);
			rmss=sqrt(ss/(length(x)-1));
			return(rmss);
		});

	persamp_rms_resid_dist=sort(persamp_rms_resid_dist, decreasing=T);
	return(persamp_rms_resid_dist);
}

persamp_rms_resid_dist=estimate_residuals(dbrda_res);
print(persamp_rms_resid_dist);

##############################################################################

used_factors=intersect(model_var, factor_names);
used_factors_df=as.data.frame(factors[,used_factors]);

factor_summary=capture.output(summary(used_factors_df));
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
cat("\n\n");
print(out_text, quote=F);
cat("\n\n");

# Output model and ANOVA table
anova_lines=capture.output(print(adonis2_res));
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
cat("\n\n");
print(out_text, quote=F);
cat("\n\n");

##############################################################################

get_clean_aov_tab=function(adns_res){
	# Just grab variables where F was calculable
	adns_tab=as.data.frame(adns_res);
	Fval=adns_tab[,"F"];
	nona=!is.na(Fval);
	clean_tab=adns_tab[nona,,drop=F];
	return(clean_tab);
}
clean_permanova_tab=get_clean_aov_tab(adonis2_res);

cat("Clean PERMANOVA Tab:\n");
print(clean_permanova_tab);

##############################################################################
# Plot SS barplots

plot_sumsqr_barplot=function(clean_tab){

	#print(res);
	#print(names(res));
	#print(res[["aov.tab"]]);
	#print(rownames(res[["aov.tab"]]));
	#print(colnames(res[["aov.tab"]]));

	cat("Plotting SumSqrs Barplot...\n");

	# Extract out variables we need
	pval=clean_tab[,"Pr(>F)"];
	r2=clean_tab[,"R2"];
	varnames=rownames(clean_tab);
	names(pval)=varnames;
	names(r2)=varnames;
	num_var=nrow(clean_tab);
	signf_varnames=varnames[pval<0.1];

	unexplained_r2=1-sum(r2);

	r2_sorted=sort(r2, decreasing=T);
	r2_out=c(unexplained_r2, r2_sorted);
	names(r2_out)=c("\"Unexplained\"", names(r2_sorted));

	pval_out=c(0, pval[names(r2_sorted)]);

	#----------------------------------------------------------------------
	# Generate plot

	orig_par=par(no.readonly=T);

	par(mar=c(10,5,5,10));

	num_bars=length(r2_out);
	barcol=rep("grey", num_bars);
	textcol=rep("grey33", num_bars);
	names(barcol)=names(r2_out);
	names(textcol)=names(r2_out);

	# Color the unexplained differently
	barcol[signf_varnames]="blue";
	barcol["\"Unexplained\""]="red";
	textcol[signf_varnames]="black";
	textcol["\"Unexplained\""]="darkred";

	# Generate bars	
	mids=barplot(r2_out, main="R^2 By Factor", xlab="", las=2, col=barcol,
		ylim=c(0,1.1), ylab="Proportion of Sum of Squares (SS)",
		names.arg=""
		);

	# Label the R2 values
	text(mids, r2_out, sprintf("%2.3f", r2_out), pos=3, cex=.7, col=textcol);

	# Label the variable names below
        bar_width=mids[2]-mids[1];
        plot_range=par()$usr;
        label_size=min(c(1,.7*bar_width/par()$cxy[1]));
        text(
		mids-par()$cxy[1]/2, 
		rep(-par()$cxy[2]/2, length(r2_out)), 
		names(r2_out), srt=-45, xpd=T, pos=4, cex=label_size, col=textcol
	);

	legend(max(mids)*3/4, 1, 
		fill=c("blue", "grey"), 
		legend=c("p-value < 0.1", "Not Significant"));

	par(orig_par);
}

plot_sumsqr_barplot(clean_permanova_tab);

##############################################################################
# Pre-Compute PCA and MDS

# Precompute nMDS and MDS:
# Remember, you can do PCA on distance matrices, you have do do PCoA
# cmdscale() is classical metric MDS PCoA
# metaMDS() is non-metric MDS

metricMDS=cmdscale(distmat, eig=T);
# Returns: "points" "eig"    "x"      "ac"     "GOF"
metMDS1=metricMDS$points[,1];
metMDS2=metricMDS$points[,2];
names(metMDS1)=rownames(metricMDS$points);

eig_val=metricMDS$eig;
valid_eig_val=eig_val[eig_val>0];
sum_eig=sum(valid_eig_val);
eig_prop=valid_eig_val/sum_eig;

num_eig_gt_1pct=eig_prop>=0.01;
top_eig=eig_prop[num_eig_gt_1pct];

out_info_text=c(
	paste("Total Eigen Values: ", length(eig_val)),
	paste("Total Positive Eig Val:", length(valid_eig_val)),
	paste("Total Eig Val > 0.01: ", length(top_eig)),
	"",
	"(Eigen Values can become negative when distance matrices are not Euclidian.)");

PC_contributions=top_eig;

plot_top_dist_eign=function(topeig){
	#print(topeig);
	remaining=1-(sum(topeig));
	out_bar_hts=c(remaining, topeig);
	num_top_pcs=length(topeig);

	# Color the unexplained differently
	barcol="blue";

	mids=barplot(out_bar_hts,
		names.arg=c("Rem", 1:num_top_pcs),
		col=c("grey", rep("blue", num_top_pcs)),
		xlab="PCs", ylab="Proportion of Variance",
		ylim=c(0, max(out_bar_hts)+.05),
		main="Top metric MDS / PCoA Eigen Values of Distance Matrix");
	text(mids, out_bar_hts, sprintf("%2.3f", out_bar_hts), pos=3, cex=.6, col="dark red")

	title(main=out_info_text, line=-2, cex.main=.8, font.main=3);
	
}

plot_top_dist_eign(top_eig);
	
#------------------------------------------------------------------------------

nonMetricMDS=metaMDS(distmat, k=2);
print(names(nonMetricMDS));
nonMetMDS1=nonMetricMDS$points[,1];
nonMetMDS2=nonMetricMDS$points[,2];
names(nonMetMDS1)=rownames(nonMetricMDS$points);
#plot(nonMetMDS1, nonMetMDS2, xlim=c(-3,3), ylim=c(-1,1));

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

fitted_preds=rownames(clean_permanova_tab);
num_fitted_preds=length(fitted_preds);
cat("Fitted Predictors:\n");
print(fitted_preds);

###############################################################################

bin_continuous_values=function(values, num_bins=10){

	h=hist(values, breaks=num_bins, plot=FALSE);
	nbin=length(h$counts);

	out_names=c();
	for(i in 1:nbin){
		catnam=paste("[", 
			sprintf("% g", h$breaks[i]), 
			" - ", 
			sprintf("% g", h$breaks[i+1]), 
			"]", sep="");
		out_names=c(out_names, catnam);
	}

	#print(out_names);

	out_bin=rep(out_names[1], length(values));
	for(i in 1:nbin){
		below_ix=values>=h$breaks[i];
		out_bin[below_ix]=out_names[i];
	}

	#print(out_bin);

	return(out_bin);
	
}

#------------------------------------------------------------------------------

flatten_terms=function(term_name, factor_df, num_target_bins=10){
	# This function will look into the term/predictors values and
	#   determine how to create bins for the values.

	if(!is.data.frame(factor_df)){
		cat("Error Input Matrix must be a data frame.\n");
		quit();
	}

	cat("Flattening: ", term_name, "\n", sep="");
	components=strsplit(term_name, ":")[[1]];
	cat("Components: \n");
	print(components);
	num_components=length(components);
	num_samples=nrow(factor_df);

	# For some reason this the drop=F doesn't work.
	if(num_components>1){
		used_df=factor_df[,components];
	}else{
		used_df=as.data.frame(factor_df[,components,drop=F]);
		colnames(used_df)=components;
		rownames(used_df)=rownames(factor_df);
	}

	#print(used_df);

	# Quantize the categories if component is continuous
	used_bins=max(2, floor((num_target_bins)^(1/num_components)));

	quantized=used_df;
	is_continuous=F;
	for(i in 1:num_components){

		val=used_df[,i];
		levels=unique(val);
		num_levels=length(levels);

		if(
			is.ordered(val)||
			is.factor(val)||
			num_levels<=num_target_bins){
				quantized[,i]=val;
		}else{
			quantized[,i]=bin_continuous_values(val, used_bins);
			is_continuous=T;
		}

	}

	# Combine components together to make a string
	mapping=apply(quantized, 1, function(x){
		paste(x, collapse=":")});

	unique_combos=sort(unique(mapping));	
	num_unique_combos=length(unique_combos);
	level_ids=1:num_unique_combos;
	names(level_ids)=unique_combos;

	id_mapping=numeric(length(mapping));
	names(id_mapping)=names(mapping);
	for(i in 1:length(mapping)){
		id_mapping[i]=level_ids[mapping[i]];
	}

	results=list();
	results[["mapping"]]=mapping;
	results[["level_names"]]=unique_combos;
	results[["num_levels"]]=num_unique_combos;
	results[["level_ids"]]=level_ids;
	results[["id_mapping"]]=id_mapping;
	results[["continuous"]]=is_continuous;

	#print(mapping);
	return(results);

}


par(oma=c(0,0,4,0));

for(pred_ix in 1:num_fitted_preds){

	pred_name=fitted_preds[pred_ix];

	cat("---------------------------------------------------------\n");
	cat("Working on: ", pred_name, "\n");

	# Get ANOVA information
	df=clean_permanova_tab[pred_name, "Df"];
	Fstat=clean_permanova_tab[pred_name,"F"];
	R2=clean_permanova_tab[pred_name,"R2"];
	pval=clean_permanova_tab[pred_name,"Pr(>F)"];

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

		
	if(R2 < 0.01){
		effect_size = "Very Small";
	}else if(R2 < 0.035){
		effect_size = "Small";
	}else if(R2 < 0.06){
		effect_size = "Medium-Small";
	}else if(R2 < 0.10){
		effect_size = "Medium";
	}else if(R2 < 0.14){
		effect_size = "Medium-Large";
	}else if(R2 < 0.20){
		effect_size = "Large";
	}else{
		effect_size = "Very Large";
	}

	
	cat("df: ", df, "\n");
	cat("F: ", Fstat, "\n");
	cat("pval: ", pval, "\n");
	cat("Sgnf: ", signf_char, "\n");
	cat("R2: ", R2, "\n");
	cat("EffSize: ", effect_size, "\n");
	cat("\n");

	# Get Predictor info
	flattened_terms_res=flatten_terms(pred_name, factors);
	#print(flattened_terms_res);

	num_levels=flattened_terms_res[["num_levels"]];
	# allocate/assign colors to palette

	dark_rainbow=function(n){
		basic_rb=rev(rainbow(n, start=0, end=4/6));
		hsv = rgb2hsv(col2rgb(basic_rb))
		yellow = hsv["h", ] > 0.10 & hsv["h", ] < 0.25
		hsv["v", yellow] = hsv["v", yellow] * 0.85
		darkened_rb = hsv(
		    hsv["h", ],
		    hsv["s", ],
		    hsv["v", ]
		);
		return(darkened_rb);
	}
	
	palette(dark_rainbow(num_levels));

	#----------------------------------------------------------------------
	# Set up layout
	layout(layout_mat);

	XPAD=0.15;
	YPAD=0.05
	
	samp_cols=flattened_terms_res[["id_mapping"]];
	level_names=flattened_terms_res[["level_names"]];
	is_continuous=flattened_terms_res[["continuous"]];

	#----------------------------------------------------------------------
	# Plot nonMetric MDS 

	xrange=range(nonMetMDS1); xspan=abs(diff(xrange));
	yrange=range(nonMetMDS2); yspan=abs(diff(yrange));
	sample_names=names(nonMetMDS1);
	plot(nonMetMDS1, nonMetMDS2, type="n", 
		xlab="Dim 1",
		ylab="Dim 2",
		main="non-Metric MDS",
		xlim=c(xrange[1]-XPAD*xspan, xrange[2]+XPAD*xspan),
		ylim=c(yrange[1]-YPAD*yspan, yrange[2]+YPAD*yspan)
	);

	text(nonMetMDS1, nonMetMDS2, labels=sample_names, cex=.7, col=samp_cols[sample_names]);

	#----------------------------------------------------------------------
	# Plot metric MDS

	xrange=range(metMDS1); xspan=abs(diff(xrange));
	yrange=range(metMDS2); yspan=abs(diff(yrange));
	sample_names=names(metMDS1);
	plot(metMDS1, metMDS2, type="n",
		xlab=sprintf("Dim 1 (%3.1f%%)", PC_contributions[1]*100), 
		ylab=sprintf("Dim 2 (%3.1f%%)", PC_contributions[2]*100),
		main=sprintf("Metric MDS / PCoA: (%3.1f%%)", (PC_contributions[1]+PC_contributions[2])*100),
		xlim=c(xrange[1]-XPAD*xspan, xrange[2]+XPAD*xspan),
		ylim=c(yrange[1]-YPAD*yspan, yrange[2]+YPAD*yspan)
	);
	text(metMDS1, metMDS2, labels=sample_names, cex=.7, col=samp_cols[sample_names]);

	#----------------------------------------------------------------------
	# Plot Legend

	mar=par()$mar;
	par(mar=c(0,0,0,0));
	plot(0,0, type="n", xlim=c(0,10), ylim=c(0,10), ylab="", xlab="", xaxt="n", yaxt="n", bty="n");
	legend(0,9, legend=level_names, fill=1:num_levels, bty="n", title=pred_name);
	text(5,2, sprintf("df = %i\nF = %5.4f\np-value = %5.4f%s\nR^2 = eta^2 = %5.4f\nEffect Size = %s", 
		df, Fstat, pval, signf_char, R2, effect_size));
	par(mar=mar);

	mtext(pred_name, side=3, outer=T, line=1.5, cex=1.2, font=2);

	###############################################################################
	###############################################################################
	# Plot Reoriented MDS with centroids

	# Reorient points so that second factor level is on the the right of the first factor level
	# samp_col is also the grouping by factor levels
	reoriented=orient_points_by_centroid(metMDS1, metMDS2, samp_cols);
	mds1_reori=reoriented$x;
	mds2_reori=reoriented$y;
	mds1_centoid=reoriented$x_centroids;
	mds2_centoid=reoriented$y_centroids;

	#----------------------------------------------------------------------

	# Plot oriented with labels
	plot(mds1_reori, mds2_reori, type="n",
		xlab="Dim 1",
		ylab="Dim 2",
		main="Rotated metric MDS: Samples Labeled",
		xlim=xrange,
		ylim=yrange
	);
	text(mds1_reori, mds2_reori, labels=sample_names, cex=.7, col=samp_cols[sample_names]);

	#----------------------------------------------------------------------

	# Plot reoriented with glyphs
	plot(mds1_reori, mds2_reori, type="n",
		xlab="Dim 1",
		ylab="Dim 2",
		main="Rotated metric MDS: Centroids Labeled",
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
	points(mds1_reori, mds2_reori, cex=pt_size, col=samp_cols[sample_names]);

	# bull eye
	points(mds1_centoid, mds2_centoid, cex=1.9, col=1:num_levels, pch=19);
	points(mds1_centoid, mds2_centoid, cex=1.9, col="black", pch=21);

	if(!is_continuous){
		text(mds1_centoid, mds2_centoid, labels=level_names, cex=1.1, font=2, pos=1);
	}

	#----------------------------------------------------------------------

	# Plot Legend
	mar=par()$mar;
	par(mar=c(0,0,0,0));
	plot(0,0, type="n", xlim=c(0,10), ylim=c(0,10), ylab="", xlab="", xaxt="n", yaxt="n", bty="n");
	legend(0,9, legend=level_names, fill=1:num_levels, bty="n", title=pred_name);
	text(5,2, sprintf("df = %i\nF = %5.4f\np-value = %5.4f%s\nR^2 = eta^2 = %5.4f\nEffect Size = %s", 
		df, Fstat, pval, signf_char, R2, effect_size));
	par(mar=mar);

	###############################################################################
	###############################################################################

	# Plot "residuals"
	layout(variations_layout_mat);
	par(mar=c(10, 4.1, 2.1, 2.1));

	persamp_rms_resid_dist_decr=sort(persamp_rms_resid_dist, decreasing=T);
	names_by_decr_variation=names(persamp_rms_resid_dist_decr);

	barplot(persamp_rms_resid_dist_decr, names=names_by_decr_variation, 
		col=samp_cols[names_by_decr_variation], las=2, cex.names=.6,
		main="Unexplained (Residuals) Distances");

	# Plot Legend
	mar=par()$mar;
	par(mar=c(0,0,0,0));
	plot(0,0, type="n", xlim=c(0,10), ylim=c(0,10), ylab="", xlab="", xaxt="n", yaxt="n", bty="n");
	legend(0,9, legend=level_names, fill=1:num_levels, bty="n", title=pred_name);
	par(mar=mar);

	###############################################################################

	# Plot dispersion analyses
	factor_dispersion=compute_dispersion(persamp_rms_resid_dist_decr, samp_cols, level_names);
	#print(factor_dispersion);	
	layout(variation_comparison_layout_mat);
	ymax=max(persamp_rms_resid_dist_decr);
		
	par(mar=c(10, 4.1, 4.1, 2.1));
	boxplot(factor_dispersion$points, col=1:num_levels, 
		main=paste("Dispersion Ranges by Factor Level:\n", pred_name, "\n", sep=""),
		ylim=c(0, ymax),
		ylab="Unexplained (Residual) Distances",
		xaxt="n",
		cex=0
	);
	stripchart(factor_dispersion$points, vertical=T, method="jitter", add=T, pch=1, col="grey40");
	abline(h=0, col="grey");
	
	# Label levels under boxplot
	fd_num_levels=length(factor_dispersion$points);
	fd_level_names=names(factor_dispersion$points);
	for(i in 1:num_levels){
                text(i, -ymax*.1, level_names[i], pos=4,
                        srt=-45, xpd=T, cex=min(c(1, (35/num_levels)), pos=4));
        }

	# Plot the heat map
	plot_pval_heatmap(factor_dispersion$pvals, "Differences in Dispersion, P-values");

	###############################################################################

}


##############################################################################
# Output pvalue 

output_results=function(rootfn, tag_name, anova_tab){

	print(anova_tab);

	pred_names=rownames(anova_tab);
	pvals=anova_tab[,"Pr(>F)"];

	if(tag_name==""){
		tag_name=rootfn;
	}

	# As lines
	outfn=paste(rootfn, ".perm.pval.rows.tsv", sep="");
	fh=file(outfn, "w");
	cat(file=fh, "#", paste(c("AnalysisName", pred_names), collapse="\t"), "\n", sep="");
	cat(file=fh, paste(c(tag_name, sprintf("%5.4g", pvals)), collapse="\t"), "\n", sep="");
	close(fh);

	# As columns
	outfn=paste(rootfn, ".perm.pval.cols.tsv", sep="");

	signf_char=sapply(pvals, sig_char);
	out_tab=cbind(rownames(anova_tab), sprintf("%5.4g", pvals), signf_char);
	colnames(out_tab)=c(tag_name, "p-value", "signf_char");

	write.table(x=out_tab, file=outfn, quote=F, sep="\t",
		row.names=F, col.names=T);
}
output_results(OutputFnameRoot, TagName, clean_permanova_tab);

#------------------------------------------------------------------------------

output_anova_tab=function(rootfn, tag_name, adns_res){
	# Whole table
	outtxt=capture.output(print(adns_res, quotes=F));
	outfn=paste(rootfn, ".perm.anova_tab.txt", sep="");
	fh=file(outfn, "w");
	cat(file=fh, paste(c(tag_name, "", outtxt,""), collapse="\n"));
	close(fh);
}
output_anova_tab(OutputFnameRoot, TagName, adonis2_res);

##############################################################################

cat("--------------------------------------------------------------------------\n");
cat("Done.\n");
dev.off();

print(warnings());
q(status=0);
