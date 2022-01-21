#!/usr/bin/env Rscript

###############################################################################

library(MASS)
library('getopt');
library('vegan');

library(glm2);

source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');

DEF_DISTTYPE="man";

params=c(
	"input_summary_table", "s", 1, "character",
	"mapping_file", "m", 1, "character",
	"a_name", "a", 1, "character",
	"b_name", "b", 1, "character",
	"factor_file", "f", 1, "character",

	"factor_target_variables", "c", 1, "character",
	"collection_time_cname", "l", 1, "character",
	"which_static_value", "w", 2, "character",

	"output_filename_root", "o", 1, "character",

	"dist_type", "d", 2, "character",
	"tag_name", "t", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\n\nUsage:\n", script_name, "\n",
	"\n",
	"	-s <input summary_table.tsv file>\n",
	"	[-d <euc/wrd/man/bray/horn/bin/gow/tyc/minkp5/minkp3, default =", DEF_DISTTYPE, ">]\n",
	"\n",
	"	-m <mapping file, from A to B Sample IDs>\n",
	"	-a <column name in map file for A>\n",
	"	-b <column name in map file for B>\n",
	"\n",
	"	-f <factor file>\n",
	"	-c <factor target variables list filename>\n",
	"	-l <column name for time of sample, e.g. collection_date as integer>\n",
	"	[-w <which value to use for from sample: A or B, default=A>]\n",
	"\n",
	"	-o <output file root name>\n",
	"\n",
	"	[-t <tag name>]\n",
	"\n",
	"\n",
	"	The contents of the mapping files should be:\n",
	"	<SubjectID>\\t<SampleID_1>\\t<SampleID_2>\\n\n",
	"\n",
	"	Samples are matched by the SubjectID in the first column of the mapping file.\n",
	"\n",
	"	1.) Load summary tables\n",
	"	2.) Load Mapping files\n",
	"	3.) Compute distance between samples A and B\n",
	"	4.) Fit distance=mu*(1-exp(-t*lamba))\n",
	"	5.) Fit distance=linear_model\n",
	"	6.) Fit distance=mu*(1-exp(-t*lamba)) + linear_model\n",
	"\n",
	"\n"
);

if(
	!length(opt$input_summary_table) || 
	!length(opt$mapping_file) || 
	!length(opt$a_name) ||
	!length(opt$b_name) ||
	!length(opt$output_filename_root)
){
	cat(usage);
	q(status=-1);
}

InputSumTab=opt$input_summary_table;
MappingFile=opt$mapping_file;
Aname=opt$a_name;
Bname=opt$b_name;

OutputFileRoot=opt$output_filename_root;

FactorFile="";
TargetVariablesFname="";
CollectionTimeColname="";
WhichStaticValue="A";


if(length(opt$factor_file)){
	FactorFile=opt$factor_file;
}
if(length(opt$factor_target_variables)){
	TargetVariablesFname=opt$factor_target_variables;
}
if(length(opt$collection_time_cname)){
	CollectionTimeColname=opt$collection_time_cname;
}
if(length(opt$which_value_to_use)){
	WhichStaticValue=opt$which_value_to_use;
}

DistType=DEF_DISTTYPE;
if(length(opt$dist_type)){
	DistType=opt$dist_type;
}

if(!any(DistType== c("wrd","man","bray","horn","bin","gow","euc","tyc","minkp3","minkp5"))){
	cat("Error: Specified distance type: ", DistType, " not recognized.\n");
	quit(status=-1);
}

OutputFileRoot=paste(OutputFileRoot, ".", DistType, sep="");

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

script_options=capture.output({
cat("\n");
cat("Input Summary Table:", InputSumTab, "\n");
cat("Mapping File       :", MappingFile, "\n");
cat(" A Col             :", Aname, "\n");
cat(" B Col             :", Bname, "\n");
cat("Output File Root  : ", OutputFileRoot, "\n");
cat("\n");
cat("Distance Type      :", DistType, "\n");
cat("\n");
cat("Factor File        : ", FactorFile,"\n");
cat("  Target Variables File            : " , TargetVariablesFname, "\n");
cat("  Collection Time Colname          : " , CollectionTimeColname, "\n");
cat("  Which sample to use metadata from: " , WhichStaticValue, "\n");
});
print(script_options);

###############################################################################
# See http://www.mothur.org/wiki/Thetayc for formula

tyc_fun=function(v1, v2){
	sum_intersect=sum(v1*v2);
	sum_sqrd_diff=sum((v1-v2)^2);
	denominator=sum_sqrd_diff + sum_intersect;
	tyc=1-(sum_intersect/denominator);
	return(tyc);
}

thetaYC=function(matrix){
	
	nsamples=nrow(matrix);
	ycdist=matrix(0, ncol=nsamples, nrow=nsamples);
	for(i in 1:nsamples){
		for(j in 1:nsamples){
			if(i<j){
				ycdist[i,j]=tyc_fun(matrix[i,], matrix[j,]);
			}else{
				ycdist[i,j]=ycdist[j,i];			
			}
		}
	}
	
	as.dist(return(ycdist));
}

weight_rank_dist_opt=function(M, deg){
        NumSamples=nrow(M);
        order_matrix=matrix(0, nrow=nrow(M), ncol=ncol(M));
        for(i in 1:NumSamples){
                order_matrix[i,]=rank(M[i,], ties.method="average");
        }

        dist_mat=matrix(0, nrow=NumSamples, ncol=NumSamples);
        colnames(dist_mat)=rownames(M);
        rownames(dist_mat)=rownames(M);
        for(i in 1:NumSamples){
                for(j in 1:i){
                        dist_mat[i,j]=
                                sqrt(sum((
                                        (order_matrix[i,]-order_matrix[j,])^2)*
                                        (((M[i,]+M[j,])/2)^deg)
                                        )
                                );
                }
        }
        return(as.dist(dist_mat));

}


###############################################################################

load_summary_table=function(st_fname){
	inmat=as.matrix(read.delim(st_fname, sep="\t", header=TRUE, row.names=1, 
		check.names=FALSE, comment.char="", quote=""))

	num_categories=ncol(inmat)-1;
	num_samples=nrow(inmat);

	cat("Loaded Summary Table: ", st_fname, "\n", sep="");
	cat("  Num Categories: ", num_categories, "\n", sep="");
	cat("  Num Samples: ", num_samples, "\n", sep="");

	countsmat=inmat[,2:(num_categories+1)];

	return(countsmat);
}

#------------------------------------------------------------------------------

load_mapping_file=function(mp_fname, keep_ids, a_colname, b_colname){

	inmat=as.matrix(read.delim(mp_fname, sep="\t", header=TRUE, check.names=F, comment.char="", quote=""));

	if(ncol(inmat)==3){
		cat("\nWarning: 3 Columns Found in Mapping File.  Assuming First Column is Subject ID\n");
		cat("As Seen:\n");
		print(inmat);
		subject_ids=inmat[,1];
		inmat=inmat[,c(2,3)];
		cat("\nAs Interpretted:\n");
		print(inmat);
		cat("\n");
	}else{
		subject_ids=character(0);
	}

	# Remove unpaired
	keep_ix=apply(inmat, 1, function(x){ all(!is.na(x))});
	inmat=inmat[keep_ix, ,drop=F];
	
	if(length(subject_ids)){
		subject_ids=subject_ids[keep_ix];
	}

	# Make A the first column, B the second column
	inmat=inmat[, c(a_colname, b_colname),drop=F];

	# Keep Entry if record is in both lists
	keep_ix=c();
	orig_mat_rows=nrow(inmat);
	cat("Number of Mapping Entries Read: ", orig_mat_rows, "\n");
	for(i in 1:orig_mat_rows){
		if(any(inmat[i,1]==keep_ids) && any(inmat[i,2]==keep_ids)){
			keep_ix=c(keep_ix, i);
		}
	}
	inmat=inmat[keep_ix,];
	num_kept_matrows=nrow(inmat);
	cat("Number of Mapping Entries Kept: ", num_kept_matrows, "\n");

	ba_mapping=as.list(x=inmat[,1]);
	names(ba_mapping)=inmat[,2];
	ab_mapping=as.list(x=inmat[,2]);
	names(ab_mapping)=inmat[,1];

	coln=colnames(inmat);
	ab_name_mapping=list();
	ab_name_mapping[[coln[1]]]="a_id";
	ab_name_mapping[[coln[2]]]="b_id";

	pair_mat=matrix(NA, nrow=num_kept_matrows, ncol=2);
	pair_mat[,1]=inmat[,1];
	pair_mat[,2]=inmat[,2];
	rownames(pair_mat)=subject_ids;
	colnames(pair_mat)=c(coln[1],coln[2]);

	map_info=list();
	map_info[["subject_ids"]]=subject_ids;
	map_info[["ab_map"]]=ab_mapping;
	map_info[["ba_map"]]=ba_mapping;
	map_info[["a"]]=coln[1];
	map_info[["b"]]=coln[2];
	map_info[["name_to_ab"]]=ab_name_mapping;
	map_info[["a_id"]]=inmat[,1];
	map_info[["b_id"]]=inmat[,2];
	map_info[["num_pairs"]]=num_kept_matrows;
	map_info[["pair_matrix"]]=pair_mat;

	return(map_info);	
}

#------------------------------------------------------------------------------

normalize=function(st){
	num_samples=nrow(st);
	num_categories=ncol(st);

	normalized=matrix(0, nrow=num_samples, ncol=num_categories);
	colnames(normalized)=colnames(st);
	rownames(normalized)=rownames(st);

	sample_counts=apply(st, 1, sum);
	for(i in 1:num_samples){
		normalized[i,]=st[i,]/sample_counts[i];
	}
	return(normalized);
}

#------------------------------------------------------------------------------

compute_dist=function(norm_st, type){

	if(type=="euc"){
		dist_mat=dist(norm_st);
	}else if (type=="wrd"){
		dist_mat=weight_rank_dist_opt(norm_st, deg=4);
	}else if (type=="man"){
		dist_mat=vegdist(norm_st, method="manhattan");
	}else if (type=="bray"){
		dist_mat=vegdist(norm_st, method="bray");
	}else if (type=="horn"){
		dist_mat=vegdist(norm_st, method="horn");
	}else if (type=="bin"){
		dist_mat=vegdist(norm_st, method="bin");
	}else if (type=="gow"){
		dist_mat=vegdist(norm_st, method="gower");
	}else if (type=="tyc"){
		dist_mat=thetaYC(norm_st);
	}else if (type=="minkp3"){
		dist_mat=dist(norm_st, method="minkowski", p=1/3);
	}else if (type=="minkp5"){
		dist_mat=dist(norm_st, method="minkowski", p=1/2);
	}

	dist_mat[dist_mat==0]=1e-323;

	return(dist_mat);
}

plot_text=function(strings){

        orig_par=par(no.readonly=T);
	options(width=100);

        par(family="Courier");
        par(mar=rep(0,4));

        num_lines=length(strings);

        top=max(as.integer(num_lines), 50);

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
	
	if(orig_par$fig[1]<0){
		orig_par$fig[1]=0;
	}

	#print(orig_par);
        par(orig_par);

}

load_list=function(fn){
	cat("Loading ", fn, " as list...\n");
	data=read.delim(fn, header=F, sep="\t", quote="", comment.char="#", as.is=T);
	res=data[,1];
	cat(" Number of ID's loaded: ", length(res), "\n");
	return(res);
}

load_factors=function(fname){
	cat("Loading ", fname, "as Factor file\n");
        factors=data.frame(read.table(fname,  sep="\t", header=TRUE, row.names=1, check.names=FALSE));
        factor_names=colnames(factors);

        ignore_idx=grep("^IGNORE\\.", factor_names);

	nr=nrow(factors);
	nc=ncol(factors);
	
	cat("Rows: ", nr, " x Cols: ", nc, "\n");

        if(length(ignore_idx)!=0){
                return(factors[-ignore_idx]);
        }else{
                return(factors);
        }
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

##############################################################################

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

mask_matrix=function(val_mat, mask_mat, mask_thres, mask_val){
        masked_matrix=val_mat;
        masked_matrix[mask_mat>mask_thres]=mask_val;
        return(masked_matrix);
}

###############################################################################

pdf(paste(OutputFileRoot, ".tcpd.pdf", sep=""), height=8.5, width=8);
plot_text(script_options);

###############################################################################
# Load summary table, factor file and pairing file, and reconcile

cat("Loading summary table:", InputSumTab, "\n");
counts_mat=load_summary_table(InputSumTab);
sumtab_sample_ids=rownames(counts_mat);
num_orig_samples=nrow(counts_mat);

cat("Loading Mapping file:", MappingFile, "\n");
map_info=load_mapping_file(MappingFile, sumtab_sample_ids, Aname, Bname);
#print(map_info);

cat("Number of complete pairs: ", map_info[["num_pairs"]], "\n");
sample_ids_pairable=c(map_info[["a_id"]], map_info[["b_id"]]);

##############################################################################

cat("Extracting shared and pairable samples from counts table (summary table)...\n");
cat("Removing samples without complete mappings...\n");
counts_mat=counts_mat[sample_ids_pairable,,drop=F];

##############################################################################

data_info=capture.output({
	cat("Summary Table: ", InputSumTab, "\n");
	cat("  Num Categories: ", ncol(counts_mat), "\n");
	cat("  Num Orig Samples: ", num_orig_samples, "\n");
	cat("\n");
	cat("Num Pairable Samples: ", map_info[["num_pairs"]], "\n");
});

plot_text(data_info);
print(data_info);

##############################################################################

#load_factors
if(FactorFile!=""){
	factor_matrix=load_factors(FactorFile);
}else{
	cat("Factor File not specified.\n");
	quit();
}

if(CollectionTimeColname==""){
	cat("Collection Date/Time needs to be specified.\n");
	quit();
}

cat("Collection Time/Dates Column Name: ", CollectionTimeColname, "\n");


if(TargetVariablesFname!=""){
	target_variables_array=load_list(TargetVariablesFname);
	cat("Target Variables:\n");
	print(target_variables_array);
}else{
	cat("Target List needs to be specified.\n");
	quit();
}

if(WhichStaticValue=="A"){
	selected_factors_mat=
		factor_matrix[map_info[["a_id"]], c(CollectionTimeColname, target_variables_array), drop=F];
}else if(WhichStaticValue=="B"){
	selected_factors_mat=
		factor_matrix[map_info[["b_id"]], c(CollectionTimeColname, target_variables_array), drop=F];
}


##############################################################################

# Remove subjects with NAs in them
samples_woNA_ix=apply(selected_factors_mat, 1, function(x){ !any(is.na(x))});
selected_factors_woNA_mat=selected_factors_mat[samples_woNA_ix,,drop=F];

num_samp_wNA=nrow(selected_factors_mat);
num_samp_woNA=nrow(selected_factors_woNA_mat);

targeted_nonNA_factor_sample_ids=rownames(selected_factors_woNA_mat);
targeted_nonNA_factor_sample_ids_mate=character(num_samp_woNA);

print(targeted_nonNA_factor_sample_ids);

if(WhichStaticValue=="A"){
	cat("Using A's metdata...\n");
	for(i in 1:num_samp_woNA){
		targeted_nonNA_factor_sample_ids_mate[i]=
			map_info[["ab_map"]][[targeted_nonNA_factor_sample_ids[i]]];
	}
}else{
	cat("Using B's metdata...\n");
	for(i in 1:num_samp_woNA){
		targeted_nonNA_factor_sample_ids_mate[i]=
			map_info[["ba_map"]][[targeted_nonNA_factor_sample_ids[i]]];
	}
}


targeted_nonNA_factor_samples_comb=c(targeted_nonNA_factor_sample_ids, targeted_nonNA_factor_sample_ids_mate);

cat("Num Samples with NAs in Factor Matrix: ", num_samp_wNA, "\n");
cat("Num Samples after Subjects w/ NAs Removed:", num_samp_woNA, "\n");

# Build model string/formula
pred_string=paste(target_variables_array, collapse=" + ");

cat("Model Predictor String: \n");
print(pred_string);

sample_IDs_woNAs=rownames(selected_factors_woNA_mat);

factors_info=capture.output({
	cat("Factor File: ", FactorFile, "\n");
	cat("Collection Time Column Name: ", CollectionTimeColname, "\n");
	cat("Target Variables: \n");
	print(target_variables_array);
	cat("\n");
	cat("Original Factor Matrix: ", nrow(factor_matrix), " x ", ncol(factor_matrix), "\n");
	cat("No NA Factor Matrix for Targeted Variables: ", 
		nrow(selected_factors_mat), " x ", ncol(selected_factors_mat), "\n");
});
print(factors_info);

##############################################################################

if(WhichStaticValue=="A"){
	final_sample_pairing=cbind(targeted_nonNA_factor_sample_ids, targeted_nonNA_factor_sample_ids_mate);
}else if(WhichStaticValue=="B"){
	final_sample_pairing=cbind(targeted_nonNA_factor_sample_ids_mate, targeted_nonNA_factor_sample_ids);
}
colnames(final_sample_pairing)=c(Aname, Bname);

#print(selected_factors_woNA_mat[final_sample_pairing[,1], CollectionTimeColname] );
#print(selected_factors_woNA_mat[final_sample_pairing[,2], CollectionTimeColname] );

final_selected_factors_mat=selected_factors_woNA_mat;
final_timespans=
	factor_matrix[final_sample_pairing[,2], CollectionTimeColname] - 
	factor_matrix[final_sample_pairing[,1], CollectionTimeColname];
final_sample_pair_count=nrow(final_selected_factors_mat);

print(final_timespans);

plot_text(c(
	data_info,
	"", 
	factors_info,
	"",
	paste("Final Sample Count: ", final_sample_pair_count)
));


##############################################################################

# Normalize counts
cat("Normalizing counts...\n");
cat("Counts dimensions: \n");
print(dim(counts_mat));
norm_mat=normalize(counts_mat);

# Compute full distances
cat("Computing distances...\n");
dist_mat=compute_dist(norm_mat, DistType);
sqr_distmat=as.matrix(dist_mat);
cat("Distmat dimensions:\n");
print(dim(sqr_distmat));

# Extract paired distances
final_paired_distances=numeric(final_sample_pair_count);
for(i in 1:final_sample_pair_count){
	final_paired_distances[i]=sqr_distmat[final_sample_pairing[i,1], final_sample_pairing[i,2]];
}
#print(final_paired_distances);

##############################################################################

mds=cmdscale(dist_mat);

mds_min_x=min(mds[,1]);
mds_max_x=max(mds[,1]);
mds_min_y=min(mds[,2]);
mds_max_y=max(mds[,2]);

mds_range_x=mds_max_x-mds_min_x;
mds_range_y=mds_max_y-mds_min_y;
mds_pad_x=mds_range_x*.05;
mds_pad_y=mds_range_y*.05;

plot(0, type="n", 
	xlim=c(mds_min_x-mds_pad_x, mds_max_x+mds_pad_x),
	ylim=c(mds_min_y-mds_pad_y, mds_max_y+mds_pad_y),
	xlab="Dimension 1", ylab="Dimension 2");

for(i in 1:final_sample_pair_count){
	points(
		c(mds[final_sample_pairing[i,1],1], mds[final_sample_pairing[i,2],1]),
		c(mds[final_sample_pairing[i,1],2], mds[final_sample_pairing[i,2],2]),
		type="l",
		col="grey");
}

points(mds[final_sample_pairing[,1],1], mds[final_sample_pairing[,1],2], col="blue");
points(mds[final_sample_pairing[,2],1], mds[final_sample_pairing[,2],2], col="green");

##############################################################################

hist(final_paired_distances, main="Distribution of Paired Distances",
	xlab="Paired Distance");

##############################################################################

max_distance=max(final_paired_distances);
max_time=max(final_timespans);
plot(final_timespans, final_paired_distances,
	ylim=c(0, max_distance*1.1),
	xlim=c(0, max_time *1.1),
	xlab="Time", ylab="Paired Distance", main="Relationship between Time and Distance");

# Fit loess with 0/0
distance_w0=c(0,final_paired_distances);
timespan_w0=c(0,final_timespans);
loess_res=loess(distance_w0~timespan_w0);
loess_y=loess_res[["fitted"]];
loess_x=loess_res[["x"]];
xix=order(loess_x);
points(loess_x[xix], loess_y[xix], col="red", type="l", lwd=2, lty="dashed");

# Fit line
lmfit=lm(final_paired_distances~final_timespans);
abline(lmfit,col="green");

quit();

##############################################################################

fit_distance_rate_model=function(dist, time, factor_mat){

	cat("Inside Fit Distance Rate Model...\n");

	#print(dist);
	#print(factor_mat);

	if(length(dist)!=nrow(factor_mat)){
		cat("Error Num Distances doesn't match Num Rows in Factor Matrix.\n");
	}

	if(length(dist)!=length(time)){
		cat("Error Num Distances doesn't match Num Times.\n");
	}

	num_samples=nrow(factor_mat);
	num_factors=ncol(factor_mat);
	max_time=max(time);
	mean_dist=mean(dist);

	prev_SSD=Inf;
	obs=NULL;
	pred=NULL;

	calc_ssd=function(dist, time, factors, parameters){

		num_param=length(parameters);
		mu=parameters[1];
		lambda=parameters[2];

		prod=sweep(factors, parameters[3:num_param], MARGIN=2, FUN="*");
		linear=apply(prod, 1, sum);

		obs=dist;
		pred=mu*(1-exp(-time*lambda))+linear;

		res=list();
		res[["ssd"]]=sum((obs-pred)^2);
		res[["obs"]]=obs;
		res[["pred"]]=pred;
		return(res);		
	}
	
	# Fit rate and factors

	res_factor_mat=factor_mat;
	res_dist=dist;
	res_time=time;
	res_seq=NA;

	initial_conditions=rep(0, num_factors+2);
	initial_conditions[1]=mean_dist;
	initial_conditions[2]=mean_dist/max_time;

	print(num_factors);
	num_bootstraps=80;

	parameters=matrix(0, nrow=num_bootstraps, ncol=num_factors+2);
	colnames(parameters)=c("mu", "lambda", colnames(factor_mat));

	cat("Running Simulated Annealing for Initial Conditions:\n");

	sann_objective_fn_wfactors=function(b){
		res=calc_ssd(dist, time, factor_mat, b);
		return(res$ssd);		
	}

	sann_opt_res=optim(initial_conditions, sann_objective_fn_wfactors, method="SANN", control=list(reltol=1));
	
	for(bsit in 1:num_bootstraps){
		cat("Running BFGS for Final Estimate: [", bsit, "]\n");
		resamp_ix=sample(num_samples, replace=T);
		res_seq=resamp_ix;

		res_time=time[resamp_ix];
		res_factor_mat=factor_mat[resamp_ix,];
		res_dist=dist[resamp_ix];

		objective_fn_wfactors=function(b){
			res=calc_ssd(res_dist, res_time, res_factor_mat, b);
			return(res$ssd);		
		}

		#opt_res=optim(sann_opt_res$par, objective_fn_wfactors, method="BFGS", control=list(reltol=1e3));
		opt_res=optim(sann_opt_res$par, objective_fn_wfactors, method="BFGS", control=list(reltol=1e-4));
		print(opt_res); 
		parameters[bsit,]=opt_res$par;

	
		res=calc_ssd(res_dist, res_time, res_factor_mat, opt_res$par);

		#plot(res$obs, res$pred, type="n", xlab="observed", ylab="predicted", 
		#	main=opt_res$value, xlim=c(0,2), ylim=c(0,2));
		#text(res$obs, res$pred, res_seq);
	}

	print(parameters);

	means=apply(parameters, 2, mean);
	sds=apply(parameters, 2, sd);
	lb95=apply(parameters, 2, function(x){quantile(x, 0.025)});
	ub95=apply(parameters, 2, function(x){quantile(x, 0.975)});
	pval=apply(parameters, 2, 
		function(x){
			gtz=sum(x>0); 
			n=length(x);
			pr_nz=min(gtz/n, 1-gtz/n);
			return(min(1,pr_nz*2));
		}
	);
			
			

	cat("Means:\n");
	print(means);

	cat("StDev:\n");
	print(sds);

	print(lb95);

	print(ub95);
	
	regr_table=matrix(NA, nrow=num_factors+2, ncol=4);
	rownames(regr_table)=c("mu", "lambda", colnames(factor_mat));
	colnames(regr_table)=c("mean", "lb95", "ub95", "pval");

	regr_table[,"mean"]=means;
	regr_table[,"lb95"]=lb95;
	regr_table[,"ub95"]=ub95;
	regr_table[,"pval"]=pval;

	print(regr_table);

	plot_text(capture.output(print(regr_table)));
}

analyze_distance_time=function(distance, timespan, factors, name){

	#print(distance);
	#print(factors);
	rate=distance/timespan;

	max_time=max(timespan);
	max_dist=max(distance);
	plot(timespan, distance, xlab="Time Span", ylab="Distance",
		xlim=c(0, max_time),
		ylim=c(0, max_dist),
		main=paste("Type ", name, ": Effect of Time on Distance",  sep="")
		);
	
	# Fit loess
	loess_res=loess(distance~timespan);
	loess_y=loess_res[["fitted"]];
	loess_x=loess_res[["x"]];
	xix=order(loess_x);
	points(loess_x[xix], loess_y[xix], col="blue", type="l", lwd=4);


	#x=1:1000;
	#r=.01;
	#points(x, 1-exp(-x*r), col="pink");
	#r=.02;
	#points(x, 1-exp(-x*r), col="brown");
	#r=.03;
	#points(x, 1-exp(-x*r), col="orange");


	# Fit 1-exp(t*r), to estimate r.
	objective_fn=function(x){
		
		r=x[1];
		mu=x[2];

		obs=distance;
		pred=mu*(1-exp(-timespan*r));

		#cat("r: ", r, "\n");
		#cat("mu: ", mu, "\n");
		#cat("Obs:\n");
		#print(obs);
		#cat("Pred:\n");
		#print(pred);

		sse=sum((obs-pred)^2);
		#cat("SSE: ", sse, "\n\n");

		return(sse);
	}

	# Fit max-dist model without any metadata first.
	opt_res=optim(c(mean(distance)/max_time, mean(distance)), objective_fn, control=list(reltol=1e-60) );
	print(opt_res);
	x=1:max_time;
	r=opt_res$par[1];
	mu=opt_res$par[2];
	cat("Rate: ", r, " Mean: ", mu, "\n");
	points(x, mu*(1-exp(-x*r)), col="pink", type="b");

	hist(rate, main=paste("Type ", name, ": Distribution of Rate = Distance/Timespan"));

	

	result=fit_distance_rate_model(distance, timespan, factors);

}

analyze_distance_time(
	pairedXY_dist[subject_IDs_woNAs, "X"], 
	pairedXY_timespan[subject_IDs_woNAs, "X"],
	selected_factors_woNA_mat,
	TypenameX
	);

###############################################################################

cat("\nDone.\n")
dev.off();
warn=warnings();
if(length(warn)){
        print(warn);
}
q(status=0);

