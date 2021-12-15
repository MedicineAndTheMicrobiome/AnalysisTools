#!/usr/bin/env Rscript

###############################################################################

cat("\n\n");

library(MASS)
library('getopt');
library('vegan');

library(car);
library(glm2);

source('~/git/AnalysisTools/Metadata/RemoveNAs/Remove_NAs.r');

DEF_DISTTYPE="man";

params=c(
	"input_summary_table_X", "s", 1, "character",
	"mapping_file_X", "m", 1, "character",
	"a_name_X", "a", 1, "character",
	"b_name_X", "b", 1, "character",
	"typename_X", "n", 1, "character",

	"input_summary_table_Y", "S", 1, "character",
	"mapping_file_Y", "M", 1, "character",
	"a_name_Y", "A", 1, "character",
	"b_name_Y", "B", 1, "character",
	"typename_Y", "N", 1, "character",

	"factor_file_X", "f", 2, "character",
	"factor_file_Y", "F", 2, "character",
	"factor_target_variables", "c", 2, "character",
	"collection_time_cname", "l", 2, "character",
	"which_static_value", "w", 2, "character",


	"output_filename_root", "o", 1, "character",

	"dist_type", "d", 2, "character",
	"num_clust_cuts", "k", 2, "character",
	"tag_name", "t", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\n\nUsage:\n", script_name, "\n",
	"\n",
	"	[Sample Type X Parameters]\n",
	"	-s <input summary_table.tsv Type X file>\n",
	"	-m <mapping file, from A to B Sample IDs for Type X>\n",
	"	-a <column name in map file for A for Type X>\n",
	"	-b <column name in map file for B for Type X>\n",
	"	-n <sample type name X for display>\n",
	"\n",
	"	[Sample Type Y Parameters]\n",
	"	-S <input summary_table.tsv Type Y file>\n",
	"	-M <mapping file, from A to B Sample IDs for Type Y>\n",
	"	-A <column name in map file for A for Type Y>\n",
	"	-B <column name in map file for B for Type Y>\n",
	"	-N <sample type name Y for display>\n",
	"\n",
	"	-o <output file root name>\n",
	"\n",
	"	Options:\n",
	"	[-f <factor file for X sample IDs>]\n",
	"	[-F <factor file for Y sample IDs>]\n",
	"	[-c <factor target variables list filename>]\n",
	"	[-l <column name for time of sample, e.g. collection_date as integer]\n",
	"	[-w <which value to use for static values: AX, AY, BX, or BY>]\n",
	"\n",
	"	[-d <euc/wrd/man/bray/horn/bin/gow/tyc/minkp5/minkp3, default =", DEF_DISTTYPE, ">]\n",
	"	[-t <tag name>]\n",
	"\n",
	"\n",
	"\n",
	"\n",
	"	The contents of the mapping files should be:\n",
	"	<SubjectID>\\t<SampleID_1>\\t<SampleID_2>\\n\n",
	"\n",
	"	The -a/-b/-A/-B options will specify whether SampleID_1 or SampleID_2 should be\n",
	"	considered A or B.\n",
	"\n",
	"	The sample summary table can be specified twice for -s and -S\n",
	"\n",
	"	Samples are matched by the SubjectID in the first column of the mapping file.\n",
	"\n",
	"	1.) Load summary tables\n",
	"	2.) Load Mapping files\n",
	"	3.) Compute distance between samples A and B\n",
	"	4.) Calculate correlation between X and Y distances\n",
	"\n",
	"\n"
);

if(
	!length(opt$input_summary_table_X) || 
	!length(opt$mapping_file_X) || 
	!length(opt$a_name_X) ||
	!length(opt$b_name_X) ||
	!length(opt$typename_X) ||

	!length(opt$input_summary_table_Y) || 
	!length(opt$mapping_file_Y) || 
	!length(opt$a_name_Y) ||
	!length(opt$b_name_Y) ||
	!length(opt$typename_Y) ||

	!length(opt$output_filename_root)
){
	cat(usage);
	q(status=-1);
}

InputSumTabX=opt$input_summary_table_X;
MappingFileX=opt$mapping_file_X;
AnameX=opt$a_name_X;
BnameX=opt$b_name_X;
TypenameX=opt$typename_X;

InputSumTabY=opt$input_summary_table_Y;
MappingFileY=opt$mapping_file_Y;
AnameY=opt$a_name_Y;
BnameY=opt$b_name_Y;
TypenameY=opt$typename_Y;

OutputFileRoot=opt$output_filename_root;


FactorFileX="";
FactorFileY="";
TargetVariablesFname="";
CollectionTimeColname="";
WhichStaticValue="AX";


if(length(opt$factor_file_X)){
	FactorFileX=opt$factor_file_X;
}
if(length(opt$factor_file_Y)){
	FactorFileY=opt$factor_file_Y;
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

cat("\n");
cat("Input Summary Table X   :", InputSumTabX, "\n");
cat("Mapping File X         :", MappingFileX, "\n");
cat(" A Col Name X          :", AnameX, "\n");
cat(" B Col Name X          :", BnameX, "\n");
cat("\n");
cat("Input Summary Table Y   :", InputSumTabY, "\n");
cat("Mapping File Y         :", MappingFileY, "\n");
cat(" A Col Name Y          :", AnameY, "\n");
cat(" B Col Name Y          :", BnameY, "\n");
cat("Output File           :", OutputFileRoot, "\n");
cat("\n");
#cat("Factor File           :", FactorFile, "\n");
#cat("   (Map) Column Name  :", FactorMapColname, "\n");
#cat("Target Var File       :", TargetVariablesFile, "\n");
cat("Distance Type         :", DistType, "\n");
cat("\n");

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

find_height_at_k=function(hclust, k){
# Computes the height on the dendrogram for a particular k

        heights=hclust$height;
        num_heights=length(heights);
        num_clust=numeric(num_heights);
        for(i in 1:num_heights){
                num_clust[i]=length(unique(cutree(hclust, h=heights[i])));
        }
        height_idx=which(num_clust==k);
        midpoint=(heights[height_idx+1]+heights[height_idx])/2;
        return(midpoint);
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

get_middle_of_groups=function(clustered_leaf_names, group_asgn){
# Finds middle of each group in the plot
        num_leaves=length(group_asgn);
        groups=sort(unique(group_asgn));
        num_groups=length(groups);

        reord_grps=numeric(num_leaves);
        names(reord_grps)=clustered_leaf_names;
        reord_grps[clustered_leaf_names]=group_asgn[clustered_leaf_names];

        mids=numeric(num_groups);
        names(mids)=1:num_groups;
        for(i in 1:num_groups){
                mids[i]=mean(range(which(reord_grps==i)));
        }
        return(mids);

}

reorder_member_ids=function(members_cut, dendr_names){

	grp_mids=get_middle_of_groups(dendr_names, members_cut);

        # Reorder cluster assignments to match dendrogram left/right
        plot_order=order(grp_mids);
        mem_tmp=numeric(num_samples);
	num_cl=length(unique(members_cut));
        for(gr_ix in 1:num_cl){
		old_id=(members_cut==plot_order[gr_ix]);
		mem_tmp[old_id]=gr_ix;
        }
        names(mem_tmp)=names(members_cut);
        members_cut=mem_tmp;
	return(members_cut);
} 

remap_coord=function(x, sbeg, send, dbeg, dend){
	srang=send-sbeg;
	norm=(x-sbeg)/srang;
	drang=dend-dbeg;
	return(norm*drang+dbeg);
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


##############################################################################

output_fname_root = paste(OutputFileRoot, ".", DistType, sep="");

###############################################################################
# Load summary table, factor file and pairing file, and reconcile

cat("\n");

cat("Loading summary table X:", InputSumTabX, "\n");
counts_mat_X=load_summary_table(InputSumTabX);
sumtab_sample_ids_X=rownames(counts_mat_X);

cat("Loading Mapping file X:", MappingFileX, "\n");
map_info_X=load_mapping_file(MappingFileX, sumtab_sample_ids_X, AnameX, BnameX);

cat("\n");

cat("Loading summary table Y:", InputSumTabY, "\n");
counts_mat_Y=load_summary_table(InputSumTabY);
sumtab_sample_ids_Y=rownames(counts_mat_Y);

cat("Loading Mapping file Y:", MappingFileY, "\n");
map_info_Y=load_mapping_file(MappingFileY, sumtab_sample_ids_Y, AnameY, BnameY);

##############################################################################


cat("Sample Type X:", map_info_X[["num_pairs"]], " complete pairs.\n");
cat("Sample Type Y:", map_info_Y[["num_pairs"]], " complete pairs.\n");
cat("\n");

samples_shared_between_XY=intersect(map_info_X[["subject_ids"]], map_info_Y[["subject_ids"]]);
num_shared_XY=length(samples_shared_between_XY);
cat("Num Subjects Shared between X and Y: ", num_shared_XY, " (complete: paired pairs)\n");
print(samples_shared_between_XY);

##############################################################################

cat("Extracting pairable sample IDs based on shared subject IDs...\n");
sample_ids_pairable_X=map_info_X[["pair_matrix"]][samples_shared_between_XY,];
sample_ids_pairable_Y=map_info_Y[["pair_matrix"]][samples_shared_between_XY,];

print(sample_ids_pairable_X);
print(sample_ids_pairable_Y);

##############################################################################

cat("Extracting shared and pairable samples from counts table (summary table)...\n");

cat("Removing samples without complete mappings...\n");
counts_mat_X=counts_mat_X[sample_ids_pairable_X,,drop=F];

cat("Removing samples without complete mappings...\n");
counts_mat_Y=counts_mat_Y[sample_ids_pairable_Y,,drop=F];

##############################################################################

# Normalize counts
cat("Normalizing counts...\n");
norm_mat_X=normalize(counts_mat_X);
norm_mat_Y=normalize(counts_mat_Y);

# Compute full distances
cat("Computing distances...\n");
dist_mat_X=compute_dist(norm_mat_X, DistType);
dist_mat_Y=compute_dist(norm_mat_Y, DistType);

sqr_distmat_X=as.matrix(dist_mat_X);
sqr_distmat_Y=as.matrix(dist_mat_Y);

pairedXY_dist=matrix(NA, nrow=num_shared_XY, ncol=2);
colnames(pairedXY_dist)=c("X", "Y");
rownames(pairedXY_dist)=samples_shared_between_XY;

for(i in 1:num_shared_XY){
	a_id=sample_ids_pairable_X[i,1];
	b_id=sample_ids_pairable_X[i,2];
	pairedXY_dist[i, 1]=sqr_distmat_X[a_id, b_id];

	a_id=sample_ids_pairable_Y[i,1];
	b_id=sample_ids_pairable_Y[i,2];
	pairedXY_dist[i, 2]=sqr_distmat_Y[a_id, b_id];
}

print(pairedXY_dist);

cor_res=capture.output(print(cor.test(pairedXY_dist[,1], pairedXY_dist[,2])));

##############################################################################

pdf(paste(OutputFileRoot, ".dbl_pair_diff.pdf", sep=""), height=8.5, width=8);

max_dist=max(pairedXY_dist);

plot(pairedXY_dist[,1], pairedXY_dist[,2], xlab=TypenameX, ylab=TypenameY,
	xlim=c(0, max_dist), ylim=c(0, max_dist)
	);

plot_text(cor_res);


##############################################################################


#load_factors
if(FactorFileX!=""){
	factor_matrix_X=load_factors(FactorFileX);
}else{
	cat("Factor File for X not specified.\n");
}

if(FactorFileY!=""){
	factor_matrix_Y=load_factors(FactorFileY);
}else{
	cat("Factor File for Y not specified.\n");
}

if(TargetVariablesFname!=""){
	target_variables_array=load_list(TargetVariablesFname);
	cat("Target Variables:\n");
	print(target_variables_array);
}else{
	cat("Target List needs to be specified.\n");
}

wtype=substr(WhichStaticValue, 1,1);
wpair=substr(WhichStaticValue, 2,2);

cat("Get Static Data from Type:", wtype, "\n");
cat("Get Static Data from Member: ", wpair, "\n");

cat("Collection Time/Dates Column Name: ", CollectionTimeColname, "\n");

##############################################################################

# Calculate days/time between samples
collection_times_X=factor_matrix_X[,CollectionTimeColname,drop=F];
collection_times_Y=factor_matrix_Y[,CollectionTimeColname,drop=F];

print(collection_times_X);
print(collection_times_Y);

pairedXY_timespan=matrix(NA, nrow=num_shared_XY, ncol=2);
colnames(pairedXY_timespan)=c("X", "Y");
rownames(pairedXY_timespan)=samples_shared_between_XY;

for(i in 1:num_shared_XY){
        a_id=sample_ids_pairable_X[i,1];
        b_id=sample_ids_pairable_X[i,2];
        pairedXY_timespan[i, 1]=
		abs(collection_times_X[b_id, CollectionTimeColname]-
		collection_times_X[a_id, CollectionTimeColname]);

        a_id=sample_ids_pairable_Y[i,1];
        b_id=sample_ids_pairable_Y[i,2];
        pairedXY_timespan[i, 2]=
		abs(collection_times_Y[b_id, CollectionTimeColname]-
		collection_times_Y[a_id, CollectionTimeColname]);
}

print(pairedXY_timespan);



max_timespan=max(pairedXY_timespan);

plot(pairedXY_timespan[,1], pairedXY_timespan[,2], xlab=TypenameX, ylab=TypenameY,
	xlim=c(0, max_timespan), ylim=c(0, max_timespan),
	main="Time Spans Between Samples from Same Subject"
	);

###############################################################################

pairedXY_rates=pairedXY_dist/pairedXY_timespan;
print(pairedXY_rates);
max_rate=max(pairedXY_rates);

plot(pairedXY_rates[,1], pairedXY_rates[,2], xlab=TypenameX, ylab=TypenameY,
	xlim=c(0, max_rate), ylim=c(0, max_rate),
	main="Rate=Dist/Timespan"
	);

print(cor.test(pairedXY_rates[,1], pairedXY_rates[,2]));

cat("\nDone.\n")
dev.off();
warn=warnings();
if(length(warn)){
        print(warn);
}
q(status=0);

quit();

###############################################################################

# Export sample IDs that were pairable 
kept_sample_ids=rownames(counts_mat);
write(kept_sample_ids, file=paste(OutputFileRoot, ".used_sample_ids.txt", sep=""), 
	ncolumns=1, sep="\n");


###############################################################################


cat("Hierarchical clustering...\n");
hcl=hclust(dist_mat, method="ward.D2");

###############################################################################

param_summary=capture.output({
	cat("Input Summary Table   :", InputSumTab, "\n");
	cat("Mapping File          :", MappingFile, "\n");
	cat("Output File           :", OutputFileRoot, "\n");
	cat("Factor File           :", FactorFile, "\n");
	cat("Target Var File       :", TargetVariablesFile, "\n");
	cat("Number of Cluster Cuts:", NumClusterCuts, "\n");
	cat("Distance Type         :", DistType, "\n");
	cat("\n");
	cat("Num pairs before Factor File Reconc:", num_pairs_pre_factor, "\n");
	cat("Num pairs after Factor File Reconc:", num_pairs_post_factor, "\n");
});
plot_text(param_summary);

###############################################################################


