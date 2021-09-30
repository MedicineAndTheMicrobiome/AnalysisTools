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
MAX_CLUSTER_CUTS=8;

params=c(
	"input_summary_table", "s", 1, "character",
	"mapping_file", "m", 1, "character",
	"output_filename_root", "o", 1, "character",
	"a_name", "A", 1, "character",
	"b_name", "B", 1, "character",

	"factor_file", "f", 2, "character",
	"factor_mapcolname", "", 2, "character",
	"factor_target_variables", "c", 2, "character",

	"dist_type", "d", 2, "character",
	"num_clust_cuts", "k", 2, "character",
	"tag_name", "t", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\n\nUsage:\n", script_name, "\n",
	"	-s <input summary_table.tsv file>\n",
	"	-m <mapping file, from A to B Identifiers>\n",
	"	-o <output file root name>\n",
	"	-A <column name in map file for A>\n",
	"	-B <column name in map file for B>\n",
	"\n",
	"	Options:\n",
	"	[-f <factor file>]\n",
	"	[--factor_mapcolname=<column name in map file (i.e. A or B) to use sample ids from>]\n",
	"	[-c <target variables>]\n",
	"\n",
	"	[-d <euc/wrd/man/bray/horn/bin/gow/tyc/minkp5/minkp3, default =", DEF_DISTTYPE, ">]\n",
	"	[-k <max clusters to cut, default=", MAX_CLUSTER_CUTS, "\n",
	"	[-t <tag name>]\n",
	"\n",
	"\n",
	"	1.) Read in summary table\n",
	"	2.) Read in factor file, keep data for 'A' samples\n",
	"	3.) Find shared samples\n",
	"	4.) Calculate distance matrix across all samples\n",
	"	5.) Cluster, and cut from 2:", MAX_CLUSTER_CUTS, "\n",
	"\n",
	"\n",
	"The format of the sample mapping file:\n",
	"<sample group A name>\\t<sample group B name>\n",
	"<id.A.1>\\t<id.B.1>\n",
	"<id.A.2>\\t<id.B.2>\n",
	"<id.A.3>\\t<id.B.3>\n",
	"...\n",
	"<id.A.N>\\t<id.B.N>\n",
	"\n",
	"\n",
	"For the distance types:\n",
	" minkp5 is the minkowski with p=.5, i.e. sum((x_i-y_i)^1/2)^2\n",
	" minkp3 is the minkowski with p=.3, i.e. sum((x_i-y_i)^1/3)^3\n",
	"\n");

if(
	!length(opt$input_summary_table) || 
	!length(opt$mapping_file) || 
	!length(opt$output_filename_root) ||
	!length(opt$a_name) ||
	!length(opt$b_name) 
){
	cat(usage);
	q(status=-1);
}

InputSumTab=opt$input_summary_table;
MappingFile=opt$mapping_file;
OutputFileRoot=opt$output_filename_root;
Aname=opt$a_name;
Bname=opt$b_name;

if(length(opt$factor_file)){
	FactorFile=opt$factor_file;
	FactorMapColname=opt$factor_mapcolname;

	if(!length(FactorMapColname)){
		cat("\n");
		cat("If you specify a factor file, then you need to specify\n");
		cat("the column name in the map file to use for the metadata.\n");
		cat("The column name for A or B, needs to be specified.\n");
		quit(status=-1);
	}

	TargetVariablesFilename=opt$factor_target_variables;
	if(!length(TargetVariablesFilename)){
		cat("\n");
		cat("You should specify a target list of variables to include\n");
		cat("if you have specified a factor/metadata file.\n");
		quit(status=-1);
	}

}else{
	FactorFile="";
	FactorMapColname="";
}


TargetVariablesFile=ifelse(length(opt$factor_target_variables), opt$factor_target_variables, "");

DistType=DEF_DISTTYPE;
if(length(opt$dist_type)){
	DistType=opt$dist_type;
}

if(!any(DistType== c("wrd","man","bray","horn","bin","gow","euc","tyc","minkp3","minkp5"))){
	cat("Error: Specified distance type: ", DistType, " not recognized.\n");
	quit(status=-1);
}

NumClusterCuts=opt$num_clust_cuts;


OutputFileRoot_nodist=OutputFileRoot;
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

cat("Input Summary Table   :", InputSumTab, "\n");
cat("Mapping File          :", MappingFile, "\n");
cat(" A Col Name           :", Aname, "\n");
cat(" B Col Name           :", Bname, "\n");
cat("Output File           :", OutputFileRoot, "\n");
cat("Factor File           :", FactorFile, "\n");
cat("   (Map) Column Name  :", FactorMapColname, "\n");
cat("Target Var File       :", TargetVariablesFile, "\n");
cat("Number of Cluster Cuts:", NumClusterCuts, "\n");
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
		inmat=inmat[,c(2,3)];
		cat("\nAs Interpretted:\n");
		print(inmat);
		cat("\n");
	}

	# Remove unpaired
	keep_ix=apply(inmat, 1, function(x){ all(!is.na(x))});
	inmat=inmat[keep_ix, ,drop=F];

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

	map_info=list();
	map_info[["ab_map"]]=ab_mapping;
	map_info[["ba_map"]]=ba_mapping;
	map_info[["a"]]=coln[1];
	map_info[["b"]]=coln[2];
	map_info[["name_to_ab"]]=ab_name_mapping;
	map_info[["a_id"]]=inmat[,1];
	map_info[["b_id"]]=inmat[,2];
	map_info[["num_pairs"]]=num_kept_matrows;

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

###############################################################################
# Load summary table, factor file and pairing file, and reconcile

output_fname_root = paste(OutputFileRoot, ".", DistType, sep="");

cat("\n");
cat("Loading summary table:", InputSumTab, "\n");
counts_mat=load_summary_table(InputSumTab);
sumtab_sample_ids=rownames(counts_mat);

cat("Loading Mapping file:", MappingFile, "\n");
map_info=load_mapping_file(MappingFile, sumtab_sample_ids, Aname, Bname);

sample_ids_pairable=c(map_info[["a_id"]], map_info[["b_id"]]);

cat("Removing samples without complete mappings...\n");
counts_mat=counts_mat[sample_ids_pairable,];

num_samples=nrow(counts_mat);
cat("Num usable samples: ", num_samples, "\n");

num_pairs_pre_factor=num_samples/2;

if(FactorFile!=""){
	factors_matrix=load_factors(FactorFile);
	factors_samp_ids=rownames(factors_matrix);

	# Grab sample ID's from A or B list
	a_or_b=map_info[["name_to_ab"]][[FactorMapColname]];
	a_or_b_samp_ids=map_info[[a_or_b]];
	cat("Using samples IDs from: ", FactorMapColname, "(", a_or_b, ")\n", sep="");

	# Keep rows with sample IDs in A or B list
	shared_ids=intersect(a_or_b_samp_ids, factors_samp_ids);
	num_shared_ids=length(shared_ids);
	factors_matrix=factors_matrix[shared_ids,,drop=F];

	# Insert sample ID's from b_or_a
	a_sample_ids=character(num_shared_ids);
	b_sample_ids=character(num_shared_ids);

	if(a_or_b == "a_id"){
		a_sample_ids=shared_ids;
		for(i in 1:num_shared_ids){
			b_sample_ids[i]=map_info[["ab_map"]][[shared_ids[i]]];
		}
	}else{
		b_sample_ids=shared_ids;
		for(i in 1:num_shared_ids){
			a_sample_ids[i]=map_info[["ba_map"]][[shared_ids[i]]];
		}
	}

	factors_matrix=cbind(factors_matrix, a_sample_ids, b_sample_ids);

	# Keep columns in target list
	target_variables=load_list(TargetVariablesFilename);
	factors_matrix=factors_matrix[,c(target_variables, "a_sample_ids", "b_sample_ids"),drop=F];

	# NA removal
	remove_na_res=remove_sample_or_factors_wNA_parallel(factors_matrix, required=target_variables,
        num_trials=640000, num_cores=64, outfile=OutputFileRoot);
	factors_matrix=remove_na_res$factors;
	target_variables=colnames(factors_matrix);
	samp_wo_nas=rownames(factors_matrix);	

	# Remove samples from counts_mat
	paired_samp_ids=c();
	for(samp_id in samp_wo_nas){
		paired_samp_ids= c(paired_samp_ids, samp_id, 
			map_info$ab_map[[samp_id]], map_info$ba_map[[samp_id]]);
	}

	map_info$a_id=intersect(paired_samp_ids, map_info$a_id);
	bid=c();
	for(aid in map_info$a_id){
		bid=c(bid, map_info$ab_map[[aid]]);
	}
	map_info$b_id=bid;
	map_info$num_pairs=length(map_info$a_id);

	counts_mat=counts_mat[paired_samp_ids,,drop=F];
	num_samples=nrow(counts_mat);

	cat("Num samples after NA removal:", num_samples, "\n");

}else{
	factors_matrix=NULL;
}

num_pairs_post_factor=num_samples/2;

###############################################################################

# Export sample IDs that were pairable and had metadata
kept_sample_ids=rownames(counts_mat);
write(kept_sample_ids, file=paste(OutputFileRoot, ".used_sample_ids.txt", sep=""), 
	ncolumns=1, sep="\n");

# Normalize counts
cat("Normalizing counts...\n");
norm_mat=normalize(counts_mat);

###############################################################################

# Compute full distances
cat("Computing distances...\n");
dist_mat=compute_dist(norm_mat, DistType);

cat("Hierarchical clustering...\n");
hcl=hclust(dist_mat, method="ward.D2");

###############################################################################

pdf(paste(OutputFileRoot, ".clust_trans.pdf", sep=""), height=8.5, width=8);

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

###############################################################################

plot_dendro_colored_by_pairtype=function(
	clus_mem, inden, aname, bname, a_ids, b_ids, k, h, grp_sizes, grp_col,
	option_name, option_args
){

	orig_par=par(no.readonly=T);

	num_pairs=length(a_ids);
	num_samples=num_pairs*2;
	sample_to_color_map=clus_mem;
	group_lines=cumsum(grp_sizes);
	group_centers=(c(0, group_lines)+c(grp_sizes/2,0))[1:k];


	label_scale=30/(num_pairs);

	color_denfun_bySample=function(n, target_ids){
		if(is.leaf(n)){
			leaf_attr=attributes(n);
			leaf_name=leaf_attr$label;

			if(any(leaf_name==target_ids)){
				ind_color=sample_to_color_map[leaf_name];
			}else{
				ind_color="grey"
			}

			if(is.null(ind_color)){
				ind_color="white";
			}

			attr(n, "nodePar") = c(leaf_attr$nodePar,
							list(lab.col=ind_color));
		}
		return(n);
	}

	text_scale_denfun=function(n){
		if(is.leaf(n)){
			leaf_attr=attributes(n);
			leaf_name=leaf_attr$label;
			attr(n, "nodePar") = c(leaf_attr$nodePar,
						cex=0,
						lab.cex=label_scale);
		}
		return(n);
	}
	
	dendro_mod=dendrapply(inden, text_scale_denfun);
	dendra_mod=dendrapply(dendro_mod, function(n){color_denfun_bySample(n, a_ids)});
	dendrb_mod=dendrapply(dendro_mod, function(n){color_denfun_bySample(n, b_ids)});

	table_sp=5;
	layout_mat=matrix(c(
		1,rep(2, table_sp),
		rep(c(3,rep(4, table_sp)), table_sp)),
		nrow=table_sp+1, byrow=T);
	layout(layout_mat);

	title_spc=3;
	samp_labl_spc=4;	
	dendro_plot_range=c(0, num_samples+1);
	palette(grp_col);

	bottom_right_margins=1;

	par(oma=c(0,0,4,0));

	# Top-Left
	par(mar=c(0,0,0,0));
	plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n",
		 main="", xlab="", ylab="", xaxt="n", yaxt="n", bty="n");
	text(.5, .5, paste("k = ", k, sep=""), cex=2, font=2);

	# Top Margin
	par(mar=c(samp_labl_spc, 0, title_spc, bottom_right_margins));
	plot(dendrb_mod, main="", horiz=F, xaxt="n", yaxt="n", xlab="", ylab="", xlim=dendro_plot_range);
	mtext(bname, side=3, font=2, line=.5);
	abline(h=h, col="blue", lty=2, lwd=.7);
	#abline(v=c(1,num_samples));

	# Right Margin
	par(mar=c(bottom_right_margins, title_spc,0, samp_labl_spc));
	plot(dendra_mod, main="", horiz=T, xaxt="n", yaxt="n", xlab="", ylab="", ylim=dendro_plot_range);
	mtext(aname, side=2, font=2, line=.5)
	abline(v=h, col="blue", lty=2, lwd=.7);
	#abline(h=c(1,num_samples));

	# Center Field
	par(mar=c(bottom_right_margins,0,0,bottom_right_margins));
	plot(0,0, xlim=dendro_plot_range, ylim=dendro_plot_range, 
		type="n", main="", xlab="", ylab="", xaxt="n", yaxt="n", bty="n");
	#abline(h=dendro_plot_range);
	#abline(v=dendro_plot_range);
	abline(h=c(0, group_lines)+.5, lwd=.5);
	abline(v=c(0, group_lines)+.5, lwd=.5);

	
	topmost=par()$usr[4];
	leftmost=par()$usr[1];
	rightmost=par()$usr[2];
	bottommost=par()$usr[3];

	# Label cluster numbers at top and left margin
	text(group_centers, rep(topmost, k), 1:k, cex=2, font=2, pos=1);
	text(rep(leftmost, k), group_centers, 1:k, cex=2, font=2, pos=4);

	if(option_name=="drawMatchingDots"){

		points(option_args[,1], option_args[,2], pch=22, col=option_args[,3]);
		mtext("Sample Pairing Scatter Plot", side=3, outer=T);

	}else if(option_name=="count.contingency"){
	
		for(rowix in 1:k){
			for(colix in 1:k){
				text(group_centers[colix], group_centers[rowix], 
					option_args[[option_name]][rowix, colix]);
			}
		}

		# bottom margin
		text(group_centers, rep(bottommost, k), option_args[["count.marginal_B"]], 
			cex=1, font=4, pos=3);
		# right margin
		text(rep(rightmost, k), group_centers, option_args[["count.marginal_A"]],
			cex=1, font=4, pos=2);

		mtext("Contingency Table (Counts)", side=3, outer=T);

	}else if(any(option_name==c("prob.contingency", "prob.BgivenA"))){

		for(rowix in 1:k){
			for(colix in 1:k){
				text(group_centers[colix], group_centers[rowix], 
					paste(round(option_args[[option_name]][rowix, colix], 3)));
			}
		}

		if(option_name=="prob.contingency"){

			mtext(paste("Joint Prob: Pr(", aname, " & ", bname, ")", sep=""), 
				side=3, outer=T);

			# bottom margin
			text(group_centers, rep(bottommost, k), 
				sprintf("%4.3f", option_args[["prob.marginal_B"]]), 
				cex=1, font=4, pos=3);
			# right margin
			text(rep(rightmost, k), group_centers, 
				sprintf("%4.3f", option_args[["prob.marginal_A"]]),
				cex=1, font=4, pos=3, srt=90);

		}else if(option_name=="prob.BgivenA"){

			# right margin
			text(rep(rightmost, k), group_centers, 
				apply(option_args[["prob.BgivenA"]], 1, sum),
				cex=1, font=4, pos=2);

			mtext(paste("Conditional Prob: Pr(", bname, " | ", aname, ")",  sep=""), 
				side=3, outer=T);

		}

	}

	par(orig_par);
}


###############################################################################

analyze_cut=function(k, num_pairs, a_clus_mem, b_clus_mem, a_ids, b_ids){
	#print(a_clus_mem);
	#print(b_clus_mem);
	#print(a_ids);
	#print(b_ids);

	ctng_cnt=matrix(0, nrow=k, ncol=k);
	anames=paste("a_cl_", 1:k, sep="");
	bnames=paste("b_cl_", 1:k, sep="");

	rownames(ctng_cnt)=anames;
	colnames(ctng_cnt)=bnames;
	
	# Sum of contingency table
	for(i in 1:num_pairs){
		ctng_cnt[
			a_clus_mem[a_ids[i]],
			b_clus_mem[b_ids[i]]
		]=
			ctng_cnt[
				a_clus_mem[a_ids[i]],
				b_clus_mem[b_ids[i]]
			]+1;
	}

	cat("\n*********************************************************************\n");
	#cat("Contigency Table:\n");
	#print(ctng_cnt[k:1,]);
	total=sum(ctng_cnt);

	#cat("\n");
	#cat("Totals: \n");
	#print(total);

	ctng_prb=ctng_cnt/total;
	#print(ctng_prb);

	margin_A_cnt=apply(ctng_cnt, 1, sum);
	margin_B_cnt=apply(ctng_cnt, 2, sum);
	margin_A_prb=margin_A_cnt/total;
	margin_B_prb=margin_B_cnt/total;

	logdiff=log(margin_B_prb/margin_A_prb);
	propchange=(margin_B_prb-margin_A_prb)/margin_A_prb;
	prop_no_change=sum(diag(ctng_prb));	

	chisq.homogeneity.res=chisq.test(rbind(margin_A_cnt, margin_B_cnt));
	chisq.homogeneity.pval=chisq.homogeneity.res$p.value;

	#cat("\n");
	#cat("Marginal A: \n");
	#print(margin_A_cnt);
	#print(margin_A_prb);

	#cat("\n");
	#cat("Marginal B: \n");
	#print(margin_B_cnt);
	#print(margin_B_prb);

	cond_prb_BgA=matrix(0, nrow=k, ncol=k);
	cond_prb_AgB=matrix(0, nrow=k, ncol=k);
	colnames(cond_prb_BgA)=bnames;
	rownames(cond_prb_BgA)=anames;
	colnames(cond_prb_AgB)=bnames;
	rownames(cond_prb_AgB)=anames;

	# Conditional Prob: p(B|A)
	# P(A|B) = P(A^B)/P(B)
	# P(B|A) = P(A^B)/P(A)

	for(aix in 1:k){
		for(bix in 1:k){
			cond_prb_AgB[aix, bix]=ctng_prb[aix, bix] / margin_B_prb[bix];
			cond_prb_BgA[aix, bix]=ctng_prb[aix, bix] / margin_A_prb[aix];
		}
	}

	#cat("P(A|B):\n");
	#print(cond_prb_AgB);

	#cat("P(B|A):\n");
	#print(cond_prb_BgA);


	results=list();
	results[["k"]]=k;
	results[["count.totals"]]=total;
	results[["count.contingency"]]=ctng_cnt;
	results[["count.marginal_A"]]=margin_A_cnt;
	results[["count.marginal_B"]]=margin_B_cnt;
	results[["count.chisq.homo.pval"]]=chisq.homogeneity.pval;

	results[["count.remaining"]]=diag(ctng_cnt)
	results[["count.departures"]]=results[["count.marginal_A"]]-results[["count.remaining"]];
	results[["count.arrivals"]]=results[["count.marginal_B"]]-results[["count.remaining"]];
	
	results[["prob.contingency"]]=ctng_prb;
	results[["prob.marginal_A"]]=margin_A_prb;
	results[["prob.marginal_B"]]=margin_B_prb;
	results[["prob.AgivenB"]]=cond_prb_AgB;
	results[["prob.BgivenA"]]=cond_prb_BgA;
	
	results[["log_diff_perCl"]]=logdiff;
	results[["prop_chang_perCl"]]=propchange;

	results[["prop_nochange"]]=prop_no_change;

	return(results);

}

print_cut_info=function(cinfo, aname, bname){
	
	k=cinfo[["k"]];

	aclnames=c(paste(aname, "_", 1:k, sep=""), "", "Totals");
	bclnames=c(paste(bname, "_", 1:k, sep=""), "", "Totals");

	cont_tab=cinfo[["count.contingency"]];
	cont_tab=cbind(cont_tab, "", cinfo[["count.marginal_A"]]);
	cont_tab=rbind(cont_tab, "", c(cinfo[["count.marginal_B"]], "", cinfo[["count.totals"]]));
	rownames(cont_tab)=aclnames;
	colnames(cont_tab)=bclnames;

	#print(cont_tab);

	prob_tab=apply(cinfo[["prob.contingency"]], 1:2, function(x){sprintf("%3.3f", x)});
	prob_tab=cbind(prob_tab, "", sprintf("%3.3f", cinfo[["prob.marginal_A"]]));
	prob_tab=rbind(prob_tab, "", c(sprintf("%3.3f", cinfo[["prob.marginal_B"]]), "", 1));
	rownames(prob_tab)=aclnames;
	colnames(prob_tab)=bclnames;

	#print(prob_tab);

	cont_bga=apply(cinfo[["prob.BgivenA"]], 1:2, function(x){sprintf("%3.3f", x)});
	cont_bga=cbind(cont_bga, "", 1);
	rownames(cont_bga)=paste(aname, "_", 1:k, sep="");
	colnames(cont_bga)=bclnames;

	prob_staying_same=sprintf("%3.3f", diag(cinfo[["prob.BgivenA"]]));
	names(prob_staying_same)=paste("cl_", 1:k, sep="");
	
	#print(cont_bga);
	lr=cinfo[["log_diff_perCl"]];
	names(lr)=paste("cl_", 1:k, sep="");
	
	pd=cinfo[["prop_chang_perCl"]];
	names(pd)=paste("cl_", 1:k, sep="");

	out=c(
		paste("k = ", cinfo[["k"]], sep=""),
		paste("Total Pairs = ", cinfo[["count.totals"]], sep=""),
		paste("Proportion in Same Cluster = ", 
			sprintf("%3.3f", cinfo[["prop_nochange"]]), sep=""),
		paste("Chi^2 Test of Homogeneity, p-val = ", 
			sprintf("%3.3f", cinfo[["count.chisq.homo.pval"]]), sep=""),
		"",	
		"",
		"Contigency Table (Counts):",
		capture.output(print(cont_tab, quote=F)),
		"",
		"",
		"Contigency Table (Probabilities):",
		capture.output(print(prob_tab, quote=F)),
		"",
		"",
		paste("Conditional Probability: P(", bname, "|", aname, "):", sep=""),
		capture.output(print(cont_bga, quote=F)),
		"",
		"Cond. Prob. of Staying in Same Cluster:",	
		capture.output(print(prob_staying_same, quote=F)),
		"",
		"",
		paste("Difference (Positive Number is Increase of Cluster Size in ", bname, "):", sep=""),
		"",
		paste("Log Ratio, Ln(", bname, "/", aname, "):", sep=""),
		capture.output(print(sprintf("%3.3f", lr), quote=F)),
		"",
		paste("Prop Diff, (", bname, "-", aname, ")/", aname, ":", sep=""),
		capture.output(print(sprintf("%3.3f", pd), quote=F))
	);

	print(out, quote=F);

	plot_text(
		out
	);
}


add_memberships_to_factors=function(pair_map_info, factors_matrix, memberships){

	cat("Adding membership information to factor matrix...\n");
	cat("Factor Matrix Dimensions:\n");
	
	a_samp_ids=as.character(factors_matrix[,"a_sample_ids"]);
	a_memb_str=paste("cl_", memberships[a_samp_ids], sep="");
	a_cluster_id=as.factor(a_memb_str);

	b_samp_ids=as.character(factors_matrix[,"b_sample_ids"]);
	b_memb_str=paste("cl_", memberships[b_samp_ids], sep="");
	b_cluster_id=as.factor(b_memb_str);
	
	combined=cbind(as.data.frame(factors_matrix), a_cluster_id, b_cluster_id);

	return(combined);
}


fit_arrivers_logistic_regression=function(num_clusters, metadata, targ_pred, memship){
# For each cluster fit logistic regression with, 
#	post membership as response
#	pre membership and covariates as predictor

#	print(metadata);
#	print(memship);

	predictors=setdiff(targ_pred, c("a_sample_ids", "b_sample_ids", "a_cluster_id", "b_cluster_id"));
	
	regr_formla_str=paste("b_cluster_in ~ ", 
		paste(predictors, collapse= " + "), " + a_cluster_id",
		sep="");

	cat("Formula:\n");
	print(regr_formla_str);

	summ_arr_list=list();
	coeff_names=c();
	cluster_names=c();

	for(kix in 1:num_clusters){
		cat("\nFitting cluster: ", kix, " / ", num_clusters, "\n");
		target_cluster_id=paste("cl_", kix, sep="");
		cluster_names=c(cluster_names, target_cluster_id);
		b_cluster_in=as.character(metadata[,"b_cluster_id"])==target_cluster_id;
		
		cat("Percent in cluster: ", round(100*sum(b_cluster_in)/length(b_cluster_in), 2), "%\n");		
		current_df=cbind(b_cluster_in, metadata);

		resp_var=var(b_cluster_in);
		if(resp_var==0){
			cat("There is no variance in the b cluster / response.\n");
			coef_pval_mat[,"Estimate"]=0;
			coef_pval_mat[,"Pr(>|z|)"]=1;
		}else{
			values=current_df[,"a_cluster_id"];
			avail_levels=levels(values);
			if(any(target_cluster_id==avail_levels)){

				# Relevel cluster id so reference is the current cluster of interest
				current_df[,"a_cluster_id"]=relevel(values, target_cluster_id);

				#mm=model.matrix(as.formula(regr_formla_str), data=current_df);
				#num_pred_wdummies=ncol(mm);

				glm_res=glm2(as.formula(regr_formla_str), family=binomial, 
					#start=rep(0, num_pred_wdummies),
					data=current_df, singular.ok=F, trace=T);
				glm_sum=summary(glm_res);

				cat("GLM Results:\n");
				print(glm_sum);

				cat("VIF Results:\n");
				print(vif(glm_res));

				coef_pval_mat=glm_sum$coefficients[,c("Estimate", "Pr(>|z|)"), drop=F];

				if(glm_res$converged==FALSE){
					cat("Regression did not converge!!!\n");
					coef_pval_mat[,"Estimate"]=0;
					coef_pval_mat[,"Pr(>|z|)"]=1;
				}else{
					cat("Regression converged...\n");
					na_ix=is.na(coef_pval_mat[,"Pr(>|z|)"]);
					coef_pval_mat[na_ix,"Estimate"]=0;
					coef_pval_mat[na_ix,"Pr(>|z|)"]=1;
				}

			}else{
				cat("There is no ending / response members in cluster: ", target_cluster_id, "\n");
				coef_pval_mat[,"Estimate"]=0;
				coef_pval_mat[,"Pr(>|z|)"]=1;
			}
		}

		summ_arr_list[[target_cluster_id]]=coef_pval_mat;
		coeff_names=c(coeff_names, rownames(summ_arr_list[[target_cluster_id]]));

		cat("\n\n");
	}

	# Allocate matrices
	cluster_coeff_names=paste("a_cluster_id", cluster_names, sep="");
	coeff_names=sort(setdiff(unique(coeff_names), c("(Intercept)", cluster_coeff_names)));

	num_pred=length(coeff_names)+num_clusters;

	coef_matrix=matrix(0, nrow=num_pred, ncol=num_clusters);
	colnames(coef_matrix)=cluster_names;
	rownames(coef_matrix)=c(cluster_names, coeff_names);

	pval_matrix=matrix(1, nrow=num_pred, ncol=num_clusters);
	colnames(pval_matrix)=cluster_names;
	rownames(pval_matrix)=c(cluster_names, coeff_names);

	for(kix in 1:num_clusters){
		
		# Copy covariates estim/pval
		avail_coef=intersect(coeff_names, rownames(summ_arr_list[[cluster_names[kix]]]));
		for(cn in avail_coef){
			coef_matrix[cn, kix]=summ_arr_list[[cluster_names[kix]]][cn, "Estimate"];
			pval_matrix[cn, kix]=summ_arr_list[[cluster_names[kix]]][cn, "Pr(>|z|)"];
		}

		# Copy cluster estim/pval
		for(clix in 1:num_clusters){
			if(clix==kix){
				next;
			}else{
				src_clid=paste("a_cluster_id", "cl_", clix, sep="");
				dst_clid=paste("cl_", clix, sep="");

				rwn=rownames(summ_arr_list[[cluster_names[kix]]]);
				if(any(src_clid==rwn)){
					coef_matrix[dst_clid, kix]=
						summ_arr_list[[cluster_names[kix]]][src_clid, "Estimate"];
					pval_matrix[dst_clid, kix]=
						summ_arr_list[[cluster_names[kix]]][src_clid, "Pr(>|z|)"];
				}else{
					coef_matrix[dst_clid, kix]=0;
					pval_matrix[dst_clid, kix]=1;
				}
			}
		}

	}

	results=list();
	results[["coef"]]=coef_matrix;
	results[["pval"]]=pval_matrix;
	return(results);


}

fit_departers_logistic_regression=function(num_clusters, metadata, targ_pred, memship){
# For each a cluster fit logistic regression with, 
#	only it's pre-members
#	pre membership and covariates as predictor

#	print(metadata);
#	print(memship);

	print(colnames(metadata));

	predictors=setdiff(targ_pred, c("a_sample_ids", "b_sample_ids"));
	regr_formla_str=paste("b_cluster_out ~ ", 
		paste(predictors, collapse= " + "),
		sep="");

	cat("\nFormula:\n");
	print(regr_formla_str);

	summ_arr_list=list();
	coeff_names=c();
	cluster_names=c();

	for(kix in 1:num_clusters){

		cat("\nFitting cluster: ", kix, "\n");
		
		target_cluster_id=paste("cl_", kix, sep="");		
		cluster_names=c(cluster_names, target_cluster_id);

		cl_members_ix=as.character(metadata[,"a_cluster_id"])==target_cluster_id;
		cl_metadata=metadata[cl_members_ix,,drop=F];
		num_cl_members=nrow(cl_metadata);		

		if(num_cl_members>1){

			b_cluster_out=cl_metadata[,"b_cluster_id"]!=target_cluster_id;

			current_df=cbind(cl_metadata, b_cluster_out);

			glm_res=glm(as.formula(regr_formla_str), family=binomial, data=current_df);
			glm_sum=summary(glm_res);

			coef_pval_mat=glm_sum$coefficients[,c("Estimate", "Pr(>|z|)"), drop=F];
			if(glm_res$converged==FALSE){
				cat("Regression did not converge!!!\n");
				coef_pval_mat[,"Estimate"]=0;
				coef_pval_mat[,"Pr(>|z|)"]=1;
			}else{
				cat("Regression converged...\n");
				na_ix=is.na(coef_pval_mat[,"Pr(>|z|)"]);
				coef_pval_mat[na_ix,"Estimate"]=0;
				coef_pval_mat[na_ix,"Pr(>|z|)"]=1;
			}
		}else{
			cat("No starting members...\n");
			coef_pval_mat[,"Estimate"]=0;
			coef_pval_mat[,"Pr(>|z|)"]=1;
		}

		summ_arr_list[[target_cluster_id]]=coef_pval_mat;
		coeff_names=c(coeff_names, rownames(summ_arr_list[[target_cluster_id]]));

		cat("\n\n");


	}

	# Allocate matrices
	coeff_names=sort(setdiff(unique(coeff_names), "(Intercept)"));
	num_pred=length(coeff_names);

	coef_matrix=matrix(0, nrow=num_pred, ncol=num_clusters);
	colnames(coef_matrix)=cluster_names;
	rownames(coef_matrix)=coeff_names;

	pval_matrix=matrix(1, nrow=num_pred, ncol=num_clusters);
	colnames(pval_matrix)=cluster_names;
	rownames(pval_matrix)=coeff_names;


	for(kix in 1:num_clusters){
		
		# Copy covariates estim/pval
		avail_coef=setdiff(rownames(summ_arr_list[[cluster_names[kix]]]), "(Intercept)");
		for(cn in avail_coef){
			coef_matrix[cn, kix]=summ_arr_list[[cluster_names[kix]]][cn, "Estimate"];
			pval_matrix[cn, kix]=summ_arr_list[[cluster_names[kix]]][cn, "Pr(>|z|)"];
		}

	}

	#print(summ_arr_list);
	results=list();
	results[["coef"]]=coef_matrix;
	results[["pval"]]=pval_matrix;
	return(results);

}

###############################################################################

draw_arrivers_diagram=function(cinfo, mpinfo, arr_res=NULL, title="", pval_cutoff=1){

	cat("Drawing arrivers diagram...\n");
	k=cinfo$k;

	a_pos=0;
	b_pos=.5;
	margin=.05;
	top=1;
	bottom=0;
	var_space=.5;


	par(mar=c(.5,.5,6,.5));
	plot(0,0, type="n", bty="n", 
		xaxt="n", yaxt="n", 
		xlab="", ylab="",
		main=title,
		xlim=c(a_pos-margin, b_pos+var_space+margin), 
		ylim=c(bottom-margin, top+margin));

	#abline(h=c(0,1));
	#abline(v=c(0,1));

	# Draw A clusters
	clus_cent_y=seq(bottom, top, length.out=k+2)[2:(k+1)];

	b_rect_radius=(clus_cent_y[2]-clus_cent_y[1])*.95/2;
	a_rect_radius=min(b_rect_radius*.75, margin);

	# Draw lines from A to B
	for(aix in 1:k){
		for(bix in 1:k){
			lwd=10*cinfo$prob.BgivenA[aix, bix];
			if(length(lwd) && !is.nan(lwd)){
				if(lwd>0){
					points(
						c(a_pos+a_rect_radius, b_pos-b_rect_radius),
						c(clus_cent_y[aix], clus_cent_y[bix]),
						lwd=lwd,
						type="l"
					);
				}
			}
		}
	}



	#points(rep(0,k), clus_cent_y);
	
	rect(
		xleft=a_pos-a_rect_radius, 
		ybottom=clus_cent_y-a_rect_radius,
		xright=a_pos+a_rect_radius,
		ytop=clus_cent_y+a_rect_radius,
		col="white"
	);

	rect(
		xleft=b_pos-b_rect_radius, 
		ybottom=clus_cent_y-b_rect_radius,
		xright=b_pos+b_rect_radius,
		ytop=clus_cent_y+b_rect_radius,
		col="white"
	);

	for(i in 1:k){
		text(a_pos, clus_cent_y[i], i, font=2);
		text(b_pos, clus_cent_y[i], i, font=2);

		text(a_pos, clus_cent_y[i], pos=1,
			paste("[n=", cinfo$count.marginal_A[i], "]", sep=""), cex=.5);
		text(b_pos, clus_cent_y[i], pos=1,
			paste("\n(orig=", cinfo$count.remaining[i], ", new=", cinfo$count.arrivals[i], ")\n", 
			"[n=", cinfo$count.marginal_B[i], "]",
		sep=""), cex=.5);

	}

	mtext(mpinfo[["a"]], side=3, line=-2, at=a_pos, cex=2, font=2);
	mtext(mpinfo[["b"]], side=3, line=-2, at=b_pos, cex=2, font=2);

	if(!is.null(arr_res)){

		num_var_wcl=nrow(arr_res$coef);

		coef=arr_res$coef[(k+1):num_var_wcl,, drop=F];
		pval=arr_res$pval[(k+1):num_var_wcl,, drop=F];
		varnames=rownames(coef);
		num_var=length(varnames);

		#print(coef);
		#print(pval);
		
		for(clix in 1:k){

			signf_list=c();
			color_list=c();
			for(varix in 1:num_var){
				if(pval[varix, clix]<=pval_cutoff){

					if(coef[varix, clix]>0){
						dir="+";
						color_list=c(color_list, "darkgreen");
					}else{
						dir="-";
						color_list=c(color_list, "darkred");
					}

					info=paste(
						"[", dir, "] ", varnames[varix],
						" (p-val=", sprintf("%5.3f", pval[varix, clix]), ")",
						sep=""
					);
					signf_list=c(signf_list, info);
				}
			}

			cat("Signif for: ", clix, "\n");
			num_signf=length(signf_list);
			print(signf_list);

			spacing=max(num_signf, 7);
			offset=2*b_rect_radius/spacing;

			kscaling=min(4/k, 1);

			for(varix in 1:num_signf){
				text(b_pos+b_rect_radius, 
					clus_cent_y[clix]+b_rect_radius-offset*(varix-1)-b_rect_radius/8,
					signf_list[varix],
					pos=4, cex=.75*kscaling, col=color_list[varix]
					);
			}

		}

	}

	#print(cinfo);
	#print(arr_res);

	#quit();
}

###############################################################################

draw_departers_diagram=function(cinfo, mpinfo, dep_res=NULL, title="", pval_cutoff=1){

	cat("Drawing departers diagram...\n");
	k=cinfo$k;

	a_pos=0;
	margin=.5;
	top=1;
	bottom=0;
	var_space=.5;

	par(mar=c(.5,.5,6,.5));
	plot(0,0, type="n", bty="n", 
		xaxt="n", yaxt="n", 
		xlab="", ylab="",
		main=title,
		xlim=c(a_pos-margin, a_pos+var_space+margin), 
		ylim=c(bottom, top));

	# Label
	mtext(mpinfo[["a"]], side=3, line=-2, at=a_pos, cex=2, font=2);

	# Draw A clusters
	clus_cent_y=seq(bottom, top, length.out=k+2)[2:(k+1)];
	a_rect_radius=min(margin, (clus_cent_y[2]-clus_cent_y[1])*.95/2);

	cl_stay_same=diag(cinfo$prob.BgivenA);
	cl_change=1-cl_stay_same;

	# Draw lines returning to cluster
	large_arrow=.15;
	small_arrow=.5*large_arrow;
	for(aix in 1:k){


		top_y=clus_cent_y[aix]+.5*a_rect_radius;
		left_x=a_pos-2*a_rect_radius;
		bottom_y=clus_cent_y[aix]-.5*a_rect_radius;
		right_x=a_pos-a_rect_radius;

		arhd_w=a_rect_radius*.4;
		arhd_l=arhd_w*1.75;

		if(cl_stay_same[aix]>0 && !is.nan(cl_stay_same[aix])){

			stay_lwd=cl_stay_same[aix]*10;

			points(
				c(right_x, left_x, left_x, right_x-arhd_l),
				c(top_y, top_y, bottom_y, bottom_y),
				type="l", lwd=stay_lwd, col="grey50");
				
		

			# Draw arrow heads
			polygon(x=c(
					right_x,
					right_x-arhd_l,
					right_x-arhd_l,
					right_x),
				y=c(
					bottom_y,
					bottom_y+arhd_w/2,
					bottom_y-arhd_w/2,
					bottom_y),
				col="grey50", border="grey50", lwd=1
			);


		}

		text(left_x, clus_cent_y[aix], 
			sprintf("Remained:\n%3.1f%%\n[n=%i]", 100*cl_stay_same[aix], cinfo[["count.remaining"]][aix]),
			cex=.75, pos=2, font=3
		);

		left_x=a_pos+a_rect_radius;
		right_x=a_pos+2*a_rect_radius;
		bottom_y=clus_cent_y[aix]+.5*a_rect_radius;
		top_y=bottom_y;

		if(cl_change[aix]>0 && !is.nan(cl_change[aix])){

			leave_lwd=cl_change[aix]*10;

			#arrows(left_x, bottom_y-arhd_l, right_x, top_y, length=0, 
			#	angle=15, col="grey50", lwd=leave_lwd);

			points(c(left_x, right_x-arhd_l), c(bottom_y, top_y), 
				type="l", lwd=leave_lwd, col="grey50");

			# Draw arrow heads
			polygon(x=c(
					right_x,
					right_x-arhd_l,
					right_x-arhd_l,
					right_x),
				y=c(
					bottom_y,
					bottom_y+arhd_w/2,
					bottom_y-arhd_w/2,
					bottom_y),
				col="grey50", border="grey50", lwd=1
			);


		}
		text(right_x, clus_cent_y[aix]+a_rect_radius, 
			sprintf("Departures: %3.1f%% [n=%i]", 100*cl_change[aix], cinfo[["count.departures"]][aix]),
			cex=.75, pos=4, font=3
		);

	}

	# Draw rectangles for each cluster
	rect(
		xleft=a_pos-a_rect_radius, 
		ybottom=clus_cent_y-a_rect_radius,
		xright=a_pos+a_rect_radius,
		ytop=clus_cent_y+a_rect_radius,
		col="white"
	);

	# Label each rectangle
	for(i in 1:k){
		text(a_pos, clus_cent_y[i], i, font=2);
		text(a_pos, clus_cent_y[i], pos=1,
			paste("[n=", cinfo$count.marginal_A[i], "]", sep=""), cex=.5);
	}


	if(!is.null(dep_res)){

		num_var_wcl=nrow(dep_res$coef);

		coef=dep_res$coef;
		pval=dep_res$pval;

		varnames=rownames(coef);
		num_var=length(varnames);
		
		for(clix in 1:k){

			signf_list=c();
			color_list=c();
			for(varix in 1:num_var){
				if(pval[varix, clix]<=pval_cutoff){

					if(coef[varix, clix]>0){
						dir="+";
						color_list=c(color_list, "darkgreen");
					}else{
						dir="-";
						color_list=c(color_list, "darkred");
					}

					info=paste(
						"[", dir, "] ", varnames[varix],
						" (p-val=", sprintf("%5.3f", pval[varix, clix]), ")",
						sep=""
					);
					signf_list=c(signf_list, info);
				}
			}

			cat("Signif for: ", clix, "\n");
			num_signf=length(signf_list);
			print(signf_list);

			spacing=max(num_signf, 10);
			offset=2*(a_rect_radius*.9)/spacing;

			kscaling=min(4/k, 1);

			for(varix in 1:num_signf){
				text(a_pos+2*a_rect_radius, 
					(-a_rect_radius*.05)+clus_cent_y[clix]+
						a_rect_radius-offset*(varix-1)-a_rect_radius/7,
					signf_list[varix],
					pos=4, cex=.75*kscaling, col=color_list[varix]
					);
			}

		}

	}

	#print(cinfo);
	#print(arr_res);

	#quit();
}

###############################################################################
###############################################################################

test=F;

max_cuts=min(log2(num_samples), NumClusterCuts);

palette_col=c("red", "green", "blue", "cyan", "magenta", "orange", 
	"pink", "black", "purple", "brown", "aquamarine");
num_pref_col=length(palette_col);
num_clus=length(unique(max_cuts));	

if(max_cuts>num_pref_col){
	palette_col=rainbow(n=max_cuts, start=0, end=4/6);
}


orig_dendro=as.dendrogram(hcl);
lf_names=get_clstrd_leaf_names(orig_dendro);
print(lf_names);

num_pairs=map_info[["num_pairs"]];
pos_map=matrix(0, nrow=num_pairs, ncol=3);
a_id=map_info[["a_id"]];
b_id=map_info[["b_id"]];
colnames(pos_map)=c(map_info[["b"]], map_info[["a"]], "same_clus");

prop_change=rep(0,max_cuts);
chisq_homo_pval=rep(0,max_cuts);

arrivers_list=list();
departers_list=list();

#for(k in 10){
for(k in 2:max_cuts){

	
	cat("Cutting at: ", k, "\n");
	memberships=cutree(hcl, k=k);

        grp_mids=get_middle_of_groups(lf_names, memberships);

        # Reorder cluster assignments to match dendrogram left/right
        plot_order=order(grp_mids);
        mem_tmp=numeric(num_samples);
        for(gr_ix in 1:k){
                old_id=(memberships==plot_order[gr_ix]);
                mem_tmp[old_id]=gr_ix;
        }
        names(mem_tmp)=names(memberships);
        memberships=mem_tmp;
        grp_mids=grp_mids[plot_order];

	# Calculate where to draw line
	h=find_height_at_k(hcl, k);

	# Calculate group sizes
	grp_sizes_all=table(memberships);

	cut_info=analyze_cut(
		k, num_pairs,
		memberships[map_info[["a_id"]]],
		memberships[map_info[["b_id"]]],
		map_info[["a_id"]],
		map_info[["b_id"]]
	);

	# Calculate positions for pairs
	for(pix in 1:num_pairs){
		a_pos=which(lf_names==a_id[pix]);
		b_pos=which(lf_names==b_id[pix]);
		same_clus=ifelse(
			memberships[a_id[pix]]==memberships[b_id[pix]], 
			memberships[a_id[pix]],
			"grey"
			);
		pos_map[pix,]=c(b_pos, a_pos, same_clus);
	}

	plot_dendro_colored_by_pairtype(
		memberships,
		orig_dendro, 
		aname=map_info[["a"]], bname=map_info[["b"]],
		a_ids=map_info[["a_id"]], b_ids=map_info[["b_id"]],
		k, h, grp_sizes_all, palette_col,
		option_name="drawMatchingDots",
		option_args=pos_map
	);

	for(opt_name in c("count.contingency", "prob.contingency", "prob.BgivenA")){
		plot_dendro_colored_by_pairtype(
			memberships,
			orig_dendro, 
			aname=map_info[["a"]], bname=map_info[["b"]],
			a_ids=map_info[["a_id"]], b_ids=map_info[["b_id"]],
			k, h, grp_sizes_all, palette_col,
			option_name=opt_name,
			option_args=cut_info
		);
	}

	aname=map_info[["a"]];
	bname=map_info[["b"]];
	rown=paste(aname, ":\n", 1:k, sep="");
	coln=paste(bname, ":\n", 1:k, sep="");

	mat=cut_info[["count.contingency"]];
	rownames(mat)=rown;
	colnames(mat)=coln;
	paint_matrix(mat[k:1,], 
		paste("k=", k, ", Contingency Table (Counts) Heat Map", sep=""), counts=T);

	mat=cut_info[["prob.contingency"]];
	rownames(mat)=rown;
	colnames(mat)=coln;
	paint_matrix(mat[k:1,], 
		paste("k=", k, ", Joint Probability Heat Map", sep=""), counts=F, deci_pts=3,
		plot_min=0, plot_max=1);

	mat=cut_info[["prob.BgivenA"]];
	rownames(mat)=rown;
	colnames(mat)=coln;
	paint_matrix(mat[k:1,], 
		paste("k=", k, ", Conditional Pr(", map_info[["b"]], " | ", map_info[["a"]], 
			") Heat Map", sep=""), 
		counts=F, deci_pts=3,
		plot_min=0, plot_max=1);

	print_cut_info(cut_info, aname, bname);

	prop_change[k]=cut_info[["prop_nochange"]];
	chisq_homo_pval[k]=cut_info[["count.chisq.homo.pval"]];

	arrivers_res=NULL;
	departers_res=NULL;
	if(!is.null(factors_matrix)){

		factor_wmemberships=add_memberships_to_factors(map_info, factors_matrix, memberships);

		a_marginal=cut_info[["count.marginal_A"]];
		b_marginal=cut_info[["count.marginal_B"]];

		###############################################################

		# Fit and plot arrivers
		arrivers_res=fit_arrivers_logistic_regression(
			k, factor_wmemberships, target_variables, memberships);

		arrivers_list[[k]]=arrivers_res;

		# Adjust the colnames
		augmented_colnames=paste(
                        bname, ": ", colnames(arrivers_res[["coef"]]),
                        " [n=", b_marginal, "]",
                        sep="");

		colnames(arrivers_res[["coef"]])=augmented_colnames;
		colnames(arrivers_res[["pval"]])=augmented_colnames;

		# Adjust the rownames
		augmented_rownames=rownames(arrivers_res[["coef"]]);
		clix=grep("^cl_\\d", augmented_rownames);
		augmented_rownames[clix]=paste(aname, ": ", augmented_rownames[1:k], 
			" [n=", a_marginal, "]", sep="");
		rownames(arrivers_res[["coef"]])=augmented_rownames;
		rownames(arrivers_res[["pval"]])=augmented_rownames;

		# Generate all the heatmaps
		paint_matrix(arrivers_res[["coef"]], 
			title=paste("k = ", k, ", Arrivers: All Coefficients", sep=""), 
			high_is_hot=T, deci_pts=3);

		paint_matrix(arrivers_res[["pval"]], 
			title=paste("k = ", k, ", Arrivers: All P-values", sep=""), 
			high_is_hot=F, deci_pts=3, plot_min=0, plot_max=1);

		msk=mask_matrix(arrivers_res[["coef"]], arrivers_res[["pval"]], .1, 0);
		paint_matrix(msk, 
			title=paste("k = ", k, ", Arrivers: Coefficients w/ p-val < 0.10", sep=""), 
			label_zeros=F, deci_pts=3);

		msk=mask_matrix(arrivers_res[["coef"]], arrivers_res[["pval"]], .05, 0);
		paint_matrix(msk, 
			title=paste("k = ", k, ", Arrivers: Coefficients w/ p-val < 0.05", sep=""), 
			label_zeros=F, deci_pts=3);

		msk=mask_matrix(arrivers_res[["coef"]], arrivers_res[["pval"]], .01, 0);
		paint_matrix(msk, 
			title=paste("k = ", k, ", Arrivers: Coefficients w/ p-val < 0.01", sep=""), 
			label_zeros=F, deci_pts=3);


		###############################################################

		# Fit and plot leavers

		departers_res=fit_departers_logistic_regression(
			k, factor_wmemberships, target_variables, memberships);

		departers_list[[k]]=departers_res;

		# Adjust the column names
		augmented_colnames=paste(
			aname, ": ", colnames(departers_res[["coef"]]),
			" [n=", a_marginal, "]",
			sep="");
		colnames(departers_res[["coef"]])=augmented_colnames;
		colnames(departers_res[["pval"]])=augmented_colnames;
			
		# Generate all the heatmaps
		paint_matrix(departers_res[["coef"]], 
			title=paste("k = ", k, ", Departers: All Coefficients", sep=""), 
			high_is_hot=T, deci_pts=3);

		paint_matrix(departers_res[["pval"]], 
			title=paste("k = ", k, ", Departers: All P-values", sep=""), 
			high_is_hot=F, deci_pts=3, plot_min=0, plot_max=1);

		msk=mask_matrix(departers_res[["coef"]], departers_res[["pval"]], .1, 0);
		paint_matrix(msk, 
			title=paste("k = ", k, ", Departers: Coefficients w/ p-val < 0.10", sep=""), 
			label_zeros=F, deci_pts=3);

		msk=mask_matrix(departers_res[["coef"]], departers_res[["pval"]], .05, 0);
		paint_matrix(msk, 
			title=paste("k = ", k, ", Departers: Coefficients w/ p-val < 0.05", sep=""), 
			label_zeros=F, deci_pts=3);

		msk=mask_matrix(departers_res[["coef"]], departers_res[["pval"]], .01, 0);
		paint_matrix(msk, 
			title=paste("k = ", k, ", Departers: Coefficients w/ p-val < 0.01", sep=""), 
			label_zeros=F, deci_pts=3);

	}


	# Plot diagram
	draw_arrivers_diagram(cut_info, map_info, arrivers_res, title="Arrivers: P-val<0.10", pval_cutoff=.1);
	draw_arrivers_diagram(cut_info, map_info, arrivers_res, title="Arrivers: P-val<0.05", pval_cutoff=.05);
	draw_arrivers_diagram(cut_info, map_info, arrivers_res, title="Arrivers: P-val<0.01", pval_cutoff=.01);
	
	draw_departers_diagram(cut_info, map_info, departers_res, title="Departers: P-value<0.10", pval_cutoff=.1);
	draw_departers_diagram(cut_info, map_info, departers_res, title="Departers: P-value<0.05", pval_cutoff=.05);
	draw_departers_diagram(cut_info, map_info, departers_res, title="Departers: P-value<0.01", pval_cutoff=.01);

}

# Proportion not changing clusters
par(mar=c(4,4,4,1));
plot(2:k, prop_change[2:k], ylim=c(0,1), main="Proportion Not Changing Clusters", 
	ylab="Proportion", xlab="Cuts (k)", type="b");
abline(h=c(0,1), lty=2, col="blue");

# Chisq test of homogeneity 
signf_levels=c(1, .1, .05, .01, .001, .0001);
neglog_signf_levels=-log10(signf_levels);
neglog_homo_pval=-log10(chisq_homo_pval[2:k]);

par(mar=c(5, 6, 5, 6));
plot(2:k, -log10(chisq_homo_pval[2:k]), 
	ylim=range(neglog_signf_levels, neglog_homo_pval), 
	main="Chi^2 Test of Homogeneity across cuts", 
	yaxt="n",
	ylab="", xlab="Cuts (k)", type="b");
mtext("p-value", side=4, line=4);
mtext("-log10(p-value)", side=2, line=4);
abline(h=neglog_signf_levels, lty=2, col="blue");
axis(side=2, at=neglog_signf_levels, labels=sprintf("%3.2f", neglog_signf_levels), las=2);
axis(side=4, at=neglog_signf_levels, labels=signf_levels, las=2);

##############################################################################

if(!is.null(factors_matrix)){

	#print(arrivers_list);
	#print(departers_list);

	pred_names=setdiff(rownames(arrivers_list[[2]]$coef), c("cl_1", "cl_2"));
	print(pred_names);

	num_pred=length(pred_names);
	arr_mspval_mat=matrix(1, nrow=num_pred, ncol=k);
	colnames(arr_mspval_mat)=paste("cl_", 1:k, sep="");
	rownames(arr_mspval_mat)=pred_names;

	dep_mspval_mat=matrix(1, nrow=num_pred, ncol=k);
	colnames(dep_mspval_mat)=paste("cl_", 1:k, sep="");
	rownames(dep_mspval_mat)=pred_names;

	for(kix in 2:k){

		avail=intersect(pred_names, rownames(arrivers_list[[kix]]$pval));
		min_pval=apply(arrivers_list[[kix]]$pval, 1, min);
		arr_mspval_mat[avail, kix]=min_pval[avail];

		avail=intersect(pred_names, rownames(departers_list[[kix]]$pval));
		min_pval=apply(departers_list[[kix]]$pval, 1, min);
		dep_mspval_mat[avail, kix]=min_pval[avail];

	}

	cat("Arrivers:\n");
	print(arr_mspval_mat);

	cat("Departers:\n");
	print(dep_mspval_mat);

	# Function to plot most significant variables over k
	plot_pval_over_k=function(mat, title, pval_cutoff=1){
		num_var=nrow(mat);
		k=ncol(mat);

		cols=rainbow(num_var);

		# Reference lines
		pval_ref=c(0.1, 0.05, 0.01, 0.001);
		neglog10_pval_ref=-log10(pval_ref);

		neglog10_pval_mat=apply(mat, 1:2, function(x){-log10(x);});

		neglog10_pval_mat[neglog10_pval_mat==Inf]=-log10(min(pval_ref)/10);

		neglog10_pval_range=range(c(neglog10_pval_ref, neglog10_pval_mat));

		par(mar=c(4,4,5,1));
		plot(0, type="n", 
			xlim=c(1.5,k),
			ylim=neglog10_pval_range,
			main=title, xlab="k", ylab="-log10(p-value)");

		# Draw significance reference lines
		abline(h=neglog10_pval_ref, col="grey", lwd=.5, lty=2);
		text(1.5, neglog10_pval_ref, pval_ref, pos=4, font=3, cex=.75);

		# Calculate transformed cutoff so we can apply it in transformed space
		neglog10_pval_cutoff=-log10(pval_cutoff);

		for(i in 1:num_var){

			# Find most significant p-value across the k's
			most_sigf=max(neglog10_pval_mat[i,2:k]);
			most_sigf_k=min(which(most_sigf==neglog10_pval_mat[i,2:k]))+1;

			if(most_sigf>neglog10_pval_cutoff){

				points(2:k, neglog10_pval_mat[i,2:k], type="l", 
					col=cols[i], lwd=2);
				points(2:k, neglog10_pval_mat[i,2:k], type="l", 
					col="black", lwd=.5);

				points(most_sigf_k, most_sigf, col=cols[i], cex=1, pch=16);
				points(most_sigf_k, most_sigf, col="black", cex=1.5, pch=1);

				# Label with variable names
				text(most_sigf_k, most_sigf, pred_names[i], font=2, cex=.75, 
					adj=c((most_sigf_k-2)/(k-2), -.5)
					);
			}
		}
		
		

	}


	plot_pval_over_k(arr_mspval_mat, title="Arrivers (All)");
	plot_pval_over_k(arr_mspval_mat, title="Arrivers (p-val<0.10)", pval_cutoff=0.1);
	plot_pval_over_k(arr_mspval_mat, title="Arrivers (p-val<0.05)", pval_cutoff=0.05);
	plot_pval_over_k(arr_mspval_mat, title="Arrivers (p-val<0.01)", pval_cutoff=0.01);

	plot_pval_over_k(dep_mspval_mat, title="Departers (All)");
	plot_pval_over_k(dep_mspval_mat, title="Departers (p-val<0.10)", pval_cutoff=0.1);
	plot_pval_over_k(dep_mspval_mat, title="Departers (p-val<0.05)", pval_cutoff=0.05);
	plot_pval_over_k(dep_mspval_mat, title="Departers (p-val<0.01)", pval_cutoff=0.01);

}

##############################################################################

cat("\nDone.\n")
dev.off();
warn=warnings();
if(length(warn)){
        print(warn);
}
q(status=0);
