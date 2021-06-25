#!/usr/bin/env Rscript

###############################################################################

cat("\n\n");

library(MASS)
library('getopt');
library('vegan');
library('plotrix');


DEF_DISTTYPE="man";
MAX_CLUSTER_CUTS=8;

params=c(
	"input_summary_table", "s", 1, "character",
	"mapping_file", "m", 1, "character",
	"output_filename_root", "o", 1, "character",
	"a_name", "A", 1, "character",
	"b_name", "B", 1, "character",

	"factor_file", "f", 2, "character",
	"target_variables", "c", 2, "character",
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

FactorFile=ifelse(length(opt$factor_file), opt$factor_file, "");
TargetVariablesFile=ifelse(length(opt$target_variables), opt$target_variables, "");

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

	mapping=as.list(x=inmat[,1]);
	names(mapping)=inmat[,2];

	coln=colnames(inmat);
	map_info=list();
	map_info[["map"]]=mapping;
	map_info[["a"]]=coln[1];
	map_info[["b"]]=coln[2];
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
	data=read.delim(fn, header=F, sep="\t", quote="", as.is=T);
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


##############################################################################

###############################################################################
# Load summary table, factor file and pairing file, and reconcile

output_fname_root = paste(OutputFileRoot, ".", DistType, sep="");

cat("\n");
cat("Loading summary table:", InputSumTab, "\n");
counts_mat=load_summary_table(InputSumTab);

sumtab_sample_ids=rownames(counts_mat);

if(FactorFile!=""){
	factors_matrix=load_factors(FactorFile);
}else{
	factors_matrix=NA;
}

cat("Loading Mapping file:", MappingFile, "\n");
map_info=load_mapping_file(MappingFile, sumtab_sample_ids, Aname, Bname);
print(map_info);

sample_ids_pairable=c(map_info[["a_id"]], map_info[["b_id"]]);

cat("Removing samples without complete mappings...\n");
counts_mat=counts_mat[sample_ids_pairable,];

num_samples=nrow(counts_mat);
cat("Num usable samples: ", num_samples, "\n");

###############################################################################

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
});
plot_text(param_summary);

###############################################################################

###############################################################################

plot_dendro_colored_by_pairtype=function(
	clus_mem, inden, aname, bname, a_ids, b_ids, k, h, grp_sizes, grp_col,
	option_name, option_args
){

	num_pairs=length(a_ids);
	num_samples=num_pairs*2;
	sample_to_color_map=clus_mem;
	names(sample_to_color_map)=c(a_ids, b_ids);
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

	# Top-Left
	par(mar=c(0,0,0,0));
	plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n",
		 main="", xlab="", ylab="", xaxt="n", yaxt="n", bty="n");
	text(.5, .5, paste("k = ", k, sep=""), cex=2, font=2);

	# Top Margin
	par(mar=c(samp_labl_spc, 0, title_spc,0));
	plot(dendra_mod, main="", horiz=F, xaxt="n", yaxt="n", xlab="", ylab="", xlim=dendro_plot_range);
	mtext(aname, side=3, font=2, line=.5);
	abline(h=h, col="blue", lty=2, lwd=.7);
	#abline(v=c(1,num_samples));

	# Right Margin
	par(mar=c(0, title_spc,0, samp_labl_spc));
	plot(dendrb_mod, main="", horiz=T, xaxt="n", yaxt="n", xlab="", ylab="", ylim=dendro_plot_range);
	mtext(bname, side=2, font=2, line=.5)
	abline(v=h, col="blue", lty=2, lwd=.7);
	#abline(h=c(1,num_samples));

	# Center Field
	par(mar=c(0,0,0,0));
	plot(0,0, xlim=dendro_plot_range, ylim=dendro_plot_range, 
		type="n", main="", xlab="", ylab="", xaxt="n", yaxt="n", bty="n");
	#abline(h=dendro_plot_range);
	#abline(v=dendro_plot_range);
	abline(h=c(0, group_lines)+.5, lwd=.5);
	abline(v=c(0, group_lines)+.5, lwd=.5);

	topmost=par()$usr[4];
	leftmost=par()$usr[1];

	text(
		group_centers, rep(topmost, k), 1:k, cex=2, font=2, pos=1
	);
	text(
		rep(leftmost, k), group_centers, 1:k, cex=2, font=2, pos=4
	);

	if(option_name=="drawMatchingDots"){
		
		
		points(option_args[,1], option_args[,2], pch=22, col=option_args[,3]);

	}

}


###############################################################################

analyze_cut=function(k, num_pairs, a_clus_mem, b_clus_mem, a_ids, b_ids){
	#print(a_clus_mem);
	#print(b_clus_mem);
	#print(a_ids);
	#print(b_ids);

	ctng_counts=matrix(0, nrow=k, ncol=k);
	rownames(ctng_counts)=paste("a_cl_", 1:k, sep="");
	colnames(ctng_counts)=paste("b_cl_", 1:k, sep="");

	for(i in 1:num_pairs){
		ctng_counts[
			a_clus_mem[a_ids[i]],
			b_clus_mem[b_ids[i]]
		]=
			ctng_counts[
				a_clus_mem[a_ids[i]],
				b_clus_mem[b_ids[i]]
			]+1;
	}

	print(ctng_counts);

}

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

grp_a_col="blue";
grp_b_col="green";

lf_names=get_clstrd_leaf_names(orig_dendro);
print(lf_names);

num_pairs=map_info[["num_pairs"]];
pos_map=matrix(0, nrow=num_pairs, ncol=3);
a_id=map_info[["a_id"]];
b_id=map_info[["b_id"]];
colnames(pos_map)=c(map_info[["a"]], map_info[["b"]], "same_clus");


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
		pos_map[pix,]=c(a_pos, b_pos, same_clus);
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

}


quit();

plot_dendro_contigency=function(hclA, hclB, acuts, bcuts, namea, nameb, idsb){

	color_denfun_bySample=function(n){
		if(is.leaf(n)){
			leaf_attr=attributes(n);
			leaf_name=leaf_attr$label;
			ind_color=sample_to_color_map[leaf_name];
			if(is.null(ind_color)){
				ind_color="black";
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

	# Compute
	cat("Working on: ", nameb, ": ", bcuts, " x ", namea, ": ", acuts, "\n", sep="");

	dendra=as.dendrogram(hclA);
	dendrb=as.dendrogram(hclB);

	memb_byA=cutree(hclA, k=acuts);
	memb_byB=cutree(hclB, k=bcuts);

	dendr_names_a=get_clstrd_leaf_names(dendra);
	dendr_names_b=get_clstrd_leaf_names(dendrb);

	memb_byA=reorder_member_ids(memb_byA, dendr_names_a);
	memb_byB=reorder_member_ids(memb_byB, dendr_names_b);

	dend_mids_a=get_middle_of_groups(dendr_names_a, memb_byA);
	dend_mids_b=get_middle_of_groups(dendr_names_b, memb_byB);

	num_members=length(memb_byA);

	grp_cnts_a=(table(memb_byA)[1:acuts]);
	grp_cnts_b=(table(memb_byB)[1:bcuts]);;
	
	grp_prop_a=grp_cnts_a/num_members;
	grp_prop_b=grp_cnts_b/num_members;;

	cat("Group Proportions A:\n");
	print(grp_prop_a);

	cat("Group Proportions B:\n");
	print(grp_prop_b);

	# Count up observed
	ab_cnts_obs_mat=matrix(0, nrow=bcuts, ncol=acuts);
	for(i in 1:num_members){
		ma=memb_byA[i];
		mb=memb_byB[i];
		ab_cnts_obs_mat[mb, ma]=ab_cnts_obs_mat[mb, ma]+1;
	}

	ab_prop_obs_mat=ab_cnts_obs_mat/num_members;
	
	# Calculate expected
	ab_prop_exp_mat=grp_prop_b %*% t(grp_prop_a);
	ab_cnts_exp_mat=ab_prop_exp_mat * num_members;

	cat("Observed Counts\n");
	print(ab_cnts_obs_mat);

	cat("Observed Proportions\n");
	print(ab_prop_obs_mat);

	cat("Expected Counts:\n");
	print(ab_cnts_exp_mat);

	cat("Expected Proportions:\n");
	print(ab_prop_exp_mat);
	
	# Compute pvalue for contingency table
	cst=chisq.test(ab_cnts_obs_mat);
	ct_cst_pval=cst$p.value;

	#print(cst);
	#print(cst$observed);
	#print(cst$expected);

	fish_exact_mat=matrix(0, nrow=bcuts, ncol=acuts);

	for(rix in 1:bcuts){
		for(cix in 1:acuts){

			twobytwo=matrix(0, nrow=2, ncol=2);
			twobytwo[1,1]=ab_cnts_obs_mat[rix, cix];
			twobytwo[1,2]=sum(ab_cnts_obs_mat[-rix, cix]);
			twobytwo[2,1]=sum(ab_cnts_obs_mat[rix, -cix]);
			twobytwo[2,2]=sum(ab_cnts_obs_mat[-rix, -cix]);
			ind_cst_res=fisher.test(twobytwo);
			fish_exact_mat[rix, cix]=ind_cst_res$p.value;
		}
	}
	#print(fish_exact_mat);

	# Calculate conditional probabilities
	cnd_pr_agb=matrix(0, nrow=bcuts, ncol=acuts);
	cnd_pr_bga=matrix(0, nrow=bcuts, ncol=acuts);

	for(rix in 1:bcuts){
		cnd_pr_agb[rix,]=ab_prop_obs_mat[rix,]/grp_prop_b[rix];
	}

	for(cix in 1:acuts){
		cnd_pr_bga[,cix]=ab_prop_obs_mat[,cix]/grp_prop_a[cix];
	}

	# Calculate cumulative probability for highest probability mapping
	cum_pr_bga=sum(apply(ab_prop_obs_mat, 2, max));
	cum_pr_agb=sum(apply(ab_prop_obs_mat, 1, max));
	

	##########################################
	# Plot

	orig_par=par(no.readonly=T);

	par(oma=c(1,1,1,1));

	table_sp=5;
	layout_mat=matrix(c(
		1,rep(2, table_sp),
		rep(c(3,rep(4, table_sp)), table_sp)),
		nrow=table_sp+1, byrow=T);
	#print(layout_mat);
	layout(layout_mat);

	# plot top/left spacer
	par(mar=c(0,0,0,0));
	plot(0,0,type="n", bty="n", xlab="", ylab="", main="", xaxt="n", yaxt="n");
	text(0,0, paste
		(nameb, ": ", bcuts, "\n x \n", namea, ": ", acuts, 
		"\n\nX^2 Test p-value:\n", sprintf("%1.3g", ct_cst_pval),
		"\n\nCumulative Top:\nPr(", nameb, "|", namea, ")=\n", round(cum_pr_bga, 3), 
		"\nPr(", namea, "|", nameb, ")=\n", round(cum_pr_agb, 3),
		 sep=""), cex=.8, font=2);



	# Scale leaf sample IDs
	label_scale=.2;
	dendra=dendrapply(dendra, text_scale_denfun);
	dendrb=dendrapply(dendrb, text_scale_denfun);
	
	# Color both dendrograms by A clustering
	sample_to_color_map=memb_byA;
	dendra=dendrapply(dendra, color_denfun_bySample);
	names(sample_to_color_map)=idsb;
	dendrb=dendrapply(dendrb, color_denfun_bySample);
	
	# Find height where clusters separate
	acutheight=find_height_at_k(hclA, acuts);
	bcutheight=find_height_at_k(hclB, bcuts);

	top_label_spc=4;
	left_label_spc=4;
	title_spc=2;

	# Plot A Dendrogram
	par(mar=c(5,left_label_spc,title_spc,0));
	plot(dendra, main=namea, horiz=F, yaxt="n", xaxt="n", xlab="", ylab="", xlim=c(-1,num_members+1));
	abline(h=acutheight, col="blue", lty=2, lwd=.7);
	abline(v=c(0,cumsum(grp_cnts_a)+.5), col="grey75", lwd=.5);
	trans_dend_mids_a=remap_coord(dend_mids_a, 0, num_members, 0, 1);

	# Plot B Dendrogram
	par(mar=c(0,title_spc,top_label_spc,5));
	plot(dendrb, main=nameb, horiz=T, xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(-1,num_members+1));
	abline(v=bcutheight, col="blue", lty=2, lwd=.7);
	abline(h=c(0, cumsum(grp_cnts_b)+.5), col="grey75", lwd=.5);
	trans_dend_mids_b=remap_coord(dend_mids_b, 0, num_members, 0, 1);

	# Plot shared statistics
	par(mar=c(0,left_label_spc,top_label_spc,0));
	plot(0,0, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n");
	#points(c(0,0,1,1), c(0,1,0,1));
	axis(3, at=trans_dend_mids_a, 1:acuts, tick=F, line=NA, font=2, cex.axis=2);
	axis(3, at=trans_dend_mids_a, grp_cnts_a, tick=F, line=-1, font=2, cex.axis=1);
	axis(2, at=trans_dend_mids_b, 1:bcuts, tick=F, line=NA, font=2, cex.axis=2);
	axis(2, at=trans_dend_mids_b, grp_cnts_b, tick=F, line=-1, font=2, cex.axis=1);

	cellab_size=min(1, 4/sqrt(acuts^2+bcuts^2));

	minpval=min(fish_exact_mat);
	cell_bounds_x=c(0,cumsum(grp_cnts_a)+.5)/num_members;
	cell_bounds_y=c(0,cumsum(grp_cnts_b)+.5)/num_members;

	for(colx in 1:acuts){
		for(rowx in 1:bcuts){

			cur_pval=fish_exact_mat[rowx, colx];

			cell_info=paste(
				"ob ct: ", ab_cnts_obs_mat[rowx, colx], "\n",
				"ob pr: ", round(ab_prop_obs_mat[rowx, colx], 3), "\n",
				"ex ct: ", round(ab_cnts_exp_mat[rowx, colx], 1), "\n",
				"ex pr: ", round(ab_prop_exp_mat[rowx, colx], 3), "\n",
				"fe pv: ", sprintf("%3.3g", cur_pval), "\n",
				"\n",
				"Pr(", namea, "|", nameb, ")=\n", sprintf("%3.3g", cnd_pr_agb[rowx, colx]), "\n",
				"Pr(", nameb, "|", namea, ")=\n", sprintf("%3.3g", cnd_pr_bga[rowx, colx]), "\n",
				sep="");

			
			# Highlight significant cells labels
			col="grey";
			font=1;
			if(cur_pval<.05){
				col="blue";
			}
			if(cur_pval<=minpval){
				font=2;
				points(c(cell_bounds_x[colx], cell_bounds_x[colx], cell_bounds_x[colx+1], 
					cell_bounds_x[colx+1], cell_bounds_x[colx]),
					c(cell_bounds_y[rowx], cell_bounds_y[rowx+1], cell_bounds_y[rowx+1], 
					cell_bounds_y[rowx], cell_bounds_y[rowx]), 
					col="cornflowerblue", type="l");
			}
				
			text(trans_dend_mids_a[colx], trans_dend_mids_b[rowx], cell_info, cex=cellab_size,
				font=font, col=col);
		}
	}

	par(orig_par);

	result=list();
	result[["chisqr_test_pval"]]=ct_cst_pval;
	result[["cumulative_cond_prob_agb"]]=cum_pr_agb;
	result[["cumulative_cond_prob_bga"]]=cum_pr_bga;

	return(result);

}

##############################################################################

plot_dendro_group_compare=function(hclA, hclB, acuts, bcuts, namea, nameb, idsb){

	color_denfun_bySample=function(n){
		if(is.leaf(n)){
			leaf_attr=attributes(n);
			leaf_name=leaf_attr$label;
			ind_color=sample_to_color_map[leaf_name];
			if(is.null(ind_color)){
				ind_color="black";
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

	# Compute
	cat("Working on: ", nameb, ": ", bcuts, " x ", namea, ": ", acuts, "\n", sep="");

	dendra=as.dendrogram(hclA);
	dendrb=as.dendrogram(hclB);

	memb_byA=cutree(hclA, k=acuts);
	memb_byB=cutree(hclB, k=bcuts);

	dendr_names_a=get_clstrd_leaf_names(dendra);
	dendr_names_b=get_clstrd_leaf_names(dendrb);

	memb_byA=reorder_member_ids(memb_byA, dendr_names_a);
	memb_byB=reorder_member_ids(memb_byB, dendr_names_b);

	dend_mids_a=get_middle_of_groups(dendr_names_a, memb_byA);
	dend_mids_b=get_middle_of_groups(dendr_names_b, memb_byB);

	num_members=length(memb_byA);

	grp_cnts_a=(table(memb_byA)[1:acuts]);
	grp_cnts_b=(table(memb_byB)[1:bcuts]);;
	
	grp_prop_a=grp_cnts_a/num_members;
	grp_prop_b=grp_cnts_b/num_members;;

	cat("Group Proportions A:\n");
	print(grp_prop_a);

	cat("Group Proportions B:\n");
	print(grp_prop_b);

	##########################################

	orig_par=par(no.readonly=T);
	#par(mfrow=c(4,1));
	plsz=3;
	layout_mat=matrix(c(
		rep(1, plsz),
		rep(2, plsz),
		3,
		rep(4, plsz),
		rep(5, plsz)
		), byrow=T, ncol=1);
	layout(layout_mat);

	# Find height where clusters separate
	acutheight=find_height_at_k(hclA, acuts);
	bcutheight=find_height_at_k(hclB, bcuts);

	top_label_spc=4;
	left_label_spc=1;
	title_spc=2;

	#-----------------------------------------------------------------------------

	# Scale leaf sample IDs
	label_scale=.2;
	dendraA=dendrapply(dendra, text_scale_denfun);
	dendrbA=dendrapply(dendrb, text_scale_denfun);
	
	# Color both dendrograms by A clustering
	sample_to_color_map=memb_byA;
	dendraA=dendrapply(dendraA, color_denfun_bySample);
	idsa=names(sample_to_color_map);
	names(sample_to_color_map)=idsb;
	dendrbA=dendrapply(dendrbA, color_denfun_bySample);

	# Plot A Dendrogram
	par(mar=c(5,left_label_spc,title_spc,0));
	plot(dendraA, main=paste(namea, ": ", acuts, " cuts", sep=""), 
		horiz=F, yaxt="n", xaxt="n", xlab="", ylab="", xlim=c(-1,num_members+1));
	abline(h=acutheight, col="blue", lty=2, lwd=.7);
	abline(v=c(0,cumsum(grp_cnts_a)+.5), col="grey75", lwd=.5);

	# Plot B Dendrogram
	par(mar=c(5,left_label_spc,title_spc,0));
	plot(dendrbA, main=paste(nameb, ": colored by ", namea, sep=""), 
		horiz=F, xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(-1,num_members+1));
	abline(h=bcutheight, col="blue", lty=2, lwd=.7);
	abline(v=c(0, cumsum(grp_cnts_b)+.5), col="grey75", lwd=.5);

	#-----------------------------------------------------------------------------
	par(mar=c(0,0,0,0));
	plot(0,0, type="n", main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n");
	abline(h=0, col="grey", lwd=5);

	#-----------------------------------------------------------------------------

	# Color both dendrograms by B clustering
	label_scale=.2;
	dendraB=dendrapply(dendra, text_scale_denfun);
	dendrbB=dendrapply(dendrb, text_scale_denfun);

	sample_to_color_map=memb_byB;
	dendrbB=dendrapply(dendrbB, color_denfun_bySample);
	names(sample_to_color_map)=idsa;
	dendraB=dendrapply(dendraB, color_denfun_bySample);

	# Plot B Dendrogram
	par(mar=c(5,left_label_spc,title_spc,0));
	plot(dendrbB, main=paste(nameb, ": ", bcuts, " cuts", sep=""),
		horiz=F, xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(-1,num_members+1));
	abline(h=bcutheight, col="blue", lty=2, lwd=.7);
	abline(v=c(0, cumsum(grp_cnts_b)+.5), col="grey75", lwd=.5);

	# Plot A Dendrogram
	par(mar=c(5,left_label_spc,title_spc,0));
	plot(dendraB, main=paste(namea, ": colored by ", nameb, sep=""),
		horiz=F, yaxt="n", xaxt="n", xlab="", ylab="", xlim=c(-1,num_members+1));
	abline(h=acutheight, col="blue", lty=2, lwd=.7);
	abline(v=c(0,cumsum(grp_cnts_a)+.5), col="grey75", lwd=.5);

	par(orig_par);


}

##############################################################################

classic_mds_pts_A=matrix(0, nrow=num_samples, ncol=2); 
nonparm_mds_pts_A=matrix(0, nrow=num_samples, ncol=2); 
classic_mds_pts_B=matrix(0, nrow=num_samples, ncol=2); 
nonparm_mds_pts_B=matrix(0, nrow=num_samples, ncol=2); 

class_mdsA_res=cmdscale(dist_mat_A, k=2);
class_mdsB_res=cmdscale(dist_mat_B, k=2);
classic_mds_pts_A=class_mdsA_res;
classic_mds_pts_B=class_mdsB_res;

nonpr_mdsA_res=isoMDS(dist_mat_A, k=2);
nonpr_mdsB_res=isoMDS(dist_mat_B, k=2);
nonparm_mds_pts_A=nonpr_mdsA_res$points;
nonparm_mds_pts_B=nonpr_mdsB_res$points;

hcl_A=hclust(dist_mat_A, method="ward.D2");
hcl_B=hclust(dist_mat_B, method="ward.D2");

clus4_A=cutree(hcl_A, k=max_cuts);
clus4_B=cutree(hcl_B, k=max_cuts);

cat("Comparing MDS plots:\n");
par(mfrow=c(2,2));
compare_mds(classic_mds_pts_A, classic_mds_pts_B, "Classical MDS", 
	clus4_A, clus4_B, map_info[["a"]], map_info[["b"]]);
compare_mds(nonparm_mds_pts_A, nonparm_mds_pts_B, "NonMetric MDS", 
	clus4_A, clus4_B, map_info[["a"]], map_info[["b"]]);

##############################################################################

analyze_dendro_cont=T;

if(analyze_dendro_cont){

	probmat=matrix(0, nrow=max_cuts, ncol=max_cuts);
	rownames(probmat)=c(paste(map_info[["b"]], 1:max_cuts));
	colnames(probmat)=c(paste(map_info[["a"]], 1:max_cuts));

	pval_mat=probmat;
	cum_agb_mat=probmat;
	cum_bga_mat=probmat;

	for(acuts in 2:max_cuts){
		for(bcuts in 2:max_cuts){

			cont_res=plot_dendro_contigency(hcl_A, hcl_B, acuts, bcuts, 
				map_info[["a"]], map_info[["b"]], map_info[["b_id"]]);

			plot_dendro_group_compare(hcl_A, hcl_B, acuts, bcuts,
                                map_info[["a"]], map_info[["b"]], map_info[["b_id"]]);

			pval_mat[bcuts, acuts]=cont_res[["chisqr_test_pval"]];
			cum_agb_mat[bcuts, acuts]=cont_res[["cumulative_cond_prob_agb"]];
			cum_bga_mat[bcuts, acuts]=cont_res[["cumulative_cond_prob_bga"]];
		}
	}

	pval_mat=pval_mat[2:max_cuts, 2:max_cuts, drop=F];
	cum_agb_mat=cum_agb_mat[2:max_cuts, 2:max_cuts, drop=F];
	cum_bga_mat=cum_bga_mat[2:max_cuts, 2:max_cuts, drop=F];

	logodds=log2(cum_agb_mat/cum_bga_mat);


	print(pval_mat);
	min_cont_pval=min(pval_mat);
	min_idx=which(pval_mat==min_cont_pval, arr.ind=T);

	anames=colnames(pval_mat);
	bnames=rownames(pval_mat);

	paint_matrix(-log10(pval_mat), title="Contigency Table Dimension Log10(P-Values)");

	par(mfrow=c(1,1));
	plot_text(c(
		"Contingency Table Chi-Squared Tests by Num Clusters, p-value:",
		"",
		capture.output(print(signif(pval_mat, 3))),
		"",
		"",
		"",
		"Contingency Table Chi-Squared Tests by Num Clusters, -log10(p-value):",
		"",
		capture.output(print(-log10(pval_mat))),
		"",
		"",
		paste("Min P-Value: ", sprintf("%3.3g", min_cont_pval), 
			" at (", bnames[min_idx[1]], ", ", anames[min_idx[2]], ")", sep="")
	));

	plot_text(c(
		paste("Given a sample from ", map_info[["b"]], 
			" what's the probability classifying it in ",  map_info[["a"]], "?", sep=""),
		paste("Cumulative Pr(", map_info[["a"]], "|", map_info[["b"]], "):", sep=""),
		"",
		capture.output(print(round(cum_agb_mat,3))),
		"",
		"",
		paste("Given a sample from ", map_info[["a"]], 
			" what's the probability classifying it in ",  map_info[["b"]], "?", sep=""),
		paste("Cumulative Pr(", map_info[["b"]], "|", map_info[["a"]], "):", sep=""),
		"",
		capture.output(print(round(cum_bga_mat,3))),
		"",
		"",
		paste("Log(Pr(", map_info[["a"]], "|", map_info[["b"]], ")/Pr(", 
			map_info[["b"]], "|", map_info[["a"]], ")):", sep=""),
		paste("  Positive Log Prob Ratio implies ", map_info[["b"]], 
			" predicts ", map_info[["a"]], " better than vice versa.", sep=""),
		"",
		capture.output(print(round(logodds, 2)))
	));

	plot_dendro_contigency(hcl_A, hcl_B, min_idx[2]+1, min_idx[1]+1, 
		map_info[["a"]], map_info[["b"]], map_info[["b_id"]]);

	plot_dendro_group_compare(hcl_A, hcl_B, min_idx[2]+1, min_idx[1]+1,
                map_info[["a"]], map_info[["b"]], map_info[["b_id"]]);

	# Write optimal cuts to file
	cnt_fh=file(paste(OutputFileRoot, ".", map_info[["a"]], ".cuts", sep=""), "w");
	cat(file=cnt_fh, min_idx[2]+1, "\n", sep="");
	close(cnt_fh);

	cnt_fh=file(paste(OutputFileRoot, ".", map_info[["b"]], ".cuts", sep=""), "w");
	cat(file=cnt_fh, min_idx[1]+1, "\n", sep="");
	close(cnt_fh);
}


##############################################################################

cat("\nDone.\n")
dev.off();
warn=warnings();
if(length(warn)){
        print(warn);
}
q(status=0);
