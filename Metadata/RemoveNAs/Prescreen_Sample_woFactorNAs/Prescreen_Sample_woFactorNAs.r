#!/usr/bin/env Rscript

###############################################################################

library('getopt');

params=c(
        "summary_table_A", "s", 1, "character",
        "factor_file_A", "f", 1, "character",
	"inclusion_list_A", "l", 2, "character",
	"inclusion_list_B", "L", 2, "character",
        "summary_table_B", "S", 2, "character",
        "map_file", "m", 2, "character",
        "output_dir", "d", 1, "character",
        "target_var", "t", 2, "character",
        "required_var", "q", 2, "character",
	"num_iterations", "n", 2, "numeric"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

NUM_ITERATIONS=50000;

usage = paste(
        "\nUsage:\n", script_name, "\n",
	"\n",
	"	Input Parameters:\n",
	"	-s <data with A's IDs, i.e. summary_table>\n",
        "       -f <factor file with A's IDs>\n",
	"\n",
	"	[-l <list of A IDs to include>]\n",
	"	[-L <list of B IDs to include>]\n",
	"\n",
	"	Paired Screening parameters:\n",
	"	[-S <data with B's IDs, i.e. summary_table>\n",
	"	[-m <map file, mapping between A and B>\n",
	"\n",
	"	Output Parameters:\n",
	"	-d <output directory>\n",
	"\n",
	"	Factor constraints:\n",
	"	[-t <target variables file, default is all>]\n",
	"	[-q <required variables file, default is none>]\n",
	"\n",	
	"	Search Options:\n",
	"	[-n <number of search iterations, default=", NUM_ITERATIONS, ">\n",
        "\n",
	"This script will determine which samples can be used across\n",
	"the metadata and the sample files.  \n",
	"  The following criteria will be used:\n",
	"   1.) Samples not pairing up: -S and -m\n",
	"   2.) Samples not in the inclusion list: -l\n",
	"   3.) Samples with factors with NAs: -f (-t, -q)\n",
	"\n",
	"The output will be the following:\n",
	"	If 1 file specified:\n",
	"		<output directory>/<data A>.prescr.tsv\n",
	"		<output directory>/<factor file>.prescr.tsv\n",
	"\n",
	"	If 2 file specified:\n",
	"		<output directory>/<data A>.<A_name>.prescr.tsv\n",
	"		<output directory>/<factor file>.<A_name>.prescr.tsv\n",
	"		<output directory>/<data B>.<B_name>.prescr.tsv\n",
	"		<output directory>/<factor file>.<B_name>.prescr.tsv\n",
	"\n",
        "\n");

if(
	!length(opt$factor_file_A) || 
	!length(opt$summary_table_A) || 
	!length(opt$output_dir)
){
        cat(usage);
        q(status=-1);
}

FactorFile=opt$factor_file;
SummaryTableA=opt$summary_table_A;
OutputDir=opt$output_dir;

if(length(opt$summary_table_B)){
	SummaryTableB=opt$summary_table_B;
}else{
	SummaryTableB="";
}

if(length(opt$map_file)){
	MapFile=opt$map_file;
}else{
	MapFile="";
}

if(length(opt$target_var)){
	TargetVarFile=opt$target_var;	
}else{
	TargetVarFile="";
}

if(length(opt$inclusion_list_A)){
	InclusionListA=opt$inclusion_list_A;
}else{
	InclusionListA="";
}

if(length(opt$inclusion_list_B)){
	InclusionListB=opt$inclusion_list_B;
}else{
	InclusionListB="";
}

if(length(opt$required_var)){
	RequiredVarFile=opt$required_var;	
}else{
	RequiredVarFile="";
}

if(length(opt$num_iterations)){
	NumIterations=opt$num_iterations;	
}else{
	NumIterations=NUM_ITERATIONS;
}

cat("\n");
cat("Required:\n");
cat(" Summary Table (A): ", SummaryTableA, "\n");
cat(" Factor File: ", FactorFile, "\n");
cat(" Output Directory: ", OutputDir, "\n");
cat("\n");
cat("Optional:\n");
cat(" Inclusion List A: ", InclusionListA, "\n");
cat(" Inclusion List B: ", InclusionListB, "\n");
cat(" Summary Table B: ", SummaryTableB, "\n");
cat(" Map File: ", MapFile, "\n");
cat(" Targeted Variables: ", TargetVarFile, "\n");
cat(" Required Variables: ", RequiredVarFile, "\n");
cat(" Num Iterations: ", NumIterations, "\n");
cat("\n");

###############################################################################

load_factor_file=function(fn){
        inmat=read.delim(fn, sep="\t", header=TRUE, row.names=1, check.names=F, comment.char="", quote="");

        # Changes spaces to underscore
        var_names=colnames(inmat);
        var_names=gsub(" ", "_", var_names);
        colnames(inmat)=var_names;

        cat("  Num Factors: ", ncol(inmat), "\n", sep="");
        cat("  Num Samples: ", nrow(inmat), "\n", sep="");
        return(inmat);
}

write_factors=function(fname, table){

	cat("Writing: ", fname, " as factor file.\n", sep="");
        dimen=dim(factors);
        cat("Rows Exporting: ", dimen[1], "\n");
        cat("Cols Exporting: ", dimen[2], "\n");

        write.table(table, fname, quote=F, row.names=T, col.names=NA, sep="\t");

}

load_summary_table=function(st_fname){
        inmat=as.matrix(read.delim(st_fname, sep="\t", header=TRUE, row.names=1, check.names=FALSE, comment.char="", quote=""))

        num_categories=ncol(inmat)-1;
        num_samples=nrow(inmat);

        cat("Loaded Summary Table: ", st_fname, "\n", sep="");
        cat("  Num Categories: ", num_categories, "\n", sep="");
        cat("  Num Samples: ", num_samples, "\n", sep="");

        countsmat=inmat[,2:(num_categories+1)];

        return(countsmat);
}

write_summary_table=function(out_mat, fname){

	cat("Writing: ", fname, " as summary table.\n");
        fc=file(fname, "w");
        cat(file=fc, paste("sample_id\ttotal", paste(colnames(out_mat), collapse="\t"), sep="\t"));
        cat(file=fc, "\n");
        sample_names=rownames(out_mat);
        num_samples=nrow(out_mat);
        for(samp_idx in 1:num_samples){
                total=sum(out_mat[samp_idx,]);
                outline=paste(sample_names[samp_idx], total,
                        paste(out_mat[samp_idx,], collapse="\t"), sep="\t");
                cat(file=fc, outline);
                cat(file=fc, "\n");
        }
        close(fc);
}

load_list=function(list_fname){
        list=read.delim(list_fname, sep="\t", header=F, row.names=NULL, as.is=T, check.names=F, comment.char="#", quote="");
        return(list[,1]);
}

load_mapping_file=function(mp_fname, keep_a_ids, keep_b_ids){

        num_keep_a=length(keep_a_ids);
        num_keep_b=length(keep_b_ids);
        cat("Num A's IDs to keep: ", num_keep_a, "\n");
        cat("Num B's IDs to keep: ", num_keep_b, "\n");

        inmat=as.matrix(read.delim(mp_fname, sep="\t", header=TRUE, check.names=F, comment.char="", quote=""));

        # Keep Entry if record is in both lists
        keep_ix=c();
        orig_mat_rows=nrow(inmat);
        cat("Number of Mapping Entries Read: ", orig_mat_rows, "\n");
        for(i in 1:orig_mat_rows){
                if(any(inmat[i,1]==keep_a_ids) && any(inmat[i,2]==keep_b_ids)){
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

        return(map_info);
}

###############################################################################

cat("Loading Summary Table (A): ", SummaryTableA, "\n");
stA=load_summary_table(SummaryTableA);

cat("Loading Factor File: ", FactorFile, "\n");
factors=load_factor_file(FactorFile);

if(InclusionListA!=""){
	cat("Loading Inclusion List A: ", InclusionListA, "\n");
	incl_id_A=load_list(InclusionListA);
}else{
	cat("Inclusion list for A not specified.\n");
	incl_id_A=c();
}

if(InclusionListB!=""){
	cat("Loading Inclusion List B: ", InclusionListB, "\n");
	incl_id_B=load_list(InclusionListB);
}else{
	cat("Inclusion list for B not specified.\n");
	incl_id_B=c();
}

if(SummaryTableB!=""){
	cat("Loading Summary Table B: ", SummaryTableB, "\n");
	stB=load_summary_table(SummaryTableB);
}else{
	cat("Summary Table B not specified.\n");
	stB=c();
}

if(TargetVarFile!=""){
	cat("Loading Targeted Variables List: ", TargetVarFile, "\n");
	target_var=load_list(TargetVarFile);
}else{
	cat("No targeted variables specified.\n");
	target_var=c();
}

if(RequiredVarFile!=""){
	cat("Loading Required Variables List: ", RequiredVarFile, "\n");
	required_var=load_list(RequiredVarFile);
}else{
	cat("No required variables specified.\n");
}

###############################################################################

cat("\n");

stA_sampid=rownames(stA);
avail_id=stA_sampid;
cat("Example of Summary Table A IDs:\n");
print(head(stA_sampid));
cat("Available A IDs: ", length(avail_id), "\n");

if(InclusionListA!=""){
	avail_id=intersect(avail_id, incl_id_A);	
	cat("** Intersecting A's IDs with inclusion list.\n");
	cat("Available IDs: ", length(avail_id), "\n");
}

if(SummaryTableB!=""){
	stB_sampid=rownames(stB);
	cat("Example of Summary Table B IDs:\n");
	cat("Available B IDs: ", length(stB_sampid), "\n");
	print(head(stB_sampid));

	if(InclusionListB!=""){
		stB_sampid=intersect(stB_sampid, incl_id_B);	
		cat("Intersecting B's IDs with inclusion list.\n");
		cat("Available B IDs: ", length(stB_sampid), "\n");
	}
}

cat("\n");

if(SummaryTableB!=""){
	if(MapFile!=""){
		cat("Loading Mapping Info: ", MapFile, "\n");
		map_info=load_mapping_file(MapFile, stA_sampid, stB_sampid);
		avail_id=map_info[["a_id"]];
		cat("** Matching A/B IDs.\n");
	}else{
		cat("** Since no mapping file specified, intersecting summary tables A and B.\n");
		avail_id=intersect(avail_id, stB_sampid);
	}

	cat("Summary Tables A/B IDs:", length(avail_id));
}

cat("\n");
cat("\n");

###############################################################################

factors_IDs=rownames(factors);

cat("Example of factor IDs:\n");
print(head(factors_IDs));

cat("Intersecting Factor IDs with Available IDs.\n");
avail_id=intersect(factors_IDs, avail_id);
cat("Matching IDs: ", length(avail_id), "\n");

factors=factors[avail_id,,drop=F];

factor_var=colnames(factors);

cat("Checking targeted variables.\n");
if(length(target_var)){
	target_var_avail=intersect(factor_var, target_var);
}else{
	target_var_avail=factor_var;
}

if(length(target_var)>0){
	if(length(target_var_avail)!=length(target_var)){
		cat("WARNING: Targeted variables missing from available factors.\n");
		print(setdiff(target_var, target_var_avail));
	}
}
cat("\nTargeted Variables:\n");
print(target_var_avail);

cat("Subsetting factors to targeted varaibles.\n");
factors=factors[,target_var_avail,drop=F];


cat("Checking required variables.\n");
required_var_avail=intersect(target_var_avail, required_var);
if(length(required_var_avail)!=length(required_var)){
	cat("WARNING: Required variables missing from targeted factors.\n");
	print(setdiff(required_var, required_var_avail));
}
cat("\nRequired Variables:\n");
print(required_var_avail);

if(any(is.na(factors))){
        cat("NAs's found in factors...\n\n");

        script_path=paste(head(strsplit(script_name, "/")[[1]], -1), collapse="/");
	cat("Path of this script: ", script_path, "\n");
        source(paste(script_path, "/../Remove_NAs.r", sep=""));
        na_rm_res=remove_sample_or_factors_wNA_parallel(
		factors=factors, 
		required_variables=required_var_avail,
		num_trials=NumIterations, 
		num_cores=64
		);

	factors=na_rm_res$factors;

}else{
	cat("No NA's found in factors.  Exiting with no new output.\n");
}

###############################################################################

cat("Writing output to directory: ", OutputDir, "\n");

if(!file.exists(OutputDir)){
	dir.create(OutputDir);
}

if(!file.exists(OutputDir)){
	cat("Error, could not create or find directory: ", OutputDir, "\n");
	quit(-1);
}


if(MapFile!=""){
	aname_ext=paste(".", map_info[["a"]], sep="");
	bname_ext=paste(".", map_info[["b"]], sep="");
}else{
	aname_ext="";
	bname_ext="";
}
cat("Output File Extensions: ", aname_ext, " and ", bname_ext, "\n", sep="");

# Output Summary Table A
cat("Writing Presreened Summary Table A...\n");
stA_root=tail(strsplit(SummaryTableA, "/")[[1]],1);
stA_root=gsub("\\.summary_table\\.tsv$", "", stA_root);

recon_samp_ids=rownames(factors);
stA_recon=stA[recon_samp_ids,,drop=F];
write_summary_table(stA_recon, paste(OutputDir, "/", stA_root, aname_ext,".prescr.summary_table.tsv", sep=""));



# Output Factor File A
cat("Writing Factor File with A IDs as the primary key...\n");
factorf_root=tail(strsplit(FactorFile, "/")[[1]],1);
factorf_root=gsub("\\.tsv$", "", factorf_root);
write_factors(paste(OutputDir, "/", factorf_root, aname_ext, ".prescr.tsv", sep=""), factors);


if(SummaryTableB!=""){

	# Output Summary Table B
	cat("\n");
	cat("Writing Presreened Summary Table B...\n");
	stB_root=tail(strsplit(SummaryTableB, "/")[[1]],1);
	stB_root=gsub("\\.summary_table\\.tsv$", "", stB_root);

	if(MapFile!=""){
		cat("Remapping A's IDs to B's IDs before writing screened Summary Table B.\n");
		b_ids=map_info[["b_id"]]		
		names(b_ids)=map_info[["a_id"]];
		recon_samp_ids=b_ids[recon_samp_ids];
	}
	stB_recon=stB[recon_samp_ids,,drop=F];
	write_summary_table(stB_recon, paste(OutputDir, "/", stB_root, bname_ext, ".prescr.summary_table.tsv", sep=""));

	if(MapFile!=""){
		factor_a_ids=rownames(factors);
		factor_b_ids=b_ids[factor_a_ids];
		rownames(factors)=factor_b_ids;
	}	

	# Output Factor File B
	cat("\n");
	cat("Writing Factor File with B IDs as the primary key...\n");
	factorf_root=tail(strsplit(FactorFile, "/")[[1]],1);
	factorf_root=gsub("\\.tsv$", "", factorf_root);
	write_factors(paste(OutputDir, "/", factorf_root, bname_ext, ".prescr.tsv", sep=""), factors);
}




cat("done.\n");

