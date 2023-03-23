#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);
options(width=120);

params=c(
	"input", "i", 1, "character",
	"samp_id_colname", "s", 1, "character",
	"subj_id_colname", "S", 2, "character",
	"paired_map", "p", 1, "character",
	"a_colname", "A", 1, "character",
	"b_colname", "B", 1, "character",
	"percdiff_fname", "c", 2, "character",
	"log2diff_fname", "l", 2, "character",
	"nomndiff_fname", "n", 2, "character",
	"export_a", "a", 2, "logical",
	"export_b", "b", 2, "logical",
	"output", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input tab-separated column metadata/factors file>\n",
	"	-s <sample ID column name in metadata/factors file, default is first column>\n",
	"\n",
	"	-p <paired samples map>\n",
	"	-A <A Mapping ColumnName>\n",
	"	-B <B Mapping ColumnName>\n",
	"\n",
	"	[-c <Percent Difference, Target Variables Filename>]\n",
	"	[-l <Log2 Ratio, Target Variables Filename>]\n",
	"	[-n <Nominal Difference, Target Variables Filename>]\n",
	"\n",
	"	[-a (Export A values)]\n",
	"	[-b (Export B values)]\n",
	"	[-S <subject ID column name, for export>]\n",
	"	-o <output tab_separated column data file>\n",
	"\n",
	"This script will read in the metadata and paired map and compute the\n",
	"paired differences between the samples.\n",
	"\n",
	"You must choose 1 of the -c, -l, or -n options.\n",
	"\n",
	"\n");

if(!length(opt$input) || 
	!length(opt$samp_id_colname) || 
	!length(opt$paired_map) || 
	!length(opt$a_colname) || 
	!length(opt$b_colname) || 
	!length(opt$output)){
	cat(usage);
	q(status=-1);
}

InputFName=opt$input;
SampleIDColname=opt$samp_id_colname;

PairingsFile=opt$paired_map;
AColname=opt$a_colname;
BColname=opt$b_colname;


PercDiffFname="";
Log2DiffFname="";
NomnDiffFname="";

if(length(opt$percdiff_fname)){
	PercDiffFname=opt$percdiff_fname;
}
if(length(opt$log2diff_fname)){
	Log2DiffFname=opt$log2diff_fname;
}
if(length(opt$nomndiff_fname)){
	NomnDiffFname=opt$nomndiff_fname;
}


ExportA=F;
ExportB=F;

if(length(opt$export_a)){
	ExportA=T;
}
if(length(opt$export_b)){
	ExportB=T;
}

SubjectIDColname="";
if(length(opt$subj_id_colname)){
	SubjectIDColname=opt$subj_id_colname;
}

OutputFName=opt$output;

cat("Input Filename: ", InputFName, "\n");
cat("Sample ID Colname: ", SampleIDColname, "\n");
cat("Subject ID Colname: ", SubjectIDColname, "\n");
cat("Pairings Map: ", PairingsFile, "\n");
cat("A Colname: ", AColname, "\n");
cat("B Colname: ", BColname, "\n");
cat("\n");
cat("Perc Diff Target Filename: ", PercDiffFname, "\n");
cat("Log2 Ratio Target Filename: ", Log2DiffFname, "\n");
cat("Nominal Diff Filename: ", NomnDiffFname, "\n");
cat("\n");
cat("Export A: ", ExportA, "\n");
cat("Export B: ", ExportA, "\n");
cat("Output Filename: ", OutputFName, "\n");

##############################################################################

load_factors=function(fname){
	factors=data.frame(read.table(fname,  header=TRUE, check.names=FALSE, 
		as.is=T, comment.char="", quote="", sep="\t"));

	#print(factors);

	dimen=dim(factors);
	cat("Rows Loaded: ", dimen[1], "\n");
	cat("Cols Loaded: ", dimen[2], "\n");

	return(factors);
}

write_factors=function(fname, table){

	dimen=dim(factors);
	cat("Rows Exporting: ", dimen[1], "\n");
	cat("Cols Exporting: ", dimen[2], "\n");
	
	write.table(table, fname, quote=F, row.names=F, sep="\t");

}

load_list=function(fname){
	cat("Reading in: ", fname, "\n");
	arr=read.delim(fname, header=F, sep="\t", stringsAsFactors=F)[,1];
	cat("Num variables loaded: ", length(arr), "\n");
	return(arr);
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

##############################################################################
# Additional Functions
pick_A=function(val){
	return(val$a);
}

pick_B=function(val){
	return(val$b);
}

##############################################################################

filenames=list();

if(NomnDiffFname!=""){
	filenames[["nomn"]]=NomnDiffFname;
}
if(PercDiffFname!=""){
	filenames[["perc"]]=PercDiffFname;
}
if(Log2DiffFname!=""){
	filenames[["log2"]]=Log2DiffFname;
}

target_types=names(filenames);

# Load variable list and how to calculate differences
target_variables_list=list();
all_diffed_variables=character();
for(type in target_types){
	target_variables_list[[type]]=sort(load_list(filenames[[type]]));
	all_diffed_variables=c(all_diffed_variables, target_variables_list[[type]]);
}
all_diffed_variables=sort(all_diffed_variables);

cat("\n");
cat("Variables for each difference calculation type:\n");
print(target_variables_list);

##############################################################################

# Load factors
cat("Loading Factors/Metadata File: ", InputFName, "\n");

factors=load_factors(InputFName);
if(SampleIDColname!=""){
	rownames(factors)=factors[,SampleIDColname];
}else{
	rownames(factors)=factors[,1];
}
factor_rownames=rownames(factors);

num_samples=nrow(factors);

# Load Pairings
cat("\n");
cat("Loading Pairings File: ", PairingsFile, "\n");
all_pairings_map=load_mapping(PairingsFile, AColname, BColname);
num_pairs=nrow(all_pairings_map);
print(all_pairings_map);

##############################################################################

calc_diff=function(a_mat, b_mat, diff_type){

	num_samp=nrow(a_mat);
	num_var=ncol(a_mat);

	cat("Num Rows / Cols : ", num_samp, " / ", num_var, "\n");
	
	out_mat=matrix(NA, nrow=num_samp, ncol=num_var);
	rownames(out_mat)=rownames(a_mat);
	colnames(out_mat)=colnames(a_mat);

	for(rix in 1:num_samp){
		arow=a_mat[rix,];
		brow=b_mat[rix,];

		if(diff_type == "nomn"){
			res=brow-arow;
		}else if (diff_type == "perc"){
			res=100*(brow-arow)/arow;
		}else if (diff_type == "log2"){
			res=log2(brow/arow);
		}else{
			cat("Error, invalid difference type: ", diff_type, "\n");
			quit(status=-1);
		}

		out_mat[rix, ]=unlist(res[1,]);

	}

	return(out_mat);
}

output_matrix=numeric();

if(SubjectIDColname!=""){
	output_matrix=as.matrix(factors[all_pairings_map[, AColname, drop=F], SubjectIDColname], ncol=1);
	colnames(output_matrix)=SubjectIDColname;
}

for(type in target_types){

	A_matrix=factors[all_pairings_map[, AColname, drop=F], target_variables_list[[type]]];
	B_matrix=factors[all_pairings_map[, BColname, drop=F], target_variables_list[[type]]];

	diff_matrix=calc_diff(A_matrix, B_matrix, type);

	varnames=colnames(diff_matrix);
	prefixed_varnames=paste(type, "_diff_", varnames, sep="");
	colnames(diff_matrix)=prefixed_varnames;

	output_matrix=cbind(output_matrix, diff_matrix);
}

if(ExportA){
	Avals=factors[all_pairings_map[, AColname], all_diffed_variables];
	varnames=colnames(Avals);
	prefixed_varnames=paste(AColname, "_", varnames, sep="");
	colnames(Avals)=prefixed_varnames;
	output_matrix=cbind(output_matrix, Avals);
}

if(ExportB){
	Bvals=factors[all_pairings_map[, BColname], all_diffed_variables];
	varnames=colnames(Bvals);
	prefixed_varnames=paste(BColname, "_", varnames, sep="");
	colnames(Bvals)=prefixed_varnames;
	output_matrix=cbind(output_matrix, Bvals);
}

##############################################################################

write_factors(paste(OutputFName, ".tsv", sep=""), output_matrix);

##############################################################################

cat("\nDone.\n");

print(warnings());
q(status=0);
