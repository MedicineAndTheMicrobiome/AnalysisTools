#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);
options(width=120);

params=c(
	"input", "i", 1, "character",
	"samp_id_colname", "S", 2, "character",
	"paired_map", "p", 1, "character",
	"a_colname", "A", 1, "character",
	"b_colname", "B", 1, "character",
	"formulas", "f", 1, "character",
	"output", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input tab-separated column data file>\n",
	"	[-S <sample ID column name, default is first column>]\n",
	"	-p <paired samples map>\n",
	"	-A <A Mapping ColumnName>\n",
	"	-B <B Mapping ColumnName>\n",
	"	-f <formulas file>\n",
	"	-o <output tab_separated column data file>\n",
	"\n",
	"This script will read in a metadata file and paired mapping file\n",
	"and run through the list of formulas to add columns of\n",
	"derived data to the table.\n",
	"\n",
	"Example formulas file entries:\n",
	"	nomn_score=nomn_diff(score);\n",
	"	pctd_counts=perc_diff(counts);\n",
	"	lrat_percentage=log2_ratio(percentage);\n",
	"\n",
	"Note: The formula you want will depend on what linearizes the differences.\n",
	"The nominal difference is just the difference between the two values:\n",
	"	For example between day 50 and day 30, would be 20 days.\n",
	"	The nominal difference is useful if 0 has no meaning.\n",
	"\n",
	"The percent difference is 100*(B-A)/A, which is useful if A and B\n",
	"	have positive values and 0 is absolute.\n",
	"	For example the difference between $80 and $90 is 12.5%\n",
	"\n",
	"The log2 ratio is is log2(B/A).\n",
	"	This is useful for proportions/percentages.\n",
	"\n",
	"There are also 2 special functions: pick_A(x) and pick_B(x)\n",
	"	These two functions will return the values from A or B.\n",
	"	This will allow you to explicity keep the values from A or B.\n",
	"	This is useful if you aren't sure what pairs of sample values\n",
	"	go together or you want to explicitly rename a variable name.\n",
	"	For example:  As_x=pick_A(x);\n",
	"	              Bs_x=pick_A(x);\n",
	"\n",
	"\n");

if(!length(opt$input) || 
	!length(opt$paired_map) || 
	!length(opt$a_colname) || 
	!length(opt$b_colname) || 
	!length(opt$formulas) || 
	!length(opt$output)){
	cat(usage);
	q(status=-1);
}

InputFName=opt$input;
Formulas=opt$formulas;
OutputFName=opt$output;
PairingsFile=opt$paired_map;
AColname=opt$a_colname;
BColname=opt$b_colname;
SampIDColname="";


cat("Input Filename: ", InputFName, "\n");
cat("Sample ID Colname: ", SampIDColname, "\n");
cat("Formulas: ", Formulas, "\n");
cat("Output Filename: ", OutputFName, "\n");
cat("Pairings Map: ", PairingsFile, "\n");
cat("A Colname: ", AColname, "\n");
cat("B Colname: ", BColname, "\n");

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

load_commands=function(fname){
	commands=readLines(fname);	
	return(commands);
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

nomn_diff=function(val){
	return(val$b-val$a);
}

perc_diff=function(val){
	return(100*(val$b-val$a)/val$a);
}

log2_ratio=function(val){
	return(log2(val$b/val$a));
}

pick_A=function(val){
	return(val$a);
}

pick_B=function(val){
	return(val$b);
}


##############################################################################

# Load factors
factors=load_factors(InputFName);
if(SampIDColname!=""){
	rownames(factors)=factors[,SampIDColname];
}else{
	rownames(factors)=factors[,1];
}
factor_rownames=rownames(factors);

num_samples=nrow(factors);

# Load Pairings
all_pairings_map=load_mapping(PairingsFile, AColname, BColname);
num_pairs=nrow(all_pairings_map);
print(all_pairings_map);

# Load commands
commands=load_commands(Formulas);
num_commands=length(commands);

cat("Number of Commands: ", num_commands, "\n");

print(commands);
for(cmd in commands){
	cat("\n");

	cmd=strsplit(cmd, "#")[[1]][1]; # Strip off comments
	if(is.na(cmd) || length(cmd)==0 || cmd == ""){
		next;
	}

	cmd=gsub("^\\s+", "", cmd);
	cmd=gsub("\\s+$", "", cmd);

	print(cmd);

	cat("Working on: ", cmd, "\n");

	# Add variable to factors
	#cmd=gsub("\\s+", "", cmd);
	lhs=strsplit(cmd, "=")[[1]][1];
	cat("LHS: ", lhs, "\n");

	rhs=strsplit(cmd, "=")[[1]][2];
	cat("Function Call:\n");
	print(rhs);

	param=gsub(".+\\(", "", rhs);
	param=gsub("\\);*$", "", param);
	cat("Param:", param, "\n");

	fname=gsub("\\(.+\\);*$", "", rhs);
	cat("Function Name: ", fname, "\n");
	funct=get(fname);

	# For each pair, gather metadata for that parameter
	# and perform the difference calculation

	result_arr=rep(NA, num_samples);
	names(result_arr)=factor_rownames;
	
	for(pair_ix in 1:num_pairs){

		pair=all_pairings_map[pair_ix, ];
		#cat("Working on Pair:\n");

		pair_factorsA=factors[pair[1], param];
		pair_factorsB=factors[pair[2], param];

		valpair=list();
		valpair$a=pair_factorsA;
		valpair$b=pair_factorsB;
		#print(valpair);

		diff_res=funct(valpair);
	
		result_arr[pair[1]]=diff_res;
		result_arr[pair[2]]=diff_res;
	}

	cnames=colnames(factors);
	factors=cbind(factors, result_arr);
	colnames(factors)=c(cnames, lhs);

	print(factors);
	
	cat("ok.\n");
}

cat("\n");

# Report any columns with no information
num_factors=ncol(factors);
factor_names=colnames(factors);
for(i in 1:num_factors){
	na_ix=is.na(factors[,i]);
	uniqs=unique(factors[!na_ix,i]);
	
	if(length(uniqs)<=1){
		cat("Warning: \"", factor_names[i], "\" has no information. (All \"", uniqs, "\")\n", sep="");
	}	
}
cat("\n");

write_factors(OutputFName, factors);

##############################################################################

cat("\nDone.\n");

print(warnings());
q(status=0);
