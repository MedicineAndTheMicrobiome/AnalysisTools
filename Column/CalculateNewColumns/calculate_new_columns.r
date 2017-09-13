#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);
options(width=120);

params=c(
	"input", "i", 1, "character",
	"formulas", "f", 1, "character",
	"output", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input tab-separated column data file>\n",
	"	-f <formulas file>\n",
	"	-o <output tab_separated column data file>\n",
	"	[-r <replace/overwrite columns>]\n", 
	"\n",
	"This script will read in a metadata file\n",
	"and run through the list of formulas to add\n",
	"derived data to the table.\n",
	"\n",
	"Example formulas file entries:\n",
	"	fev1.over.fvc=bdfev1/bdfvc\n",
	"	age=(dob-baldat)/365\n",
	"	smoke=(env.cigrete==\"Yes\" || env.cigar==\"Yes\" || env.pipe==\"Yes\")\n",
	"	delete env.cigrete\n",
	"	delete env.cigar\n",
	"	delete env.pipe\n",
	"\n",
	"Additional columns will be added to the end\n",
	"of the existing columns.  If the formula line is \"delete\"\n",
	"Then the column will be deleted.\n",
	"\n",
	"In addition, the following non-standard R functions have been implemented:\n",
	"	remap(x, key, value):  This will remap the keys in x to the corresponding values.\n",
	"\n",
	"	  example:\n",
	"	     gpa=remap(grade, c(\"A\", \"B\", \"C\", \"D\", \"E\"), c(4, 3, 2, 1, 0))\n",
	"\n",
	"	min.list(list(x, y, ...), na.rm=T):  This will take the minimum value between the pairs x, y, ...\n",
	"	max.list(list(x, y, ...), na.rm=T):  This will take the maximum value between the pairs x, y, ...\n",
	"	  * Note that this is necessary because R's min and max builtin function will return a scalar for min(x, y)\n",
	"\n",
	"	  example:\n",
	"	     max_long_weekend_smokes=min.list(list(saturday, sunday, monday)); \n",
	"\n",
	"\n");

if(!length(opt$input) || !length(opt$formulas) || !length(opt$output)){
	cat(usage);
	q(status=-1);
}

InputFName=opt$input;
Formulas=opt$formulas;
OutputFName=opt$output;

cat("Input Filename: ", InputFName, "\n");
cat("Formulas: ", Formulas, "\n");
cat("Output Filename: ", OutputFName, "\n");

##############################################################################

load_factors=function(fname){
	factors=data.frame(read.table(fname,  header=TRUE, check.names=FALSE, sep="\t"));

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

##############################################################################
# Additional Functions

remap=function(x, key, value){
	len=length(x);
	if(length(key)!=length(value)){
		cat("Error!  Key/Value lengths are not the same for 'remap'\n");
		quit(-1);
	}
	new=numeric();
	for(i in 1:len){
		ix=which(x[i]==key);
		if(length(ix)==0){
			new[i]=NA;
		}else{
			new[i]=value[ix];
		}
	}
	return(new);
}

min.list=function(arglist, na.rm=T){
	m=matrix(unlist(arglist), byrow=F, ncol=length(arglist));
	out=apply(m, 1, function(x){min(x, na.rm=na.rm)});	
	return(out);
}

max.list=function(arglist, na.rm=T){
	m=matrix(unlist(arglist), byrow=F, ncol=length(arglist));
	out=apply(m, 1, function(x){max(x, na.rm=na.rm)});	
	return(out);
}


##############################################################################

# Load factors
factors=load_factors(InputFName);

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

	if(length(grep("^delete ", cmd))==1){
		# Delete variable from factors
		var=strsplit(cmd, "\\s+")[[1]][2];
		cat("Deleting: ", var, "\n"); 
		cnames=colnames(factors);
		cnames=setdiff(cnames, var);
		factors=factors[,cnames];
	}else{
		# Add variable to factors
		cmd=gsub("\\s+", "", cmd);
		lhs=strsplit(cmd, "=")[[1]][1];
		cat("LHS: ", lhs, "\n");
		results=eval(parse(text=cmd), envir=factors);
		print(results);
		cnames=colnames(factors);
		factors=cbind(factors, results);
		colnames(factors)=c(cnames, lhs);
	}
}

cat("\n");

write_factors(OutputFName, factors);

##############################################################################

cat("\nDone.\n");

print(warnings());
q(status=0);
