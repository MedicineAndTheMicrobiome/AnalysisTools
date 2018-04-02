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
	"	good.ids=!is.na(subj_ids)\n",
	"	keep good.ids\n",
	"	make_key subj_ids\n",
	"\n",
	"Additional columns will be added to the end\n",
	"of the existing columns.  If the formula line is \"delete\"\n",
	"Then the column will be deleted.\n",
	"\n",
	"Rows/Samples can be conditionally kept/removed by using\n",
	"the \"keep\" command.  First create/identify a variable/column that is T/F\n",
	"then specify that variable as parameter for the keep command.  Rows that\n",
	"are F will be excluded.\n",
	"\n",
	"A column can be moved to the first column position perhaps making it a\n",
	"primary key or sample id for the matrix with the \"make_key\" command.\n",
	"\n",
	"In addition, the following non-standard R functions have been implemented:\n",
	"\n",
	"	bin(x, range_name, breaks): \n",
	"		This will assign the range name to the values that fall between and outside the breaks\n",
	"\n",	
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
	"	to.bool(x):  This will convert x to upper case and then to 0/1 so R won't treat them as factors:\n",
	"			FALSE / TRUE\n",
	"			F / T\n",
	"			NO / YES\n",
	"			OFF / ON\n",
	"		For other less commonly used boolean pairs, use the remap function.\n",
	"\n",
	"	String Manipulation Functions:\n",
	"		left/right will extract the number of characters from the string.\n",
	"		chop_left/chop_right will remove the number of characters from the string.\n",
	"\n",
	"		str.left(x, n)\n",
	"		str.right(x, n)\n",
	"		str.chop_left(x, n)\n",
	"		str.chop_right(x, n)\n",
	"\n",
	"	Indexing:\n",
	"		index(start=1) will assign a number from start to number of rows.\n",
	"\n",
	"	days_from_date(x):\n",
	"		This will convert a date string, e.g. 2014-05-16, into the number of days\n",
	"		since Thursday, 1 1970.  If it is before then, the number will be negative\n",
	"		The function is a wrapper for as.Date, so you can check that to see what\n",	
	"		formats it can handle.  The value may not be useful by itself, but should\n",
	"		be used when calculating days between two events.\n",
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
	factors=data.frame(read.table(fname,  header=TRUE, check.names=FALSE, comment.char="", quote="", sep="\t"));

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

bin=function(x, range_names, breaks){
	len=length(x);
	num_ranges=length(range_names);
	num_breaks=length(breaks);

	if(num_ranges!=(num_breaks+1)){
		cat("Error! Range names and breaks don't agree. (num_breaks + 1) = num_ranges\n");
		cat("Num Ranges: ", num_ranges, "\n");
		print(range_names);
		cat("Num Breaks: ", num_breaks, "\n");
		print(breaks);
		quit(-1);
	}

	if(!all(breaks==sort(breaks))){
		cat("Error! Your breaks are out of order.  Are your range names ok?\n");
		quit(-1);
	}

	new=character(len);
	for(i in 1:len){
		if(!is.na(x[i])){
			if(x[i]>=breaks[num_breaks]){
				ix=num_ranges;
			}else{
				ix=min(which(x[i]<breaks));
			}
			new[i]=range_names[ix];
		}else{
			new[i]=NA;
		}
	}
	print(new);
	return(new);
}


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

to.bool=function(x){
	uppered=toupper(x);

	bools=rep(NA, length(x));

	ones=grep("YES", uppered);
	ones=c(ones, grep("TRUE", uppered));
	ones=c(ones, grep("ON", uppered));
	ones=c(ones, grep("T", uppered));

	zeros=grep("NO", uppered);
	zeros=c(zeros, grep("FALSE", uppered));
	zeros=c(zeros, grep("OFF", uppered));
	zeros=c(zeros, grep("F", uppered));

	bools[zeros]=0;
	bools[ones]=1;

	return(bools);

}

#------------------------------------------------------------------------------

str.left=function(x, n){
	x=as.character(x);
	substr(x, 1, n);	
}

str.right=function(x, n){
	x=as.character(x);
	last=nchar(x);
	substr(x, last-n+1, last);	
}

str.chop_left=function(x, n){
	x=as.character(x);
	last=nchar(x);
	substr(x, n+1, last);	
}

str.chop_right=function(x, n){
	x=as.character(x);
	last=nchar(x);
	substr(x, 1, last-n);	
}

#------------------------------------------------------------------------------

index=function(start=1){
	# This will generate an index for each row in the input file
	return(seq(start, num_samples));
}

#------------------------------------------------------------------------------

days_from_date=function(x){
	# This will convert a date string to the number of days since Jan 1, 1970.
	return(as.numeric(as.Date(x)));
}

##############################################################################

# Load factors
factors=load_factors(InputFName);
num_samples=nrow(factors);

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
		factors=factors[,cnames, drop=F];
	}else if(length(grep("^keep ", cmd))==1){
		# Remove rows with factors that are F
		var=strsplit(cmd, "\\s+")[[1]][2];
                cat("Keeping Rows where ", var, " is TRUE\n", sep="");
		keep_ix=factors[,var];
		factors=factors[keep_ix,, drop=F];
		rowcol=dim(factors);
		cat("Rows: ", rowcol[1], " x Cols: ", rowcol[2], "\n", sep="");
	}else if(length(grep("^make_key ", cmd))==1){
		# move column to first position
		var=strsplit(cmd, "\\s+")[[1]][2];
		cat("Making ", var, " key column.\n", sep="");
		key_col_val=factors[, var, drop=F];
		cnames=colnames(factors);
		cnames=setdiff(cnames, var);
		factors=cbind(key_col_val, factors[,cnames,drop=F]);
	}else{
		# Add variable to factors
		#cmd=gsub("\\s+", "", cmd);
		lhs=strsplit(cmd, "=")[[1]][1];
		cat("LHS: ", lhs, "\n");
		results=eval(parse(text=cmd), envir=factors);
		print(results);
		cnames=colnames(factors);
		factors=cbind(factors, results);
		colnames(factors)=c(cnames, lhs);
	}
	
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
