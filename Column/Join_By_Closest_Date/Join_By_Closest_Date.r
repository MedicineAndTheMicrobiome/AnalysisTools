#!/usr/bin/env Rscript

###############################################################################

library('getopt');

options(useFancyQuotes=F);

DEFAULT_INTERP_PREFIX="interp_";
DEFAULT_CLOSEST_PREFIX="closest_";

params=c(
	"main_filen", "F", 1, "character",
	"main_keycn", "K", 1, "character",
	"main_offcn", "D", 1, "character",
	"main_dateformat", "T", 1, "character",
	"aux_filen", "f", 1, "character",
	"aux_keycn", "k", 1, "character",
	"aux_offcn", "d", 1, "character",
	"aux_dateformat", "t", 1, "character",
	"aux_valuecn", "v", 2, "character",
	"aux_cn_fn_list", "l", 2, "character", 
	"aux_choose_closest", "h", 2, "character",
	"outputfn", "o", 1, "character",
	"output_cname_prefix", "p", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"\n",
	"	--main_filen <factor file name (to add to)>\n",
	"	--main_keycn <column name with key (e.g. subject id)>\n",
	"	--main_offcn <column name with Date/Offset>\n",
	"	--main_dateformat <eg. \"%Y-%m-%d\">\n",
	"\n",
	"	--aux_filen <factor file name (to add from>\n",
	"	--aux_keycn <column name with key (e.g. subject id)>\n",
	"	--aux_offcn <column name with Date/Offset>\n",
	"	--aux_dateformat <eg. \"%Y-%m-%d\">\n",
	"\n",
	"	One of: \n",
	"	[--aux_valuecn <column name with value to link>]\n",
	"	[--aux_cn_fn_list <filename with list of column names to link>]\n",
	"\n",
	"	Use this option for categorical or non-ordinal variables:\n",
	"	[--choose_closest]  (default: weighted average of bracketing measurements, i.e. interpolate)\n",
	"\n",
	"	-o <output map file name>\n",
	"	[-p <appended values column prefix>]\n",
	"	  interpolated default =", DEFAULT_INTERP_PREFIX, "\n",
	"         closest default      =", DEFAULT_CLOSEST_PREFIX, "\n",
	"\n",
	"This script will read in the main file and using the key (subject_id)\n",
	"and then link it with specified date with the value in the auxiliary\n",
	"file with the closest matching date.\n",
	"\n", sep="");

if(
	!length(opt$main_filen) || 
	!length(opt$main_keycn) || 
	!length(opt$main_offcn) || 
	!length(opt$main_dateformat) || 
	!length(opt$aux_filen) || 
	!length(opt$aux_keycn) || 
	!length(opt$aux_offcn) || 
	!length(opt$aux_dateformat) || 
	!length(opt$outputfn)
){
	cat(usage);
	q(status=-1);
}

MainFilename=opt$main_filen;
MainKeyColname=opt$main_keycn;
MainOffsetColname=opt$main_offcn;
MainDateFormat=opt$main_dateformat;
AuxFilename=opt$aux_filen;
AuxKeyColname=opt$aux_keycn;
AuxOffsetColname=opt$aux_offcn;
AuxDateFormat=opt$aux_dateformat;
OutputFilename=opt$outputfn;
OutputColumnPrefix=opt$output_cname_prefix;

AuxValueColname=opt$aux_valuecn;
AuxFilenameColnameList=opt$aux_cn_fn_list;
AuxChooseClosest=opt$aux_choose_closest;

if(length(AuxChooseClosest)){
	AuxChooseClosest=T;
}else{
	AuxChooseClosest=F;
}

if(!length(opt$aux_valuecn)){
	AuxValueColname="";
}

if(length(opt$output_cname_prefix)){
	OutputColNamePrefix=opt$output_cname_prefix;
}else{
	OutputColNamePrefix=ifelse(AuxChooseClosest, DEFAULT_CLOSEST_PREFIX, DEFAULT_INTERP_PREFIX);
}

cat("\n");
cat("Main:\n");
cat("	Filename:       ", MainFilename, "\n", sep="");
cat("	Key Colname:    ", MainKeyColname, "\n", sep="");
cat("	Offset Colname: ", MainOffsetColname, "\n", sep="");
cat("	Date Format:    ", MainDateFormat , "\n", sep="");
cat("\n");
cat("Auxiliary:\n");
cat("	Filename:       ", AuxFilename, "\n", sep="");
cat("	Key Colname:    ", AuxKeyColname, "\n", sep="");
cat("	Offset Colname: ", AuxOffsetColname, "\n", sep="");
cat("	Date Format:    ", AuxDateFormat, "\n", sep="");
cat("	Fname Colname List: ", AuxFilenameColnameList, "\n", sep="");
cat("\n");
cat("Choose Closest?: ", AuxChooseClosest, "\n", sep="");
cat("\n");
cat("Output Filename:   ", OutputFilename, "\n", sep="");
cat("	ColnamePrefix:  ", OutputColNamePrefix, "\n", sep="");
cat("\n");

##############################################################################
##############################################################################

load_factors=function(fname){
	factors=data.frame(read.table(fname,  sep="\t", header=TRUE, 
		as.is=T, check.names=FALSE, comment.char=""));
	factor_names=colnames(factors);

	ignore_idx=grep("^IGNORE\\.", factor_names);

	if(length(ignore_idx)!=0){
		return(factors[-ignore_idx]);
	}else{
		return(factors);
	}
}

write_factors=function(fname, table){

        dimen=dim(table);
        cat("Rows Exporting: ", dimen[1], "\n");
        cat("Cols Exporting: ", dimen[2], "\n");

        write.table(table, fname, quote=F, row.names=F, sep="\t");

}

load_list=function(fname){
	arr=as.matrix(read.delim(fname, header=F, as.is=T, stringsAsFactors=F, comment.char="#"))[,1];
	return(arr);
}


##############################################################################

target_variables_arr=c();
cat("\nIdentification of Targeted Values to Join:\n");
if(AuxFilenameColnameList!=""){
	cat("Loading Target File List: ", AuxFilenameColnameList, "\n");
	target_variables_arr=load_list(AuxFilenameColnameList);
}else if(AuxValueColname!=""){
	cat("Using single specified column name: ", AuxValueColname, "\n");
	target_variables_arr=AuxValueColname;
}else{
	cat("Error: Neither a list (via Filename) nor column name, was specified for joining.\n");
}
cat("Targeted Variable Name(s) for joining:\n");
print(target_variables_arr);
cat("\n\n");

##############################################################################

cat("Loading Main File:\n");
main_matrix=load_factors(MainFilename);
#print(main_matrix);
main_data=main_matrix[, c(MainKeyColname, MainOffsetColname), drop=F];
#print(main_data);
main_nrows=nrow(main_data);
cat("Main: Num Rows to Match:",  main_nrows, "\n");

cat("Loading Aux File:\n");
aux_matrix=load_factors(AuxFilename);
aux_data=aux_matrix[, c(AuxKeyColname, AuxOffsetColname, target_variables_arr), drop=F];
aux_nrows=nrow(aux_matrix);
cat("Aux: Num Rows to Select From:", aux_nrows, "\n");
cat("Excerpt:\n");
print(head(aux_data));
cat("\n");


###############################################################################
# Convert dates to offsets

if(MainOffsetColname!=""){
	cat("Converting Main Dates to Offsets...\n");

	main_date_conv=apply(main_data[, MainOffsetColname, drop=F], 1, 
		function(x){
			as.numeric(as.Date(x, MainDateFormat));
		}
	);

	main_data=cbind(main_data,main_date_conv);
	#print(main_data);
}

if(AuxOffsetColname!=""){
	cat("Converting Aux Dates to Offsets...\n");

	aux_date_conv=apply(aux_data[, AuxOffsetColname, drop=F], 1, 
		function(x){
			as.numeric(as.Date(x, AuxDateFormat));
		}
	);

	aux_data=cbind(aux_data,aux_date_conv);
	#print(aux_data);
}

###############################################################################

choose_values=function(values, offsets, target_offset, interpolate=T){

	#print(values);
	#print(offsets);
	#print(target_offset);

	# Sort by offset, just in case 
	offsets_sort_ix=order(offsets);
	values=values[offsets_sort_ix];
	offsets=offsets[offsets_sort_ix];

	# Find target position:
	#   - is target is after sample
	#   + is target is before sample
	relative_pos=offsets-target_offset;

	# Return exact match
	exact_matches=which(relative_pos==0);
	num_exact=length(exact_matches);
	if(num_exact>0){

		matching_values=values[exact_matches];
		matching_values=matching_values[!is.na(matching_values)];
		mean_mv=mean(matching_values);

		# If we can't take the mean across the matching values, then return first one.
		if(is.na(mean_mv)){
			matching_value=matching_values[1];
		}else{
			matching_value=mean_mv;
		}

		cat("Exact Offset/Date Match Found: ", matching_value, "\n");
		return(matching_value);
	}

	#print(relative_pos);

	before_target_ix=relative_pos<0;
	after_target_ix=relative_pos>0;

	before_ix=which(before_target_ix);
	after_ix=which(after_target_ix);

	before_bracket_ix=ifelse(length(before_ix), max(before_ix), NA);
	after_bracket_ix=ifelse(length(after_ix), min(after_ix), NA);

	cat("Bracket ix: ", before_bracket_ix, " / ", after_bracket_ix, "\n");

	before_bracket_offset=unique(offsets[before_bracket_ix]);
	after_bracket_offset=unique(offsets[after_bracket_ix]);

	cat("Offsets ix: ", before_bracket_offset, " / ", after_bracket_offset, "\n");

	before_bracket_value=unique(values[before_bracket_ix]);
	after_bracket_value=unique(values[after_bracket_ix]);

	cat("Values: ", before_bracket_value, " / ", after_bracket_value, "\n");

	cat("Bracketing Offsets: [", before_bracket_offset, ", ", after_bracket_offset, "]\n", sep="");
	cat("Bracketing Values: [", before_bracket_value, ", ", after_bracket_value, "]\n", sep="");

	if(is.na(before_bracket_ix)){
		cat("Using closest after value: ", after_bracket_value, "\n");
		return(after_bracket_value);
	}else if(is.na(after_bracket_ix)){
		cat("Using closest before value: ", before_bracket_value, "\n");
		return(before_bracket_value);
	}else{

		# If one or the other bracketed value is NA, return non-NA value.
		#   returns NA if both bracketed values are NA, else continues to interp/closest
		if(is.na(before_bracket_value)){
			cat("Before value is NA.\n")
			return(after_bracket_value);
		}else if(is.na(after_bracket_value)){
			cat("After value is NA.\n")
			return(before_bracket_value);
		}

		# Calculate weighted average
		prop_before=1-(target_offset-before_bracket_offset)/
			(after_bracket_offset-before_bracket_offset);

		cat("Prop of Before to use: ", prop_before, "\n");

		if(interpolate){

			interpolated_value=prop_before*before_bracket_value + 
				(1-prop_before)*after_bracket_value;

			cat("Interpolated Value: ", interpolated_value, "\n");

			return(interpolated_value);

		}else{

			closest_value=ifelse(prop_before<.5, before_bracket_value, after_bracket_value);

			cat("Closest Value:", closest_value, "\n");

			return(closest_value);

		}
	}

}

# Allocate empty matrix to accumulated selected values
num_targeted_columns=length(target_variables_arr);
new_columns=matrix(NA, nrow=main_nrows, ncol=num_targeted_columns);

colnames(new_columns)=paste(OutputColNamePrefix, target_variables_arr, sep="");

# Loop through targeted columns
col_ix=1;
for(aux_targeted_colname in target_variables_arr){

	cat("Working on Column: ", aux_targeted_colname, "\n");

	for(i in 1:main_nrows){
		cur_key=main_data[i, MainKeyColname];
		cur_date=main_data[i, MainOffsetColname];
		cur_offset=main_data[i, "main_date_conv"];
		cat("[", aux_targeted_colname, "] Working on Subject: ", cur_key, 
			" [Offset: ", cur_date, " / ", cur_offset, "]\n", sep="");

		# Select out cur_key / subject_id
		per_key_aux_ix=aux_data[,AuxKeyColname]==cur_key;
		per_key_aux=aux_data[per_key_aux_ix,];

		# Get available offsets and values
		offsets=per_key_aux[,"aux_date_conv"];
		values=per_key_aux[, aux_targeted_colname];

		new_columns[i, col_ix]=choose_values(values, offsets, cur_offset, interpolate=!AuxChooseClosest);

		cat("\n");
	}

	col_ix=col_ix+1;

}

###############################################################################

if(AuxChooseClosest){
	# Don't round if choosing closest, since it could be categorical
	out_cols=new_columns;
}else{
	out_cols=round(new_columns, 5);
}

main_appended=cbind(main_matrix, out_cols);

write_factors(OutputFilename, main_appended);

###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
