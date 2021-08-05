#!/usr/bin/env Rscript

###############################################################################

library('getopt');

options(useFancyQuotes=F);

params=c(
	"main_filen", "F", 1, "character",
	"main_keycn", "K", 1, "character",
	"main_offcn", "D", 1, "character",
	"main_dateformat", "T", 1, "character",
	"aux_filen", "f", 1, "character",
	"aux_keycn", "k", 1, "character",
	"aux_offcn", "d", 1, "character",
	"aux_dateformat", "t", 1, "character",
	"aux_valuecn", "v", 1, "character",
	"outputfn", "o", 1, "character",
	"output_cname", "c", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"\n",
	" (All parameters are required.)\n",
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
	"	--aux_valuecn <column name with value to link>\n",
	"\n",
	"	-o <output map file name>\n",
	"	[-c <appended values column name, default=\"interp_<aux_valuecn>\">]\n",
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
	!length(opt$aux_valuecn) || 
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
AuxValueColname=opt$aux_valuecn;
OutputFilename=opt$outputfn;

if(length(opt$output_cname)){
	OutputColName=opt$output_cname;
}else{
	OutputColName=paste("interp_", AuxValueColname, sep="");
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
cat("	Values:         ", AuxValueColname, "\n", sep="");
cat("\n");
cat("Output Filename:   ", OutputFilename, "\n", sep="");
cat("	Colname:        ", OutputColName, "\n", sep="");
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
#print(aux_matrix);
aux_data=aux_matrix[, c(AuxKeyColname, AuxOffsetColname, AuxValueColname), drop=F];
aux_nrows=nrow(aux_matrix);
cat("Aux: Num Rows to Select From:", aux_nrows, "\n");
#print(aux_data);

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
	print(aux_data);
}

###############################################################################

interpolate_values=function(values, offsets, target_offset){

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
		return(mean(values[exact_matches]));
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
		# Calculate weighted average
		prop_before=(target_offset-before_bracket_offset)/
			(after_bracket_offset-before_bracket_offset);

		interpolated_value=prop_before*before_bracket_value + 
			(1-prop_before)*after_bracket_value;

		cat("Prop of Before to use: ", prop_before, "\n");
		cat("Interpolated Value: ", interpolated_value, "\n");

		return(interpolated_value);
	}

}

# Keep values in same order as main input file
selected_values=numeric(main_nrows);

for(i in 1:main_nrows){
	cur_key=main_data[i, MainKeyColname];
	cur_date=main_data[i, MainOffsetColname];
	cur_offset=main_data[i, "main_date_conv"];
	cat("Working on: ", cur_key, " [Offset: ", cur_date, " / ", cur_offset, "]\n", sep="");

	# Select out cur_key / subject_id
	per_key_aux_ix=aux_data[,AuxKeyColname]==cur_key;
	per_key_aux=aux_data[per_key_aux_ix,];

	# Get available offsets and values
	offsets=per_key_aux[,"aux_date_conv"];
	values=per_key_aux[, AuxValueColname];

	#print(per_key_aux);
	selected_values[i]=interpolate_values(values, offsets, cur_offset)

	cat("\n");
}

###############################################################################

if(OutputColName!=""){
	appended_cname=paste("interp_", AuxValueColname, sep="");
}else{
	appended_cname=OutputColName;
}

orig_cnames=colnames(main_matrix);
new_cnames=c(orig_cnames, appended_cname);

cat("Output Colnames:\n");
print(new_cnames);

main_appended=cbind(main_matrix, selected_values);
colnames(main_appended)=new_cnames;
print(main_appended);

write_factors(OutputFilename, main_appended);

###############################################################################

cat("Done.\n");
#dev.off();
print(warnings());
q(status=0);
