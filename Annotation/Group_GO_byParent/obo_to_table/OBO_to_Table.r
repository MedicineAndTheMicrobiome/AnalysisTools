#!/usr/bin/env Rscript

###############################################################################

library(getopt);

options(useFancyQuotes=F);
options(width=120);

params=c(
	"input", "i", 1, "character",
	"output", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	-i <input .OBO file>\n",
	"	-o <mapping file>\n",
	"\n");

if(
	!length(opt$input) |
	!length(opt$output)
){
	cat(usage);
	q(status=-1);
}

InputOBOFile=opt$input;
OutputOBOMap=opt$output;

##############################################################################

load_obo_file=function(obofn){

	#----------------------------------------------------------------------
	# Read in file

	fh=file(obofn, "r");
	lines=readLines(fh);
	num_lines=length(lines);
	close(fh);
	cat("Num Lines Read: ", num_lines, "\n");

	#----------------------------------------------------------------------
	# Count number of records
	
	num_recs_exp=0;
	for(i in 1:num_lines){
		if(lines[i]=="[Term]"){
			num_recs_exp=num_recs_exp+1;
		}
	}

	#----------------------------------------------------------------------
	
	cnames=c("id", "name", "namespace", "is_a", "def");
	out_mat=matrix(NA, ncol=length(cnames), nrow=num_recs_exp);
	colnames(out_mat)=cnames;

	row_tmpl=character(length(cnames));
	names(row_tmpl)=cnames;
	
	num_entries=0;
	num_records_checked=0;
	num_obsol=0;

	parse_rec=function(rec_arr){
		#cat("--------------------------------\n");
		#print(rec_arr);

		num_recs=length(rec_arr);
		splits=strsplit(rec_arr, ": ");

		cur_row=row_tmpl;
		obs=F;

		for(rix in 1:num_recs){

			sp_arr=splits[[rix]];	
			fldnm=sp_arr[1];
			data=paste(sp_arr[2:length(sp_arr)], collapse=":");

			# If record is obsolete, ignore it
			if(fldnm=="is_obsolete"){
				obs=T;
				break;
			} 

			# If there are any fields we don't want, skip the field
			if(any(fldnm==
				c("alt_id", "synonym", "replaced_by", "comment", "xref", 
				"relationship", "subset")
			)){
				next;
			}

			# Clean up the field data if it's a is_a record
			if(fldnm=="is_a"){
				is_a_split=strsplit(data, " ! ")[[1]];
				data=is_a_split[1];
			}

			# If there are many of the same field name, concatenate together
			if(nchar(cur_row[fldnm])){
				cur_row[fldnm]=paste(cur_row[fldnm], ";", data, sep="");	
			}else{
				cur_row[fldnm]=data;	
			}
		}
		
		# If the record is obsolete, don't store it in the out matrix
		if(!obs){
			out_mat[num_entries+1,] <<- cur_row;
			num_entries <<- num_entries+1;
		}else{
			num_obsol <<- num_obsol+1;
		}

		num_records_checked <<- num_records_checked+1;
	}

	#----------------------------------------------------------------------

	# Skip until we seen first entry
	i=1;
	while(i <= num_lines){
		if(lines[i]=="[Term]"){
			record_arr=c();
			eor=F;
			i=i+1;
			#cat("Line: ", i, "\n");
			while(i <= num_lines && !eor){
				if(lines[i]==""){	
					eor=T;
					parse_rec(record_arr);
				}else{
					record_arr=c(record_arr, lines[i]);
					i=i+1;
				}
			}
		}
		i=i+1;
	}


	cat("Num Records Expected: ", num_recs_exp, "\n");
	cat("Num Records Kept: ", num_entries, "\n");
	cat("Num Records Obsol: ", num_obsol, "\n");
	cat("Num Records Checked: ", num_records_checked, "\n");

	if((num_entries+num_obsol) != num_records_checked){
		cat("WARNING: Num Entries Kept + Obsolete don't equal Records Checked.\n");
	}else{
		cat("GOOD: Kept and Obsolete Record Counts add up to Records Check.\n");
	}

	out_mat=out_mat[1:num_entries,,drop=F];
	return(out_mat);

}

##############################################################################

obo_mat=load_obo_file(InputOBOFile);
write.table(obo_mat, OutputOBOMap, quote=F, sep="\t", row.names=F);


cat("\nDone.\n");

print(warnings());
q(status=0);
