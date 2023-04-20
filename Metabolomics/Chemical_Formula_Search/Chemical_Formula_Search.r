#!/usr/bin/env Rscript

###############################################################################

library('getopt');
options(useFancyQuotes=F);
options(width=220);

params=c(
	"input_file", "i", 1, "character",
	"formula_colname", "c", 1, "character",
	"input_formulas_text", "d", 1, "character",
	"input_formulas_matrix", "m", 1, "character",
	"output_root", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste(
	"\nUsage:\n", script_name, "\n",
	"	Use one of the input options:\n",
	"	-i <input, with formulas>\n",
	"	-c <column name with formulas>\n",
	"	[-d <input, unparsed database with names>]",
	"	[-m <input, matrix database with names>]",
	"	-o <output, filename root>\n",
	"\n",
	"\n");

if(
	!length(opt$input_file) ||
	!length(opt$formula_colname) ||
	(!length(opt$input_formulas_text) && !length(opt$input_formulas_matrix)) ||
	!length(opt$output_root)
){
	cat(usage);
	q(status=-1);
}

InputFile=opt$input_file;
FormulaColname=opt$formula_colname;
OutputFName=opt$output_root;

InputDatabaseList="";
if(length(opt$input_formulas_text)){
	InputDatabaseList=opt$input_formulas_text;
}

InputDatabaseMatrix="";
if(length(opt$input_formulas_matrix)){
	InputDatabaseMatrix=opt$input_formulas_matrix;
}

##############################################################################

parse_formula=function(formula){

	numbers=as.character(0:9);
	any_letters=c(LETTERS, letters);

	cl_form=formula;
	cl_form=gsub("\\+", "", cl_form);
	cl_form=gsub("\\-", "", cl_form);
	cl_form=gsub("\\.", "", cl_form);

	chars=strsplit(cl_form, "")[[1]];

	# Place space after lower case
	i=1;
	while(i < length(chars)){
		if(any(chars[i]==letters)){
			len=length(chars);
			chars=c(chars[1:i], " ", chars[(i+1):len]);
		}
		i=i+1;
	}

	# Place space between consequtive upper case
	i=1;
	while(i < length(chars)){
		if(any(chars[i]==LETTERS) && any(chars[i+1]==LETTERS)){
			len=length(chars);
			chars=c(chars[1:i], " ", chars[(i+1):len]);
		}
		i=i+1;
	}

	# Place space before numbers
	i=1;
	while(i < length(chars)){
		if(any(chars[i]==any_letters) && any(chars[i+1]==numbers)){
			len=length(chars);
			chars=c(chars[1:i], " ", chars[(i+1):len]);
		}
		i=i+1;
	}

	# Place space after numbers
	i=1;
	while(i < length(chars)){
		if(any(chars[i]==numbers) && any(chars[i+1]==any_letters)){
			len=length(chars);
			chars=c(chars[1:i], " ", chars[(i+1):len]);
		}
		i=i+1;
	}
	
	spaced=paste(chars, collapse="");

	resplit=strsplit(spaced, " ")[[1]];

	if(cl_form != gsub(" ", "", spaced)){
		cat("Error: parsing formula:\n");
		print(cl_form);
	}
	
	return(resplit);

}

count_atoms=function(split_form){

	#cat("Split Formula:\n");
	#print(split_form);

	num_tokens=length(split_form);

	suppressWarnings({values=as.numeric(split_form)});

	i=1;
	results=list();
	while(i<=num_tokens){
		if(is.na(values[i])){
			if(!is.na(values[i+1])){
				results[[split_form[i]]]=values[i+1];
			}else{
				results[[split_form[i]]]=1;
			}
		}
		i=i+1;
	}

	#print(results);
	#cat("--------------------------------------\n");
	return(results);

}

load_formula_database=function(fname){

	cat("Loading Formula Database: ", fname, "\n", sep="");

	data=data.frame(read.table(fname,  header=F, check.names=FALSE, as.is=T, 
		comment.char="", quote="", sep="\t"));

	colnames(data)=c("Name", "Formula", "CAS");
	dimen=dim(data);
	cat("Rows Loaded: ", dimen[1], "\n");
	cat("Cols Loaded: ", dimen[2], "\n");

	nrows=dimen[1];

	atom_counts_list=list();
	molecule_names_list=list();

	unique_atoms=c();
	ticks=nrows/100;
	for(i in 1:nrows){
		formula=data[i,"Formula"];	
		if(formula==""){next;}
		if(length(grep("\\[", formula))){ next; }
		if(length(grep("\\]", formula))){ next; }
		if(length(grep("\\(", formula))){ next; }
		if(length(grep("\\)", formula))){ next; }
		if(length(grep(">", formula))){ next; }
		if(length(grep("<", formula))){ next; }

		split_formula=parse_formula(formula);
		ac=count_atoms(split_formula);

		atom_counts_list[[formula]]=ac;

		if(is.null(molecule_names_list[[formula]])){
			molecule_names_list[[formula]]=data[i,"Name"];
		}else{
			molecule_names_list[[formula]]=c(
				molecule_names_list[[formula]], data[i,"Name"]);
		}

		unique_atoms=unique(c(unique_atoms, names(ac)));

		if(!(i%%ticks)){
			cat(".");
		}
	}
	cat("\n");

	unique_atoms=sort(unique_atoms);
	num_unique_atoms=length(unique_atoms);
	cat("Unique 'atoms' found in ", fname, "\n");
	print(unique_atoms);


	num_entries=length(atom_counts_list);
	kept_formulas=names(atom_counts_list);
	entry_names=character(num_entries);
	for(i in 1:num_entries){
		formula=kept_formulas[i];
		entry_names[i]=paste(
			formula, ": ",
			paste(molecule_names_list[[formula]], collapse=";"), sep=""
			);
	}
	
	atom_matrix=matrix(0, nrow=num_entries, ncol=num_unique_atoms);
	rownames(atom_matrix)=entry_names;
	colnames(atom_matrix)=unique_atoms;

	cat("Loading Matrix with Parsed Counts...\n");
	for(i in 1:num_entries){
		formula=kept_formulas[i];
		counts=atom_counts_list[[formula]];
		atom_names=names(counts);
		atom_matrix[i, atom_names]=unlist(counts[atom_names]);
	}

	return(atom_matrix);
}

###############################################################################

lookup=function(qry_form, ref_db_mat, num_top_hits=5){

	# Parse the query string
	qry_formula=gsub(" ", "", qry_formula);
	qry_split_formula=parse_formula(qry_formula);
	qry_ac=count_atoms(qry_split_formula);
	#print(qry_ac);

	# Get DB metadata
	atom_colnames=colnames(ref_db_mat);
	num_db_col=ncol(ref_db_mat);
	db_rownames=rownames(ref_db_mat);
	num_db_row=nrow(ref_db_mat);

	# Make sure the query doesn't have atoms missing from database
	extra_atoms=setdiff(names(qry_ac), atom_colnames);
	if(length(extra_atoms)){
		extra_col=matrix(0, nrow=num_db_row, ncol=length(extra_atoms));
		colnames(extra_col)=extra_atoms;
		ref_db_mat=cbind(ref_db_mat, extra_col);
		cat("Zero padded atoms missing in DB:\n");
		print(extra_atoms);
		atom_colnames=colnames(ref_db_mat);
		num_db_col=ncol(ref_db_mat);
	}

	# Create query vector that matches DB
	search_vector=rep(0, num_db_col);
	names(search_vector)=atom_colnames;
	for(atom_name in names(qry_ac)){
		search_vector[atom_name]=qry_ac[[atom_name]];
	}
	print(search_vector);

	# Calculated edit distance across all the matrix rows
	diffs=apply(ref_db_mat, 1, function(ref){
			diff=search_vector-ref;
			subtractions=sum(diff[diff<0]);
			additions=sum(diff[diff>0]);
			tot_edits=additions+abs(subtractions);
			num_atoms=sum(ref);			
			net_prop_change=round(tot_edits/num_atoms, 5);	
			return(c(tot_edits, additions, subtractions, net_prop_change));
		});

	diff_tab=t(diffs);
	colnames(diff_tab)=c("Edits", "Adds", "Subs", "dProp");

	# Sort
	order_ix=order(diff_tab[,"Edits"]);
	diff_tab=diff_tab[order_ix,];

	res=list();
	res[["diffs"]]=diff_tab[1:num_top_hits,,drop=F];
	res[["hits"]]=rownames(diff_tab)[1:num_top_hits];

	return(res);

}
	
###############################################################################
# Main program

# Load database either from unparsed list, or previously parsed matrix
if(InputDatabaseMatrix!=""){
	cat("Reading atoms matrix from: ", InputDatabaseMatrix, "\n");
	db_mat=read.table(file=InputDatabaseMatrix, header=T, sep=" ", quote="\"");
}else{
	cat("Reading atoms list from: ", InputDatabaseList, "\n");
	db_mat=load_formula_database(InputDatabaseList);
	db_mat_fn=paste(InputDatabaseList, ".mat", sep="");
	write.table(db_mat, file=db_mat_fn, row.names=T);
	cat("Wrote ", nrow(db_mat), " entries to ", db_mat_fn, "\n");
}

num_entries=nrow(db_mat);
cat("Num Entries Read: ", num_entries, "\n");

# Sanity check results
if(0){
	for(i in 1:num_entries){
		row=db_mat[i,,drop=F];
		nozero=row>0;
		print(rownames(row));
		print(row[1, nozero]);
		cat("\n\n");
	}
}

#------------------------------------------------------------------------------

# Read in file that contains query chemical formulas to annotate
in_file=read.table(file=InputFile, header=T, sep="\t", stringsAsFactors=F);

num_queries=nrow(in_file);
cat("Num Queries/Lines: ", num_queries, "\n");

matches=character(num_queries);
additional_match_info=list();

for(i in 1:num_queries){

	qry_formula=in_file[i, FormulaColname];
	cat("--------------------------------------------------\n");
	cat("Looking up: ", qry_formula, "\n");
	if(qry_formula =="" || is.na(qry_formula)){
		cat("Not query-able...\n");
		next;
	}

	results=lookup(qry_formula, db_mat);

	top_match=results[["diffs"]][1,,drop=F];
	top_name=rownames(top_match);
	
	# Save best match
	matches[i]=paste("[", top_match[,"Edits"], "] ", top_name, sep="");

	cat("\n");
	cat("Top match: ", top_match[, "Edits"], " edits\n", sep="");
	cat("\t", top_name, "\n");
	
	additional_match_info[[qry_formula]]=results;

	cat("\n\n");
}

out_mat=cbind(in_file, matches);
last_col=ncol(out_mat);
colnames(out_mat)=c(colnames(out_mat)[1:(last_col-1)], "[edits] DBChemName");

#------------------------------------------------------------------------------
# Output query file with addition column of annotation

outfn=paste(OutputFName, ".chem_matched.best.tsv", sep="");
write.table(out_mat, outfn, sep="\t", row.names=F, quote=F);

#------------------------------------------------------------------------------
# Output additional matches by query formula

outfn=paste(OutputFName, ".chem_matched.additional_hits.tsv", sep="");

fh=file(outfn, "w");
for(qry_form in names(additional_match_info)){
	cat(file=fh, qry_form, "\n", sep="");	
	hits=additional_match_info[[qry_form]]$diffs;

	num_hits=nrow(hits);
	hit_names=rownames(hits);
	rownames(hits)=1:num_hits;

	cat(file=fh, "\t", paste(c(colnames(hits), "Name"), collapse="\t"), "\n");
	for(i in 1:num_hits){
		cat(file=fh, "\t", paste(c(hits[i,], hit_names[i]), collapse="\t"), "\n");
	}

	cat(file=fh, "\n\n");

}
close(fh);

###############################################################################

cat("\nDone.\n");

print(warnings());
q(status=0);
