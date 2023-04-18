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

	qry_formula=gsub(" ", "", qry_formula);
	qry_split_formula=parse_formula(qry_formula);
	qry_ac=count_atoms(qry_split_formula);

	atom_colnames=colnames(ref_db_mat);
	num_db_col=ncol(ref_db_mat);
	num_db_row=nrow(ref_db_mat);
	db_rownames=rownames(ref_db_mat);

	search_vector=rep(0, num_db_col);
	names(search_vector)=atom_colnames;

	for(atom_name in names(qry_ac)){
		search_vector[atom_name]=qry_ac[[atom_name]];
	}
	#print(search_vector);

	#diffs=numeric(num_db_row);
	#for(i in 1:num_db_row){
	#	cat(i, "\n");
	#	diffs[i]=sum(abs(search_vector-ref_db_mat[i,]));
	#}

	diffs=apply(ref_db_mat, 1, function(x){
			sum(abs(x-search_vector));
		});

	order_ix=order(diffs);
	top_hits=order_ix[1:num_top_hits];

	#for(i in 1:5){
	#	cat(diffs[order_ix[i]], db_rownames[order_ix[i]], "\n");
	#}

	res=list();
	res[["diffs"]]=diffs[top_hits];
	res[["hits"]]=db_rownames[top_hits];
	
	return(res);

}
	


#------------------------------------------------------------------------------

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

in_file=read.table(file=InputFile, header=T, sep="\t", stringsAsFactors=F);

num_queries=nrow(in_file);

cat("Num Queries/Lines: ", num_queries, "\n");

matches=character(num_queries);
for(i in 1:num_queries){

	qry_formula=in_file[i, FormulaColname];
	cat("Looking up: ", qry_formula, "\n");
	if(qry_formula =="" || is.na(qry_formula)){
		cat("Not query-able...\n");
		next;
	}

	results=lookup(qry_formula, db_mat);
	
	matches[i]=paste("[", results[["diffs"]][1], "] ", results[["hits"]][1], sep="");

}

out_mat=cbind(in_file, matches);

#------------------------------------------------------------------------------

outfn=paste(OutputFName, ".chem_matched.tsv", sep="");
write.table(out_mat, outfn, sep="\t", row.names=F, quote=F);

#------------------------------------------------------------------------------


###############################################################################

cat("\nDone.\n");

print(warnings());
q(status=0);
