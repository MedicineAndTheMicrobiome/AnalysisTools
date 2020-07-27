
library(xml2);


###############################################################################
# Model Rec Functions

ModelRec.init=function(){
	ModelRec=list();
	ModelRec[["Excluded"]]=c();
	ModelRec[["Available"]]=c();
	ModelRec[["Covariates"]]=c();
	ModelRec[["Groups"]]=list();
	ModelRec<<-ModelRec;
	return;
}

ModelRec.write_model=function(filename){

	# Build Excluded, Available and Covariates XML
	ml=list();
	for(vartype in c("Excluded", "Available", "Covariates")){
		ml[[vartype]]=list();
		
		ix=0;
		for(vname in ModelRec[[vartype]]){
			ix=ix+1;
			ml[[vartype]][[paste("variables",ix,sep="")]]=list(vname);
		}
		names(ml[[vartype]])=rep("name", ix);
		
	}
	
	# Build Groups Lists XML
	gr=list();
	grp_names=names(ModelRec[["Groups"]]);
	for(gname in grp_names){
		gr[[gname]]=structure(list(), id=gname);
		
		ix=0;
		for(vname in ModelRec[["Groups"]][[gname]]){
			ix=ix+1;
			gr[[gname]][[paste("variables",ix,sep="")]]=list(vname);
		}
		names(gr[[gname]])=rep("name", ix);
		
	}
	names(gr)=rep("Groups", length(grp_names));
	
	# Combine two XML records 
	ml=c(ml, gr);
	
	# Create a root
	root=list();
	root[["Model"]]=ml;
	
	xml_doc=as_xml_document(root);
	write_xml(xml_doc, filename);

	return;
}

ModelRec.read_model=function(filename){

	xml_doc=read_xml(filename);
	xml_list=as_list(xml_doc);
	
	ModelRec.init();
	
	for(vartype in c("Excluded", "Available", "Covariates")){
		var_list=xml_list[["Model"]][[vartype]];
		num_var=length(var_list);
		for(i in 1:num_var){
			ModelRec[[vartype]]=c(ModelRec[[vartype]], var_list[[i]][[1]][1]);
		}
		
		# Delete 
		xml_list[["Model"]][[vartype]]=c();
	}
	
	# Only remaining should be "Groups"
	ModelRec[["Groups"]]=list();
	num_groups=length(xml_list[["Model"]]);
	
	for(g in 1:num_groups){
	
		group_name=attr(xml_list[["Model"]][[g]], "id");
		num_var=length(xml_list[["Model"]][[g]]);
		
		for(v in 1:num_var){
			ModelRec[["Groups"]][[group_name]]=c(
				ModelRec[["Groups"]][[group_name]], xml_list[["Model"]][[g]][[v]][[1]]);
		}

	}
	
	ModelRec<<-ModelRec;
	return;
}

ModelRec.set_list=function(varcategory, name=NULL, varlist){

	if(varcategory=="Excluded" || varcategory=="Available" || varcategory=="Covariates"){
		ModelRec[[varcategory]]=varlist;
	}else if(varcategory=="Groups"){
		ModelRec[[varcategory]][[name]]=varlist;
	}
	ModelRec<<-ModelRec;
	return;
}

ModelRec.get_list=function(varcategory, name=NULL){

	if(varcategory=="Excluded" || varcategory=="Available" || varcategory=="Covariates"){
		varlist=ModelRec[[varcategory]];
	}else if(varcategory=="Groups"){
		varlist=ModelRec[[varcategory]][[name]];
	}

	return(varlist);
}

ModelRec.remove_list=function(varcategory, name=NULL){
	if(is.null(name)){
		ModelRec[[varcategory]]=c();
	}else{
		ModelRec[[varcategory]][[name]]=c();
	}
	ModelRec<<-ModelRec;
	return;
}

#------------------------------------------------------------------------------

if(test){
	ModelRec.init();
	ModelRec.set_list("Excluded", varlist=c("banned", "barred", "blocked", "ignored"));
	ModelRec.set_list("Available", varlist=c("accessible", "applicable", "free", "usable"));
	ModelRec.set_list("Covariates", varlist=c("age", "sex", "ethnicity"));
	ModelRec.set_list("Groups", "fruits", c("apples", "bananas", "cantaloupe", "eggplant"));
	ModelRec.set_list("Groups", "animals", c("bird", "frog", "dog", "cat"));
	ModelRec.set_list("Groups", "colors", c("red", "orange", "yellow", "green"));
	ModelRec.set_list("Groups", "numbers", c("one-hundred", "ninety-three", "eighty-six"));
	ModelRec.set_list("Groups", "shapes", c("circle", "square", "trapezoid", "cone"));
	ModelRec.set_list("Groups", "books", c("biography", "history", "computer", "self-help"));
	ModelRec.write_model("test");
	ModelRec.read_model("test");
	ModelRec.remove_list("Groups", "books");
	ModelRec.get_list("Excluded");
}

###############################################################################
# Metadata Rec Functions

MetadataRec.init=function(){
	MetadataRec=c();
}

MetadataRec.write_metadata=function(filename){
	write.table(MetadataRec, file=filename, quote=F, sep="\t", col.names=T, row.names=F);
}

MetadataRec.read_metadata=function(filename){

	cat("File Path: '", filename, "'\n", sep="");
	lc=tolower(filename);
	
	if(length(grep("\\.tsv$", lc))){
		cat("Reading as Tab-Separated File.\n");
		separator="\t";
	}else if(length(grep("\\.csv$", lc))){
		cat("Reading as Comma-Separated File.\n");
		separator=",";
	}else if(length(grep("\\.txt$", lc))){
		cat("Generic Text File specified. Trying to Autodetect.\n");
		fh=file(filename);
		line_buffer=readLines(filename);
		close(fh);
		
		if(length(grep("\\t", line_buffer))){
			separator="\t";
		}else if(length(grep(",", line_buffer))){
			separator=",";
		}else{
			cat("Error: Could not detect file separator.\n");
			return(-1);
		}
		
		cat("Detected: ", separator, "\n");
		
	}
	
	MetadataRec=read.table(file=filename, header=T, quote="", sep=separator);
	MetadataRec<<-MetadataRec;
	return(1);
}

MetadataRec.print=function(){
	print(MetadataRec);
}


###############################################################################
# Study Rec Functions

StudyRec.init=function(){

	StudyRec=list();
	
	StudyRec[["CrossSection"]]=list();
	StudyRec[["CrossSection"]][["SubjectID"]]="";
	
	StudyRec[["Longitudinal"]]=list();
	StudyRec[["Longitudinal"]][["SubjectID"]]="";
	StudyRec[["Longitudinal"]][["TimeOffsets"]]="";
	
	StudyRec[["Paired"]]=list();
	StudyRec[["Paired"]][["SubjectID"]]="";
	StudyRec[["Paired"]][["PairingCriteria"]]="";
	StudyRec[["Confirmed"]]="";
	
	StudyRec<<-StudyRec;
	return;
}

StudyRec.write_study=function(filename){
	
	rname=c("CrossSection", "Longitudinal", "Paired");
	cname=c("SubjectID", "TimeOffsets", "PairingCriteria", "Confirmed");
	
	tmp_mat=matrix("", nrow=length(rname), ncol=length(cname));
	rownames(tmp_mat)=rname;
	colnames(tmp_mat)=cname;
	
	cat("Study as Matrix:\n");
	print(tmp_mat);
	
	write.table(tmp_mat, file=filename, quote=F, sep="\t");
	cat("Table written to: ", filename, "\n");
	
	return;
}

StudyRec.read_study=function(filename){

	cat("Reading: ", filename, "\n");
	tmp_mat=read.table(file=filename, header=T, quote="", sep="\t");
	
	cat("Study as Matrix:\n");
	print(tmp_mat);
	
	StudyRec.init();
	StudyRec[["CrossSection"]][["SubjectID"]]=tmp_mat["CrossSection", "SubjectID"];
	StudyRec[["Longitudinal"]][["SubjectID"]]=tmp_mat["Longitudinal", "SubjectID"];
	StudyRec[["Longitudinal"]][["TimeOffsets"]]=tmp_mat["Longitudinal", "TimeOffsets"];
	StudyRec[["Paired"]][["SubjectID"]]=tmp_mat["Paired", "SubjectID"];
	StudyRec[["Paired"]][["PairingCriteria"]]=tmp_mat["Paired", "PairingCriteria"];
	StudyRec[["Confirmed"]]=rownames(tmp_mat)[which(tmp_mat[,"Confirmed"]==1)];
	
	StudyRec<<-StudyRec;
	return;
}

StudyRec.set_study=function(study_type, attribute_type, value, confirmed){

	if(!exists("StudyRec")){
		StudyRec.init();
	}

	StudyRec[[study_type]][[attribute_type]]=value;
	if(confirmed){
		StudyRec[["Confirmed"]]=study_type;
	}
	
	StudyRec<<-StudyRec;
	return;
}

###############################################################################
# Model Rec Functions




