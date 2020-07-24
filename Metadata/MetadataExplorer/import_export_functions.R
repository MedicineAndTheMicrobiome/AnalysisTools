
###############################################################################
# Model Rec Functions


ModelRec.init=function(){

}

ModelRec.write_model=function(filename){
	ModelRec=list();

}

ModelRec.read_model=function(filename){

	return(model_info);
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




