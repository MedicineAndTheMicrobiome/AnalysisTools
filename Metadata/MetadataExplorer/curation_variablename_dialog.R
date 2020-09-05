

BadVariableName_DialogBox=function(variable, mesg=NULL){
	modalDialog(
		fluidPage(
			tags$h4(paste("The variable \"", variable, "\" is invalid.", sep="")),
			tags$h5(ifelse(is.null(mesg), "", mesg)),
			tags$br(),
			tags$h5("Please try again."),
		),
		footer=fluidRow(
			column(2, actionButton("BadVN.okButton", label="OK")),
		)
	);
}

check_variable_name=function(varname, existing_varnames){
# This function is just for a quick check for egregiously bad names

	if(length(grep("^[0-9]", varname))){
		return("The variable name may not start with a number.");
	}
	
	if(length(grep("\\s", varname))){
		return("The variable name may not contain any spaces.");
	}
	
	if(any(varname==existing_varnames)){
		return("The variable name is already being used.");
	}
	
	bad_chars=gsub("[a-zA-Z0-9_\\.]", "", varname);
	if(nchar(bad_chars)>0){
		return(paste(
			"The variable name contains the following bad characters:\n",
				bad_chars, sep=""));
	}
	
	return("");
}


RenameDialogBox=function(old_name, suggested_name, mesg=NULL){
	cat(file=stderr(), old_name, suggested_name, mesg, "\n", sep=",");
	modalDialog(
		fluidPage(
			ifelse(is.null(mesg), "", tags$h4(mesg)),
			tags$b("Current Variable Name"),
			tags$h5(old_name),
			textInput("ReNmDB.new_variable_name", label="New Variable Name", value=suggested_name)
		),
		footer=fluidRow(
			column(2, actionButton("ReNmDB.okButton", label="OK")),
			column(3, actionButton("ReNmDB.cancelButton", label="Cancel"))
		)
	);
}

BadVariableName_DialogBox=function(variable, mesg=NULL){
	modalDialog(
		fluidPage(
			tags$h4(paste("The variable \"", variable, "\" is invalid.", sep="")),
			tags$h5(ifelse(is.null(mesg), "", mesg)),
			tags$br(),
			tags$h5("Please try again."),
		),
		footer=fluidRow(
			column(2, actionButton("BadVN.okButton", label="OK")),
		)
	);
}

CantSave_DialogBox=function(mesg=NULL){
	modalDialog(
		fluidPage(
			tags$h4(paste("Cannot save changes: ", mesg, sep=""))
		),
		footer=fluidRow(
			column(2, actionButton("CS.okButton", label="OK"))
		)
	);
}


observe_VariableNameDialogBoxEvents=function(input, output, session){


}
