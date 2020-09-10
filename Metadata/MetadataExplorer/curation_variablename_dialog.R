
RenameDialogBoxUI=function(id, old_name, suggested_name, mesg){

	ns=NS(id);

	modalDialog(
		fluidPage(
			tags$h4(mesg),
			tags$hr(),
			tags$b("Current Variable Name:"),
			tags$p(tags$h5(old_name)),
			textInput(ns("ReNmDB.new_variable_name"), label="New Variable Name:", value=suggested_name)
		),
		footer=fluidRow(
			column(2, actionButton(ns("ReNmDB.okButton"), label="OK")),
			column(3, actionButton(ns("ReNmDB.cancelButton"), label="Cancel"))
		)
	);
}

RenameDialogBoxServer=function(id, old_name, suggested_name, current_variable_names, mesg){
	
	moduleServer(
	
		id,
		
		function(input, output, session){
		
			cat("RenameDialogBoxServer started with id: ", id, "\n");
		
			showModal(RenameDialogBoxUI(id, old_name, suggested_name, mesg));
			
			observeEvent(input$ReNmDB.cancelButton,{
				cat("ReNmDB.cancelButton pressed.\n");
				removeModal();
				#showModal(DataTransformationDialogBoxUI(id, in_values, in_varname), session);
				#updateSelectInput(session, "DTDB.transformation_select", selected=input$DTDB.transformation_select);
			});
		  
			observeEvent(input$ReNmDB.okButton, {
				cat("ReNmDB.okButton pressed.\n");
				removeModal();
				
				check_variable_name_msg=check_variable_name(input$ReNmDB.new_variable_name, current_variable_names);
				
				cur_new_variable_name=input$ReNmDB.new_variable_name;
				
				if(check_variable_name_msg!=""){
					showModal(
						BadVariableName_DialogBox(input$ReNmDB.new_variable_name, 
							check_variable_name_msg)
					);
				}else{
					removeModal();
					return(input$ReNmDB.new_variable_name);
				}
			});
			
			#--------------------------------------------------------------------------
			# BadVN (Bad variable name) Events
			observeEvent(input$BadVN.okButton, {
				removeModal();
				showModal(RenameDialogBox(in_varname));
				updateTextInput(session, "ReNmDB.new_variable_name", value=input$ReNmDB.new_variable_name);
			});
			
		}
	);
}



BadVariableName_DialogBoxUI=function(id, variable, mesg=NULL){

	ns=NS(id);

	modalDialog(
		fluidPage(
			tags$h4(paste("The variable \"", variable, "\" is invalid.", sep="")),
			tags$h5(ifelse(is.null(mesg), "", mesg)),
			tags$br(),
			tags$h5("Please try again."),
		),
		footer=fluidRow(
			column(2, actionButton(ns("BadVN.okButton"), label="OK")),
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

