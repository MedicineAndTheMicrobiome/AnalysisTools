
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

RenameDialogBoxServer=function(id, in_varname, in_availname, save_transf_data_callback){
	
	moduleServer(
	
		id,
		
		function(input, output, session){
		
			ns=NS(id);
		
			cat("RenameDialogBoxServer called:\n");
			cat("  Session Namespace=", session$ns(""), " id=", id, " ns() prepender=", ns(""), "\n");
	
			# I think there is a bug in the namespce.  It's either not prepending ns of calling function to id,
			# or if the ns() is supposed to be placed around the id when called, it's double appending the ns.
			id=gsub("-$", "", session$ns(""));
			cat("  Used id=", id, "\n");
			
			#ret_val=reactiveValues();
		
			#showModal(RenameDialogBoxUI(id, old_name, suggested_name, mesg));
			
			observeEvent(input$ReNmDB.cancelButton,{
				cat("ReNmDB.cancelButton pressed.\n");
				removeModal();
				
				#return(
				#	list(
				#		new_name=reactive({input$ReNmDB.new_variable_name}), 
				#		cancel=reactive({input$ReNmDB.cancelButton})
				#));
				#showModal(DataTransformationDialogBoxUI(id, in_values, in_varname), session);
				#updateSelectInput(session, "DTDB.transformation_select", selected=input$DTDB.transformation_select);
			});
		  
			observeEvent(input$ReNmDB.okButton, {
				cat("ReNmDB.okButton pressed.\n");
				
				check_variable_name_msg=check_variable_name(input$ReNmDB.new_variable_name, session$userData[["in_availname"]]);
				
				cat(date(), "ok button: new_variable_name: ", input$ReNmDB.new_variable_name, "\n");
				
				cur_new_variable_name=input$ReNmDB.new_variable_name;
				
				if(check_variable_name_msg!=""){
					
					removeModal();
					showModal(
						BadVariableName_DialogBox(input$ReNmDB.new_variable_name, 
							check_variable_name_msg)
					);
				}else{
					cat("Acceptable new name: ", input$ReNmDB.new_variable_name, "\n");
					
					save_transf_data_callback(
						new_col_name=input$ReNmDB.new_variable_name, 
						session=session);
					
					removeModal();
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

