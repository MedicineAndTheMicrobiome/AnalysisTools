
library(shiny)
library(tableHTML)

###############################################################################
# All the dialog boxes and UI defined here

BooleanConversionDialogBoxUI=function(id, in_values, in_varname){		

	ns=NS(id);
	
	uniq_values=sort(unique(in_values));
		
	modalDialog(
		fluidPage(
			 titlePanel(paste("Boolean Conversion: ", in_varname, sep="")),
			 helpText("Please convert your 2 level variable into a Boolean variable."),
			 fluidRow(
				column(3,
					radioButtons(ns("BCDB.reference_radio"), 
								label = "Pick a Reference:",
								choices = uniq_values
							  )
				),
				column(5,
					tableOutput(ns("bool_stats_table"))
				),
				column(4,
					textInput(ns("BCDB.new_variable_name_textInput"), "Example New Variable Name:", value="")
				)
			)
		),
		footer=fluidRow(
				column(1, actionButton(ns("BCDB.helpButton"), label="Help")),
				column(7),
				column(2, actionButton(ns("BCDB.saveButton"), label="Save")),
				column(2, actionButton(ns("BCDB.cancelButton"), label="Cancel"))
			)
		
	)
}


CantSave_DialogBox=function(id, mesg=NULL){

	ns=NS(id);

	modalDialog(
		fluidPage(
			tags$h4(paste("Cannot save changes: ", mesg, sep=""))
		),
		footer=fluidRow(
			column(2, actionButton(ns("CS.okButton"), label="OK"))
		)
	);
}

BCDB.HelpDialog=function(id){

	ns=NS(id);
	
	modalDialog(title="Boolean Conversion", 
		BCDB.help_txt, 
		footer=actionButton(ns("BCDB.dismissHelp"), label="OK"));
		
}

###############################################################################
# All the responses/event handlers here
# Do not place modalDialogs inline, unless they are really simple.

BooleanConversionDialogBoxServer=function(id, invalname, invarname, transformname, save_transf_data_callback){

	moduleServer(
	
		id,
		
		function(input, output, session){
			
			observeEvent(input$BCDB.reference_radio, {
				cat("radio button pressed.\n");
				
				in_varname=session$userData[[invarname]];
				in_values=session$userData[[invalname]];
				nas=is.na(in_values);
				nona_values=in_values[!nas];
				
				bool_levels=sort(unique(nona_values));
				
				alternate_level=setdiff(bool_levels, input$BCDB.reference_radio);
				
				suggested_variable_name=paste(
					get_prefix(in_varname), "_", in_varname, "_", alternate_level, sep="");

				updateTextInput(session, "BCDB.new_variable_name_textInput", value=suggested_variable_name);
				
				# Build/Refresh Matrix
				mat=matrix(character(), nrow=3, ncol=3);
				colnames(mat)=c(in_varname, "n", "Boolean");
				mat[,in_varname]=c(bool_levels, "NA");
				mat[,"n"]=c(sum(nona_values==bool_levels[1]), sum(nona_values==bool_levels[2]), sum(nas));
				mat[,"Boolean"]=c(input$BCDB.reference_radio!=bool_levels, "NA");	
				output$bool_stats_table=renderTable({mat}, colnames=T);
			});
			
			observeEvent(input$BCDB.new_variable_name_textInput, {
				cat("User modified text input.\n");
			});
			
			observeEvent(input$BCDB.helpButton, {
				cat("Dialog Help Pushed.\n");
				removeModal();
				showModal(BCDB.HelpDialog(id));
			});
			
			observeEvent(input$BCDB.cancelButton, {
				cat("Dialog Cancel Pushed.\n");
				removeModal();
			});
			
			observeEvent(input$BCDB.saveButton, {
				cat("Dialog Save Pushed.\n");
				
				removeModal();
				
				showModal(RenameDialogBoxUI(
					id="bool_rename", 
					old_name=session$userData[[invarname]], 
					suggested_name=input$BCDB.new_variable_name_textInput, 
					mesg="Please specify a new name for your boolean converted data."));
					
				session$userData[["curate.cancel_call_back"]]=function(){
						updateTextInput(session, "BCDB.new_variable_name_textInput", value=input$BCDB.new_variable_name_textInput);
						updateSelectInput(session, "BCDB.reference_radio", selected=input$BCDB.reference_radio);
					}

			});
			
			# Help Events
			observeEvent(input$BCDB.dismissHelp, {
				removeModal();
				showModal(BooleanConversionDialogBoxUI(id, session$userData[[invalname]], session$userData[[invarname]]));
				updateTextInput(session, "BCDB.new_variable_name_textInput", value=input$BCDB.new_variable_name_textInput);
				updateSelectInput(session, "BCDB.reference_radio", selected=input$BCDB.reference_radio);
			});

			# CS, Can't save Events
			observeEvent(input$CS.okButton, {
				showModal(BooleanConversionDialogBoxUI(id, session$userData[[invalname]], session$userData[[invarname]]));
				updateSelectInput(session, "BCDB.reference_radio", selected=input$BCDB.reference_radio);
			});
			
			# BadVN (Bad variable name) Events
			observeEvent(input$BadVN.okButton, {
				removeModal();
				showModal(RenameDialogBox(session$userData[[invarname]]));
				updateTextInput(session, "ReNmDB.new_variable_name", value=input$ReNmDB.new_variable_name);
			});
		
		}
	);
}

###############################################################################
# Supporting constants and functions specific to this group of UIs

get_prefix=function(name){
	# Suggest prefix to variable name
	if(length(grep("s$", name))){
		prefix="are";
	}else{
		prefix="is";
	}
	return(prefix);
}

###############################################################################
# Use this function so we can see what is being returned.
# As we share code, we will look for this call to quickly identify what is being
# returned by this set of functions.  This code below will not be executed
# in the main application.

SetReturn=function(varname, value, session){
	updateTextInput(session, varname, value=value); 
}

###############################################################################
# Static variables
BCDB_static=list();
BCDB_static[["selected_reference"]]="";
BCDB_static[["new_variable_name"]]="";
BCDB_static[["levels"]]="";
BCDB_static[["variable_name"]]="";

BCDB.help_txt=tags$html(
	tags$h4(tags$b("Pick a Reference:")),
		tags$p(
			"Select the data value to use as the reference. ",
			"The reference will be assigned a value of 0 (False)", 
			"and the alternate value will be assigned a value of 1 (True)"),
	tags$h4(tags$b("Statistics:")),
		tags$ul(
			tags$li(tags$b("Categories"), ": The factor levels/values of the variable."),
			tags$li(tags$b("n"), ": The number of samples with that data value."),
			tags$li(tags$b("Boolean"), ": The Boolean representation of the data value.")),				
	tags$h4(tags$b("New Variable Name:")),
		"A new variable name will be suggested, but you can modify this."
);

if(!exists("integration")){

	###############################################################################
	# This will be provided, for now hard coded for testing
	if(1){
		in_values=c("Apple", "Orange", "Orange", "Orange", "Apple", "Apple", "Orange");
		in_varname="FruitType";
	}else{
		in_values=c(8, 8, 16, 8, 16, 16, 16, 8, 8);
		in_varname="Ounces";
	}

	print(in_varname);

	in_current_variable_names=c("apple_pie", "peach_cobbler", "portuguese_tart", "shoofly_pie");

	###############################################################################
	# Boiling plate for unit testing

	ui = fluidPage(
		actionButton("startButton", label="Start"),
		textInput("Selected_Reference", label="Selected Reference:"),
		textInput("Selected_Variable_Name", label="Selected Variable Name:")
	);

	###############################################################################

	server = function(input, output, session) {

		observeEvent(input$startButton,{
			showModal(BooleanConversionDialogBox(in_values, in_varname), session);
		});

		observe_BooleanConversionDialogBoxEvents(input, output, session, in_values, in_varname, in_current_variable_names);
	}

	###############################################################################

	# Launch
	shinyApp(ui, server);
}