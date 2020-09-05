# Do not include library that are not necessary, but list all
# used libraries at the top

library(shiny)
library(tableHTML)

###############################################################################
# All the dialog boxes and UI defined here

DateFormatConversionDialogBoxUI=function(id, in_date, in_varname){

	ns=NS(id);

	format_strings=c(
		"%d%b%Y", "%m/%d/%y", "%Y-%d-%m", "%Y.%d", "%d-%b-%Y");
		
	modalDialog(
		fluidPage(
			 titlePanel(paste("Date Format Conversion: ", in_varname, sep="")),
			 helpText("Please convert your date into the standard ISO 8601 format."),
			 fluidRow(
				column(4,
					selectInput(ns("DFCD.input_dates_select"), 
							  label = "Original Date String",
							  choices = in_date,
							  selected = "",
							  size=16, selectize=F
							  )
				),
				column(4,
					selectInput(ns("DFCD.date_format_select"), 
							  label = "Select from Common Format Strings",
							  choices = format_strings,
							  selected = "",
							  size=10, selectize=F
							  ),
					textInput("DFCD.date_format_textInput", "Format String", value="")
				),
				column(4,
					selectInput(ns("DFCD.converted_dates_select"), 
							  label = "'ISO 8601' Converted",
							  choices = c(),
							  selected = "",
							  size=16, selectize=F
							  )
				)
			)
		),
		footer=fluidRow(
				column(1, actionButton(ns("DFCD.helpButton"), label="Help")),
				column(7),
				column(2, actionButton(ns("DFCD.saveAsButton"), label="Save As...")),
				column(2, actionButton(ns("DFCD.cancelButton"), label="Cancel"))
			)
		
	)
}

###############################################################################
# All the responses/event handlers here
# Do not place modalDialogs inline, unless they are really simple.

DateFormatConversionDialogBoxServer=function(id, in_date, in_varname){

	moduleServer(

		id,
		
		function(input, output, session){
		
			showModal(DateFormatConversionDialogBoxUI(id, in_date, in_varname));
		
			# DFCD Events:
			observeEvent(input$DFCD.helpButton, {
				cat("Dialog Help Pushed.\n");
				removeModal();
				showModal(modalDialog(title="Date Conversion Formats", 
					help_txt, 
					footer=actionButton("dismissHelp", label="OK")));
			});
			
			observeEvent(input$DFCD.cancelButton, {
				cat("Dialog Cancel Pushed.\n");
				removeModal();
			});
			
			observeEvent(input$DFCD.saveAsButton, {
				cat("Dialog SaveAs Pushed.\n");
				if(input$DFCD.date_format_textInput!=""){
					removeModal();
					showModal(RenameDialogBox(in_varname));
				}else{
					showModal(CantSave_DialogBox("No format string specified"));
				}
			});
			
			observeEvent(input$DFCD.date_format_select,{
				cat("observeEvent: input$DFCD.date_format_select.\n");
				fmt_str=input$DFCD.date_format_select;
				cat("Selected format string: '", fmt_str, "'\n", sep="");
				if(length(fmt_str)){
					updateTextInput(session, "DFCD.date_format_textInput", value=fmt_str);
				}
			});

			observeEvent(input$DFCD.date_format_textInput, {
				cat("Format text box changed.\n");
				txt_fmt_str=input$DFCD.date_format_textInput;
				if(length(txt_fmt_str)){
					if(txt_fmt_str==""){
						conv_date=c();
					}else{
						conv_date=tryCatch({
								as.Date(in_date, txt_fmt_str);
							});
					}
					updateSelectInput(session, "DFCD.converted_dates_select", choices=conv_date);
				}
			});
			
			# Help Events
			observeEvent(input$dismissHelp, {
				cat("Dismissed Help.\n");
				removeModal();
				showModal(DateFormatConversionDialogBox(in_date, in_varname), session);
				updateTextInput(session, "DFCD.date_format_textInput", value=input$DFCD.date_format_textInput);
			});

			# CS, Can't save Events
			observeEvent(input$CS.okButton, {
					showModal(DateFormatConversionDialogBox(in_date, in_varname), session);
			});

			# ReNmDB (Rename variable) Events
			observeEvent(input$ReNmDB.cancelButton,{
				removeModal();
				showModal(DateFormatConversionDialogBox(in_date, in_varname), session);
				updateTextInput(session, "DFCD.date_format_textInput", value=input$DFCD.date_format_textInput);
			});
		  
			observeEvent(input$ReNmDB.okButton, {
				removeModal();
				
				check_variable_name_msg=check_variable_name(input$ReNmDB.new_variable_name, current_variable_names);
				
				cur_new_variable_name=input$ReNmDB.new_variable_name;
				
				if(check_variable_name_msg!=""){
					showModal(
						BadVariableName_DialogBox(input$ReNmDB.new_variable_name, 
							check_variable_name_msg)
					);
				}else{
					SetReturn("Selected_Variable_Name", input$ReNmDB.new_variable_name, session);
					SetReturn("Selected_Date_Format", input$DFCD.date_format_textInput, session);
				}
			});
			
			# BadVN (Bad variable name) Events
			observeEvent(input$BadVN.okButton, {
				removeModal();
				showModal(RenameDialogBox(in_varname));
				updateTextInput(session, "ReNmDB.new_variable_name", value=input$ReNmDB.new_variable_name);
			});
 
		}
	);
 
}

###############################################################################
# Supporting constants and functions specific to this group of UIs

# Help textInput
date_conversion_help_table=matrix(
	paste("&nbsp", c(
		"%d", "day as a number (0-31)", "01-31",
		"%a", "abbreviated weekday", "Mon",
		"%A", "unabbreviated weekday", "Monday", 	
		"%m", "month (00-12)", "00-12",
		"%b", "abbreviated month", "Jan",
		"%B", "unabbreviated month", "January",
		"%y", "2-digit year", "07",
		"%Y", "4-digit year", "2007"
)), byrow=T, ncol=3);

rownames(date_conversion_help_table)=rep("",8);
colnames(date_conversion_help_table)=paste("&nbsp", c("Symbol", "Meaning", "Example"));

help_txt=tableHTML(date_conversion_help_table, rownames=F, widths=c(90,200,100));

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

if(!exists("integration")){

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
	conv_date=c();
	cur_new_variable_name="";

	# This will be provided, for now hard coded for testing
	in_date=c("03-may-1974", "30-oct_2013", "10-may-1492", "1-jan-2000");
	in_varname="Birthdate";
	in_current_variable_names=c("apple_pie", "peach_cobbler", "portuguese_tart", "shoofly_pie");

	###############################################################################
	# Boiling plate for unit testing

	ui = fluidPage(
		actionButton("startButton", label="Start"),
		textInput("Selected_Date_Format", label="Selected Date Format:"),
		textInput("Selected_Variable_Name", label="Selected Variable Name:")
	);

	###############################################################################

	server = function(input, output, session) {

		observeEvent(input$startButton,{
			showModal(DateFormatConversionDialogBox(in_date, in_varname), session);
		});

		observe_DateFormatConversionDialogBoxEvents(input, output, session, in_current_variable_names);
	}

	###############################################################################
	# Launch
	shinyApp(ui, server);

}
