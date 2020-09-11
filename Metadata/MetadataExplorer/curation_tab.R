library(shiny)
library(DT)

source("D:\\work_git\\AnalysisTools\\Metadata\\MetadataExplorer\\curation_table_functions.R");

source("D:\\work_git\\AnalysisTools\\Metadata\\MetadataExplorer\\curation_variablename_dialog.R");
source("D:\\work_git\\AnalysisTools\\Metadata\\MetadataExplorer\\curation_data_transformation.R");
source("D:\\work_git\\AnalysisTools\\Metadata\\MetadataExplorer\\curation_date_reformat_dialog.R");
source("D:\\work_git\\AnalysisTools\\Metadata\\MetadataExplorer\\curation_convert_to_boolean.R");

CurationTab=function(){

	tabPanel("Curation",
	
		tags$h1(tags$b("Curation")),
		
		"Let us check your metadata for common issues:  ",
		actionButton(
			inputId="CurationTab.curateButton", 
			label="Review", icon=icon("tasks")
		),

		tags$hr(),
		
		fluidPage(
	
			DTOutput("CurationTab.warnings_table")
		),
		
	) #tabPanel
}

observe_CurationTabEvents=function(input, output, session){

	observeEvent(input$CurationTab.curateButton, {
	
		metadata=session$userData[["Metadata"]]
		variable_info=VariableInfo.build(metadata);
	
		warnings_table=WarningsTable.GenerateFromVariableInfo(variable_info);
		session$userData[["Curation"]][["WarningsTable"]]=warnings_table;
		
		output$CurationTab.warnings_table=
			renderDT(
				warnings_table,
				selection='single'
			);
		
		updateActionButton(session, inputId="CurationTab.curateButton", label="Refresh", icon=NULL);
			
	});
	
	observeEvent(input$CurationTab.warnings_table_rows_selected, {
		selected_row=input$CurationTab.warnings_table_rows_selected;
		cat("Row Selected: ", selected_row, "\n");
		
		warnings_table=session$userData[["Curation"]][["WarningsTable"]];
		warning_line=warnings_table[selected_row,];
		varname=warning_line["VariableName"];
		warning_code=warning_line["WarningCode"];
		
		metadata=session$userData[["Metadata"]];
		used_names=colnames(metadata);
		values=metadata[,varname]
		
		print(warning_code);
		print(values);
		print(varname);
		
		dt_proxy=dataTableProxy("CurationTab.warnings_table");
		
		switch(warning_code,
			"INV_VAR_NAME"={
				suggested_name=gsub("[^a-zA-Z0-9_\\.]", "", varname);
				RenameDialogBoxServer("bad_var_name", varname, suggested_name, used_names, mesg=NULL);
			},
			"NA_GT50"={
			},
			"NA_GT25"={
			},
			"NA_GT10"={
				suggested_name=gsub("[^a-zA-Z0-9_\\.]", "", varname);
				RenameDialogBoxServer("nagt10", varname, suggested_name, used_names, mesg=NULL);
			},
			"DICH_NOTBOOL"={
				BooleanConversionDialogBoxServer("bool_conv", values, varname);
			},
			"ALL_IDENT"={
			},
			"REC_SQRT_TRANS"={
				DataTransformationDialogBoxServer("sqrt_trans", values, varname, default_trans="sqrt(x)");
			},
			"REC_LOG_TRANS"={
				DataTransformationDialogBoxServer("log_trans", values, varname, default_trans="ln(x)");
			},
			"REC_LOGIT_TRANS"={
				DataTransformationDialogBoxServer("logit_trans", values, varname, default_trans="logit(x)");
			},
			"NON_ISODATE"={
				DateFormatConversionDialogBoxServer("date_format", values, varname, used_names);
			},
			"INV_LVL_NAMES"={
				#showModal();
			},
			"EXCESS_CATS"={
			},	
			"UNDERREP_CATS"={
			}
		);
		
		selectRows(dt_proxy, NULL);
		
	});
	
}

###############################################################################

if(!exists("integration")){

###############################################################################

	ui = fluidPage(
		mainPanel(
			tabsetPanel(
				tabPanel("Import/Export", ""),
				tabPanel("Data", ""),
				CurationTab(),
				tabPanel("Study",""),
				tabPanel("Model Explorer", "")
			)
		)
	);

	ImportExportTab.data_matrix=c();

	###############################################################################

	server = function(input, output, session) {
		observe_CurationTabEvents(input, output, session);
	}

	###############################################################################
	# Launch
	shinyApp(ui, server);

}