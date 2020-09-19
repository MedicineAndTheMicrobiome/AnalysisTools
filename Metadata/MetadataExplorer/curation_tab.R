library(shiny)
library(DT)

source("D:\\work_git\\AnalysisTools\\Metadata\\MetadataExplorer\\curation_table_functions.R");

source("D:\\work_git\\AnalysisTools\\Metadata\\MetadataExplorer\\curation_variablename_dialog.R");
source("D:\\work_git\\AnalysisTools\\Metadata\\MetadataExplorer\\curation_data_transformation.R");
source("D:\\work_git\\AnalysisTools\\Metadata\\MetadataExplorer\\curation_date_reformat_dialog.R");
source("D:\\work_git\\AnalysisTools\\Metadata\\MetadataExplorer\\curation_convert_to_boolean.R");


CurationTab.add_transformed_column=function(new_col_name, session){

	cat("Adding: ", new_col_name, "\n");
	cat("Using Transform: ", session$userData[["current_transform_type"]], "\n");
	
}


CurationTab.get_existing_variable_names=function(session){
	return(colnames(session$userData[["Metadata"]]));
}

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
		
		session$userData[["curation_values"]]=values;
		session$userData[["curation_varname"]]=varname;
		session$userData[["avail_varnames"]]=colnames(session$userData[["Metadata"]]);
		
		
		switch(warning_code,
			"INV_VAR_NAME"={
				suggested_name=gsub("[^a-zA-Z0-9_\\.]", "", varname);
				#RenameDialogBoxServer("bad_var_name", varname, suggested_name, used_names, mesg=NULL);
			},
			"NA_GT50"={
			},
			"NA_GT25"={
			},
			"NA_GT10"={
				suggested_name=gsub("[^a-zA-Z0-9_\\.]", "", varname);
				
				#RenameDialogBoxServer("nagt10", varname, suggested_name, used_names, mesg=NULL);
			},
			"DICH_NOTBOOL"={
				#BooleanConversionDialogBoxServer("bool_conv", values, varname);
			},
			"ALL_IDENT"={
			},
			"REC_SQRT_TRANS"={
				session$userData[["id_name"]]="data_trans";
				showModal(DataTransformationDialogBoxUI("data_trans", values, varname, default_trans="sqrt(x)"));
			},
			"REC_LOG_TRANS"={
				session$userData[["id_name"]]="data_trans";
				showModal(DataTransformationDialogBoxUI("data_trans", values, varname, default_trans="ln(x)"));
			},
			"REC_LOGIT_TRANS"={
				session$userData[["id_name"]]="data_trans";
				showModal(DataTransformationDialogBoxUI("data_trans", values, varname, default_trans="logit(x)"));
			},
			"NON_ISODATE"={
				session$userData[["id_name"]]="date_format";
				showModal(DateFormatConversionDialogBoxUI("date_format", values, varname));
			},
			"INV_LVL_NAMES"={
			},
			"EXCESS_CATS"={
			},	
			"UNDERREP_CATS"={
			}
		);
		
		selectRows(dt_proxy, NULL);
		
	});
	
	
	
	RenameDialogBoxServer(id="rename", in_values="curation_values", in_varname="curation_varname", 
		in_availname="avail_varnames", 
		return_ui=NULL,
		save_transf_data_callback=CurationTab.add_transformed_column);
		
	RenameDialogBoxServer(id="data_rename",  in_values="curation_values", in_varname="curation_varname", 
		in_availname="avail_varnames", 
		return_ui=DataTransformationDialogBoxUI,
		save_transf_data_callback=CurationTab.add_transformed_column);

	RenameDialogBoxServer(id="date_rename", in_values="curation_values", in_varname="curation_varname", 
		in_availname="avail_varnames", 
		return_ui=DateFormatConversionDialogBoxUI,
		save_transf_data_callback=CurationTab.add_transformed_column);
		
		
		
	
	DateFormatConversionDialogBoxServer(id="date_format", 
		invalname="curation_values", invarname="curation_varname", transformname="current_transform_type",
		save_transf_data_callback=CurationTab.add_transformed_column);
		
	DataTransformationDialogBoxServer("data_trans", 
		invalname="curation_values", invarname="curation_varname", transformname="current_transform_type",
		save_transf_data_callback=CurationTab.add_transformed_column);
		
	#BooleanConversionDialogBoxServer("bool_conv", save_transf_data_callback=CurationTab.add_transformed_column);
	
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