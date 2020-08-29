library(shiny)
library(DT)

source("D:\\work_git\\AnalysisTools\\Metadata\\MetadataExplorer\\curation_table_functions.R");

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
		print(variable_info);

		warnings_table=WarningsTable.GenerateFromVariableInfo(variable_info);
		print(warnings_table);
		
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