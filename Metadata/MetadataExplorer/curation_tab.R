library(shiny)
library(DT)

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

# <input type="button" name="hello">

#reviewed_checkbox=matrix(sprintf("<input type=\"checkbox\" id=\"checkbox_%i\"/>", 1:10), nrow=10,ncol=1);
#review_button=matrix(sprintf("<input type=\"button\" value=\"review\" id=\"button_%i\"/>", 1:10), nrow=10, ncol=1);
#testmat=cbind(reviewed_checkbox, review_button);

testmat=matrix(runif(150), ncol=3);

print(testmat);

observe_CurationTabEvents=function(input, output, session){

	observeEvent(input$CurationTab.curateButton, {
	
		print(input);
	
		output$CurationTab.warnings_table=
			renderDT(
				testmat,
				selection='single'
			);
	});
	
	
	observeEvent(input$CurationTab.warnings_table_rows_selected, {
		selected_row=input$CurationTab.warnings_table_rows_selected;
		cat("Row Selected: ", selected_row, "\n");
	});
	
}

###############################################################################


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
