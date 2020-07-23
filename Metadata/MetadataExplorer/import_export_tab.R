library(shiny)

ImportExportTab=function(){

	tabPanel("Import/Export",
		
		#tags$h5("Center for"),
		#tags$h4("Medicine and the Microbiome"),
		#tags$h6("University of Pittsburgh"),
		
		img(src="microbiome-logo-option.png", align="right"),
		
		tags$h1(tags$b("Metadata Explorer")),
		
		tags$h4("Welcome to the CMM's Metadata Explorer application!"),
		
		tags$hr(),
		
		fixedPage(
		
			fixedRow(
				column(4,
					tags$b("Upload your Metadata File:"),
					tags$p("If this is a first or returning session, please upload your metadata file here."),
					fileInput(
						inputId="ImportExportTab.metadata_load",
						label=NULL, #label="Upload Your Metadata File:",
						accept=c(
							".csv",
							".tsv"))
						),
				column(5,
					tags$p(tags$b("Download your Metadata File:")),
					tags$p("You may download your new metadata, if you have modified it."),
					fixedRow(
						column(5,
							downloadButton(
								outputId="ImportExportTab.metadata_download",
								label="Download Metadata"
								)),
						column(1,
							textOutput("ImportExportTab.metadata_ready"))
					)
				)
			),#fixedRow
			
			tags$hr(),
			
			fixedRow(
				column(4,
					tags$b("Upload your Study File:"),
					tags$p("You may upload a previously saved Study File (.stu) description here."),
					fileInput(
						inputId="ImportExportTab.study_load",
						label=NULL, #label="Upload Your Study File:",
						accept=c(
							".txt",
							".stu"
							))
					),
				column(5,
					tags$p(tags$b("Download your Study File:")),
					
					tags$p("You may download your Study File (.stu), after you have specified it."),
					fixedRow(
						column(5, 
							downloadButton(
								outputId="ImportExportTab.study_download",
								label="Download Study"
								)),
						column(1,
							textOutput("ImportExportTab.study_ready"))
					)
				),
			),#fixedRow
			
			tags$hr(),
			
			fixedRow(
				column(4,
					tags$b("Upload your Model File:"),
					tags$p("You may upload a previously saved Model File (.mod) description here."),
					fileInput(
						inputId="ImportExportTab.model_load",
						label=NULL, #label="Upload your Model File:",
						accept=c(
							".txt",
							".mod"))
					),
				column(5,
					tags$p(tags$b("Download your Model File:")),
					tags$p("You may download your Model File (.mod), after you have specified it."),
					fixedRow(
						column(5,
							downloadButton(
								outputId="ImportExportTab.model_download",
								label="Download Model"
								)
						),
						column(1,	
							textOutput("ImportExportTab.model_ready"))
					)
				)
			) #fixedRow
		) #fixedPage
	) #tabPanel
}

observe_ImportExportTabEvents=function(input, output, session){

	output$ImportExportTab.metadata_ready=renderText("Not Ready")
	output$ImportExportTab.study_ready=renderText("Not Ready")
	output$ImportExportTab.model_ready=renderText("Not Ready")

	observeEvent(input$ImportExportTab.file_load,{
		print(input$ImportExportTab.file_load);
	});
	
}

###############################################################################


###############################################################################

ui = fluidPage(
	mainPanel(
		tabsetPanel(
			ImportExportTab(),
			tabPanel("Data", ""),
			tabPanel("Curation", ""),
			tabPanel("Study",""),
			tabPanel("Model Explorer", "")
		)
	)
);

ImportExportTab.data_matrix=c();

###############################################################################

server = function(input, output, session) {
	observe_ImportExportTabEvents(input, output, session);
}

###############################################################################
# Launch
shinyApp(ui, server);
