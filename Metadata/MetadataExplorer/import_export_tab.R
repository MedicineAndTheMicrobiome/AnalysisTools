library(shiny)

source("D:\\work_git\\AnalysisTools\\Metadata\\MetadataExplorer\\import_export_functions.R");

ImportExportTab=function(){

	tabPanel("Import/Export",
		
		img(src="microbiome-logo-option.png", align="right"),
		
		tags$h1(tags$b("Metadata Explorer")),
		
		tags$h4("Welcome to the CMM's Metadata Explorer application!"),
		tags$h5("If you are starting afresh with a new metadata file, please go ahead and upload it below."),
		
		tags$hr(),
		
		fixedPage(
		
			fixedRow(
				column(4,
					tags$b("Upload your Metadata File:"),
					tags$style(".fa-flag {color:#00AA00}"),
					actionLink("ImportExportTab.start_here", NULL, icon=icon("flag")),
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
						column(2,
							textOutput("ImportExportTab.metadata_ready"))
					)
				)
			),#fixedRow
			
			tags$hr(),
			
			fixedRow(
				column(4,
					tags$b("Upload your Transformation File:"),
					tags$p("You may upload a previously saved Transformation File (.trn) here."),
					fileInput(
						inputId="ImportExportTab.transformation_load",
						label=NULL, #label="Upload Your Metadata File:",
						accept=c(
							".txt",
							".trn"))
						),
				column(5,
					tags$p(tags$b("Download your Transformation File:")),
					tags$p("You may download your new transformations, if you have modified it."),
					fixedRow(
						column(5,
							downloadButton(
								outputId="ImportExportTab.transformation_download",
								label="Download Transformation"
								)),
						column(2,
							textOutput("ImportExportTab.transformation_ready"))
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
						column(2,
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
						column(2,	
							textOutput("ImportExportTab.model_ready"))
					)
				)
			) #fixedRow
		) #fixedPage
	) #tabPanel
}

observe_ImportExportTabEvents=function(input, output, session){

	output$ImportExportTab.metadata_ready=renderText("Not Ready")
	output$ImportExportTab.transformation_ready=renderText("Not Ready")
	output$ImportExportTab.study_ready=renderText("Not Ready")
	output$ImportExportTab.model_ready=renderText("Not Ready")
	
	# Import handlers
	
	observeEvent(input$ImportExportTab.metadata_load, {
		cat("Metadata Load File: ", input$ImportExportTab.metadata_load$datapath, "\n");
		metadata=MetadataRec.read_metadata(input$ImportExportTab.metadata_load$datapath);
		
		#######################################################################
		
		# Update data across tabs
		output$DataTab.table=renderDT(metadata);
		
		available_variables=sort(colnames(metadata));
		updateSelectInput(session, inputId="ModelBuilderTab.available_selector",
			choices=available_variables);
			
		session$userData[["Metadata"]]=metadata;	
		session$userData[["ModelBuilderSets"]][["ModelBuilderTab.available_selector"]]=available_variables;
	
		
		#
		#######################################################################
		
		
		metadata_colnames=colnames(metadata);
		updateSelectInput(session, "DataTab.disp_col_selector", 
			choices=metadata_colnames, selected=metadata_colnames);
		
		cat("Leaving metadata load\n");
	});
	
	#observeEvent(input$ImportExportTab.transformation_load, {
	#	cat("Study Load File: ", input$ImportExportTab.study_load$datapath, "\n");
	#});

	observeEvent(input$ImportExportTab.study_load, {
		cat("Study Load File: ", input$ImportExportTab.study_load$datapath, "\n");
	});

	observeEvent(input$ImportExportTab.model_load, {
		cat("Model Load File: ", input$ImportExportTab.model_load$datapath, "\n");
	});

	# Download handlers
	
	output$ImportExportTab.metadata_download=downloadHandler(
			filename=function(){
				paste("MyMetadata_", format(Sys.time(), "%y%m%d_%H%M"), ".tsv", sep="");
			},
			content=function(filename){
				MetadataRec.write_metadata(filename);
			}
		);
		
	#output$ImportExportTab.transformation_download=downloadHandler(
	#		filename=function(){
	#			paste("MyTransformations_", format(Sys.time(), "%y%m%d_%H%M"), ".tsv", sep="");
	#		},
	#		content=function(filename){
	#			TransformationRec.write_transformations(filename);
	#		}
	#	);
		
	output$ImportExportTab.study_download=downloadHandler(
			filename=function(){
				paste("MyStudy_", format(Sys.time(), "%y%m%d_%H%M"), ".stu", sep="");
			},
			content=function(filename){
				StudyRec.write_study(filename);
			}
		);
		
		
	output$ImportExportTab.model_download=downloadHandler(
			filename=function(){
				paste("MyModel_", format(Sys.time(), "%y%m%d_%H%M"), ".mod", sep="");
			},
			content=function(filename){
				ModelRec.write_model(filename);
			}
		);
	
	cat("Leaving observe_ImportExportTabEvents\n");
}

###############################################################################

if(!exists("integration")){

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

}
